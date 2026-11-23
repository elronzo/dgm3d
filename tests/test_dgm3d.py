"""Offline tests for dgm3d — synthetic terrain, no network access to BEV needed.

Run:  python3 -m pytest -q
"""
import http.server, math, os, re, struct, sys, threading
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, osr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import dgm3d  # noqa: E402

gdal.UseExceptions()

# A 20 × 20 km synthetic tile on the BEV grid, near the Großglockner (EPSG:3035).
TILE_E, TILE_N = 4500000, 2650000
PX = 25.0
N = 800


def synthetic_terrain(seed=3):
    rng = np.random.default_rng(seed)
    f = np.fft.fftfreq(N)
    fx, fy = np.meshgrid(f, f)
    k = np.hypot(fx, fy)
    k[0, 0] = 1
    spec = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / k ** 2 * np.exp(-(k / 0.1) ** 2)
    spec[0, 0] = 0
    z = np.real(np.fft.ifft2(spec))
    z = (z - z.min()) / (z.max() - z.min())
    return (800 + 2800 * z ** 1.6).astype("f4")


def write_tile(path, z, nodata=None, overviews=False):
    ds = gdal.GetDriverByName("GTiff").Create(str(path), z.shape[1], z.shape[0], 1, gdal.GDT_Float32,
                                              ["TILED=YES", "COMPRESS=DEFLATE"])
    # tile occupies the upper-left 20 km of the 50 km cell
    ds.SetGeoTransform((TILE_E, PX, 0, TILE_N + 50000, 0, -PX))
    s = osr.SpatialReference()
    s.ImportFromEPSG(3035)
    ds.SetProjection(s.ExportToWkt())
    b = ds.GetRasterBand(1)
    b.WriteArray(z)
    if nodata is not None:
        b.SetNoDataValue(nodata)
    if overviews:
        ds.BuildOverviews("AVERAGE", [2, 4, 8])
    ds = None
    return path


def to_lonlat(x, y):
    t = osr.CoordinateTransformation(dgm3d.srs(3035), dgm3d.srs(4326))
    lon, lat, _ = t.TransformPoint(x, y)
    return lon, lat


def read_stl(path):
    with open(path, "rb") as f:
        f.read(80)
        n = struct.unpack("<I", f.read(4))[0]
        rec = np.fromfile(f, dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")], count=n)
    assert len(rec) == n
    return rec["v"].astype(np.float64)


def assert_closed(tris):
    """Every directed edge exactly once and its reverse exactly once → closed, consistently oriented,
    edge-manifold. Also checks for positive volume (normals point outwards)."""
    assert len(tris)
    _, idx = np.unique(np.round(tris.reshape(-1, 3), 4), axis=0, return_inverse=True)
    t = idx.reshape(-1, 3)
    e = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    e = e[e[:, 0] != e[:, 1]]
    key = e[:, 0].astype(np.int64) * (idx.max() + 1) + e[:, 1]
    rkey = e[:, 1].astype(np.int64) * (idx.max() + 1) + e[:, 0]
    u, c = np.unique(key, return_counts=True)
    assert c.max() == 1, "an edge is used twice in the same direction (non-manifold or flipped)"
    assert np.isin(rkey, u).all(), "open edge (hole in the mesh)"
    vol = np.einsum("ij,ij->i", tris[:, 0], np.cross(tris[:, 1], tris[:, 2])).sum() / 6
    assert vol > 0, "mesh is inside out"
    return vol


# ───────────────────────────── mesh ──────────────────────────────

def test_grid_mesh_rectangle_closed():
    rng = np.random.default_rng(0)
    top = 2 + rng.random((40, 60)) * 5
    vol = assert_closed(dgm3d.grid_mesh(top, 0.5, 0.5).astype(np.float64))
    assert vol == pytest.approx(top[:-1, :-1].mean() * 59 * 0.5 * 39 * 0.5, rel=0.05)


@pytest.mark.parametrize("seed", range(5))
def test_grid_mesh_random_mask_closed(seed):
    rng = np.random.default_rng(seed)
    top = 2 + rng.random((50, 50)) * 3
    mask = rng.random((49, 49)) > 0.45            # lots of diagonal contacts
    mask = dgm3d.no_diagonals(mask)
    assert_closed(dgm3d.grid_mesh(top, 0.3, 0.3, mask).astype(np.float64))


def test_no_diagonals_removes_corner_contacts():
    m = np.array([[1, 0], [0, 1]], bool)
    out = dgm3d.no_diagonals(m)
    assert out.sum() == 3


# ───────────────────────────── stl ──────────────────────────────

@pytest.fixture(scope="module")
def dem_tif(tmp_path_factory):
    return write_tile(tmp_path_factory.mktemp("dem") / "dem.tif", synthetic_terrain())


def test_stl_dimensions_and_closed(dem_tif, tmp_path):
    out = tmp_path / "a.stl"
    dgm3d.dem_to_stl(str(dem_tif), str(out), width=100, exag=1.0, base=2, res=0.5)
    tris = read_stl(out)
    assert_closed(tris)
    ext = tris.reshape(-1, 3).max(0) - tris.reshape(-1, 3).min(0)
    assert ext[0] == pytest.approx(100, abs=0.6)
    # 2800 m of relief at 1:200,000 = 14 mm, plus 2 mm base
    assert ext[2] == pytest.approx(2 + 2800 / 200_000 * 1000, rel=0.08)


def test_gamma_keeps_peak_and_lifts_lowlands(dem_tif, tmp_path):
    a, b = tmp_path / "g1.stl", tmp_path / "g06.stl"
    dgm3d.dem_to_stl(str(dem_tif), str(a), width=80, exag=3, base=2, res=0.8, gamma=1.0)
    dgm3d.dem_to_stl(str(dem_tif), str(b), width=80, exag=3, base=2, res=0.8, gamma=0.6)
    za, zb = read_stl(a)[..., 2], read_stl(b)[..., 2]
    assert zb.max() == pytest.approx(za.max(), abs=1e-3)
    assert zb[zb > 2].mean() > za[za > 2].mean()


def test_crop_follows_data_outline(tmp_path):
    z = synthetic_terrain()
    yy, xx = np.mgrid[0:N, 0:N] / N
    z[xx < 0.3 + 0.1 * np.sin(yy * 12)] = -9999          # wiggly "border"
    z[(xx > 0.8) & (yy > 0.8)] = -9999
    tif = write_tile(tmp_path / "border.tif", z, nodata=-9999)
    out = tmp_path / "c.stl"
    dgm3d.dem_to_stl(str(tif), str(out), width=100, exag=1.5, base=2, res=0.5, crop=True)
    tris = read_stl(out)
    assert_closed(tris)
    xs = tris[..., 0]
    assert xs.min() > 15          # the western ~20–40 % is gone


def test_bands_printed(dem_tif, tmp_path, capsys):
    dgm3d.dem_to_stl(str(dem_tif), str(tmp_path / "b.stl"), width=80, exag=2, base=2, res=1.0,
                     bands=[1500, 9000], layer=0.1)
    err = capsys.readouterr().err
    lines = [l for l in err.splitlines() if "→" in l and "mm" in l and " m " in l]
    assert len(lines) == 2
    mm = float(re.search(r"→\s+([\d.]+) mm", lines[0]).group(1))
    assert mm == pytest.approx(round(mm, 1))
    assert "above the highest point" in lines[1]


# ───────────────────────────── view ──────────────────────────────

def test_view_all_modes(dem_tif, tmp_path):
    # standpoint: on a slope 6 km west of the highest summit, looking straight at it
    z = synthetic_terrain()
    r, c = np.unravel_index(np.argmax(z[150:650, 350:700]), (500, 350))
    r, c = r + 150, c + 350
    x, y = TILE_E + (c + 0.5) * PX - 6000, TILE_N + 50000 - (r + 0.5) * PX
    lon, lat = to_lonlat(x, y)
    prefix = tmp_path / "v"
    dgm3d.main(["view", "--dem", str(dem_tif), "--lat", str(lat), "--lon", str(lon),
                "--az", "90", "--fov", "60", "--dist", "8", "--mode", "all", "--out", str(prefix),
                "--width-wedge", "60", "--res", "0.6", "--width", "60", "--res-view", "0.4",
                "--steps", "1500"])
    for part in ("wedge", "visible", "hidden", "visible_engraved", "relief", "litho"):
        assert_closed(read_stl(f"{prefix}_{part}.stl"))
    for png in ("wedge", "visible", "view", "relief", "litho"):
        assert os.path.getsize(f"{prefix}_{png}.png") > 200
    # visible + hidden make up the whole wedge
    v = sum(assert_closed(read_stl(f"{prefix}_{p}.stl")) for p in ("visible", "hidden"))
    assert v == pytest.approx(assert_closed(read_stl(f"{prefix}_wedge.stl")), rel=1e-3)


def test_sector_bbox_contains_standpoint():
    w, s, e, n = dgm3d.sector_bbox_lonlat(47.07, 12.75, 268, 60, 12000)
    assert w < 12.75 < e and s < 47.07 < n
    assert (12.75 - w) > (e - 12.75)          # looking west → bbox extends west


# ───────────────────────────── fetch ──────────────────────────────

class RangeHandler(http.server.SimpleHTTPRequestHandler):
    """Minimal HTTP server with Range support, as GDAL's /vsicurl/ needs."""

    def log_message(self, *a):
        pass

    def send_head(self):
        path = self.translate_path(self.path)
        if not os.path.isfile(path):
            self.send_error(404)
            return None
        size = os.path.getsize(path)
        rng = self.headers.get("Range")
        f = open(path, "rb")
        if rng:
            a, b = rng.replace("bytes=", "").split("-")
            a, b = int(a), min(int(b) if b else size - 1, size - 1)
            f.seek(a)
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {a}-{b}/{size}")
            self.send_header("Content-Length", str(b - a + 1))
            self.end_headers()
            self._remaining = b - a + 1
        else:
            self.send_response(200)
            self.send_header("Content-Length", str(size))
            self.send_header("Accept-Ranges", "bytes")
            self.end_headers()
            self._remaining = size
        return f

    def copyfile(self, src, dst):
        dst.write(src.read(self._remaining))


@pytest.fixture
def tile_server(tmp_path, monkeypatch):
    srv_dir = tmp_path / "srv"
    srv_dir.mkdir()
    write_tile(srv_dir / f"ALS_DTM_CRS3035RES50000mN{TILE_N}E{TILE_E}.tif", synthetic_terrain(),
               nodata=-9999, overviews=True)
    handler = lambda *a, **k: RangeHandler(*a, directory=str(srv_dir), **k)  # noqa: E731
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    for var in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("NO_PROXY", "127.0.0.1")
    monkeypatch.setattr(dgm3d, "BASE_URL", f"http://127.0.0.1:{httpd.server_port}")
    yield
    httpd.shutdown()


def test_fetch_to_stl(tile_server, tmp_path):
    x0, y0 = TILE_E + 5000, TILE_N + 50000 - 15000
    w, s = to_lonlat(x0, y0)
    e, n = to_lonlat(x0 + 6000, y0 + 4000)
    out = tmp_path / "f.stl"
    dgm3d.main(["fetch", str(w), str(s), str(e), str(n), str(out),
                "--dem-res", "25", "--width", "90", "--res", "0.5"])
    assert (tmp_path / "f.tif").exists()
    assert_closed(read_stl(out))
