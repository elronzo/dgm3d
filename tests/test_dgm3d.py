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


# ───────────────────────────── buildings ──────────────────────────────

# a 1 km² "village" in the middle of the synthetic tile, 1 m resolution
VX0, VY0 = TILE_E + 10000, TILE_N + 50000 - 10000          # upper-left corner (EPSG:3035)
VN = 1000
HOUSES = [(200, 300, 30, 20, 9.0), (500, 500, 60, 40, 18.0), (700, 250, 15, 15, 6.0)]  # col,row,w,h,height
TREES = [(350, 700), (360, 720), (800, 800)]


def village(tmp, dsm_offset=0.0):
    """Return (dtm_1m, dsm_1m) arrays for the village: DSM = terrain + houses + trees."""
    dtm_src = write_tile(tmp / "dtm25.tif", synthetic_terrain())
    ds = gdal.Warp("", str(dtm_src), format="MEM", outputBounds=(VX0, VY0 - VN, VX0 + VN, VY0),
                   xRes=1, yRes=1, resampleAlg="bilinear")
    dtm = ds.GetRasterBand(1).ReadAsArray().astype("f4")
    dsm = dtm.copy()
    for c, r, w, h, ht in HOUSES:
        dsm[r:r + h, c:c + w] = dtm[r:r + h, c:c + w] + ht
    yy, xx = np.mgrid[0:VN, 0:VN]
    for c, r in TREES:
        dsm = np.maximum(dsm, dtm + 15 * np.exp(-((xx - c) ** 2 + (yy - r) ** 2) / 30))
    return dtm, dsm + dsm_offset


def write_raster(path, arr, x0, y0, px, epsg=3035, nodata=-9999):
    ds = gdal.GetDriverByName("GTiff").Create(str(path), arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
    ds.SetGeoTransform((x0, px, 0, y0, 0, -px))
    if epsg:
        s = osr.SpatialReference(); s.ImportFromEPSG(epsg); ds.SetProjection(s.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(arr)
    ds.GetRasterBand(1).SetNoDataValue(nodata)
    ds = None
    return path


def write_footprints(path, epsg=31287):
    """House polygons in Lambert Austria (≠ EPSG:3035, so reprojection is exercised) + a line layer."""
    from osgeo import ogr
    ogr.UseExceptions()
    drv = ogr.GetDriverByName("GPKG")
    ds = drv.CreateDataSource(str(path))
    lam = osr.SpatialReference(); lam.ImportFromEPSG(epsg)
    lam.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(dgm3d.srs(3035), lam)
    # like the real BEV DLM: a coverage polygon over everything + the building layer
    cov = ds.CreateLayer("BWK_8100_AKTUALITAETSSTAND", lam, ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for dx, dy in ((-500, 500), (1500, 500), (1500, -1500), (-500, -1500), (-500, 500)):
        ring.AddPoint_2D(*tr.TransformPoint(VX0 + dx, VY0 + dy)[:2])
    poly = ogr.Geometry(ogr.wkbPolygon); poly.AddGeometry(ring)
    f = ogr.Feature(cov.GetLayerDefn()); f.SetGeometry(poly); cov.CreateFeature(f)
    lyr = ds.CreateLayer("BWK_8100_BAUWERK_F", lam, ogr.wkbPolygon)
    for c, r, w, h, _ in HOUSES:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for dx, dy in ((0, 0), (w, 0), (w, h), (0, h), (0, 0)):
            x, y, _ = tr.TransformPoint(VX0 + c + dx, VY0 - r - dy)
            ring.AddPoint_2D(x, y)
        poly = ogr.Geometry(ogr.wkbPolygon); poly.AddGeometry(ring)
        f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(poly); lyr.CreateFeature(f)
    lines = ds.CreateLayer("bauwerk_linie", lam, ogr.wkbLineString)
    ln = ogr.Geometry(ogr.wkbLineString)
    ln.AddPoint_2D(*tr.TransformPoint(VX0, VY0 - 10)[:2]); ln.AddPoint_2D(*tr.TransformPoint(VX0 + 900, VY0 - 10)[:2])
    f = ogr.Feature(lines.GetLayerDefn()); f.SetGeometry(ln); lines.CreateFeature(f)
    ds = None
    return path


def mask_at(path, c, r):
    """Value of a mask raster at village pixel (col c, row r), via the raster's georeference."""
    d = dgm3d.DEM(str(path))
    x, y = VX0 + c + 0.5, VY0 - r - 0.5
    return d.z[int((y - d.gt[3]) / d.gt[5]), int((x - d.gt[0]) / d.gt[1])] > 0.5


def village_bbox():
    w, s = to_lonlat(VX0 + 100, VY0 - VN + 100)
    e, n = to_lonlat(VX0 + VN - 100, VY0 - 100)
    return w, s, e, n


@pytest.fixture
def village_server(tmp_path, monkeypatch):
    srv = tmp_path / "srv"; srv.mkdir()
    dtm, dsm = village(tmp_path)
    write_raster(srv / f"ALS_DTM_CRS3035RES50000mN{TILE_N}E{TILE_E}.tif", dtm, VX0, VY0, 1)
    write_raster(srv / f"ALS_DSM_CRS3035RES50000mN{TILE_N}E{TILE_E}.tif", dsm, VX0, VY0, 1)
    write_footprints(srv / "bauwerke.gpkg")
    handler = lambda *a, **k: RangeHandler(*a, directory=str(srv), **k)  # noqa: E731
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    for var in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("NO_PROXY", "127.0.0.1")
    base = f"http://127.0.0.1:{httpd.server_port}"
    monkeypatch.setattr(dgm3d, "BASE_URL", base)
    yield base, dtm, dsm
    httpd.shutdown()


def test_buildings_from_bev(village_server, tmp_path, capsys):
    base, dtm, dsm = village_server
    out = tmp_path / "dorf.stl"
    dgm3d.main(["fetch", *map(str, village_bbox()), str(out), "--dem-res", "1", "--width", "100",
                "--res", "0.25", "--exag", "1", "--base", "2", "--buildings",
                "--footprints", f"{base}/bauwerke.gpkg"])
    err = capsys.readouterr().err
    assert "layer BWK_8100_BAUWERK_F: 3 polygons" in err
    assert "AKTUALITAETSSTAND" not in err
    whole = assert_closed(read_stl(out))
    terrain = assert_closed(read_stl(tmp_path / "dorf_terrain.stl"))
    houses = assert_closed(read_stl(tmp_path / "dorf_buildings.stl"))
    assert terrain + houses == pytest.approx(whole, rel=1e-3)
    # the tallest house (18 m above ground) shows up …
    assert abs(float(re.search(r"buildings: up to (\d+) m", err).group(1)) - 18) <= 1
    mask = dgm3d.DEM(str(tmp_path / "dorf_bldg.tif")).z > 0.5
    assert mask.sum() == pytest.approx(sum(w * h for _, _, w, h, _ in HOUSES), rel=0.25)
    # … the 15 m trees do not; the house centres do
    for c, r in TREES:
        assert not mask_at(tmp_path / "dorf_bldg.tif", c, r)
    for c, r, w, h, _ in HOUSES:
        assert mask_at(tmp_path / "dorf_bldg.tif", c + w // 2, r + h // 2)


def test_vienna_local_dsm_with_offset(village_server, tmp_path, capsys):
    """Vienna-style DSM: other CRS (EPSG:31256), 0.5 m, heights above 'Wiener Null' (−156.68 m),
    file without CRS → --dsm-srs; the tool must find the offset itself."""
    base, dtm, dsm = village_server
    # resample the true DSM to 0.5 m in EPSG:31256, shift by the datum difference
    src = write_raster(tmp_path / "dsm3035.tif", dsm, VX0, VY0, 1)
    gk = gdal.Warp("", str(src), format="MEM", dstSRS="EPSG:31256", xRes=0.5, yRes=0.5,
                   resampleAlg="bilinear", dstNodata=-9999)
    arr = gk.GetRasterBand(1).ReadAsArray()
    arr = np.where(arr == -9999, -9999, arr - 156.68).astype("f4")
    gt = gk.GetGeoTransform()
    tiles = []
    h = arr.shape[0] // 2                    # split in two "sheets" like the Vienna download
    for i, part in enumerate((arr[:h], arr[h:])):
        tiles.append(str(write_raster(tmp_path / f"wien_dom_{i}.tif", part, gt[0], gt[3] - i * h * 0.5,
                                      0.5, epsg=None)))
    out = tmp_path / "wien.stl"
    dgm3d.main(["fetch", *map(str, village_bbox()), str(out), "--dem-res", "1", "--width", "100",
                "--res", "0.25", "--exag", "1", "--buildings", "--dsm-files", *tiles,
                "--dsm-srs", "EPSG:31256", "--footprints", f"{base}/bauwerke.gpkg"])
    err = capsys.readouterr().err
    off = float(re.search(r"height offset: ([+-][\d.]+) m", err).group(1))
    assert off == pytest.approx(156.68, abs=0.3)
    assert abs(float(re.search(r"buildings: up to (\d+) m", err).group(1)) - 18) <= 1
    for part in ("", "_terrain", "_buildings"):
        assert_closed(read_stl(tmp_path / f"wien{part}.stl"))


def test_stl_reuses_building_files(village_server, tmp_path):
    base, _, _ = village_server
    tif = tmp_path / "d.tif"
    dgm3d.main(["fetch", *map(str, village_bbox()), str(tif), "--dem-res", "1", "--buildings",
                "--footprints", f"{base}/bauwerke.gpkg"])
    dgm3d.main(["stl", str(tif), str(tmp_path / "d2.stl"), "--width", "80", "--res", "0.4", "--buildings"])
    assert_closed(read_stl(tmp_path / "d2_buildings.stl"))


def test_footprints_none_includes_trees(village_server, tmp_path):
    dgm3d.main(["fetch", *map(str, village_bbox()), str(tmp_path / "n.stl"), "--dem-res", "1",
                "--width", "80", "--res", "0.4", "--buildings", "--footprints", "none"])
    for c, r in TREES:
        assert mask_at(tmp_path / "n_bldg.tif", c, r)    # trees count as "buildings" here
    assert_closed(read_stl(tmp_path / "n_buildings.stl"))


def test_dsm_files_missing_and_zip(tmp_path):
    with pytest.raises(SystemExit, match="Geodatenviewer"):
        dgm3d.resolve_raster_files([str(tmp_path / "wien_dom/*.tif")])
    import zipfile
    t = write_raster(tmp_path / "sheet.tif", np.ones((10, 10), "f4"), 0, 10, 1, epsg=None)
    with zipfile.ZipFile(tmp_path / "sheet.zip", "w") as z:
        z.write(t, "dom/sheet.tif")
    files = dgm3d.resolve_raster_files([str(tmp_path / "*.zip")])
    assert files == [f"/vsizip/{tmp_path}/sheet.zip/dom/sheet.tif"]
    assert gdal.Open(files[0]).RasterXSize == 10
