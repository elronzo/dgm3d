#!/usr/bin/env python3
"""dgm3d.py — printable 3D terrain from Austria's 1 m ALS terrain model (BEV DGM).

Subcommands
  fetch  W S E N out.tif|out.stl   cut a WGS84 bbox from the BEV tiles (remote, no full download);
                                   an .stl output also builds the mesh
  stl    in.tif out.stl            any projected DEM GeoTIFF (metres) -> watertight STL
  view   --lat --lon --az ...      prints of the view from a standpoint, four modes:
           wedge     pie-slice of terrain from the standpoint along the view direction
           visible   same wedge split into visible/hidden parts (two STLs for multi-colour)
                     plus a single-colour version with hidden terrain sunk in
           relief    upright bas-relief plaque of the panorama (near ridges stand out more)
           litho     lithophane of the rendered panorama (hold against a window)
           all       all of the above

Examples
  ./dgm3d.py fetch 16.30 48.25 16.38 48.29 kahlenberg.stl --dem-res 2 --width 180 --exag 1.5
  ./dgm3d.py fetch 9.53 46.37 17.17 49.02 austria.stl --dem-res 250 --width 175 --exag 4 \
             --gamma 0.7 --crop --bands 1000,2500
  ./dgm3d.py view --lat 47.0748 --lon 12.7506 --az 268 --fov 60 --dist 12 --mode all --out glockner

Requires numpy and the GDAL Python bindings (Fedora: sudo dnf install python3-gdal python3-numpy);
scipy is optional (removes tiny islands with --crop).

Code:  PolyForm Noncommercial 1.0.0 — https://github.com/elronzo/dgm3d
Data:  © BEV – Bundesamt für Eich- und Vermessungswesen, CC BY 4.0 — https://data.bev.gv.at
"""
import argparse, math, os, struct, sys
import numpy as np
from osgeo import gdal, osr

__version__ = "0.1.0"

gdal.UseExceptions()

STAND = os.environ.get("STAND", "20230915")
BASE_URL = os.environ.get("BASE", f"https://data.bev.gv.at/download/ALS/DTM/{STAND}")
TILE = 50_000
R_EARTH, REFRACTION = 6_371_000.0, 0.13

for k, v in {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
             "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
             "VSI_CACHE": "TRUE", "VSI_CACHE_SIZE": "200000000",
             "GDAL_HTTP_MULTIRANGE": "YES",
             "GDAL_HTTP_MERGE_CONSECUTIVE_RANGES": "YES"}.items():
    gdal.SetConfigOption(k, v)


def log(*a):
    print(*a, file=sys.stderr)


# ───────────────────────────── coordinates ──────────────────────────────

def srs(defn):
    s = osr.SpatialReference()
    if isinstance(defn, int):
        s.ImportFromEPSG(defn)
    else:
        s.ImportFromWkt(defn)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)   # always x/lon first
    return s


def lonlat_to(target, pts):
    tr = osr.CoordinateTransformation(srs(4326), srs(target))
    return np.array([tr.TransformPoint(float(lon), float(lat))[:2] for lon, lat in pts])


# ───────────────────────────── fetch ──────────────────────────────

def fetch(w, s, e, n, out_tif, res=1.0):
    """Clip a WGS84 bbox out of the BEV 1 m DTM tiles into a GeoTIFF (EPSG:3035)."""
    xy = lonlat_to(3035, [(w, s), (w, n), (e, s), (e, n)])
    xmin, ymin = np.floor(xy.min(0)).astype(int)
    xmax, ymax = np.ceil(xy.max(0)).astype(int)
    log(f"bbox EPSG:3035  E {xmin}–{xmax}  N {ymin}–{ymax}  ({xmax-xmin} × {ymax-ymin} m)")

    # probe one ring of neighbouring 50 km tiles (works for either corner convention)
    inputs = []
    for ty in range(ymin // TILE * TILE - TILE, ymax // TILE * TILE + TILE + 1, TILE):
        for tx in range(xmin // TILE * TILE - TILE, xmax // TILE * TILE + TILE + 1, TILE):
            for name in (f"ALS_DTM_CRS3035RES50000mN{ty}E{tx}.tif", f"CRS3035RES50000mN{ty}E{tx}.tif"):
                url = f"/vsicurl/{BASE_URL}/{name}"
                try:
                    ds = gdal.Open(url)
                except RuntimeError:
                    continue
                gt = ds.GetGeoTransform()
                txmin, tymax = gt[0], gt[3]
                txmax, tymin = txmin + gt[1] * ds.RasterXSize, tymax + gt[5] * ds.RasterYSize
                if txmin < xmax and txmax > xmin and tymin < ymax and tymax > ymin:
                    ovr = ds.GetRasterBand(1).GetOverviewCount()
                    log(f"  tile {name}  ({ovr} overviews)")
                    if res >= 20 and ovr == 0:
                        log("  warning: no overviews — a coarse fetch must read every 1 m pixel; this will be slow")
                    inputs.append(url)
                break
    if not inputs:
        sys.exit(f"No BEV tile found — bbox outside Austria, or no release at {BASE_URL}/ "
                 "(set STAND=YYYYMMDD for a newer one)")

    gdal.Warp(out_tif, inputs, dstSRS="EPSG:3035", dstNodata=-9999,
              outputBounds=(xmin, ymin, xmax, ymax), xRes=res, yRes=res, targetAlignedPixels=True,
              resampleAlg="bilinear" if res <= 1 else "average", multithread=True,
              warpOptions=["NUM_THREADS=ALL_CPUS"],
              creationOptions=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"])
    log(f"GeoTIFF → {out_tif}")
    return out_tif


# ───────────────────────────── DEM ──────────────────────────────

class DEM:
    def __init__(self, path, max_px=None):
        ds = gdal.Open(path)
        if max_px and max(ds.RasterXSize, ds.RasterYSize) > max_px:
            f = max_px / max(ds.RasterXSize, ds.RasterYSize)
            ds = gdal.Translate("", ds, format="MEM", resampleAlg="average",
                                width=max(2, round(ds.RasterXSize * f)),
                                height=max(2, round(ds.RasterYSize * f)))
        b = ds.GetRasterBand(1)
        z = b.ReadAsArray().astype(np.float64)
        nd = b.GetNoDataValue()
        if nd is not None:
            z[np.isclose(z, nd)] = np.nan
        z[(z < -1000) | (z > 9000)] = np.nan
        self.z, self.gt, self.wkt = z, ds.GetGeoTransform(), ds.GetProjection()
        self.rows, self.cols = z.shape
        self.px = abs(self.gt[1])

    @property
    def extent_m(self):
        return abs(self.gt[1]) * self.cols, abs(self.gt[5]) * self.rows

    def sample(self, x, y, grid=None):
        """Bilinear sample at map coords (NaN outside)."""
        g = self.z if grid is None else grid
        fc = (np.asarray(x) - self.gt[0]) / self.gt[1] - 0.5
        fr = (np.asarray(y) - self.gt[3]) / self.gt[5] - 0.5
        c0, r0 = np.floor(fc).astype(np.int64), np.floor(fr).astype(np.int64)
        ok = (c0 >= 0) & (r0 >= 0) & (c0 < self.cols - 1) & (r0 < self.rows - 1)
        c0c, r0c = np.clip(c0, 0, self.cols - 2), np.clip(r0, 0, self.rows - 2)
        tx, ty = fc - c0, fr - r0
        v = (g[r0c, c0c] * (1 - tx) * (1 - ty) + g[r0c, c0c + 1] * tx * (1 - ty) +
             g[r0c + 1, c0c] * (1 - tx) * ty + g[r0c + 1, c0c + 1] * tx * ty)
        return np.where(ok, v, np.nan)

    def hillshade(self, sun_az, sun_alt):
        z = np.where(np.isnan(self.z), np.nanmin(self.z), self.z)
        dzdx = np.gradient(z, axis=1) / self.gt[1]
        dzdy = np.gradient(z, axis=0) / self.gt[5]          # northward derivative
        n = np.stack([-dzdx, -dzdy, np.ones_like(z)])
        n /= np.linalg.norm(n, axis=0)
        a, h = math.radians(sun_az), math.radians(sun_alt)
        L = np.array([math.sin(a) * math.cos(h), math.cos(a) * math.cos(h), math.sin(h)])
        return np.clip(np.tensordot(L, n, 1), 0, 1)


# ───────────────────────────── mesh + STL ──────────────────────────────

def grid_mesh(top, sx, sy, cellmask=None):
    """Watertight solid from a vertex height grid (mm). cellmask selects cells (rows-1, cols-1);
    None = whole rectangle. Bottom at z=0, north (row 0) = +Y."""
    rows, cols = top.shape
    X, Y = np.meshgrid(np.arange(cols) * sx, (rows - 1 - np.arange(rows)) * sy)
    T = np.stack([X, Y, top], -1)
    B = np.stack([X, Y, np.zeros_like(top)], -1)
    tris = []

    if cellmask is None or cellmask.all():
        a, b, c, d = T[:-1, :-1], T[:-1, 1:], T[1:, :-1], T[1:, 1:]
        tris += [np.stack([a, c, b], -2).reshape(-1, 3, 3), np.stack([b, c, d], -2).reshape(-1, 3, 3)]
        r, q = rows - 1, cols - 1
        loop = np.array([(r, j) for j in range(q)] + [(i, q) for i in range(r, 0, -1)] +
                        [(0, j) for j in range(q, 0, -1)] + [(i, 0) for i in range(r)])
        nxt = np.roll(loop, -1, 0)
        t0, t1 = T[loop[:, 0], loop[:, 1]], T[nxt[:, 0], nxt[:, 1]]
        b0, b1 = B[loop[:, 0], loop[:, 1]], B[nxt[:, 0], nxt[:, 1]]
        tris += [np.stack([b0, b1, t1], 1), np.stack([b0, t1, t0], 1)]
        centre = np.broadcast_to([X.max() / 2, Y.max() / 2, 0.0], b0.shape)
        tris.append(np.stack([centre, b1, b0], 1))
        return np.concatenate(tris).astype(np.float32)

    r, c = np.nonzero(cellmask)
    a, b, cc, d = T[r, c], T[r, c + 1], T[r + 1, c], T[r + 1, c + 1]
    tris += [np.stack([a, cc, b], 1), np.stack([b, cc, d], 1)]
    A, Bb, C, D = B[r, c], B[r, c + 1], B[r + 1, c], B[r + 1, c + 1]
    tris += [np.stack([A, Bb, C], 1), np.stack([Bb, D, C], 1)]
    P = np.pad(cellmask, 1)
    # (neighbour present?, edge start, edge end) — edges run counter-clockwise seen from above
    for nb, p0, p1 in ((P[r, c + 1],     (r, c + 1),     (r, c)),           # north
                       (P[r + 2, c + 1], (r + 1, c),     (r + 1, c + 1)),   # south
                       (P[r + 1, c + 2], (r + 1, c + 1), (r, c + 1)),       # east
                       (P[r + 1, c],     (r, c),         (r + 1, c))):      # west
        o = ~nb
        t0, t1 = T[p0[0][o], p0[1][o]], T[p1[0][o], p1[1][o]]
        b0, b1 = B[p0[0][o], p0[1][o]], B[p1[0][o], p1[1][o]]
        tris += [np.stack([b0, b1, t1], 1), np.stack([b0, t1, t0], 1)]
    return np.concatenate(tris).astype(np.float32)


def no_diagonals(m, keep=None):
    """Remove cells touching only at a corner (would give non-manifold edges) by filling one
    off-diagonal cell; 'keep' limits where cells may be added."""
    m = m.copy()
    for _ in range(50):
        a, b, c, d = m[:-1, :-1], m[:-1, 1:], m[1:, :-1], m[1:, 1:]
        p1, p2 = a & d & ~b & ~c, b & c & ~a & ~d
        if keep is not None:
            kb, ka = keep[:-1, 1:], keep[:-1, :-1]
            kc, kd = keep[1:, :-1], keep[1:, 1:]
        else:
            ka = kb = kc = kd = np.ones_like(p1)
        if not (p1.any() or p2.any()):
            break
        f = np.zeros_like(m)
        f[:-1, 1:] |= p1 & kb; f[1:, :-1] |= p1 & ~kb & kc
        f[:-1, :-1] |= p2 & ka; f[1:, 1:] |= p2 & ~ka & kd
        if not f.any():
            break
        m |= f
    return m


def drop_specks(m, min_cells=50):
    """Remove small disconnected islands (e.g. enclaves of data along the border)."""
    try:
        from scipy import ndimage
    except ImportError:
        return m
    lab, n = ndimage.label(m)
    if n <= 1:
        return m
    sizes = np.bincount(lab.ravel())
    keep = sizes >= min_cells
    keep[0] = False
    return keep[lab]


def write_stl(path, tris):
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.divide(n, ln, out=np.zeros_like(n), where=ln > 0)
    rec = np.zeros(len(tris), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, tris
    with open(path, "wb") as f:
        f.write(b"dgm3d".ljust(80, b" "))
        f.write(struct.pack("<I", len(tris)))
        rec.tofile(f)
    log(f"  STL  {len(tris):,} triangles, {os.path.getsize(path)/1e6:.1f} MB → {path}")


def write_png(path, img):
    """img: HxW float 0..1 (grey) or HxWx3."""
    a = (np.clip(np.nan_to_num(img, nan=1.0), 0, 1) * 255).astype(np.uint8)
    if a.ndim == 2:
        a = a[..., None]
    mem = gdal.GetDriverByName("MEM").Create("", a.shape[1], a.shape[0], a.shape[2], gdal.GDT_Byte)
    for i in range(a.shape[2]):
        mem.GetRasterBand(i + 1).WriteArray(a[..., i])
    gdal.GetDriverByName("PNG").CreateCopy(path, mem)
    for ext in (".aux.xml",):
        if os.path.exists(path + ext):
            os.remove(path + ext)
    log(f"  PNG  → {path}")


def shade_heights(z, mask=None, az=315, alt=45, scale=1.0):
    dzdx, dzdy = np.gradient(z * scale)
    n = np.stack([-dzdx, dzdy, np.ones_like(z)])
    n /= np.linalg.norm(n, axis=0)
    a, h = math.radians(az), math.radians(alt)
    L = np.array([math.sin(a) * math.cos(h), math.cos(a) * math.cos(h), math.sin(h)])
    s = 0.25 + 0.75 * np.clip(np.tensordot(L, n, 1), 0, 1)
    return s if mask is None else np.where(mask, s, 1.0)


# ───────────────────────────── stl ──────────────────────────────

def dem_to_stl(dem_path, out, width=150, exag=1.5, base=3, res=0.25, zmin=None, crop=False,
               bands=None, layer=0.1, gamma=1.0):
    dem = DEM(dem_path, max_px=int(width / res) + 1)
    z = dem.z
    valid = ~np.isnan(z)
    if not valid.any():
        sys.exit("DEM contains only nodata")
    cellmask = None
    if crop:
        # keep cells whose four corners all have data → model follows the data outline (e.g. the border)
        cellmask = valid[:-1, :-1] & valid[1:, 1:] & valid[:-1, 1:] & valid[1:, :-1]
        cellmask = drop_specks(no_diagonals(cellmask))
        log(f"  cropped to data outline: {100 * cellmask.mean():.0f} % of the bbox")
    elif (~valid).mean() > 0.01:
        log(f"  warning: {100*(~valid).mean():.0f} % nodata (outside Austria?) — printed flat at the lowest "
            "height; use --crop to cut along the data outline")
    zmin_data = np.nanmin(z)
    z = np.where(valid, z, zmin_data)
    wm, hm = dem.extent_m
    mm_per_m = width / max(wm, hm)
    zmin = z.min() if zmin is None else zmin
    span = max(z.max() - zmin, 1e-6)

    def lift(h):
        """Height above zmin in m; gamma < 1 lifts low terrain, the highest point stays put."""
        rel = np.clip((np.asarray(h, float) - zmin) / span, 0, 1)
        return span * rel ** gamma

    top = base + lift(z) * mm_per_m * exag
    sx, sy = wm * mm_per_m / (dem.cols - 1), hm * mm_per_m / (dem.rows - 1)
    log(f"  terrain {wm/1000:.2f} × {hm/1000:.2f} km, {z.min():.0f}–{z.max():.0f} m")
    log(f"  print   {wm*mm_per_m:.1f} × {hm*mm_per_m:.1f} × {top.max():.1f} mm "
        f"(1:{1000/mm_per_m:,.0f}, exag {exag}×" + (f", gamma {gamma}" if gamma != 1 else "") + ")")
    write_stl(out, grid_mesh(top, sx, sy, cellmask))
    if bands:
        log(f"  filament changes (layer {layer} mm) — add in the slicer at these heights:")
        for b in bands:
            h = base + float(lift(b)) * mm_per_m * exag
            note = "  (above the highest point — skip)" if b > z.max() else ""
            log(f"    {b:6.0f} m  →  {round(h / layer) * layer:6.2f} mm{note}")


# ───────────────────────────── view ──────────────────────────────

def sector_bbox_lonlat(lat, lon, az, fov, dist_m, margin=0.03):
    """WGS84 bbox of the view sector (for fetching)."""
    angs = np.radians(az + np.linspace(-fov / 2, fov / 2, 90))
    dx = np.concatenate([[0], dist_m * np.sin(angs)])
    dy = np.concatenate([[0], dist_m * np.cos(angs)])
    pad = dist_m * margin
    m_lat, m_lon = 111_320.0, 111_320.0 * math.cos(math.radians(lat))
    return (lon + (dx.min() - pad) / m_lon, lat + (dy.min() - pad) / m_lat,
            lon + (dx.max() + pad) / m_lon, lat + (dy.max() + pad) / m_lat)


def cast(dem, x0, y0, eye_z, az, fov, near, dist, ncols, nsteps):
    """March rays across the field of view. Returns per-sample distances, ray azimuths,
    sample coordinates and the running max of elevation angle (the 'horizon so far')."""
    azr = np.radians(az + np.linspace(-fov / 2, fov / 2, ncols))
    d = np.geomspace(near, dist, nsteps)
    X = x0 + d[:, None] * np.sin(azr)[None, :]
    Y = y0 + d[:, None] * np.cos(azr)[None, :]
    h = dem.sample(X, Y)
    drop = d ** 2 / (2 * R_EARTH) * (1 - REFRACTION)            # curvature + refraction
    ang = np.arctan2(h - drop[:, None] - eye_z, d[:, None])
    ang = np.where(np.isnan(ang), -np.pi / 2, ang)
    return d, X, Y, np.maximum.accumulate(ang, axis=0), ang


def panorama(dem, x0, y0, eye_z, a):
    """Render the view: depth (m, inf = sky) and shading per pixel."""
    ncols = int(a.width / a.res_view)
    near = max(dem.px, 10.0)
    d, X, Y, cm, ang = cast(dem, x0, y0, eye_z, a.az, a.fov, near, a.dist * 1000, ncols, a.steps)
    top = cm[-1].max()
    if a.bottom is None:            # crop the foreground: horizon reached after ~8 % of the distance
        k = min(int(np.searchsorted(d, 0.08 * d[-1])), len(d) - 1)
        bottom = cm[k].min() - 0.03 * (top - cm[k].min())
    else:
        bottom = math.radians(a.bottom)
    top += (top - bottom) * a.sky
    vfrac = (top - bottom) / math.radians(a.fov)
    if a.vstretch is None:                                       # aim for a ~3:1 plaque
        a.vstretch = round(float(np.clip(0.33 / vfrac, 1.0, 4.0)), 2)
    nrows = max(8, int(round(ncols * vfrac * a.vstretch)))
    rowang = np.linspace(bottom, top, nrows)                     # ascending
    idx = np.empty((nrows, ncols), np.int64)
    for c in range(ncols):
        idx[:, c] = np.searchsorted(cm[:, c], rowang)
    idx = idx[::-1]                                              # row 0 = top
    sky = idx >= len(d)
    ic = np.minimum(idx, len(d) - 1)
    ip = np.maximum(ic - 1, 0)
    cols = np.broadcast_to(np.arange(ncols), idx.shape)
    # interpolate the hit between the last sample below the ray and the first one above it
    ra = np.broadcast_to(rowang[::-1, None], idx.shape)
    a0, a1 = ang[ip, cols], ang[ic, cols]
    t = np.where(ic > 0, np.clip((ra - a0) / np.maximum(a1 - a0, 1e-12), 0, 1), 1.0)
    depth = np.where(sky, np.inf, d[ip] + t * (d[ic] - d[ip]))
    hx = X[ip, cols] + t * (X[ic, cols] - X[ip, cols])
    hy = Y[ip, cols] + t * (Y[ic, cols] - Y[ip, cols])
    hs = dem.hillshade(a.az + a.sun_az, a.sun_alt)
    shade = dem.sample(hx, hy, grid=hs)
    shade = np.where(sky, 1.0, 0.2 + 0.8 * np.nan_to_num(shade))
    log(f"  view  {ncols} × {nrows} px, {a.fov:.0f}° × {math.degrees(top-bottom):.1f}° "
        f"(vertical stretch {a.vstretch}×)")
    return depth, shade, sky


def view_image(depth, shade, sky, dist_m):
    haze = 1 - np.exp(-np.where(sky, 0, depth) / (0.6 * dist_m))
    return np.where(sky, 1.0, shade * (1 - haze) + 0.88 * haze)


def wedge_grid(dem, x0, y0, a):
    """Resample the sector's bbox to the print grid; return heights (m), cell mask, cell coords."""
    dist = a.dist * 1000
    angs = np.radians(a.az + np.linspace(-a.fov / 2, a.fov / 2, 181))
    px = np.concatenate([[x0], x0 + dist * np.sin(angs)])
    py = np.concatenate([[y0], y0 + dist * np.cos(angs)])
    xmin, xmax, ymin, ymax = px.min(), px.max(), py.min(), py.max()
    mm_per_m = a.width / max(xmax - xmin, ymax - ymin)
    cell = a.res / mm_per_m
    xs = np.arange(xmin, xmax + cell, cell)
    ys = np.arange(ymax, ymin - cell, -cell)
    GX, GY = np.meshgrid(xs, ys)
    h = dem.sample(GX, GY)
    cx, cy = (GX[:-1, :-1] + GX[1:, 1:]) / 2, (GY[:-1, :-1] + GY[1:, 1:]) / 2
    dx, dy = cx - x0, cy - y0
    dc = np.hypot(dx, dy)
    diff = (np.degrees(np.arctan2(dx, dy)) - a.az + 180) % 360 - 180
    mask = (dc <= dist) & (np.abs(diff) <= a.fov / 2)
    mask |= dc <= cell * 1.5                                     # keep the apex solid
    return h, mask, cell, mm_per_m, (cx, cy, dc, diff)


def viewshed(dem, x0, y0, eye_z, a, cells, cell_m):
    cx, cy, dc, diff = cells
    dist = a.dist * 1000
    ncols = int(min(3000, math.radians(a.fov) * dist / cell_m * 2 + 2))
    nsteps = int(min(4000, dist / cell_m * 2))
    d, _, _, cm, _ = cast(dem, x0, y0, eye_z, a.az, a.fov, max(cell_m / 2, 1.0), dist, ncols, nsteps)
    col = np.clip(np.round((diff + a.fov / 2) / a.fov * (ncols - 1)).astype(int), 0, ncols - 1)
    k = np.searchsorted(d, dc)
    prev = np.where(k > 0, cm[np.maximum(k - 1, 0), col], -np.pi / 2)
    hc = dem.sample(cx, cy)
    drop = dc ** 2 / (2 * R_EARTH) * (1 - REFRACTION)
    angc = np.arctan2(hc - drop - eye_z, np.maximum(dc, 1e-3))
    return (angc >= prev - a.tolerance / np.maximum(dc, 1.0)) | (dc < cell_m * 2)


def cmd_view(a):
    modes = ["wedge", "visible", "relief", "litho"] if a.mode == "all" else [a.mode]
    dist_m = a.dist * 1000

    if a.dem is None:
        bbox = sector_bbox_lonlat(a.lat, a.lon, a.az, a.fov, dist_m)
        res = a.dem_res or max(2.0, round(dist_m / 2500))
        a.dem = f"{a.out}_dem.tif"
        log(f"fetching DEM for the view sector at {res:g} m …")
        fetch(*bbox, a.dem, res)
    dem = DEM(a.dem, max_px=a.max_dem_px)
    x0, y0 = lonlat_to(dem.wkt, [(a.lon, a.lat)])[0]
    ground = float(dem.sample(np.array([x0]), np.array([y0]))[0])
    if math.isnan(ground):
        sys.exit("Standpoint is outside the DEM")
    eye_z = ground + a.eye
    log(f"standpoint {a.lat:.5f}, {a.lon:.5f}: ground {ground:.0f} m, eye {eye_z:.0f} m, "
        f"looking {a.az:.0f}° ± {a.fov/2:.0f}°, {a.dist:g} km")

    if {"wedge", "visible"} & set(modes):
        width = a.width
        a.width = a.width_wedge or 150
        h, mask, cell, mm_per_m, cells = wedge_grid(dem, x0, y0, a)
        mask = no_diagonals(mask)
        hin = h[:-1, :-1][mask]
        hmin = np.nanmin(hin)
        hfill = np.where(np.isnan(h), hmin, h)
        top = a.base + np.clip(hfill - hmin, 0, None) * mm_per_m * a.exag
        sx = sy = cell * mm_per_m
        log(f"  wedge {mask.shape[1]*sx:.0f} × {mask.shape[0]*sy:.0f} mm, 1:{1000/mm_per_m:,.0f}, "
            f"{np.nanmin(hin):.0f}–{np.nanmax(hin):.0f} m")
        hs = shade_heights(top, scale=1 / sx)[:-1, :-1]
        if "wedge" in modes:
            write_stl(f"{a.out}_wedge.stl", grid_mesh(top, sx, sy, mask))
            write_png(f"{a.out}_wedge.png", np.where(mask, hs, 1.0))
        if "visible" in modes:
            vis = viewshed(dem, x0, y0, eye_z, a, cells, cell) & mask
            for _ in range(20):        # make both parts manifold
                vis = no_diagonals(vis, keep=mask)
                hid = no_diagonals(mask & ~vis, keep=mask)
                if not (hid & vis).any():
                    break
                vis &= ~hid
            log(f"  visible {100 * (vis & mask).sum() / mask.sum():.0f} % of the wedge")
            write_stl(f"{a.out}_visible.stl", grid_mesh(top, sx, sy, mask & vis))
            write_stl(f"{a.out}_hidden.stl", grid_mesh(top, sx, sy, mask & ~vis))
            # single-colour variant: hidden terrain sunk by --sink mm
            vgrid = np.zeros_like(top, bool)
            vgrid[:-1, :-1] |= vis; vgrid[1:, 1:] |= vis; vgrid[:-1, 1:] |= vis; vgrid[1:, :-1] |= vis
            sunk = np.where(vgrid, top, np.maximum(top - a.sink, a.base * 0.5))
            write_stl(f"{a.out}_visible_engraved.stl", grid_mesh(sunk, sx, sy, mask))
            col = np.where(vis[..., None], [1.0, 0.62, 0.25], [0.45, 0.55, 0.7]) * hs[..., None]
            write_png(f"{a.out}_visible.png", np.where(mask[..., None], col, 1.0))
        a.width = width

    if {"relief", "litho"} & set(modes):
        a.width = a.width or 200
        depth, shade, sky = panorama(dem, x0, y0, eye_z, a)
        img = view_image(depth, shade, sky, dist_m)
        write_png(f"{a.out}_view.png", img)

        if "relief" in modes:
            # depth layers: near = high, far = low, sky = flat; shading adds surface texture
            t = (np.log(np.where(sky, dist_m, depth)) - math.log(a.near)) / (math.log(dist_m) - math.log(a.near))
            t = np.clip(t, 0, 1)
            z = a.base + a.relief * (0.12 + 0.73 * (1 - t) + 0.15 * shade)
            z = np.where(sky, a.base, z)
            px = a.width / (z.shape[1] - 1)
            write_stl(f"{a.out}_relief.stl", grid_mesh(z, px, px))
            write_png(f"{a.out}_relief.png", shade_heights(z, az=300, alt=35, scale=1 / px))
            log(f"  relief plaque {a.width:.0f} × {z.shape[0]*px:.0f} mm, print face-up, stand upright")

        if "litho" in modes:
            # dark = thick; reads correctly from the textured side (flat side towards the light)
            thick = a.tmin + (1 - img) * (a.tmax - a.tmin)
            b = max(1, int(round(a.border / a.res_view)))
            thick = np.pad(thick, b, constant_values=a.tmax + 0.4)
            px = a.width / (thick.shape[1] - 1)
            write_stl(f"{a.out}_litho.stl", grid_mesh(thick, px, px))
            backlit = np.exp(-1.5 * thick)
            write_png(f"{a.out}_litho.png", (backlit - backlit.min()) / (backlit.max() - backlit.min()))
            log(f"  lithophane {thick.shape[1]*px:.0f} × {thick.shape[0]*px:.0f} mm, "
                f"{a.tmin}–{a.tmax} mm, print upright with 100 % infill, white PLA")


# ───────────────────────────── CLI ──────────────────────────────

def parse_bands(text):
    return [float(x) for x in text.split(",")] if text else None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"dgm3d {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    def stl_opts(sp, width=150):
        sp.add_argument("--width", type=float, default=width, help="long side in mm")
        sp.add_argument("--exag", type=float, default=1.5, help="vertical exaggeration")
        sp.add_argument("--base", type=float, default=3, help="base thickness in mm")
        sp.add_argument("--res", type=float, default=0.25, help="mesh spacing in mm")

    f = sub.add_parser("fetch", help="cut a bbox from the BEV DGM (→ .tif or .stl)")
    for c in ("west", "south", "east", "north"):
        f.add_argument(c, type=float)
    f.add_argument("out")
    f.add_argument("--dem-res", type=float, default=1, help="DEM resolution in m (default 1)")
    stl_opts(f)
    f.add_argument("--zmin", type=float, help="height (m) at the top of the base; default lowest point")
    f.add_argument("--crop", action="store_true", help="cut the model along the data outline (Austria's border)")
    f.add_argument("--bands", help="elevations in m for filament changes, e.g. 1000,2500")
    f.add_argument("--layer", type=float, default=0.1, help="layer height in mm (for --bands)")
    f.add_argument("--gamma", type=float, default=1.0,
                   help="<1 lifts lowlands relative to peaks (e.g. 0.6); peak height unchanged")

    s = sub.add_parser("stl", help="DEM GeoTIFF → STL")
    s.add_argument("dem"); s.add_argument("out")
    stl_opts(s)
    s.add_argument("--zmin", type=float)
    s.add_argument("--crop", action="store_true", help="cut the model along the data outline")
    s.add_argument("--bands", help="elevations in m for filament changes, e.g. 1000,2500")
    s.add_argument("--layer", type=float, default=0.1, help="layer height in mm (for --bands)")
    s.add_argument("--gamma", type=float, default=1.0,
                   help="<1 lifts lowlands relative to peaks (e.g. 0.6); peak height unchanged")

    v = sub.add_parser("view", help="prints of the view from a standpoint")
    v.add_argument("--lat", type=float, required=True)
    v.add_argument("--lon", type=float, required=True)
    v.add_argument("--az", type=float, required=True, help="view direction, degrees (0=N, 90=E)")
    v.add_argument("--fov", type=float, default=60, help="horizontal field of view, degrees")
    v.add_argument("--dist", type=float, default=30, help="how far to look, km")
    v.add_argument("--eye", type=float, default=1.7, help="eye height above ground, m")
    v.add_argument("--mode", choices=["wedge", "visible", "relief", "litho", "all"], default="all")
    v.add_argument("--out", default="view", help="output file prefix")
    v.add_argument("--dem", help="use this GeoTIFF instead of fetching from BEV")
    v.add_argument("--dem-res", type=float, help="fetch resolution in m (default: dist/2500, ≥2)")
    v.add_argument("--max-dem-px", type=int, default=6000)
    g = v.add_argument_group("wedge / visible")
    g.add_argument("--width-wedge", type=float, default=150, help="long side in mm")
    g.add_argument("--exag", type=float, default=1.5)
    g.add_argument("--base", type=float, default=3)
    g.add_argument("--res", type=float, default=0.3, help="mesh spacing in mm")
    g.add_argument("--sink", type=float, default=0.8, help="how deep hidden terrain is sunk, mm")
    g.add_argument("--tolerance", type=float, default=3, help="visibility tolerance in m")
    g = v.add_argument_group("relief / litho")
    g.add_argument("--width", type=float, default=200, help="plaque width in mm")
    g.add_argument("--res-view", type=float, default=0.2, help="pixel size in mm")
    g.add_argument("--near", type=float, default=400, help="relief: everything closer counts as foreground, m")
    g.add_argument("--vstretch", type=float, help="vertical stretch of the panorama (default auto, 1–4×)")
    g.add_argument("--bottom", type=float, help="lower edge of the view, degrees (default auto)")
    g.add_argument("--sky", type=float, default=0.12, help="sky margin above the skyline (fraction)")
    g.add_argument("--steps", type=int, default=5000, help="samples per ray")
    g.add_argument("--relief", type=float, default=4, help="relief depth, mm")
    g.add_argument("--tmin", type=float, default=0.8); g.add_argument("--tmax", type=float, default=3.2)
    g.add_argument("--border", type=float, default=3, help="lithophane frame, mm")
    g.add_argument("--sun-az", type=float, default=-70, help="light direction relative to view, degrees")
    g.add_argument("--sun-alt", type=float, default=30)

    a = p.parse_args(argv)
    if a.cmd == "fetch":
        if a.out.endswith(".stl"):
            tif = a.out[:-4] + ".tif"
            fetch(a.west, a.south, a.east, a.north, tif, a.dem_res)
            dem_to_stl(tif, a.out, a.width, a.exag, a.base, a.res, a.zmin, a.crop, parse_bands(a.bands), a.layer, a.gamma)
        else:
            fetch(a.west, a.south, a.east, a.north, a.out, a.dem_res)
    elif a.cmd == "stl":
        dem_to_stl(a.dem, a.out, a.width, a.exag, a.base, a.res, a.zmin, a.crop, parse_bands(a.bands), a.layer, a.gamma)
    else:
        cmd_view(a)


if __name__ == "__main__":
    main()
