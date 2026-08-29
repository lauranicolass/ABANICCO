"""
ABANICCO: AB ANgular Illustrative Classification of COlor
Python port of AnalyzeColorImage_final2_short.m

Reference: https://www.mdpi.com/1424-8220/23/6/3338
"""

import json
import os
import numpy as np
from skimage import color as skcolor
from skimage import morphology as skmorph
from sklearn.cluster import KMeans

_DATA_PATH = os.path.join(os.path.dirname(__file__), "abanicco_data.json")


# ----------------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------------

class Region:
    """One color category: an angular (and radial) region of the AB plane."""

    def __init__(self, d):
        self.name = d["name"]
        self.code = d["code"]                 # hex color for display
        self.theta1 = d["theta1"]             # radians (bounds of the wedge)
        self.theta2 = d["theta2"]
        self.rho1 = d["rho1"]                 # display-ring radii (100..145)
        self.rho2 = d["rho2"]
        self.base = d["base"]                 # central angle of the wedge
        # Membership lookup: points are rows of [a+127, b+127] (integers)
        pts = np.asarray(d["points"], dtype=int)
        self.lut = np.zeros((256, 256), dtype=bool)
        ok = (pts[:, 0] >= 0) & (pts[:, 0] < 256) & (pts[:, 1] >= 0) & (pts[:, 1] < 256)
        self.lut[pts[ok, 0], pts[ok, 1]] = True
        # results (filled by analyze)
        self.mask = None
        self.segmented = None
        self.mask_percentage = 0.0
        self.shades = []


def load_regions(path=_DATA_PATH):
    """Load the default ABANICCO region definitions and mixing intervals."""
    with open(path) as f:
        data = json.load(f)
    regions = [Region(d) for d in data["regions"]]
    lines = np.asarray(data["lines"])     # (10,2) boundary line angles
    lines2 = np.asarray(data["lines2"])   # (10,2) radial extents (x98)
    intervals = data["intervals"]         # angular mixing intervals
    intervals_rad = data["intervals_rad"]  # radial (brown/pure) intervals
    return regions, lines, lines2, intervals, intervals_rad


# ----------------------------------------------------------------------------
# Core analysis
# ----------------------------------------------------------------------------

def _classify(lab, regions):
    """Assign each pixel to every region whose AB lookup table contains it."""
    a = lab[:, :, 1]
    b = lab[:, :, 2]
    ai = np.clip(np.floor(a + 127).astype(int), 0, 255)
    bi = np.clip(np.floor(b + 127).astype(int), 0, 255)
    masks = {}
    for r in regions:
        masks[r.name] = r.lut[ai, bi]
    return masks


def _remove_small(mask, min_size):
    """MATLAB bwareaopen: remove components with fewer than min_size pixels."""
    try:  # skimage >= 0.26
        return skmorph.remove_small_objects(mask, max_size=min_size - 1,
                                            connectivity=2)
    except TypeError:  # older skimage
        return skmorph.remove_small_objects(mask, min_size=min_size,
                                            connectivity=2)


def _clean_mask(mask, morpho):
    """bwareaopen + hole filling, as in the MATLAB code."""
    if morpho <= 0:
        return mask
    m = _remove_small(mask, morpho)
    inv = _remove_small(~m, 2 * morpho)
    return ~inv


def _quantize_shades(rgb, mask, numshades, seed=0):
    """Quantize masked pixels into `numshades` clusters (MATLAB rgb2ind analog).

    Returns (label_image, palette) where label 0 = background and labels
    1..k index rows of `palette` (float RGB in 0..1).
    """
    labels = np.zeros(mask.shape, dtype=int)
    pix = rgb[mask].astype(float)
    if rgb.dtype == np.uint8:
        pix /= 255.0
    n = pix.shape[0]
    if n == 0:
        return labels, np.zeros((0, 3))
    k = int(min(numshades, max(1, len(np.unique((pix * 255).astype(np.uint8), axis=0)))))
    k = min(k, n)
    km = KMeans(n_clusters=k, n_init=4, random_state=seed).fit(pix)
    labels[mask] = km.labels_ + 1
    return labels, km.cluster_centers_


def wrap_to_360(deg):
    d = np.mod(deg, 360.0)
    return d


def analyze(an_im, regions=None, intervals=None, intervals_rad=None,
            numshades=15, morpho=10, min_shade_percent=5.0, seed=0):
    """Run the full ABANICCO analysis.

    Parameters
    ----------
    an_im : (H,W,3) uint8 RGB image
    regions : list[Region] or None (default regions loaded from data file)
    numshades : shades detected per category
    morpho : min object size / half min hole size, in pixels
    min_shade_percent : keep shades covering more than this % of the category

    Returns
    -------
    regions : list of Region with .mask, .segmented, .mask_percentage, .shades
    shade_table : list of dicts (one row per detected shade), with keys
        theta, rho, L, deg, percent, rgb, region, code, base, shade_name
    """
    if regions is None:
        regions, _, _, dflt_int, dflt_rad = load_regions()
        intervals = intervals or dflt_int
        intervals_rad = intervals_rad or dflt_rad

    an_im = np.asarray(an_im)
    if an_im.ndim == 2:
        an_im = np.stack([an_im] * 3, axis=-1)
    if an_im.shape[2] == 4:
        an_im = an_im[:, :, :3]

    lab = skcolor.rgb2lab(an_im)
    total_pix = an_im.shape[0] * an_im.shape[1]

    masks = _classify(lab, regions)

    shade_rows = []
    for reg in regions:
        m = _clean_mask(masks[reg.name], morpho)
        reg.mask = m
        reg.segmented = an_im * m[:, :, None].astype(an_im.dtype)
        reg.mask_percentage = 100.0 * m.sum() / total_pix
        reg.shades = []

        if m.sum() == 0:
            continue

        labels, palette = _quantize_shades(an_im, m, numshades, seed=seed)
        region_pix = (labels > 0).sum()
        for j in range(palette.shape[0]):
            mj = labels == (j + 1)
            pct_of_region = 100.0 * mj.sum() / max(region_pix, 1)
            if pct_of_region <= min_shade_percent:
                continue
            rgb_shade = np.clip(palette[j], 0, 1)
            lab_shade = skcolor.rgb2lab(rgb_shade.reshape(1, 1, 3)).ravel()
            L, A, B = lab_shade
            theta = np.arctan2(A, B)  # cart2pol(b, a) -> atan2(a, b)
            rho = np.hypot(A, B)
            row = dict(
                theta=float(theta),
                rho=float(rho),
                L=float(L),
                deg=float(wrap_to_360(np.degrees(theta))),
                percent=float(100.0 * mj.sum() / total_pix),
                rgb=rgb_shade.tolist(),
                region=reg.name,
                code=reg.code,
                base=reg.base,
                mask=mj,
            )
            reg.shades.append(row)
            shade_rows.append(row)

    if intervals is not None:
        _name_shades(shade_rows, intervals, intervals_rad)

    return regions, shade_rows


# ----------------------------------------------------------------------------
# Shade naming (mixing percentages), port of the "Show mixing" section
# ----------------------------------------------------------------------------

def _name_shades(rows, intervals, intervals_rad):
    # Find the pair of radial intervals sharing identical bounds (Brown/Pure)
    bounds = [(iv["x1"], iv["x2"]) for iv in intervals_rad]
    dup_idx = []
    for i, b in enumerate(bounds):
        if bounds.count(b) > 1:
            dup_idx.append(i)
    Z = dup_idx[:2] if len(dup_idx) >= 2 else None

    for row in rows:
        th = wrap_to_360(np.degrees(row["theta"]))
        R = row["rho"]
        belongs = [k for k, iv in enumerate(intervals)
                   if iv["x1"] < th < iv["x2"]]
        if not belongs:
            row["shade_name"] = row["region"]
            row["angle"] = float(th)
            continue

        idxs = sorted(belongs, reverse=True)
        if intervals[idxs[0]]["Color"] != row["region"]:
            idxs = sorted(belongs)
        iv1 = intervals[idxs[0]]
        val1 = th * iv1["p1"] + iv1["p2"]
        totcolors = 1
        iv2 = None
        val2 = None
        if len(idxs) > 1:
            iv2 = intervals[idxs[1]]
            val2 = th * iv2["p1"] + iv2["p2"]
            totcolors = 2

        # Radial brown mixing (only for 0 < theta < 90 degrees)
        iswithbrown = False
        val_pure = val_brown = None
        if 0 < th < 90 and Z is not None:
            z1, z2 = intervals_rad[Z[0]], intervals_rad[Z[1]]
            if z1["x1"] < R < z1["x2"]:
                iswithbrown = True
                v1 = R * z1["p1"] + z1["p2"]
                v2 = R * z2["p1"] + z2["p2"]
                if z1["Color"] == "Pure":
                    val_pure, val_brown = v1, v2
                else:
                    val_pure, val_brown = v2, v1

        if not iswithbrown:
            if totcolors == 1:
                name = f"{val1 * 100:.1f}% of {iv1['Color']}"
            else:
                name = (f"{val1 * 100:.1f}% of {iv1['Color']} and "
                        f"{val2 * 100:.1f}% of {iv2['Color']}")
        else:
            if totcolors == 1:
                name = (f"{val_brown * 100:.1f}% of Brown, and "
                        f"{val_pure * val1 * 100:.1f}% of {iv1['Color']}")
            else:
                name = (f"{val_brown * 100:.1f}% of Brown, "
                        f"{val_pure * val1 * 100:.1f}% of {iv1['Color']}, and "
                        f"{val_pure * val2 * 100:.1f}% of {iv2['Color']}")
            row["val_pure"] = float(val_pure)
            row["val_brown"] = float(val_brown)

        row["val1"] = float(val1)
        row["val1_name"] = iv1["Color"]
        if totcolors == 2:
            row["val2"] = float(val2)
            row["val2_name"] = iv2["Color"]
        row["totcolors"] = totcolors
        row["shade_name"] = name
        row["angle"] = float(th)


# ----------------------------------------------------------------------------
# Post-processing utilities (cluster creation and custom boundaries)
# ----------------------------------------------------------------------------

def make_cluster(an_im, regions, names, cluster_name="Cluster"):
    """Merge the masks of several categories into one cluster segmentation."""
    an_im = np.asarray(an_im)
    mask = np.zeros(an_im.shape[:2], dtype=bool)
    for nm in names:
        for r in regions:
            if r.name == nm and r.mask is not None:
                mask |= r.mask
    seg = an_im * mask[:, :, None].astype(an_im.dtype)
    return mask, seg


def custom_region(an_im, theta1_deg, theta2_deg, rho_min=0.0, rho_max=1e9):
    """Segment pixels whose AB polar coordinates fall in a custom region.

    Angles in degrees on the wrapped 0-360 convention (theta measured with
    cart2pol(b, a), i.e. atan2(a, b)).
    """
    an_im = np.asarray(an_im)
    lab = skcolor.rgb2lab(an_im[:, :, :3])
    a = lab[:, :, 1]
    b = lab[:, :, 2]
    theta = wrap_to_360(np.degrees(np.arctan2(a, b)))
    rho = np.hypot(a, b)
    if theta1_deg <= theta2_deg:
        m = (theta > theta1_deg) & (theta < theta2_deg)
    else:  # wrap-around region
        m = (theta > theta1_deg) | (theta < theta2_deg)
    m &= (rho > rho_min) & (rho < rho_max)
    seg = an_im * m[:, :, None].astype(an_im.dtype)
    return m, seg


# ----------------------------------------------------------------------------
# Export helpers
# ----------------------------------------------------------------------------

def shades_to_records(rows):
    """Shade rows -> plain records (no mask arrays) for CSV/JSON export."""
    recs = []
    for r in rows:
        recs.append({
            "region": r["region"],
            "shade_name": r.get("shade_name", ""),
            "theta_rad": r["theta"],
            "angle_deg": r.get("angle", wrap_to_360(np.degrees(r["theta"]))),
            "rho": r["rho"],
            "L": r["L"],
            "percent_of_image": r["percent"],
            "R": r["rgb"][0], "G": r["rgb"][1], "B": r["rgb"][2],
            "hex": "#%02x%02x%02x" % tuple(int(round(c * 255)) for c in r["rgb"]),
        })
    return recs


# ----------------------------------------------------------------------------
# Automatic clustering of shades on the AB plane
# ----------------------------------------------------------------------------

CLUSTER_METHODS = [
    "Affinity Propagation",
    "Agglomerative Clustering",
    "BIRCH",
    "DBSCAN",
    "K-Means",
    "Mini-Batch K-Means",
    "Mean Shift",
    "OPTICS",
    "Spectral Clustering",
    "Mixture of Gaussians",
]

# Methods that determine the number of clusters themselves
AUTO_K_METHODS = {"Affinity Propagation", "DBSCAN", "Mean Shift", "OPTICS"}


def cluster_shades(rows, method="K-Means", n_clusters=4, **kwargs):
    """Automatically cluster the detected shades on the AB plane.

    The polar coordinates (theta, rho) of each shade are converted back to
    Cartesian AB coordinates (b = rho*cos(theta), a = rho*sin(theta)) before
    clustering, so distances are Euclidean in the AB plane and there is no
    wrap-around problem at 0/360 degrees.

    Parameters
    ----------
    rows : list of shade dicts (as returned by ``analyze``)
    method : one of ``CLUSTER_METHODS``
    n_clusters : requested number of clusters. Ignored by the methods in
        ``AUTO_K_METHODS`` (Affinity Propagation, DBSCAN, Mean Shift, OPTICS),
        which estimate it from the data.
    **kwargs : forwarded to the underlying scikit-learn estimator
        (e.g. ``eps=`` for DBSCAN, ``min_samples=`` for OPTICS).

    Returns
    -------
    labels : (n,) int array; -1 marks noise points (DBSCAN / OPTICS).
    """
    from sklearn import cluster as skc
    from sklearn.mixture import GaussianMixture

    theta = np.array([r["theta"] for r in rows])
    rho = np.array([r["rho"] for r in rows])
    X = np.column_stack([rho * np.cos(theta), rho * np.sin(theta)])  # (b, a)

    n_clusters = int(max(1, min(n_clusters, len(rows))))
    m = method.strip().lower()

    if m == "affinity propagation":
        model = skc.AffinityPropagation(random_state=0, **kwargs)
    elif m == "agglomerative clustering":
        model = skc.AgglomerativeClustering(n_clusters=n_clusters, **kwargs)
    elif m == "birch":
        model = skc.Birch(n_clusters=n_clusters, **kwargs)
    elif m == "dbscan":
        kwargs.setdefault("eps", 12.0)
        kwargs.setdefault("min_samples", 3)
        model = skc.DBSCAN(**kwargs)
    elif m == "k-means":
        model = skc.KMeans(n_clusters=n_clusters, n_init=10, random_state=0,
                           **kwargs)
    elif m == "mini-batch k-means":
        model = skc.MiniBatchKMeans(n_clusters=n_clusters, n_init=10,
                                    random_state=0, **kwargs)
    elif m == "mean shift":
        model = skc.MeanShift(**kwargs)
    elif m == "optics":
        kwargs.setdefault("min_samples", 3)
        model = skc.OPTICS(**kwargs)
    elif m == "spectral clustering":
        kwargs.setdefault("assign_labels", "kmeans")
        model = skc.SpectralClustering(n_clusters=n_clusters, random_state=0,
                                       **kwargs)
    elif m in ("mixture of gaussians", "gaussian mixture"):
        gm = GaussianMixture(n_components=n_clusters, random_state=0, **kwargs)
        return gm.fit_predict(X)
    else:
        raise ValueError(f"Unknown clustering method: {method!r}. "
                         f"Choose one of {CLUSTER_METHODS}")
    return model.fit_predict(X)


def cluster_summary(rows, labels):
    """Per-cluster summary records: size, % of image, dominant categories."""
    out = []
    labels = np.asarray(labels)
    for lb in sorted(set(labels.tolist())):
        idx = np.where(labels == lb)[0]
        members = [rows[i] for i in idx]
        per = sum(r["percent"] for r in members)
        cats = {}
        for r in members:
            cats[r["region"]] = cats.get(r["region"], 0) + r["percent"]
        cats = sorted(cats.items(), key=lambda kv: -kv[1])
        out.append({
            "cluster": "noise" if lb == -1 else int(lb),
            "n_shades": len(members),
            "percent_of_image": round(per, 2),
            "categories": ", ".join(f"{k} ({v:.1f}%)" for k, v in cats[:3]),
        })
    return out
