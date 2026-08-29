# ABANICCO — Python

Python implementation of ABANICCO with both:

- a reusable Python package (`abanicco`), and
- an interactive Jupyter notebook (`notebooks/ABANICCO_interactive.ipynb`).

This is the recommended option for scripting, batch analysis, integration with other Python image-processing tools, and reproducible notebook workflows.

## Requirements

- Python 3.9 or newer
- NumPy
- scikit-image
- scikit-learn
- Matplotlib
- Plotly
- pandas
- Pillow
- Jupyter + ipywidgets for the interactive notebook

## Install

From the repository root:

```bash
cd python
python -m pip install -e ".[notebook]"
```

For package-only use without the notebook extras:

```bash
python -m pip install -e .
```

## Interactive notebook

Launch:

```bash
jupyter lab notebooks/ABANICCO_interactive.ipynb
```

Running all cells provides four sections:

1. **Input and analysis** — choose an example image or enter a path, then set the number of shades and morphology size.
2. **Visualization** — inspect category segmentations, the ABANICCO polar representation, detected shades, and the shade table.
3. **Post-processing** — combine color categories, redefine AB boundaries, or automatically cluster detected shades in the AB plane.
4. **Export** — save polar plots, segmentations, individual shade masks, a CSV table, the modification log, and reproducible Python code.

## Minimal Python API example

```python
import numpy as np
from PIL import Image

from abanicco import analyze, show_segmentations, polar_description

im = np.asarray(Image.open("../MATLAB/example_images/fish2.jpg").convert("RGB"))

regions, shades = analyze(
    im,
    numshades=15,
    morpho=10,
)

show_segmentations(im, regions)
polar_description(regions, shades)
```

`analyze()` returns:

- `regions`: the 12 ABANICCO color-category objects with masks, segmentations, percentages, and detected shades;
- `shades`: one record per retained representative shade.

## Load the underlying ABANICCO definitions

Most users do **not** need to call `load_regions()` directly because `analyze()` loads the bundled defaults automatically. For advanced use:

```python
from abanicco import load_regions

regions, lines, lines2, intervals, intervals_rad = load_regions()
```

## Export shade information

```python
import pandas as pd
from abanicco import shades_to_records

pd.DataFrame(shades_to_records(shades)).to_csv("shades.csv", index=False)
```

## Combine color categories

```python
from abanicco import make_cluster

mask, segmented = make_cluster(
    im,
    regions,
    ["Blue", "Ultramarine", "Teal"],
    cluster_name="Cold colours",
)
```

## Define a custom AB region

```python
from abanicco import custom_region

mask, segmented = custom_region(
    im,
    theta1_deg=90,
    theta2_deg=180,
    rho_min=30,
)
```

## Automatically cluster the detected shades

```python
from abanicco import cluster_shades, cluster_polar, cluster_summary

labels = cluster_shades(shades, method="K-Means", n_clusters=4)
cluster_polar(regions, shades, labels).show()
summary = cluster_summary(shades, labels)
```

Available methods include Affinity Propagation, Agglomerative Clustering, BIRCH, DBSCAN, K-Means, Mini-Batch K-Means, Mean Shift, OPTICS, Spectral Clustering, and Gaussian Mixture Models.

## Implementation note

The Python port uses the same ABANICCO color-region definitions and shade-naming intervals as the MATLAB workflow. Representative shade extraction uses k-means rather than MATLAB's `rgb2ind`, so exact shade centroids can differ slightly from the MATLAB output.

The bundled `abanicco_data.json` makes the Python package self-contained; the original MATLAB `.mat` files in `MATLAB/data/` are not required when using the Python package. The example images used by the notebook are read from `MATLAB/example_images/`.
