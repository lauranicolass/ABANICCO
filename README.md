# ABANICCO

![ABANICCO logo](https://user-images.githubusercontent.com/49556605/201992997-359d29d3-8c88-4d22-8bba-1517db572d01.png)

**ABANICCO** (**AB AN**gular **I**llustrative **C**lassification of **CO**lor) is an automatic method for color classification, naming, visualization, and color-based image segmentation in the CIELAB *a–b* plane.

The method was introduced in:

> Nicolás-Sáenz L, Ledezma A, Pascau J, Muñoz-Barrutia A. **ABANICCO: A New Color Space for Multi-Label Pixel Classification and Color Analysis.** *Sensors*. 2023;23(6):3338. https://doi.org/10.3390/s23063338

ABANICCO is now available in **three implementations**. They follow the same underlying color-space model, but provide different ways of working with it.

## Choose an implementation

| Option | Implementation | Recommended for | Requirements |
|---|---|---|---|
| **1** | **[MATLAB — original/reference implementation](docs/MATLAB.md)** | Reproducing the original published workflow and MATLAB-based analysis | MATLAB R2021 or newer; Image Processing Toolbox |
| **2** | **[Python — package + Jupyter notebook](python/README.md)** | Reproducible scripting, notebooks, batch processing, and extending the method in Python | Python 3.9+ |
| **3** | **[Fiji / ImageJ — graphical plugin](fiji/README.md)** | Interactive image analysis without writing code | Fiji or ImageJ |

### Which one should I use?

- Use **MATLAB** if you want the implementation closest to the code released with the 2023 paper.
- Use **Python** if you want a fully scriptable open-source workflow, an interactive notebook, or integration with a larger Python image-analysis pipeline.
- Use **Fiji** if you want to run ABANICCO directly from a graphical image-analysis environment.

## What ABANICCO does

For an RGB image, ABANICCO:

1. converts the image to CIELAB;
2. classifies pixels into 12 predefined color categories using their position in the *a–b* plane;
3. optionally removes small objects and fills small holes;
4. identifies representative shades within each category;
5. describes each shade using its angle, radius, lightness, and a mixture-based color name;
6. visualizes the detected colors in the ABANICCO polar representation;
7. allows color categories to be inspected, combined, or redefined for color-based segmentation.

## Quick start

### Option 1 — MATLAB

The original/reference MATLAB implementation is contained in the `MATLAB/` folder. From the repository root:

```matlab
cd MATLAB
examplerun
```

It demonstrates the main analysis, shade descriptions, category retrieval, custom clusters, and boundary redefinition. See **[docs/MATLAB.md](docs/MATLAB.md)** for details.

If you do not have a local MATLAB installation, the repository can also be opened in MATLAB Online:

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=lauranicolass/ABANICCO)

In MATLAB Online, open the `MATLAB/` folder and run `examplerun.m`.

### Option 2 — Python

```bash
cd python
python -m pip install -e ".[notebook]"
jupyter lab notebooks/ABANICCO_interactive.ipynb
```

The Python folder also exposes ABANICCO as a normal Python package for scripts and batch processing. See **[python/README.md](python/README.md)**.

### Option 3 — Fiji / ImageJ

1. Download `fiji/ABANICCO-1.0.0.jar`.
2. Copy it into `Fiji.app/plugins/` (or the `plugins/` folder of ImageJ).
3. Restart Fiji/ImageJ.
4. Open an RGB image.
5. Run **Plugins ▸ Analyze ▸ ABANICCO Color Analysis**.

See **[fiji/README.md](fiji/README.md)** for usage and build instructions.

## Implementation note

The **color-category model is shared across the implementations**, but shade quantization is not numerically identical in every port. The original MATLAB implementation uses `rgb2ind` minimum-variance quantization, whereas the Python and Fiji implementations use k-means. Therefore, category-level segmentations should be directly comparable, while the exact number and RGB values of representative shades can differ slightly between implementations.

For analyses where exact reproduction of previously published MATLAB shade outputs is required, use the MATLAB implementation.

## Repository layout

```text
ABANICCO/
├── README.md
├── MATLAB/                             # original/reference MATLAB implementation
│   ├── examplerun.m                    # recommended MATLAB entry point
│   ├── AnalyzeColorImage_final2.m
│   ├── AnalyzeColorImage_final2_short.m
│   ├── hex2rgb.m
│   ├── data/                           # MATLAB region definitions
│   └── example_images/                 # example images used by MATLAB and the Python notebook
├── python/                             # Python package + notebook
├── fiji/                               # Fiji/ImageJ plugin + source
└── docs/
    └── MATLAB.md                       # MATLAB guide
```

## Multiplex images

For multiplex images, it is usually more informative to run ABANICCO on the individual channels and compare those results with the merged image. Brightness/contrast adjustment can make visualization clearer, but it can also amplify noise; ABANICCO itself does not require contrast enhancement as a preprocessing step.

## Citation

If ABANICCO is useful in your work, please cite the original paper:

```bibtex
@article{nicolas-saenz2023abanicco,
  title   = {ABANICCO: A New Color Space for Multi-Label Pixel Classification and Color Analysis},
  author  = {Nicolás-Sáenz, Laura and Ledezma, Agapito and Pascau, Javier and Muñoz-Barrutia, Arrate},
  journal = {Sensors},
  volume  = {23},
  number  = {6},
  pages   = {3338},
  year    = {2023},
  doi     = {10.3390/s23063338}
}
```

## Feedback and issues

The Python and Fiji versions are newer ports of the original MATLAB implementation. If you find an inconsistency between implementations, please open a GitHub issue and include the input image, parameters used, and the implementation/version.
