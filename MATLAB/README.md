# ABANICCO — MATLAB

This is the **original/reference implementation** released with the 2023 ABANICCO paper. All MATLAB code, region definitions, and example images are contained in this folder.

## Requirements

- MATLAB R2021 or newer
- Image Processing Toolbox

The code uses functionality including `rgb2lab`, `rgb2ind`, `bwareaopen`, polar plots, and bubble charts.

## Files

| File | Purpose |
|---|---|
| `examplerun.m` | Recommended entry point; demonstrates the full workflow |
| `AnalyzeColorImage_final2_short.m` | Main analysis with a user-defined number of shades per color category |
| `AnalyzeColorImage_final2.m` | Full/original shade-analysis version |
| `hex2rgb.m` | Utility used by the MATLAB visualization code |
| `data/Regions_Angles0.mat` | ABANICCO color-region definitions |
| `data/Regions_Angles1.mat` | Shade naming / interval information |

## Run the example

Make this folder the MATLAB current folder and run:

```matlab
examplerun
```

The example loads `example_images/totoro.jpg` by default and defines three main parameters:

```matlab
wantfigure = 1;      % set to 0 to suppress optional figures
numberofshades = 15; % number of representative shades per category
morpho = 10;         % minimum morphology size in pixels
```

The recommended configurable call is:

```matlab
[Results, Regions_Angles1a, subTablesa] = ...
    AnalyzeColorImage_final2_short( ...
        an_im, Regions_Angles, Regions_Angles1, ...
        numberofshades, morpho, wantfigure);
```

## What `examplerun.m` demonstrates

The example script shows how to:

1. load the ABANICCO region definitions;
2. analyze an RGB image;
3. inspect the detected shades and their ABANICCO names;
4. retrieve a particular color-category segmentation;
5. create a custom cluster by combining categories (for example Blue + Ultramarine + Teal);
6. redefine angular/radial boundaries to create a custom color segmentation.

## MATLAB Online

The repository can also be opened directly in MATLAB Online:

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=lauranicolass/ABANICCO)

After opening the repository, enter this `MATLAB/` folder and run `examplerun.m`.

## Reproducibility

Use this implementation when you need the closest reproduction of the shade outputs from the original published workflow. MATLAB uses `rgb2ind` for shade quantization; the newer Python and Fiji ports use k-means, so representative shades may differ slightly even when the underlying category classification agrees.
