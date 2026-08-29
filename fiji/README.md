# ABANICCO — Fiji / ImageJ

Fiji/ImageJ implementation of ABANICCO as a graphical ImageJ1 plugin. This is the recommended option for users who want interactive color analysis and segmentation without writing MATLAB or Python code.

## Install the prebuilt plugin

1. Download `ABANICCO-1.0.0.jar` from this folder.
2. Copy it into:
   - `Fiji.app/plugins/` for Fiji, or
   - the `plugins/` directory of ImageJ.
3. Restart Fiji/ImageJ.
4. Open an RGB image.
5. Run **Plugins ▸ Analyze ▸ ABANICCO Color Analysis**.

The plugin bundles the ABANICCO region definitions and naming intervals, so no additional data files are required.

## Basic workflow

When the plugin starts, choose:

- the number of representative shades per color category;
- the morphology size used to remove small objects and fill small holes.

The main interface then provides:

1. **Original image** — the image being analyzed.
2. **Color-category segmentations** — the 12 ABANICCO categories and their detected content.
3. **Polar color representation** — the ABANICCO fan together with detected shades.
4. **Interactive inspection** — select category regions or shade markers to inspect percentages and shade information.
5. **Post-processing** — redefine boundaries, combine categories into clusters, or rerun with different shade/morphology parameters.
6. **Export** — save the polar representation, category segmentations, individual shade masks, shade measurements, modification history, and reproducibility code.

## Polar display

The interface supports:

- a **bubble view**, where marker size reflects the percentage of the image occupied by a shade;
- a **scatter view**, with fixed-size markers for easier inspection when many shades are close together.

Clicking a category reports its percentage of the image. Clicking a shade reports its angular/radial CIELAB position, lightness, percentage, and ABANICCO mixture-based color name.

## Build from source

The source is under `src/main/java/abanicco/` and uses Java 8-compatible bytecode.

With Maven:

```bash
mvn package
```

The resulting JAR will be created under `target/`.

The project depends on ImageJ1 (`net.imagej:ij`) with `provided` scope because Fiji/ImageJ supplies that dependency at runtime.

## Code layout

| Component | Purpose |
|---|---|
| `ABANICCO_.java` | ImageJ plugin entry point and graphical workflow |
| `AbaniccoEngine.java` | Core ABANICCO analysis logic |
| `Morphology.java` | Small-object removal and hole filling |
| `KMeans.java` | Shade quantization |
| `PolarPanel.java` | Interactive polar fan display |
| `src/main/resources/abanicco/` | Bundled ABANICCO regions, lookup table, boundaries, and naming intervals |
| `src/main/resources/plugins.config` | Registers **Plugins ▸ Analyze ▸ ABANICCO Color Analysis** |

## Implementation note

The Fiji port follows the same ABANICCO color-category definitions as the original MATLAB implementation. Shade quantization uses k-means rather than MATLAB's `rgb2ind`, so representative shades can differ slightly from MATLAB even when the category segmentation is comparable.
