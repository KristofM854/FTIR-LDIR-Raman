# Multi-Instrument Microplastic Particle Matching Pipeline

Spatial alignment and material cross-validation of microplastic particles detected by **FTIR**, **Raman**, and **LDIR** microspectroscopy on the same filter sample.

## What it does

When the same filter is analysed by multiple instruments, each produces its own particle list with coordinates, sizes, and material identifications. This pipeline:

1. **Ingests** data from all three instruments (CSV/Excel)
2. **Aligns** the coordinate systems using a tiered approach (landmark -> RANSAC -> ICP)
3. **Matches** particles spatially across instruments
4. **Compares** material identifications to assess inter-instrument agreement
5. **Exports** matched pairs, agreement statistics, and diagnostic plots

The LDIR instrument does not export spatial coordinates, so its comparison is limited to bulk material composition.

## Quick start

### Prerequisites

R (>= 4.0) with the following packages:

```r
install.packages(c("readxl", "ggplot2", "RANN", "png"))
```

### Running the pipeline

**Option A — Hardcoded paths** (non-interactive):

```r
ftir_file  <- "Comparstic Spotlight F2Ba Au 240926.csv"
raman_file <- "Comparstic Raman F2Ba Au IAEA 240930.csv"
ldir_file  <- "Comparstic LDIR F2Ba Au MP2 240925.xlsx"          # optional
ftir_image <- "Average Abs.( Comparstic Spotlight F2Ba Au 240926 ).png"  # optional
source("main.R")
```

**Option B — Interactive file picker**:

```r
input_mode <- "interactive"
source("main.R")
```

Results are written to a timestamped subfolder under `output/`.

## Input data formats

| Instrument | Format | Key columns |
|------------|--------|-------------|
| **FTIR** (PerkinElmer Spotlight) | CSV (comma-separated) | `Coord. [um]` (bracketed X;Y), `Group` (material), `Max AAU score` |
| **Raman** (Horiba / WITec) | CSV (semicolon-separated) | `Visual Center Point X/Y [um]`, `Material`, `HQI`, `Feret Max [um]` |
| **LDIR** (Agilent 8700) | XLSX | `Identification` (material), `Quality`, `Width/Height [um]` — no X/Y coordinates |
| **FTIR image** | PNG | Chemical absorption map exported from the FTIR instrument |

The ingest module auto-detects CSV delimiters and handles encoding differences (Latin-1, UTF-8).

## Pipeline architecture

```
main.R                          # Orchestrator
R/
  00_config.R                   # All tuneable parameters
  00b_file_input.R              # Interactive file picker
  01_ingest.R                   # FTIR + Raman CSV/Excel parsing
  01b_ingest_image.R            # Image-based particle extraction
  01c_ingest_ldir.R             # LDIR Excel parsing
  02_prefilter.R                # Quality/size pre-filtering
  03_normalize.R                # Coordinate centering
  03b_landmark_align.R          # Tier 1: large-particle landmark alignment
  04_ransac.R                   # Tier 2: RANSAC with rotation grid search
  05_transform.R                # Affine transform helpers
  06_icp_refine.R               # Iterative Closest Point refinement
  07_match.R                    # Nearest-neighbor spatial matching
  08_agreement.R                # Material agreement scoring
  08b_material_map.R            # Material name equivalence mapping
  09_diagnostics.R              # Overlay plots, histograms, confusion matrix
  10_export.R                   # CSV/PDF/TXT export
  utils.R                       # Logging, geometry, transform utilities
```

### Alignment strategy

The coordinate systems from FTIR and Raman differ in origin, rotation (often ~180 degrees), and sometimes scale. The pipeline finds the spatial transform in three tiers:

1. **Tier 1 — Landmarks**: Particles >= 100 um and fibers (aspect ratio >= 3) are matched between instruments. If enough landmarks agree (>= 50% inliers, residual < 50 um), RANSAC is skipped.

2. **Tier 2 — Material-anchored RANSAC**: PET and PP particles (identified by both instruments with HQI >= 70 on the Raman side) serve as anchor points. A coarse rotation grid search (1-degree steps, including mirror check) finds the best angle, then RANSAC refines the similarity transform.

3. **ICP refinement**: Iterative Closest Point polishes the transform using all particles >= 20 um, with reciprocal nearest-neighbour filtering and 10% trimming.

### Raman particle filtering by pipeline step

| Step | Size filter | Material filter | HQI filter |
|------|-------------|-----------------|------------|
| Spatial transform (landmarks, ICP) | >= 20 um | none | none |
| Material-based alignment (RANSAC) | >= 20 um | PET / PP only | >= 70 |
| Spatial matching | all | none | none |
| Agreement scoring | all | none | >= 70 |

### Image-based particle extraction

When an FTIR chemical image (PNG) is provided, particles are extracted using:

- **Grayscale**: `mean(R, G, B)` — colour-agnostic, no channel bias
- **Adaptive threshold**: Each pixel is compared to its local mean (box filter via integral image). Foreground = pixels deviating by more than the offset from local mean.
- **Size filter**: Blobs smaller than 20 um (physical) are removed
- **Sanity check**: If extracted count exceeds 1.5x the tabular count, only the largest blobs are kept

## Output files

Each run creates a timestamped folder (e.g., `output/2026-02-16_6/`) containing:

| File | Description |
|------|-------------|
| `matched_particles.csv` | All matched FTIR-Raman pairs with coordinates, materials, distances |
| `unmatched_ftir.csv` | FTIR particles with no Raman match |
| `unmatched_raman.csv` | Raman particles with no FTIR match |
| `agreement_summary.csv` | Per-material agreement rates |
| `agreement_pairwise.csv` | Per-pair material comparison |
| `transform_params.txt` | Final spatial transform (scale, rotation, translation) |
| `match_statistics.txt` | Summary statistics |
| `ldir_material_distribution.csv` | LDIR material counts (if LDIR provided) |
| `plots/overlay.png` | Spatial overlay of aligned point clouds |
| `plots/confusion.png` | Material confusion matrix |
| `plots/distance_hist.png` | Match distance distribution |
| `plots/all_diagnostics.pdf` | Combined PDF of all plots |

## Configuration

All parameters are in `R/00_config.R`. Key settings:

```r
# Alignment anchors
align_ftir_materials      = c("PET", "Polypro")
align_raman_materials     = c("Polyethylene terephtalate", "Polypropylene")
align_raman_min_size_um   = 20     # Raman particles below this are excluded from alignment

# Quality thresholds
raman_hqi_threshold       = 70     # HQI cutoff for material agreement scoring

# Landmark alignment
landmark_min_size_um      = 100    # particles >= this are landmarks
landmark_fiber_aspect_ratio = 3.0  # fibers detected by aspect ratio

# Matching
match_dist_threshold_um   = 100    # max distance for a spatial match (um)

# RANSAC
ransac_inlier_dist_um     = 200    # inlier threshold (um)
ransac_allow_mirror        = TRUE  # search over reflections too
```

## Limitations

- **FTIR feret size**: The PerkinElmer Spotlight CSV does not export a dedicated Feret Max column. The pipeline uses `Major dim` as a proxy, which may overestimate actual Feret dimensions for non-convex particles.

- **Raman minimum size**: The Raman instrument detects particles down to ~1 um, but FTIR resolution is typically >= 20 um. Raman particles below 20 um have no FTIR counterpart and inflate the "unmatched Raman" count.

- **Material name mapping**: FTIR and Raman use different spectral libraries with different naming conventions (e.g., FTIR: "Polypro", Raman: "Polypropylene (PP)", LDIR: "Polypropylene"). The pipeline normalises common abbreviations but may miss unusual names. Extend the mapping in `08b_material_map.R` or `00_config.R` for your datasets.

- **LDIR spatial matching**: The Agilent 8700 LDIR does not export particle coordinates, so LDIR comparison is limited to bulk material composition. Spatial matching is only possible between FTIR and Raman.

- **Image extraction quality**: The adaptive threshold works well on high-contrast images but may over- or under-segment low-contrast or noisy images. The sanity check helps, but better results come from providing B&W pre-processed images from the instrument software.

- **Single-filter assumption**: The pipeline assumes all instruments analysed the same physical filter. If scan areas differ significantly, the alignment may fail or produce a poor transform.

## License

Internal research tool. Not yet licensed for public distribution.
