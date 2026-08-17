# Multi-Instrument Microplastic Particle Matching Pipeline

Spatial alignment and material cross-validation of microplastic particles detected by **FTIR**, **Raman**, and **LDIR** microspectroscopy on the same filter sample.

## What it does

When the same filter is analysed by multiple instruments, each produces its own particle list with coordinates, sizes, and material identifications. This pipeline:

1. **Ingests** tabular data and companion images from all three instruments
2. **Extracts spatial coordinates** for LDIR via a Python-based image processing backend (scipy/numpy sigma-clipping + auto-tuned threshold)
3. **Aligns** coordinate systems using a tiered approach (Landmark → RANSAC → ICP)
4. **Matches** particles spatially across all three instrument pairs (FTIR↔Raman, LDIR↔Raman, LDIR↔FTIR)
5. **Identifies triplets** — particles detected by all three instruments simultaneously, with material-agreement quality scoring
6. **Compares** material identifications to assess inter-instrument agreement at Exact / Family / Disagree tiers
7. **Exports** matched pairs, agreement statistics, diagnostic plots, and a Shiny viewer

## Quick start

### Prerequisites

**R (≥ 4.1)** — the code uses the native `|>` pipe. Install the full package set:
```r
install.packages(c(
  "clue", "dplyr", "ggplot2", "ggrepel", "jpeg", "jsonlite", "magick",
  "png", "RANN", "readxl", "rintrojs", "shiny",
  "reticulate"   # optional — enables the preferred Python LDIR backend
))
```

**Python 3** with these packages (optional — for the preferred LDIR image backend):
```bash
pip3 install numpy scipy Pillow matplotlib
```

The Python backend is invoked automatically via `reticulate`. If Python or its
packages are unavailable, the pipeline falls back to an R-based saturation
segmentation approach. To pin a specific interpreter, set the
`RETICULATE_PYTHON` environment variable (e.g. in `~/.Renviron`) before
launching R; when unset, `reticulate` auto-discovers Python on `PATH`.

#### Reproducible environment (recommended for the deposit)

The exact dependency versions are captured with [`renv`](https://rstudio.github.io/renv/).
A `DESCRIPTION` file declares the dependency surface; run once on your machine to
freeze versions into a committed `renv.lock`:
```r
install.packages("renv")
renv::init(bare = TRUE)              # generates renv/activate.R
renv::snapshot(type = "explicit")   # writes renv.lock from DESCRIPTION
```
Anyone re-running the analysis then restores the identical environment with:
```r
renv::restore()
```

### Running the pipeline

**Option A — Hardcoded paths** (non-interactive):

```r
ftir_file  <- "Comparstic Spotlight F2Ba Au 240926.csv"
raman_file <- "Comparstic Raman F2Ba Au IAEA 240930.csv"
ldir_file  <- "Comparstic LDIR F2Ba Au MP2 240925.xlsx"          # optional
ftir_image <- "Average Abs.( Comparstic Spotlight F2Ba Au 240926 ).png"  # optional
ldir_image <- "LDIR_particle_map.png"                             # optional
source("main.R")
```

**Option B — Interactive file picker**:

```r
input_mode <- "interactive"
source("main.R")
```

Results are written to a timestamped subfolder under `output/` (e.g. `output/2026-02-19_1/`).

## Input data formats

| Instrument | Format | Key columns |
|------------|--------|-------------|
| **FTIR** (PerkinElmer Spotlight) | CSV (comma-separated) | `Coord. [um]` (bracketed X;Y), `Group` (material), `Max AAU score` |
| **Raman** (Horiba / WITec) | CSV (semicolon-separated) | `Visual Center Point X/Y [um]`, `Material`, `HQI`, `Feret Max [um]` |
| **LDIR** (Agilent 8700) | XLSX | `Identification` (material), `Quality`, `Width/Height [um]` — **no X/Y coordinates** |
| **FTIR image** | PNG | Chemical absorption map exported from the FTIR instrument |
| **LDIR image** | PNG | Colored particle map exported from Agilent LDIR software |

The ingest module auto-detects CSV delimiters and handles encoding differences (Latin-1, UTF-8).

## Pipeline architecture

```
main.R                          # Orchestrator
R/
  00_config.R                   # All tunable parameters
  00b_file_input.R              # Interactive file picker
  01_ingest.R                   # FTIR + Raman CSV/Excel parsing
  01b_ingest_image.R            # Image-based particle extraction (FTIR)
  01c_ingest_ldir.R             # LDIR Excel parsing + image coordinate extraction
  02_prefilter.R                # Quality/size pre-filtering
  03_normalize.R                # Coordinate centering
  03b_landmark_align.R          # Tier 1: large-particle landmark alignment
  04_ransac.R                   # Tier 2: RANSAC with rotation grid search
  05_transform.R                # Affine transform helpers
  06_icp_refine.R               # Iterative Closest Point refinement
  07_match.R                    # Nearest-neighbor spatial matching
  08_agreement.R                # Material agreement scoring
  08b_material_map.R            # Polymer family classification & equivalence mapping
  09_diagnostics.R              # Overlay plots, histograms, confusion matrix
  10_export.R                   # CSV/PDF/TXT export
  utils.R                       # Logging, geometry, transform utilities
  utils_python.R                # Reticulate bridge for Python LDIR detector
inst/
  python/
    particle_detector.py        # Python LDIR particle detection backend
shiny_app/
  app.R                         # Interactive viewer
  global.R                      # Data loading for Shiny
```

### Alignment strategy (FTIR ↔ Raman)

The coordinate systems from FTIR and Raman differ in origin, rotation (often ~180°), and sometimes scale. The pipeline finds the spatial transform in three tiers:

Alignment is **geometric** — it uses position, size, and shape, never material identification. (Material is compared only afterwards, in agreement scoring, so the geometry never depends on the identifications it is used to validate.)

1. **Tier 1 — Landmarks**: Particles ≥ 100 µm and fibers (aspect ratio ≥ 3) are matched between instruments. If no particle clears 100 µm, the 12 largest (≥ 30 µm) are used instead — landmarks are conspicuous *relative to the sample*, which matters on field samples where nothing is 100 µm. If enough landmarks agree (≥ 50% inliers, residual < 50 µm), Tier 2 is skipped.

2. **Tier 2 — RANSAC + global registration**: Two aligners run and the one pairing more particles wins. The coarse RANSAC does a rotation grid search (1° steps, including mirror check) and refines a similarity transform. Global registration sweeps rotation × scale and recovers translation by voting over all pairwise offsets, scoring one-to-one — more robust when the two clouds overlap only partially.

3. **ICP refinement**: Iterative Closest Point polishes the transform using all particles ≥ 20 µm, with reciprocal nearest-neighbour filtering and 10% trimming of worst pairs.

The same RANSAC + ICP pipeline is applied separately to align **LDIR → Raman**.

`align_ftir_materials` / `align_raman_materials` optionally narrow the Tier 2 anchor set to polymers both instruments should agree on. This is a **hint for spiked lab samples, not a requirement**: the filter is honoured only while it leaves at least `align_min_anchor_count` particles, otherwise the full cloud is used and the log says so. Set them to `NULL` for field samples.

### LDIR coordinate extraction

Since the Agilent 8700 LDIR does not export particle coordinates, the pipeline extracts them from the companion PNG image:

**Primary method (Python, preferred)**:
- Per-tile iterative sigma-clipping background correction (4×4 grid, scipy Gaussian)
- Global thresholding after background subtraction (not per-tile, because ~3% of particles span tile boundaries)
- Auto-tuned threshold — binary search to match the expected count from the Excel file
- Achieves ~99.6% join rate between image-detected and Excel particles on tested datasets

**Fallback (R, no Python)**:
- HSV saturation segmentation for colored particle maps
- Adaptive threshold for grayscale images

After detection, image pixel centroids are mapped to physical µm coordinates using the configured scan diameter. An optional scan-order validation (Kendall τ correlation) checks that particle IDs follow the expected raster pattern.

### Raman particle filtering by pipeline step

| Step | Size filter | Material filter | HQI filter |
|------|-------------|-----------------|------------|
| Spatial transform (landmarks, ICP) | ≥ 20 µm | none | none |
| Tier 2 alignment anchors | ≥ 20 µm | PET / PP if ≥ 4 available, else none | ≥ 70 if ≥ 4 available, else none |
| Spatial matching | all | none | none |
| Agreement scoring | all | none | ≥ 70 |

### Three-way triplets

Particles matched by all three instruments are identified by finding the intersection of FTIR↔Raman pairs and LDIR↔Raman pairs that share the same Raman particle ID. Each triplet includes:

- `n_instrument_agreement`: count of pairwise family agreements (0–3)
- `material_consensus`: human-readable agreement level ("Full agreement", "FTIR+Raman agree", etc.)
- `ftir_family`, `raman_family`, `ldir_family`: polymer family classifications for each instrument

## Output files

Each run creates a timestamped folder (e.g., `output/2026-02-19_1/`) containing:

| File | Description |
|------|-------------|
| `matched_particles.csv` | All matched FTIR↔Raman pairs with coordinates, materials, distances |
| `unmatched_ftir.csv` | FTIR particles with no Raman match |
| `unmatched_raman.csv` | Raman particles with no FTIR match |
| `unmatched_ldir.csv` | LDIR particles with no Raman match |
| `agreement_summary.csv` | Per-material agreement rates (Exact / Family / Disagree) |
| `agreement_pairwise.csv` | Per-pair material comparison with tier scoring |
| `ldir_raman_matched.csv` | LDIR↔Raman matched pairs |
| `ldir_ftir_matched.csv` | LDIR↔FTIR matched pairs |
| `ldir_raman_agreement.csv` | LDIR↔Raman material agreement |
| `triplets_3way.csv` | Particles matched across all 3 instruments with material consensus |
| `ldir_image_extracted.csv` | Raw image-extracted LDIR coordinates (before Excel join) |
| `ldir_coord_quality.csv` | Join quality metrics and scan-order validation |
| `composite_matches.csv` | 1:many FTIR→Raman composite matches (fragmented particles) |
| `transform_params.txt` | FTIR→Raman 3×3 transform matrix + ICP stats |
| `ldir_transform_params.txt` | LDIR→Raman transform + ICP stats |
| `match_statistics.txt` | Summary statistics |
| `tps_assessment.txt` | Thin-plate spline distortion assessment |
| `triage_top_pairs.csv` | Highest-confidence matches for manual review |
| `plots/overlay.png` | Spatial overlay of aligned FTIR + Raman |
| `plots/ldir_overlay.png` | Spatial overlay of aligned LDIR + Raman |
| `plots/confusion.png` | FTIR↔Raman material confusion matrix |
| `plots/ldir_confusion.png` | LDIR↔Raman material confusion matrix |
| `plots/distance_hist.png` | Match distance distribution |
| `plots/icp_convergence.png` | ICP RMS history |
| `plots/residual_scatter.png` | Spatial residuals per particle |
| `plots/residual_vectors.png` | Residual vector field |
| `plots/size_comparison.png` | Matched pair size comparison |
| `plots/tiered_agreement.png` | Tiered agreement bar chart |
| `plots/ambiguity.png` | Ambiguity score distribution |
| `plots/all_diagnostics.pdf` | Combined PDF of all plots |

## Shiny viewer

Launch the interactive viewer after running the pipeline:
```r
shiny::runApp("shiny_app/")
```

The viewer includes:
- **FTIR tab**: Native coordinate view with chemical absorption image overlay
- **Raman tab**: Native coordinate view
- **LDIR tab**: Native coordinate view with raw / processed / extracted-points image layers
- **Overlay tab**: All instruments in the shared Raman coordinate frame, with layers for individual instruments, matched pairs, and triple matches (gold rings)

### PDF report

The **Summary** tab has a **Download PDF report** button that writes one
multi-page PDF containing:

| Page | Content |
|---|---|
| 1 | Title page — run ID, generation timestamp, particle scope, per-instrument counts |
| 2 | Material comparison bar chart for the selected family |
| 3 | Plastics by Instrument (family × device counts, with totals) |
| 4 | Material breakdown pie charts, 2×2 |
| 5 | Size distributions per instrument |
| 6 | Size statistics table |
| 7+ | Each instrument view — image with its particles overlaid (FTIR PerkinElmer, FTIR Bruker, Raman, LDIR) |
| last | Overlay — all instruments in the shared Raman frame |

The report is a **snapshot of what the viewer is currently showing**: the
Summary tab's scope toggle, the quality / size / material / match filters set on
each instrument tab, the pie display modes and the view rotations. Every figure
carries a caption recording the filters behind it, so a page lifted out of the
PDF still says what it is showing. Instruments with no data loaded are skipped.

Figures are the same ggplot objects the app renders and the tables come from the
same builders as the on-screen HTML (`report_plastics_table()` /
`report_size_stats_table()` in `shiny_app/global.R`), so the report cannot drift
from the viewer. Output is written with `grDevices::cairo_pdf()` — no pandoc or
LaTeX needed.

## Configuration

All parameters are in `R/00_config.R`. Key settings:

```r
# Alignment anchors — a HINT, not a requirement. Alignment is geometric; these
# only narrow the anchor set when the sample happens to contain these polymers.
# Each filter is honoured only while it leaves >= align_min_anchor_count
# particles, otherwise the full cloud is used. Set to NULL to disable material
# anchoring entirely (recommended for field samples).
align_ftir_materials      = c("PET", "Polypro")
align_raman_materials     = c("Polyethylene terephtalate", "Polypropylene")
align_ldir_materials      = c("Polyethylene terephthalate", "Polypropylene", "Polycarbonate")
align_min_anchor_count    = 4      # below this, the material filter is dropped
align_raman_min_size_um   = 20     # Raman particles below this are excluded from alignment

# Quality thresholds
raman_hqi_threshold       = 70     # HQI cutoff for material agreement scoring
ldir_quality_threshold    = 0.6    # Agilent quality score threshold

# LDIR scan geometry — CRITICAL: must match actual scan area
ldir_scan_diameter_um     = 13000  # Physical extent of the LDIR scan area (µm)
                                   # Increase/decrease if ICP RMS > 100 µm in log

# Landmark alignment
landmark_min_size_um      = 100    # Particles >= this are landmark candidates
landmark_fiber_aspect_ratio = 3.0  # Fibers detected by this aspect ratio threshold
landmark_adaptive_size    = TRUE   # If < landmark_min_count particles clear the
                                   # absolute threshold, take the largest ones
                                   # instead — landmarks are conspicuous
                                   # relative to the sample, not in absolute µm
landmark_target_count     = 12     # How many to take in that fallback
landmark_min_size_floor_um = 30    # ...but never anything below this

# Matching
match_dist_threshold_um   = 100    # Max distance for a valid spatial match (µm)

# RANSAC
ransac_inlier_dist_um     = 200    # Inlier distance threshold (µm)
ransac_allow_mirror        = TRUE  # Search reflections (needed for 180° rotations)
ftir_use_global_register  = TRUE   # Tier 2 also runs global registration
                                   # (rotation x scale sweep + translation
                                   # voting); the transform pairing more
                                   # particles wins
```

## Known limitations

- **LDIR scan area**: `ldir_scan_diameter_um` defaults to 13,000 µm (13 mm filter). If the LDIR software scans a smaller or differently-shaped area, adjust this setting. The pipeline emits a WARN log when LDIR ICP RMS > 100 µm, which is the primary diagnostic for a mismatched scan area.

- **LDIR material anchor count**: The LDIR RANSAC alignment uses PET, PP, and PC particles as anchors. If fewer than 4 such particles exist in the dataset, all particles are used, which can produce a poor initial transform for the alignment step.

- **FTIR feret size**: The PerkinElmer Spotlight CSV does not export a dedicated Feret Max column. The pipeline uses `Major dim` as a proxy, which may overestimate actual Feret dimensions for non-convex particles.

- **Raman minimum size**: The Raman instrument detects particles down to ~1 µm, but FTIR resolution is typically ≥ 20 µm. Raman particles below 20 µm have no FTIR counterpart and inflate the "unmatched Raman" count.

- **Material name mapping**: FTIR, Raman, and LDIR use different spectral libraries with different naming conventions (e.g., FTIR "Polypro" ↔ Raman "Polypropylene (PP)" ↔ LDIR "Polypropylene"). The pipeline normalises common abbreviations and uses polymer family classification. Extend the mapping in `R/08b_material_map.R` for dataset-specific names.

- **"Polyamide (naturally occurring)" from LDIR**: This LDIR material name refers to biological polyamides (silk, protein-like), not synthetic nylon (PA6, PA12). It is classified into the "Natural" polymer family (not "PA") to avoid false material-agreement matches between synthetic and biological materials.

- **Single-filter assumption**: The pipeline assumes all instruments analysed the same physical filter. If scan areas differ significantly, the alignment may fail or produce a poor transform.

## License

Released under the [MIT License](LICENSE) © 2026 Kristof Moeller.
