# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Entries that
change results on existing data are marked **(results)**.

## [1.0.0] — first archived release

Pre-release audit and portability refactor (branch `audit/pre-release`).

### Added
- **Data contract**: every reader returns one canonical particle table
  (`R/particle_contract.R`), checked by `validate_particles()`; images are
  placed by one pixel ↔ µm mapping, recorded per run in
  `00_manifest/images.json`.
- **Column-mapping reader**: any particle table (CSV/Excel), with units,
  y orientation and decimal commas, via `config$<role>_column_map` (R list or
  JSON). An LDIR table with coordinates needs no image.
- **Per-instrument material tables** `inst/materials/<role>.csv` to translate
  instrument labels without editing code; `tools/material_vocabulary.R`.
- **`run_pipeline(config)`** and `load_pipeline()`: the whole run as a function
  of its config, without file dialogs; `config_overrides` when sourcing
  `main.R`.
- **WITec Raman image placement file** (`<image>.txt`) used automatically;
  large Raman images are downsized for the run and originals over 100 MB are
  no longer copied into every run.
- **LDIR join check** in every run (log line, `07_diagnostics/ldir_join_check.csv`).
- **Validation tools** (`tools/validation.R`, `examples/validate_example_A.R`):
  chance matches, gate and alignment sensitivity; `docs/validation.md`.
- **Example dataset A** (one filter, all four instruments) with a regression
  test, a synthetic end-to-end test with known transforms, a determinism
  test, and tests for the Python detector.
- Provenance: full parameter dump, R session and package versions, git commit
  (also outside a git checkout), input MD5s.
- `docs/input_format.md`, `docs/validation.md`, `docs/ai_assistance.md`,
  `CITATION.cff`, `.zenodo.json`, this changelog.

### Changed
- **(results)** Material agreement: family agreement is primary, exact polymer
  secondary; only identifications above each instrument's cutoff (FTIR AAU
  0.70, Raman HQI 70, LDIR 0.85) and not "unknown" are scored; unresolvable
  name pairs form an *Unresolved* tier; rates are over resolved pairs.
  Triplets apply the same gate.
- **(results)** LDIR distance gate 150 µm (was 250 µm).
- **(results)** Acrylic polymers (incl. polyacrylonitrile) are *Synthetic*.
- The PerkinElmer and Bruker FTIR paths share one `align_and_match()`;
  `main.R` is file input and config only.
- Python is required for LDIR images (no R fallback); dependencies declared
  in `DESCRIPTION` and `inst/python/requirements.txt`.
- The LDIR coordinate-swap file is found under any name containing
  coord/cord and swap.
- **(results on other data)** The landmark tier is trusted only with at
  least 8 one-to-one landmark pairs (`landmark_min_pairs`); the race between
  RANSAC and global registration is decided on one-to-one pairs. Dataset A
  unchanged.
- **(results)** Composite (1:many) matches no longer share a Raman fragment.
- Material anchors for Tier 2 alignment are off by default
  (`align_*_materials = NULL`); alignment is purely geometric unless set.
- Triplets carry `direct_ldir_ftir_match`; approximate image placements are
  labelled in the viewer and report; debug trace IDs no longer default to
  dataset A's particles.

### Fixed
- **(results)** Images read through magick had scrambled pixels (wrong array
  layout), which misplaced the FTIR background image.
- **(results)** Coarse RANSAC could drop scale 1 from its candidates when a
  few particles lay outside the shared field, then align 4% off scale.
- Fallback translation of the coarse pose; deterministic subsampling in
  alignment and diagnostics.
- Readers fail on inputs that would be read into wrong numbers (decimal
  commas, pixel instead of µm coordinates, duplicate IDs, empty files).
- The run folder is byte-reproducible (PNG timestamps removed).

### Removed
- The LDIR scan-order test (could not detect a bad join), the R fallback for
  LDIR detection, the unused `material_map_*` config keys, the rotation
  warning, `normalize_coordinates()`.

## [0.1.0] — development, February–October 2026

FTIR (PerkinElmer, Bruker) / Raman / LDIR ingestion, LDIR particle detection
from images, tiered alignment (landmarks, RANSAC, global registration, ICP,
TPS), one-to-one matching, tiered material agreement with family rules and an
optional spectral-library index, diagnostics and report, Shiny viewer with
multi-run comparison and reproducibility tools. History in git (pull requests
#1–#207).
