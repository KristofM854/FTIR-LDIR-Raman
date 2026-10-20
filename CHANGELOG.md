# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Entries that
change results on existing data are marked **(results)**.

## [1.0.0] — first archived release

Pre-release audit and portability refactor (branch `audit/pre-release`).

### Added
- **Any combination of two to four instruments**: the shared frame is Raman
  when present, else FTIR (PerkinElmer), else FTIR (Bruker); every other
  instrument is aligned into it and every pair is matched. Output files keep
  their names and columns (`matched_<a>_<b>.csv`; new: `matched_ldir_ftir_bruker.csv`).
  `04_alignment/aligned_<instrument>.csv` holds every instrument in the shared
  frame; the viewer and the report show any combination. The file picker
  requires no particular instrument, only two data files.
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
- **(results) Field-sample quality control** (dense Raman fields: 12 000+
  particles put a chance partner within 150 um of every LDIR particle):
  - only particles passing their instrument's quality cutoff are matched
    (`match_quality_filter`; the agreement cutoffs);
  - distance gates 50 um (was 100), 75 um with LDIR (was 150);
  - Tier 2 alignment anchors are the quality-passing plastics
    (`align_anchor_plastics`);
  - landmarks are trusted also when their pairs are far above chance
    (Poisson p < 1e-8, `landmark_significance_p`), not only at >= 50% inliers;
  - RANSAC and global registration search scales 0.8-1.25 only (the span
    ratio was read as a scale: up to 3.8 for an FTIR scan of part of a
    filter, giving wrong poses and most of a field-sample run's time);
  - every run logs the expected share of chance pairs per instrument pair
    (`07_diagnostics/chance_pairs.csv`);
  - only particles >= 20 um (Feret max) are matched (`min_particle_size_um`,
    was 0) and a pair whose sizes differ more than 5x is not made
    (`match_max_size_ratio`). Dataset A: pairs 31/20/30 -> 31/18/29, Bruker
    with PerkinElmer 17 -> 16, with LDIR 19 -> 18; the five lost pairs are a
    Bruker particle sized 18 um (110-120 um on the other three) and two
    Raman particles sized 6-7 um next to a 23 um and a 318 um particle.
  Dataset A: transforms unchanged; pairs 32/29/34 -> 31/20/30 (FTIR PE,
  FTIR Bruker, LDIR with Raman), scored pairs 32/21/31 -> 31/20/30, 100%
  agreement; chance share 2-5%.
- Landmarks need 5 one-to-one pairs (was 8) besides the chance test: 8 rejected
  7 pairs at p = 2e-13 on a reference filter (Tier 2 then found the same pose).
- The LDIR "tiled/repeated" image warning compares the quadrants' detail
  (correlation after removing smooth shading) instead of their mean and SD,
  which matched on every evenly lit filter image and warned on every run.
- The WITec placement file is also found by its content (any .txt next to
  the image with Width, Height, X, Y; of several, the one that puts the most
  particles on the image).
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
- **(results on large samples)** Matching is optimal (Hungarian) for any
  number of particles; it used to switch to greedy above 2,000 particles on
  the smaller side.
- Explicit anchor polymers for Tier 2 alignment (`align_*_materials`) are
  off by default (`NULL`); see the field-sample quality control above for the
  default Tier 2 anchors.
- Triplets carry `direct_ldir_ftir_match`; approximate image placements are
  labelled in the viewer and report; debug trace IDs no longer default to
  dataset A's particles.

### Fixed
- Optimal matching built one dense n × n cost matrix (n = the larger
  particle set): 1.7 GB and most of the run time for a 14,500-particle
  Raman field sample. It now solves each group of particles that can pair
  within the gate on its own; the result is identical (same pairs and
  scores; checked against the dense solve), matching takes seconds.
- Viewer: run images were read at full resolution; a field-sample
  micrograph (or the original file, when the run's downsized copy could not
  be found) took ~10 GB of memory per 400 MB image and minutes to draw.
  Images are now resized by ImageMagick to at most 3000 px before they are
  read into R (same placement); a moved run folder uses its own `inputs/`
  copy before the original. Sample B (12000 px): 11 GB -> 2 GB, viewer ready
  in 16 s instead of 96 s.
- Raman images of 1-2 GB were not downsized: ImageMagick refused them under
  its default resource limits ("cache resources exhausted"), the run went on
  with only a WARN, and the viewer then loaded the original. Large BMPs are
  now downsized strip by strip from disk (`read_bmp_downsized()`, < 1 GB of
  memory, ~30 s for 1.4 GB), in the pipeline and in the viewer.
- Viewer: "Ignoring unknown parameters: label.size" with ggplot2 >= 3.5.
- Report: a "Matching Quality per Instrument Pair" page (gate, pairs,
  distance median / p90, pairs expected by chance, scored pairs, family
  agreement); FTIR (Bruker) in the size distributions and statistics; the
  title page no longer runs off the page; the input-file list fits.
- The Raman image auto-calibration cross-check (with a WITec placement file)
  is off by default (`raman_placement_cross_check`): 7 min on a field sample,
  and it does not change the placement.
- Viewer: the "Highlight particle" lists sent every particle ID to the
  browser, which froze it on field samples; they now load on demand
  (type to search).
- Viewer and report: the instrument views said "matched to Raman / FTIR"
  whatever the run's instruments; they name the run's partners.
- An Excel export whose sheet is not called `Long_Table` failed ("Sheet
  'Long_Table' not found"); a workbook with a single sheet is now read from
  it, and with several sheets the error lists them.
- **(results)** The PerkinElmer library label "PolyEthylen" was classified
  "Other polymer"; `inst/materials/ftir.csv` now translates it to
  polyethylene (PE).
- Viewer: a missing table was silently replaced by another one whose name
  starts the same (`$` partial matching on R lists): e.g. with no unmatched
  FTIR (PerkinElmer) particles, `unmatched_ftir` returned the unmatched FTIR
  (Bruker) particles, shown as PerkinElmer ones. Every table now has an
  explicit entry.
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
