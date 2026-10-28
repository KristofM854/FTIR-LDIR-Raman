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
  - distance gate 50 um (was 100); for LDIR see below;
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
- Landmarks need 5 one-to-one pairs and the chance test (p < 1e-8);
  the share of landmarks that pair (>= 50%) no longer suffices on its own: 5
  of 6 random landmarks on a dense reference pair (p = 0.02). An interim
  minimum of 8 pairs rejected 7 pairs at p = 2e-13 on a reference filter.
- **(results)** LDIR gate 20 um (was 250): ~3 sigma of the LDIR positional
  error, which every run now fits and logs (sigma 7.1 um dataset A, 7.3 um
  field sample A; `position_error_window_um`, `chance_pairs.csv`, report).
  On field sample A the 75 um gate gave ~41% wrong pairs, 20 um ~5%. Pairs of
  two aligned instruments take the wider gate (50 um). Dataset A: LDIR-Raman
  29 -> 28, FTIR (Bruker)-LDIR 18 -> 14, triplets 29 -> 28.
- **Pair confidence** `p_true` on every matched pair (matched and agreement
  files): the probability that the pair is the true partner and not a chance
  neighbour, from the positional-error fit (larger particles: sigma at least
  0.05 x major axis). Agreement is also reported over pairs with p_true >= 0.9
  (`pair_confidence_high`). The viewer can colour matched particles by it.
  Synthetic field sample: expected wrong pairs 2.5, actual 2.
- **95% intervals** (Wilson) on every agreement rate: log, match_statistics,
  report.
- **Detection rates by size** (`07_diagnostics/detection_rates.csv`, two
  report pages, a Summary-tab plot): for every instrument pair, both ways, the
  share of one instrument's particles (quality-passing, >= 20 um, inside the
  other's covered area) that the other also found, per size class 20-50,
  50-100, 100-200, > 200 um, all particles and plastics, with 95% intervals
  and a chance-corrected rate (sum of p_true).
- Test: synthetic dense field sample end to end (`test-pipeline-field.R`).
- **(results on field samples)** The LDIR gate is set per run: 3 x the
  positional error fitted to that alignment, within 10-50 um
  (`match_gate_from_sigma`, `match_gate_sigma_k`, `match_gate_min_um`,
  `match_gate_max_um`); the fixed 20 um when there is no fit. With the
  analysed LDIR image sigma is ~3 um (field sample A), so the gate is 10 um and
  the chance pairs drop to about a quarter. Dataset A unchanged (21.2 um).
- **(results on field samples)** LDIR analysed image: the scale stays the scan
  diameter's unless the Excel areas imply one more than 1.5x off
  (`ldir_processed_scale_tolerance`; the overlay then covers part of the scan,
  dataset A 2.7x). The area estimate was 13% off on field sample A, where the
  alignment had to correct it. A difference above 5% is a WARN.
- Detection rates in two senses: *identified* (paired with a particle passing
  the other's quality cutoff, as before) and *detected* (any particle of the
  other instrument within the pair's gate, of any quality, with its
  chance-corrected rate); in `detection_rates.csv`, the report (four pages)
  and the Summary plot.
- Segmented PerkinElmer exports (temporary, `R/ftir_segments.R`): when the
  selected FTIR table is one of six segments, the run stitches the tables and
  images (2 x 3 grid column by column, 2 um/px) into `inputs/ftir_stitched/`
  and continues with them; `tools/stitch_ftir_segments.R` does it on its own.
  Checked on field sample A (15 plastics land on particles of the stitched
  image). The stitched table numbers its particles P_1 .. P_n (the export labels
  every non-plastic MP_-1, so its identifiers repeat) and keeps the exported
  identifier and the segment in `Identifier (export)` and `Segment`. A single
  PerkinElmer export with repeated identifiers is renumbered the same way on
  reading (`particle_id` P_1 .. P_n, the original in `identifier_export`).
  It writes the image placement next to the
  image (`<image>.txt`), which the pipeline reads for PerkinElmer and Bruker
  images when the config sets no extent.
- Viewer: hovering (and clicking) on the overlay finds only the particles the
  overlay draws -- its material, quality and size filters -- not every particle
  of the run, and a hovered particle that a filter change hides is cleared, on
  the overlay and on each instrument's view.
- FTIR (and any other non-LDIR) pairs against the reference take the per-run
  gate (3 x the fitted positional error, 10-50 um) when the fixed 50 um gate
  would let in more than `chance_share_warn` (20%) chance pairs, and keep
  50 um otherwise: FTIR has a tail of true pairs at 4-6 sigma (dataset A: 3 of
  49 at 41-48 um), which a 3 sigma gate on a sparse filter would lose for
  nothing. Field sample A, all six PerkinElmer segments unfiltered: sigma
  12.4 um, 39% chance pairs at 50 um, so ~37 um. LDIR keeps its per-run gate.
  Aligning FTIR to LDIR instead of Raman (open decision) is dropped: with the
  six segments the FTIR-Raman landmarks are confident (35 pairs, p = 4e-40).
- Agreement: a pair in which either identification is a reference chemical
  (family "Other compound": a drug, solvent, oil, reagent from the spectral
  libraries) is Unresolved, not Disagree, and two different such chemicals no
  longer count as a family hit. The Raman libraries hold almost no
  natural-fibre references, so natural particles take a chemical's spectrum:
  field sample A, FTIR (PerkinElmer)-Raman family agreement 26% (13-45%) of
  27 resolved pairs instead of 7% of 107, with 81 unresolved. Without the
  licensed library index no name is "Other compound" and nothing changes.
- Viewer overlay, multi-instrument option: a particle is ringed when the
  shown pairs -- with Raman or direct (FTIR (Bruker) <-> FTIR (PerkinElmer),
  ... <-> LDIR) -- link it into a group of three or more instruments
  (`link_groups()`). Before, only PerkinElmer + Raman + LDIR triplets were
  ringed, never the Bruker particle, and triplets of large particles whose
  LDIR leg lay in the gate widened for their size lost the ring (dataset A:
  4 of 28). Hover and click also find particles paired only directly (e.g.
  FTIR + LDIR without Raman) and show the whole group; Clear also empties the
  per-instrument highlight selections.
- Viewer: the Summary detection plot calls `shiny::validate()` explicitly; with
  jsonlite attached in the session its `validate()` masked Shiny's and the
  page failed ("is.character(txt) is not TRUE").
- `docs/roadmap.md` (open work and open identification questions) replaces
  the internal pre-release to-do list; `docs/ai_assistance.md` re-measured.
- A Tier 2 (anchor) alignment is "not determined" unless, after ICP, more
  anchors pair within the gate than the pose search finds by chance
  (`align_determined_p` 0.01); its pairs are written but not scored, pairs
  through it (LDIR-FTIR, triplets, Bruker-PerkinElmer) are skipped, and the
  report says so. Field sample A: both FTIR exports (pre-filtered to 15 and
  5 plastics) cannot fix a pose.
- The LDIR "tiled/repeated" image warning compares the quadrants' detail
  (correlation after removing smooth shading) instead of their mean and SD,
  which matched on every evenly lit filter image and warned on every run.
- Viewer: material filters have "plastics" (synthetic + semi-synthetic, as
  in the Summary tab) and "synthetic" quick selections; the boxes no longer
  appear unticked on load; no Raman image-scale warning when the image is
  placed from the scan's own WITec placement file (the particle-area estimate
  reads ~3x too coarse on a downsized field-sample image).
- The WITec placement file is also found by its content (any .txt next to
  the image with Width, Height, X, Y; of several, the one that puts the most
  particles on the image).
- **(results)** Material agreement: family agreement is primary, exact polymer
  secondary; only identifications above each instrument's cutoff (FTIR AAU
  0.70, Raman HQI 70, LDIR 0.85) and not "unknown" are scored; unresolvable
  name pairs form an *Unresolved* tier; rates are over resolved pairs.
  Triplets apply the same gate.
- **(results)** Acrylic polymers (incl. polyacrylonitrile) are *Synthetic*.
- The PerkinElmer and Bruker FTIR paths share one `align_and_match()`;
  `main.R` is file input and config only.
- Python is required for LDIR images (no R fallback); dependencies declared
  in `DESCRIPTION` and `inst/python/requirements.txt`.
- The LDIR coordinate-swap file is found under any name containing
  coord/cord and swap.
- **(results on other data)** The landmark tier counts one-to-one landmark
  pairs (`landmark_min_pairs`, see above); the race between RANSAC and global
  registration is decided on one-to-one pairs. Dataset A unchanged.
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
