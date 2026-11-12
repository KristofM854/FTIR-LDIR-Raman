# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Entries that
change results on existing data are marked **(results)**.

## [1.0.0] — first archived release

Pre-release audit and portability refactor (branch `audit/pre-release`).

### Fixed (audit, branch `audit/fable-fixes`; see AUDIT.md)
- M-01: `detection_rates()` no longer stops ("missing value where TRUE/FALSE
  needed") when `min_particle_size_um` is below 20 and a particle is smaller
  than the first size class; such particles count in "all sizes" only.
- M-02: a run in which every particle of one instrument is below its quality
  cutoff (or no particle is left to match) no longer stops in the matcher
  ("no points in query!"); the pair is written with 0 matches.

### Added
- **Run Info: supporting files.** A table of the side files the pipeline looks
  for next to each instrument's data (FTIR scan segments and placement .txt,
  Bruker placement / OPUS values, Raman WITec placement .txt, LDIR analyzed
  image, .sam project file, coordinate swaps): how to name them, whether this
  run found them, and what they do (or the fallback without them). Runs record
  the ones they used in `00_manifest/supporting_files.csv`.
- Overlay: the four instrument filter cards are one height, with a fixed-height
  material list, so their controls line up whether or not a list scrolls; the
  ID-pattern box and its Apply button sit on one line (also on the instrument
  tabs).
- **Viewer reorganised around one question per tab.** The instrument tabs show
  what each instrument found; the Overlay and the Summary tab show how the
  instruments agree.
  - Instrument tabs: match status removed (Match Status boxes, "Show all
    detected", pair-confidence colouring, LDIR match-score slider, "Show Raman
    partners", "Hide unmatched"). Its "matched" meant matched to Raman (to FTIR
    PerkinElmer on the Raman tab) and so differed from tab to tab. Points are
    now filled above the instrument's quality cutoff and open below it; the
    hover table says which.
  - Summary tab, in three parts: (1) **alignment check** per instrument
    (method, scale, rotation, mirror, residual, positional error, pose by
    chance, verdict; new pipeline output `04_alignment/alignment_summary.csv`);
    (2) what each instrument found (as before; size statistics now count
    particles above / below the quality cutoff instead of matched / unmatched);
    (3) how the instruments agree: **pair matrix** (identified / detected, each
    direction, with 95% intervals), matching quality per pair, **why particles
    are not matched** (from `partner_status.csv`), **particle groups** (which
    combinations of instruments found the same particle, as in the overlay's
    symbol view) and **material concordance** per instrument pair; then the
    detection rates by size.
  - Report: the same alignment and agreement pages, the concordance tables as
    an appendix.
  - Removed: the Upload Data tab (loaded the old CSV bundle from before run
    folders) and the overlay's "Multi-instrument (3+/4)" option (the symbol
    view's combinations show it).
- **Partner status, independent of quality** (`05_matches/partner_status.csv`,
  `partner_status()` in `R/detection_rates.R`): for every particle and every
  other instrument, what that instrument has at its position whatever the
  quality: matched, partner below its quality cutoff, this particle below
  its cutoff, both below, partner paired elsewhere, nothing found, or outside
  the partner's area; with the partner's ID, material, quality, distance and
  the chance of an unrelated particle within the gate. Matching itself is
  unchanged (quality-passing particles only). The run log counts the statuses
  per pair. In the viewer's overlay, the "Match Details" table fills an
  instrument with no match from this file (in italics, with a note such as
  "below the HQI cutoff: not confidently identified" and the distance), and
  unmatched particles get a ring coloured by status in both views
  (checkbox, on by default).
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
- Viewer overlay: a symbol view (default) draws each group of linked
  particles once, as a circle in four fixed quarters (top-left FTIR
  (PerkinElmer), top-right FTIR (Bruker), bottom-left Raman, bottom-right
  LDIR), filled where that instrument found the particle; colours from the
  Okabe-Ito palette, checked for deuteranopia, protanopia and tritanopia.
  A combination list with counts filters the groups (e.g. only
  PerkinElmer + LDIR). The former view stays as "Raw positions".
- Viewer: one colour-blind-safe instrument palette (Okabe-Ito) across the whole
  app and the report -- FTIR (PerkinElmer) yellow, FTIR (Bruker) reddish
  purple, Raman blue, LDIR vermillion -- one exact colour per instrument
  everywhere; labels are badges filled with it (black or white text). Particles without a partner are grey on
  every instrument (formerly a second colour each, e.g. green vs red), the
  pair-confidence gradient is viridis (formerly red-orange-green), and the
  size histograms take the instrument colour. Overlay: the per-instrument
  filters sit side by side below the image, under the Match Details box;
  background image and offsets at
  the bottom of the sidebar. The overlay's LDIR match-gate slider is
  removed: LDIR pairs are counted with the run's own gate (the per-run gate
  the pipeline matched with). The Match Details table has fixed rows (Particle
  ID, Material, Quality, Feret Max, Position) and a column for each of the
  four instruments, dashes where one has no particle, so it keeps its height. A global
  microplastics choice in the sidebar (all particles / plastics / synthetic
  only) sets the four instruments' material lists at once. Material lists:
  checkbox and label centred on each other, and a short list centred in its
  box (a long one scrolls from the top).
- Viewer: the run list names each run's input folder (after "My Documents"),
  e.g. `2026-10-28_1 — Automatisations/MPs/.../Bora2025/A`.
- Precision study after ISO 5725-2 (`R/precision_study.R`,
  `tools/precision_study.R`, viewer tab "Precision"): the same filter
  measured n times on each instrument; particle counts per polymer family and
  size class; repeatability, between-instrument and reproducibility SDs with
  95% intervals, r / R limits, CVs, Mandel's h and k with indicator levels
  (checked against their null distributions), and the replicate variance
  against Poisson counting noise. CSVs and a PDF report.
- Replicate tool (`R/reproducibility.R`, `tools/reproducibility.R`), four
  fixes: run-to-run matching uses the pipeline's matcher (per group of nearby
  particles; was a dense n x n matrix, ~1.25 GB per pair for 12 500 Raman
  particles); registration starts from the pipeline's landmark search on the
  200 largest particles (scale 1, no mirror), so a remounted filter at any
  angle aligns (was: start at 0 deg); the gate is 3x the positional error
  fitted to each pair of runs, 10-75 um (was a fixed 75 um), with the chance
  share and p_true per link; links are merged most confident first and never
  join two particles of one run (was: union-find kept the first of two
  silently). Synthetic dense test (12 000 particles, 3 replicates at 0 / 90 /
  180 deg): rotations exact, 0.6% of groups impure, detected-in-all 8 671 vs
  8 691 true, 27 s. `--design design.csv` runs every instrument of a
  precision-study design.
- Precision study and replicate tool from pipeline runs: one run per
  replicate, `--runs <run folders>` for `tools/precision_study.R` and
  `tools/reproducibility.R`, and in the viewer's Precision tab (select runs,
  Compute). The runs' instrument tables are used (stitched PerkinElmer
  segments, LDIR joined to its image coordinates), the counting cutoffs come
  from the first run's config; no raw files to name.
- Precision study views: a group reported by one instrument only is marked
  "*" (an identification difference, not a precision) and left out of the
  Mandel h / k plots, which state their ceilings ((p - 1) / sqrt(p) and
  sqrt(p)); a short explanation under every table and figure, in the viewer's
  Precision tab and the PDF report; the counts-per-replicate plot is
  interactive (plotly: hover a dot for its replicate, a bar for mean and SD)
  with the Summary tab's switches (all particles / plastics only, synthetic
  only / synthetic + semi-synthetic); new total group "Synthetic".
- Precision report: the first page opens with a short results paragraph for
  the synthetic polymers (`precision_summary_text()`), written from the
  numbers: mean count, CV_r / CV_L / CV_R with r and R, the instrument share
  of the reproducibility variance, the instrument range, Mandel flags,
  counting-noise check and the size classes with the highest / lowest CV_R.
  Study folders written before the "Synthetic" group get it on reading,
  summed from the synthetic family counts (the same count) with its
  statistics, so their summary and tables are complete.
- Precision tab: the study selector names each study like the run list
  does: study folder, its runs (or number of export files) and the input
  folder of the data (`precision_study_labels()`).
- File input (`main.R`): a list of instruments first, so only their dialogs
  open (Cancel there ends the run), and a confirmation of the chosen files
  before the run starts (No ends it) -- no more cancelling through every
  remaining slot after a wrong file. The confirmation is a Tk window as wide
  as the longest file name (Windows message box when Tk is unavailable).
- Viewer: the Precision and Replicates tabs resolve their results folders
  against the project folder (`project_path()`), so a study computed in the
  viewer lands in the pipeline's `output/precision/`, not in
  `shiny_app/output/precision/`.
- FTIR / Bruker reader: coordinates written as "[x,y]" (no semicolon, one
  comma; an occasional Lumos export glitch) are read with the comma as the
  x/y separator, with a warning naming the file and the number of rows,
  instead of failing with "could not be parsed as [x;y]".
- Images: a 32-bit BMP whose alpha channel is all zero (ImageMagick reads it
  as fully transparent; suspected for the Bruker Lumos background image) is
  drawn opaque, in the viewer and in the pipeline's reader
  (`opaque_if_blank()`). Before, the image loaded but showed nothing.
- Report: input file names on the title page in full (on their own line,
  broken at a space when too long; a stitched PerkinElmer table also names
  its segments); wide tables (detection rates) are fitted to the page (long
  cells on two lines, then a smaller font); stacked bars no longer drop the
  natural and inorganic families -- they had no colour, and an NA fill is
  drawn as nothing. One family palette (`FAMILY_COLOURS`, `family_colour()`)
  replaces four copies, in the report and the viewer.
- FTIR (Bruker) image placement from the OPUS values of the scan
  (`opus_ftir_bruker_image`: pixel size, bitmap position, image size, map
  origin; pre-filled for Faanui SW1L R4, 2026-07-27). The Lumos export writes
  no metadata, so the values are copied from OPUS; they are checked against
  the image size and the particles and ignored with a warning when they belong
  to another scan. The particle-extent fallback drew this image 25% too small.
- Alignment pose limits per instrument (`align_pose_limits`): FTIR (Bruker)
  and FTIR (PerkinElmer) are searched without mirror and at scale 0.97-1.03,
  in every tier; the pose
  check counts mirrored poses only when they were searched.
- Alignment on image spots: when an FTIR-type particle list cannot fix its
  pose ("not determined"), the pose is fitted on the particle-like spots of
  its exactly placed image (background-flattened, within the particle extent)
  and kept when it passes the same chance test (`align_use_image_spots`,
  `image_spot_cloud()`, `align_and_match(image_cloud)`). The spots are
  fitted to the reference's largest particles of any material
  (`align_image_spots_ref_max`), not to its plastic anchors: a field sample's
  few plastics (52 Raman anchors) left the spots without partners.
- Image-spot alignment, checked on the Faanui field sample (Lumos photo,
  OPUS values, Raman): the pose search on spots is now a rotation sweep with
  translation voting on the 300 largest spots against the 1500 largest
  reference particles (`vote_align()`, `fit_pose_spots()`); RANSAC and the
  global registration (random 150-point subsamples, 200 um tolerance) found
  no pose in the dense field. Result: rotation -4.79 deg, scale 0.999, 231
  of 300 spots paired (p = 7e-128; LDIR -> Raman on the same filter -4.76
  deg).
- Lumos table fitted to its own image (`fit_table_to_image()`,
  `ftir_table_to_image_fit`): the table's y runs ~0.7% short of the image
  (it converts both axes with the map's x step) plus a small shift; per-axis
  scale and shift are fitted to the particles visible in the image and the
  table corrected into the image frame. Faanui: median distance of a Bruker
  particle to the nearest Raman particle 44 um -> 9 um; Bruker-Raman pairs
  6 -> 12. Confirmed on the photo: the table's y points up from the map
  origin.
- Viewer overlay: "Number all particles" puts each instrument's ID in its
  corner of the marker (the symbol view's corners) in the instrument's
  colour, so a matched pair's IDs no longer overprint (`add_overlay_labels()`).
- LDIR particle outlines: with the exact `.sam` join the run writes
  `02_joined/ldir_outlines.csv` (each particle's outline, native LDIR frame;
  `read_ldir_sam_particles()` keeps the vertices, `join_ldir_by_sam()`
  converts them with the overlay's pixel -> um map). The viewer draws them as
  filled shapes under the points in the LDIR tab (native and aligned) and in
  both overlay views (`ldir_outline_polygons()`: in the aligned frame each
  outline is anchored on its particle and turned with the native -> aligned
  map). Manual LDIR join corrections are not applied after the `.sam` join.
- Viewer: the four instrument tabs and the overlay open with only the
  synthetic material families ticked (all families when the data has no
  synthetic one); the overlay's microplastics choice defaults to "Synthetic
  only". The "all" / "plastics" links and the radio buttons widen it.
- LDIR join from the Clarity project file (`.sam`): the file (decoded from
  Clarity 1.7.27) stores each particle's ID, outline (pixels of the IR scan
  image = the analyzed overlay's grid, y from the bottom) and area / width /
  height. `read_ldir_sam_particles()` scans it for the particle records;
  `join_ldir_by_sam()` picks the particle group whose IDs and areas match
  the Excel table and gives each row the overlay blob its outline covers
  (`coord_source = "sam"`). Used automatically when the only `.sam` next to
  the LDIR Excel (or `ldir_sam_path`) fits; else the size / shape / colour
  join. Faanui: A1's outline 337-467 x 1884-1920 px = 650 x 180 um, area
  39300, the Excel values, on the overlay particle at the same pixels.
- LDIR join: colour classes of the analyzed image. The LDIR "analyzed"
  export colours every particle by its identification class (palette
  reused when there are more classes than colours). `ldir_colour_classes()`
  matches the colours to the table's classes by count (no legend needed;
  left-over classes allowed in the colours with surplus blobs; equal counts
  share) and the join adds a cost (`ldir_join_weight_colour_class`, 5) for a
  row on a blob of another class, in every pass. Faanui field sample: 10
  colours for 14 classes, counts exact (Undefined 946 = rows below quality
  0.85); without the term only 54% of the rows sat on a blob of their own
  class, with it 99.9%; LDIR-Raman family agreement where both identify
  9% (1 of 11) -> 41% (9 of 22); LDIR polyethylene on Raman LLDPE/paraffin
  8 of 8. The extraction records each blob's modal colour (`mode_rgb`).
- LDIR join check: new circularity statistic (LDIR table vs Raman partner;
  circularity is not used by the join, so unlike size it tests whether a row
  sits on the right particle, not just one of the right size; the Raman
  reader keeps the exported circularity). The family statistic reports the
  highest agreement a perfect join could reach (`max_possible`) and is
  informative only when that is at least 0.2 above the random join, so a
  sample where the instruments name different families (Raman reference
  chemicals against LDIR natural materials) no longer counts as a test of the
  join. A warning when sizes agree but circularity does not. Dataset A:
  circularity 0.75 (p < 0.001); Faanui field sample: 0.20 (p 0.04).
- Viewer: the FTIR (PerkinElmer) and FTIR (Bruker) tabs frame the whole
  background image as well as the particles (an exactly placed image is
  larger than the particles' extent and was cropped).
- Run: the Bruker image object is recorded for an OPUS placement (it failed
  with "subscript out of bounds").
- Viewer overlay, symbol view: the "Found by" legend failed ("`breaks` and
  `labels` must have the same length") once an instrument was hidden; it
  now labels only the shown instruments.
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
