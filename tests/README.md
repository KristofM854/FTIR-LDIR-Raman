# Test suite

The project is run by sourcing modules (`R/load.R`, `main.R`), not as an
installed package, so the tests source the `R/` modules they need directly.

## Running

```r
# from the repository root
Rscript run_tests.R
```

About 5 minutes. Needs the packages in `DESCRIPTION` (incl. `testthat`,
`withr`, `magick`) and Python with `inst/python/requirements.txt`; tests that
need Python skip without it, and `SKIP_E2E=1` skips the dataset A run.
`.github/workflows/tests.yml` runs the suite on pushes to `main` and on pull
requests.

The Python detector's own tests also run on their own:

```sh
python -m unittest discover -s inst/python -p "test_*.py"
```

## What is covered

**End to end**

| File | What it checks |
|---|---|
| `test-pipeline-synthetic.R` | `run_pipeline()` on a synthetic filter seen by all four instruments, each read through the column-mapping reader with a *known* transform (178° rotation in mm; 3° with y-down; −90° with scale 1.03; particles outside the shared field). Recovered transforms, every matched pair a true correspondence, the exact agreement tiers and exclusions designed in, triplets, output files. Then the same run twice: the output folder is byte-identical apart from timestamps, the run path and PDF creation dates. |
| `test-pipeline-field.R` | A synthetic *dense* field sample (12 000 Raman particles on a 13 mm filter, LDIR with a known 7 µm positional error, five unrelated FTIR particles): LDIR pairs correct and found at the 20 µm gate, the positional error recovered, a low chance share, pair confidence calibrated (expected vs actual wrong pairs), detection rates written, and the short FTIR list reported as not determined. |
| `test-example-dataset-A.R` | The bundled real dataset (`examples/dataset_A`) through the full pipeline in a subprocess, compared with a stored fingerprint (row counts, every matched pair and its distance, transforms, agreement tier counts). |

**Data contract and readers**

| File | What it checks |
|---|---|
| `test-particle-contract.R` | Canonical particle table and image object; `validate_particles()` / `validate_image()` list every problem; pixel ↔ µm mapping in both y orientations. |
| `test-readers.R` | Column-mapping reader (units, y orientation, decimal comma, JSON and Excel), material tables, `read_instrument()` for every role, and that vendor column names appear only in the reader files. |
| `test-ingest-contract.R` | Built-in readers: decimal commas, required columns, empty files, duplicate IDs. |
| `test-raman-image-metadata.R` | WITec placement file (corner → centre, rotated / incomplete files ignored, discovery next to the image), large images downsized and not copied, the viewer honouring a stage-frame centre. |
| `test-image-read.R` | Image reading (magick bitmap layout). |

**Alignment and matching**

| File | What it checks |
|---|---|
| `test-align-and-match.R` | Raman reference frame, `align_and_match()` recovering a known rotation + shift, cross-pair matching, the triplet quality gate, the module list in `R/load.R`. |
| `test-utils-transform.R`, `test-align-helpers.R`, `test-alignment-anchors.R`, `test-alignment-collapse.R` | Transform primitives, pose search, anchor selection, collapse guards. |
| `test-determinism.R` | Aligners reproducible run to run and leave the caller's RNG untouched. |
| `test-matching.R`, `test-ldir-gate.R` | Hungarian one-to-one matching, distance gates (50 µm; LDIR 20 µm), size-ratio limit, pair confidence. |
| `test-ldir-processed-image.R` | LDIR analyzed-overlay extraction and the image join. |
| `test-python-detector.R` → `inst/python/test_particle_detector.py` | Python LDIR detector on synthetic particle maps: every particle found at sub-pixel position, area, threshold auto-tuning, scan-circle mask, determinism, scan-circle centre. |

**Materials, agreement, validation**

| File | What it checks |
|---|---|
| `test-agreement-tiers.R`, `test-material-families.R`, `test-material-filter.R`, `test-reference-material-menu.R` | Exact / Family / Filler / Unresolved / Disagree tiers, quality and identification exclusions, polymer names and families. |
| `test-validation-tools.R` | Null model (few chance matches in a sparse field, many in a dense one; RNG restored), gate sensitivity, LDIR join check (correct join significant, scrambled not, single-polymer sample flagged uninformative) and its log line. |

**Viewer, report, provenance**

`test-ftir-image-placement.R`, `test-multirun-placement.R`, `test-plot-bounds.R`,
`test-view-orientation.R`, `test-image-scale-check.R`,
`test-raman-calibration-resolve.R`, `test-report.R`, `test-manifest.R`,
`test-reproducibility.R`, `test-repro-quality-gate.R`, `test-file-picker-cancel.R`,
`test-try-or.R`.

## Known gaps

- **Scan-circle radius.** `detect_scan_circle()` overestimates the radius by
  ~12% on a synthetic disk (recorded as an expected failure in
  `test_particle_detector.py`), and on dataset A's optical image (a crop with
  no visible circle) it returns a wrong circle that the R side rejects. Needs
  a real optical LDIR export with a visible scan circle before changing the
  detector. The pipeline uses the analyzed overlay when present.
- **Shiny viewer interaction** (no `shinytest2` test; the viewer's helper
  functions are tested individually).
- **One real dataset in the regression test** (dataset A). Further samples
  are listed in [docs/roadmap.md](../docs/roadmap.md).
