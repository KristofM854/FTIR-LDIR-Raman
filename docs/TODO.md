# To do before the release

Open items from the pre-release audit (branch `audit/pre-release`). Items are
grouped by what they are waiting for.

## Waiting for example datasets

The audit had only one FTIR export and synthetic data. With the full
datasets (FTIR + image, Raman + micrograph + WITec values, LDIR .xlsx + image,
Bruker if available):

- [ ] Before/after comparison of all outputs on every dataset (Phase 2 refactor
      must reproduce them or explain each difference).
- [ ] Null model: expected chance matches at the chosen gates, per instrument
      pair (method validated on synthetic data: ~28 chance matches of 130 at
      160/158 particles, 100 µm gate).
- [ ] Parameter sensitivity: match_dist_threshold_um, ransac_inlier_dist_um,
      landmark sizes, size prefilters, HQI and LDIR quality cutoffs.
- [ ] LDIR join check against Raman (same-family neighbour rate vs chance),
      replacing the scan-order test, which cannot detect a bad join.
- [x] Actual scales on dataset A (all four instruments): FTIR→Raman 0.9965,
      Bruker→Raman 1.0021, LDIR→Raman 1.0381; the ICP clamp never engaged, so
      the 0.8-1.25 band stays for every pair. Re-check on the other datasets.
- [ ] Material vocabulary: counts of names falling to "Unknown" / "Other
      polymer" / Unresolved per instrument (`tools/material_vocabulary.R`;
      dataset A: all identified labels resolve, PET only).
- [ ] Regression test on one real dataset (needs your OK to commit a small one).
- [ ] docs/validation.md with all numbers per dataset.

## Dataset A (all four instruments) — numbers so far

| Pair | gate | n src / Raman | matched | chance matches (shift / rotation / Poisson) | match distance median / p90 / max |
|---|---|---|---|---|---|
| FTIR (PerkinElmer)–Raman | 100 µm | 33 / 43 | 32 | 1.6 / 2.7 / 2.1 | 8.7 / 21.1 / 61.5 µm |
| FTIR (Bruker)–Raman | 100 µm | 42 / 43 | 29 | 1.9 / 2.8 / 2.6 | 15.8 / 41.3 / 56.8 µm |
| LDIR–Raman | 250 µm | 40 / 43 | 34 | 8.5 / 10.7 / 13.4 | 8.5 / 40.9 / 113.0 µm |

Matched pairs are identical before and after the Phase 1.5 changes (all six
pair files). The manual LDIR swap file corrects 10 of 40 image joins.

## Waiting for a decision

- [x] ICP / alignment scale band (audit Q6): keep 0.8-1.25 for every pair
      (dataset A, see above).
- [x] Spectral-library index (audit Q7): option (a), use it and disclose;
      internal use for now. The run records whether it was used.
- [x] "Acrylate" family (incl. polyacrylonitrile): Synthetic.
- [x] Triplet consensus applies the agreement quality gate (build_triplets();
      "Not scored" when fewer than two sides pass).

## Manual steps on your machine

- [ ] Re-snapshot `renv.lock` (it lacks plotly, reticulate, gridExtra,
      shinyFiles, zip, tiff, pdftools, testthat) and commit `renv/activate.R`:
      `renv::init(bare = TRUE); renv::snapshot(type = "explicit")`.
- [ ] Replace the versions in `inst/python/requirements.txt` with the ones you
      ran the published analysis with (`pip freeze`).
- [ ] Check that the GitHub workflow is green after the push (it now pins
      R 4.6.1 and installs everything DESCRIPTION declares; not run on GitHub
      from this branch yet).
- [ ] ORCID, co-authors and DOI in CITATION.cff / .zenodo.json (Phase 4).

## Deferred audit findings (Phase 2 or later)

- M2  Transform selection counts many-to-one inliers (RANSAC vs global race).
- M3  Landmark "confident" has no chance baseline (as few as 4 inliers).
- M4  Matching cost is effectively distance-only; document in the paper.
- M5  Composite (1:many) matches can reuse a Raman fragment.
- M8  Raman/LDIR image-placement fallbacks still derive extent from particles
      (display only).
- M10 Triplets: no consistency check against the direct LDIR↔FTIR match;
      `ldir_alignment_method` records only landmark/ransac.
- M11 Dataset-specific defaults (LDIR −90° pre-rotation, debug IDs A3/MP_11,
      anchor material vocabularies).
- [x] M12/M13 Duplicated PerkinElmer/Bruker pipeline; 2,074-line `main.R`:
      now `run_pipeline(config)` + `align_and_match()`; main.R is file input
      and config only.
- [x] M16 FTIR "size" is the Major dim, Raman's is Feret max: stated in
      docs/input_format.md.
- Python detector has no unit tests (Phase 3).

## Decided

- Git history is not rewritten. The only items in history are a 49 KB
  processed FTIR export (Feb–Mar 2026, your own measurement) and Windows
  user paths in old file versions; no licensed library files, credentials or
  large binaries. A rewrite would invalidate every existing clone, fork and
  PR reference for no real protection, since the Zenodo deposit is a snapshot
  of the release, not of the history.
