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
- [ ] Actual LDIR scale on your data (is it inside the ICP band 0.8-1.25?).
- [ ] Material vocabulary: counts of names falling to "Unknown" / "Other
      polymer" / Unresolved per instrument.
- [ ] Regression test on one real dataset (needs your OK to commit a small one).
- [ ] docs/validation.md with all numbers per dataset.

## Waiting for a decision

- [ ] ICP / alignment scale band (audit Q6) — see the explanation in the
      Phase 1.5 report. Currently 0.8-1.25 for every instrument pair, now a
      declared parameter.
- [ ] Spectral-library index for published numbers (audit Q7). The run now
      records whether it was used (manifest `spectral_library_index`).
- [ ] "Acrylate" family (incl. polyacrylonitrile) is categorised
      Semi-synthetic (`R/08b_material_map.R`, `semi_synthetic_families`).
- [ ] Triplet consensus does not yet apply the agreement quality gate
      (it already uses the Unknown/catch-all rule). Planned with Phase 2b.

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
- M12/M13 Duplicated PerkinElmer/Bruker pipeline; 2,074-line `main.R`
      (→ `run_pipeline(config)` + `align_and_match()`, Phase 2b).
- M16 FTIR "size" is the Major dim, Raman's is Feret max — state in the
      input contract.
- Python detector has no unit tests (Phase 3).

## Decided

- Git history is not rewritten. The only items in history are a 49 KB
  processed FTIR export (Feb–Mar 2026, your own measurement) and Windows
  user paths in old file versions; no licensed library files, credentials or
  large binaries. A rewrite would invalidate every existing clone, fork and
  PR reference for no real protection, since the Zenodo deposit is a snapshot
  of the release, not of the history.
