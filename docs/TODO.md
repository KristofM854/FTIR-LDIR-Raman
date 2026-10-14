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
      pair: `tools/validation.R` (`examples/validate_example_A.R`); dataset A
      below. Repeat on the other datasets.
- [ ] Parameter sensitivity (gate, landmark tier on/off, landmark size,
      Raman size filter, TPS, RANSAC inlier distance): `tools/validation.R`;
      dataset A below. Quality-cutoff sweep still to add.
- [x] LDIR join check against Raman: `ldir_join_check()` (R/ldir_join_check.R;
      family agreement and LDIR-table vs Raman size correlation, against
      randomly re-assigned LDIR rows). Logged by every run (WARN when no better
      than random), written to 07_diagnostics/ldir_join_check.csv and
      ldir_join_pairs.csv; replaces the scan-order test. Dataset A below.
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

Validation (`examples/validate_example_A.R`, 200 displacements per null):

| Pair | observed | chance: rotation / shift / Poisson | gate 0.5x → 2x: matched |
|---|---|---|---|
| FTIR (PerkinElmer)–Raman | 32 | 2.9 / 1.7 / 2.1 | 31 → 33 |
| FTIR (Bruker)–Raman | 29 | 3.1 / 2.0 / 2.6 | 27 → 31 |
| LDIR–Raman (150 µm) | 34 | 5.8 / 4.1 / 5.4 | 33 → 34 |

- Alignment: every sweep value (landmark tier on/off, landmark size 50–200 µm,
  Raman size filter 10–40 µm, TPS on/off, RANSAC inlier distance 100–300 µm
  with Tier 2 forced) keeps 100% of the baseline pairs. Different coarse
  starts (landmark vs RANSAC, rotations differing by up to 0.4°) converge to
  the same final transform after ICP.
- LDIR join: size Spearman LDIR table vs Raman 0.78 (random join 0.00, 95%
  0.28; p < 0.001; biased upward by the matcher's size term). Family
  agreement is uninformative here (single-polymer sample: random join 91%).

- [x] WITec image placement file (sample Faanui SW1L R4 B): X/Y = top-left
      corner, stage frame; 100% of 14,562 particles on the image, brightness
      peak at zero offset and scale 1.00. Used automatically when found.
      The current auto-calibration found no confident fit on this sample.

## Waiting for a decision

- [x] ICP / alignment scale band (audit Q6): keep 0.8-1.25 for every pair
      (dataset A, see above).
- [x] FTIR<->Raman keeps a fitted scale (not fixed to 1): dataset A scales
      0.9965 / 1.0021.
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
- [ ] CITATION.cff / .zenodo.json: ORCID, affiliation, co-authors,
      date-released, and the DOI once Zenodo has minted it (then the
      method paper as preferred-citation).
- [ ] docs/ai_assistance.md: confirm whether any commits under your name
      (e.g. 2026-02-27) contain code from other AI tools.

## Deferred audit findings

- [x] M2 Tier 2 race scored on one-to-one pairs (tier2_pairs()).
- [x] M3 Landmarks trusted only with >= 8 one-to-one landmark pairs
      (landmark_min_pairs); dataset A has 18 / 17 / 16.
- [ ] M4 Matching cost is effectively distance-only; state in the paper.
- [x] M5 Composite matches no longer share a Raman fragment.
- [x] M8 Image placements without metadata / calibration / registration are
      labelled "approximate image placement" in the viewer and report.
- [x] M10 Triplets: direct_ldir_ftir_match flags LDIR and FTIR particles that
      are not each other's direct match (dataset A: 0 of 32).
- [ ] M11 Anchor material lists (PET/PP) as default: open question.
      LDIR -90 deg pre-rotation kept as the Agilent export convention;
      debug trace IDs removed.
- [x] M12/M13 run_pipeline(config) + align_and_match().
- [x] M16 Size definitions stated in docs/input_format.md.
- [x] L3 unmatched_ftir / unmatched_raman aliases removed for non-FTIR pairs.
- [x] Python detector unit tests.

## Fixed during Phase 3

- Coarse RANSAC dropped scale 1 from its candidates when the clouds' span
  ratio came out just below 1 (a few particles outside the shared field,
  e.g. 0.962), and then aligned 4% off scale: synthetic test with 4 of 104
  FTIR particles 20 mm away matched 4 instead of 100 pairs. Scale 1 is now
  always a candidate (R/04_ransac.R). Dataset A unchanged.
- PNG files written by the run carried their write time (tIME chunk), so two
  identical runs differed byte-wise; now written without it.

## Decided

- Git history is not rewritten. The only items in history are a 49 KB
  processed FTIR export (Feb–Mar 2026, your own measurement) and Windows
  user paths in old file versions; no licensed library files, credentials or
  large binaries. A rewrite would invalidate every existing clone, fork and
  PR reference for no real protection, since the Zenodo deposit is a snapshot
  of the release, not of the history.
