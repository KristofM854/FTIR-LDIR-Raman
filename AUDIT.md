# Pre-release audit

Base: `audit/pre-release` at `5eaf53e`. Fixes: branch `audit/fable-fixes`.
Environment of the audit: Ubuntu 24.04, R 4.3.3 (Ubuntu packages; `renv.lock` is R 4.6.1),
Python 3.13 with numpy 2.5.3 / scipy 1.18.1 / Pillow 12.3.0, C locale, no display.
`rintrojs` could not be installed (no CRAN access); a stand-in exporting `introjsUI`,
`introjs`, `introBox` was used for the viewer check.

## 1. Summary

> **Update after the author's review.** H-01, H-02 and H-03 were implemented on
> request (commits `44f5149`, `c7bef98`, `56ba751`; they change results and are marked
> **(results)** in the CHANGELOG). The README now states the limits of the LDIR join
> without a `.sam` file (`d09a467`). Both open questions (paraffin, tin(II) chloride)
> stay as they are: unresolved. The sections below describe the state at the base
> commit; "Follow-up" at the end gives the numbers after these changes. H-04 is not
> done: the file is no longer in the tree, and removing it from the history means
> rewriting the history of `main`.

The pipeline is in good shape for an archive. The baseline test suite passes
(1468 passed, 0 failed, 2 skipped for a missing package). Dataset A runs end to end,
two runs give byte-identical tables, the transforms reproduce the aligned coordinates
from the raw files to 4e-5 µm, and the log, the CSVs, the report and the viewer show the
same pair counts and rates. On an R version three minor releases older than the lockfile
the regression fingerprint still matches exactly. No `$` partial-matching warning comes
from the repository's code in a full run, and no cached viewer output reads an input
that is missing from its cache key.

The alignments of dataset A are plausible and far above chance, and the matched pairs
are real. The weak points are not in the geometry but in three places where a number
or a label says more than the code supports.

The three most important problems (all open, for the author to decide):

1. **H-01. Exact-polymer scoring ignores bracketed qualifiers.** LDIR
   "Polyamide (naturally occurring)" is scored as an *Exact* match with PA, the opposite
   of what the README promises; "Polyethylene (chlorinated)" likewise with PE. Dataset A
   is not affected (all scored pairs are PET), field samples are.
2. **H-02. The LDIR alignment is never tested for being determined.** When the LDIR
   landmark tier fails, the Tier 2 pose is used without the chance test FTIR gets, the
   run records `determined = TRUE` with the failed landmark p-value, and the viewer
   prints "fixed by the data". Shown with 12 unrelated particles: as FTIR they are
   "NOT DETERMINED", as LDIR they are accepted.
3. **H-03. The pair confidence `p_true` is too pessimistic on real data.** On
   dataset A it implies 1.4 / 3.0 / 1.0 wrong pairs where the density model of the same
   run expects 0.6 / 0.5 / 0.1. The "chance-corrected" identified rate built on it drops
   FTIR (Bruker) → Raman from 75% to 62% while the chance-corrected detected rate on the
   same line stays at 79%.

Also: one real data file is in the git history (H-04, report only).

Six Medium and three Low findings were fixed on `audit/fable-fixes`; every table of a
dataset A run after the fixes is byte-identical to the baseline except the new `md5`
column of `supporting_files.csv`.

## 2. Dataset A assessment

Command: `Rscript examples/run_example_A.R <tmpdir>`, run twice.

**Alignment** (`04_alignment/alignment_summary.csv`)

| Instrument → Raman | Method | Scale | Rotation | Mirror | ICP RMS | Landmark pairs, p | σ fitted | Gate |
|---|---|---|---|---|---|---|---|---|
| FTIR (PerkinElmer) | landmark | 0.9965 | −178.69° | no | 14.2 µm | 18, 2.9e-25 | 8.0 µm | 50 µm fixed |
| FTIR (Bruker) | landmark | 1.0021 | 1.62° | no | 20.1 µm | 17, 7.4e-25 | 12.3 µm | 50 µm fixed |
| LDIR | landmark | 1.0384 | −178.45° | no | 14.1 µm | 16, 4.1e-22 | 7.1 µm | 21.2 µm (3σ) |

All three are plausible: scales within 0.4% of 1 for both FTIR, rotations of ~180° /
~0° as expected for a remounted filter, no mirror, landmark pairs 50 times above the
chance expectation (0.3). The LDIR scale of 1.038 is a real 4% scale error of the LDIR
coordinates, which come from an image whose µm/px is calibrated from the summed particle
areas; it is absorbed by the fit and the same value appears against FTIR (1.041 in an
FTIR-reference run). The 50 µm FTIR gates are 6σ and 4σ of the fitted errors: wide, but
with 38 reference particles the chance of an unrelated partner is 2%, so nothing is
gained by narrowing them here.

**Matches**

| Pair | Pairs | Distance median / p90 / max (µm) | Feret ratio median (range) | Size Spearman | Chance pairs expected |
|---|---|---|---|---|---|
| FTIR (PE) – Raman | 31 | 8.6 / 20.6 / 61.5 | 1.15 (0.90–2.05) | 0.98 | 0.6 |
| FTIR (Bruker) – Raman | 18 | 16.6 / 34.4 / 43.3 | 0.76 (0.32–1.68) | 0.74 | 0.5 |
| LDIR – Raman | 28 | 7.3 / 20.7 / 58.4 | 1.05 (0.89–2.60) | 0.97 | 0.1 |
| LDIR – FTIR (PE) | 31 | 6.6 / 11.1 / 21.9 | 0.93 (0.66–2.27) | 0.96 | 0.5 |
| FTIR (Bruker) – FTIR (PE) | 16 | 20.6 / 32.9 / 49.7 | 0.65 (0.30–0.83) | 0.87 | 0.4 |
| FTIR (Bruker) – LDIR | 14 | 16.9 / 27.7 / 42.1 | 0.74 (0.33–0.95) | 0.77 | 0.3 |

The pairs are real: sizes of partners correlate strongly, all scored pairs are PET–PET,
all 28 triplets have the LDIR and FTIR particle as each other's direct match, and the
LDIR–FTIR (PE) pairs sit within 22 µm although neither is the reference. The Bruker
table is the outlier: it sizes particles about a third smaller than the other three
instruments and has the largest positional error.

**Rates and partner status agree with each other.** Checked for every direction:
`n_found` equals the matched count, `n_base` equals the partner-status rows that pass
and are not `outside_area`, and `n_detected − n_found` equals the passing rows with
`partner_below_cutoff` or `paired_elsewhere`. Example, Raman particles found by
FTIR (Bruker): 38 particles ≥ 20 µm = 18 matched + 10 not found + 3 outside + 7 partner
below cutoff; base 35; identified 18/35 = 51%; detected 25/35 = 71%. The asymmetries fit
the instruments: the Bruker list has 16 of 41 particles below AAU 0.70, so Bruker
*identifies* about half of what the others found but *detects* about 70%; PerkinElmer and
LDIR see nearly the same particles (97–100%).

**Particles traced by hand** (raw table → normalised → transform file → aligned CSV →
match and partner-status files). FTIR (PE) aligned coordinates recomputed from the raw
CSV and `transform_params_ftir_perkin_raman.txt`: maximum deviation 4e-5 µm. Raman:
raw minus the centroid of all 43 particles, 2e-5 µm.

Pairs:
- FTIR (PE) `MP_25` (448 µm) ↔ Raman `42` (218 µm) at 61.5 µm, inside the widened gate
  (0.15 × 448 = 67 µm). LDIR `A2` (428 µm) sits 6 µm from `MP_25` and pairs with the same
  Raman particle at 58 µm. Three instruments agree on the place; Raman reports half the
  length. Real, a Raman fragment or a different outline.
- FTIR (PE) `MP_8` (251 µm) ↔ Raman `97` (217 µm) at 48.2 µm, just inside 50 µm. LDIR
  `A13` (223 µm) is 2.1 µm from `MP_8` but 46 µm from Raman `97`, outside its own gate
  (33 µm): the same particle is matched FTIR–Raman and unmatched LDIR–Raman.
- FTIR (Bruker) `MP_20` (38 µm) ↔ Raman `8` (119 µm) at 43 µm. FTIR (PE) `MP_15`
  (127 µm) and LDIR `A22` (117 µm) are at the Raman position; Bruker reports a third of
  the size, 28–31 µm off. Same particle, sized and centred differently by Bruker.
- LDIR `A22` ↔ Raman `8` at 17.5 µm (gate 21.2), ↔ `MP_15` at 2.9 µm. Clean.

Unmatched:
- FTIR (PE) `MP_10` (374 µm) and LDIR `A1` (349 µm) are 3.2 µm apart and both unmatched
  to Raman. Raman `1` (341 µm) is 112 µm away and is the same particle by size; Bruker
  lists it as two particles (`MP_26`, `MP_27`), one at each end.
- FTIR (PE) `MP_1` (321 µm) and LDIR `A4` (318 µm), 5.8 µm apart, unmatched to Raman;
  Raman `80` (300 µm) is 81 µm away.
- Raman `4` (25.6 µm): nothing of any instrument within 120 µm. A small particle only
  Raman found.
- Raman `90` (205 µm): FTIR (Bruker) `MP_7` at 9.9 µm with AAU 0.687, correctly
  `partner_below_cutoff`.
- FTIR (Bruker) `MP_30`: status `paired_elsewhere` with Raman `44`, which is 7.3 µm in
  size and paired with nothing (see M-08).

**Suspicious or worth stating in the paper**

- *Large particles next to each other but unmatched* (M-12). For particles above
  ~200 µm the centres reported by Raman differ from FTIR/LDIR by 50–110 µm. The
  \> 200 µm class therefore has the lowest rates of all classes (64–85% against Raman,
  100% between FTIR (PE) and LDIR). That dip is a centroid effect, not missed detection.
- *The LDIR result rests on a manual file* (M-11). `LDIR_coord_swaps.csv` moves 16 of
  the 40 LDIR positions. Without it the run gives 29 instead of 28 LDIR–Raman pairs and
  the join check reads the same (size Spearman 0.97 with, 0.97 without): the check
  cannot see the permutations among similar-sized particles that the file corrects.
- *`p_true`* marks the six largest-distance pairs (all large particles, all with a
  third instrument confirming them) as probably wrong (H-03).
- The LDIR gate is 21.2 µm in the run and "≤20 µm" in the viewer (fixed, M-03).
- Nothing else: no count differs between log, CSVs and report.

## 3. Findings

Status "fixed" gives the commit on `audit/fable-fixes`.

| ID | Sev. | Area | File:line (base) | Summary | Evidence | Status |
|---|---|---|---|---|---|---|
| H-01 | High | Material | `R/08b_material_map.R:682-693`, `:772-773` | A trailing bracket is dropped before the exact-polymer comparison, and Exact overrides the family: "Polyamide (naturally occurring)" = PA, "Polyethylene (chlorinated)" = PE | CONFIRMED | fixed `44f5149` |
| H-02 | High | Alignment | `R/run_pipeline.R:643-644`, `R/ldir_stage.R:346-429` | LDIR Tier 2 pose has no chance test; recorded as determined with the failed landmark p | CONFIRMED | fixed `c7bef98` |
| H-03 | High | Statistics | `R/07_match.R:120-131`, `:206-229`; `R/detection_rates.R:110-118` | `p_true` and the "chance-corrected" identified rate count heavy-tailed true pairs as chance | CONFIRMED | fixed `56ba751` |
| H-04 | High | Hygiene | git history `9dec7f6` | `processed_data_ftir_20260211_102404.xlsx` (150 particles, samples UW1M1–5) added and later deleted | CONFIRMED | open (report only) |
| M-01 | Medium | Statistics | `R/detection_rates.R:105` | Run stops when `min_particle_size_um` < 20 and a particle is below 20 µm | CONFIRMED | fixed `4d920a2` |
| M-02 | Medium | Matching | `R/07_match.R:276-277` | Run stops ("no points in query!") when all particles of an instrument are below its cutoff | CONFIRMED | fixed `c952773`, `ef36576` |
| M-03 | Medium | Viewer | `shiny_app/global.R:542-556` | Viewer re-gates LDIR–Raman pairs with the configured 20 µm, not the gate the run used | CONFIRMED (function level; no pair lost on dataset A) | fixed `df209bd` |
| M-04 | Medium | Viewer | `shiny_app/global.R:3978` | Positional error missing from the alignment check for runs without Raman | CONFIRMED | fixed `24a964d` |
| M-05 | Medium | Provenance | `R/utils.R:1588-1601` | Supporting files that change results recorded without a hash | CONFIRMED | fixed `b32317f` |
| M-06 | Medium | Report | `R/10_export.R:105` | `all_diagnostics.pdf`: dashes and µ become dots outside a UTF-8 locale | CONFIRMED | fixed `e040384` |
| M-07 | Medium | Alignment | `R/ldir_stage.R:327`, `:478` | LDIR landmark tier and ICP accept a mirrored pose although `ldir_allow_reflection = FALSE` | CONFIRMED | open (comment corrected, L-03) |
| M-08 | Medium | Statistics | `R/detection_rates.R:192` | `paired_elsewhere` also labels partners that are paired with nothing (below the size limit, size ratio > 5) | CONFIRMED | open |
| M-09 | Medium | Material | `R/08b_material_map.R:145-149`, `:233` | Every "Cellulos…" label is Semi-synthetic and so in the "plastics" group, including natural cellulose | CONFIRMED (classification); judgement | open |
| M-10 | Medium | Statistics | `R/03b_landmark_align.R:205`, `R/align_helpers.R:225`, `shiny_app/global.R:4007` | `pose_p` holds two different quantities under one name and one verdict threshold | CONFIRMED | open |
| M-11 | Medium | LDIR join | `R/ldir_join_check.R`, `examples/dataset_A/LDIR_coord_swaps.csv` | Join check blind to swaps among similar-sized particles; dataset A needs 10 manual rows | CONFIRMED | documented in the README `d09a467`; check itself open |
| M-12 | Medium | Matching | `R/07_match.R:268-273` | Centroid-based gate fails for particles > 200 µm; lowers the rates of that class | CONFIRMED | open |
| M-13 | Medium | Viewer | `shiny_app/app.R:5063-5065` | Overlay "unmatched" counts use all particles; the pipeline's use the matchable ones (12 vs 7 Raman, 24 vs 7 Bruker) | CONFIRMED | open |
| M-14 | Medium | Provenance | `R/utils.R:925` | Python package versions not recorded in the run | CONFIRMED | open |
| M-15 | Medium | Hygiene | `R/00_config.R:163-167`; `tests/testthat/test-ldir-gate.R:177` | Default config carries the OPUS values of a field sample; field sample names and a personal folder path in committed files | CONFIRMED | open (report only) |
| M-16 | Medium | Tests | `tools/regression_fingerprint.R` | Fingerprint does not cover detection rates, partner status, chance pairs, `p_true`, LDIR join positions | CONFIRMED | open |
| L-01 | Low | Tests | `DESCRIPTION:40-44` | `writexl` undeclared: two reader tests skip in CI | CONFIRMED | fixed `eba49d8` |
| L-02 | Low | Docs | `README.md` | Tab list, output table and module list behind the code | CONFIRMED | fixed `389ec3a`, `3c3c218` |
| L-03 | Low | Docs | `R/00_config.R:546`, `:593`; `README.md` | Comments on the LDIR gate and the mirror search contradict the code | CONFIRMED | fixed `81b8a13` |
| L-04 | Low | Code | `shiny_app/app.R:805`, `:1050` | `uploaded_data` is always NULL (18 references); stale "Tab 7: Data Upload" comment | CONFIRMED | open |
| L-05 | Low | Provenance | `04_alignment/transform_params_ldir_raman.txt` | LDIR transform file has no matrix, no normalisation, no −90° pre-rotation; method names differ ("landmark" vs "landmark_ransac") | CONFIRMED | open |
| L-06 | Low | Tests | several | Tests that cut code out of `main.R` / `app.R` by text markers skip when the marker moves; `.sam` tests use a writer built from the reader's own spec | CONFIRMED | open |
| L-07 | Low | Readers | `R/01_ingest.R:76-101`, `:312` | Column lookup by substring: "Area" takes the first column containing "area" | PLAUSIBLE | open |
| L-08 | Low | Statistics | `R/02_prefilter.R:78`, `R/detection_rates.R:67` | A particle without a size is matched but not counted in the rate base | PLAUSIBLE | open |
| L-09 | Low | Provenance | manifest | Absolute input paths and the user name in `manifest.json`; `manifest.json.bak` left behind | CONFIRMED | open |
| L-10 | Low | Archive | `CITATION.cff`, `.zenodo.json` | No release date, ORCID or DOI; Zenodo licence lists MIT only although the data is CC BY 4.0 | CONFIRMED | open |
| L-11 | Low | Logs | several | Bruker pair logged as "ftir -> raman"; "Raman landmarks" for any reference; viewer says "placed from WITec extent" for an auto-calibrated image; a literal BOM character in `R/01c_ingest_ldir.R:2147` | CONFIRMED | open |

## 4. Details per finding

### H-01 Exact-polymer scoring ignores bracketed qualifiers

`canonical_polymer_vec()` tries the whole name, then the name *without its trailing
bracket*. That is right for an abbreviation ("Polycarbonate (PC)") and wrong for a
qualifier. `score_tiered_agreement_vec()` then sets `result[exact] <- "Exact"` last, over
whatever the families said.

```r
a <- c("Polyamide (naturally occurring)", "Polyethylene (chlorinated)", "Paraffin (PE)")
b <- c("Polyamide (PA)", "PE", "PE")
data.frame(pol_a = canonical_polymer_vec(a), fam_a = classify_family_vec(a),
           fam_b = classify_family_vec(b), tier = score_tiered_agreement_vec(a, b))
#>  pol_a         fam_a fam_b  tier
#>     PA       Natural    PA Exact
#>     PE Other polymer    PE Exact
#>     PE       Unknown    PE Exact
```

The README states that "Polyamide (naturally occurring)" is kept out of PA "to avoid
false material-agreement matches". A related case goes the other way: the same LDIR
label against Raman "Protein", "Wool" or "Silk" is *Disagree* (family Natural vs
Protein), although natural polyamide is protein.

Affects: `tier`, `agree_exact`, `agree_family` in `agreement_pairwise_*.csv`, every
family and exact rate, `material_polymer` in the particle tables, triplet consensus.
Dataset A: none (scored pairs are PET–PET).

Proposed fix: (1) in `canonical_polymer_vec()` use the base name only when the bracket
is an abbreviation (maps to the same canonical polymer, or is not a word such as
"naturally occurring", "chlorinated", "natural"); (2) in `score_tiered_agreement_vec()`
require the families to be equal for Exact; (3) decide whether LDIR natural polyamide
belongs to Protein. Risk: changes agreement numbers on any run with these labels, which
is the point; needs the vocabulary of the field samples to be re-checked with
`tools/material_vocabulary.R`.

### H-02 LDIR alignment is never tested for being determined

`align_and_match()` tests every Tier 2 pose with `pose_significance()` and withholds
scoring when p ≥ 0.01. `run_ldir_stage()` has no such test, and `run_pipeline()` writes
`determined = TRUE, pose_p = <landmark p>` for LDIR whatever tier was used.

Shown with a synthetic dense reference (12 000 particles) and 12 particles with no
counterpart, given once as LDIR and once as FTIR (Bruker):

```
 instrument  method          scale  icp_rms_um  pose_p   determined  n_matched
 ftir_bruker global_register 1.03   72.3        1.0e+00  FALSE       3
 ldir        ransac_icp      0.80   65.7        3.6e-05  TRUE        0
```

The LDIR scale ran into the 0.8 clamp, `pose_p` is the p of the *rejected* landmark
attempt ("not confident: p = 3.6e-05; 4 < 5 pairs"), and the viewer's rule
(`pose_p < 1e-3` → "fixed by the data") approves it. Here the 20 µm gate left no pair;
with a longer unrelated list, or a wrong LDIR image scale, chance pairs would be scored
and LDIR ↔ FTIR pairs and triplets built on them.

Proposed fix: after the LDIR ICP, when the landmark tier was not used, compute
`pose_significance()` on the LDIR anchors as in `align_and_match()`; store its p and
`determined`; when not determined, score nothing and skip LDIR ↔ FTIR and triplets.
Risk: low for good data (dataset A uses landmarks, p = 4e-22); a field run that
currently passes through Tier 2 could be reported as not determined, which needs a look
at those runs before release.

### H-03 `p_true` and the "chance-corrected" identified rate

`.position_error()` fits a two-component mixture to the nearest-neighbour distances
within 75 µm: a Rayleigh (true partners, 2-D normal error) and the chance
nearest-neighbour law at the reference density, with a free mixing weight w. True
errors are not normal for large particles (section 2), so the fit explains the tail by
lowering w, that is by calling those pairs chance:

```
              sigma   w      chance neighbours in the window
                             fitted (1-w)·n   expected from the density
ftir_perkin   8.0     0.93   2.3              0.10
ftir_bruker   12.3    0.71   6.0              0.19
ldir          7.1     0.88   3.5              0.14
```

The density is known (2.6e-6 per µm², chance of an unrelated neighbour within 75 µm:
4.6%), so the fitted chance share is 20–30 times too high. Consequences on dataset A:

- `sum(1 − p_true)` = 1.4 / 3.0 / 1.0 "wrong" pairs against `chance_pairs_expected`
  0.6 / 0.5 / 0.1 in `chance_pairs.csv` of the same run;
- `rate_true_pct`: FTIR (Bruker) → Raman 75.0% → 62.4%, FTIR (PE) → Raman 93.9% → 89.6%,
  while `detected_corrected_pct` in the same row moves by 0.1–0.5 points;
- "pairs with p_true ≥ 0.9" drops `MP_25`–`42`, `A2`–`42`, `MP_8`–`97` and the Bruker
  pairs traced in section 2, all confirmed by a third instrument.

The synthetic field test passes because its errors are drawn from a normal distribution.

Proposed fix, in order of preference: (a) do not fit w; set the chance share from the
density, `1 − w = n_unpaired · P(chance within W) / n`, and fit σ only (robustly, e.g.
from the median); (b) keep the fit but describe `p_true` as a distance-based confidence
and stop calling `rate_true_pct` chance-corrected; (c) drop `rate_true_pct` and keep the
density-based correction used for the detected rate. Risk: (a) changes `p_true`, the
high-confidence agreement and `rate_true_pct` in every run; it should be re-validated on
the synthetic field sample, where both approaches must agree.

### H-04 A real data file is in the git history

```
$ git log --all --format="%h %ad %s" --date=short -- processed_data_ftir_20260211_102404.xlsx
c569646 2026-03-06 Delete processed_data_ftir_20260211_102404.xlsx
9dec7f6 2026-02-13 Add files via upload
```

49 kB, sheet `Long_Table` with 150 FTIR particles from source files `UW1M1.csv` …
`UW1M5.csv` plus count summaries. It is reachable from `main`. Apart from this, the full
history (730 commits) only ever added `examples/dataset_A/*`, `inst/materials/*.csv` and
the fingerprint. A Zenodo archive of a release contains the tagged tree, not the
history, so the deposit itself would be clean; the public GitHub history is not.
Nothing was rewritten. Options for the author: accept it, or publish the release from a
fresh repository / squashed history.

### M-01 Detection rates stop on particles below the first size class (fixed)

`cut()` gives NA below 20 µm, `sum(s)` became NA and `if (n > 0)` failed. Reproduced
with `config$min_particle_size_um <- 0` on dataset A: "Error in if (n > 0) ... missing
value where TRUE/FALSE needed" after all other outputs, before the report. Fix: such
particles count in "all sizes" only. Risk: none at the default limit (all tables
identical). Test in `test-align-helpers.R`.

### M-02 Run stops when one instrument has nothing to match (fixed)

Reproduced with `config$ftir_agreement_min_quality <- 0.999` on dataset A:
`Error in RANN::nn2(...): no points in query!`. Fix: `.match_hungarian()` and
`.match_greedy()` return an empty result for an empty side; the run completes with
0 pairs for that instrument (`chance_pairs.csv` shows `n_matched 0`). Risk: low; only
the empty case takes the new path. Test in `test-matching.R`. The first commit
(`c952773`) carried a broken test fixture, corrected in `ef36576`.

### M-03 Viewer re-gates LDIR pairs with the configured gate (fixed)

`ldir_acceptance_gate()` read `match_dist_threshold_ldir_um` (20) from the manifest;
the run's gate is fitted (dataset A 21.2 µm, up to 50 µm). `annotate_ldir_gate()` with
a pair at 20.8 µm and a 30 µm particle returned `within_gate = FALSE`. On dataset A no
pair falls in between (closest: `A22` at 17.5 µm), so the visible effect was the label
"28/40 LDIR↔Raman matched (≤20µm)". Fix: read `gate_um` of `ldir_raman` from
`chance_pairs.csv`, then the config, then 20. Risk: low, viewer and report only.

### M-04 Positional error missing for runs without Raman (fixed)

`paste0(instrument, "_", tolower(reference))` gives "ldir_ftir (perkinelmer)"; the file
has "ldir_ftir_perkin". Fix: map the display name to the file name. Risk: none.

### M-05 Supporting files recorded without a hash (fixed)

`supporting_files.csv` gets an `md5` column. Risk: none (additive; the viewer selects
columns by name).

### M-06 Non-ASCII titles in `all_diagnostics.pdf` (fixed)

`pdftotext` of the baseline file: "FTIR...Raman Particle Overlay"; about 50 "conversion
failure" warnings per run in a C locale. Fix: `cairo_pdf` when available, as the report
uses; device closed in a `finally`. After the fix: "FTIR–Raman Particle Overlay", 11
pages, tables identical. `particle_report.pdf` was already clean: its tables are ASCII
(`.report_ascii()`), figures carry "—", "–" and "µ", which `cairo_pdf` renders.

### M-07 LDIR accepts a mirrored pose in the landmark tier

`config$ldir_flip_y_for_alignment <- TRUE` on dataset A (a mirrored LDIR cloud):

```
Coarse search best: angle = 181 deg, mirror = TRUE ...
"ldir","Raman","landmark_ransac",1.0384,-178.45,TRUE,14.1,4.1e-22,TRUE,28
```

The reflection ban is applied to Tier 2 only (`.ldir_align_cfg`); `landmark_align()` and
`icp_refine()` receive the unmodified config with `ransac_allow_mirror = TRUE`. The
result is visible ("Mirrored: yes" in the alignment check), and a spurious mirrored pose
would still have to pass p < 1e-8. Proposed fix: pass `.ldir_align_cfg` to the landmark
tier and ICP. Not applied: it changes behaviour for any run that relies on it, and
whether a mirrored LDIR export should be accepted or refused is the author's call.
Dataset A would not change (the unmirrored pose wins).

### M-08 `paired_elsewhere` for partners that are paired with nothing

The last branch of the status rule is "both pass" and is named `paired_elsewhere`.
FTIR (Bruker) `MP_30` → Raman `44`: the Raman particle is 7.3 µm, below the 20 µm limit,
and appears in no match file. The same branch catches partners rejected by the 5× size
ratio. Proposed fix: test whether the partner is in the pair's match table; otherwise
use a new status ("partner not matchable: below the size limit or size ratio"). Needs a
new category in the viewer's legend, hence not applied. It also means "detected" counts
partners of any size, which is what the documentation says.

### M-09 Natural cellulose counts as semi-synthetic

Rule `"Cellulos"` puts "Cellulose", "Cellulosic" (the LDIR class for natural cellulose)
and cellulose acetate into family Cellulose, category Semi-synthetic.
`detection_rates()` builds the "plastics" group from synthetic + semi-synthetic, so on a
field sample the plastics rows include paper and cotton fibres labelled cellulose, and
"Cellulose" vs "Cotton" is a Disagree. Proposed: plain "Cellulose"/"Cellulosic" →
Natural; only regenerated or derivatised names (rayon, viscose, lyocell, acetate, CAB) →
Semi-synthetic. Changes results on field samples; a scientific decision.

### M-10 One column, two p-values

Landmark poses store the raw Poisson tail (accepted below 1e-8, to allow for ~40 000
hypotheses); Tier 2 poses store a Bonferroni-corrected p (accepted below 0.01). Both go
to `pose_p`, shown as "Pose by chance (p)" with one rule (< 1e-3 "fixed by the data").
A Tier 2 pose accepted at p = 0.005 reads "check"; a landmark p of 1e-5 that the
pipeline *rejected* reads "fixed by the data" (H-02). Proposed: store the corrected p
for both, or add a column naming the test, and derive the verdict from `determined`.

### M-11 The join check does not see what the manual swaps correct

See section 2. The check compares size ranks of LDIR rows and their Raman partners;
swapping two particles of similar size changes neither. Without a `.sam` file the size /
shape / colour join mis-places 16 of 40 rows on the dataset that ships as the example.
Proposed: say so in the README and the paper (the join is exact only with the `.sam`
file; otherwise it needs manual checking), and report in the log how many rows a swap
file moved as a share of the table.

### M-12 Centroid matching and large particles

The gate widens with the *source* particle's major axis only (0.15×). `MP_10`/`A1` vs
Raman `1` (112 µm apart, all ~350 µm long) and `MP_1`/`A4` vs Raman `80` (81 µm) stay
unmatched; `A13` is matched to FTIR but not to Raman. Options: widen by the larger of
the two major axes; or accept a pair when one centroid lies inside the other particle's
extent; or state the limitation and report the > 200 µm class separately. Any change
alters the matches of dataset A.

### M-13 "Unmatched" in the overlay summary

Headless viewer on the dataset A run: "2 FTIR unmatched | 12 Raman unmatched | 24 Bruker
unmatched". The log: 7 unmatched Raman (of 38 matchable), 7 unmatched Bruker (of 25).
The viewer counts every particle without a pair, including those below the size limit
or quality cutoff that the pipeline never tries to match. Proposed: count over the
matchable particles, or name the line "without a pair (any quality)".

### M-14 Python package versions

`session_info.txt` records the interpreter path only. The LDIR positions come from
numpy / scipy / Pillow. Proposed: write their versions (one `reticulate` call in
`setup_python_detector()`) into `session_info.txt`.

### M-15 Field-sample values and names in committed files

`make_config()` ships `opus_ftir_bruker_image` "Pre-filled: Faanui SW1L R4": per-scan
calibration of a field sample as a default for every run (it is checked against the
image and ignored with a warning when it does not fit). The sample name appears 23
times (config, comments, docs, tests), and `test-ldir-gate.R:177` contains a path with
an organisation's OneDrive folder and a field campaign name. None of this is data.
Proposed: default `opus_ftir_bruker_image = NULL` with the values as a commented
example; neutral names in tests and comments.

### M-16 What the regression fingerprint covers

Row counts of the ingested / prefiltered / match tables, every pair with its distance,
scale / rotation / translation, agreement tier counts. Not covered: `detection_rates.csv`,
`partner_status.csv`, `chance_pairs.csv` (gate, σ), `p_true`, triplet consensus,
`alignment_summary.csv`, the LDIR joined positions. A change to any statistic passes the
regression. It is stable (identical on R 4.3.3 and across two runs).

### Low findings (open)

- **L-04** `uploaded_data <- reactiveVal(NULL)` is never set to anything else; 18
  references, among them cache keys and report captions for "uploaded data". Removing
  it touches all of them in a 6 400-line file; proposed as its own change with a viewer
  check.
- **L-05** `transform_params_ldir_raman.txt` has scale, rotation and translation only;
  the LDIR centroid and the −90° pre-rotation are in `debug/ldir_raman_transform_audit.json`.
  After a TPS refinement the aligned coordinates no longer follow from the transform.
- **L-06** `test-file-picker-cancel.R`, `test-plot-bounds.R`, `test-image-scale-check.R`
  and `test-reference-material-menu.R` extract code by text markers and skip when a
  marker is not found. One test warns "no DISPLAY variable so Tk is not available" on a
  headless machine and falls back to the console prompt.
- **L-07** `find_column()` falls back to substring matching; "Area" resolves to the
  first column whose name contains "area". Correct for the WITec column order, silent
  if the order changes.
- **L-08** The prefilters keep particles with a missing size; the rate base drops them.
- **L-09** `manifest.json` stores absolute input paths and `Sys.info()` user. Worth
  knowing before run folders are shared.
- **L-10**, **L-11** as in the table.

### Open questions from the brief

- **Raman paraffin → PE family?** Today "Paraffin" is family Unknown (with the library
  index: Other compound), so a paraffin–PE pair is *Unresolved*: counted, flagged,
  excluded from the rates. That is the defensible default. Paraffin and PE spectra are
  close, and mapping paraffin to PE would turn a library ambiguity into agreement.
  Recommendation: keep it unresolved and report the number of paraffin–PE pairs
  separately. Note "Paraffin (PE)" currently scores Exact with PE through H-01.
- **Tin(II) chloride → unresolved?** Yes, already: family Unknown, tier Unresolved
  against any other name. Two instruments reporting the identical unknown name score
  Exact; that cannot happen between different libraries in practice.

## 5. Statistical review

**Wilson intervals.** Implemented correctly, vectorised, NA for n = 0, inside 0–1 for
0% and 100% (checked: 0/1 → 0–79%, 5/5 → 57–100%). They are binomial intervals for the
particles of one filter. They do not include the uncertainty of the alignment, of the
gate, or between filters, and the particles of one filter are not independent draws of a
population. For a methods section: "95% Wilson score intervals for the proportion within
one measurement".

**Denominators.** A rate of "A's particles found by B" is over A's particles that pass
A's quality cutoff, are at least 20 µm, and lie inside the convex hull of *all* of B's
particles (any quality, any size), plus any matched particle outside that hull.
*Identified* = paired one-to-one with a particle of B that passes B's cutoff and is at
least 20 µm. *Detected* = any particle of B within the gate, or paired. The convex hull
of B's particles is a proxy for B's scan area: it is smaller than the true area when B
found little near an edge, which removes edge particles of A from the base and raises
the rate slightly. The code, the log, the CSV, the report captions and the viewer
describe these definitions consistently.

**Chance pairs.** For each pair of instruments: reference density ρ = n / hull area;
chance of an unrelated reference particle within the gate r is p = 1 − exp(−ρπr²)
(Poisson, uniform); expected chance pairs = p × source particles inside the hull. This
is the number of pairs a random arrangement would produce, which is how the report words
it. It is an upper bound on wrong pairs in so far as a particle with a true partner
rarely takes a chance one, and a lower bound in so far as the widened gate of large
particles is ignored (the log says so). Uniformity is the main assumption; particles on
real filters cluster, which raises the local chance rate.

**Chance-corrected detected rate.** (k − e)/(n − e) with e the summed per-particle chance
probability: the standard correction, correctly applied and clamped to 0–100.

**Chance-corrected identified rate and `p_true`.** Not sound on real data as it stands
(H-03). On dataset A it removes 4–13 percentage points where the density model of the
same run supports about 1. Until it is changed, a paper should quote the raw identified
rate with its Wilson interval and the expected chance pairs from `chance_pairs.csv`, and
not `rate_true_pct` or "agreement over p_true ≥ 0.9" as chance-free values.

**Positional error σ.** A useful per-run number, and the 3σ LDIR gate works on dataset A
(28 pairs, 0.1 expected by chance). Two caveats: σ is fitted to the same nearest-neighbour
distances the gate is then applied to, and with a free mixing weight it describes the core
of the distribution only. The statement "99% of true pairs lie within 3σ" in the log
assumes a normal error and is contradicted by the run's own pairs (2 of 31 FTIR (PerkinElmer)
pairs beyond 3σ = 24 µm, 6% where a normal error gives 1%).

**Pose significance.** Tier 2: a Poisson tail on the anchor pairs beyond the two that fix
a pose, Bonferroni-corrected by (n_src · n_ref)² × 2 hypotheses. Conservative and
reasonable. Landmarks: an uncorrected Poisson tail with a threshold of 1e-8. Both are
defensible; they must not share a column (M-10), and LDIR needs the Tier 2 test (H-02).

**Agreement rates.** Family agreement = (Exact + Family) / resolved pairs, where a pair
enters only if both identifications pass their cutoffs and neither is "unknown", and
unresolved pairs are counted but excluded. Clearly defined and consistently reported.
The exclusion of unresolved pairs raises the rate when one instrument often returns a
reference chemical; the unresolved share is printed next to it, and should be quoted with
it. H-01 affects which pairs count as Exact.

**Detection by size class.** Classes [20, 50), [50, 100), [100, 200), ≥ 200 µm of the
*base* instrument's Feret maximum. With 3–15 particles per class on dataset A the
intervals span 40–60 points; the classes support no size trend on a single filter. The
≥ 200 µm class is biased low by M-12.

## 6. Test gaps

Added with the fixes: size classes below 20 µm (`test-align-helpers.R`), empty matcher
input for both methods (`test-matching.R`), recorded LDIR gate (`test-ldir-gate.R`),
alignment table without Raman (`test-report.R`), supporting-file MD5 (`test-manifest.R`).

Worth adding, in priority order:

1. Material tiers for bracketed qualifiers: "Polyamide (naturally occurring)" vs PA must
   not be Exact; "Polyethylene (chlorinated)" vs PE; "Polycarbonate (PC)" vs PC must stay
   Exact (H-01).
2. An LDIR table with no counterpart in a dense reference must be reported as not
   determined (H-02); the generator of `test-pipeline-field.R` can be reused.
3. `p_true` calibration with heavy-tailed errors (a share of pairs at 4–6σ), not only
   normal errors (H-03).
4. `partner_status()`: one test per status, including a partner below the size limit
   (M-08), n = 1 base particle, and an `other` with fewer than 3 particles.
5. `detection_rates()` against a hand-built case with known hull membership, including a
   matched particle outside the hull and a zero base.
6. Fingerprint extended to `detection_rates.csv`, `partner_status.csv` status counts,
   `chance_pairs.csv` and `02_joined/ldir_joined_raw.csv` positions (M-16).
7. A mirrored LDIR cloud, asserting whichever behaviour is decided (M-07).
8. `pose_significance()` and `chance_partner_rate()` edge cases: fewer than 3 reference
   points, collinear points (zero hull area), n_pairs ≤ 2.
9. A real `.sam` fragment as a fixture, if one can be published; today the decoder is
   tested only against bytes written from its own description.
10. A full run with `options(warnPartialMatchDollar = TRUE)` turned into a test that
    fails on a warning from repository code (currently clean).

## 7. Archive checklist

Before the Zenodo release:

- Decide H-01, H-02 and H-03; they change numbers a paper would quote.
- Decide what to do about the data file in the history (H-04).
- Remove the field-sample OPUS values from the default config and the personal path
  from the test (M-15).
- State in the README and the paper that the LDIR join without a `.sam` file needed 10
  manual corrections on dataset A (M-11) and that large particles can stay unmatched
  (M-12).
- `CITATION.cff`: add `date-released`, ORCID and, after the first deposit, the DOI;
  `.zenodo.json`: mention the CC BY 4.0 licence of the example data.
- Regenerate `renv.lock` if `writexl` should be available under `renv` (L-01).
- Tag the release from the commit whose hash the example run's manifest will show.

Can wait:

- M-08, M-10, M-13, M-14, M-16 and all Low findings.
- Removing `uploaded_data` (L-04).
- A decision on cellulose (M-09) can wait only if the paper reports no "plastics" rates
  from field samples.

## 8. What could not be checked, and why

- **The `.sam` decoder on a real file.** The repository has none (dataset A has no
  `.sam`); the reader was checked against synthetic records and for robustness to
  truncated and stray records only. The frame conversion (0-based column, y from the
  bottom → 1-based column / row from the top) is consistent with the code that consumes
  it, but unverified against Clarity output.
- **Field samples, the PC reference filter, segmented PerkinElmer exports, WITec
  placement files, OPUS placement, optical LDIR images without the analyzed overlay.**
  No real files; covered by synthetic tests only. Statements about field behaviour above
  are from code and synthetic data.
- **The spectral-library index.** The PDFs are licensed and absent; material results
  above are from the name rules alone.
- **Windows, file dialogs, Tk.** Linux, headless.
- **The viewer in a browser.** Checked headlessly (`shiny::testServer`, run loading,
  overlay summary, the Summary builders in `global.R`) and by static analysis of the
  cache keys; no rendering, no clicks, and the guided tour was replaced by a stand-in.
- **`renv::restore()` and R 4.6.1.** No CRAN access; R 4.3.3 from the distribution.
- **`docs/validation.md` numbers** other than those the example run reproduces, and
  `examples/validate_example_A.R`, were not re-run.
- **The two Excel reader tests** (`writexl` not installable here).
- **Coordinate frames of the image placements** were checked through the transform
  round trip and the image objects in `images.json` (for the LDIR overlay the origin is
  pixel 501 of 1002, half a pixel = 2.4 µm from the centre of a 1-based grid, absorbed
  by the alignment translation); the viewer's rotate / flip code was not traced line by
  line beyond its existing tests.

## Verification of the fixes

- Baseline (`5eaf53e`): full suite 1468 passed, 0 failed, 2 skipped (`writexl`), 2 warnings
  (no display for Tk; a missing file in a label test).
- After the fixes: every CSV, JSON and text table of a dataset A run is byte-identical
  to the baseline run, apart from timestamps, run paths and the new `md5` column of
  `supporting_files.csv` (49 files identical; compared file by file, which is stricter
  than the fingerprint).
- The targeted test files of each fix pass at its commit.
- Full suite (`testthat::test_dir("tests/testthat")`, including the dataset A
  fingerprint test) run on a checkout of each fix commit: `4d920a2`, `ef36576`,
  `df209bd`, `24a964d`, `b32317f`, `e040384` and `81b8a13` all finish with 0 failures;
  the only skips are the two `writexl` tests. `c952773` was not run on its own: its new
  test had a broken fixture, corrected in `ef36576`. `eba49d8`, `389ec3a` and `3c3c218`
  change no code and are covered by the run of the next commit.

## Follow-up: H-01 to H-03 implemented

- **H-01** (`44f5149`). Only an abbreviation in a trailing bracket (upper case / digits,
  "(PC)") is dropped from the exact-polymer name; Exact needs equal families.
  "Polyamide (naturally occurring)" vs PA and "Polyethylene (chlorinated)" vs PE are now
  *Disagree*; "Paraffin (PE)" vs PE is *Unresolved*. Not changed: natural polyamide
  against Protein / Wool / Silk is still *Disagree* (family Natural vs Protein).
- **H-02** (`c7bef98`). The LDIR Tier 2 pose is tested with `pose_significance()` after
  ICP. The 12 unrelated particles are now `determined = FALSE`, `pose_p = 1` as LDIR as
  well; nothing is scored or paired through them. Accepted landmark poses are unchanged.
- **H-03** (`56ba751`). The share of true pairs in `p_true` comes from the reference
  density; the fitted σ, the gates and the pairs are unchanged. Dataset A:

  | Pair | Σ(1 − p_true) before → after | Chance estimate | Pairs with p_true ≥ 0.9 | Identified, chance-corrected (raw) |
  |---|---|---|---|---|
  | FTIR (PE) – Raman | 1.44 → 0.25 | 0.6 | 29 → 30 of 31 | 89.6% → 93.2% (93.9%) |
  | FTIR (Bruker) – Raman | 3.03 → 0.35 | 0.5 | 14 → 16 of 18 | 62.4% → 73.6% (75.0%) |
  | LDIR – Raman | 0.96 → 0.05 | 0.1 | 26 → 28 of 28 | 87.2% → 90.2% (90.3%) |

  The three pairs still below 0.9 are `MP_8`–`97` (0.78) and the Bruker pairs `MP_20`–`8`
  and `MP_4`–`42`: distant pairs of large or mis-sized particles (M-12).
- Dataset A after the three changes: every pair, distance, transform and agreement tier
  is identical to the baseline; only the `p_true` columns (match and agreement files) and
  `n_found_true` / `rate_true_pct` in `detection_rates.csv` differ.
- Statistical review, updated: the chance-corrected identified rate is now consistent
  with `chance_pairs.csv` and can be quoted; `p_true` still assumes a normal error for
  the distance term, so it remains conservative for individual distant pairs.

## Follow-up 2: decisions on the open findings

All decisions were the author's (answers to 20 questions); the commits are on
`audit/fable-fixes`. Entries marked **(results)** in the CHANGELOG change numbers.

| ID | Decision | Commit | Changes results |
|---|---|---|---|
| M-12 | Extent rule (one centre inside the other particle) and additional pairs (several particles of one instrument on one of another); dataset A may change | `327242d` | yes |
| M-08 | New status `partner_not_matchable` | `41a11e2`, `e998b66` | partner status only |
| M-09 | Native cellulose natural, "Regenerated cellulose" semi-synthetic; natural polyamide is Protein | `843a2de` | yes (field samples) |
| M-15 | `opus_ftir_bruker_image = NULL`; sample names stay | `8dd0281` | runs that relied on the default |
| M-14 | Python package versions in `session_info.txt` | `d5a6b3e` | no |
| M-07 | Mirrored LDIR pose accepted and flagged (WARN) | `d5b49ea` | no |
| M-11 | WARN for a join without `.sam`; join check unchanged | `d5b49ea` | no |
| M-10 | One search-corrected `pose_p`; verdict from `determined` | `3074fe3` | `pose_p` values |
| M-13 | Overlay counts unmatched like the log | `8673f90` | viewer only |
| L-04 | `uploaded_data` removed | `5cb06d9` | viewer only |
| L-05 | Full LDIR transform file; one method name | `e3e27cb` | no |
| L-09 | Manifest paths relative, no user name | `8c79e57` | no |
| M-16 | Fingerprint extended to 918 values | `87843d4` | no |
| L-06, L-07, L-08, L-11 | Fixed | `8008eee`, `4f5fc41`, `31c6f84`, `7a7383e` | no on dataset A |
| L-10 | Deferred by the author (ORCID, release date) | | |
| H-04 | The author removes the file from the history | | |

### M-12 in detail

- **Rule.** A candidate outside the distance gate is accepted when the distance is at
  most 0.5 × the larger major axis (one centre lies inside the other particle), the
  Feret max agree within a factor 2, and the Poisson chance of an unrelated particle of
  comparable size within that radius is at most 20%. The assignment stays one-to-one.
- **Additional pairs.** After the assignment, a particle without a partner whose centre
  lies within the extent of an already paired particle (same chance limit, size ratio
  ≤ 5) is attached to it. Both particles stay separate rows; `pair_type` is
  `additional`. Agreement and triplets use the primary pairs; detection rates and
  partner status count every paired particle as found.
- **The chance limit was added after a test.** Without it, two unrelated random clouds
  of 5–150 particles gained up to 61 pairs from the rule. With it, a dense synthetic
  field of 120 + 120 large particles gains none.
- **Dataset A.**

  | | before | after |
  |---|---|---|
  | FTIR (PE) – Raman pairs | 31 | 33 (`MP_10`–`1` at 112 µm, `MP_1`–`80` at 81 µm) |
  | LDIR – Raman pairs | 28 | 31 (`A1`–`1`, `A4`–`80`, `A13`–`97`) |
  | FTIR (Bruker) – Raman pairs | 18 | 21 (two by extent; one additional: `MP_26` on Raman `1`, which already has `MP_27`) |
  | Triplets | 28 | 31, all with full agreement |
  | FTIR (PE) particles > 200 µm found by Raman | 11 of 13 | 13 of 13 |
  | LDIR particles > 200 µm found by Raman | 7 of 10 | 10 of 10 |
  | FTIR (PE) particles found by Raman, all sizes | 93.9% | 100% |
  | Raman particles found by LDIR, all sizes | 77.8% | 86.1% |

  These are exactly the particles traced by hand in section 2. `p_true` of the new
  pairs is 0.98–0.99 (0.70 for the additional pair).
- **Limits.** The extent is approximated by half the major axis around the centre, not
  the real outline. With the `.sam` outlines an LDIR pair could be tested by
  point-in-outline instead; that needs a real `.sam` file to build and test against.
  The greedy matcher (`match_method = "greedy"`) does not use the rule.
  `docs/validation.md` still quotes the dataset A pair counts from before this change in
  places; re-run `examples/validate_example_A.R` to refresh it.

### Notes on the other items

- **M-10.** The landmark acceptance still uses the uncorrected p < 1e-8. With very many
  landmarks on both sides the corrected p of an accepted pose can exceed 0.01; the
  verdict then still follows `determined`.
- **L-09.** A file chosen with the file picker outside the working directory keeps its
  full path in the manifest: the multi-run tools read raw files through it.
- **L-11.** An escaped `"\ufeff"` in `sub()` fails in a C locale ("input string 1 is
  invalid"); the BOM is matched by bytes instead.
- Three commits carried a test whose fixture was wrong and were followed by a correcting
  commit (`c952773` → `ef36576`, `41a11e2` → `e998b66`); history was not amended.
- **Not verified in a browser:** the viewer after L-04, M-08, M-13 and L-09. Checked
  headlessly: a dataset A run and a moved copy of it load, the overlay summary reads
  "33 FTIR↔Raman pairs | 0 FTIR unmatched | 5 Raman unmatched | 21 Bruker↔Raman pairs |
  4 Bruker unmatched | 31/40 LDIR↔Raman matched", which matches the log.

Full suite on the branch head (`7a7383e`): result added in the next commit.
