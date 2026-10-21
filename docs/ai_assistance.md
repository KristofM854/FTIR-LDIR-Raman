# Use of AI assistance

This software was developed by Kristof Moeller with extensive use of an AI
coding assistant: Anthropic's Claude, through Claude Code (GitHub
integration and cloud sessions), from February to October 2026. This page
records where it was used, how that was measured, and which decisions were
the author's.

## How it was measured

From the git history (`git log`, `git blame` on the release commit):

- Commits authored by "Claude" are AI-written: 277 of the 631 commits (21
  October 2026). The author's 354 commits are pull-request merges of those
  changes (206) and manual edits and file uploads (148).
- For every tracked text file (images, CSV data and `renv.lock` excluded),
  each current line is attributed to the commit that last changed it.

Limits of this measure: a line last touched by an AI commit may have been
written or dictated by the author before (and vice versa), and some commits
under the author's name (for example 2026-02-27, Shiny image provenance
fixes) contain code produced with AI tools outside this record, so the AI
share is if anything higher than shown.

## Where

Share of current lines last changed by an AI-authored commit:

| Area | lines | AI-authored | share |
|---|---|---|---|
| Pipeline (`R/`) | 14,880 | 14,795 | 99.4% |
| Viewer (`shiny_app/`) | 9,126 | 9,053 | 99.2% |
| Python LDIR detector and its tests (`inst/python/`) | 813 | 813 | 100% |
| Tests (`tests/`) | 4,968 | 4,968 | 100% |
| Tools and examples | 2,436 | 2,435 | 100% |
| Documentation and metadata | 1,831 | 1,829 | 99.9% |
| **Total** | **34,054** | **33,893** | **99.5%** |

Files with author-written lines: `shiny_app/global.R` (71 lines), `R/utils.R`
(66), `R/00_config.R` (14), `R/10_export.R` (4), `.gitignore` (2),
`shiny_app/app.R` (2), `R/00b_file_input.R` (1), `tools/reproducibility.R` (1).

In short: essentially all code, tests and documentation text were written by
the assistant; the author directed the work, reviewed and merged every change,
ran the pipeline on the instruments' data and supplied the datasets.

## Contributions

| Role (CRediT) | Who |
|---|---|
| Conceptualisation: the idea of matching particles across FTIR, Raman and LDIR on the same filter, and the comparison it serves | Author |
| Methodology: which steps the pipeline takes, the agreement rules, cutoffs and gates, what counts as validation | Author, with the assistant proposing options and the author deciding |
| Data: filters, measurements, instrument exports | Author |
| Software: implementation, tests, documentation | AI assistant, directed and reviewed by the author |
| Validation: running the checks, interpreting the results | Author, with the checks implemented by the assistant |

## What the assistant did

- Wrote the pipeline: readers, LDIR image analysis (R and Python), alignment
  (landmarks, RANSAC, global registration, ICP, TPS), matching, agreement
  scoring, material classification, export, report and the Shiny viewer.
- Wrote the tests, the validation tools and this documentation.
- Pre-release audit (October 2026, branch `audit/pre-release`): code review
  against the input data, the data contract and readers, the
  `run_pipeline()` refactor, the validation tools, the test suite, the
  release documents; it ran the pipeline on the example dataset and reported
  every change in the outputs.

## What the author decided

Scientific and design choices were put to the author as explicit questions
and decided by the author; the assistant implemented them. From the pre-release
audit:

- Low-quality and unidentified particles never enter the agreement
  denominator; cutoffs FTIR AAU 0.70, Raman HQI 70, LDIR 0.85.
- Forced LDIR joins are kept with a quality warning; unresolvable name pairs
  form a separate tier for human inspection.
- Agreement is reported at two levels: family (primary) and exact polymer
  (secondary).
- Alignment scale band 0.8–1.25 kept; FTIR↔Raman keeps a fitted scale.
- LDIR size discrepancy < 5% accepted without a manual scan-diameter entry.
- Acrylic polymers (incl. PAN) counted as synthetic.
- Python required for LDIR images (no R fallback); git history not rewritten.
- The spectral-library index is used and disclosed.
- The LDIR join check is logged in every run; the quality-cutoff sweep was
  dropped.
- Which datasets are published (dataset A).

From the first field samples (dense filters, thousands of Raman particles):

- Only particles passing their instrument's quality cutoff and >= 20 µm are
  matched; pairs whose sizes differ more than 5x are rejected.
- Distance gates 50 µm, and 20 µm for LDIR against the reference (about three
  times its fitted positional error; the earlier 150 µm gate let in ~40%
  chance pairs on a field sample).
- Landmarks are trusted only when far above chance; an alignment that the
  data cannot fix (a short, pre-filtered list) is reported as not determined.
- Every pair gets a probability of being the true partner; agreement and
  detection rates are reported with 95% intervals.
- Open identification questions are referred to the plastics group
  ([roadmap.md](roadmap.md)), not decided in code.

## Verification

AI-written code was checked by running it, not only by reading it: the test
suite (`tests/README.md`), the regression test on dataset A, synthetic
end-to-end runs with known answers, and the validation in
[validation.md](validation.md). The author inspected the results in the
viewer.
