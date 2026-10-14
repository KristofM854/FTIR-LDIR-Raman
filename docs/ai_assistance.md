# Use of AI assistance

This software was developed by Kristof Moeller with extensive use of an AI
coding assistant: Anthropic's Claude, through Claude Code (GitHub
integration and cloud sessions), from February to October 2026. This page
records where it was used, how that was measured, and which decisions were
the author's.

## How it was measured

From the git history (`git log`, `git blame` on the release commit):

- Commits authored by "Claude" are AI-written. 257 of the 611 commits are;
  the author's 354 commits are almost all pull-request merges of those
  changes (319), plus small manual edits and file uploads.
- For every tracked text file, each current line is attributed to the commit
  that last changed it.

Limits of this measure: a line last touched by an AI commit may have been
written or dictated by the author before (and vice versa), and commits under
the author's name (for example 2026-02-27, Shiny image provenance fixes) may
contain code produced with AI tools outside this record. The author should
confirm those.

## Where

Share of current lines last changed by an AI-authored commit:

| Area | lines | AI-authored | share |
|---|---|---|---|
| Pipeline (`R/`) | 13,680 | 13,595 | 99.4% |
| Viewer (`shiny_app/`) | 8,538 | 8,463 | 99.1% |
| Python LDIR detector and its tests (`inst/python/`) | 813 | 813 | 100% |
| Tests (`tests/`) | 4,634 | 4,634 | 100% |
| Tools and examples | 2,566 | 2,565 | 100% |
| Documentation and metadata | 2,744 | 2,742 | 99.9% |
| **Total** | **32,975** | **32,812** | **99.5%** |

Files with author-written lines: `R/utils.R` (66 lines), `shiny_app/global.R`
(73), `R/00_config.R` (14), `R/10_export.R` (4), `.gitignore` (2).

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
- Distance gates 100 µm, LDIR 150 µm; LDIR size discrepancy < 5% accepted
  without a manual scan-diameter entry.
- Acrylic polymers (incl. PAN) counted as synthetic.
- Python required for LDIR images (no R fallback); git history not rewritten.
- The spectral-library index is used and disclosed.
- The LDIR join check is logged in every run; the quality-cutoff sweep was
  dropped.
- Which datasets are published (dataset A).

## Verification

AI-written code was checked by running it, not only by reading it: the test
suite (`tests/README.md`), the regression test on dataset A, synthetic
end-to-end runs with known answers, and the validation in
[validation.md](validation.md). The author inspected the results in the
viewer.
