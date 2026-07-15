# Test suite

Regression tests for the pipeline's critical, deterministic paths. The project
is run by sourcing `main.R` (it is not an installed package), so the tests
source the relevant `R/` modules directly — no build/install step.

## Running

```r
# from the repo root
Rscript run_tests.R
```

Requires `testthat`, `RANN`, `clue`, `jsonlite`, `withr` (the packages the
tested modules touch). They install as binaries on CI (`use-public-rspm`), and
`.github/workflows/tests.yml` runs the suite on every push to `main` and on
pull requests.

## What is covered

| File | Focus |
|------|-------|
| `test-utils-transform.R` | Similarity-transform round-trip, identity, reflection detection — the coordinate primitives every alignment path uses. |
| `test-determinism.R` | **B3 regression:** `ransac_align` / `global_register_align` are reproducible run-to-run *and* leave the caller's global RNG untouched. |
| `test-matching.R` | Hungarian 1-to-1 assignment on known clouds; the dedicated looser LDIR distance gate. |
| `test-manifest.R` | Provenance: `write_manifest` records config + environment and round-trips through JSON. |
| `test-try-or.R` | **B4 regression:** `try_or` logs on failure instead of swallowing silently. |

Fixtures are synthetic and self-contained (`tests/testthat/helper-setup.R`):
point clouds transformed by a *known* similarity, so recovery/match counts have
a ground truth without shipping instrument data. Determinism tests deliberately
shrink the RANSAC grid/iteration count — reproducibility is independent of how
exhaustive the search is.

## Not yet covered (tracked follow-ups)

- **Real-data golden pipeline test.** An end-to-end run over a small frozen
  3-instrument dataset with known-truth ingest/alignment/match/agreement
  outputs. The synthetic determinism + matching tests already guard against
  match-count drift; a real-data baseline additionally guards the ingest and
  agreement stages. Drop a minimal dataset under `tests/golden/` and assert the
  staged outputs against frozen CSVs.
- **Agreement scoring** unit tests (`R/08_agreement.R`).
- **Shiny viewer** smoke test (`shinytest2`): app boots, overlay renders, a
  filter change updates the point count.
