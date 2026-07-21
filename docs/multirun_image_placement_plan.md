# Plan v3: Multi-Run image↔coordinate placement — current status

**Status: FIXED. The lead hypothesis was confirmed on the user's machine
(see "v3 diagnostic findings" below) and the fix is implemented — see "v3 fix
implemented" at the end of this section. This section supersedes v2's §2-4
(kept below as historical context) now that its central open question has
been answered by production log evidence. v1 is Appendix A.**

## v3 update — the "P3 fallback" theory is DISPROVEN; new lead hypothesis

**New evidence (from the user's console log on the single Raman tab):**

```
[Particle Viewer] Raman image placed from WITec extent (Center Y as reported; 100% of particles inside).
```

This message can only be emitted from `raman_native_image_info`'s **Priority 1**
branch (`shiny_app/app.R:1580`, confirmed by reading the code). It fires
**twice** in the log — once on run load, once after the user additionally
uploaded `raman_image_canonical.png` through the tab's own upload control
(downsized 8956×8828 → 1791×1765 there, same as everywhere else) — and both
times reports **100% of run 1's particles inside the WITec-defined box**,
using the center-Y value **as configured, not negated**. That is a real,
correctly-calibrated placement — not a coincidental fallback. **v2's §3
hypothesis (single tab silently using P3/fit-to-particles) is wrong; discard
it.**

**The zoom-state theory is also ruled out** — the user confirmed the single
Raman tab is not zoomed and is displaying the full (downsampled) image at
100%, yet still visually looks different in scale from the Multi-Run tab
showing "the same" image.

**So: P1 (WITec) is confirmed correct and active on the single tab. The
question is now why Multi-Run — also supposedly using P1 with WITec values —
produces a different, wrong-looking scale for what should be the identical
placement math.**

### Lead hypothesis (check this FIRST — cheap, concrete, high-probability)

`R/00_config.R` explicitly documents `raman_image_width_um` /
`_height_um` / `_center_x_um` / `_center_y_um` as **per-dataset values**
("PER-DATASET value — read them from WITec's Particle Scout for each new
Raman scan and update them here, or the viewer falls back to heuristic
placement"). This means:

- The single Raman tab reads its WITec numbers from **`active_manifest()$config_snapshot`**
  — the config that was actually in effect **when that specific pipeline run
  was processed** (baked into that run's `manifest.json` at run time).
- The reproducibility tool's `CONFIG$raman_image_*` fields were filled in by
  the user **by hand**, presumably by copying whatever is **currently** in
  `R/00_config.R` — which may have since been edited for a *different* scan,
  or may simply not be the same values that were in effect for the specific
  raw Raman file(s) used as the three reproducibility replicates.
- **Quantitative smoking gun:** the WITec width in the current `R/00_config.R`
  is `12471.2`. The visually-apparent scale mismatch has consistently looked
  like ~2.5×. `12471.2 / ~5000 ≈ 2.5`. This is exactly the ratio you'd expect
  if Multi-Run is using stale/wrong-scan WITec values roughly 2.5× too large
  for the actual replicate dataset's true calibration.
- It is also possible the pipeline run currently selected in the single
  Raman tab (`2026-07-16_1`, the newest) is not even the **same physical scan**
  as the reproducibility tool's "run 1" raw file — in which case their WITec
  values could legitimately differ and comparing them was apples-to-oranges
  from the start. This needs to be checked, not assumed either way.

### Sharpened diagnostic step (replaces v2 §4 as the first thing to do)

1. Identify exactly which raw Raman file was used as
   `CONFIG$runs[[1]]$file` in `tools/reproducibility.R` for this comparison
   (its filename, and ideally its MD5 — `tools::md5sum(path)`).
2. Find the pipeline `output/<run_dir>/` whose inputs correspond to that same
   raw file: open each candidate run's `manifest.json` and check
   `inputs.raman.basename` (and/or `inputs.raman.md5`, which is more
   reliable if the file was ever renamed/copied) — confirmed field names,
   written by `write_manifest()` in `R/utils.R:783-799`. **Do not assume**
   it's "the newest run" (`2026-07-16_1`) or "whatever's currently selected
   in the single tab" without actually matching the file.
3. Open that specific run's `manifest.json` and read
   `config_snapshot.raman_image_width_um`, `_height_um`, `_center_x_um`,
   `_center_y_um` — these are the numbers that produced the "100% inside,
   correct" placement.
4. Diff those against what's currently entered in `tools/reproducibility.R`'s
   `CONFIG$raman_image_*` fields.
5. **If they differ:** that's almost certainly the entire bug. Fix is trivial
   — enter the exact manifest values into the reproducibility tool's CONFIG
   (or, better, have the tool read them automatically from the corresponding
   run's manifest instead of requiring hand entry — this removes the
   copy-paste/staleness failure mode entirely for future datasets). Re-run
   the tool, reload Multi-Run, compare.
6. **If they match exactly** and the problem still reproduces, THEN fall back
   to v2's original §4 (add temporary diagnostic logging to both code paths
   and diff the actual computed extents/tier numbers) and §5 (extract the
   native placement logic into a function shared by both tabs) — the
   duplication risk documented in v2 is still real and still worth fixing
   regardless of whether it's the cause of this specific symptom.

### Design implication worth carrying into the fix

Regardless of what the diagnostic finds, **requiring the user to hand-copy
per-dataset WITec values into a second place (the reproducibility tool's
CONFIG) is itself a latent bug generator** — it's exactly the kind of value
that drifts silently. If feasible, the tool should instead resolve these
values automatically from the pipeline run corresponding to each replicate
file (matching on source filename against existing `output/*/manifest.json`
entries), falling back to manual CONFIG entry only when no matching run is
found. This is a stronger, more permanent fix than "remember to keep two
config blocks in sync."

### v3 diagnostic findings (remote session 2026-07-21)

The decisive numeric comparison (steps 1–4 above) **cannot be completed from
the repository alone** — all three artifacts it needs exist only on the
user's machine:

1. The locally-edited `tools/reproducibility.R` CONFIG (the committed copy
   has all `raman_image_*` fields `NULL` and placeholder run paths).
2. The pipeline run for the problem dataset (e.g. `2026-07-16_1`) — the
   committed `output/` holds only `2026-03-06_2` and `2026-03-09_3`, which are
   a **different dataset** (Raman input is a CSV with an 892×864 PNG image,
   not the 8956×8828 BMP), and whose `config_snapshot` contains **no
   `raman_image_*` keys at all** (they were `NULL` at run time; `write_manifest`
   drops NULLs).
3. Any `output/reproducibility/*/reproducibility_meta.csv`.

What the repo/history **does** establish, all supporting the lead hypothesis:

- **The per-dataset drift mechanism is real and recent.** `R/00_config.R` git
  history shows two different WITec calibration sets within two days:
  `12569.2097 / 12153.0107 / -272.4737 / +7277.1328` (until 2026-07-09,
  removed in `21a401b`, which set the fields NULL with `raman_um_per_px = 0.38`)
  → `12471.2 / 12313.0 / 1487.6 / -4394.0` (activated 2026-07-10 in `6857109`,
  still current). Note: **neither historical set is ~5000 µm wide** — if the
  observed mismatch really is ~2.5× in scale, the replicate scan's true
  calibration has never been committed anywhere; it exists only in WITec's
  Particle Scout panel and possibly the local run's manifest.
- **"100% of particles inside" does not validate the scale.**
  `raman_image_extent_from_config()` accepts any extent with
  `frac_inside >= 0.5`; an oversized box trivially contains 100% of the
  particles, so the containment check can reject too-small/mispositioned
  values but *structurally cannot reject too-large ones*. Both tabs would
  silently accept a 2.5×-too-big calibration.
- **The P1 math and the framing genuinely cannot diverge.** Both the single
  tab (`app.R:1575`) and Multi-Run (`place_image_raman_meta`,
  `global.R:1095`) call the same `raman_image_extent_from_config()`, and both
  frame the plot on the image extent + 200 µm pad (`app.R:2761` vs
  `app.R:4534`). The Multi-Run manual Width/Height override defaults to `NA`
  (inactive). So with identical WITec values and identical run-1 particles the
  two tabs mathematically render the same scale — the only remaining degrees
  of freedom are (a) the values themselves (meta vs manifest snapshot — the
  lead hypothesis), (b) the particle sets (the replicate raw file may not be
  the same physical scan the pipeline run ingested), or (c) a leftover manual
  Width entry in the Multi-Run sidebar during testing.
- **Traceability gap:** `reproducibility_meta.csv` records no source
  filenames/MD5s, so a repro output cannot be traced back to its raw files.
  The auto-resolve fix should also record `run<i>_file` / `run<i>_md5` in the
  meta for future diagnosability.

**New tool committed for the decisive check:**
`tools/diagnose_multirun_placement.R` (read-only) automates steps 1–4 on the
user's machine:

```
Rscript tools/diagnose_multirun_placement.R <path/to/run1-raman-file> [repro_output_dir]
```

It prints (A) the WITec values Multi-Run actually uses (newest
`reproducibility_meta.csv`), (B) every pipeline run's
`inputs.raman.basename/md5` + `config_snapshot.raman_image_*`, flagging the
run whose MD5 matches the replicate file, (C) the current `R/00_config.R`
values, and a per-field MATCH/MISMATCH verdict that maps directly onto step 5
(values fix) vs step 6 (fall back to v2 §4/§5).

### User-run diagnostic result: MISMATCH confirmed (2026-07-21)

The user ran `tools/diagnose_multirun_placement.R` against their real
`output/` (17 pipeline runs, not the 2 committed to this repo) with the PET-A
replicate file (`Comparstic2026 PET A Kev 260417.csv`). Findings:

- **The raw file was reprocessed repeatedly with drifting calibration.**
  MD5-matching runs for this exact file show at least three distinct WITec
  calibrations over six weeks: `<NA>` (2026-05-20 to 2026-07-09, P3 fallback)
  → `5133.146/5068.031/2296.34/2355.731` (2026-07-09_1) →
  `12471.2/12313.0/1487.6/-4394.0` (2026-07-09_4, matches current
  `R/00_config.R`) → `5237.904/5171.46/2314.045/2354.009` (2026-07-10_1
  onward, stable through the most recent run, 2026-07-15_1).
- **Multi-Run's `reproducibility_meta.csv` had `12471.2/12313.0/1487.6/-4394.0`**
  — a straight copy of whatever was in `R/00_config.R` when the user last
  edited `tools/reproducibility.R`'s CONFIG. That is the value from
  2026-07-09_4, six days stale relative to the run the single tab is actually
  showing (2026-07-15_1, `5237.904/...`) — a **2.38x width ratio**
  (12471.2 / 5237.904), matching the observed compression almost exactly.
- **Verdict: MISMATCH CONFIRMED.** The lead hypothesis (manual hand-copy of a
  per-dataset value going stale) is the entire root cause. No further
  investigation into v2 §4/§5 (instrumentation, shared-function extraction)
  is needed for this symptom — that work remains valid future hardening but
  is no longer required to fix this bug.

### v3 fix implemented (2026-07-21)

Per step 5's stronger recommendation, `tools/reproducibility.R` now
**auto-resolves** the Raman calibration instead of requiring hand entry:

- `resolve_raman_calibration_from_manifests()` (`R/utils.R`) takes the raw
  Raman file used as `CONFIG$runs[[1]]$file`, MD5-matches it against every
  `output/*/manifest.json`'s `inputs.raman.md5`, and returns the
  `config_snapshot.raman_image_*` values from the **latest** matching run
  (latest = most recent processing of that exact file, i.e. the user's most
  recent belief about the correct calibration for that scan). Returns `NULL`
  if no run matches or the matching run's snapshot is incomplete.
- `tools/reproducibility.R` calls this for `CONFIG$instrument == "raman"`
  before writing output. Priority order: **auto-resolved from manifest** →
  **manual `CONFIG$raman_image_*` fields** (now a fallback, not the primary
  path) → none (particle-extent fit). The resolution outcome is logged and
  recorded in `reproducibility_meta.csv` as `raman_calibration_source`
  (`"manifest:<run_id>"` / `"manual"` / `"none"`).
- `reproducibility_meta.csv` now also always records `run1_file`/`run1_md5`
  (any instrument) so a reproducibility output can be traced back to its
  exact raw file — closing the traceability gap noted in the original v3
  lead hypothesis write-up.
- `tools/diagnose_multirun_placement.R` uses the new `run1_md5` field to
  auto-detect the replicate file on a re-run, so it needs no arguments once a
  reproducibility output has been produced by the fixed tool.
- Regression coverage: `tests/testthat/test-raman-calibration-resolve.R`
  (MD5 match + latest-wins, no-match, incomplete-snapshot, and
  nonexistent-file/dir cases). Full suite (`Rscript run_tests.R`) passes.

**Remaining step for the user:** re-run `tools/reproducibility.R` for the
Raman instrument (no CONFIG edits needed — the calibration now resolves
automatically) and confirm the Multi-Run tab's image scale now matches the
single Raman tab.

---

# Plan v2 (superseded — kept for its still-valid architectural analysis, now Phase 2 material)

**Original status when written: plan only. Its central diagnostic question
(§3/§4 below) has been answered — see the v3 update above. v2's §5
(extract-and-share architecture) and its non-goals/guardrails remain valid
and should still be done if the lead hypothesis above doesn't fully resolve
the symptom, or as follow-up hardening regardless.**

**Why v2 exists:** v1 was implemented in good faith but iterated on *screenshots*
without ever positively confirming what the single-instrument tab actually
computes for the problem dataset. Four sequential "fixes" each looked plausible
and each turned out to be wrong or incomplete:
1. Frame on particle extent → cropped to a sub-region of the image.
2. Frame on image extent, width-only placement → stretched (non-square pixels).
3. Width + height (WITec values) → placement scale still ~2.5x too large,
   particles compressed into the image centre.
4. Add the "P2 µm-per-pixel" tier → **did not fix it**, because for a BMP source
   `extract_tiff_um_per_px()` is TIFF-only and returns `NULL` — P2 is dead code
   for this dataset on *both* tabs. Uploading the canonical PNG (same pixel
   dimensions as the original BMP, confirmed via the pipeline's own manifest
   metadata table) **still shows the same compression**, which rules out image
   *decoding* as the cause and points back at the placement *math* or at stale
   config being reused across runs of the tool.

This pattern — plausible fix, still wrong — means the remaining work needs to
stop guessing from renders and instead **instrument both code paths and
compare actual numbers**, then **fix the duplication that let them diverge in
the first place.**

---

## 1. What we know for certain (verified by reading the code, not by inference)

| Fact | Where | Confirmed how |
|---|---|---|
| Raman source image for this dataset is BMP, 8956×8828 px | user-reported, from the pipeline's manifest metadata panel | direct user report |
| `canonical.png` is a full-resolution re-encode (same 8956×8828 dims) | `R/utils.R:canonicalize_instrument_image` writes `<inst>_image_canonical.png` at full res | code read this session |
| `preview.png` is a downscaled copy (~22% in this case) | same function, writes `<inst>_image_preview.png` | code read this session |
| `extract_tiff_um_per_px()` returns `NULL` for any non-TIFF file, and also needs `magick` | `shiny_app/global.R:990-993` — `if (sniff_image_type(path) != "TIFF") return(NULL)` | code read this session |
| ⇒ Raman placement **Priority 2** (µm-per-px) is structurally unreachable for this BMP-sourced dataset, on **both** the single Raman tab and the Multi-Run tab | derived from the above | logical consequence, not yet empirically confirmed inside a running app |
| `raman_native_image_info` (single tab) is a 3-tier cascade: P1 WITec extent → P2 µm/px (TIFF only) → P3 `compute_image_bounds` (aspect-preserving fit to the particle extent, 300 µm pad) | `shiny_app/app.R:1554-1605` (native reactive) | code read this session |
| ⇒ **For this dataset specifically, if no WITec metadata is configured in the pipeline run's config, the single Raman tab is almost certainly falling through to P3** (fit-to-particles), not P1 or P2 | derived from the two facts above | **not yet confirmed** — needs direct instrumentation (§4) |
| The Multi-Run tab's Raman placement (`place_image_raman_meta` → `place_image_raman_umpx` → caller's `compute_image_bounds` fallback) is a **separate reimplementation**, not a call into `raman_native_image_info` | `shiny_app/global.R` (helpers added this session) | this session's own edits |
| `options(shiny.maxRequestSize = 50 * 1024^2)` — 50 MB upload cap is already set | `shiny_app/global.R:17` | code read this session |
| The reproducibility tool's `reproducibility_meta.csv` is **overwritten wholesale on every run** (`write.csv(meta, ...)`, no merge) | `tools/reproducibility.R` | code read this session — so stale values are NOT the mechanism *if* the user re-ran the tool after editing CONFIG. **This has not been confirmed with the user** — ask directly whether they re-ran the tool (not just reloaded the viewer) after clearing/changing the WITec fields. |

## 2. Where v1's reasoning broke down

v1 assumed the single Raman tab's displayed extent (~5000 µm, read off the
screenshot axes) meant P1 or P2 must be active with a value near 5000. That
was never verified — it was inferred backwards from a rendered plot. The
actual mechanism is very likely P3 (fit-to-particle-extent with padding),
which *coincidentally* also produces a plot that looks like "the image roughly
matches the particle spread" — because that's literally what P3 does by
construction (it force-fits the image to wherever the particles are, with a
7300 padding). If that's correct, then:

- The **single Raman tab's apparent alignment is not evidence of a WITec or
  µm/px value being correct** — it's the fallback doing its job (fitting the
  image to the particles), which trivially "looks right" but carries **no
  real spatial information** — the image could be scaled arbitrarily and P3
  would still make it "fit."
- If that's the mechanism, then the Multi-Run tab reaching for real physical
  metadata (WITec width/height) was **never going to visually match** P3's
  fit-to-particles behaviour, because they encode fundamentally different
  claims about the image (true physical scale vs. "make it fit
  cosmetically"). One of these is right and one is wrong for this dataset —
  determining which is priority zero (§3).

This is the central open question the next session must resolve before
writing any more placement code.

## 3. The open question that gates everything else

**Does this Raman dataset actually have a reliable physical scale (WITec
metadata, or a known µm/px), or does the single Raman tab currently show a
cosmetically-fit image with no real spatial meaning (P3 fallback)?**

- If **P3 is what's "correct"** (i.e., the user has always visually accepted
  a fit-to-particles image, and there genuinely is no reliable calibration
  for this scan) → the fix is simple: make Multi-Run call P3 with the exact
  same inputs (particle set, padding) as the single tab, and stop trying to
  apply WITec/µm-per-px for datasets that don't have real calibration. No
  further physical-accuracy work needed; this was over-engineered.
- If **there IS a real WITec calibration** for this scan that the pipeline
  run's config doesn't currently have entered (i.e. the single Raman tab is
  *unknowingly* falling back to P3 when it shouldn't, and the "true" image is
  actually offset/mis-scaled in the single tab too) → that's a **separate,
  pre-existing bug** in the single Raman tab, not something to reproduce in
  Multi-Run. Surface it to the user rather than chase it further in
  Multi-Run.
- Either way, Multi-Run should track whatever the single tab **actually
  does**, not what we assume it does.

## 4. Mandatory first step: instrument, don't guess

Before writing any placement code, add temporary diagnostic output (console
`log_message()` / `message()`, or a debug text panel in the UI) to
**both** `raman_native_image_info` (single tab) and the Multi-Run placement
call, printing, for the exact same dataset/run:

- which tier fired (P1 / P2 / P3)
- the resolved `um_per_px` value (if P2)
- the WITec config values seen (if P1) — confirm whether they're `NULL`/absent
  for this run's pipeline config (this determines whether P1 is even reachable
  on the single tab)
- the final `xmin/xmax/ymin/ymax` extent
- the raster's pixel dimensions as loaded

Run both tabs against the identical run-1 dataset and **diff the printed
numbers**, not the rendered plots. This single step will conclusively answer
§3 and should take under 15 minutes — do it first, before touching any
placement logic.

Also worth checking directly: open the pipeline run's `manifest.json` for this
Raman dataset and read `config_snapshot.raman_image_width_um` /
`_height_um` / `_center_x_um` / `_center_y_um` — if all four are `NULL`/absent,
P1 is confirmed dead on the single tab for this run, which is strong evidence
for the P3 hypothesis.

## 5. Architectural fix — eliminate the duplication that caused this

Regardless of what §4 finds, **the reimplementation itself is the root
process failure**: v1 built parallel logic (`place_image_raman_meta`,
`place_image_raman_umpx`, a separate `compute_image_bounds` fallback call) that
*approximates* `raman_native_image_info` instead of calling it. Every
iteration this session added a "tier" to the copy that already existed in the
original — and still didn't match, because copies drift.

**Do this instead: extract the native reactive's *body* into a plain
function, and call that same function from both places.**

Concretely, for Raman:

```r
# Pure function - no reactive context, callable from anywhere.
compute_raman_image_info <- function(raw, x_orig, y_orig, cfg, um_per_px_config,
                                     raman_image_path, ox = 0, oy = 0) {
  # ...move the existing 3-tier body of raman_native_image_info here verbatim...
}
```

Then:
- `raman_native_image_info <- reactive({ compute_raman_image_info(raman_raw_image(), ...) })`
  — unchanged behaviour, just delegates.
- The Multi-Run tab calls `compute_raman_image_info(run1_raster, run1_x, run1_y, cfg_from_meta, ...)`
  with run 1's data and whatever config the reproducibility tool recorded.

Same treatment for FTIR/Bruker (`place_image_particle_extent` should become
the literal body of `ftir_native_image_info`, not a hand-written
reimplementation — check whether it already matches exactly; it's simpler so
it's lower risk, but verify) and LDIR (`ldir_native_image_info`'s circle
branch and its symmetric fallback both need to be the literal shared function,
including the fallback — v1 only ported the circle branch).

This is strictly more work than v1's approach but is the only way to
**guarantee** the two tabs can't diverge again — every future change to
Raman/FTIR/LDIR placement logic then automatically applies to Multi-Run too.

## 6. Also verify: is `reproducibility_meta.csv` actually fresh?

Confirm with the user whether they re-ran `tools/reproducibility.R` (not just
reloaded the Shiny app / re-selected the run folder) after each config change.
The tool always creates a **new timestamped subfolder** per run (from an
earlier fix this session) — if the user has been re-loading an **old** run
folder in the viewer's dropdown after editing CONFIG without re-running the
tool, they'd be looking at stale metadata from before the WITec fields were
set. Cheap to rule out: check the `run` dropdown's selected timestamp against
when they last edited the tool config.

## 7. Step-by-step execution order for the next session

1. **Instrument and diff** (§4) — confirm which tier fires on the single tab
   for this dataset, with exact numbers. Check `manifest.json`'s
   `config_snapshot` directly too.
2. **Confirm fresh metadata** (§6) — rule out stale-run confusion.
3. Based on §1's finding, decide: is P3 "correct" for this dataset, or is
   there a genuine calibration bug in the single tab? Get the user's read on
   whether the single tab's current image placement actually looks
   *physically* correct to them (features aligned with real particle
   positions they can independently verify) or just *cosmetically parked*
   over the particles.
4. **Do the architectural extraction** (§5) for whichever instrument(s) are in
   scope — start with Raman since it's the one in front of us, but the same
   duplication risk exists for FTIR/Bruker/LDIR helpers added this session;
   audit and convert all four.
5. Re-run the reproducibility tool, reload the Multi-Run tab, and compare
   against the single tab **using the diagnostic output from step 1**, not
   just visually — confirm the numbers match exactly, then confirm visually.
6. Remove or gate the temporary diagnostic logging behind `config$debug` /
   an existing debug flag before committing (don't ship noisy logging).
7. Update/extend `tests/testthat/test-multirun-placement.R` to assert the
   shared function returns identical output whether called from the "native
   tab" code path or the "Multi-Run" code path, using a shared fixture (this
   test would have caught v1's drift immediately).

## 8. Non-goals / guardrails for whoever picks this up

- **Do not** re-touch the run-directory image downsampling revert from
  earlier this session (`d0243a1` reverted to full-res loading for
  Overlay/instrument tabs) — that was a separate, already-resolved bug
  (block-average distortion), unrelated to this placement issue. Full-res
  loading must stay.
- **Do not** re-derive the non-square-pixel fix (explicit width *and* height
  fields) — that part was correct and should be preserved once the underlying
  metadata question (§3) is resolved.
- **Do not** iterate on more "tiers" or "fallbacks" bolted onto the existing
  duplicated helpers. If §5's extraction isn't done, any further placement fix
  is liable to repeat this exact failure mode a fifth time.
- The manual Width/Height/Offset/Rotate controls in the Multi-Run sidebar
  should remain as an override path (useful when there's genuinely no
  metadata, e.g. a dataset with no WITec info and the user wants to eyeball a
  fit) — just make sure they sit clearly *after* the shared-function result in
  the priority order, and are visually labeled as a manual override, not the
  default.

## 9. File/function map (for orientation)

| File | Symbol | Role |
|---|---|---|
| `shiny_app/app.R` | `raman_native_image_info` (reactive, ~line 1554) | Single Raman tab's placement — the ground truth to extract from |
| `shiny_app/app.R` | `ftir_native_image_info`, `ftir_bruker_native_image_info` (~line 1519, 1535) | FTIR/Bruker placement — extract similarly |
| `shiny_app/app.R` | `ldir_native_image_info` (~line 1678) | LDIR placement — extract similarly, including its fallback branch |
| `shiny_app/app.R` | `output$repro_plot` (Multi-Run tab renderPlot) | Where the extracted functions should be called from, replacing the current `place_image_*` calls |
| `shiny_app/global.R` | `raman_image_extent_from_config`, `place_image_*` (added this session) | Existing partial extraction — audit, complete, or replace per §5 |
| `shiny_app/global.R` | `compute_image_bounds` | The P3 fallback fit — used by both tabs already, lowest risk |
| `shiny_app/global.R` | `extract_tiff_um_per_px` | Confirmed TIFF-only; dead for this dataset |
| `tools/reproducibility.R` | `CONFIG$raman_image_*`, meta-writing block | Where Raman metadata is entered/recorded — revisit once §3 is resolved |
| `R/utils.R` | `canonicalize_instrument_image` | Produces `canonical.png` (full-res) / `preview.png` (downscaled) — already confirmed correct, not a suspect |

## 10. Acceptance criteria

- A single, shared placement function per instrument, called identically by
  the native tab and the Multi-Run tab (verified by a unit test, not just eyes).
- For the problem Raman dataset: diagnostic output confirms which tier fires
  on the single tab, and the Multi-Run tab's rendered image visually matches
  it — including whether that match is a "true calibration" match or a
  "both fell back to P3" match (either is fine, as long as it's understood
  and both tabs agree).
- No regression to Overlay/instrument-tab image quality (full-res loading) or
  to the non-square-pixel width/height fix.

---

# Appendix A — v1 plan (superseded, kept for history)

Status when written: plan only. Decisions locked with the user at the time:
1. Metadata source = **(a)** entered in the reproducibility tool's config.
2. **Run 1 is the reference** image; runs 2–3 register to it (residual = jitter).
3. Must work for **every instrument** (Raman, FTIR PerkinElmer, FTIR Bruker, LDIR),
   not just Raman.

## Problem

The single-instrument viewers place the background image exactly, but the
Multi-Run tab guesses the placement (particle-mean centre + manual width/height/
offset), producing a large systematic offset. We want the Multi-Run tab to
reproduce each instrument's *native* image↔coordinate relationship.

## Root cause — each instrument places its image differently

From `shiny_app/app.R`:

| Instrument | Native reactive | Placement logic | Metadata needed |
|---|---|---|---|
| FTIR (PerkinElmer) | `ftir_native_image_info` | image spans the **raw particle extent** `min/max(x_orig,y_orig)` + offset | none (just run-1 coords) |
| FTIR (Bruker) | `ftir_bruker_native_image_info` | same particle-extent placement | none |
| Raman | `raman_native_image_info` | 3-tier: **(1) WITec extent** (center + width/height, Y auto-detect) → (2) `um_per_px` scale → (3) `compute_image_bounds` | `raman_image_width_um`, `height`, `center_x`, `center_y` |
| LDIR | `ldir_native_image_info` | **scan-circle calibration** (`cx_px`, `cy_px`, `scale_um_per_px`, image px dims) → fallback: symmetric ±scan_diameter/2 | circle calibration (or scan diameter) |

*(v1's error, discovered in v2 §2: this table is accurate as a description of
the code, but v1 wrongly assumed P1 or P2 was firing on the single Raman tab
for the problem dataset, when P3 is far more likely — see v2 §3-4.)*

**Rest of v1 elided here — its implementation pieces (metadata recording in
the tool, per-instrument placement branches, control demotion) are subsumed by
v2 §5's "extract-and-share" approach, which is stricter than v1's "write a
parallel implementation that's supposed to match" approach.**
