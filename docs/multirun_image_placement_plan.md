# Plan: exact image↔coordinate placement in the Multi-Run tab (all instruments)

Status: **plan only — not implemented.** Decisions locked with the user:
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

The Multi-Run tab currently uses `compute_image_bounds` (aspect-fit, particle-mean
centre) for all of them — which matches none exactly.

## Unifying principle

**For each instrument, the Multi-Run tab reuses that instrument's native
placement logic, fed with run 1's coordinates plus the metadata the tool
records for that instrument.**

Coordinate frame is consistent by construction: in the reproducibility engine
run 1 is the reference, so `x_aligned == x_um == x_orig` (raw). That is the same
frame each native `*_image_info` places into, so:

- FTIR/Bruker: `min/max(run1 x,y)` reproduces the native extent directly.
- Raman: `raman_image_extent_from_config(cfg, run1_x, run1_y)` returns the extent
  in the raw stage frame; the native tab's centroid shift `mean(x_orig - x_norm)`
  is **zero** here (run 1 isn't normalized), so no shift is needed. Y-up/Y-down
  auto-detection comes for free from the reused function.
- LDIR: the circle formula produces coords centred at 0, which is exactly the
  joined LDIR frame run 1 already sits in.

## Implementation

### Piece 1 — Tool records per-instrument placement metadata
File: `tools/reproducibility.R` (+ maybe a small helper in `R/reproducibility.R`).

Extend `reproducibility_meta.csv` with an instrument-appropriate block, driven by
the tool's CONFIG (option **a**):

- **Raman**: add `raman_image_width_um`, `raman_image_height_um`,
  `raman_image_center_x_um`, `raman_image_center_y_um` to CONFIG; write them to meta.
- **LDIR**: the tool already calls `extract_ldir_image_coords()`, which returns
  `circle_info`. Write its calibration (`cx_px`, `cy_px`, `scale_um_per_px`,
  `image_width_px`, `image_height_px`) — and `ldir_scan_diameter_um` for the
  fallback — to meta.
- **FTIR / Bruker**: nothing to record (particle-extent placement is
  self-contained).

Keep `bg_image` (already recorded) as the copied run-1 image.

### Piece 2 — Viewer places per instrument, reusing native logic
File: `shiny_app/app.R`, the Multi-Run `output$repro_plot` image block.

Replace the single `compute_image_bounds` call with a branch on
`meta$instrument`:

```
raman        -> raman_image_extent_from_config(cfg_from_meta, run1_x, run1_y)
                with um_per_px / compute_image_bounds as the same 2nd/3rd fallbacks
ftir_perkin  -> list(xmin=min(x), xmax=max(x), ymin=min(y), ymax=max(y))   # run 1
ftir_bruker  -> same particle-extent placement
ldir         -> circle-calibration extent from meta (mirror ldir_native_image_info);
                fallback symmetric ±scan_diameter/2
```

Best done by extracting the four native placement bodies into small shared
helpers in `global.R` (e.g. `place_image_ftir_extent()`,
`place_image_raman_config()`, `place_image_ldir_circle()`) and calling the same
helper from both the native tab and the Multi-Run tab — guarantees they can't
drift apart. `raman_image_extent_from_config` already exists and is reusable
as-is.

### Piece 3 — Simplify the Multi-Run controls
- The manual **Width / Height / X-offset / Y-offset / rotate** inputs become a
  **fallback / fine-tune** only (used when metadata is absent, or for a nudge).
- Default behaviour with metadata present = exact placement, no manual entry.

## Coordinate-frame checks to verify during implementation
1. Confirm `raman_image_extent_from_config` returns extent in the *raw* (x_orig)
   frame when fed run-1 coords (centroid shift resolves to 0). If it instead
   assumes normalized coords, either pass `x_orig` for both args or apply the
   zero-shift explicitly.
2. Confirm LDIR run-1 `x_aligned` equals the circle-calibrated joined coords
   (centred at 0) that `ldir_native_image_info` expects.
3. FTIR/Bruker: confirm the native tab plots particles in `x_orig` (raw) so the
   `min/max` extent aligns — the Multi-Run run-1 frame is raw, so they should match.

## Verification
- No automated test possible for placement (needs the running app), so verify
  visually: load each instrument's 3-run set and confirm points sit on image
  features **identically** to that instrument's single-viewer tab (the Raman
  side-by-side the user already used is the acceptance test).
- Add a lightweight unit test only for any pure helper that computes an extent
  from metadata (e.g. the LDIR circle formula and the FTIR min/max helper).

## Effort
Small–moderate and contained: one metadata block per instrument in the tool,
one placement branch (ideally via 3 shared helpers) in the viewer, and a UI
demotion of the manual controls to fallback. No pipeline changes. No new deps.

## Open follow-up (not blocking)
- LDIR: each run has its own circle calibration; using run 1's is consistent with
  decision 2, but if run-to-run circle calibration differs materially, the
  registration already absorbs it — worth a sanity check on real data.
