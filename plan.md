# Plan: Fix LDIR↔Raman 90° rotation mismatch

## Files changed
| File | Change |
|---|---|
| `R/00_config.R` | Add `ldir_rotate_deg_for_alignment = -90` |
| `R/utils.R` | Add `rotate_coords_90()` helper |
| `R/03c_procrustes_align.R` | Add `rotate_deg` param to `normalize_coords_ldir()` |
| `main.R` | Pass `rotate_deg` to `normalize_coords_ldir()` call |

---

## Step 0 — `R/00_config.R`

Insert after the `ldir_flip_y_for_alignment` line (line 58):

```r
    ldir_rotate_deg_for_alignment = -90,  # clockwise rotation before alignment;
                                           # corrects instrument export convention vs Raman
                                           # must be one of: 0, 90, -90, 180
```

---

## Step 1 — `R/utils.R`: add `rotate_coords_90()`

Insert a small helper near the other coordinate transform utilities (after `extract_transform_params()`, ~line 249):

```r
#' Apply a multiple-of-90° rotation to centred (x, y) coordinate vectors
#'
#' Only exact multiples of 90 are accepted. Positive = counter-clockwise,
#' negative = clockwise (standard mathematical convention with y-up).
#'
#' @param x,y Numeric vectors (centred, same length).
#' @param deg Integer. Must be in {0, 90, -90, 180}.
#' @return List with rotated x and y vectors.
rotate_coords_90 <- function(x, y, deg) {
  deg <- as.integer(round(deg)) %% 360L
  if (deg < 0L) deg <- deg + 360L
  switch(as.character(deg),
    "0"   = list(x = x,  y = y),
    "90"  = list(x = -y, y = x),   # counter-clockwise
    "180" = list(x = -x, y = -y),
    "270" = list(x = y,  y = -x),  # same as -90 clockwise
    stop("rotate_coords_90: deg must be a multiple of 90, got: ", deg)
  )
}
```

Note: the modulo arithmetic normalises −90 → 270, so `switch("270")` is the -90/CW case — mathematically equivalent.

---

## Step 2 — `R/03c_procrustes_align.R`: update `normalize_coords_ldir()`

### Signature change

```r
normalize_coords_ldir <- function(df, flip_y = TRUE, scale_coords = FALSE,
                                   rotate_deg = 0, debug_dir = NULL) {
```

### Early-exit guard (fewer than 2 valid points) — update norm_params to include rotate_deg

```r
    return(list(df = df,
                norm_params = list(centroid_x = 0, centroid_y = 0,
                                   scale_factor = 1, y_flip_applied = flip_y,
                                   rotate_deg_applied = rotate_deg)))
```

### After centering, before flip_y — insert rotation block

Replace:
```r
  df$x_norm <- x_c / sf
  df$y_norm  <- if (flip_y) -(y_c / sf) else (y_c / sf)
```

With:
```r
  # Validate rotate_deg
  if (!rotate_deg %in% c(0L, 90L, -90L, 180L)) {
    stop("normalize_coords_ldir: rotate_deg must be one of {0, 90, -90, 180}, got: ",
         rotate_deg)
  }

  # Apply rotation (after centering, before flip_y, before scale)
  rot <- rotate_coords_90(x_c, y_c, rotate_deg)
  x_r <- rot$x
  y_r <- rot$y

  df$x_norm <- x_r / sf
  df$y_norm  <- if (flip_y) -(y_r / sf) else (y_r / sf)
```

### Update `norm_params` to include `rotate_deg_applied`

```r
  norm_params <- list(
    centroid_x         = cx,
    centroid_y         = cy,
    scale_factor       = sf,
    y_flip_applied     = flip_y,
    rotate_deg_applied = rotate_deg,
    n_valid            = sum(valid),
    x_norm_range       = range(df$x_norm[valid]),
    y_norm_range       = range(df$y_norm[valid])
  )
```

### Update log message

```r
  log_message("  LDIR normalize: centroid=(", round(cx, 1), ", ", round(cy, 1), ")",
              ", scale=", round(sf, 4),
              ", y_flip=", flip_y,
              ", rotate_deg=", rotate_deg,
              ", n=", sum(valid))
```

### Update JSON output (add `rotate_deg_applied` before `n_valid`)

```r
      json_lines <- c(
        "{",
        paste0('  "centroid_x": ',         round(cx, 6),    ","),
        paste0('  "centroid_y": ',         round(cy, 6),    ","),
        paste0('  "scale_factor": ',       round(sf, 6),    ","),
        paste0('  "y_flip_applied": ',     tolower(as.character(flip_y)), ","),
        paste0('  "rotate_deg_applied": ', rotate_deg,      ","),
        paste0('  "n_valid": ',            sum(valid),      ","),
        paste0('  "x_norm_min": ',         round(min(df$x_norm[valid]), 2), ","),
        paste0('  "x_norm_max": ',         round(max(df$x_norm[valid]), 2), ","),
        paste0('  "y_norm_min": ',         round(min(df$y_norm[valid]), 2), ","),
        paste0('  "y_norm_max": ',         round(max(df$y_norm[valid]), 2)),
        "}"
      )
```

---

## Step 3 — `main.R`: wire the new argument

In the `normalize_coords_ldir()` call (line ~615), add one argument:

```r
    ldir_norm_result <- normalize_coords_ldir(
      ldir_with_coords,
      flip_y       = isTRUE(config$ldir_flip_y_for_alignment),
      scale_coords = isTRUE(config$normalize_scale),
      rotate_deg   = config$ldir_rotate_deg_for_alignment %||% 0,
      debug_dir    = if (isTRUE(config$debug)) config$debug_dir else NULL
    )
```

---

## What is NOT changed

- Circle detection (untouched)
- `map_pixels_to_um_circle()` (untouched)
- Raw `x_um` / `y_um` columns — rotation only applies to `x_norm` / `y_norm`
- LDIR image display / Shiny viewer (raw coords still correct there)
- `icp_max_rotation_deg` guardrail in config (already 90°, no change needed)

---

## Rotation math reference (y-up convention)

| `deg` argument | normalised to | formula | description |
|---|---|---|---|
| `0` | 0° | `(x, y)` | no-op |
| `90` | 90° | `(-y, x)` | CCW |
| `-90` | 270° | `(y, -x)` | CW (the fix) |
| `180` | 180° | `(-x, -y)` | half-turn |
