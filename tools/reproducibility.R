# =============================================================================
# reproducibility.R — compare repeat runs of one filter on one instrument
# =============================================================================
# Quantifies intra-instrument reproducibility (and, given the known reference
# plastic, accuracy) across N repeat measurements of the SAME filter on the
# SAME instrument. See R/reproducibility.R for the method.
#
# Usage — edit the CONFIG block below, then:
#   Rscript tools/reproducibility.R
#
# Or pass an instrument and run files on the command line:
#   Rscript tools/reproducibility.R <instrument> <run1> <run2> <run3> [...]
#   (instrument = ftir_perkin | ftir_bruker | raman | ldir)
#   For ldir, pass image files interleaved:  ldir <excel1> <img1> <excel2> <img2> ...
# =============================================================================

# --- locate repo root and load modules --------------------------------------
.script_path <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1) normalizePath(f) else NA_character_
}, error = function(e) NA_character_)
REPO_ROOT <- if (!is.na(.script_path)) dirname(dirname(.script_path)) else getwd()

for (m in c("utils.R", "00_config.R", "01_ingest.R", "01c_ingest_ldir.R",
            "07_match.R", "08b_material_map.R", "reproducibility.R")) {
  source(file.path(REPO_ROOT, "R", m))
}

# =============================================================================
# CONFIG — edit these, or override via command-line arguments
# =============================================================================
CONFIG <- list(
  instrument = "ldir",           # ftir_perkin | ftir_bruker | raman | ldir
  # One entry per run. For ldir, set image = path to the companion image.
  runs = list(
    list(file = "C:\\Users\\moellerkr\\OneDrive - IAEA\\My Documents\\Automatisations\\MPs\\test-data\\PET\\A\\Comparstic2026 PET A Kev 260401.xlsx"
         ,
         image = "C:\\Users\\moellerkr\\OneDrive - IAEA\\My Documents\\Automatisations\\MPs\\test-data\\PET\\A\\Comparstic2026 PET A Kev 260401.png"
         ),
    list(file = "C:\\Users\\moellerkr\\OneDrive - IAEA\\My Documents\\Automatisations\\MPs\\test-data\\PET\\B\\Comparstic2026 PET B Kev 260401.xlsx"
         ,
         image = "C:\\Users\\moellerkr\\OneDrive - IAEA\\My Documents\\Automatisations\\MPs\\test-data\\PET\\B\\Comparstic2026 PET B Kev 260401.png"
    ),
    list(file = "C:\\Users\\moellerkr\\OneDrive - IAEA\\My Documents\\Automatisations\\MPs\\test-data\\PET\\C\\Comparstic2026 PET C Kev 260401.xlsx"
         ,
         image = "C:\\Users\\moellerkr\\OneDrive - IAEA\\My Documents\\Automatisations\\MPs\\test-data\\PET\\C\\Comparstic2026 PET C Kev 260401.png"
    )
  ),
  reference_material = "Polyethylene terephthalate",  # NULL to skip accuracy
  match_gate_um  = 75,     # tight: same-instrument localization is precise
  align_gate_um  = 300,    # ICP correspondence gate (absorbs a slight re-seat)
  output_dir     = "output/reproducibility",

  # Raman background-image placement (WITec metadata) — recorded to the meta so
  # the Multi-Run viewer places the image EXACTLY like the single Raman tab.
  #
  # You normally do NOT need to fill these in: the script auto-resolves the
  # correct per-dataset calibration by finding the pipeline run whose
  # inputs$raman file matches runs[[1]]$file (by MD5) and reading that run's
  # OWN manifest.json config_snapshot — the exact values that were in effect
  # when that scan was processed, not whatever happens to be in R/00_config.R
  # right now. This is what used to go stale: hand-copying "whatever's
  # currently in 00_config.R" silently drifted from the actual calibration
  # once R/00_config.R was edited for a different scan (see
  # docs/multirun_image_placement_plan.md, plan v3).
  #
  # These fields are now only a MANUAL OVERRIDE/fallback, used only when no
  # pipeline run's manifest matches runs[[1]]$file's content. Leave them NULL
  # unless you need to force specific values (e.g. no matching run exists yet).
  raman_image_width_um    = NULL,
  raman_image_height_um   = NULL,
  raman_image_center_x_um = NULL,
  raman_image_center_y_um = NULL,
  # If your dataset has NO WITec extent (the single Raman tab shows the image at
  # the µm-per-pixel scale, not the WITec width), leave the four fields above
  # blank and either upload the original TIFF in the viewer (its DPI is read
  # automatically) or set the µm-per-pixel here:
  raman_um_per_px         = NULL
)

# --- command-line override ---------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 3) {
  CONFIG$instrument <- .args[1]
  rest <- .args[-1]
  if (identical(CONFIG$instrument, "ldir")) {
    if (length(rest) %% 2 != 0)
      stop("ldir needs interleaved <excel> <image> pairs")
    CONFIG$runs <- lapply(seq(1, length(rest), by = 2),
                          function(i) list(file = rest[i], image = rest[i + 1]))
  } else {
    CONFIG$runs <- lapply(rest, function(f) list(file = f, image = NULL))
  }
}

# --- instrument-agnostic particle acquisition -------------------------------
# Returns a standardized particle frame (x_um/y_um/material/feret...) for one
# run, dispatching on instrument. LDIR coordinates come from the image pipeline.
get_run_particles <- function(instrument, file, image = NULL, config) {
  ldir_circle <- NULL
  df <- switch(instrument,
    ftir_perkin = ingest_ftir(file),
    ftir_bruker = ingest_ftir_bruker(file),
    raman       = ingest_raman(file),
    ldir        = {
      raw <- ingest_ldir(file)
      if (is.null(image) || !nzchar(image) || !file.exists(image))
        stop("LDIR run needs a companion image with spatial coordinates: ", file)
      diam   <- config$ldir_scan_diameter_um %||% 13000
      bounds <- list(x_min = 0, x_max = diam, y_min = 0, y_max = diam)
      ext <- extract_ldir_image_coords(image, scan_bounds = bounds,
                                       expected_count = nrow(raw), config = config)
      ldir_circle <<- ext$circle_info            # captured for image placement
      join_ldir_coords(raw, ext$particles, config = config)
    },
    stop("Unknown instrument: ", instrument))
  if (!all(c("x_um", "y_um") %in% names(df)))
    stop("Run ingestion did not yield x_um/y_um for ", instrument, " (", file, ")")
  df <- df[is.finite(df$x_um) & is.finite(df$y_um), , drop = FALSE]
  attr(df, "ldir_circle") <- ldir_circle
  df
}

# --- run ---------------------------------------------------------------------
cfg <- make_config()
if (length(CONFIG$runs) < 2) stop("Need at least 2 runs to compare.")

# Raman calibration: auto-resolve from the pipeline run whose raw Raman file
# matches runs[[1]]$file's content, falling back to the manual CONFIG fields
# above only when no matching run is found. See CONFIG comment above and
# resolve_raman_calibration_from_manifests() in R/utils.R for the rationale.
raman_calibration_source <- "none"
if (identical(CONFIG$instrument, "raman")) {
  run1_file <- CONFIG$runs[[1]]$file
  resolved <- tryCatch(
    resolve_raman_calibration_from_manifests(run1_file, file.path(REPO_ROOT, "output")),
    error = function(e) NULL)
  if (!is.null(resolved)) {
    log_message("Raman calibration: auto-resolved from pipeline run ", resolved$run_id,
                " (manifest timestamp ", resolved$timestamp, ") for ",
                basename(run1_file), " — width=", resolved$width_um,
                " height=", resolved$height_um, " center=(", resolved$center_x_um,
                ", ", resolved$center_y_um, ")")
    CONFIG$raman_image_width_um    <- resolved$width_um
    CONFIG$raman_image_height_um   <- resolved$height_um
    CONFIG$raman_image_center_x_um <- resolved$center_x_um
    CONFIG$raman_image_center_y_um <- resolved$center_y_um
    raman_calibration_source <- paste0("manifest:", resolved$run_id)
  } else if (!is.null(CONFIG$raman_image_width_um)) {
    log_message("Raman calibration: no pipeline run's manifest matches ",
                basename(run1_file), " (by MD5) — using the manually-entered ",
                "CONFIG$raman_image_* values. Verify these are correct for ",
                "THIS scan, not a leftover from a previous dataset.")
    raman_calibration_source <- "manual"
  } else {
    log_message("Raman calibration: no matching pipeline run and no manual ",
                "CONFIG$raman_image_* values — the Multi-Run viewer will fall ",
                "back to a particle-extent fit for the background image.")
  }
}

log_message("Reproducibility: ", CONFIG$instrument, " — ", length(CONFIG$runs), " runs")
runs <- lapply(CONFIG$runs, function(r)
  get_run_particles(CONFIG$instrument, r$file, r$image, cfg))
for (i in seq_along(runs))
  log_message("  run", i, ": ", nrow(runs[[i]]), " particles (", CONFIG$runs[[i]]$file, ")")

ref_fam <- if (!is.null(CONFIG$reference_material))
             classify_family_vec(CONFIG$reference_material) else NULL

res <- run_reproducibility(runs,
                           gate = CONFIG$match_gate_um,
                           align_gate = CONFIG$align_gate_um,
                           reference_family = ref_fam)

# --- write outputs -----------------------------------------------------------
# Each run goes in its own timestamped subfolder so results are never
# overwritten and the Shiny Multi-Run tab can list them side by side.
run_id <- paste0(CONFIG$instrument, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
out <- file.path(CONFIG$output_dir, run_id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

write.csv(res$consensus, file.path(out, "reproducibility_particles.csv"),
          row.names = FALSE)

s <- res$summary
summary_df <- data.frame(
  instrument            = CONFIG$instrument,
  n_runs                = s$n_runs,
  count_mean            = round(s$count_mean, 2),
  count_cv              = round(s$count_cv, 4),
  n_consensus           = s$n_consensus,
  detected_in_all       = s$detected_in_all,
  detected_in_all_frac  = round(s$detected_in_all_frac, 4),
  material_concordance  = round(s$material_concordance, 4),
  accuracy_vs_reference = if (is.na(s$accuracy_vs_reference)) NA else round(s$accuracy_vs_reference, 4),
  median_feret_cv       = round(s$median_feret_cv, 4),
  median_pos_jitter_um  = round(s$median_pos_jitter_um, 2),
  stringsAsFactors = FALSE)
write.csv(summary_df, file.path(out, "reproducibility_summary.csv"), row.names = FALSE)

# Viewer-ready long table (per run × consensus particle, aligned frame).
long <- repro_long_table(res)
write.csv(long, file.path(out, "reproducibility_points.csv"), row.names = FALSE)

# Backdrop image for the Multi-Run view: copy run 1's image (if supplied) into
# the output dir. The viewer places it at the point extent (aspect-preserving),
# so no calibration is replicated here.
bg_image <- NA_character_
img1 <- CONFIG$runs[[1]]$image
if (!is.null(img1) && nzchar(img1) && file.exists(img1)) {
  bg_image <- paste0("background", tools::file_ext(img1) |> (\(e) if (nzchar(e)) paste0(".", e) else ".png")())
  file.copy(img1, file.path(out, bg_image), overwrite = TRUE)
}
# Per-instrument image-placement metadata so the Multi-Run viewer reproduces
# the single-instrument tab's exact image<->coordinate relationship.
# run1_file/run1_md5 are recorded for every instrument so a reproducibility
# output can always be traced back to the exact raw file it was built from
# (tools/diagnose_multirun_placement.R reads these).
meta <- data.frame(instrument = CONFIG$instrument, n_runs = length(runs),
                   bg_image = bg_image,
                   run1_file = basename(CONFIG$runs[[1]]$file),
                   run1_md5  = file_md5(CONFIG$runs[[1]]$file) %||% NA_character_,
                   stringsAsFactors = FALSE)
if (identical(CONFIG$instrument, "raman")) {
  meta$raman_image_width_um    <- CONFIG$raman_image_width_um    %||% NA_real_
  meta$raman_image_height_um   <- CONFIG$raman_image_height_um   %||% NA_real_
  meta$raman_image_center_x_um <- CONFIG$raman_image_center_x_um %||% NA_real_
  meta$raman_image_center_y_um <- CONFIG$raman_image_center_y_um %||% NA_real_
  meta$raman_um_per_px         <- CONFIG$raman_um_per_px         %||% NA_real_
  meta$raman_calibration_source <- raman_calibration_source
}
if (identical(CONFIG$instrument, "ldir")) {
  ci <- attr(runs[[1]], "ldir_circle")
  if (!is.null(ci)) {
    meta$ldir_cx_px           <- ci$cx_px           %||% NA_real_
    meta$ldir_cy_px           <- ci$cy_px           %||% NA_real_
    meta$ldir_scale_um_per_px <- ci$scale_um_per_px %||% NA_real_
    meta$ldir_image_width_px  <- ci$width           %||% NA_real_
    meta$ldir_image_height_px <- ci$height          %||% NA_real_
  }
}
write.csv(meta, file.path(out, "reproducibility_meta.csv"), row.names = FALSE)

plots <- repro_plots(res, out, title_prefix = CONFIG$instrument)

# --- console overview --------------------------------------------------------
log_message(strrep("=", 60))
log_message("Reproducibility summary (", CONFIG$instrument, ")")
log_message("  Runs / counts:        ", paste(s$counts_per_run, collapse = ", "),
            "  (CV ", round(100 * s$count_cv, 1), "%)")
log_message("  Physical particles:   ", s$n_consensus)
log_message("  Detected in all runs: ", s$detected_in_all,
            " (", round(100 * s$detected_in_all_frac, 1), "%)")
log_message("  Material concordance: ", round(100 * s$material_concordance, 1), "%")
if (!is.na(s$accuracy_vs_reference))
  log_message("  Accuracy vs '", CONFIG$reference_material, "': ",
              round(100 * s$accuracy_vs_reference, 1), "%")
log_message("  Median Feret CV:      ", round(100 * s$median_feret_cv, 1), "%")
log_message("  Median jitter:        ", round(s$median_pos_jitter_um, 1), " um")
log_message("  Wrote: reproducibility_particles.csv, reproducibility_summary.csv, ",
            length(plots), " plots -> ", out)
