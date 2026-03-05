# =============================================================================
# 03c_procrustes_align.R — Explicit normalization + SVD-Procrustes alignment
# =============================================================================
#
# Provides:
#   normalize_coords_ldir()        — explicit centering/Y-flip, returns norm_params
#   fit_similarity_from_landmarks() — SVD Procrustes on named correspondences
#
# This is the "hard-constraint" alternative to RANSAC + ICP when the user
# can supply explicit LDIR→Raman particle ID correspondences via
# config$ldir_landmark_map.  The SVD fit guarantees landmark residuals are
# minimized exactly, so A3 cannot drift due to ICP re-weighting.
# =============================================================================


#' Normalize LDIR coordinates with an explicit, auditable norm_params object
#'
#' Replaces the ad-hoc inline normalization in main.R so that the exact
#' centroid, scale factor, and Y-flip applied are captured and written to
#' a JSON file for debugging.
#'
#' @param df Data frame with x_um, y_um (must have at least some non-NA rows)
#' @param flip_y Logical. If TRUE negate y after centering (undo image Y-flip).
#' @param scale_coords Logical. If TRUE divide by RMS distance from centroid.
#' @param rotate_deg Integer. In-plane rotation applied after centering, before
#'   flip_y. Must be one of {0, 90, -90, 180}. Use -90 to correct the
#'   LDIR instrument export convention (90° CW mismatch vs Raman).
#' @param debug_dir Optional path. If set, writes ldir_norm_params.json there.
#' @return List with:
#'   df          — original df with x_norm, y_norm columns added/replaced
#'   norm_params — list(centroid_x, centroid_y, scale_factor, y_flip_applied,
#'                      rotate_deg_applied)
normalize_coords_ldir <- function(df, flip_y = TRUE, scale_coords = FALSE,
                                   rotate_deg = 0, debug_dir = NULL) {

  valid <- !is.na(df$x_um) & !is.na(df$y_um)
  if (sum(valid) < 2) {
    warning("normalize_coords_ldir: fewer than 2 valid particles")
    df$x_norm <- df$x_um
    df$y_norm <- df$y_um
    return(list(df = df,
                norm_params = list(centroid_x = 0, centroid_y = 0,
                                   scale_factor = 1, y_flip_applied = flip_y,
                                   rotate_deg_applied = rotate_deg)))
  }

  # Guardrail: only exact multiples of 90 are permitted
  if (!rotate_deg %in% c(0L, 90L, -90L, 180L)) {
    stop("normalize_coords_ldir: rotate_deg must be one of {0, 90, -90, 180}, got: ",
         rotate_deg)
  }

  cx <- mean(df$x_um[valid])
  cy <- mean(df$y_um[valid])

  x_c <- df$x_um - cx
  y_c <- df$y_um - cy

  # Apply rotation after centering, before flip_y and scale
  rot <- rotate_coords_90(x_c, y_c, rotate_deg)
  x_r <- rot$x
  y_r <- rot$y

  # Scale
  sf <- 1.0
  if (scale_coords) {
    rms <- sqrt(mean(x_r[valid]^2 + y_r[valid]^2))
    if (rms > 0) sf <- rms
  }

  df$x_norm <- x_r / sf
  df$y_norm  <- if (flip_y) -(y_r / sf) else (y_r / sf)

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

  log_message("  LDIR normalize: centroid=(", round(cx, 1), ", ", round(cy, 1), ")",
              ", scale=", round(sf, 4),
              ", y_flip=", flip_y,
              ", rotate_deg=", rotate_deg,
              ", n=", sum(valid))

  # Write JSON if debug enabled
  if (!is.null(debug_dir) && dir.exists(debug_dir)) {
    json_path <- file.path(debug_dir, "ldir_norm_params.json")
    tryCatch({
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
      writeLines(json_lines, json_path)
      stopifnot(file.exists(json_path))
      log_message("  Wrote ldir_norm_params.json")
    }, error = function(e) {
      log_message("  WARN: could not write ldir_norm_params.json: ", e$message,
                  level = "WARN")
    })
  }

  list(df = df, norm_params = norm_params)
}


#' Fit a 2D similarity transform from named landmark correspondences (SVD)
#'
#' Given explicit LDIR→Raman particle ID pairs, extract the paired coordinates
#' and solve the optimal similarity transform (scale + rotation + translation)
#' via SVD.  This is an exact, closed-form solution — no iteration, no
#' convergence issues.
#'
#' The transform maps LDIR normalized space → Raman normalized space so that
#' the landmark residuals are globally minimized.
#'
#' @param src_df LDIR data frame with x_norm, y_norm, particle_id
#' @param tgt_df Raman data frame with x_norm, y_norm, particle_id
#' @param landmark_map Named character vector. Names = LDIR particle_id,
#'   values = Raman particle_id.
#'   Example: c("A3" = "A3", "MP_11" = "Raman_190")
#' @param debug_dir Optional path to debug directory. If set, saves
#'   landmark_pairs.csv, landmark_residuals_pre.csv, transform_similarity.json
#' @return List with:
#'   success          — logical
#'   matrix           — 3x3 homogeneous transform (compatible with apply_transform_points)
#'   params           — list(scale, rotation_deg, tx, ty, reflected)
#'   n_pairs          — number of valid landmark pairs used
#'   residual_rms     — RMS landmark residual after transform (µm in norm space)
#'   landmark_pairs   — data frame with per-landmark residuals
#'   message          — human-readable summary
fit_similarity_from_landmarks <- function(src_df, tgt_df, landmark_map,
                                           debug_dir = NULL) {

  if (is.null(landmark_map) || length(landmark_map) == 0) {
    return(list(success = FALSE, matrix = NULL, params = NULL,
                n_pairs = 0, residual_rms = NA_real_,
                landmark_pairs = data.frame(),
                message = "No landmark map provided"))
  }

  ldir_ids  <- names(landmark_map)
  raman_ids <- unname(landmark_map)

  # Resolve rows
  src_rows <- match(ldir_ids,  src_df$particle_id)
  tgt_rows <- match(raman_ids, tgt_df$particle_id)

  ok_mask <- !is.na(src_rows) & !is.na(tgt_rows)

  # Also require non-NA coordinates
  ok_src_coord <- ok_mask
  ok_tgt_coord <- ok_mask
  for (k in seq_along(src_rows)) {
    if (!ok_mask[k]) next
    if (is.na(src_df$x_norm[src_rows[k]]) || is.na(src_df$y_norm[src_rows[k]])) {
      ok_src_coord[k] <- FALSE
    }
    if (is.na(tgt_df$x_norm[tgt_rows[k]]) || is.na(tgt_df$y_norm[tgt_rows[k]])) {
      ok_tgt_coord[k] <- FALSE
    }
  }
  ok_mask <- ok_src_coord & ok_tgt_coord

  n_valid <- sum(ok_mask)

  log_message("  Procrustes landmark map: ", length(ldir_ids), " requested, ",
              n_valid, " resolved")

  # Log unresolved IDs for diagnosis
  for (k in seq_along(ldir_ids)) {
    status <- if (!ok_mask[k]) {
      src_found <- !is.na(src_rows[k])
      tgt_found <- !is.na(tgt_rows[k])
      if (!src_found && !tgt_found) "LDIR+Raman ID not found"
      else if (!src_found) paste0("LDIR '", ldir_ids[k], "' not found")
      else if (!tgt_found) paste0("Raman '", raman_ids[k], "' not found")
      else "coord NA"
    } else "OK"
    if (status != "OK") {
      log_message("    Landmark '", ldir_ids[k], "' → '", raman_ids[k],
                  "': ", status, level = "WARN")
    }
  }

  if (n_valid < 2) {
    msg <- paste0("Only ", n_valid, " valid landmark pairs (need >= 2)")
    log_message("  Procrustes: ", msg, level = "WARN")
    return(list(success = FALSE, matrix = NULL, params = NULL,
                n_pairs = n_valid, residual_rms = NA_real_,
                landmark_pairs = data.frame(), message = msg))
  }

  # Extract coordinate pairs
  valid_idx <- which(ok_mask)
  X_src <- cbind(src_df$x_norm[src_rows[valid_idx]],
                 src_df$y_norm[src_rows[valid_idx]])  # n×2
  Y_tgt <- cbind(tgt_df$x_norm[tgt_rows[valid_idx]],
                 tgt_df$y_norm[tgt_rows[valid_idx]])  # n×2

  # --- SVD 2D similarity ---
  # min_{s,R,t} || s R X + t·1^T - Y ||^2_F
  mu_X <- colMeans(X_src)
  mu_Y <- colMeans(Y_tgt)
  Xc   <- sweep(X_src, 2, mu_X)
  Yc   <- sweep(Y_tgt, 2, mu_Y)

  M     <- t(Yc) %*% Xc      # 2×2
  svd_M <- svd(M)

  # Handle reflection: ensure det(R) = +1 unless reflection is needed
  det_UV  <- det(svd_M$u %*% t(svd_M$v))
  det_sgn <- if (det_UV >= 0) 1 else -1
  D_diag  <- c(1, det_sgn)   # flip last singular vector if needed

  R_2d <- svd_M$u %*% diag(D_diag) %*% t(svd_M$v)  # 2×2 rotation (or improper)

  # Scale: s = trace(D S) / trace(Xc^T Xc)
  denom <- sum(Xc^2)
  s     <- if (denom > 0) sum(svd_M$d * D_diag) / denom else 1.0

  # Translation (in norm space)
  t_vec <- mu_Y - s * (R_2d %*% mu_X)

  # Build 3×3 homogeneous matrix compatible with apply_transform_points:
  #   result = M %*% [x; y; 1]
  #
  # Columns filled column-major in R:
  #   col1=[s*R[1,1], s*R[2,1], 0],  col2=[s*R[1,2], s*R[2,2], 0],  col3=[tx, ty, 1]
  sR <- s * R_2d
  M_3x3 <- matrix(
    c(sR[1,1], sR[2,1], 0,
      sR[1,2], sR[2,2], 0,
      t_vec[1], t_vec[2], 1),
    nrow = 3, byrow = FALSE
  )

  # --- Compute per-landmark residuals ---
  transformed_lm <- apply_transform_points(X_src[,1], X_src[,2], M_3x3)
  dx <- transformed_lm$x_transformed - Y_tgt[,1]
  dy <- transformed_lm$y_transformed - Y_tgt[,2]
  resid <- sqrt(dx^2 + dy^2)
  rms   <- sqrt(mean(resid^2))

  reflected  <- det_sgn < 0
  angle_deg  <- atan2(R_2d[2,1], R_2d[1,1]) * 180 / pi

  params <- list(
    scale        = s,
    rotation_deg = angle_deg,
    tx           = t_vec[1],
    ty           = t_vec[2],
    reflected    = reflected
  )

  landmark_pairs_df <- data.frame(
    ldir_id      = ldir_ids[valid_idx],
    raman_id     = raman_ids[valid_idx],
    ldir_x_norm  = X_src[,1],
    ldir_y_norm  = X_src[,2],
    raman_x_norm = Y_tgt[,1],
    raman_y_norm = Y_tgt[,2],
    aligned_x    = transformed_lm$x_transformed,
    aligned_y    = transformed_lm$y_transformed,
    residual_um  = resid,
    stringsAsFactors = FALSE
  )

  msg <- paste0("SVD Procrustes: n=", n_valid,
                ", scale=", round(s, 4),
                ", rot=", round(angle_deg, 2), "deg",
                ", reflected=", reflected,
                ", RMS=", round(rms, 2), " µm")
  log_message("  ", msg)

  # --- Debug artifacts ---
  if (!is.null(debug_dir) && dir.exists(debug_dir)) {
    tryCatch({
      # landmark_pairs.csv
      lp_path <- file.path(debug_dir, "landmark_pairs.csv")
      write.csv(landmark_pairs_df, lp_path, row.names = FALSE)
      stopifnot(file.exists(lp_path))

      # landmark_residuals_pre.csv
      lr_path <- file.path(debug_dir, "landmark_residuals_pre.csv")
      write.csv(landmark_pairs_df[, c("ldir_id","raman_id","residual_um",
                                       "aligned_x","aligned_y",
                                       "raman_x_norm","raman_y_norm")],
                lr_path, row.names = FALSE)
      stopifnot(file.exists(lr_path))

      # transform_similarity.json
      tf_path <- file.path(debug_dir, "transform_similarity.json")
      tf_lines <- c(
        "{",
        paste0('  "n_pairs": ',        n_valid,             ","),
        paste0('  "scale": ',          round(s, 6),         ","),
        paste0('  "rotation_deg": ',   round(angle_deg, 4), ","),
        paste0('  "tx": ',             round(t_vec[1], 4),  ","),
        paste0('  "ty": ',             round(t_vec[2], 4),  ","),
        paste0('  "reflected": ',      tolower(as.character(reflected)), ","),
        paste0('  "residual_rms": ',   round(rms, 4)),
        "}"
      )
      writeLines(tf_lines, tf_path)
      stopifnot(file.exists(tf_path))

      log_message("  Procrustes debug artifacts saved to ", debug_dir)
    }, error = function(e) {
      log_message("  WARN: Procrustes debug save failed: ", e$message, level = "WARN")
    })
  }

  list(
    success        = TRUE,
    matrix         = M_3x3,
    params         = params,
    n_pairs        = n_valid,
    residual_rms   = rms,
    landmark_pairs = landmark_pairs_df,
    message        = msg
  )
}
