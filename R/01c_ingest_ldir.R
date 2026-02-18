# =============================================================================
# 01c_ingest_ldir.R — LDIR data ingestion + image-based coordinate extraction
# =============================================================================
#
# The Agilent 8700 LDIR exports particle data (size, shape, material) in an
# Excel file with two sheets:
#   - "Info"      — sample name, instrument serial, etc.
#   - "Particles" — one row per particle with morphological + spectral data
#
# Critically, LDIR does NOT export particle X/Y centroid coordinates.
# This module:
#   1. Ingests the tabular data (particle sizes, materials, quality)
#   2. Extracts coordinates from the companion LDIR PNG image
#   3. Joins image-derived coordinates with Excel particles by size matching
# =============================================================================


#' Ingest LDIR particle data
#'
#' Reads the Agilent 8700 LDIR export and standardizes columns.
#'
#' @param filepath Path to the LDIR Excel file (.xlsx)
#' @param sheet    Sheet name to read (default "Particles")
#' @return Data frame with standardized columns:
#'   particle_id, x_um, y_um (NA), area_um2, major_um, minor_um,
#'   feret_min_um, feret_max_um, material, quality, source_file
ingest_ldir <- function(filepath, sheet = "Particles") {
  log_message("Reading LDIR data from: ", filepath)

  raw <- readxl::read_excel(filepath, sheet = sheet)
  log_message("  Raw LDIR data: ", nrow(raw), " rows, ", ncol(raw), " columns")
  log_message("  Columns: ", paste(names(raw), collapse = ", "))

  # Try to read sample info from "Info" sheet
  info_sheet <- tryCatch({
    readxl::read_excel(filepath, sheet = "Info")
  }, error = function(e) NULL)

  if (!is.null(info_sheet)) {
    sample_name <- tryCatch({
      info_sheet[[2]][info_sheet[[1]] == "Sample Name"]
    }, error = function(e) NA_character_)
    instrument <- tryCatch({
      info_sheet[[2]][info_sheet[[1]] == "Instrument Name"]
    }, error = function(e) NA_character_)
    log_message("  Sample: ", sample_name)
    log_message("  Instrument: ", instrument)
  }

  # Find columns by fuzzy matching
  id_col     <- find_column(raw, c("Id", "Identifier", "#"))
  width_col  <- find_column(raw, c("Width"))
  height_col <- find_column(raw, c("Height"))
  diam_col   <- find_column(raw, c("Diameter"))
  area_col   <- find_column(raw, c("Area"))
  perim_col  <- find_column(raw, c("Perimeter"))
  mat_col    <- find_column(raw, c("Identification", "Material"))
  qual_col   <- find_column(raw, c("Quality"))
  valid_col  <- find_column(raw, c("Is Valid"))
  aspect_col <- find_column(raw, c("Aspect Ratio"))

  # Particle IDs
  particle_ids <- if (!is.null(id_col)) {
    as.character(raw[[id_col]])
  } else {
    paste0("LDIR_", seq_len(nrow(raw)))
  }

  # Sizes
  width_um  <- safe_col_numeric(raw, width_col)
  height_um <- safe_col_numeric(raw, height_col)
  diam_um   <- safe_col_numeric(raw, diam_col)

  # For feret max/min: use width and height as proxies
  feret_max <- pmax(width_um, height_um, na.rm = TRUE)
  feret_min <- pmin(width_um, height_um, na.rm = TRUE)

  # Material
  material <- if (!is.null(mat_col)) trimws(as.character(raw[[mat_col]])) else NA_character_

  # Quality score
  quality <- safe_col_numeric(raw, qual_col)

  # Is Valid filter
  if (!is.null(valid_col)) {
    valid <- raw[[valid_col]]
    n_invalid <- sum(tolower(valid) != "true", na.rm = TRUE)
    if (n_invalid > 0) {
      log_message("  Flagged ", n_invalid, " invalid particles (keeping all for now)")
    }
  }

  # Build standardized data frame
  df <- data.frame(
    particle_id  = particle_ids,
    x_um         = NA_real_,    # LDIR does not export coordinates
    y_um         = NA_real_,    # will be populated by image extraction
    area_um2     = safe_col_numeric(raw, area_col),
    major_um     = feret_max,
    minor_um     = feret_min,
    feret_min_um = feret_min,
    feret_max_um = feret_max,
    material     = material,
    quality      = quality,
    source_file  = basename(filepath),
    # Keep LDIR-specific columns
    diameter_um  = diam_um,
    aspect_ratio = safe_col_numeric(raw, aspect_col),
    stringsAsFactors = FALSE
  )

  log_message("  Parsed ", nrow(df), " LDIR particles")
  log_message("  Coordinates: NOT available (LDIR does not export X/Y)")
  log_message("  Size range (feret max): ",
              round(min(df$feret_max_um, na.rm = TRUE), 1), " – ",
              round(max(df$feret_max_um, na.rm = TRUE), 1), " µm")
  log_message("  Materials: ", paste(names(head(sort(table(df$material), decreasing = TRUE), 10)),
                                     collapse = ", "))

  df
}


#' Extract LDIR particle coordinates from the companion PNG image
#'
#' The Agilent 8700 LDIR exports a particle map PNG where particles are
#' rendered as colored markers on a background. This function detects those
#' markers and extracts centroids.
#'
#' @param image_path Path to LDIR PNG image file
#' @param scan_bounds Physical scan bounds in µm (list with x_min, x_max, y_min, y_max)
#' @param expected_count Expected number of particles (from Excel data)
#' @return Data frame with particle_id, x_um, y_um, area_um2, etc.
extract_ldir_image_coords <- function(image_path,
                                      scan_bounds = NULL,
                                      expected_count = NULL) {
  log_message("Extracting LDIR particle coordinates from image")

  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' required. Install with: install.packages('png')")
  }

  # Read image
  img <- png::readPNG(image_path)
  h_full <- nrow(img)
  w_full <- ncol(img)
  n_ch <- if (length(dim(img)) == 3) dim(img)[3] else 1
  log_message("  LDIR image: ", w_full, " x ", h_full, " px, ", n_ch, " channels")

  # Determine extraction method based on image characteristics
  if (n_ch >= 3) {
    r <- img[,,1]; g <- img[,,2]; b <- img[,,3]

    # Check for colored markers (high saturation)
    mx <- pmax(r, g, b)
    mn <- pmin(r, g, b)
    sat <- ifelse(mx > 0, (mx - mn) / mx, 0)
    mean_sat <- mean(sat)
    n_colored <- sum(sat > 0.3)
    frac_colored <- n_colored / (h_full * w_full)

    log_message("  Mean saturation: ", round(mean_sat, 3),
                ", colored pixels (sat>0.3): ", round(frac_colored * 100, 1), "%")

    # Strategy decision: if significant colored pixels exist, use
    # color-channel segmentation. Otherwise, fall back to grayscale adaptive.
    if (frac_colored > 0.001 && frac_colored < 0.5) {
      log_message("  Using color-marker extraction strategy")
      result <- .extract_ldir_color_markers(img, h_full, w_full, scan_bounds,
                                             expected_count)
    } else {
      log_message("  Using adaptive-threshold extraction strategy")
      result <- extract_particles_from_image(
        image_path,
        scan_bounds     = scan_bounds,
        adaptive_radius = 50L,
        adaptive_offset = 0.03,
        min_pixels      = 3,
        min_size_um     = 10,
        downsample      = 2L,
        expected_count  = expected_count,
        instrument      = "LDIR"
      )
    }
  } else {
    # Grayscale image — use standard extraction
    result <- extract_particles_from_image(
      image_path,
      scan_bounds     = scan_bounds,
      adaptive_radius = 50L,
      adaptive_offset = 0.03,
      min_pixels      = 3,
      min_size_um     = 10,
      downsample      = 2L,
      expected_count  = expected_count,
      instrument      = "LDIR"
    )
  }

  log_message("  Extracted ", nrow(result), " particles from LDIR image")
  result
}


#' Color-marker extraction for LDIR particle map images
#'
#' Detects colored markers (non-gray pixels) against a uniform background.
#' Works well for Agilent LDIR exports where particles are colored dots.
#'
#' @param img 3D array (h x w x channels)
#' @param h Image height in pixels
#' @param w Image width in pixels
#' @param scan_bounds Physical scan bounds
#' @param expected_count Expected particle count
#' @return Data frame with standard particle columns
.extract_ldir_color_markers <- function(img, h, w, scan_bounds, expected_count) {
  r <- img[,,1]; g <- img[,,2]; b <- img[,,3]

  # Detect colored (saturated) pixels
  mx <- pmax(r, g, b)
  mn <- pmin(r, g, b)
  sat <- ifelse(mx > 0, (mx - mn) / mx, 0)

  # Also detect bright non-background pixels
  gray <- (r + g + b) / 3

  # Background is typically the most common value — find mode
  gray_hist <- hist(gray, breaks = 256, plot = FALSE)
  bg_val <- gray_hist$mids[which.max(gray_hist$counts)]

  # Foreground: either colored OR different from background
  binary <- sat > 0.2 | abs(gray - bg_val) > 0.15

  # Morphological cleanup: remove very small isolated pixels (noise)
  # Simple approach: just rely on min_pixels in CCL

  # Connected components
  cc <- .two_pass_ccl(binary, h, w)
  lab_mat <- cc$labels
  n_components <- cc$n_components

  if (n_components == 0) return(.empty_image_df())

  min_pixels <- 3
  fg_idx    <- which(binary, arr.ind = TRUE)
  fg_labels <- lab_mat[binary]
  tab       <- tabulate(fg_labels, nbins = n_components)
  keep_ids  <- which(tab >= min_pixels)

  if (length(keep_ids) == 0) return(.empty_image_df())

  label_factor <- factor(fg_labels, levels = keep_ids)
  valid <- !is.na(label_factor)
  fg_rows <- fg_idx[valid, 1]
  fg_cols <- fg_idx[valid, 2]
  lf      <- droplevels(label_factor[valid])

  # Centroids and properties
  cy <- tapply(fg_rows, lf, mean)
  cx <- tapply(fg_cols, lf, mean)
  row_min <- tapply(fg_rows, lf, min)
  row_max <- tapply(fg_rows, lf, max)
  col_min <- tapply(fg_cols, lf, min)
  col_max <- tapply(fg_cols, lf, max)
  bb_h <- row_max - row_min + 1
  bb_w <- col_max - col_min + 1
  areas <- tab[keep_ids]

  # Scale to physical coordinates (Cartesian: y increases upward)
  # Image pixels have y=0 at top (row 1), but physical/Raman convention
  # is y increasing upward.  Flip y so row 0 (top) → y_max and
  # row h (bottom) → y_min, matching Raman's Cartesian frame.
  cx_full <- as.numeric(cx) - 0.5
  cy_full <- as.numeric(cy) - 0.5

  if (!is.null(scan_bounds)) {
    x_scale <- (scan_bounds$x_max - scan_bounds$x_min) / w
    y_scale <- (scan_bounds$y_max - scan_bounds$y_min) / h
    x_um <- scan_bounds$x_min + cx_full * x_scale
    y_um <- scan_bounds$y_max - cy_full * y_scale
    size_scale <- max(x_scale, y_scale)
    area_scale <- x_scale * y_scale
  } else {
    x_um <- cx_full
    y_um <- cy_full
    size_scale <- 1
    area_scale <- 1
  }

  n_kept <- length(keep_ids)

  df <- data.frame(
    particle_id  = paste0("LDIR_IMG_", seq_len(n_kept)),
    x_um         = x_um,
    y_um         = y_um,
    area_um2     = as.numeric(areas) * area_scale,
    major_um     = as.numeric(pmax(bb_w, bb_h)) * size_scale,
    minor_um     = as.numeric(pmin(bb_w, bb_h)) * size_scale,
    feret_min_um = as.numeric(pmin(bb_w, bb_h)) * size_scale,
    feret_max_um = as.numeric(pmax(bb_w, bb_h)) * size_scale,
    material     = NA_character_,
    quality      = NA_real_,
    source_file  = "LDIR_image",
    stringsAsFactors = FALSE
  )

  # Auto-trim if way too many components
  if (!is.null(expected_count) && expected_count > 0 && nrow(df) > expected_count * 1.5) {
    target_n <- round(expected_count * 1.3)
    if (nrow(df) > target_n) {
      size_cutoff <- sort(df$feret_max_um, decreasing = TRUE)[min(target_n, nrow(df))]
      df <- df[df$feret_max_um >= size_cutoff, ]
      log_message("  Auto-trimmed to ", nrow(df), " particles (kept largest)")
    }
  }

  df
}


#' Join image-extracted coordinates with Excel particle table
#'
#' Since the LDIR Excel has IDs + sizes but no coordinates, and the image
#' has coordinates + sizes but no IDs, we use size-based optimal matching
#' (Hungarian algorithm) to pair them.
#'
#' @param excel_df Data frame from ingest_ldir() (has particle_id, sizes, material)
#' @param image_df Data frame from extract_ldir_image_coords() (has x_um, y_um, sizes)
#' @return excel_df with x_um and y_um populated, plus coord_match_cost column
join_ldir_coords <- function(excel_df, image_df) {
  log_message("Joining LDIR image coordinates with Excel particle table")

  n_excel <- nrow(excel_df)
  n_image <- nrow(image_df)
  log_message("  Excel: ", n_excel, " particles, Image: ", n_image, " particles")

  if (n_image == 0) {
    log_message("  No image particles to join — coordinates remain NA", level = "WARN")
    excel_df$coord_match_cost <- NA_real_
    excel_df$coord_source <- "none"
    return(excel_df)
  }

  # Build cost matrix using morphological similarity
  n_max <- max(n_excel, n_image)

  # Pre-compute size features
  excel_area  <- excel_df$area_um2
  image_area  <- image_df$area_um2
  excel_feret <- excel_df$feret_max_um
  image_feret <- image_df$feret_max_um

  # Vectorized cost computation (sparse — only fill relevant pairs)
  BIG <- 1e9
  cost <- matrix(BIG, nrow = n_max, ncol = n_max)

  for (i in seq_len(n_excel)) {
    for (j in seq_len(n_image)) {
      c_area <- 0
      c_feret <- 0

      if (!is.na(excel_area[i]) && !is.na(image_area[j]) &&
          excel_area[i] > 0 && image_area[j] > 0) {
        c_area <- abs(log(excel_area[i] / image_area[j]))
      }
      if (!is.na(excel_feret[i]) && !is.na(image_feret[j]) &&
          excel_feret[i] > 0 && image_feret[j] > 0) {
        c_feret <- abs(log(excel_feret[i] / image_feret[j]))
      }

      cost[i, j] <- c_area + 0.5 * c_feret
    }
  }

  # Hungarian assignment
  if (requireNamespace("clue", quietly = TRUE)) {
    assignment <- as.integer(clue::solve_LSAP(cost, maximum = FALSE))
  } else {
    # Greedy fallback: for each excel row, pick best remaining image
    assignment <- rep(NA_integer_, n_excel)
    taken <- logical(n_image)
    for (i in order(apply(cost[seq_len(n_excel), seq_len(n_image), drop = FALSE], 1, min))) {
      best_j <- which.min(cost[i, seq_len(n_image)] + ifelse(taken, BIG, 0))
      if (cost[i, best_j] < BIG && !taken[best_j]) {
        assignment[i] <- best_j
        taken[best_j] <- TRUE
      }
    }
  }

  # Apply coordinates (reject matches above cost threshold)
  max_cost <- 2.0
  excel_df$coord_match_cost <- NA_real_
  excel_df$coord_source <- "none"
  n_joined <- 0
  n_rejected <- 0

  for (i in seq_len(n_excel)) {
    j <- assignment[i]
    if (is.na(j) || j > n_image || cost[i, j] >= BIG) next

    if (cost[i, j] > max_cost) {
      n_rejected <- n_rejected + 1
      next
    }

    excel_df$x_um[i] <- image_df$x_um[j]
    excel_df$y_um[i] <- image_df$y_um[j]
    excel_df$coord_match_cost[i] <- cost[i, j]
    excel_df$coord_source[i] <- "image"
    n_joined <- n_joined + 1
  }

  n_missing <- n_excel - n_joined
  log_message("  Joined coordinates: ", n_joined, " of ", n_excel, " particles")
  if (n_rejected > 0) {
    log_message("  Rejected ", n_rejected, " joins with cost > ", max_cost,
                " (poor size match)")
  }
  if (n_missing > 0) {
    log_message("  Missing coordinates: ", n_missing, " particles (no image match)")
  }

  # Quality diagnostics
  costs <- excel_df$coord_match_cost[!is.na(excel_df$coord_match_cost)]
  if (length(costs) > 0) {
    log_message("  Join quality — cost median: ", round(median(costs), 3),
                ", mean: ", round(mean(costs), 3),
                ", max: ", round(max(costs), 3))
  }

  excel_df
}


#' Validate LDIR coordinate join by checking scan-order correlation
#'
#' If LDIR IDs follow a raster scan order, the ID sequence should correlate
#' with the spatial position of joined coordinates.
#'
#' @param df Data frame with particle_id and x_um, y_um (joined)
#' @return List with tau correlation coefficient and assessment
validate_ldir_scan_order <- function(df) {
  # Only use particles with valid coordinates
  valid <- !is.na(df$x_um) & !is.na(df$y_um)
  if (sum(valid) < 10) {
    return(list(tau = NA_real_, scan_order_consistent = NA,
                message = "Too few coordinated particles to test"))
  }

  df_valid <- df[valid, ]

  # Extract numeric part of ID for ordering
  id_num <- as.numeric(gsub("[^0-9]", "", df_valid$particle_id))
  if (all(is.na(id_num))) {
    return(list(tau = NA_real_, scan_order_consistent = NA,
                message = "Non-numeric IDs — cannot test scan order"))
  }

  # Test raster-scan correlation: sort by (y descending, x ascending)
  spatial_rank <- rank(-df_valid$y_um * 1e6 + df_valid$x_um)
  id_rank <- rank(id_num, na.last = "keep")

  tau <- cor(spatial_rank, id_rank, method = "kendall", use = "complete.obs")

  list(
    tau = tau,
    scan_order_consistent = !is.na(tau) && abs(tau) > 0.5,
    message = if (!is.na(tau)) {
      paste0("Kendall tau = ", round(tau, 3),
             if (abs(tau) > 0.7) " (strong)" else if (abs(tau) > 0.5) " (moderate)" else " (weak)")
    } else "Could not compute"
  )
}


#' Pre-filter LDIR particle data
#'
#' @param df Data frame from ingest_ldir()
#' @param min_quality Minimum quality score (0 = keep all)
#' @param min_size_um Minimum particle size in µm (applied to feret_max_um)
#' @param remove_invalid Logical, remove "Is Valid" = FALSE particles
#' @return Filtered data frame
prefilter_ldir <- function(df, min_quality = 0, min_size_um = 0,
                           remove_invalid = FALSE) {
  n_start <- nrow(df)
  log_message("Pre-filtering LDIR: ", n_start, " particles")

  # Quality filter
  if (min_quality > 0 && any(!is.na(df$quality))) {
    keep <- is.na(df$quality) | df$quality >= min_quality
    df <- df[keep, ]
    log_message("  Quality filter (>= ", min_quality, "): kept ", nrow(df),
                " of ", n_start, " particles")
  }

  # Size filter
  if (min_size_um > 0 && any(!is.na(df$feret_max_um))) {
    keep <- is.na(df$feret_max_um) | df$feret_max_um >= min_size_um
    df <- df[keep, ]
    log_message("  Size filter (>= ", min_size_um, " µm): kept ", nrow(df), " particles")
  }

  log_message("  LDIR after filtering: ", nrow(df), " particles")
  df
}
