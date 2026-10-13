# =============================================================================
# 01b_ingest_image.R — Connected-component helpers for image segmentation
# =============================================================================
#
# Shared by the LDIR processed-image path (R/01c_ingest_ldir.R). The former
# generic R particle extractor (extract_particles_from_image) was only reachable
# from the LDIR R fallback, which was removed: LDIR image analysis requires the
# Python backend (inst/python/particle_detector.py).
# =============================================================================


#' Two-pass connected component labeling (4-connected)
#'
#' Efficient O(n) algorithm using union-find. Works on binary matrix.
#' @param binary Logical matrix (TRUE = foreground)
#' @param h Number of rows
#' @param w Number of columns
#' @return List with labels (integer matrix) and n_components (int)
.two_pass_ccl <- function(binary, h, w) {
  lab <- matrix(0L, nrow = h, ncol = w)

  # Union-find data structure (dynamic, grows as needed)
  uf_parent <- integer(0)
  next_label <- 1L

  # --- Pass 1: assign provisional labels ---
  for (r in seq_len(h)) {
    for (cc in seq_len(w)) {
      if (!binary[r, cc]) next

      up   <- if (r > 1 && binary[r - 1, cc]) lab[r - 1, cc] else 0L
      left <- if (cc > 1 && binary[r, cc - 1]) lab[r, cc - 1] else 0L

      if (up == 0L && left == 0L) {
        uf_parent <- c(uf_parent, next_label)
        lab[r, cc] <- next_label
        next_label <- next_label + 1L
      } else if (up != 0L && left == 0L) {
        lab[r, cc] <- up
      } else if (up == 0L && left != 0L) {
        lab[r, cc] <- left
      } else {
        lab[r, cc] <- min(up, left)
        if (up != left) {
          # Union
          ra <- up; while (uf_parent[ra] != ra) ra <- uf_parent[ra]
          rb <- left; while (uf_parent[rb] != rb) rb <- uf_parent[rb]
          if (ra != rb) uf_parent[max(ra, rb)] <- min(ra, rb)
        }
      }
    }
  }

  if (next_label == 1L) {
    return(list(labels = lab, n_components = 0L))
  }

  # --- Pass 2: resolve equivalences ---
  n_prov <- next_label - 1L
  roots <- integer(n_prov)
  for (i in seq_len(n_prov)) {
    x <- i
    while (uf_parent[x] != x) x <- uf_parent[x]
    roots[i] <- x
  }

  unique_roots <- unique(roots)
  remap <- integer(n_prov)
  remap[unique_roots] <- seq_along(unique_roots)
  final_map <- remap[roots]

  fg <- lab > 0L
  lab[fg] <- final_map[lab[fg]]

  list(labels = lab, n_components = length(unique_roots))
}


#' Return an empty data frame with standard image particle columns
.empty_image_df <- function() {
  data.frame(
    particle_id  = character(),
    x_um         = numeric(),
    y_um         = numeric(),
    area_um2     = numeric(),
    major_um     = numeric(),
    minor_um     = numeric(),
    feret_min_um = numeric(),
    feret_max_um = numeric(),
    material     = character(),
    quality      = numeric(),
    source_file  = character(),
    stringsAsFactors = FALSE
  )
}

