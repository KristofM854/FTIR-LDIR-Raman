# LDIR<->Raman acceptance-gate classification (viewer).
#
# Regression guard for the "40/40 LDIR matched" miscount: with
# ldir_force_complete_match = TRUE the pipeline pairs EVERY LDIR particle with a
# Raman particle regardless of distance, so the raw matched CSV always looks
# fully matched. The viewer must reclassify forced over-gate pairs (match_distance
# above the acceptance gate) as unmatched so the summary, overlay and tables agree.
#
# The functions under test live in shiny_app/global.R, which normally boots the
# Shiny app (library()/relative source()). We load only the pure helpers by
# neutralising those side effects and stubbing the one external dependency.

.load_global_helpers <- function() {
  gpath <- file.path(REPO_ROOT, "shiny_app", "global.R")
  skip_if_not(file.exists(gpath), "shiny_app/global.R not found")
  env <- new.env(parent = globalenv())
  env$classify_family_vec <- function(x) x
  # Deterministic manifest so the gate resolves to a known value.
  env$load_run_manifest <- function(run_dir)
    list(config_snapshot = list(match_dist_threshold_ldir_um = 250))
  txt <- paste(readLines(gpath, warn = FALSE), collapse = "\n")
  txt <- gsub("library\\([^)]*\\)", "invisible(NULL)", txt)
  txt <- gsub("source\\(file\\.path[^\n]*\\)", "invisible(NULL)", txt)
  eval(parse(text = txt), envir = env)
  env
}

.synthetic_ldir_match <- function(distances) {
  n <- length(distances)
  data.frame(
    ldir_particle_id  = paste0("L", seq_len(n)),
    raman_particle_id = paste0("R", seq_len(n)),
    ldir_x_aligned = seq_len(n) * 100, ldir_y_aligned = seq_len(n) * 100,
    ldir_x_um = seq_len(n), ldir_y_um = seq_len(n),
    ldir_area_um2 = 1, ldir_major_um = 1, ldir_minor_um = 1, ldir_feret_max_um = 1,
    ldir_material = "PET", ldir_quality = 1,
    match_id = seq_len(n), match_score = 1,
    match_distance = distances,
    stringsAsFactors = FALSE
  )
}

test_that("annotate_ldir_gate flags only within-gate pairs as genuine", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(5, 12, 14, 19, 1, 40, 60, 90))
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")

  expect_equal(data$ldir_match_gate_um, 20)
  expect_equal(data$ldir_raman_matched$within_gate,
               c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("build_instrument_dfs reports over-gate forced pairs as unmatched LDIR", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(5, 12, 14, 19, 1, 40, 60, 90))
  data <- env$annotate_ldir_gate(
    list(ldir_raman_matched = m, unmatched_ldir = NULL,
         matched = NULL, unmatched_raman = NULL), "dummy")

  res <- env$build_instrument_dfs(data)
  # 5 within gate -> matched; 3 over gate -> unmatched; all 8 still present.
  expect_equal(sum(res$ldir$match_status == "matched"), 5L)
  expect_equal(sum(res$ldir$match_status == "unmatched"), 3L)
  expect_equal(nrow(res$ldir), 8L)
  # Over-gate rows must not claim a Raman partner.
  expect_false(any(res$ldir$matched_to_raman[res$ldir$match_status == "unmatched"]))
})

test_that("genuine-pair helpers gate the Raman->LDIR flag", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(5, 40))  # R1 genuine, R2 over-gate
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")

  # Raw unmatched_raman CSV schema (x_norm/x_um/feret_max_um), as loaded.
  raman <- data.frame(particle_id = c("R1", "R2", "R3"),
                      x_norm = 0, y_norm = 0, x_um = 0, y_um = 0,
                      area_um2 = 1, major_um = 1, minor_um = 1, feret_max_um = 1,
                      material = "PET", quality = 1, stringsAsFactors = FALSE)
  data$matched <- NULL
  data$unmatched_raman <- raman
  res <- env$build_instrument_dfs(data)
  flag <- res$raman$matched_to_ldir[match(c("R1", "R2", "R3"), res$raman$particle_id)]
  expect_equal(flag, c(TRUE, FALSE, FALSE))
})

test_that("legacy match frames without match_distance stay fully matched", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(5, 40))
  m$match_distance <- NULL  # simulate a pre-gate export
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")
  expect_true(all(data$ldir_raman_matched$within_gate))
})

test_that("regate_ldir retunes within_gate to an explicit gate", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(50, 120, 200, 240, 10, 400, 600, 900))
  d <- list(ldir_raman_matched = m)

  expect_equal(sum(env$regate_ldir(d, 100)$ldir_raman_matched$within_gate), 2L)
  expect_equal(sum(env$regate_ldir(d, 250)$ldir_raman_matched$within_gate), 5L)
  expect_equal(sum(env$regate_ldir(d, 1000)$ldir_raman_matched$within_gate), 8L)

  # Stored gate follows the override.
  expect_equal(env$regate_ldir(d, 500)$ldir_match_gate_um, 500)

  # Invalid gate is a no-op (no within_gate written, data untouched).
  out <- env$regate_ldir(d, NA)
  expect_false("within_gate" %in% names(out$ldir_raman_matched))
})

test_that("ldir_acceptance_gate falls back to the 20 um default when unrecorded", {
  env <- .load_global_helpers()
  # Override the manifest loader to return no gate keys.
  env$load_run_manifest <- function(run_dir) list(config_snapshot = list())
  expect_equal(env$ldir_acceptance_gate("dummy"), 20)
})

# AUDIT M-03: the pipeline's LDIR gate is fitted per run (3 sigma, 10-50 um)
# and recorded in chance_pairs.csv; the viewer re-gated with the configured
# 20 um and would drop a pair the run made at e.g. 20.8 um.
test_that("ldir_acceptance_gate uses the gate the run recorded in chance_pairs.csv", {
  env <- .load_global_helpers()
  env$load_run_manifest <- function(run_dir)
    list(config_snapshot = list(match_dist_threshold_ldir_um = 20))
  run <- withr::local_tempdir()
  expect_equal(env$ldir_acceptance_gate(run), 20)           # no file: the configured gate
  dir.create(file.path(run, "07_diagnostics"))
  utils::write.csv(data.frame(pair = c("ftir_perkin_raman", "ldir_raman"), gate_um = c(50, 21.2)),
                   file.path(run, "07_diagnostics", "chance_pairs.csv"), row.names = FALSE)
  expect_equal(env$ldir_acceptance_gate(run), 21.2)
  m <- .synthetic_ldir_match(c(20.8, 25))
  m$ldir_major_um <- c(30, 30)
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), run)
  expect_equal(data$ldir_raman_matched$within_gate, c(TRUE, FALSE))
  expect_equal(data$ldir_match_gate_um, 21.2)
})

test_that("pairs made by the matcher's extent rule stay matched in the viewer", {
  env <- .load_global_helpers()
  m <- .synthetic_ldir_match(c(113, 113))
  m$ldir_major_um <- c(349, 349); m$match_rule <- c("extent", "distance")
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")
  expect_equal(data$ldir_raman_matched$within_gate, c(TRUE, FALSE))
})

test_that("large LDIR particles keep the matcher's widened gate in the viewer", {
  env <- .load_global_helpers()
  # Dataset A: a 428 um particle matched at 58 um (gate 0.15 x 428 = 64 um)
  m <- .synthetic_ldir_match(c(58, 58))
  m$ldir_major_um <- c(428, 30)
  data <- env$annotate_ldir_gate(list(ldir_raman_matched = m), "dummy")
  expect_equal(data$ldir_raman_matched$within_gate, c(TRUE, FALSE))
})

test_that("link_groups joins particles through any pair, Raman or direct", {
  env <- .load_global_helpers()
  ram <- list(inst = "raman", id = "raman_particle_id", x = "raman_x_norm", y = "raman_y_norm")
  pe  <- list(inst = "ftir_pe", id = "ftir_particle_id", x = "ftir_x_aligned", y = "ftir_y_aligned")
  ld  <- list(inst = "ldir", id = "ldir_particle_id", x = "ldir_x_aligned", y = "ldir_y_aligned")
  br  <- list(inst = "ftir_bruker", id = "ftir_bruker_particle_id", x = "ftir_bruker_x_aligned",
              y = "ftir_bruker_y_aligned")
  ld2 <- list(inst = "ldir", id = "ldir_particle_id", x = "ldir_x_aligned", y = "ldir_y_aligned")
  legs <- list(
    # P_1 pairs with Raman R1 only
    list(df = data.frame(ftir_particle_id = c("P_1", "P_2"), ftir_x_aligned = c(0, 50),
                         ftir_y_aligned = 0, raman_particle_id = c("R1", "R2"),
                         raman_x_norm = c(1, 51), raman_y_norm = 0), a = pe, b = ram),
    # Bruker B1 <-> LDIR L1 and LDIR L1 <-> PE P_1: a four-instrument group
    list(df = data.frame(ftir_bruker_particle_id = "B1", ftir_bruker_x_aligned = 2,
                         ftir_bruker_y_aligned = 0, ldir_particle_id = "L1",
                         ldir_x_aligned = 3, ldir_y_aligned = 0), a = br, b = ld2),
    list(df = data.frame(ldir_particle_id = c("L1", "L9"), ldir_x_aligned = c(3, 900),
                         ldir_y_aligned = 0, ftir_particle_id = c("P_1", "P_9"),
                         ftir_x_aligned = c(0, 901), ftir_y_aligned = 0), a = ld, b = pe))
  g <- env$link_groups(legs)
  n <- setNames(g$n_inst, paste(g$inst, g$id))
  expect_equal(unname(n[c("ftir_pe P_1", "raman R1", "ldir L1", "ftir_bruker B1")]), rep(4L, 4))
  expect_equal(unname(n[c("ftir_pe P_2", "raman R2")]), c(2L, 2L))
  expect_equal(unname(n[c("ldir L9", "ftir_pe P_9")]), c(2L, 2L))      # FTIR + LDIR only
  expect_null(env$link_groups(list()))
})

test_that("overlay symbols: one per group, unpaired particles as their own", {
  env <- .load_global_helpers()
  grp <- data.frame(inst = c("ftir_pe", "raman", "ldir"), id = c("P_1", "R1", "L1"),
                    x = c(0, 2, 4), y = 0, group = 1L, n_inst = 3L, stringsAsFactors = FALSE)
  singles <- list(raman = data.frame(particle_id = c("R1", "R2"), x = c(2, 500), y = 0),
                  ldir = data.frame(particle_id = "L1", x = 4, y = 0))
  st <- env$overlay_symbol_table(grp, singles)
  expect_equal(nrow(st$sym), 2)
  expect_equal(st$sym$combo, c("ftir_pe+raman+ldir", "raman"))
  expect_equal(st$sym$x, c(2, 500))
  g <- env$overlay_symbol_geometry(st$sym, r = 10, n_arc = 4L)
  expect_setequal(unique(g$quarters$inst), c("ftir_pe", "raman", "ldir"))
  # top-left quarter for FTIR (PerkinElmer): x <= centre, y >= centre
  q <- g$quarters[g$quarters$inst == "ftir_pe", ]
  expect_true(all(q$x <= 2 + 1e-9) && all(q$y >= -1e-9))
  expect_null(env$overlay_symbol_table(NULL, list()))
})

test_that("run labels name the input folder after (My) Documents", {
  env <- .load_global_helpers()
  m <- list(inputs = list(
    ftir = list(path = "C:/runs/out/2026-10-28_1/inputs/ftir_stitched/A_stitched.csv"),
    ftir_image = list(path = "C:\\Users\\me\\OneDrive - IAEA\\My Documents\\Automatisations\\MPs\\Field_samples\\Bora2025\\A\\sample_perkin.bmp")))
  expect_equal(env$run_input_folder(m, "C:/runs/out/2026-10-28_1"),
               "Automatisations/MPs/Field_samples/Bora2025/A")
  m2 <- list(inputs = list(raman = list(path = "/data/lab/projects/x/y/z/Raman.csv")))
  expect_equal(env$run_input_folder(m2), "projects/x/y/z")
  expect_equal(env$run_input_folder(list()), "")
})

test_that("instrument colours resolve from keys and display names", {
  env <- .load_global_helpers()
  expect_identical(env$instrument_colour(c("FTIR (PerkinElmer)", "ftir_perkin", "FTIR (Bruker)",
                                           "Raman", "LDIR", "LD-IR", "something")),
                   c("#F0E442", "#F0E442", "#CC79A7", "#0072B2", "#D55E00", "#D55E00", "grey40"))
  b <- as.character(env$instrument_badge("raman"))
  expect_match(b, "background:#0072B2", fixed = TRUE)
  expect_match(b, "color:#FFFFFF", fixed = TRUE)
})

test_that("precision view: one-instrument groups marked and left out of Mandel; plastics switch", {
  env <- .load_global_helpers()
  sys.source(file.path(REPO_ROOT, "R", "08b_material_map.R"), envir = env)
  inst <- c("ftir", "ftir_bruker", "raman", "ldir")
  grp <- c("All particles", "Plastics", "Synthetic", "PC", "Cellulose", "Protein")
  counts <- expand.grid(instrument = inst, replicate = 1:3, group = grp, size_class = "all sizes",
                        stringsAsFactors = FALSE)
  counts$count <- ifelse(counts$group == "Protein", ifelse(counts$instrument == "ftir", 6, 0),
                  ifelse(counts$group == "Cellulose", ifelse(counts$instrument %in% c("ftir", "ldir"), 1, 0),
                         20 + counts$replicate))
  ins <- stats::aggregate(count ~ instrument + group + size_class, counts, mean)
  names(ins)[4] <- "mean"; ins$h <- 0.5; ins$k <- 1
  res <- list(counts = counts, instruments = ins,
              iso = data.frame(group = grp, size_class = "all sizes", mean = 1, s_r = 1, s_r_low = 0.5,
                               s_r_high = 2, s_L = 1, s_R = 1.4, s_R_low = 0.8, s_R_high = 3,
                               r_limit = 2.8, R_limit = 4, cv_r_pct = 1, cv_L_pct = 1, cv_R_pct = 1),
              critical = list(h_5 = 1.42, h_1 = 1.49, k_5 = 1.59, k_1 = 1.77))
  expect_equal(env$precision_single_groups(res), "Protein")
  expect_equal(env$precision_groups(res, mode = "plastics", category = "both"), c("Plastics", "PC", "Cellulose"))
  expect_equal(env$precision_groups(res, mode = "plastics", category = "synthetic"), c("Synthetic", "PC"))
  expect_true("Protein *" %in% env$precision_iso_table(res)$Group)
  h <- env$precision_mandel_plot(res, "h", groups = grp[-3])
  expect_false("Protein" %in% as.character(h$data$group))
  expect_s3_class(env$precision_counts_plotly(res), "plotly")
  ins$dispersion_p <- 0.5; ins$h_flag <- ""; ins$k_flag <- ""
  ins$k_flag[ins$group == "Synthetic" & ins$instrument == "ldir"] <- "5%"
  res$instruments <- ins
  txt <- env$precision_summary_text(res)
  expect_match(txt, "Synthetic particles (all sizes): mean 1.0 per filter over 4 instruments", fixed = TRUE)
  expect_match(txt, "Mandel flags: k for LDIR (5%)", fixed = TRUE)
  expect_match(env$precision_summary_text(res, "Protein"), "Only one instrument")
  expect_match(env$precision_summary_text(res, "PET"), "No PET particles")
})

test_that("precision results written before the Synthetic group get it from the families", {
  env <- .load_global_helpers()
  sys.source(file.path(REPO_ROOT, "R", "08b_material_map.R"), envir = env)
  inst <- c("ftir", "raman")
  counts <- expand.grid(instrument = inst, replicate = 1:3, group = c("All particles", "PC", "PET", "Cellulose"),
                        size_class = c("all sizes", "20-50"), stringsAsFactors = FALSE)
  counts$count <- 5 + seq_len(nrow(counts)) %% 7
  dir <- withr::local_tempdir()
  utils::write.csv(counts, file.path(dir, "precision_counts.csv"), row.names = FALSE)
  utils::write.csv(data.frame(group = "All particles", size_class = "all sizes", mean = 1),
                   file.path(dir, "precision_iso5725.csv"), row.names = FALSE)
  utils::write.csv(data.frame(group = "All particles", size_class = "all sizes", instrument = inst, mean = 1),
                   file.path(dir, "precision_instruments.csv"), row.names = FALSE)
  r <- env$read_precision_results(dir)
  syn <- r$counts[r$counts$group == "Synthetic", ]
  want <- stats::aggregate(count ~ instrument + replicate + size_class,
                           counts[counts$group %in% c("PC", "PET"), ], sum)
  m <- merge(syn, want, by = c("instrument", "replicate", "size_class"))
  expect_equal(nrow(m), 12)
  expect_equal(m$count.x, m$count.y)
  expect_equal(r$iso$group[1:2], c("All particles", "Synthetic"))
  expect_false(grepl("No Synthetic", env$precision_summary_text(r)))
})

test_that("precision study labels name the runs and their input folder", {
  env <- .load_global_helpers()
  env$load_run_manifest <- function(run_dir)
    list(is_missing = FALSE, inputs = list(raman = list(path = "C:/Users/k/My Documents/PC_A_ref/raman/R.csv")))
  root <- withr::local_tempdir()
  a <- file.path(root, "2026-10-29_085018"); b <- file.path(root, "2026-10-30_101010")
  dir.create(a); dir.create(b)
  utils::write.csv(data.frame(instrument = "ftir", replicate = c("2026-10-28_1", "2026-10-29_1", "2026-10-29_2"),
                              file = "x.csv", source = "run", run_dir = file.path(root, "run")),
                   file.path(a, "precision_design.csv"), row.names = FALSE)
  utils::write.csv(data.frame(instrument = "ldir", replicate = 1:2,
                              file = "/home/k/Documents/lab/rep/ldir.csv"),
                   file.path(b, "precision_design.csv"), row.names = FALSE)
  lab <- env$precision_study_labels(c(a, b, file.path(root, "none")))
  expect_equal(lab[1], "2026-10-29_085018 \u2014 3 runs: 2026-10-28_1, 2026-10-29_1, 2026-10-29_2 \u2014 PC_A_ref/raman")
  expect_equal(lab[2], "2026-10-30_101010 \u2014 2 files, 2 replicates \u2014 lab/rep")
  expect_equal(lab[3], "none")
})

test_that("viewer result folders resolve against the project folder, not shiny_app/", {
  env <- .load_global_helpers()
  root <- normalizePath(REPO_ROOT, winslash = "/")
  expect_equal(normalizePath(env$project_path("output/precision"), winslash = "/", mustWork = FALSE),
               file.path(root, "output", "precision"))
  expect_equal(env$project_path("C:/data/prec"), "C:/data/prec")
  expect_equal(env$project_path("C:\\data\\prec"), "C:\\data\\prec")
  expect_equal(env$project_path("/srv/prec"), "/srv/prec")
})

test_that("LDIR outlines follow their particle into the aligned frame", {
  env <- .load_global_helpers()
  sq <- function(id, cx, cy, r) data.frame(particle_id = id, x_um = cx + c(-r, r, r, -r),
                                           y_um = cy + c(-r, -r, r, r))
  ol <- rbind(sq("A1", 0, 0, 10), sq("A2", 500, 200, 20), sq("A3", -300, 400, 5))
  # aligned = rotate 90 deg + shift, plus a small local offset on A2 (as a TPS would)
  df <- data.frame(particle_id = c("A1", "A2", "A3", "A4"), x_orig = c(0, 500, -300, 100),
                   y_orig = c(0, 200, 400, -600))
  df$x <- -df$y_orig + 1000; df$y <- df$x_orig + 50
  df$x[2] <- df$x[2] + 3
  nat <- env$ldir_outline_polygons(ol, df[1:2, ], aligned = FALSE)
  expect_equal(nrow(nat), 8)
  al <- env$ldir_outline_polygons(ol, df, aligned = TRUE)
  a2 <- al[al$particle_id == "A2", ]
  expect_equal(mean(a2$x), df$x[2], tolerance = 0.5)            # anchored on its own position
  expect_equal(diff(range(a2$x)), 40, tolerance = 1)               # size kept, turned 90 deg
  a1 <- al[al$particle_id == "A1", ]
  expect_equal(a1$x[2] - a1$x[1], 0, tolerance = 0.5)             # x offset turned into y
  expect_equal(a1$y[2] - a1$y[1], 20, tolerance = 0.5)
  expect_null(env$ldir_outline_polygons(NULL, df))
})


test_that("partner status: hover fill and colour classes", {
  g <- .load_global_helpers()
  ps <- data.frame(
    of = c("ldir", "ldir", "raman", "ldir", "ldir"),
    found_by = c("raman", "ftir_perkin", "ldir", "raman", "raman"),
    particle_id = c("A1", "A1", "R7", "A2", "A3"),
    passes = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    status = c("partner_below_cutoff", "not_found", "this_below_cutoff", "not_found", "matched"),
    partner_id = c("R7", NA, "A1", NA, "R9"),
    partner_material = c("PET", NA, "PET", NA, "PE"),
    partner_quality = c(55, NA, 0.9, NA, 88),
    partner_passes = c(FALSE, NA, TRUE, NA, TRUE),
    partner_distance_um = c(12, NA, 12, NA, 3), gate_um = 50, p_chance = c(0.02, NA, 0.03, NA, 0.01),
    stringsAsFactors = FALSE)
  # The user's case: LDIR PET hovered, Raman shown without a match
  f <- g$partner_fill(ps, list(ldir = "A1"), c("ldir", "raman"))
  expect_equal(f$raman$id, "R7")
  expect_false(f$raman$passes)
  n <- g$partner_fill_note(f$raman, "HQI")
  expect_true(n$warn)
  expect_match(n$text, "Below the HQI cutoff: not confidently identified")
  # FTIR not shown: nothing filled for it; nothing found: no fill
  expect_null(f$ftir_pe)
  expect_length(g$partner_fill(ps, list(ldir = "A2"), c("ldir", "raman")), 0)
  # The hovered particle's own cutoff state
  expect_false(attr(g$partner_fill(ps, list(raman = "R7"), c("ldir", "raman")), "own")$raman$passes)
  # Colour classes: best over the shown partners, matched gives none
  cl <- g$overlay_partner_class(ps, c("ldir", "raman", "ftir_pe"))
  expect_equal(cl$class[cl$inst == "ldir" & cl$particle_id == "A1"], "below_cutoff")
  expect_equal(cl$class[cl$inst == "raman" & cl$particle_id == "R7"], "below_cutoff")
  expect_equal(cl$class[cl$particle_id == "A2"], "not_found")
  expect_false("A3" %in% cl$particle_id)
  expect_null(g$overlay_partner_class(NULL, "ldir"))
})
