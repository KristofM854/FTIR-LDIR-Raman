# End-to-end on a synthetic DENSE field sample, shaped like Faanui SW1L R4 A:
# 12 000 Raman particles on a 13 mm filter (20% above HQI 70), an LDIR table of
# 1 400 particles with a known positional error (sigma 7 um), and an FTIR
# (Bruker) list of five particles with no counterpart (a list pre-filtered to
# plastics). Guards what the field sample taught: the LDIR gate, the
# positional-error fit, the chance estimate, and an alignment that a short list
# cannot fix being reported as not determined. About 40 s.

local({
  sys.source(file.path(REPO_ROOT, "R", "load.R"), envir = globalenv())
  load_pipeline(root = REPO_ROOT)
})

.field_inputs <- function(dir, sigma = 7) {
  set.seed(4242)
  n <- 12000; R <- 6400
  r <- R * sqrt(stats::runif(n)); a <- stats::runif(n, 0, 2 * pi)
  size <- pmin(round(exp(stats::rnorm(n, log(14), 0.7)), 1), 900)
  big <- sample(n, 250); size[big] <- round(stats::runif(250, 100, 600))
  t <- data.frame(k = seq_len(n), x = r * cos(a), y = r * sin(a), size = size)
  t$material <- sample(c("Calcium carbonate", "Cellulose", "Quartz", "P-DIAMINOBENZENE"),
                       n, replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))
  plast <- sample(which(t$size >= 50), 60)
  t$material[plast] <- sample(c("Polyethylene", "Polypropylene", "PET"), 60, replace = TRUE)
  hqi <- ifelse(stats::runif(n) < 0.2, 85, 40); hqi[plast] <- 92
  raman <- data.frame(ID = paste0("R", t$k), Xum = t$x, Yum = t$y, Feret = t$size,
                      Mat = t$material, HQI = hqi)
  utils::write.csv(raman, file.path(dir, "raman.csv"), row.names = FALSE)

  # LDIR: 1 400 particles >= 20 um, ~450 with quality >= 0.85; positions with
  # a 2-D normal error sigma, rotated -90 deg, shifted
  cand <- which(t$size >= 20)
  l <- t[sort(sample(cand, 1400)), ]
  q <- ifelse(stats::runif(nrow(l)) < 0.32, 0.92, 0.6)
  q[l$k %in% plast] <- 0.95
  th <- -pi / 2
  ex <- l$x + stats::rnorm(nrow(l), 0, sigma); ey <- l$y + stats::rnorm(nrow(l), 0, sigma)
  ldir <- data.frame(ID = paste0("L", l$k),
                     x = cos(th) * ex - sin(th) * ey + 500, y = sin(th) * ex + cos(th) * ey - 300,
                     size = l$size * 1.05, Identification = l$material, Quality = q)
  utils::write.csv(ldir, file.path(dir, "ldir.csv"), row.names = FALSE)

  # FTIR (Bruker): five particles with no counterpart in a 3 mm patch (a list
  # pre-filtered to plastics in the instrument software)
  bruker <- data.frame(ID = paste0("B", 1:5), x = stats::runif(5, 1000, 4000),
                       y = stats::runif(5, 4000, 7000), size = c(31, 40, 45, 61, 25),
                       material = "PET", score = 0.9)
  utils::write.csv(bruker, file.path(dir, "bruker.csv"), row.names = FALSE)
  invisible(t)
}

.field_config <- function(dir, out) {
  cfg <- make_config(raman_path = file.path(dir, "raman.csv"), output_dir = out)
  cfg$ldir_path <- file.path(dir, "ldir.csv")
  cfg$ftir_bruker_path <- file.path(dir, "bruker.csv")
  cfg$raman_column_map <- list(columns = list(particle_id = "ID", x_um = "Xum", y_um = "Yum",
    feret_max_um = "Feret", material = "Mat", quality = "HQI"), quality_metric = "HQI")
  cfg$ldir_column_map <- list(columns = list(particle_id = "ID", x_um = "x", y_um = "y",
    feret_max_um = "size", material = "Identification", quality = "Quality"))
  cfg$ftir_bruker_column_map <- list(columns = list(particle_id = "ID", x_um = "x", y_um = "y",
    feret_max_um = "size", material = "material", quality = "score"), quality_metric = "AAU")
  cfg
}

test_that("a dense field sample: LDIR pairs are right, sigma is recovered, a short list is not aligned", {
  dir <- withr::local_tempdir(); out <- withr::local_tempdir()
  t <- .field_inputs(dir)
  r <- withCallingHandlers(
    suppressMessages(run_pipeline(.field_config(dir, out), report = FALSE)),
    warning = function(w) if (grepl("conversion failure", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  run <- list.dirs(out, recursive = FALSE)

  # LDIR <-> Raman: nearly every pair is the true one, and nearly every true
  # pair that both quality cutoffs admit is found
  m <- r$ldir$ldir_ref_match$matched
  ok <- sub("^L", "", m$ldir_particle_id) == sub("^R", "", m$raman_particle_id)
  expect_gte(mean(ok), 0.95)
  ld <- utils::read.csv(file.path(dir, "ldir.csv"))
  ra <- utils::read.csv(file.path(dir, "raman.csv"))
  k  <- sub("^L", "", ld$ID[ld$Quality >= 0.85])
  possible <- sum(ra$HQI[match(paste0("R", k), ra$ID)] >= 70 &
                  ra$Feret[match(paste0("R", k), ra$ID)] >= 20)
  expect_gte(sum(ok) / possible, 0.85)

  ch <- utils::read.csv(file.path(run, "07_diagnostics", "chance_pairs.csv"))
  l <- ch[ch$pair == "ldir_raman", ]
  expect_equal(l$gate_um, 20)
  expect_lt(l$chance_share, 0.1)
  expect_gt(l$position_sigma_um, 4.5); expect_lt(l$position_sigma_um, 8.5)

  # FTIR (Bruker): five unrelated particles cannot fix a pose
  expect_false(r$ftir_bruker$alignment_determined)
  expect_false(ch$alignment_determined[ch$pair == "ftir_bruker_raman"])
  expect_false(file.exists(file.path(run, "06_agreement",
                                     "agreement_pairwise_ftir_bruker_raman.csv")))
  expect_false(file.exists(file.path(run, "05_matches", "matched_ftir_bruker_ldir.csv")))
})
