# =============================================================================
# 00b_file_input.R — recognise input files by their content
# =============================================================================
#
# main.R's default input mode opens one dialog in which all files of a sample
# are selected at once. classify_input_file() recognises each file from its
# content (table columns, Excel sheets, image headers, the first bytes), using
# the file name only where the content cannot decide; assign_input_roles()
# then gives each instrument its table, image and side files, and asks the
# operator whenever there are two candidates for one role or the content
# cannot tell (which FTIR table is the Bruker one, for instance).
#
# What decides, per file (checked on dataset A and the Faanui field sample):
#   FTIR table      columns Identifier, Coord. [pixels], Coord. [µm], AAU score.
#                   PerkinElmer and Bruker exports have identical columns, so
#                   the instrument comes from the paired image (same file name
#                   stem), the file name, or the other FTIR table.
#   segment tables  FTIR columns and a name ending in M<k> (segmented
#                   PerkinElmer scan; R/ftir_segments.R)
#   Raman table     Particle Name, HQI, Visual Center Point X (WITec)
#   LDIR table      Excel workbook with a Particles sheet (Identification,
#                   Quality, Width / Diameter); its Info sheet names the 8700
#   LDIR .sam       "Agilent.Clarity" in its first kilobytes (byte 112 on Faanui)
#   coordinate swaps  two-column table with id_clarity / id_R style headers
#   WITec placement   .txt with Width, Height, X, Y lines
#   FTIR placement    .txt starting "FTIR image placement"
#   LDIR analyzed image  few colours on black (Clarity's particle colouring)
#   BMP images      the header: Bruker Lumos 32-bit; WITec mosaics without a
#                   resolution field (px/m 0); others by name or pairing
# =============================================================================

INPUT_ROLE_LABELS <- c(ftir = "FTIR (PerkinElmer)", ftir_bruker = "FTIR (Bruker)",
                       raman = "Raman", ldir = "LDIR")

# File-name hints, used only where the content cannot decide
.INPUT_NAME_HINTS <- list(
  ftir        = "perkin|spotlight",
  ftir_bruker = "bruker|lumos|opus",
  raman       = "raman|witec",
  ldir        = "ldir|8700|clarity")

.INPUT_IMAGE_EXT <- c("png", "jpg", "jpeg", "tif", "tiff", "bmp", "webp")

# Header names the coordinate-swap loader accepts (load_ldir_coord_swaps())
.SWAP_HEADERS <- c("id_a", "from", "id1", "particle_a", "id_clarity", "clarity", "software",
                   "correct", "true", "id_b", "to", "id2", "particle_b", "id_r",
                   "id_pipeline", "pipeline", "assigned", "current")

.norm_header <- function(x) gsub("[^a-z0-9]+", "", tolower(iconv(x, "", "ASCII", sub = "")))

.name_hint <- function(path) {
  b <- tolower(basename(path))
  hit <- names(.INPUT_NAME_HINTS)[vapply(.INPUT_NAME_HINTS, grepl, logical(1), b)]
  if (length(hit) == 1) hit else NA_character_
}

.file_stem <- function(path) tolower(tools::file_path_sans_ext(basename(path)))

# Header line of a delimited text file, split on its most frequent separator
.text_header <- function(path) {
  l <- tryCatch(readLines(path, n = 1L, warn = FALSE, encoding = "latin1"),
                error = function(e) character(0))
  if (!length(l)) return(character(0))
  l <- sub("^\xef\xbb\xbf", "", l, useBytes = TRUE)
  seps <- c(";", ",", "\t")
  n <- vapply(seps, function(s) lengths(regmatches(l, gregexpr(s, l, fixed = TRUE))), 1L)
  trimws(strsplit(l, seps[which.max(n)], fixed = TRUE)[[1]])
}

.table_kind <- function(cols) {
  h <- .norm_header(cols)
  has <- function(p) any(startsWith(h, p))
  if (has("particlename") && has("hqi") && has("visualcenterpointx")) return("raman_table")
  if (has("identifier") && has("coord") && has("maxaauscore")) return("ftir_table")
  if (any(h %in% c("identification")) && any(h %in% c("quality")) &&
      (has("width") || has("diameter"))) return("ldir_table")
  if (length(h) <= 4 && sum(h %in% .norm_header(.SWAP_HEADERS)) >= 2) return("ldir_coord_swaps")
  NA_character_
}

# BMP header fields (first 54 bytes)
.bmp_header <- function(path) {
  b <- tryCatch(readBin(path, "raw", 54L), error = function(e) raw(0))
  if (length(b) < 54 || !identical(rawToChar(b[1:2]), "BM")) return(NULL)
  i32 <- function(i) readBin(b[i:(i + 3)], "integer", size = 4, endian = "little")
  list(width = abs(i32(19)), height = abs(i32(23)),
       bits = readBin(b[29:30], "integer", size = 2, signed = FALSE, endian = "little"),
       ppm_x = i32(39), ppm_y = i32(43))
}

# PNG width and height from its IHDR chunk
.png_dims <- function(path) {
  b <- tryCatch(readBin(path, "raw", 24L), error = function(e) raw(0))
  if (length(b) < 24 || !identical(b[2:4], charToRaw("PNG"))) return(NULL)
  u32 <- function(i) sum(as.numeric(b[i:(i + 3)]) * 256^(3:0))
  list(width = u32(17), height = u32(21))
}

# Clarity's analyzed image: a few flat colours on black. Counted on an evenly
# spaced sample of pixels (no resampling, so the colours stay exact).
.is_palette_image <- function(path, max_px = 40e6, max_colours = 32L) {
  d <- .png_dims(path)
  if (is.null(d) || d$width * d$height > max_px || !requireNamespace("png", quietly = TRUE))
    return(FALSE)
  a <- tryCatch(png::readPNG(path), error = function(e) NULL)
  if (is.null(a) || length(dim(a)) < 3) return(FALSE)
  n <- dim(a)[1] * dim(a)[2]
  idx <- unique(round(seq(1, n, length.out = min(n, 2e5))))
  rgb <- vapply(1:3, function(k) as.vector(a[, , k])[idx], numeric(length(idx)))
  key <- round(rgb[, 1] * 255) * 65536 + round(rgb[, 2] * 255) * 256 + round(rgb[, 3] * 255)
  length(unique(key)) <= max_colours && mean(key == 0) >= 0.5
}

#' Recognise one input file from its content.
#'
#' @param path File
#' @return list(path, kind, role, reason, info): kind one of ftir_table,
#'   ftir_segment_table, raman_table, ldir_table, ldir_coord_swaps, ldir_sam,
#'   raman_placement, ftir_placement, image, ldir_analyzed_image, unknown; role
#'   the instrument when the content (or, failing that, the name) says so, else
#'   NA; reason a short explanation for the confirmation list; info the image
#'   size where read
classify_input_file <- function(path) {
  out <- function(kind, role = NA_character_, reason = "", info = list())
    list(path = path, kind = kind, role = role, reason = reason, info = info)
  if (!file.exists(path) || dir.exists(path)) return(out("unknown", reason = "not a file"))
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("csv", "txt", "tsv")) {
    if (ext == "txt") {
      first <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) "")
      if (length(first) && grepl("FTIR image placement", first[1], ignore.case = TRUE))
        return(out("ftir_placement", reason = "starts 'FTIR image placement'"))
      w <- if (exists("read_witec_image_placement", mode = "function"))
        suppressMessages(read_witec_image_placement(path, quiet = TRUE)) else NULL
      if (!is.null(w))
        return(out("raman_placement", "raman", "Width / Height / X / Y lines (WITec)",
                   info = list(aspect = w$width_um / w$height_um)))
    }
    kind <- .table_kind(.text_header(path))
    if (!is.na(kind)) return(.classified_table(path, kind, out))
    return(out("unknown", reason = "columns not recognised"))
  }

  if (ext %in% c("xlsx", "xls")) {
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
    for (s in c(intersect(c("Particles", "Data", "Long_Table"), sheets), sheets)) {
      cols <- tryCatch(names(readxl::read_excel(path, sheet = s, n_max = 0)),
                       error = function(e) character(0))
      kind <- .table_kind(cols)
      if (!is.na(kind)) {
        r <- .classified_table(path, kind, out)
        if (kind == "ldir_table" && "Info" %in% sheets) r$reason <- paste0(r$reason, "; Info sheet")
        return(r)
      }
    }
    return(out("unknown", reason = "no sheet with known columns"))
  }

  if (ext == "sam") {
    b <- tryCatch(readBin(path, "raw", 4096L), error = function(e) raw(0))
    if (length(grepRaw("Agilent.Clarity", b, fixed = TRUE)))
      return(out("ldir_sam", "ldir", "Agilent Clarity project file"))
    return(out("unknown", reason = ".sam without the Clarity signature"))
  }

  if (ext %in% .INPUT_IMAGE_EXT) {
    if (grepl(FTIR_SEGMENT_IMAGE, basename(path), ignore.case = TRUE))
      return(out("image", "ftir", "segment image (M<k>) of a PerkinElmer scan",
                 if (ext == "bmp") .bmp_header(path) %||% list() else list()))
    if (ext == "bmp") {
      h <- .bmp_header(path)
      info <- if (is.null(h)) list() else h
      if (!is.null(h) && h$bits == 32)
        return(out("image", "ftir_bruker", "32-bit BMP, as the Lumos camera writes", info))
      if (!is.null(h) && h$ppm_x == 0 && min(h$width, h$height) >= 10000)
        return(out("image", "raman", "large BMP mosaic without a resolution field (WITec)", info))
    } else if (ext == "png") {
      info <- .png_dims(path) %||% list()
      if (.is_palette_image(path))
        return(out("ldir_analyzed_image", "ldir", "a few flat colours on black (Clarity analyzed image)", info))
    } else info <- list()
    hint <- .name_hint(path)
    return(out("image", hint, if (is.na(hint)) "" else "file name", info))
  }
  out("unknown", reason = "file type not used")
}

.classified_table <- function(path, kind, out) {
  if (kind == "ftir_table" && grepl(FTIR_SEGMENT_TABLE, basename(path), ignore.case = TRUE))
    return(out("ftir_segment_table", "ftir", "FTIR columns; segment M<k> of a PerkinElmer scan"))
  switch(kind,
    raman_table = out(kind, "raman", "WITec columns (Particle Name, HQI, Visual Center Point)"),
    ldir_table  = out(kind, "ldir", "LDIR Particles sheet (Identification, Quality)"),
    ldir_coord_swaps = out(kind, "ldir", "coordinate-swap table (id_clarity / id_R)"),
    ftir_table  = {
      hint <- .name_hint(path)
      hint <- if (hint %in% c("ftir", "ftir_bruker")) hint else NA_character_
      out(kind, hint, paste0("FTIR columns (Identifier, Coord., AAU score)",
                             if (!is.na(hint)) "; instrument from the file name" else ""))
    })
}


#' Give each instrument its files.
#'
#' @param files Paths the operator selected
#' @param ask function(question, choices) -> one of choices, called whenever
#'   the content does not decide; NULL answers "not used"
#' @return list(paths = list(ftir_file, ftir_image, ftir_bruker_file,
#'   ftir_bruker_image, raman_file, raman_image, ldir_file, ldir_image),
#'   config = settings for the side files (ldir_sam_path,
#'   raman_image_metadata, ldir_coord_swaps_path, ldir_processed_image_path),
#'   lines = one line per file for the confirmation, notes = warnings)
assign_input_roles <- function(files, ask) {
  cls <- lapply(unique(files), classify_input_file)
  kind <- vapply(cls, `[[`, "", "kind")
  role <- vapply(cls, function(x) x$role %||% NA_character_, "")
  why  <- vapply(cls, `[[`, "", "reason")
  used <- rep(NA_character_, length(cls))     # what each file became
  notes <- character(0)
  paths <- list(); cfg <- list()
  base <- function(i) basename(cls[[i]]$path)
  pick <- function(idx, question) {
    if (length(idx) <= 1) return(idx)
    a <- ask(question, vapply(idx, base, ""))
    idx[match(a, vapply(idx, base, ""))][!is.na(match(a, vapply(idx, base, "")))]
  }

  # --- Tables -----------------------------------------------------------------
  take_table <- function(k, r, slot, label) {
    idx <- pick(which(kind == k & role %in% r),
                paste0("Several ", label, " tables: which one is this run's?"))
    if (length(idx)) { paths[[slot]] <<- cls[[idx]]$path; used[idx] <<- slot }
  }
  take_table("raman_table", "raman", "raman_file", "Raman")
  take_table("ldir_table", "ldir", "ldir_file", "LDIR")

  # FTIR: a segmented scan is PerkinElmer; the content cannot tell the two
  # instruments' tables apart, so the paired image (same name stem), the
  # file name or the other table decides, else the operator
  ft <- which(kind == "ftir_table")
  stems <- vapply(cls, function(x) .file_stem(x$path), "")
  # An FTIR image recognised by its content names the table of the same stem
  for (i in ft) if (is.na(role[i])) {
    j <- which(kind == "image" & role %in% c("ftir", "ftir_bruker") & stems == stems[i])
    if (length(j) == 1) {
      role[i] <- role[j]
      why[i] <- paste0(why[i], "; same name as the ", INPUT_ROLE_LABELS[[role[j]]], " image")
    }
  }
  # Segmented PerkinElmer scans: the segments of one scan share their folder
  # and name before M<k>. One scan is used; a single PerkinElmer table among
  # the selected files competes with it, and the operator chooses.
  seg <- which(kind == "ftir_segment_table")
  if (length(seg)) {
    scan_key <- function(i, pattern) {
      p <- cls[[i]]$path
      file.path(dirname(p), sub(pattern, "", basename(p), ignore.case = TRUE))
    }
    seg_img <- which(kind == "image" & grepl(FTIR_SEGMENT_IMAGE, vapply(seq_along(cls), base, ""),
                                             ignore.case = TRUE))
    key_t <- vapply(seg, scan_key, "", FTIR_SEGMENT_TABLE)
    # Images belong to the scan whose tables share their folder: the exports
    # name them differently ("Average Abs.( <scan> )M1.bmp" next to
    # "<scan>M1.csv"), and the stitching pairs them by M<k>. Only a folder
    # with two scans needs the name before M<k> to match.
    key_i <- vapply(seg_img, function(i) {
      in_dir <- unique(key_t[dirname(key_t) == dirname(cls[[i]]$path)])
      if (length(in_dir) == 1) in_dir else scan_key(i, FTIR_SEGMENT_IMAGE)
    }, "")
    scans <- unique(key_t)
    lab <- function(k) sprintf("segmented scan %s* (%d segments, folder %s)", basename(k),
                               sum(key_t == k), basename(dirname(k)))
    keep <- scans
    if (length(scans) > 1) {
      a <- ask("Several segmented PerkinElmer scans: which one is this run's?",
               c(vapply(scans, lab, ""), "none"))
      keep <- scans[match(a, vapply(scans, lab, ""))]
    }
    single <- ft[role[ft] %in% "ftir"]
    if (!is.na(keep) && length(single)) {
      a <- ask("FTIR (PerkinElmer): which input is this run's?",
               c(lab(keep), vapply(single, base, "")))
      if (!identical(a, lab(keep))) keep <- NA_character_
    }
    drop <- c(seg[is.na(keep) | key_t != keep], seg_img[is.na(keep) | key_i != keep])
    kind[drop] <- "unused_segment"
    why[drop] <- paste0(why[drop], "; not used: ",
                        if (is.na(keep)) "another FTIR (PerkinElmer) input was chosen" else "another scan was chosen")
    seg <- setdiff(seg, drop)
    if (length(seg)) {
      first <- seg[order(vapply(seg, function(i) .segment_number(cls[[i]]$path, FTIR_SEGMENT_TABLE), 1L))][1]
      paths$ftir_file <- cls[[first]]$path
      used[seg] <- "ftir_file"; why[seg] <- paste(why[seg], "(stitched before the run)")
      role[single] <- "none"
      why[single] <- paste0(why[single], "; not used: the segmented scan is used")
    }
  }
  free <- setdiff(c("ftir", "ftir_bruker"), if (length(seg)) "ftir" else character(0))
  if (length(ft)) {
    known <- ft[!is.na(role[ft]) & role[ft] %in% free]
    # One left open and one role left: it takes that role
    open <- ft[is.na(role[ft])]
    left <- setdiff(free, role[known])
    if (length(open) == 1 && length(left) == 1 && (length(known) >= 1 || length(seg))) {
      other <- if (length(known)) role[known[1]] else "ftir"
      role[open] <- left
      why[open] <- paste0(why[open], "; the other FTIR ",
                          if (length(known)) "table" else "input (the segmented scan)", " is ",
                          INPUT_ROLE_LABELS[[other]])
    } else if (length(open) && length(left)) {
      for (i in open) {
        left <- setdiff(free, role[ft[!is.na(role[ft])]])
        if (!length(left)) break
        rest <- open[is.na(role[open])]
        if (length(rest) == 1 && length(left) == 1 && any(!is.na(role[ft]))) {
          role[i] <- left
          why[i] <- paste0(why[i], "; the other FTIR table is ", INPUT_ROLE_LABELS[[setdiff(free, left)[1]]])
          next
        }
        a <- ask(paste0("Which instrument made the FTIR table ", base(i), "?"),
                 c(unname(INPUT_ROLE_LABELS[left]), "not used"))
        r <- names(INPUT_ROLE_LABELS)[match(a, INPUT_ROLE_LABELS)]
        if (!is.na(r)) { role[i] <- r; why[i] <- paste0(why[i], "; instrument chosen by you") }
      }
    }
    for (r in free) {
      idx <- pick(ft[role[ft] %in% r], paste0("Several ", INPUT_ROLE_LABELS[[r]], " tables: which one is this run's?"))
      if (length(idx)) {
        slot <- if (r == "ftir") "ftir_file" else "ftir_bruker_file"
        paths[[slot]] <- cls[[idx]]$path; used[idx] <- slot
      }
    }
  }
  have <- c(ftir = !is.null(paths$ftir_file), ftir_bruker = !is.null(paths$ftir_bruker_file),
            raman = !is.null(paths$raman_file), ldir = !is.null(paths$ldir_file))

  # --- Images -------------------------------------------------------------------
  slot_of <- c(ftir = "ftir_image", ftir_bruker = "ftir_bruker_image",
               raman = "raman_image", ldir = "ldir_image")
  im <- which(kind == "image")
  an <- which(kind == "ldir_analyzed_image")
  if (length(an) > 1) an <- pick(an, "Several LDIR analyzed images: which one is this run's?")
  if (length(an)) {
    used[an] <- "ldir_analyzed"
    # Its raw image: the same name without the suffix, else a PNG of the same size
    a_stem <- sub("_analyzed$", "", .file_stem(cls[[an]]$path))
    raw_i <- im[vapply(im, function(i) .file_stem(cls[[i]]$path) == a_stem, logical(1))]
    if (!length(raw_i)) raw_i <- im[vapply(im, function(i) identical(
      cls[[i]]$info[c("width", "height")], cls[[an]]$info[c("width", "height")]) &&
        tolower(tools::file_ext(cls[[i]]$path)) == "png", logical(1))]
    if (length(raw_i) == 1) {
      role[raw_i] <- "ldir"; why[raw_i] <- paste0("LDIR image of ", base(an))
      if (.file_stem(cls[[raw_i]]$path) != a_stem)
        cfg$ldir_processed_image_path <- cls[[an]]$path
    }
  }
  # The WITec placement file's aspect ratio points at the Raman image
  wt <- which(kind == "raman_placement")
  if (length(wt) == 1 && !"raman" %in% role[im]) {
    asp <- cls[[wt]]$info$aspect
    fit <- im[is.na(role[im]) & vapply(im, function(i) {
      d <- cls[[i]]$info
      length(d$width) && length(d$height) && abs(d$width / d$height / asp - 1) < 0.02
    }, logical(1))]
    if (length(fit) == 1) { role[fit] <- "raman"; why[fit] <- paste0("aspect ratio of ", base(wt)) }
  }
  for (r in names(slot_of)) {
    idx <- im[role[im] %in% r]
    if (!have[[r]]) {
      why[idx] <- paste0(why[idx], "; no ", INPUT_ROLE_LABELS[[r]], " table selected")
      next
    }
    if (r == "ftir" && length(seg)) {
      # The stitched table goes with the stitched image only: the scan's own
      # segment images, never another PerkinElmer image (resolve_ftir_segments()
      # stitches both before the run)
      segi <- idx[grepl(FTIR_SEGMENT_IMAGE, vapply(idx, base, ""), ignore.case = TRUE)]
      other <- setdiff(idx, segi)
      why[other] <- paste0(why[other], "; not used: the segmented scan is placed on its own segment images")
      role[other] <- "none"
      if (length(segi)) {
        segi <- segi[order(vapply(segi, function(i) .segment_number(cls[[i]]$path, FTIR_SEGMENT_IMAGE), 1L))]
        paths$ftir_image <- cls[[segi[1]]]$path; used[segi] <- "ftir_image"
        why[segi] <- paste(why[segi], "(stitched with the tables)")
      } else {
        notes <- c(notes, paste0("FTIR (PerkinElmer): the segment images (", "...M1.bmp ...) ",
                                 "were not found with the segment tables; the run stops without them"))
      }
      next
    }
    idx <- pick(idx, paste0("Several ", INPUT_ROLE_LABELS[[r]], " images: which one is this run's?"))
    if (length(idx)) { paths[[slot_of[[r]]]] <- cls[[idx]]$path; used[idx] <- slot_of[[r]] }
  }
  # Images nobody claimed: the operator says whose they are
  for (i in im[is.na(used[im]) & is.na(role[im])]) {
    open_roles <- names(slot_of)[have & vapply(slot_of, function(s) is.null(paths[[s]]), logical(1))]
    if (!length(open_roles)) break
    a <- ask(paste0("Which instrument's image is ", base(i), "?"),
             c(unname(INPUT_ROLE_LABELS[open_roles]), "not used"))
    r <- names(INPUT_ROLE_LABELS)[match(a, INPUT_ROLE_LABELS)]
    if (!is.na(r) && r %in% open_roles) {
      paths[[slot_of[[r]]]] <- cls[[i]]$path; used[i] <- slot_of[[r]]; role[i] <- r
      why[i] <- "chosen by you"
    }
  }

  # --- Side files ---------------------------------------------------------------
  side <- function(k, key, label, needs) {
    idx <- which(kind == k)
    if (!length(idx)) return(invisible())
    if (!isTRUE(have[[needs]])) { why[idx] <<- paste0(why[idx], "; no ", INPUT_ROLE_LABELS[[needs]], " table"); return(invisible()) }
    idx <- pick(idx, paste0("Several ", label, ": which one is this run's?"))
    if (length(idx)) { cfg[[key]] <<- cls[[idx]]$path; used[idx] <<- key }
  }
  side("ldir_sam", "ldir_sam_path", "LDIR project files (.sam)", "ldir")
  side("raman_placement", "raman_image_metadata", "Raman placement files", "raman")
  side("ldir_coord_swaps", "ldir_coord_swaps_path", "coordinate-swap tables", "ldir")
  # FTIR placement files are read next to their image, by its name
  for (i in which(kind == "ftir_placement")) {
    img <- unlist(paths[c("ftir_image", "ftir_bruker_image")])
    ok <- .file_stem(cls[[i]]$path) %in% vapply(img, .file_stem, "")
    if (ok) used[i] <- "ftir placement" else
      why[i] <- paste0(why[i], "; not used: it must be named like its image (<image name>.txt)")
  }
  if (!is.null(paths$ldir_file) && is.null(paths$ldir_image))
    notes <- c(notes, "LDIR: no image selected; the LDIR particles need the image (and its _analyzed image) for their positions")
  if (length(an) && is.null(paths$ldir_image))
    notes <- c(notes, "LDIR: the analyzed image is used only together with the LDIR image it was made from")

  # --- Confirmation lines -----------------------------------------------------
  label_of <- c(ftir_file = "FTIR (PerkinElmer) table", ftir_image = "FTIR (PerkinElmer) image",
                ftir_bruker_file = "FTIR (Bruker) table", ftir_bruker_image = "FTIR (Bruker) image",
                raman_file = "Raman table", raman_image = "Raman image",
                ldir_file = "LDIR table", ldir_image = "LDIR image",
                ldir_analyzed = "LDIR analyzed image", ldir_sam_path = "LDIR project file",
                raman_image_metadata = "Raman placement", ldir_coord_swaps_path = "LDIR coordinate swaps",
                `ftir placement` = "FTIR image placement")
  ord <- order(match(used, names(label_of)), na.last = TRUE)
  lines <- vapply(ord, function(i) sprintf("%-26s %s%s",
    if (is.na(used[i])) "not used" else label_of[[used[i]]], base(i),
    if (nzchar(why[i])) paste0("  (", why[i], ")") else ""), "")
  list(paths = paths, config = cfg, lines = lines, notes = notes,
       files = data.frame(path = vapply(cls, `[[`, "", "path"), kind = kind, role = role,
                          used_as = used, reason = why, stringsAsFactors = FALSE))
}


#' Segmented PerkinElmer scans in the subfolders of the selected files' folders.
#'
#' A file dialog selects files of one folder only; the segments of a scan
#' (<name>M1.csv ... with <name>M1.bmp ...) often sit in a subfolder of it. One
#' level down from each selected file's folder, every folder with at least two
#' segment tables contributes its segment tables and images; assign_input_roles()
#' then decides whether they are used.
#' @param files The selected paths
#' @return Paths of the segment files found (not among `files`)
find_segment_scans <- function(files) {
  files <- files[file.exists(files)]
  if (!length(files)) return(character(0))
  subs <- unique(unlist(lapply(unique(dirname(files)), list.dirs, recursive = FALSE)))
  out <- character(0)
  for (d in subs) {
    f <- list.files(d, full.names = TRUE)
    tab <- f[grepl(FTIR_SEGMENT_TABLE, basename(f), ignore.case = TRUE)]
    if (length(tab) >= 2)
      out <- c(out, tab, f[grepl(FTIR_SEGMENT_IMAGE, basename(f), ignore.case = TRUE)])
  }
  norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
  out[!norm(out) %in% norm(files)]
}
