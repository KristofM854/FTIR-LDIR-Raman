# Implementation plan — Multi-device Shiny features

## Files touched
| File | Changes |
|---|---|
| `R/00b_file_input.R` | Add `FTIR_perkin` / `FTIR_bruker` patterns; update `group_files_by_instrument()` |
| `R/01_ingest.R` | Add `ingest_ftir_bruker()` thin alias |
| `main.R` | Guard FTIR against NULL; wire `ftir_bruker_file`; export bruker CSVs |
| `shiny_app/global.R` | `summarise_plastics()` helper; extend `build_instrument_dfs()` + `load_run_data()` for bruker |
| `shiny_app/app.R` | All UI/server changes (tabs, click, selection, summary) |

---

## 1 — Single-device robustness (Feature 4) — do first

### `shiny_app/app.R`

**`has_data` reactive — widen definition**

```r
has_data <- reactive({
  d <- run_data()
  !is.null(d$matched)      ||
  !is.null(d$unmatched_ftir)  ||
  !is.null(d$unmatched_raman) ||
  !is.null(d$ldir_raman_matched) ||
  !is.null(d$unmatched_ldir)
})
```

**`has_device` helper reactive**

```r
has_device <- reactive({
  dfs <- instrument_dfs()
  list(
    ftir  = !is.null(dfs$ftir)  && nrow(dfs$ftir)  > 0,
    raman = !is.null(dfs$raman) && nrow(dfs$raman) > 0,
    ldir  = !is.null(dfs$ldir)  && nrow(dfs$ldir)  > 0
  )
})
```

**Guard single-viewer plots** — wrap every `renderPlot` with an early exit that shows
a placeholder when the device has no data:

```r
output$ftir_plot <- renderPlot({
  df <- ftir_filtered()
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() +
      annotate("text", x=0, y=0, label="No FTIR data loaded",
               size=6, colour="grey50") +
      theme_void())
  }
  # ... existing plot code ...
})
```
Apply the same pattern to `raman_plot` and `ldir_plot`.

**Hide/disable overlay tab when fewer than 2 devices**

In server, observe `has_device()` and use `shinyjs::hideTab()` / `shinyjs::showTab()`:

```r
observe({
  hd <- has_device()
  n_active <- sum(unlist(hd))
  if (n_active < 2) {
    shinyjs::hideTab("main_tabs", "Overlay")
  } else {
    shinyjs::showTab("main_tabs", "Overlay")
  }
})
```

Add `library(shinyjs)` to `global.R` and `useShinyjs()` in the UI.

---

## 2 — FTIR_perkin rename + FTIR_bruker (Feature 3)

### `R/00b_file_input.R`

**`.INSTRUMENT_PATTERNS`** — replace `FTIR` entry with two entries:

```r
.INSTRUMENT_PATTERNS <- list(
  FTIR_perkin = "FTIR|Spotlight|infrared|FT-IR|PerkinElmer",
  FTIR_bruker = "Bruker|OPUS|bruker",
  Raman       = "Raman",
  LDIR        = "LDIR|8700"
)
```

**`group_files_by_instrument()`** — add `FTIR_perkin` and `FTIR_bruker` slots:

```r
group_files_by_instrument <- function(manifest) {
  result <- list(
    FTIR_perkin = list(tabular = NULL, image = NULL),
    FTIR_bruker = list(tabular = NULL, image = NULL),
    Raman       = list(tabular = NULL, image = NULL),
    LDIR        = list(tabular = NULL, image = NULL)
  )
  for (f in manifest) {
    inst <- f$instrument
    # backwards compat: bare "FTIR" maps to FTIR_perkin
    if (identical(inst, "FTIR")) inst <- "FTIR_perkin"
    if (!inst %in% names(result)) {
      log_message("  Unknown instrument type: ", inst, " — skipping", level = "WARN")
      next
    }
    result[[inst]][[f$type]] <- f$path
  }
  result
}
```

### `R/01_ingest.R`

Add at the end — a thin alias that calls `ingest_ftir()` (formats are identical initially):

```r
ingest_ftir_bruker <- function(filepath, sheet = "Long_Table") {
  log_message("Reading FTIR (Bruker) data from: ", filepath)
  df <- ingest_ftir(filepath, sheet = sheet)
  df$source_instrument <- "FTIR_bruker"
  df
}
```

### `main.R`

**Input mode block** — add bruker alongside perkin:

```r
ftir_perkin_file  <- grouped$FTIR_perkin$tabular
ftir_perkin_image <- grouped$FTIR_perkin$image
ftir_bruker_file  <- grouped$FTIR_bruker$tabular
ftir_bruker_image <- grouped$FTIR_bruker$image
```

**Config** — add:

```r
config$ftir_bruker_path  <- if (exists("ftir_bruker_file"))  ftir_bruker_file  else NULL
config$ftir_bruker_image <- if (exists("ftir_bruker_image")) ftir_bruker_image else NULL
```

**Ingestion block** — guard FTIR_bruker the same way LDIR is already guarded:

```r
ftir_bruker_raw <- NULL
if (!is.null(config$ftir_bruker_path)) {
  ftir_bruker_raw <- ingest_ftir_bruker(config$ftir_bruker_path)
  log_message("FTIR (Bruker): ", nrow(ftir_bruker_raw), " particles")
} else {
  log_message("FTIR (Bruker): not provided — skipping")
}
```

**Export** — write `unmatched_ftir_bruker.csv` when bruker data exists (no alignment yet — bruker is viewers-only in this iteration):

```r
if (!is.null(ftir_bruker_raw)) {
  write.csv(ftir_bruker_raw,
            file.path(config$output_dir, "unmatched_ftir_bruker.csv"),
            row.names = FALSE)
}
```

### `shiny_app/global.R`

**`load_run_data()`** — add bruker to the file map:

```r
file_map[["unmatched_ftir_bruker"]] <- "unmatched_ftir_bruker.csv"
```

**`build_instrument_dfs()`** — add `ftir_bruker` section after LDIR:

```r
# --- FTIR Bruker (viewer-only, unmatched CSV only) ---
if (!is.null(data$unmatched_ftir_bruker) && nrow(data$unmatched_ftir_bruker) > 0) {
  u <- data$unmatched_ftir_bruker
  result$ftir_bruker <- data.frame(
    particle_id = u$particle_id,
    x = u$x_um, y = u$y_um,
    x_orig = u$x_um, y_orig = u$y_um,
    area_um2 = u$area_um2,
    major_um = u$major_um, minor_um = u$minor_um,
    feret_max = u$feret_max_um,
    material = u$material, quality = u$quality,
    match_status = "unmatched", match_id = NA_integer_,
    stringsAsFactors = FALSE)
} else {
  result$ftir_bruker <- NULL
}
```

### `shiny_app/app.R`

**UI** — Rename "FTIR" tab → "FTIR (PerkinElmer)" and insert new tab:

```r
tabPanel("FTIR (PerkinElmer)",
  instrument_panel_ui("ftir", "AAU Quality", 0, 1, 0.01, 800)
),

tabPanel("FTIR (Bruker)",
  instrument_panel_ui("ftir_bruker", "AAU Quality", 0, 1, 0.01, 800)
),
```

**Server** — add `ftir_bruker_df_full`, debounced sliders, filtered reactive, `renderPlot`, hover, summary text — identical pattern to the FTIR tab, using the `ftir_bruker` prefix.

---

## 3 — Windows multi-file selection (Feature 1)

### `shiny_app/app.R`

Detect OS once at server start or in global.R:

```r
.is_windows <- tolower(Sys.info()[["sysname"]]) == "windows"
```

Replace every `fileInput(...)` in the "Upload Data" tab and in the instrument sidebar
background-image uploads with a wrapper:

```r
multi_file_input <- function(id, label, accept) {
  fileInput(id, label, accept = accept, multiple = .is_windows)
}
```

Replace all `fileInput(...)` calls with `multi_file_input(...)`.

**Downstream handling** — wherever `input$upload_matched$datapath` is used as a scalar,
wrap with `[1]` (first file) for non-Windows, or loop over all paths on Windows:

```r
# In the Upload tab observer:
for (path in input$upload_matched$datapath) {
  # bind / merge data
}
```

For background image uploads in single viewers (`ftir_image_upload`, etc.) — use only
`input$...$datapath[1]` (only one background image makes sense per viewer).

---

## 4 — Click-to-select in single viewers (Feature 2)

### `shiny_app/app.R` — UI

**`instrument_panel_ui()`** — add `click` to the `plotOutput` call and add a "Clear selection" button + selection detail box below the plot:

```r
plotOutput(paste0(id_prefix, "_plot"), height = "650px",
           click  = paste0(id_prefix, "_click"),         # NEW
           hover  = hoverOpts(paste0(id_prefix, "_hover"), delay = 100, delayType = "throttle"),
           brush  = brushOpts(paste0(id_prefix, "_brush"), resetOnNew = TRUE),
           dblclick = paste0(id_prefix, "_dblclick")),
tags$p(class = "text-muted",
       "Drag to zoom in. Double-click to reset zoom. Click a particle to add to selection."),
hr(),
div(class = "info-box",
    fluidRow(
      column(8, h5("Selected Particles")),
      column(4, actionButton(paste0(id_prefix, "_clear_selection"), "Clear",
                             class = "btn-sm btn-default",
                             style = "float:right; margin-top:2px;"))
    ),
    detail_table_ui(paste0(id_prefix, "_selection_info")))
```

### `shiny_app/app.R` — Server

**Global selection state** — one `reactiveVal` per instrument (accumulates across clicks):

```r
selected_ids <- reactiveValues(ftir = character(0), raman = character(0),
                                ldir = character(0), ftir_bruker = character(0))
```

**Click handler** — for each instrument (example for FTIR; replicate for raman/ldir/ftir_bruker):

```r
observeEvent(input$ftir_click, {
  click <- input$ftir_click
  if (is.null(click)) return()
  df <- ftir_filtered()
  if (is.null(df) || nrow(df) == 0) return()
  # snap distance = 5% of visible range
  vis <- if (!is.null(zoom$ftir)) zoom$ftir else
           list(x = range(df$x, na.rm=TRUE), y = range(df$y, na.rm=TRUE))
  snap <- max(diff(vis$x), diff(vis$y), 200) * 0.05
  d <- sqrt((df$x - click$x)^2 + (df$y - click$y)^2)
  idx <- which.min(d)
  if (length(idx) > 0 && d[idx] <= snap) {
    pid <- df$particle_id[idx]
    cur <- selected_ids$ftir
    selected_ids$ftir <- if (pid %in% cur) cur else c(cur, pid)
  }
})
```

**Clear button handlers:**

```r
observeEvent(input$ftir_clear_selection,  { selected_ids$ftir  <- character(0) })
observeEvent(input$raman_clear_selection, { selected_ids$raman <- character(0) })
observeEvent(input$ldir_clear_selection,  { selected_ids$ldir  <- character(0) })
```

**Selection detail renderer** — shows a table of all selected particles:

```r
output$ftir_selection_info <- renderUI({
  ids <- selected_ids$ftir
  if (length(ids) == 0)
    return(tags$p(class="text-muted", "Click particles to add to selection"))
  df <- ftir_df_full()
  if (is.null(df)) return(NULL)
  rows <- df[df$particle_id %in% ids, ]
  rows <- rows[match(ids, rows$particle_id), ]  # preserve click order
  tags$table(class = "hover-tbl",
    tags$tr(tags$th("ID"), tags$th("Material"), tags$th("Feret (\u00b5m)"),
            tags$th("Quality")),
    lapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      tags$tr(tags$td(r$particle_id), tags$td(r$material),
              tags$td(round(r$feret_max, 1)), tags$td(round(r$quality, 3)))
    })
  )
})
```

**Highlight selected particles on plot** — in `make_scatter()` or inline in each `renderPlot`,
pass `selected_ids$ftir` as an additional highlight layer (gold ring, same as overlay):

```r
# After existing highlight_id block, add:
sel <- selected_ids$ftir
if (length(sel) > 0) {
  sel_df <- if (!is.null(full_df)) full_df[full_df$particle_id %in% sel, ]
             else df[df$particle_id %in% sel, ]
  if (nrow(sel_df) > 0) {
    p <- p + geom_point(data = sel_df, aes(x = x, y = y),
                         shape = 21, size = 10, stroke = 2,
                         fill = NA, colour = "#FF6600")
  }
}
```

---

## 5 — Per-device plastics summary (Feature 5)

### `shiny_app/global.R`

Add after `compute_bounds()`:

```r
# Plastics (synthetic families) from material_map
.synthetic_families <- c("PET", "PP", "PE", "PS", "PVC", "PA", "PC",
                          "PMMA", "PU", "PTFE", "ABS", "Rubber")

#' Count particles per synthetic plastic family in a device data frame.
#' Returns data.frame(material, n) sorted by n desc, or empty frame.
summarise_plastics <- function(df, family_col = "material") {
  if (is.null(df) || nrow(df) == 0 || !family_col %in% names(df))
    return(data.frame(material = character(0), n = integer(0)))
  sub <- df[df[[family_col]] %in% .synthetic_families, ]
  if (nrow(sub) == 0) return(data.frame(material = character(0), n = integer(0)))
  tbl <- sort(table(sub[[family_col]]), decreasing = TRUE)
  data.frame(material = names(tbl), n = as.integer(tbl), stringsAsFactors = FALSE)
}
```

### `shiny_app/app.R`

**UI** — In `instrument_panel_ui()`, add a plastics summary box at the bottom of the main panel:

```r
hr(),
div(class = "info-box",
    h5("Plastics (synthetic)"),
    uiOutput(paste0(id_prefix, "_plastics_summary")))
```

**Server** — for each instrument tab (example FTIR):

```r
output$ftir_plastics_summary <- renderUI({
  df <- ftir_df_full()
  tbl <- summarise_plastics(df)
  if (nrow(tbl) == 0)
    return(tags$p(class = "text-muted", "No synthetic plastics found"))
  tags$table(class = "hover-tbl",
    tags$tr(tags$th("Material"), tags$th("Count")),
    lapply(seq_len(nrow(tbl)), function(i)
      tags$tr(tags$td(tbl$material[i]), tags$td(tbl$n[i])))
  )
})
```

Replicate for raman, ldir, ftir_bruker.

---

## 6 — Overall summary tab (Feature 6)

### `shiny_app/app.R` — UI

Insert before the "Run Info" tab:

```r
tabPanel("Summary",
  fluidRow(
    column(10, offset = 1,
      div(class = "info-box", style = "margin-top: 20px;",
        h4("Plastics by Instrument"),
        p(class = "text-muted",
          "Counts of synthetic plastic particles per device (all loaded data, no filters applied)."),
        uiOutput("summary_plastics_wide")
      )
    )
  )
)
```

### `shiny_app/app.R` — Server

```r
output$summary_plastics_wide <- renderUI({
  dfs <- instrument_dfs()
  devices <- list(
    "FTIR (PerkinElmer)" = dfs$ftir,
    "FTIR (Bruker)"      = dfs$ftir_bruker,
    Raman                = dfs$raman,
    LDIR                 = dfs$ldir
  )
  # Remove devices with no data
  devices <- Filter(function(d) !is.null(d) && nrow(d) > 0, devices)
  if (length(devices) == 0)
    return(tags$p(class = "text-muted", "No data loaded"))

  # Build per-device count vectors
  per_dev <- lapply(devices, summarise_plastics)
  all_mats <- sort(unique(unlist(lapply(per_dev, `[[`, "material"))))
  if (length(all_mats) == 0)
    return(tags$p(class = "text-muted", "No synthetic plastics found"))

  # Wide table
  dev_names <- names(devices)
  header <- tags$tr(tags$th("Material"),
                    lapply(dev_names, tags$th))
  body_rows <- lapply(all_mats, function(m) {
    cells <- lapply(per_dev, function(dt) {
      idx <- match(m, dt$material)
      if (is.na(idx)) tags$td("0") else tags$td(dt$n[idx])
    })
    tags$tr(tags$td(tags$b(m)), cells)
  })
  tags$table(class = "hover-tbl", header, body_rows)
})
```

---

## What is NOT changed

- CSV column names (`ftir_*`, `raman_*`, `ldir_*`) — no schema breaking changes
- Pipeline alignment/matching logic for FTIR — bruker is viewer-only for now
- `R/08b_material_map.R` — `synthetic_families` is re-declared as `.synthetic_families` in global.R (Shiny doesn't source the R/ pipeline files)
- Existing pipeline runs are fully backwards-compatible (new fields gracefully absent)

---

## Commit sequence

1. Feature 4 (robustness guards) — `shiny_app/app.R`, `shiny_app/global.R`
2. Feature 3 (FTIR_bruker plumbing) — `R/00b_file_input.R`, `R/01_ingest.R`, `main.R`, `shiny_app/global.R`, `shiny_app/app.R`
3. Feature 1 (Windows multi-file) — `shiny_app/app.R`, `shiny_app/global.R`
4. Feature 2 (click-to-select) — `shiny_app/app.R`
5. Features 5+6 (plastics summary) — `shiny_app/global.R`, `shiny_app/app.R`
