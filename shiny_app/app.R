# =============================================================================
# app.R — Multi-Instrument Particle Viewer (Shiny + ggplot2)
# =============================================================================

source("global.R")

# ============================================================================
# Shared UI helpers
# ============================================================================

detail_table_ui <- function(id) {
  uiOutput(id)
}

make_detail_row <- function(label, value) {
  tags$tr(tags$td(tags$b(label)), tags$td(value))
}

instrument_panel_ui <- function(id_prefix, quality_label, quality_min, quality_max,
                                 quality_step, size_max = 1200) {
  sidebarLayout(
    sidebarPanel(width = 3,
      h4(paste0(toupper(id_prefix), " Filters")),
      sliderInput(paste0(id_prefix, "_quality_range"), quality_label,
                  min = quality_min, max = quality_max,
                  value = c(quality_min, quality_max), step = quality_step),
      sliderInput(paste0(id_prefix, "_size_range"), "Feret Max (\u00b5m)",
                  min = 0, max = size_max, value = c(0, size_max), step = 5),
      selectInput(paste0(id_prefix, "_material_filter"), "Material",
                  choices = c("All"), selected = "All", multiple = TRUE),
      checkboxGroupInput(paste0(id_prefix, "_match_filter"), "Match Status",
                         choices = c("matched", "unmatched"),
                         selected = c("matched", "unmatched"), inline = TRUE),
      hr(),
      div(class = "info-box",
          h5("Summary"), textOutput(paste0(id_prefix, "_summary_text"))),
      hr(),
      fileInput(paste0(id_prefix, "_image_upload"), "Background Image",
                accept = c("image/png", "image/jpeg")),
      fluidRow(
        column(6, numericInput(paste0(id_prefix, "_img_offset_x"),
                               "Img X offset (\u00b5m)", value = 0, step = 25)),
        column(6, numericInput(paste0(id_prefix, "_img_offset_y"),
                               "Img Y offset (\u00b5m)", value = 0, step = 25))
      )
    ),
    mainPanel(width = 9,
      plotOutput(paste0(id_prefix, "_plot"), height = "650px",
                 hover = hoverOpts(paste0(id_prefix, "_hover"), delay = 100,
                                   delayType = "throttle"),
                 brush = brushOpts(paste0(id_prefix, "_brush"),
                                   resetOnNew = TRUE),
                 dblclick = paste0(id_prefix, "_dblclick")),
      tags$p(class = "text-muted",
             "Drag to zoom in. Double-click to reset zoom."),
      hr(),
      div(class = "info-box",
          h5("Particle Details (hover)"),
          detail_table_ui(paste0(id_prefix, "_hover_info")))
    )
  )
}


# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background-color: #f8f9fa; font-size: 14px; }
    .navbar { margin-bottom: 8px; }
    .info-box { background: white; border-radius: 6px; padding: 12px;
                 box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin-bottom: 8px; }
    .hover-tbl { width: 100%; font-size: 13px; border-collapse: collapse; }
    .hover-tbl th { background: #e9ecef; padding: 4px 8px; text-align: left; }
    .hover-tbl td { padding: 4px 8px; border-bottom: 1px solid #dee2e6; }
    .placeholder-msg { text-align: center; padding: 80px 20px; color: #6c757d; }
    .placeholder-msg h3 { color: #495057; }
  "))),

  navbarPage(
    title = "Multi-Instrument Particle Viewer",
    id = "main_tabs",

    # Tab 1: FTIR
    tabPanel("FTIR",
      instrument_panel_ui("ftir", "AAU Quality", 0, 1, 0.01, 800)
    ),

    # Tab 2: Raman
    tabPanel("Raman",
      instrument_panel_ui("raman", "HQI", 0, 100, 1, 1200)
    ),

    # Tab 3: LDIR
    tabPanel("LDIR",
      instrument_panel_ui("ldir", "Quality", 0, 1, 0.01, 1200)
    ),

    # Tab 4: Overlay (FTIR + Raman)
    tabPanel("Overlay",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Overlay Filters"),
          sliderInput("overlay_hqi_range", "Raman HQI",
                      min = 0, max = 100, value = c(0, 100), step = 1),
          sliderInput("overlay_quality_range", "FTIR AAU Quality",
                      min = 0, max = 1, value = c(0, 1), step = 0.01),
          sliderInput("overlay_size_range", "Feret Max (\u00b5m)",
                      min = 0, max = 1200, value = c(0, 1200), step = 5),
          sliderInput("overlay_dist_range", "Match Distance (\u00b5m)",
                      min = 0, max = 100, value = c(0, 100), step = 1),
          selectInput("overlay_material_filter", "Material (either instrument)",
                      choices = c("All"), selected = "All", multiple = TRUE),
          checkboxGroupInput("overlay_layers", "Show Layers",
                             choices = c("Matched pairs" = "matched",
                                         "Unmatched FTIR" = "unmatched_ftir",
                                         "Unmatched Raman" = "unmatched_raman",
                                         "FTIR-Raman lines" = "match_lines",
                                         "LDIR matched" = "ldir_matched",
                                         "LDIR unmatched" = "ldir_unmatched",
                                         "LDIR-Raman lines" = "ldir_lines"),
                             selected = c("matched", "unmatched_ftir",
                                          "unmatched_raman", "ldir_matched"),
                             inline = FALSE),
          hr(),
          div(class = "info-box",
              h5("Match Summary"), textOutput("overlay_summary_text")),
          hr(),
          fileInput("overlay_image_upload", "Background Image",
                    accept = c("image/png", "image/jpeg")),
          fluidRow(
            column(6, numericInput("overlay_img_offset_x",
                                   "Img X offset (\u00b5m)", value = 0, step = 25)),
            column(6, numericInput("overlay_img_offset_y",
                                   "Img Y offset (\u00b5m)", value = 0, step = 25))
          )
        ),
        mainPanel(width = 9,
          plotOutput("overlay_plot", height = "650px",
                     hover = hoverOpts("overlay_hover", delay = 100,
                                       delayType = "throttle"),
                     brush = brushOpts("overlay_brush",
                                       resetOnNew = TRUE),
                     dblclick = "overlay_dblclick"),
          tags$p(class = "text-muted",
                 "Drag to zoom in. Double-click to reset zoom."),
          hr(),
          div(class = "info-box",
              h5("Match Details (hover on matched particle)"),
              detail_table_ui("overlay_hover_info"))
        )
      )
    ),

    # Tab 5: Data Upload (fallback)
    tabPanel("Upload Data",
      fluidRow(
        column(6, offset = 3,
          div(class = "info-box", style = "margin-top: 20px;",
            h4("Upload Pipeline Output"),
            p("If the app cannot find pipeline output automatically, you can
               upload the CSV files here. At minimum, upload ",
              code("matched_particles.csv"), "."),
            hr(),
            fileInput("upload_matched", "matched_particles.csv (required)",
                      accept = ".csv"),
            fileInput("upload_unmatched_ftir", "unmatched_ftir.csv (optional)",
                      accept = ".csv"),
            fileInput("upload_unmatched_raman", "unmatched_raman.csv (optional)",
                      accept = ".csv"),
            fileInput("upload_transform", "transform_params.txt (optional)",
                      accept = ".txt"),
            hr(),
            uiOutput("upload_status")
          )
        )
      )
    )
  )
)


# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {

  # ------------------------------------------------------------------
  # Load data (graceful when no pipeline output exists)
  # ------------------------------------------------------------------
  run_info <- find_latest_run()
  if (is.null(run_info)) {
    message("[Particle Viewer] No pipeline output found in ../output/")
    message("[Particle Viewer] Working directory: ", getwd())
    message("[Particle Viewer] Use the 'Upload Data' tab to load CSV files manually.")
  } else {
    message("[Particle Viewer] Loading data from: ", run_info$dir,
            " (format: ", run_info$format, ")")
  }

  # Reactive value that can be updated by CSV uploads
  uploaded_data <- reactiveVal(NULL)

  run_data <- reactive({
    # User uploads take priority
    ud <- uploaded_data()
    if (!is.null(ud)) return(ud)
    # Then try auto-detected pipeline output
    if (is.null(run_info)) return(list())
    load_run_data(run_info)
  })

  has_data <- reactive({
    d <- run_data()
    !is.null(d$matched) || !is.null(d$unmatched_ftir) || !is.null(d$unmatched_raman)
  })

  instrument_dfs <- reactive({
    if (!has_data()) return(list(ftir = NULL, raman = NULL, ldir = NULL))
    build_instrument_dfs(run_data())
  })

  # ------------------------------------------------------------------
  # Handle CSV uploads (fallback)
  # ------------------------------------------------------------------
  observeEvent(input$upload_matched, {
    matched_path <- input$upload_matched$datapath
    uf_path <- if (!is.null(input$upload_unmatched_ftir))
                 input$upload_unmatched_ftir$datapath else NULL
    ur_path <- if (!is.null(input$upload_unmatched_raman))
                 input$upload_unmatched_raman$datapath else NULL
    tp_path <- if (!is.null(input$upload_transform))
                 input$upload_transform$datapath else NULL

    tryCatch({
      data <- load_uploaded_data(matched_path, uf_path, ur_path, tp_path)
      uploaded_data(data)
      message("[Particle Viewer] Loaded uploaded data: ",
              nrow(data$matched), " matched particles")
    }, error = function(e) {
      message("[Particle Viewer] Upload error: ", conditionMessage(e))
    })
  })

  output$upload_status <- renderUI({
    if (has_data()) {
      d <- run_data()
      n_m <- if (!is.null(d$matched)) nrow(d$matched) else 0
      n_uf <- if (!is.null(d$unmatched_ftir)) nrow(d$unmatched_ftir) else 0
      n_ur <- if (!is.null(d$unmatched_raman)) nrow(d$unmatched_raman) else 0
      has_t <- !is.null(d$transform)
      tags$div(class = "alert alert-success",
        tags$b("Data loaded successfully"),
        tags$ul(
          tags$li(paste0(n_m, " matched particles")),
          tags$li(paste0(n_uf, " unmatched FTIR")),
          tags$li(paste0(n_ur, " unmatched Raman")),
          tags$li(paste0("Transform: ", if (has_t) "available" else "not loaded"))
        )
      )
    } else {
      src <- if (!is.null(uploaded_data())) "uploaded" else "auto-detected"
      tags$div(class = "alert alert-warning",
        tags$b("No data loaded"),
        tags$p(paste0("No pipeline output was ", src, ". Upload ",
                      code("matched_particles.csv"), " to get started."))
      )
    }
  })

  # Full transform matrix (FTIR original -> Raman coords), or NULL
  M_full <- reactive({
    tr <- run_data()$transform
    if (is.null(tr) || is.null(tr$M)) return(NULL)
    build_full_transform(tr)
  })

  # FTIR image bounds in original coordinates.
  # Priority: 1) from transform_params.txt, 2) from image dims, 3) from particles.
  ftir_img_bounds <- reactive({
    # Try saved scan bounds from pipeline output
    tr <- run_data()$transform
    if (!is.null(tr$ftir_scan_bounds)) return(tr$ftir_scan_bounds)

    # Estimate from image dimensions (grid geometry)
    raw_ftir_img <- ftir_raw_image()
    dfs <- instrument_dfs()
    if (!is.null(raw_ftir_img)) {
      px <- if (!is.null(dfs$ftir) && nrow(dfs$ftir) > 0) dfs$ftir$x_orig else NULL
      py <- if (!is.null(dfs$ftir) && nrow(dfs$ftir) > 0) dfs$ftir$y_orig else NULL
      return(estimate_ftir_scan_bounds(raw_ftir_img, px, py))
    }

    # Fallback: round up particle coords to nearest 500 µm
    if (!is.null(dfs$ftir) && nrow(dfs$ftir) > 0) {
      return(list(xmin = 0,
                  xmax = ceiling(max(dfs$ftir$x_orig, na.rm = TRUE) / 500) * 500,
                  ymin = 0,
                  ymax = ceiling(max(dfs$ftir$y_orig, na.rm = TRUE) / 500) * 500))
    }
    NULL
  })

  # ------------------------------------------------------------------
  # Raw image rasters
  # ------------------------------------------------------------------
  ftir_raw_image    <- reactiveVal(NULL)   # FTIR "Average Abs" image
  raman_tab_image   <- reactiveVal(NULL)   # Raman tab: user-uploaded image only
  overlay_raw_image <- reactiveVal(NULL)   # Overlay tab: raman_resized.jpg (auto-loaded)
  ldir_raw_image    <- reactiveVal(NULL)   # LDIR particle map image

  # FTIR tab: raw image placed at native FTIR scan bounds — no transform needed.
  # Particles on the FTIR tab are shown at x_orig/y_orig (FTIR instrument frame),
  # so the image just needs to sit at [xmin, xmax] × [ymin, ymax] in that same frame.
  # User fine-tuning offsets shift the image position without changing its extent.
  ftir_native_image_info <- reactive({
    raw <- ftir_raw_image()
    if (is.null(raw)) return(NULL)
    b <- ftir_img_bounds()
    if (is.null(b)) return(NULL)
    ox <- if (!is.null(input$ftir_img_offset_x)) input$ftir_img_offset_x else 0
    oy <- if (!is.null(input$ftir_img_offset_y)) input$ftir_img_offset_y else 0
    list(raster = raw, xmin = b$xmin + ox, xmax = b$xmax + ox,
         ymin = b$ymin + oy, ymax = b$ymax + oy)
  })

  # Raman tab: user-uploaded image placed at native Raman particle bounds.
  # No image is auto-loaded for the Raman tab (raman_resized.jpg is for overlay only).
  raman_native_image_info <- reactive({
    raw <- raman_tab_image()
    if (is.null(raw)) return(NULL)
    raman_df <- instrument_dfs()$raman
    if (!is.null(raman_df) && nrow(raman_df) > 0) {
      pad <- 300
      ox <- if (!is.null(input$raman_img_offset_x)) input$raman_img_offset_x else 0
      oy <- if (!is.null(input$raman_img_offset_y)) input$raman_img_offset_y else 0
      return(list(raster = raw,
                  xmin = min(raman_df$x_orig, na.rm = TRUE) - pad + ox,
                  xmax = max(raman_df$x_orig, na.rm = TRUE) + pad + ox,
                  ymin = min(raman_df$y_orig, na.rm = TRUE) - pad + oy,
                  ymax = max(raman_df$y_orig, na.rm = TRUE) + pad + oy))
    }
    NULL
  })

  # Overlay tab: raman_resized.jpg placed at Raman-normalized bounds.
  # The image covers the full filter.  We compute bounds from FTIR + Raman
  # particles only (not LDIR, whose alignment is coarser and would stretch
  # the image far beyond the filter area).
  overlay_image_info <- reactive({
    raw <- overlay_raw_image()
    if (is.null(raw)) return(NULL)
    dfs <- instrument_dfs()
    # Use FTIR + Raman only for image placement (well-aligned instruments)
    all_x <- c(dfs$ftir$x, dfs$raman$x)
    all_y <- c(dfs$ftir$y, dfs$raman$y)
    all_x <- all_x[is.finite(all_x)]
    all_y <- all_y[is.finite(all_y)]
    if (length(all_x) == 0) return(NULL)

    # Image aspect ratio (width / height)
    img_aspect <- ncol(raw) / nrow(raw)

    # Compute bounds from actual particle extent with 20% padding
    pad_frac <- 0.20
    rx <- range(all_x); ry <- range(all_y)
    w <- diff(rx); h <- diff(ry)
    xmin <- rx[1] - pad_frac * w
    xmax <- rx[2] + pad_frac * w
    ymin <- ry[1] - pad_frac * h
    ymax <- ry[2] + pad_frac * h

    # Enforce image aspect ratio (expand the tighter dimension)
    curr_aspect <- (xmax - xmin) / (ymax - ymin)
    if (curr_aspect < img_aspect) {
      mid_x <- (xmin + xmax) / 2
      half_x <- (ymax - ymin) * img_aspect / 2
      xmin <- mid_x - half_x; xmax <- mid_x + half_x
    } else {
      mid_y <- (ymin + ymax) / 2
      half_y <- (xmax - xmin) / img_aspect / 2
      ymin <- mid_y - half_y; ymax <- mid_y + half_y
    }

    # Apply user fine-tuning offsets
    ox <- if (!is.null(input$overlay_img_offset_x)) input$overlay_img_offset_x else 0
    oy <- if (!is.null(input$overlay_img_offset_y)) input$overlay_img_offset_y else 0

    list(raster = raw,
         xmin = xmin + ox, xmax = xmax + ox,
         ymin = ymin + oy, ymax = ymax + oy)
  })

  # LDIR tab: image placed at LDIR scan area bounds.
  # LDIR coordinates use image convention: y=0 at top, increasing downward.
  # ggplot uses Cartesian: y increases upward.  annotation_raster always
  # places row 1 at ymax (top of plot) regardless of ymin/ymax ordering.
  # So we PHYSICALLY flip the raster (reverse rows) so that what was the
  # image bottom (high y_um) becomes row 1 and is placed at ymax (top of
  # the plot), matching the particles plotted at y = y_um.
  ldir_native_image_info <- reactive({
    raw <- ldir_raw_image()
    if (is.null(raw)) return(NULL)
    ldir_df <- instrument_dfs()$ldir
    if (!is.null(ldir_df) && nrow(ldir_df) > 0 &&
        any(!is.na(ldir_df$x_orig))) {
      ox <- if (!is.null(input$ldir_img_offset_x)) input$ldir_img_offset_x else 0
      oy <- if (!is.null(input$ldir_img_offset_y)) input$ldir_img_offset_y else 0
      xvals <- ldir_df$x_orig[!is.na(ldir_df$x_orig)]
      yvals <- ldir_df$y_orig[!is.na(ldir_df$y_orig)]
      # Compute scan extent: round up max coordinate to nearest 1000 µm
      extent <- max(ceiling(max(xvals) / 1000) * 1000,
                    ceiling(max(yvals) / 1000) * 1000)
      # Flip raster vertically to match LDIR y convention
      flipped <- raw[nrow(raw):1, , , drop = FALSE]
      return(list(raster = flipped,
                  xmin = 0 + ox, xmax = extent + ox,
                  ymin = 0 + oy, ymax = extent + oy))
    }
    NULL
  })

  # Auto-load default images from project root
  observe({
    default_ftir <- file.path("..", "Average Abs.( Comparstic Spotlight F2Ba Au 240926 ).png")
    if (!file.exists(default_ftir)) return()
    raw <- load_image_raster(default_ftir)
    if (!is.null(raw)) ftir_raw_image(raw)
  })

  observe({
    default_raman <- file.path("..", "raman_resized.jpg")
    if (!file.exists(default_raman)) return()
    raw <- load_image_raster(default_raman)
    if (!is.null(raw)) overlay_raw_image(raw)
  })

  observe({
    default_ldir <- file.path("..", "Comparstic LDIR F2Ba_G3B AU 240925.png")
    if (!file.exists(default_ldir)) return()
    raw <- load_image_raster(default_ldir)
    if (!is.null(raw)) ldir_raw_image(raw)
  })

  # Handle uploaded images
  observeEvent(input$ftir_image_upload, {
    raw <- load_image_raster(input$ftir_image_upload$datapath)
    if (!is.null(raw)) ftir_raw_image(raw)
  })

  observeEvent(input$raman_image_upload, {
    raw <- load_image_raster(input$raman_image_upload$datapath)
    if (!is.null(raw)) raman_tab_image(raw)
  })

  observeEvent(input$overlay_image_upload, {
    raw <- load_image_raster(input$overlay_image_upload$datapath)
    if (!is.null(raw)) overlay_raw_image(raw)
  })

  observeEvent(input$ldir_image_upload, {
    raw <- load_image_raster(input$ldir_image_upload$datapath)
    if (!is.null(raw)) ldir_raw_image(raw)
  })

  # ------------------------------------------------------------------
  # Update filter controls from data
  # ------------------------------------------------------------------
  observe({
    dfs <- instrument_dfs()

    if (!is.null(dfs$ftir) && nrow(dfs$ftir) > 0) {
      ftir <- dfs$ftir
      updateSelectInput(session, "ftir_material_filter",
                        choices = c("All", sort(unique(ftir$material))), selected = "All")
      q_range <- range(ftir$quality, na.rm = TRUE)
      updateSliderInput(session, "ftir_quality_range",
                        min = floor(q_range[1] * 100) / 100,
                        max = ceiling(q_range[2] * 100) / 100,
                        value = c(floor(q_range[1] * 100) / 100,
                                  ceiling(q_range[2] * 100) / 100))
      s_max <- ceiling(max(ftir$feret_max, na.rm = TRUE) / 10) * 10
      updateSliderInput(session, "ftir_size_range", min = 0, max = s_max,
                        value = c(0, s_max))
    }

    if (!is.null(dfs$raman) && nrow(dfs$raman) > 0) {
      raman <- dfs$raman
      updateSelectInput(session, "raman_material_filter",
                        choices = c("All", sort(unique(raman$material))), selected = "All")
      q_range <- range(raman$quality, na.rm = TRUE)
      updateSliderInput(session, "raman_quality_range",
                        min = floor(q_range[1]), max = ceiling(q_range[2]),
                        value = c(floor(q_range[1]), ceiling(q_range[2])))
      s_max <- ceiling(max(raman$feret_max, na.rm = TRUE) / 10) * 10
      updateSliderInput(session, "raman_size_range", min = 0, max = s_max,
                        value = c(0, s_max))
    }

    if (!is.null(dfs$ldir) && nrow(dfs$ldir) > 0) {
      ldir <- dfs$ldir
      updateSelectInput(session, "ldir_material_filter",
                        choices = c("All", sort(unique(ldir$material))), selected = "All")
      q_range <- range(ldir$quality, na.rm = TRUE)
      if (all(is.finite(q_range))) {
        updateSliderInput(session, "ldir_quality_range",
                          min = floor(q_range[1] * 100) / 100,
                          max = ceiling(q_range[2] * 100) / 100,
                          value = c(floor(q_range[1] * 100) / 100,
                                    ceiling(q_range[2] * 100) / 100))
      }
      s_max <- ceiling(max(ldir$feret_max, na.rm = TRUE) / 10) * 10
      if (is.finite(s_max)) {
        updateSliderInput(session, "ldir_size_range", min = 0, max = s_max,
                          value = c(0, s_max))
      }
    }

    all_mats <- sort(unique(c(dfs$ftir$material, dfs$raman$material, dfs$ldir$material)))
    updateSelectInput(session, "overlay_material_filter",
                      choices = c("All", all_mats), selected = "All")

    if (!is.null(run_data()$matched)) {
      max_dist <- ceiling(max(run_data()$matched$match_distance, na.rm = TRUE))
      updateSliderInput(session, "overlay_dist_range",
                        min = 0, max = max_dist, value = c(0, max_dist))
    }
  })

  # ==================================================================
  # Helper: generic instrument filter
  # ==================================================================
  filter_instrument <- function(df, quality_range, size_range, mat_filter,
                                 match_filter) {
    df <- df[!is.na(df$quality) &
             df$quality >= quality_range[1] &
             df$quality <= quality_range[2], ]
    df <- df[!is.na(df$feret_max) &
             df$feret_max >= size_range[1] &
             df$feret_max <= size_range[2], ]
    if (!("All" %in% mat_filter))
      df <- df[df$material %in% mat_filter, ]
    df <- df[df$match_status %in% match_filter, ]
    df
  }

  # ==================================================================
  # Helper: add image background to a ggplot
  # ==================================================================
  add_image_bg <- function(p, img_info, alpha = 0.4) {
    if (is.null(img_info)) return(p)
    p + annotation_raster(img_info$raster,
          xmin = img_info$xmin, xmax = img_info$xmax,
          ymin = img_info$ymin, ymax = img_info$ymax,
          interpolate = TRUE)
  }

  # ==================================================================
  # Helper: ggplot scatter with optional image background
  # ==================================================================
  make_scatter <- function(df, img_info, bounds, title,
                            match_colours = NULL) {

    p <- ggplot(df, aes(x = x, y = y))

    # Background image (with per-image bounds)
    p <- add_image_bg(p, img_info)

    # Points
    p <- p + geom_point(aes(colour = match_status, size = feret_max),
                         alpha = 0.7)

    if (!is.null(match_colours))
      p <- p + scale_colour_manual(values = match_colours)

    p <- p +
      scale_size_continuous(range = c(2, 12), guide = "none") +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = title, x = "X (\u00b5m)", y = "Y (\u00b5m)") +
      theme_minimal(base_size = 13) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "grey98", colour = NA),
        panel.grid       = element_line(colour = "grey90"),
        legend.position  = "bottom"
      )

    p
  }

  # ==================================================================
  # Helper: detail table HTML for single instrument
  # ==================================================================
  single_detail_html <- function(row, instrument_name, quality_label) {
    if (is.null(row)) {
      return(tags$p(class = "text-muted", "Hover over a particle to see details"))
    }
    tags$table(class = "hover-tbl",
      tags$tr(tags$th("Field"), tags$th("Value")),
      make_detail_row("Instrument", instrument_name),
      make_detail_row("Particle ID", row$particle_id),
      make_detail_row("Material", tags$b(row$material)),
      make_detail_row(quality_label, round(row$quality, 3)),
      make_detail_row("Feret Max", paste0(round(row$feret_max, 1), " \u00b5m")),
      make_detail_row("Area", paste0(round(row$area_um2, 1), " \u00b5m\u00b2")),
      make_detail_row("Major Dim", paste0(round(row$major_um, 1), " \u00b5m")),
      make_detail_row("Minor Dim", paste0(round(row$minor_um, 1), " \u00b5m")),
      make_detail_row("Position (native)",
                      paste0("(", round(row$x_orig, 1), ", ", round(row$y_orig, 1), ")")),
      make_detail_row("Position (aligned)",
                      paste0("(", round(row$x, 1), ", ", round(row$y, 1), ")")),
      make_detail_row("Match Status", row$match_status)
    )
  }


  # ==================================================================
  # Sticky hover state: stores the last successfully found particle/row
  # so the info panel doesn't flicker when the cursor drifts slightly.
  # Updated only when a NEW particle is found; keeps showing the last
  # particle when hovering over empty background.
  # ==================================================================
  last_hover <- reactiveValues(ftir = NULL, raman = NULL, ldir = NULL, overlay = NULL)

  # ==================================================================
  # Zoom state: NULL means full view, otherwise list(x=c(lo,hi), y=c(lo,hi))
  # ==================================================================
  zoom <- reactiveValues(ftir = NULL, raman = NULL, ldir = NULL, overlay = NULL)

  observeEvent(input$ftir_brush, {
    b <- input$ftir_brush
    zoom$ftir <- list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax))
  })
  observeEvent(input$ftir_dblclick, { zoom$ftir <- NULL })

  observeEvent(input$raman_brush, {
    b <- input$raman_brush
    zoom$raman <- list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax))
  })
  observeEvent(input$raman_dblclick, { zoom$raman <- NULL })

  observeEvent(input$ldir_brush, {
    b <- input$ldir_brush
    zoom$ldir <- list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax))
  })
  observeEvent(input$ldir_dblclick, { zoom$ldir <- NULL })

  observeEvent(input$overlay_brush, {
    b <- input$overlay_brush
    zoom$overlay <- list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax))
  })
  observeEvent(input$overlay_dblclick, { zoom$overlay <- NULL })

  # ==================================================================
  # FTIR TAB
  # ==================================================================

  ftir_filtered <- reactive({
    dfs <- instrument_dfs()
    if (is.null(dfs$ftir) || nrow(dfs$ftir) == 0) return(data.frame())
    filter_instrument(dfs$ftir, input$ftir_quality_range, input$ftir_size_range,
                      input$ftir_material_filter, input$ftir_match_filter)
  })

  output$ftir_plot <- renderPlot({
    df <- ftir_filtered()
    # Display in native FTIR instrument frame (x_orig, y_orig)
    df_disp <- df
    if (nrow(df_disp) > 0) { df_disp$x <- df_disp$x_orig; df_disp$y <- df_disp$y_orig }

    bounds <- if (!is.null(zoom$ftir)) zoom$ftir else {
      b <- ftir_img_bounds()
      if (!is.null(b)) list(x = c(b$xmin - 200, b$xmax + 200),
                            y = c(b$ymin - 200, b$ymax + 200))
      else compute_bounds(df_disp, NULL)
    }

    img <- ftir_native_image_info()

    if (nrow(df_disp) == 0) {
      p <- ggplot() + coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = "FTIR — no particles loaded", x = "X (\u00b5m)", y = "Y (\u00b5m)") +
        theme_minimal(base_size = 13) +
        theme(plot.background = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA))
      return(add_image_bg(p, img))
    }
    make_scatter(df_disp, img, bounds,
                 paste0("FTIR Particles (", nrow(df_disp), " shown)"),
                 match_colours = c(matched = "#2ca02c", unmatched = "#d62728"))
  })

  output$ftir_summary_text <- renderText({
    df <- ftir_filtered()
    if (nrow(df) == 0) return("No pipeline data loaded")
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched | ",
           length(unique(df$material)), " materials")
  })

  observeEvent(input$ftir_hover, {
    hover <- input$ftir_hover
    if (is.null(hover)) return()
    df <- ftir_filtered()
    if (nrow(df) == 0) return()
    # Search in native FTIR coordinate space (x_orig, y_orig)
    dists <- sqrt((df$x_orig - hover$x)^2 + (df$y_orig - hover$y)^2)
    idx   <- which.min(dists)
    threshold <- max(diff(range(df$x_orig, na.rm = TRUE)),
                     diff(range(df$y_orig, na.rm = TRUE)), 500) * 0.05
    if (dists[idx] <= threshold) last_hover$ftir <- df[idx, , drop = FALSE]
  })

  output$ftir_hover_info <- renderUI({
    row <- last_hover$ftir
    single_detail_html(row, "FTIR", "AAU Quality")
  })


  # ==================================================================
  # RAMAN TAB
  # ==================================================================

  raman_filtered <- reactive({
    dfs <- instrument_dfs()
    if (is.null(dfs$raman) || nrow(dfs$raman) == 0) return(data.frame())
    filter_instrument(dfs$raman, input$raman_quality_range, input$raman_size_range,
                      input$raman_material_filter, input$raman_match_filter)
  })

  output$raman_plot <- renderPlot({
    df <- raman_filtered()
    # Display in native Raman instrument frame (x_orig, y_orig)
    df_disp <- df
    if (nrow(df_disp) > 0) { df_disp$x <- df_disp$x_orig; df_disp$y <- df_disp$y_orig }

    bounds <- if (!is.null(zoom$raman)) zoom$raman else {
      if (nrow(df_disp) > 0) {
        pad <- 300
        list(x = c(min(df_disp$x, na.rm = TRUE) - pad, max(df_disp$x, na.rm = TRUE) + pad),
             y = c(min(df_disp$y, na.rm = TRUE) - pad, max(df_disp$y, na.rm = TRUE) + pad))
      } else list(x = c(-1000, 1000), y = c(-1000, 1000))
    }

    img <- raman_native_image_info()

    if (nrow(df_disp) == 0) {
      p <- ggplot() + coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = "Raman — no particles loaded", x = "X (\u00b5m)", y = "Y (\u00b5m)") +
        theme_minimal(base_size = 13) +
        theme(plot.background = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA))
      return(add_image_bg(p, img))
    }
    make_scatter(df_disp, img, bounds,
                 paste0("Raman Particles (", nrow(df_disp), " shown)"),
                 match_colours = c(matched = "#1f77b4", unmatched = "#ff7f0e"))
  })

  output$raman_summary_text <- renderText({
    df <- raman_filtered()
    if (nrow(df) == 0) return("No pipeline data loaded")
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched | ",
           length(unique(df$material)), " materials")
  })

  observeEvent(input$raman_hover, {
    hover <- input$raman_hover
    if (is.null(hover)) return()
    df <- raman_filtered()
    if (nrow(df) == 0) return()
    # Search in native Raman coordinate space (x_orig, y_orig)
    dists <- sqrt((df$x_orig - hover$x)^2 + (df$y_orig - hover$y)^2)
    idx   <- which.min(dists)
    threshold <- max(diff(range(df$x_orig, na.rm = TRUE)),
                     diff(range(df$y_orig, na.rm = TRUE)), 500) * 0.05
    if (dists[idx] <= threshold) last_hover$raman <- df[idx, , drop = FALSE]
  })

  output$raman_hover_info <- renderUI({
    row <- last_hover$raman
    single_detail_html(row, "Raman", "HQI")
  })


  # ==================================================================
  # LDIR TAB
  # ==================================================================

  ldir_filtered <- reactive({
    dfs <- instrument_dfs()
    if (is.null(dfs$ldir) || nrow(dfs$ldir) == 0) return(data.frame())
    filter_instrument(dfs$ldir, input$ldir_quality_range, input$ldir_size_range,
                      input$ldir_material_filter, input$ldir_match_filter)
  })

  output$ldir_plot <- renderPlot({
    df <- ldir_filtered()
    # Display in native LDIR frame (x_orig, y_orig)
    df_disp <- df
    if (nrow(df_disp) > 0) { df_disp$x <- df_disp$x_orig; df_disp$y <- df_disp$y_orig }

    bounds <- if (!is.null(zoom$ldir)) zoom$ldir else {
      if (nrow(df_disp) > 0) {
        pad <- 500
        list(x = c(min(df_disp$x, na.rm = TRUE) - pad, max(df_disp$x, na.rm = TRUE) + pad),
             y = c(min(df_disp$y, na.rm = TRUE) - pad, max(df_disp$y, na.rm = TRUE) + pad))
      } else list(x = c(-1000, 14000), y = c(-1000, 14000))
    }

    img <- ldir_native_image_info()

    if (nrow(df_disp) == 0) {
      p <- ggplot() + coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = "LDIR — no particles loaded", x = "X (\u00b5m)", y = "Y (\u00b5m)") +
        theme_minimal(base_size = 13) +
        theme(plot.background = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA))
      return(add_image_bg(p, img))
    }
    make_scatter(df_disp, img, bounds,
                 paste0("LDIR Particles (", nrow(df_disp), " shown)"),
                 match_colours = c(matched = "#17becf", unmatched = "#bcbd22"))
  })

  output$ldir_summary_text <- renderText({
    df <- ldir_filtered()
    if (nrow(df) == 0) return("No LDIR data loaded")
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched | ",
           length(unique(df$material)), " materials")
  })

  observeEvent(input$ldir_hover, {
    hover <- input$ldir_hover
    if (is.null(hover)) return()
    df <- ldir_filtered()
    if (nrow(df) == 0) return()
    # Search in native LDIR coordinate space (x_orig, y_orig)
    dists <- sqrt((df$x_orig - hover$x)^2 + (df$y_orig - hover$y)^2)
    idx   <- which.min(dists)
    threshold <- max(diff(range(df$x_orig, na.rm = TRUE)),
                     diff(range(df$y_orig, na.rm = TRUE)), 500) * 0.05
    if (dists[idx] <= threshold) last_hover$ldir <- df[idx, , drop = FALSE]
  })

  output$ldir_hover_info <- renderUI({
    row <- last_hover$ldir
    single_detail_html(row, "LDIR", "Quality")
  })


  # ==================================================================
  # OVERLAY TAB
  # ==================================================================

  overlay_matched <- reactive({
    d <- run_data()
    if (is.null(d$matched) || nrow(d$matched) == 0) return(data.frame())
    df <- d$matched

    df <- df[!is.na(df$raman_quality) &
             df$raman_quality >= input$overlay_hqi_range[1] &
             df$raman_quality <= input$overlay_hqi_range[2], ]
    df <- df[!is.na(df$ftir_quality) &
             df$ftir_quality >= input$overlay_quality_range[1] &
             df$ftir_quality <= input$overlay_quality_range[2], ]

    ftir_ok <- !is.na(df$ftir_feret_max_um) &
               df$ftir_feret_max_um >= input$overlay_size_range[1] &
               df$ftir_feret_max_um <= input$overlay_size_range[2]
    raman_ok <- !is.na(df$raman_feret_max_um) &
                df$raman_feret_max_um >= input$overlay_size_range[1] &
                df$raman_feret_max_um <= input$overlay_size_range[2]
    df <- df[ftir_ok | raman_ok, ]

    df <- df[!is.na(df$match_distance) &
             df$match_distance >= input$overlay_dist_range[1] &
             df$match_distance <= input$overlay_dist_range[2], ]

    if (!("All" %in% input$overlay_material_filter)) {
      df <- df[df$ftir_material %in% input$overlay_material_filter |
               df$raman_material %in% input$overlay_material_filter, ]
    }
    df
  })

  # LDIR-Raman matched data for overlay
  overlay_ldir_matched <- reactive({
    d <- run_data()
    if (is.null(d$ldir_raman_matched) || nrow(d$ldir_raman_matched) == 0)
      return(data.frame())
    d$ldir_raman_matched
  })

  output$overlay_plot <- renderPlot({
    dfs <- instrument_dfs()
    matched <- overlay_matched()
    ldir_m <- overlay_ldir_matched()
    layers <- input$overlay_layers
    bounds <- if (!is.null(zoom$overlay)) zoom$overlay
              else compute_bounds(dfs$ftir, dfs$raman, dfs$ldir)

    p <- ggplot() +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = "FTIR + Raman + LDIR Overlay (aligned coordinates)",
           x = "X (\u00b5m)", y = "Y (\u00b5m)") +
      theme_minimal(base_size = 13) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "grey98", colour = NA),
        panel.grid       = element_line(colour = "grey90"),
        legend.position  = "bottom"
      )

    # Background image: raman_resized.jpg placed at Raman-normalized bounds
    p <- add_image_bg(p, overlay_image_info())

    # FTIR-Raman match lines
    if ("match_lines" %in% layers && nrow(matched) > 0) {
      seg_df <- data.frame(
        x    = matched$ftir_x_aligned, y    = matched$ftir_y_aligned,
        xend = matched$raman_x_norm,     yend = matched$raman_y_norm
      )
      p <- p + geom_segment(data = seg_df,
                              aes(x = x, y = y, xend = xend, yend = yend),
                              colour = "grey60", alpha = 0.3, linewidth = 0.3)
    }

    # LDIR-Raman match lines
    if ("ldir_lines" %in% layers && nrow(ldir_m) > 0 &&
        "ldir_x_aligned" %in% names(ldir_m) && "raman_x_norm" %in% names(ldir_m)) {
      ldir_seg <- data.frame(
        x    = ldir_m$ldir_x_aligned, y    = ldir_m$ldir_y_aligned,
        xend = ldir_m$raman_x_norm,    yend = ldir_m$raman_y_norm
      )
      p <- p + geom_segment(data = ldir_seg,
                              aes(x = x, y = y, xend = xend, yend = yend),
                              colour = "#9467bd", alpha = 0.3, linewidth = 0.3)
    }

    # Unmatched FTIR
    if ("unmatched_ftir" %in% layers && !is.null(dfs$ftir) && nrow(dfs$ftir) > 0) {
      um_f <- dfs$ftir[dfs$ftir$match_status == "unmatched", ]
      if (nrow(um_f) > 0) {
        p <- p + geom_point(data = um_f, aes(x = x, y = y, size = feret_max),
                             colour = "#d62728", alpha = 0.4, shape = 4) +
          scale_size_continuous(range = c(2, 12), guide = "none")
      }
    }

    # Unmatched Raman
    if ("unmatched_raman" %in% layers && !is.null(dfs$raman) && nrow(dfs$raman) > 0) {
      um_r <- dfs$raman[dfs$raman$match_status == "unmatched", ]
      if (nrow(um_r) > 0) {
        p <- p + geom_point(data = um_r, aes(x = x, y = y, size = feret_max),
                             colour = "#ff7f0e", alpha = 0.4, shape = 5)
      }
    }

    # Unmatched LDIR
    if ("ldir_unmatched" %in% layers && !is.null(dfs$ldir) && nrow(dfs$ldir) > 0) {
      um_l <- dfs$ldir[dfs$ldir$match_status == "unmatched", ]
      if (nrow(um_l) > 0) {
        p <- p + geom_point(data = um_l, aes(x = x, y = y, size = feret_max),
                             colour = "#9467bd", alpha = 0.3, shape = 8)
      }
    }

    # Matched pairs: FTIR as triangle, Raman as circle, LDIR as diamond
    all_pts <- list()
    if ("matched" %in% layers && nrow(matched) > 0) {
      all_pts[[length(all_pts) + 1]] <- data.frame(
        x = matched$ftir_x_aligned, y = matched$ftir_y_aligned,
        feret_max = matched$ftir_feret_max_um, instrument = "FTIR"
      )
      all_pts[[length(all_pts) + 1]] <- data.frame(
        x = matched$raman_x_norm, y = matched$raman_y_norm,
        feret_max = matched$raman_feret_max_um, instrument = "Raman"
      )
    }
    if ("ldir_matched" %in% layers && nrow(ldir_m) > 0 &&
        "ldir_x_aligned" %in% names(ldir_m)) {
      all_pts[[length(all_pts) + 1]] <- data.frame(
        x = ldir_m$ldir_x_aligned, y = ldir_m$ldir_y_aligned,
        feret_max = ldir_m$ldir_feret_max_um, instrument = "LDIR"
      )
    }
    if (length(all_pts) > 0) {
      both <- do.call(rbind, all_pts)
      p <- p + geom_point(data = both,
                            aes(x = x, y = y, size = feret_max,
                                shape = instrument, colour = instrument),
                            alpha = 0.7) +
        scale_shape_manual(values = c(FTIR = 17, Raman = 16, LDIR = 18)) +
        scale_colour_manual(values = c(FTIR = "#2ca02c", Raman = "#1f77b4",
                                        LDIR = "#9467bd"))
    }

    p
  })

  output$overlay_summary_text <- renderText({
    m <- overlay_matched()
    dfs <- instrument_dfs()
    if (nrow(m) == 0 && is.null(dfs$ftir)) return("No pipeline data loaded")
    n_um_f <- if (!is.null(dfs$ftir)) sum(dfs$ftir$match_status == "unmatched") else 0
    n_um_r <- if (!is.null(dfs$raman)) sum(dfs$raman$match_status == "unmatched") else 0
    n_ldir <- if (!is.null(dfs$ldir)) nrow(dfs$ldir) else 0
    n_ldir_m <- if (!is.null(dfs$ldir)) sum(dfs$ldir$match_status == "matched") else 0
    paste0(nrow(m), " FTIR-Raman pairs | ",
           n_um_f, " unmatched FTIR | ",
           n_um_r, " unmatched Raman | ",
           n_ldir_m, "/", n_ldir, " LDIR matched")
  })

  # Overlay: sticky hover — update last_hover$overlay only when a new match is found
  observeEvent(input$overlay_hover, {
    hover <- input$overlay_hover
    if (is.null(hover)) return()

    matched <- overlay_matched()
    ldir_m <- overlay_ldir_matched()

    # Snap threshold: 5% of visible range (zoom-aware)
    dfs <- instrument_dfs()
    vis <- if (!is.null(zoom$overlay)) zoom$overlay
           else compute_bounds(dfs$ftir, dfs$raman, dfs$ldir)
    snap_dist <- max(diff(vis$x), diff(vis$y), 500) * 0.05

    best_dist <- Inf
    best_row <- NULL
    best_source <- NULL

    # Check FTIR-Raman matched
    if (!is.null(matched) && nrow(matched) > 0) {
      dist_f <- sqrt((matched$ftir_x_aligned - hover$x)^2 +
                      (matched$ftir_y_aligned - hover$y)^2)
      dist_r <- sqrt((matched$raman_x_norm - hover$x)^2 +
                      (matched$raman_y_norm - hover$y)^2)
      d <- pmin(dist_f, dist_r)
      idx <- which.min(d)
      if (d[idx] < best_dist) {
        best_dist <- d[idx]
        best_row <- matched[idx, ]
        best_source <- "ftir_raman"
      }
    }

    # Check LDIR-Raman matched
    if (nrow(ldir_m) > 0 && "ldir_x_aligned" %in% names(ldir_m)) {
      dist_l <- sqrt((ldir_m$ldir_x_aligned - hover$x)^2 +
                      (ldir_m$ldir_y_aligned - hover$y)^2)
      idx_l <- which.min(dist_l)
      if (dist_l[idx_l] < best_dist) {
        best_dist <- dist_l[idx_l]
        best_row <- ldir_m[idx_l, ]
        best_source <- "ldir_raman"
      }
    }

    if (best_dist <= snap_dist && !is.null(best_row)) {
      last_hover$overlay <- best_row
      attr(last_hover$overlay, "source") <- best_source
    }
  })

  output$overlay_hover_info <- renderUI({
    row <- last_hover$overlay
    if (is.null(row)) {
      return(tags$p(class = "text-muted",
                    "Hover over a matched particle for the comparison table"))
    }

    src <- attr(row, "source")
    if (!is.null(src) && src == "ldir_raman") {
      return(tags$table(class = "hover-tbl",
        tags$tr(tags$th(""), tags$th("LDIR"), tags$th("Raman")),
        tags$tr(tags$td(tags$b("Particle ID")),
                tags$td(row$ldir_particle_id),
                tags$td(row$raman_particle_id)),
        tags$tr(tags$td(tags$b("Material")),
                tags$td(row$ldir_material),
                tags$td(row$raman_material)),
        tags$tr(tags$td(tags$b("Quality")),
                tags$td(round(row$ldir_quality, 3)),
                tags$td(paste0("HQI ", round(row$raman_quality, 2)))),
        tags$tr(tags$td(tags$b("Feret Max")),
                tags$td(paste0(round(row$ldir_feret_max_um, 1), " \u00b5m")),
                tags$td(paste0(round(row$raman_feret_max_um, 1), " \u00b5m"))),
        tags$tr(tags$td(tags$b("Match Dist.")),
                tags$td(colspan = "2",
                        paste0(round(row$match_distance, 1), " \u00b5m")))
      ))
    }

    # Default: FTIR-Raman hover info
    tags$table(class = "hover-tbl",
      tags$tr(tags$th(""), tags$th("FTIR"), tags$th("Raman")),
      tags$tr(tags$td(tags$b("Particle ID")),
              tags$td(row$ftir_particle_id),
              tags$td(row$raman_particle_id)),
      tags$tr(tags$td(tags$b("Material")),
              tags$td(row$ftir_material),
              tags$td(row$raman_material)),
      tags$tr(tags$td(tags$b("Quality")),
              tags$td(paste0("AAU ", round(row$ftir_quality, 3))),
              tags$td(paste0("HQI ", round(row$raman_quality, 2)))),
      tags$tr(tags$td(tags$b("Feret Max")),
              tags$td(paste0(round(row$ftir_feret_max_um, 1), " \u00b5m")),
              tags$td(paste0(round(row$raman_feret_max_um, 1), " \u00b5m"))),
      tags$tr(tags$td(tags$b("Area")),
              tags$td(paste0(round(row$ftir_area_um2, 1), " \u00b5m\u00b2")),
              tags$td(paste0(round(row$raman_area_um2, 1), " \u00b5m\u00b2"))),
      tags$tr(tags$td(tags$b("Position")),
              tags$td(paste0("(", round(row$ftir_x_aligned, 1), ", ",
                              round(row$ftir_y_aligned, 1), ")")),
              tags$td(paste0("(", round(row$raman_x_um, 1), ", ",
                              round(row$raman_y_um, 1), ")"))),
      tags$tr(tags$td(tags$b("Match Dist.")),
              tags$td(colspan = "2",
                      paste0(round(row$match_distance, 1), " \u00b5m")))
    )
  })
}

# ============================================================================
# Run
# ============================================================================
shinyApp(ui = ui, server = server)
