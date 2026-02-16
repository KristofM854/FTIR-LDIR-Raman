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
                accept = c("image/png", "image/jpeg"))
    ),
    mainPanel(width = 9,
      plotOutput(paste0(id_prefix, "_plot"), height = "650px",
                 hover = hoverOpts(paste0(id_prefix, "_hover"), delay = 100,
                                   delayType = "throttle")),
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

    # Tab 3: LDIR (placeholder)
    tabPanel("LDIR",
      div(class = "placeholder-msg",
        h3("LDIR Panel"),
        p("Awaiting LDIR particle coordinates."),
        p("When LDIR coordinate data or an LDIR image becomes available,
           this panel will show an interactive spatial plot identical
           to the FTIR and Raman panels."),
        p("Currently, LDIR comparison is limited to bulk material composition
           (see pipeline output ", code("ldir_material_distribution.csv"), ").")
      )
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
                                         "Match lines" = "match_lines"),
                             selected = c("matched", "unmatched_ftir",
                                          "unmatched_raman"),
                             inline = FALSE),
          hr(),
          div(class = "info-box",
              h5("Match Summary"), textOutput("overlay_summary_text")),
          hr(),
          fileInput("overlay_image_upload", "Background Image",
                    accept = c("image/png", "image/jpeg"))
        ),
        mainPanel(width = 9,
          plotOutput("overlay_plot", height = "650px",
                     hover = hoverOpts("overlay_hover", delay = 100,
                                       delayType = "throttle")),
          hr(),
          div(class = "info-box",
              h5("Match Details (hover on matched particle)"),
              detail_table_ui("overlay_hover_info"))
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
  # Load data
  # ------------------------------------------------------------------
  run_dir <- find_latest_run()
  run_data <- reactive({
    req(run_dir)
    load_run_data(run_dir)
  })

  instrument_dfs <- reactive({
    req(run_data())
    build_instrument_dfs(run_data())
  })

  # ------------------------------------------------------------------
  # Image state
  # ------------------------------------------------------------------
  images <- reactiveValues(ftir = NULL, raman = NULL, overlay = NULL)

  observe({
    default_img <- file.path("..", "Average Abs.( Comparstic Spotlight F2Ba Au 240926 ).png")
    if (file.exists(default_img)) {
      img <- load_image_raster(default_img)
      images$ftir <- img
      images$overlay <- img
    }
  })

  observeEvent(input$ftir_image_upload, {
    images$ftir <- load_image_raster(input$ftir_image_upload$datapath)
  })
  observeEvent(input$raman_image_upload, {
    images$raman <- load_image_raster(input$raman_image_upload$datapath)
  })
  observeEvent(input$overlay_image_upload, {
    images$overlay <- load_image_raster(input$overlay_image_upload$datapath)
  })

  # ------------------------------------------------------------------
  # Update filter controls from data
  # ------------------------------------------------------------------
  observe({
    dfs <- instrument_dfs()

    if (!is.null(dfs$ftir)) {
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

    if (!is.null(dfs$raman)) {
      raman <- dfs$raman
      updateSelectInput(session, "raman_material_filter",
                        choices = c("All", sort(unique(raman$material))), selected = "All")
      q_range <- range(raman$quality, na.rm = TRUE)
      updateSliderInput(session, "raman_hqi_range",
                        min = floor(q_range[1]), max = ceiling(q_range[2]),
                        value = c(floor(q_range[1]), ceiling(q_range[2])))
      s_max <- ceiling(max(raman$feret_max, na.rm = TRUE) / 10) * 10
      updateSliderInput(session, "raman_size_range", min = 0, max = s_max,
                        value = c(0, s_max))
    }

    all_mats <- sort(unique(c(dfs$ftir$material, dfs$raman$material)))
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
  # Helper: ggplot scatter with optional image background
  # ==================================================================
  make_scatter <- function(df, img_raster, bounds, title,
                            colour_by = "match_status",
                            match_colours = NULL,
                            shape_by = NULL, shapes = NULL) {

    p <- ggplot(df, aes(x = x, y = y))

    # Background image
    if (!is.null(img_raster)) {
      p <- p + annotation_raster(img_raster,
               xmin = bounds$x[1], xmax = bounds$x[2],
               ymin = bounds$y[1], ymax = bounds$y[2],
               interpolate = TRUE)
    }

    # Points
    if (!is.null(shape_by)) {
      p <- p + geom_point(aes(colour = .data[[colour_by]],
                               shape = .data[[shape_by]],
                               size = feret_max),
                           alpha = 0.7)
      if (!is.null(shapes)) p <- p + scale_shape_manual(values = shapes)
    } else {
      p <- p + geom_point(aes(colour = .data[[colour_by]],
                               size = feret_max),
                           alpha = 0.7)
    }

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
      make_detail_row("Position (aligned)",
                      paste0("(", round(row$x, 1), ", ", round(row$y, 1), ")")),
      make_detail_row("Position (original)",
                      paste0("(", round(row$x_orig, 1), ", ", round(row$y_orig, 1), ")")),
      make_detail_row("Match Status", row$match_status)
    )
  }


  # ==================================================================
  # FTIR TAB
  # ==================================================================

  ftir_filtered <- reactive({
    dfs <- instrument_dfs()
    req(dfs$ftir)
    filter_instrument(dfs$ftir, input$ftir_quality_range, input$ftir_size_range,
                      input$ftir_material_filter, input$ftir_match_filter)
  })

  output$ftir_plot <- renderPlot({
    df <- ftir_filtered()
    req(nrow(df) > 0)
    bounds <- compute_bounds(instrument_dfs()$ftir, instrument_dfs()$raman)
    make_scatter(df, images$ftir, bounds,
                 paste0("FTIR Particles (", nrow(df), " shown)"),
                 match_colours = c(matched = "#2ca02c", unmatched = "#d62728"))
  })

  output$ftir_summary_text <- renderText({
    df <- ftir_filtered()
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched | ",
           length(unique(df$material)), " materials")
  })

  output$ftir_hover_info <- renderUI({
    hover <- input$ftir_hover
    if (is.null(hover)) return(tags$p(class = "text-muted", "Hover over a particle"))
    df <- ftir_filtered()
    row <- nearest_particle(df, hover$x, hover$y)
    single_detail_html(row, "FTIR", "AAU Quality")
  })


  # ==================================================================
  # RAMAN TAB
  # ==================================================================

  raman_filtered <- reactive({
    dfs <- instrument_dfs()
    req(dfs$raman)
    filter_instrument(dfs$raman, input$raman_hqi_range, input$raman_size_range,
                      input$raman_material_filter, input$raman_match_filter)
  })

  output$raman_plot <- renderPlot({
    df <- raman_filtered()
    req(nrow(df) > 0)
    bounds <- compute_bounds(instrument_dfs()$ftir, instrument_dfs()$raman)
    make_scatter(df, images$raman, bounds,
                 paste0("Raman Particles (", nrow(df), " shown)"),
                 match_colours = c(matched = "#1f77b4", unmatched = "#ff7f0e"))
  })

  output$raman_summary_text <- renderText({
    df <- raman_filtered()
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched | ",
           length(unique(df$material)), " materials")
  })

  output$raman_hover_info <- renderUI({
    hover <- input$raman_hover
    if (is.null(hover)) return(tags$p(class = "text-muted", "Hover over a particle"))
    df <- raman_filtered()
    row <- nearest_particle(df, hover$x, hover$y)
    single_detail_html(row, "Raman", "HQI")
  })


  # ==================================================================
  # OVERLAY TAB
  # ==================================================================

  overlay_matched <- reactive({
    d <- run_data()
    req(d$matched)
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

  output$overlay_plot <- renderPlot({
    dfs <- instrument_dfs()
    matched <- overlay_matched()
    layers <- input$overlay_layers
    bounds <- compute_bounds(dfs$ftir, dfs$raman)

    p <- ggplot() +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = "FTIR + Raman Overlay (aligned coordinates)",
           x = "X (\u00b5m)", y = "Y (\u00b5m)") +
      theme_minimal(base_size = 13) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "grey98", colour = NA),
        panel.grid       = element_line(colour = "grey90"),
        legend.position  = "bottom"
      )

    # Background image
    if (!is.null(images$overlay)) {
      p <- p + annotation_raster(images$overlay,
               xmin = bounds$x[1], xmax = bounds$x[2],
               ymin = bounds$y[1], ymax = bounds$y[2],
               interpolate = TRUE)
    }

    # Match lines
    if ("match_lines" %in% layers && nrow(matched) > 0) {
      seg_df <- data.frame(
        x    = matched$ftir_x_aligned, y    = matched$ftir_y_aligned,
        xend = matched$raman_x_um,     yend = matched$raman_y_um
      )
      p <- p + geom_segment(data = seg_df,
                              aes(x = x, y = y, xend = xend, yend = yend),
                              colour = "grey60", alpha = 0.3, linewidth = 0.3)
    }

    # Unmatched FTIR
    if ("unmatched_ftir" %in% layers) {
      um_f <- dfs$ftir[dfs$ftir$match_status == "unmatched", ]
      if (nrow(um_f) > 0) {
        p <- p + geom_point(data = um_f, aes(x = x, y = y, size = feret_max),
                             colour = "#d62728", alpha = 0.4, shape = 4) +
          scale_size_continuous(range = c(2, 12), guide = "none")
      }
    }

    # Unmatched Raman
    if ("unmatched_raman" %in% layers) {
      um_r <- dfs$raman[dfs$raman$match_status == "unmatched", ]
      if (nrow(um_r) > 0) {
        p <- p + geom_point(data = um_r, aes(x = x, y = y, size = feret_max),
                             colour = "#ff7f0e", alpha = 0.4, shape = 5)
      }
    }

    # Matched pairs: FTIR as triangle, Raman as circle
    if ("matched" %in% layers && nrow(matched) > 0) {
      ftir_pts <- data.frame(
        x = matched$ftir_x_aligned, y = matched$ftir_y_aligned,
        feret_max = matched$ftir_feret_max_um, instrument = "FTIR"
      )
      raman_pts <- data.frame(
        x = matched$raman_x_um, y = matched$raman_y_um,
        feret_max = matched$raman_feret_max_um, instrument = "Raman"
      )
      both <- rbind(ftir_pts, raman_pts)
      p <- p + geom_point(data = both,
                            aes(x = x, y = y, size = feret_max,
                                shape = instrument, colour = instrument),
                            alpha = 0.7) +
        scale_shape_manual(values = c(FTIR = 17, Raman = 16)) +
        scale_colour_manual(values = c(FTIR = "#2ca02c", Raman = "#1f77b4"))
    }

    p
  })

  output$overlay_summary_text <- renderText({
    m <- overlay_matched()
    dfs <- instrument_dfs()
    paste0(nrow(m), " matched pairs shown | ",
           sum(dfs$ftir$match_status == "unmatched"), " unmatched FTIR | ",
           sum(dfs$raman$match_status == "unmatched"), " unmatched Raman")
  })

  output$overlay_hover_info <- renderUI({
    hover <- input$overlay_hover
    if (is.null(hover))
      return(tags$p(class = "text-muted",
                    "Hover over a matched particle for the comparison table"))

    matched <- overlay_matched()
    if (nrow(matched) == 0) return(NULL)

    # Find nearest matched particle (check both FTIR and Raman positions)
    dx_f <- matched$ftir_x_aligned - hover$x
    dy_f <- matched$ftir_y_aligned - hover$y
    dist_f <- sqrt(dx_f^2 + dy_f^2)

    dx_r <- matched$raman_x_um - hover$x
    dy_r <- matched$raman_y_um - hover$y
    dist_r <- sqrt(dx_r^2 + dy_r^2)

    min_f <- min(dist_f)
    min_r <- min(dist_r)

    if (min(min_f, min_r) > 200) {
      return(tags$p(class = "text-muted",
                    "Hover over a matched particle for the comparison table"))
    }

    if (min_f <= min_r) {
      row <- matched[which.min(dist_f), ]
    } else {
      row <- matched[which.min(dist_r), ]
    }

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
