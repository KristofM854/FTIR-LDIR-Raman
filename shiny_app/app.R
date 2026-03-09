# =============================================================================
# app.R — Multi-Instrument Particle Viewer (Shiny + ggplot2)
# =============================================================================

source("global.R")

# Multiple-file selection: always enable (works on all platforms; especially
# useful on Windows where the OS Open dialog can select multiple files at once).
# Previously gated to .is_windows but Sys.info() returns the *server* OS, not
# the client browser OS — so the flag was always FALSE on Linux-hosted Shiny.
.is_windows <- TRUE

FEEDBACK_URL <- "https://forms.office.com/Pages/ResponsePage.aspx?id=kxTyotGkf0utB4Gcgk9cSupFgq19XipHgo77A48m1WhUNERLUjFMNzhNMDhHM0FDWFFWNTVUVVhOVC4u"

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
                                 quality_step, size_max = 1200,
                                 match_choices = c("matched", "unmatched")) {
  sidebarLayout(
    sidebarPanel(width = 3,
      h4(paste0(toupper(id_prefix), " Filters")),
      div(id = paste0(id_prefix, "_tour_quality"),
        sliderInput(paste0(id_prefix, "_quality_range"), quality_label,
                    min = quality_min, max = quality_max,
                    value = c(quality_min, quality_max), step = quality_step)),
      div(id = paste0(id_prefix, "_tour_size"),
        sliderInput(paste0(id_prefix, "_size_range"), "Feret Max (\u00b5m)",
                    min = 0, max = size_max, value = c(0, size_max), step = 5)),
      div(id = paste0(id_prefix, "_tour_material"),
        selectInput(paste0(id_prefix, "_material_filter"), "Material",
                    choices = c("All"), selected = "All", multiple = TRUE)),
      div(id = paste0(id_prefix, "_tour_match"),
        checkboxGroupInput(paste0(id_prefix, "_match_filter"), "Match Status",
                           choices = match_choices,
                           selected = unname(match_choices), inline = TRUE)),
      # Particle highlight: selectInput for single choice, plus text pattern box
      div(id = paste0(id_prefix, "_tour_highlight"),
        selectInput(paste0(id_prefix, "_highlight_particle"), "Highlight Particle",
                    choices = c("None"), selected = "None"),
        fluidRow(
          column(8, textInput(paste0(id_prefix, "_highlight_pattern"), NULL,
                              placeholder = "IDs: 1-10, MP_*, or MP_1,MP_5")),
          column(4, actionButton(paste0(id_prefix, "_highlight_apply"), "Apply",
                                 class = "btn-sm", style = "margin-top: 25px;"))
        )),
      hr(),
      div(class = "info-box",
          h5("Summary"), textOutput(paste0(id_prefix, "_summary_text"))),
      hr(),
      div(id = paste0(id_prefix, "_tour_image"),
        fileInput(paste0(id_prefix, "_image_upload"), "Background Image",
                  accept = c("image/png", "image/jpeg", "image/tiff",
                             ".tif", ".tiff", ".bmp", ".webp"),
                  multiple = .is_windows),
        fluidRow(
          column(6, numericInput(paste0(id_prefix, "_img_offset_x"),
                                 "Img X offset (\u00b5m)", value = 0, step = 25)),
          column(6, numericInput(paste0(id_prefix, "_img_offset_y"),
                                 "Img Y offset (\u00b5m)", value = 0, step = 25))
        ))
    ),
    mainPanel(width = 9,
      plotOutput(paste0(id_prefix, "_plot"), height = "650px",
                 click  = paste0(id_prefix, "_click"),
                 hover  = hoverOpts(paste0(id_prefix, "_hover"), delay = 100,
                                    delayType = "throttle"),
                 brush  = brushOpts(paste0(id_prefix, "_brush"),
                                    resetOnNew = TRUE),
                 dblclick = paste0(id_prefix, "_dblclick")),
      fluidRow(
        column(10, tags$p(class = "text-muted",
               "Drag to zoom in. Double-click to reset. Click a particle to add it to the selection.")),
        column(2, actionButton(paste0(id_prefix, "_reset_zoom"), "Reset Zoom",
                               class = "btn-sm btn-default",
                               style = "float:right; margin-top:2px;"))
      ),
      hr(),
      div(class = "info-box",
          h5("Particle Details (hover)"),
          detail_table_ui(paste0(id_prefix, "_hover_info"))),
      hr(),
      div(class = "info-box",
          fluidRow(
            column(8, h5("Selected Particles")),
            column(4, actionButton(paste0(id_prefix, "_clear_selection"), "Clear",
                                   class = "btn-sm btn-default",
                                   style = "float:right; margin-top:2px;"))
          ),
          detail_table_ui(paste0(id_prefix, "_selection_info"))),
      hr(),
      div(class = "info-box",
          h5("Material Summary"),
          uiOutput(paste0(id_prefix, "_plastics_summary")))
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
    title = tags$span(
      "Multi-Instrument Particle Viewer",
      actionButton("start_tutorial", "Start guided tour",
                   class = "btn-xs btn-default",
                   style = "margin-left:12px; vertical-align:middle;",
                   title = "Start guided tutorial"),
      tags$a(href = FEEDBACK_URL, target = "_blank",
             class = "btn btn-xs btn-info",
             style = "margin-left:8px; vertical-align:middle;",
             shiny::icon("comment"), "Give Feedback")
    ),
    id = "main_tabs",
    header = introjsUI(),

    # Tab 1: FTIR (PerkinElmer)
    tabPanel("FTIR (PerkinElmer)",
      div(id = "ftir_viewer",
        instrument_panel_ui("ftir", "AAU Quality", 0, 1, 0.01, 800,
          match_choices = c("Matched \u2194 Raman" = "matched",
                            "Unmatched (vs Raman)" = "unmatched")))
    ),

    # Tab 2: FTIR (Bruker) — shown only when data present
    tabPanel("FTIR (Bruker)",
      div(id = "ftir_bruker_viewer",
        instrument_panel_ui("ftir_bruker", "AAU Quality", 0, 1, 0.01, 800,
          match_choices = c("Matched \u2194 Raman" = "matched",
                            "Unmatched (vs Raman)" = "unmatched")))
    ),

    # Tab 3: Raman
    tabPanel("Raman",
      div(id = "raman_viewer",
        instrument_panel_ui("raman", "HQI", 0, 100, 1, 1200,
          match_choices = c("Matched \u2194 FTIR" = "matched",
                            "Unmatched (vs FTIR)" = "unmatched")))
    ),

    # Tab 3: LDIR
    tabPanel("LDIR",
      div(id = "ldir_viewer",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("LDIR Filters"),
          sliderInput("ldir_quality_range", "Quality",
                      min = 0, max = 1, value = c(0, 1), step = 0.01),
          sliderInput("ldir_size_range", "Feret Max (\u00b5m)",
                      min = 0, max = 1200, value = c(0, 1200), step = 5),
          selectInput("ldir_material_filter", "Material",
                      choices = c("All"), selected = "All", multiple = TRUE),
          checkboxGroupInput("ldir_match_filter", "Match Status",
                             choices = c("Matched \u2194 Raman" = "matched",
                                         "Unmatched (vs Raman)" = "unmatched"),
                             selected = c("matched", "unmatched"), inline = TRUE),
          sliderInput("ldir_score_range", "Match Score (LDIR\u2194Raman)",
                      min = 0, max = 20, value = c(0, 20), step = 0.1),
          sliderInput("ldir_coord_cost_range", "Coord Match Cost (image\u2194Excel)",
                      min = 0, max = 10, value = c(0, 10), step = 0.05),
          selectInput("ldir_highlight_particle", "Highlight Particle",
                      choices = c("None"), selected = "None"),
          hr(),
          h4("Image Overlay"),
          checkboxGroupInput("ldir_overlay_mode", "Display",
                             choices = c("Raw image" = "raw_image",
                                         "Processed image" = "processed_image",
                                         "Image-extracted particles" = "extracted_pts"),
                             selected = c("raw_image"),
                             inline = FALSE),
          hr(),
          div(class = "info-box",
              h5("Summary"), textOutput("ldir_summary_text")),
          hr(),
          fileInput("ldir_image_upload", "Background Image",
                    accept = c("image/png", "image/jpeg", "image/tiff",
                               ".tif", ".tiff", ".bmp", ".webp"),
                    multiple = .is_windows),
          fluidRow(
            column(6, numericInput("ldir_img_offset_x",
                                   "Img X offset (\u00b5m)", value = 0, step = 25)),
            column(6, numericInput("ldir_img_offset_y",
                                   "Img Y offset (\u00b5m)", value = 0, step = 25))
          )
        ),
        mainPanel(width = 9,
          plotOutput("ldir_plot", height = "650px",
                     click    = "ldir_click",
                     hover    = hoverOpts("ldir_hover", delay = 100,
                                          delayType = "throttle"),
                     brush    = brushOpts("ldir_brush",
                                          resetOnNew = TRUE),
                     dblclick = "ldir_dblclick"),
          fluidRow(
            column(10, tags$p(class = "text-muted",
                   "Drag to zoom in. Double-click to reset. Click a particle to add it to the selection.")),
            column(2, actionButton("ldir_reset_zoom", "Reset Zoom",
                                   class = "btn-sm btn-default",
                                   style = "float:right; margin-top:2px;"))
          ),
          hr(),
          div(class = "info-box",
              h5("Particle Details (hover)"),
              detail_table_ui("ldir_hover_info")),
          hr(),
          div(class = "info-box",
              fluidRow(
                column(8, h5("Selected Particles")),
                column(4, actionButton("ldir_clear_selection", "Clear",
                                       class = "btn-sm btn-default",
                                       style = "float:right; margin-top:2px;"))
              ),
              detail_table_ui("ldir_selection_info")),
          hr(),
          div(class = "info-box",
              h5("Material Summary"),
              uiOutput("ldir_plastics_summary"))
        )
      )
      ) # end div#ldir_viewer
    ),

    # Tab 4: Overlay (FTIR + Raman)
    tabPanel("Overlay",
      div(id = "overlay_panel",
      sidebarLayout(
        sidebarPanel(width = 3,
          # --- GLOBAL CONTROLS ---
          h4("Global Filters"),
          sliderInput("overlay_size_range", "Feret Max (\u00b5m)",
                      min = 0, max = 1200, value = c(0, 1200), step = 5),
          sliderInput("overlay_dist_range", "Match Distance (\u00b5m)",
                      min = 0, max = 100, value = c(0, 100), step = 1),
          hr(),

          # --- FT-IR SECTION ---
          h4("FT-IR", style = "color: #2ca02c; margin-bottom: 4px;"),
          sliderInput("overlay_ftir_quality", "AAU Quality",
                      min = 0, max = 1, value = c(0, 1), step = 0.01),
          sliderInput("overlay_ftir_size", "Feret Max (\u00b5m)",
                      min = 0, max = 800, value = c(0, 800), step = 5),
          selectizeInput("overlay_ftir_material", "Material",
                         choices = c("All"), selected = "All", multiple = TRUE),
          fluidRow(
            column(8, textInput("overlay_ftir_pattern", NULL,
                                placeholder = "Range (1-10) or pattern (MP_*)")),
            column(4, actionButton("overlay_ftir_apply_pattern", "Apply",
                                   class = "btn-sm", style = "margin-top: 25px;"))
          ),
          selectizeInput("overlay_ftir_particles", "Highlight Particles",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Select particles...",
                                        plugins = list("remove_button"))),
          hr(),

          # --- RAMAN SECTION ---
          h4("Raman", style = "color: #1f77b4; margin-bottom: 4px;"),
          sliderInput("overlay_raman_quality", "HQI",
                      min = 0, max = 100, value = c(0, 100), step = 1),
          sliderInput("overlay_raman_size", "Feret Max (\u00b5m)",
                      min = 0, max = 1200, value = c(0, 1200), step = 5),
          selectizeInput("overlay_raman_material", "Material",
                         choices = c("All"), selected = "All", multiple = TRUE),
          fluidRow(
            column(8, textInput("overlay_raman_pattern", NULL,
                                placeholder = "Range (1-10) or pattern")),
            column(4, actionButton("overlay_raman_apply_pattern", "Apply",
                                   class = "btn-sm", style = "margin-top: 25px;"))
          ),
          selectizeInput("overlay_raman_particles", "Highlight Particles",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Select particles...",
                                        plugins = list("remove_button"))),
          hr(),

          # --- LD-IR SECTION ---
          h4("LD-IR", style = "color: #d62728; margin-bottom: 4px;"),
          sliderInput("overlay_ldir_quality", "Quality",
                      min = 0, max = 1, value = c(0, 1), step = 0.01),
          sliderInput("overlay_ldir_size", "Feret Max (\u00b5m)",
                      min = 0, max = 1200, value = c(0, 1200), step = 5),
          selectizeInput("overlay_ldir_material", "Material",
                         choices = c("All"), selected = "All", multiple = TRUE),
          fluidRow(
            column(8, textInput("overlay_ldir_pattern", NULL,
                                placeholder = "Range (1-10) or pattern (A*)")),
            column(4, actionButton("overlay_ldir_apply_pattern", "Apply",
                                   class = "btn-sm", style = "margin-top: 25px;"))
          ),
          selectizeInput("overlay_ldir_particles", "Highlight Particles",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Select particles...",
                                        plugins = list("remove_button"))),
          hr(),

          # --- LAYER CHECKBOXES ---
          fluidRow(
            column(6, tags$label("Show Layers")),
            column(6, actionLink("overlay_toggle_all", "Select / Deselect All",
                                 style = "float:right; font-size:12px;"))
          ),
          checkboxGroupInput("overlay_layers", NULL,
                             choices = c("FTIR\u2194Raman matched" = "matched",
                                         "FTIR unmatched (vs Raman)" = "unmatched_ftir",
                                         "Raman unmatched (vs FTIR)" = "unmatched_raman",
                                         "FTIR\u2194Raman match lines" = "match_lines",
                                         "LDIR\u2194Raman matched" = "ldir_matched",
                                         "LDIR unmatched (vs Raman)" = "ldir_unmatched",
                                         "LDIR\u2194Raman match lines" = "ldir_lines",
                                         "Triple matches only" = "triple_only"),
                             selected = c("matched", "unmatched_ftir",
                                          "unmatched_raman", "ldir_matched"),
                             inline = FALSE),
          hr(),

          # --- SUMMARY + IMAGE ---
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
                     click = "overlay_click",
                     hover = hoverOpts("overlay_hover", delay = 200,
                                       delayType = "throttle"),
                     brush = brushOpts("overlay_brush",
                                       resetOnNew = TRUE),
                     dblclick = "overlay_dblclick"),
          fluidRow(
            column(10, tags$p(class = "text-muted",
                   "Drag to zoom. Double-click to reset. Click a particle to add it to the selection.")),
            column(2, actionButton("overlay_reset_zoom", "Reset Zoom",
                                   class = "btn-sm btn-default",
                                   style = "float:right; margin-top:2px;"))
          ),
          hr(),
          div(class = "info-box",
              fluidRow(
                column(8, h5("Match Details (hover / click / select)")),
                column(4, actionButton("overlay_clear_pin", "Clear",
                                       class = "btn-sm btn-default",
                                       style = "float:right; margin-top:2px;"))
              ),
              detail_table_ui("overlay_hover_info"))
        )
      )
      ) # end div#overlay_panel
    ),

    # Tab 5: Summary
    tabPanel("Summary",
      div(id = "summary_panel",
      fluidRow(
        column(10, offset = 1,
          div(class = "info-box", style = "margin-top: 20px;",
            h4("Material Comparison Across Instruments"),
            p(class = "text-muted",
              "Select a material family to compare counts across all instruments."),
            fluidRow(
              column(4, selectInput("summary_material_select", "Material Family",
                                    choices = c("PE"), selected = "PE")),
              column(8, plotOutput("summary_material_barplot", height = "350px"))
            )
          ),
          hr(),
          div(class = "info-box",
            fluidRow(
              column(8, h4("Plastics by Instrument")),
              column(4, checkboxInput("summary_use_filters",
                                      "Apply instrument filters",
                                      value = FALSE))
            ),
            p(class = "text-muted",
              "Material family counts per device. Toggle to apply each instrument's",
              "current quality / size / match-status filters."),
            uiOutput("summary_plastics_wide")
          ),
          hr(),
          div(class = "info-box",
            fluidRow(
              column(8, h4("Material Breakdown per Instrument (Pie Charts)")),
              column(4, radioButtons("pie_display_mode", NULL,
                                     choices = c("Absolute counts" = "abs",
                                                 "Relative (%)"    = "rel"),
                                     selected = "abs", inline = TRUE))
            ),
            p(class = "text-muted",
              "Plastic particles only (non-plastic materials excluded). ",
              "Toggle above to switch between particle counts and percentage shares."),
            fluidRow(
              column(6, plotOutput("pie_ftir",        height = "300px")),
              column(6, plotOutput("pie_raman",       height = "300px"))
            ),
            fluidRow(
              column(6, plotOutput("pie_ldir",        height = "300px")),
              column(6, plotOutput("pie_ftir_bruker", height = "300px"))
            )
          )
        )
      )
      ) # end div#summary_panel
    ),

    # Tab 6: Run Selector + Provenance
    tabPanel("Run Info",
      div(id = "runinfo_panel",
      fluidRow(
        column(8, offset = 2,
          div(class = "info-box", style = "margin-top: 20px;",
            h4("Run Selector"),
            p(class = "text-muted",
              "Select which pipeline run to view. The newest run is selected by default.",
              "Changing the run reloads all data in all tabs."),
            selectInput("run_selector", "Available Runs",
                        choices = character(0), selected = NULL, width = "100%"),
            actionButton("run_reload", "Reload Selected Run",
                         class = "btn-primary", icon = icon("refresh")),
            hr(),
            h4("Run Provenance"),
            uiOutput("run_provenance_ui")
          )
        )
      )
      ) # end div#runinfo_panel
    ),

    # Tab 7: Data Upload (fallback)
    tabPanel("Upload Data",
      div(id = "upload_panel",
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
      ) # end div#upload_panel
    )
  )
)


# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {

  # ------------------------------------------------------------------
  # Guided tutorial (rintrojs)
  # ------------------------------------------------------------------
  observeEvent(input$start_tutorial, {
    # Navigate to FTIR tab first so sidebar controls are in the DOM and visible
    updateNavbarPage(session, "main_tabs", selected = "FTIR (PerkinElmer)")
    introjs(session, options = list(
      nextLabel  = "Next",
      prevLabel  = "Back",
      doneLabel  = "Done",
      steps = list(
        # Step 1: Welcome (floating)
        list(intro = paste0(
          "Welcome to the <strong>Multi-Instrument Particle Viewer</strong>!<br><br>",
          "This app integrates particle data from FTIR (PerkinElmer &amp; Bruker), ",
          "Raman spectroscopy, and LDIR into a unified spatial explorer.<br><br>",
          "The tour starts on the <strong>FTIR tab</strong>. Use <em>Next</em> to walk through the controls.")),
        # Steps 2–9: FTIR tab link + sidebar controls (all visible on FTIR tab)
        list(element = "a[data-value='FTIR (PerkinElmer)']",
             intro   = "The <strong>FTIR (PerkinElmer)</strong> tab displays particles detected by the PerkinElmer FTIR instrument. Each dot is colour-coded by match status with Raman."),
        list(element = "#ftir_tour_quality",
             intro   = "The <strong>Quality slider</strong> filters particles by spectral match score (AAU for FTIR/LDIR, HQI for Raman). Drag either handle to set a minimum/maximum range."),
        list(element = "#ftir_tour_size",
             intro   = "The <strong>Feret Max slider</strong> filters by maximum particle diameter in \u00b5m \u2014 a proxy for particle size. Use it to isolate a specific size class."),
        list(element = "#ftir_tour_material",
             intro   = "The <strong>Material filter</strong> lets you show only particles of selected polymer types. Type or pick from the dropdown; multiple selections are supported."),
        list(element = "#ftir_tour_match",
             intro   = "The <strong>Match Status</strong> checkboxes toggle visibility of matched particles (spatially paired with Raman) and unmatched particles. Both are shown by default."),
        list(element = "#ftir_tour_highlight",
             intro   = "Use <strong>Highlight Particle</strong> to visually emphasise a particle by ID, or type a range (1\u201310) or wildcard pattern (MP_*) in the text box and click Apply. Selected particles are ringed in gold."),
        list(element = "#ftir_tour_image",
             intro   = "Optionally load a <strong>background image</strong> (filter membrane photo or microscope snapshot) to display behind the scatter plot. Use the X/Y offset fields below to align it with the particle coordinates."),
        list(element = "#ftir_plot",
             intro   = "<strong>Navigate the map</strong>: drag to zoom in on a region, double-click to reset the view. Click any particle dot to add it to the selection panel below the plot. Hover to see quick-look details."),
        # Steps 10–15: remaining tabs via always-visible navbar links
        list(element = "a[data-value='Raman']",
             intro   = "The <strong>Raman</strong> tab shows particles from Raman spectroscopy, colour-coded by their match status with FTIR. The same sidebar controls apply (quality shown as HQI)."),
        list(element = "a[data-value='LDIR']",
             intro   = "The <strong>LDIR</strong> tab displays particles from the Agilent 8700 Laser Direct Infrared instrument. Two extra sliders let you filter by LDIR\u2194Raman match score and image\u2194Excel coordinate match cost."),
        list(element = "a[data-value='Overlay']",
             intro   = "The <strong>Overlay</strong> tab superimposes FTIR, Raman, and LDIR particles on a single spatial map with connecting lines between matched pairs \u2014 useful for assessing spatial alignment quality."),
        list(element = "a[data-value='Summary']",
             intro   = "The <strong>Summary</strong> tab shows aggregate material counts and per-instrument pie charts. Toggle <em>Apply instrument filters</em> to reflect your current sidebar settings."),
        list(element = "a[data-value='Run Info']",
             intro   = "The <strong>Run Info</strong> tab shows metadata for the currently loaded pipeline run: timestamp, git commit hash, and input file provenance."),
        list(element = "a[data-value='Upload Data']",
             intro   = "The <strong>Upload Data</strong> tab lets you load a pipeline output CSV bundle directly from disk without needing a run directory on this server.")
      )
    ))
  })

  # ------------------------------------------------------------------
  # Load data (graceful when no pipeline output exists)
  # ------------------------------------------------------------------

  # Populate run selector on startup
  all_runs_available <- list_all_runs()
  if (length(all_runs_available) == 0) {
    message("[Particle Viewer] No pipeline output found in ../output/")
    message("[Particle Viewer] Working directory: ", getwd())
    message("[Particle Viewer] Use the 'Upload Data' tab to load CSV files manually.")
  } else {
    message("[Particle Viewer] Found ", length(all_runs_available), " run(s). Newest: ",
            names(all_runs_available)[1])
    updateSelectInput(session, "run_selector",
                      choices = all_runs_available,
                      selected = all_runs_available[1])
  }

  # Reactive value that can be updated by CSV uploads
  uploaded_data <- reactiveVal(NULL)

  # Reactive value: currently selected run directory (responds to selector + reload)
  selected_run_dir <- reactiveVal(
    if (length(all_runs_available) > 0) all_runs_available[[1]] else NULL
  )

  observeEvent(input$run_reload, {
    chosen <- input$run_selector
    if (!is.null(chosen) && nzchar(chosen) && dir.exists(chosen)) {
      message("[Particle Viewer] Switching to run: ", chosen)
      selected_run_dir(chosen)
      uploaded_data(NULL)   # clear any uploaded data when selecting a run
    }
  })

  # Auto-select when dropdown changes (without requiring the reload button)
  observeEvent(input$run_selector, {
    chosen <- input$run_selector
    if (!is.null(chosen) && nzchar(chosen) && dir.exists(chosen)) {
      selected_run_dir(chosen)
      uploaded_data(NULL)
    }
  }, ignoreInit = TRUE)

  run_data <- reactive({
    # User uploads take priority
    ud <- uploaded_data()
    if (!is.null(ud)) return(ud)

    # Use selected run directory
    run_dir <- selected_run_dir()
    if (is.null(run_dir) || !dir.exists(run_dir)) {
      # Fallback: try find_latest_run()
      ri <- find_latest_run()
      if (is.null(ri)) return(list())
      return(load_run_data(ri))
    }

    run_info_sel <- list(dir = run_dir, format = "subdir")
    load_run_data(run_info_sel)
  }) |> bindCache(selected_run_dir(), is.null(uploaded_data()))

  # Active manifest (changes with run selection)
  active_manifest <- reactive({
    ud <- uploaded_data()
    if (!is.null(ud)) return(list(is_missing = TRUE, run_id = "uploaded"))
    run_dir <- selected_run_dir()
    if (is.null(run_dir)) return(list(is_missing = TRUE))
    load_run_manifest(run_dir)
  })

  # ==================================================================
  # SUMMARY TAB
  # ==================================================================

  # Update material family dropdown when data changes
  observe({
    dfs <- instrument_dfs()
    all_mats <- character(0)
    for (nm in c("ftir", "ftir_bruker", "raman", "ldir")) {
      d <- dfs[[nm]]
      if (!is.null(d) && nrow(d) > 0 && "material" %in% names(d)) {
        fams <- classify_family_vec(d$material)
        all_mats <- c(all_mats, fams)
      }
    }
    all_mats <- sort(unique(all_mats[!is.na(all_mats) & all_mats != "Other"]))
    if (length(all_mats) == 0) all_mats <- "PE"
    sel <- if ("PE" %in% all_mats) "PE" else all_mats[1]
    updateSelectInput(session, "summary_material_select",
                      choices = all_mats, selected = sel)
  })

  # Pre-compute per-instrument material family counts — full data (cached)
  instrument_material_counts <- reactive({
    dfs <- instrument_dfs()
    device_keys <- c("FTIR (PerkinElmer)" = "ftir", "FTIR (Bruker)" = "ftir_bruker",
                     "Raman" = "raman", "LDIR" = "ldir")
    lapply(device_keys, function(key) {
      d <- dfs[[key]]
      if (is.null(d) || nrow(d) == 0 || !"material" %in% names(d)) return(NULL)
      table(classify_family_vec(d$material))
    })
  }) |> bindCache(selected_run_dir(), is.null(uploaded_data()))

  # Filtered counts — depends on each instrument's current filter state (not cached)
  instrument_material_counts_filtered <- reactive({
    filtered_list <- list(
      "FTIR (PerkinElmer)" = ftir_filtered(),
      "FTIR (Bruker)"      = ftir_bruker_filtered(),
      "Raman"              = raman_filtered(),
      "LDIR"               = ldir_filtered()
    )
    lapply(filtered_list, function(d) {
      if (is.null(d) || nrow(d) == 0 || !"material" %in% names(d)) return(NULL)
      table(classify_family_vec(d$material))
    })
  })

  # Resolve which counts to use based on toggle
  active_material_counts <- reactive({
    if (isTRUE(input$summary_use_filters))
      instrument_material_counts_filtered()
    else
      instrument_material_counts()
  })

  # Interactive barplot: count of selected material family across instruments
  output$summary_material_barplot <- renderPlot({
    sel_fam <- input$summary_material_select
    if (is.null(sel_fam) || !nzchar(sel_fam)) return(NULL)

    device_colors <- c("FTIR (PerkinElmer)" = "#2ca02c",
                       "FTIR (Bruker)" = "#9467bd",
                       "Raman" = "#1f77b4", "LDIR" = "#d62728")
    cts <- active_material_counts()
    counts <- vapply(names(cts), function(dev_label) {
      tbl <- cts[[dev_label]]
      if (is.null(tbl)) return(0L)
      as.integer(tbl[sel_fam] %||% 0L)
    }, integer(1))

    # Only show instruments that have data
    has_data <- !vapply(names(cts), function(dev_label) is.null(cts[[dev_label]]),
                        logical(1))
    counts <- counts[has_data]
    if (length(counts) == 0) {
      plot.new()
      text(0.5, 0.5, "No data available", cex = 1.2, col = "#6c757d")
      return(NULL)
    }

    use_filt <- isTRUE(input$summary_use_filters)
    bar_df <- data.frame(
      instrument = factor(names(counts), levels = names(counts)),
      count = as.integer(counts),
      stringsAsFactors = FALSE
    )
    ggplot2::ggplot(bar_df, ggplot2::aes(x = instrument, y = count, fill = instrument)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.3, size = 5.2) +
      ggplot2::scale_fill_manual(values = device_colors[names(counts)], guide = "none") +
      ggplot2::scale_y_continuous(limits = c(0, max(counts) * 1.50)) +
      ggplot2::labs(x = NULL, y = "Particle Count",
                    title = paste0(sel_fam, " across instruments",
                                   if (use_filt) " (filtered)" else "")) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                           margin = ggplot2::margin(b = 14)),
        plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10),
        axis.text.x = ggplot2::element_text(size = 14),
        panel.grid.major.x = ggplot2::element_blank()
      )
  })

  output$summary_plastics_wide <- renderUI({
    if (isTRUE(input$summary_use_filters)) {
      devices <- list(
        "FTIR (PerkinElmer)" = ftir_filtered(),
        "FTIR (Bruker)"      = ftir_bruker_filtered(),
        Raman                = raman_filtered(),
        LDIR                 = ldir_filtered()
      )
    } else {
      dfs <- instrument_dfs()
      devices <- list(
        "FTIR (PerkinElmer)" = dfs$ftir,
        "FTIR (Bruker)"      = dfs$ftir_bruker,
        Raman                = dfs$raman,
        LDIR                 = dfs$ldir
      )
    }
    # Remove devices with no data
    devices <- Filter(function(d) !is.null(d) && nrow(d) > 0, devices)
    if (length(devices) == 0)
      return(tags$p(class = "text-muted", "No data loaded."))

    per_dev  <- lapply(devices, summarise_plastics)
    all_fams <- unique(unlist(lapply(per_dev, `[[`, "family")))
    if (length(all_fams) == 0)
      return(tags$p(class = "text-muted", "No classified materials found."))

    # Order by category then alphabetically
    fam_cats <- classify_category_vec(all_fams)
    cat_order <- c("Synthetic", "Semi-synthetic", "Natural/Organic", "Unknown")
    fam_ord <- order(match(fam_cats, cat_order, nomatch = 99), all_fams)
    all_fams <- all_fams[fam_ord]
    fam_cats <- fam_cats[fam_ord]

    dev_names <- names(devices)
    header <- tags$tr(tags$th("Family"), tags$th("Category"),
                      lapply(dev_names, tags$th))
    cur_cat <- ""
    body_rows <- lapply(seq_along(all_fams), function(i) {
      fam <- all_fams[i]
      cat <- fam_cats[i]
      cells <- lapply(per_dev, function(dt) {
        idx <- match(fam, dt$family)
        tags$td(if (is.na(idx)) "0" else as.character(dt$n[idx]))
      })
      tags$tr(tags$td(tags$b(fam)), tags$td(cat), cells)
    })
    tags$table(class = "hover-tbl", header, body_rows)
  })

  # ------------------------------------------------------------------
  # Per-instrument pie charts (Summary tab)
  # ------------------------------------------------------------------

  # Colour palette for material families (consistent across charts)
  .pie_palette <- c(
    PE = "#e41a1c", PP = "#377eb8", PS = "#4daf4a", PET = "#984ea3",
    PVC = "#ff7f00", PA = "#a65628", PU = "#f781bf", PC = "#999999",
    PMMA = "#66c2a5", PTFE = "#fc8d62", PES = "#8da0cb", Other = "#e5c494"
  )

  # Build a pie chart for a single instrument data frame
  make_instrument_pie <- function(df, title, rel_mode) {
    if (is.null(df) || nrow(df) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::labs(title = title) +
               ggplot2::theme_void(base_size = 14) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")))
    }
    fam  <- classify_family_vec(df$material)
    # Keep only plastic families (exclude non-plastic / unknown categories)
    cat  <- classify_category_vec(fam)
    keep <- cat %in% c("Synthetic", "Semi-synthetic")
    fam  <- fam[keep]
    if (length(fam) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::labs(title = title, subtitle = "No plastic particles") +
               ggplot2::theme_void(base_size = 14) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")))
    }
    tbl <- sort(table(fam), decreasing = TRUE)
    pie_df <- data.frame(material = names(tbl), count = as.integer(tbl),
                          stringsAsFactors = FALSE)
    total  <- sum(pie_df$count)
    pie_df$pct   <- pie_df$count / total * 100
    pie_df$label <- if (rel_mode)
      paste0(round(pie_df$pct, 1), "%")
    else
      as.character(pie_df$count)
    # Assign colours; grey for unmapped families
    fam_colors <- .pie_palette[pie_df$material]
    fam_colors[is.na(fam_colors)] <- "#cccccc"
    names(fam_colors) <- pie_df$material
    pie_df$material <- factor(pie_df$material, levels = pie_df$material)

    ggplot2::ggplot(pie_df, ggplot2::aes(x = "", y = count, fill = material)) +
      ggplot2::geom_col(width = 1, colour = "white", linewidth = 0.4) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 4, colour = "white", fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = fam_colors, name = "Material") +
      ggplot2::labs(title = title,
                    subtitle = paste0("n = ", total, " plastic particles")) +
      ggplot2::theme_void(base_size = 14) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "#555555"),
        legend.position = "right",
        legend.text     = ggplot2::element_text(size = 10),
        legend.title    = ggplot2::element_text(size = 11, face = "bold")
      )
  }

  # Helper reactive: resolve per-instrument data (filtered or unfiltered)
  pie_data <- reactive({
    if (isTRUE(input$summary_use_filters)) {
      list(ftir        = ftir_filtered(),
           raman       = raman_filtered(),
           ldir        = ldir_filtered(),
           ftir_bruker = ftir_bruker_filtered())
    } else {
      dfs <- instrument_dfs()
      list(ftir        = dfs$ftir,
           raman       = dfs$raman,
           ldir        = dfs$ldir,
           ftir_bruker = dfs$ftir_bruker)
    }
  })

  output$pie_ftir <- renderPlot({
    rel <- identical(input$pie_display_mode, "rel")
    make_instrument_pie(pie_data()$ftir, "FTIR (PerkinElmer)", rel)
  })
  output$pie_raman <- renderPlot({
    rel <- identical(input$pie_display_mode, "rel")
    make_instrument_pie(pie_data()$raman, "Raman", rel)
  })
  output$pie_ldir <- renderPlot({
    rel <- identical(input$pie_display_mode, "rel")
    make_instrument_pie(pie_data()$ldir, "LDIR", rel)
  })
  output$pie_ftir_bruker <- renderPlot({
    rel <- identical(input$pie_display_mode, "rel")
    make_instrument_pie(pie_data()$ftir_bruker, "FTIR (Bruker)", rel)
  })

  # Provenance panel UI
  output$run_provenance_ui <- renderUI({
    m <- active_manifest()

    if (isTRUE(m$is_missing)) {
      warn_box <- div(class = "alert alert-warning",
        tags$b("No manifest.json found for this run."),
        tags$p("Runs created before manifest support was added will not have provenance data.",
               "Re-run the pipeline to generate a manifest.")
      )
      return(warn_box)
    }

    run_id_val  <- if (!is.null(m$run_id)) m$run_id else "unknown"
    ts_val      <- if (!is.null(m$timestamp) && !is.na(m$timestamp)) m$timestamp else "N/A"
    git_val     <- if (!is.null(m$git_commit) && !is.na(m$git_commit) && nzchar(m$git_commit))
                     m$git_commit else "N/A"
    stage_val   <- if (!is.null(m$stage)) m$stage else "unknown"

    # Input files
    input_rows <- tryCatch({
      if (is.null(m$inputs) || length(m$inputs) == 0) return(list())
      lapply(m$inputs, function(inp) {
        nm   <- if (!is.null(inp$name)) inp$name else "?"
        base <- if (!is.null(inp$basename)) inp$basename else
                  if (!is.null(inp$path)) basename(inp$path) else "N/A"
        md5  <- if (!is.null(inp$md5) && !is.na(inp$md5))
                  substr(inp$md5, 1, 12) else "N/A"
        tags$tr(tags$td(tags$b(nm)), tags$td(base), tags$td(code(md5)))
      })
    }, error = function(e) list())

    # Instrument image info — reads from manifest$images (new) or manifest$ldir_image (old)
    images_src <- if (!is.null(m$images)) m$images else {
      if (!is.null(m$ldir_image)) list(ldir = m$ldir_image) else list()
    }
    images_ui_panels <- tryCatch({
      if (length(images_src) == 0) return(list())
      lapply(names(images_src), function(instr) {
        li <- images_src[[instr]]
        if (is.null(li)) return(NULL)
        fmt       <- if (!is.null(li$detected_format)) li$detected_format else "?"
        orig_dim  <- if (!is.null(li$orig_width) && !is.na(li$orig_width))
                       paste0(li$orig_width, " x ", li$orig_height) else "?"
        canon_dim <- if (!is.null(li$canonical_width) && !is.na(li$canonical_width))
                       paste0(li$canonical_width, " x ", li$canonical_height) else "?"
        prev_sc   <- if (!is.null(li$preview_scale) && !is.na(li$preview_scale))
                       paste0(round(li$preview_scale * 100), "%") else "?"
        md5_short <- if (!is.null(li$md5) && !is.na(li$md5))
                       substr(li$md5, 1, 12) else "N/A"
        orig_bn   <- if (!is.null(li$orig_basename) && nzchar(li$orig_basename))
                       li$orig_basename else "?"
        ext_warn  <- tryCatch({
          ext <- toupper(tools::file_ext(orig_bn))
          if (nzchar(ext) && toupper(fmt) != ext && fmt != "unknown")
            tags$span(class = "label label-warning",
                      paste0("Extension mismatch: .", tolower(ext),
                             " but signature=", fmt))
          else NULL
        }, error = function(e) NULL)

        div(
          h5(paste0(toupper(instr), " Image")),
          ext_warn,
          tags$table(class = "hover-tbl",
            tags$tr(tags$th("Field"), tags$th("Value")),
            tags$tr(tags$td("Original file"),      tags$td(code(orig_bn))),
            tags$tr(tags$td("Detected format"),    tags$td(tags$b(fmt))),
            tags$tr(tags$td("Original dims"),      tags$td(orig_dim)),
            tags$tr(tags$td("Canonical PNG dims"), tags$td(canon_dim)),
            tags$tr(tags$td("Preview scale"),      tags$td(prev_sc)),
            tags$tr(tags$td("MD5 (12 chars)"),     tags$td(code(md5_short)))
          )
        )
      })
    }, error = function(e) list())

    tagList(
      div(class = if (stage_val == "export_complete") "alert alert-success"
                  else "alert alert-info",
        tags$b(paste0("Run: ", run_id_val)),
        tags$span(paste0(" | Stage: ", stage_val))
      ),
      tags$table(class = "hover-tbl",
        tags$tr(tags$th("Field"), tags$th("Value")),
        tags$tr(tags$td(tags$b("Run ID")),    tags$td(run_id_val)),
        tags$tr(tags$td(tags$b("Timestamp")), tags$td(ts_val)),
        tags$tr(tags$td(tags$b("Git commit")),tags$td(code(git_val))),
        tags$tr(tags$td(tags$b("Stage")),     tags$td(stage_val))
      ),
      br(),
      if (length(input_rows) > 0) {
        div(
          h5("Input Files"),
          tags$table(class = "hover-tbl",
            tags$tr(tags$th("Input"), tags$th("File"), tags$th("MD5 (12 chars)")),
            input_rows
          )
        )
      } else NULL,
      br(),
      images_ui_panels
    )
  })

  has_data <- reactive({
    d <- run_data()
    !is.null(d$matched)            ||
    !is.null(d$unmatched_ftir)     ||
    !is.null(d$unmatched_raman)    ||
    !is.null(d$ldir_raman_matched) ||
    !is.null(d$unmatched_ldir)     ||
    !is.null(d$unmatched_ftir_bruker)
  })

  instrument_dfs <- reactive({
    if (!has_data()) return(list(ftir = NULL, raman = NULL, ldir = NULL,
                                 ftir_bruker = NULL))
    build_instrument_dfs(run_data())
  }) |> bindCache(selected_run_dir(), is.null(uploaded_data()))

  # Per-instrument full-data reactives (avoid repeated instrument_dfs()$X calls)
  ftir_df_full        <- reactive({ instrument_dfs()$ftir })
  raman_df_full       <- reactive({ instrument_dfs()$raman })
  ldir_df_full        <- reactive({ instrument_dfs()$ldir })
  ftir_bruker_df_full <- reactive({ instrument_dfs()$ftir_bruker })

  # Which devices actually have data?
  has_device <- reactive({
    list(
      ftir        = !is.null(ftir_df_full())        && nrow(ftir_df_full())        > 0,
      raman       = !is.null(raman_df_full())       && nrow(raman_df_full())       > 0,
      ldir        = !is.null(ldir_df_full())        && nrow(ldir_df_full())        > 0,
      ftir_bruker = !is.null(ftir_bruker_df_full()) && nrow(ftir_bruker_df_full()) > 0
    )
  })

  # Hide/show FTIR Bruker tab and Overlay tab based on available data
  observe({
    hd <- has_device()
    if (hd$ftir_bruker) showTab("main_tabs", "FTIR (Bruker)")
    else                hideTab("main_tabs", "FTIR (Bruker)")
    if (sum(unlist(hd)) >= 2) showTab("main_tabs", "Overlay")
    else                      hideTab("main_tabs", "Overlay")
  })

  # Click-to-select state: accumulated particle IDs per instrument
  selected_ids <- reactiveValues(
    ftir        = character(0),
    raman       = character(0),
    ldir        = character(0),
    ftir_bruker = character(0)
  )

  # Debounced slider inputs (300ms) — prevents re-render on every pixel drag
  # Individual tabs
  ftir_quality_range_d        <- debounce(reactive(input$ftir_quality_range), 300)
  ftir_size_range_d           <- debounce(reactive(input$ftir_size_range), 300)
  raman_quality_range_d       <- debounce(reactive(input$raman_quality_range), 300)
  raman_size_range_d          <- debounce(reactive(input$raman_size_range), 300)
  ldir_quality_range_d        <- debounce(reactive(input$ldir_quality_range), 300)
  ldir_size_range_d           <- debounce(reactive(input$ldir_size_range), 300)
  ldir_score_range_d          <- debounce(reactive(input$ldir_score_range), 300)
  ldir_coord_cost_range_d     <- debounce(reactive(input$ldir_coord_cost_range), 300)
  ftir_bruker_quality_range_d <- debounce(reactive(input$ftir_bruker_quality_range), 300)
  ftir_bruker_size_range_d    <- debounce(reactive(input$ftir_bruker_size_range), 300)
  # Overlay global
  overlay_size_range_d   <- debounce(reactive(input$overlay_size_range), 300)
  overlay_dist_range_d   <- debounce(reactive(input$overlay_dist_range), 300)
  # Overlay per-instrument
  overlay_ftir_quality_d <- debounce(reactive(input$overlay_ftir_quality), 300)
  overlay_ftir_size_d    <- debounce(reactive(input$overlay_ftir_size), 300)
  overlay_raman_quality_d <- debounce(reactive(input$overlay_raman_quality), 300)
  overlay_raman_size_d   <- debounce(reactive(input$overlay_raman_size), 300)
  overlay_ldir_quality_d <- debounce(reactive(input$overlay_ldir_quality), 300)
  overlay_ldir_size_d    <- debounce(reactive(input$overlay_ldir_size), 300)

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
    ftir_d <- ftir_df_full()
    if (!is.null(raw_ftir_img)) {
      px <- if (!is.null(ftir_d) && nrow(ftir_d) > 0) ftir_d$x_orig else NULL
      py <- if (!is.null(ftir_d) && nrow(ftir_d) > 0) ftir_d$y_orig else NULL
      return(estimate_ftir_scan_bounds(raw_ftir_img, px, py))
    }

    # Fallback: round up particle coords to nearest 500 µm
    if (!is.null(ftir_d) && nrow(ftir_d) > 0) {
      return(list(xmin = 0,
                  xmax = ceiling(max(ftir_d$x_orig, na.rm = TRUE) / 500) * 500,
                  ymin = 0,
                  ymax = ceiling(max(ftir_d$y_orig, na.rm = TRUE) / 500) * 500))
    }
    NULL
  })

  # ------------------------------------------------------------------
  # Helper: write overlay contract JSON for viewer debugging
  # ------------------------------------------------------------------
  # Writes a small JSON file that records the image bounds and particle bounds
  # for a given instrument viewer, allowing off-line comparison of the two.
  # The "delta" fields show how much the particle extent deviates from the image
  # placement — ideally all deltas are < 500 µm.
  write_overlay_contract <- function(instrument, image_info, points_df, run_dir) {
    if (is.null(image_info) || is.null(run_dir)) return(invisible(NULL))
    debug_dir <- file.path(run_dir, "debug")
    if (!dir.exists(debug_dir))
      tryCatch(dir.create(debug_dir, recursive = TRUE), error = function(e) NULL)
    x <- if (!is.null(points_df)) points_df$x_orig[is.finite(points_df$x_orig)] else numeric(0)
    y <- if (!is.null(points_df)) points_df$y_orig[is.finite(points_df$y_orig)] else numeric(0)
    info <- list(
      instrument       = instrument,
      image_bounds_um  = list(xmin = image_info$xmin, xmax = image_info$xmax,
                              ymin = image_info$ymin, ymax = image_info$ymax),
      points_bounds_um = list(xmin = if (length(x) > 0) min(x) else NA_real_,
                              xmax = if (length(x) > 0) max(x) else NA_real_,
                              ymin = if (length(y) > 0) min(y) else NA_real_,
                              ymax = if (length(y) > 0) max(y) else NA_real_),
      n_points         = length(x)
    )
    info$delta_bounds_um <- list(
      xmin = info$points_bounds_um$xmin - info$image_bounds_um$xmin,
      xmax = info$points_bounds_um$xmax - info$image_bounds_um$xmax,
      ymin = info$points_bounds_um$ymin - info$image_bounds_um$ymin,
      ymax = info$points_bounds_um$ymax - info$image_bounds_um$ymax
    )
    path <- file.path(debug_dir, paste0("overlay_contract_", instrument, ".json"))
    tryCatch(
      writeLines(jsonlite::toJSON(info, pretty = TRUE, auto_unbox = TRUE, null = "null"),
                 path),
      error = function(e) NULL
    )
    invisible(path)
  }

  # ------------------------------------------------------------------
  # Raw image rasters
  # ------------------------------------------------------------------
  ftir_raw_image    <- reactiveVal(NULL)   # FTIR "Average Abs" image
  raman_image       <- reactiveVal(NULL)   # Raman microscope image (Raman tab + Overlay tab)
  raman_image_path  <- reactiveVal(NULL)   # File path of the Raman image (for TIFF metadata extraction)
  ldir_raw_image    <- reactiveVal(NULL)   # LDIR particle map image

  # FTIR tab: raw image placed at native FTIR scan bounds — no transform needed.
  # FTIR image placed at the actual particle extent (x_orig / y_orig).
  # Using particle positions for bounds is more reliable than heuristic scan-area
  # estimation from pixel counts, which produced wrong bounds → tiled appearance.
  ftir_native_image_info <- reactive({
    raw <- ftir_raw_image()
    if (is.null(raw)) return(NULL)
    ftir_d <- ftir_df_full()
    if (is.null(ftir_d) || nrow(ftir_d) == 0) return(NULL)
    x_vals <- ftir_d$x_orig[is.finite(ftir_d$x_orig)]
    y_vals <- ftir_d$y_orig[is.finite(ftir_d$y_orig)]
    if (length(x_vals) == 0) return(NULL)
    ox <- if (!is.null(input$ftir_img_offset_x)) input$ftir_img_offset_x else 0
    oy <- if (!is.null(input$ftir_img_offset_y)) input$ftir_img_offset_y else 0
    list(raster = raw,
         xmin = min(x_vals) + ox, xmax = max(x_vals) + ox,
         ymin = min(y_vals) + oy, ymax = max(y_vals) + oy)
  })

  # Raman tab: image placement with 3-tier priority cascade:
  #   1. Physical extent (WITec center + width/height in µm) — resize-invariant
  #   2. Known scale (raman_um_per_px or auto-detected from TIFF) — centroid-centred
  #   3. Fallback: aspect-ratio-preserving bounds via compute_image_bounds()
  raman_native_image_info <- reactive({
    raw <- raman_image()
    if (is.null(raw)) return(NULL)
    raman_df <- raman_df_full()
    ox <- if (!is.null(input$raman_img_offset_x)) input$raman_img_offset_x else 0
    oy <- if (!is.null(input$raman_img_offset_y)) input$raman_img_offset_y else 0

    x_vals <- if (!is.null(raman_df)) raman_df$x_orig[is.finite(raman_df$x_orig)] else numeric(0)
    y_vals <- if (!is.null(raman_df)) raman_df$y_orig[is.finite(raman_df$y_orig)] else numeric(0)

    cfg <- tryCatch(active_manifest()$config_snapshot, error = function(e) list())
    h_px <- nrow(raw); w_px <- ncol(raw)

    # --- Priority 1: Physical extent from WITec metadata (resize-invariant) ---
    w_um <- cfg$raman_image_width_um
    h_um <- cfg$raman_image_height_um
    cx   <- cfg$raman_image_center_x_um
    cy   <- cfg$raman_image_center_y_um

    if (!is.null(w_um) && !is.null(h_um) && !is.null(cx) && !is.null(cy)) {
      return(list(raster = raw,
                  xmin = cx - w_um / 2 + ox,
                  xmax = cx + w_um / 2 + ox,
                  ymin = cy - h_um / 2 + oy,
                  ymax = cy + h_um / 2 + oy))
    }

    # --- Priority 2: Known scale (config or auto-detected from TIFF DPI) ---
    um_per_px <- cfg$raman_um_per_px
    if (is.null(um_per_px) || !is.numeric(um_per_px) || um_per_px <= 0) {
      um_per_px <- extract_tiff_um_per_px(raman_image_path())
    }

    if (!is.null(um_per_px) && length(x_vals) > 0) {
      cx_um <- mean(x_vals)
      cy_um <- mean(y_vals)
      half_w <- w_px * um_per_px / 2
      half_h <- h_px * um_per_px / 2
      return(list(raster = raw,
                  xmin = cx_um - half_w + ox, xmax = cx_um + half_w + ox,
                  ymin = cy_um - half_h + oy, ymax = cy_um + half_h + oy))
    }

    # --- Priority 3: Fallback — aspect-ratio-preserving bounds from particles ---
    if (length(x_vals) == 0) return(NULL)
    b <- compute_image_bounds(raw, x_vals, y_vals, padding_um = 300)
    list(raster = raw,
         xmin = b$xmin + ox, xmax = b$xmax + ox,
         ymin = b$ymin + oy, ymax = b$ymax + oy)
  })

  # Overlay tab: Raman microscope image placed at Raman particle extent in
  # normalized (centered) coordinates.  Uses same aspect-ratio-preserving
  # logic as the Raman native tab so the image appears identical in both views.
  overlay_image_info <- reactive({
    raw <- raman_image()
    if (is.null(raw)) return(NULL)
    raman_d <- raman_df_full()
    if (is.null(raman_d) || nrow(raman_d) == 0) return(NULL)

    # Use Raman particles only — the image is a Raman microscope photo
    raman_x <- raman_d$x[is.finite(raman_d$x)]
    raman_y <- raman_d$y[is.finite(raman_d$y)]
    if (length(raman_x) == 0) return(NULL)

    # Apply user fine-tuning offsets
    ox <- if (!is.null(input$overlay_img_offset_x)) input$overlay_img_offset_x else 0
    oy <- if (!is.null(input$overlay_img_offset_y)) input$overlay_img_offset_y else 0

    b <- compute_image_bounds(raw, raman_x, raman_y, padding_um = 300)
    list(raster = raw,
         xmin = b$xmin + ox, xmax = b$xmax + ox,
         ymin = b$ymin + oy, ymax = b$ymax + oy)
  })

  # LDIR tab: image placed at the full scan-circle extent.
  # LDIR coordinates are circle-calibrated and centred at (0,0) via
  # map_pixels_to_um_circle(), so the scan area spans ±(scan_diameter/2) µm.
  # annotation_raster places row 1 at ymax.  Since row 1 of the PNG is the
  # top of the physical scan = high positive y in circle coords, no row-flip
  # is needed.
  # IMPORTANT: use the scan diameter (6500 µm half-extent) as the floor —
  # never just the particle data range, because sparse joins leave most
  # particles without coordinates, causing asymmetric bounds and a
  # stretched/clipped image.
  ldir_native_image_info <- reactive({
    raw <- ldir_raw_image()
    if (is.null(raw)) return(NULL)
    ox <- if (!is.null(input$ldir_img_offset_x)) input$ldir_img_offset_x else 0
    oy <- if (!is.null(input$ldir_img_offset_y)) input$ldir_img_offset_y else 0

    m <- tryCatch(active_manifest(), error = function(e) NULL)

    # Circle-based bounds: derive exact image placement from scan-circle calibration.
    # cx_px/cy_px are the circle centre in pixel space; scale_um_per_px converts
    # pixels to µm.  The image spans from -cx_px*scale to (w-cx_px)*scale in x
    # and from (cy_px-h)*scale to cy_px*scale in y (y upward, row 0 = ymax).
    ci <- if (!is.null(m)) m$ldir_circle else NULL
    if (!is.null(ci) && !is.null(ci$scale_um_per_px) && ci$scale_um_per_px > 0) {
      s  <- ci$scale_um_per_px
      cx <- ci$cx_px;  cy <- ci$cy_px
      w  <- ci$image_width_px;  h <- ci$image_height_px
      return(list(raster = raw,
                  xmin = -cx * s + ox,      xmax = (w - cx) * s + ox,
                  ymin = (cy - h) * s + oy, ymax = cy * s + oy))
    }

    # Fallback: symmetric ±half_um from scan diameter (old behaviour)
    scan_diam_um <- tryCatch({
      d <- m$config_snapshot$ldir_scan_diameter_um
      if (!is.null(d) && is.numeric(d) && d > 0) as.integer(d) else 13000L
    }, error = function(e) 13000L)
    half_um <- as.integer(scan_diam_um / 2L)
    # Widen if any particle coordinates actually exceed the expected half-extent
    ldir_df <- ldir_df_full()
    if (!is.null(ldir_df) && nrow(ldir_df) > 0) {
      xvals <- ldir_df$x_orig[is.finite(ldir_df$x_orig)]
      yvals <- ldir_df$y_orig[is.finite(ldir_df$y_orig)]
      if (length(xvals) > 0 && length(yvals) > 0) {
        max_abs <- max(abs(c(xvals, yvals)))
        if (max_abs > half_um)
          half_um <- ceiling(max_abs / 500) * 500
      }
    }
    list(raster = raw,
         xmin = -half_um + ox, xmax = half_um + ox,
         ymin = -half_um + oy, ymax = half_um + oy)
  })

  # Write overlay contract JSONs whenever image info or particles change.
  # These lightweight JSON files record image placement vs. particle extent,
  # making it easy to verify alignment without opening the Shiny app.
  observe({
    ii  <- raman_native_image_info()
    df  <- raman_df_full()
    run <- selected_run_dir()
    if (!is.null(ii) && !is.null(run))
      write_overlay_contract("raman", ii, df, run)
  })
  observe({
    ii  <- ldir_native_image_info()
    df  <- ldir_df_full()
    run <- selected_run_dir()
    if (!is.null(ii) && !is.null(run))
      write_overlay_contract("ldir", ii, df, run)
  })

  # Load instrument images from the run manifest (manifest-driven, no hardcoded paths)
  observeEvent(selected_run_dir(), {
    run_dir <- selected_run_dir()
    if (is.null(run_dir) || !dir.exists(run_dir)) {
      ftir_raw_image(NULL)
      raman_image(NULL)
      raman_image_path(NULL)
      ldir_raw_image(NULL)
      return()
    }
    m   <- load_run_manifest(run_dir)
    img <- get_run_image_paths(m, run_dir)

    if (!is.null(img$ftir)) {
      raw <- load_image_raster(img$ftir)
      if (!is.null(raw)) ftir_raw_image(raw)
    } else {
      ftir_raw_image(NULL)
    }

    if (!is.null(img$raman)) {
      raw <- load_image_raster(img$raman)
      if (!is.null(raw)) {
        raman_image(raw)
        raman_image_path(img$raman)
      }
    } else {
      raman_image(NULL)
      raman_image_path(NULL)
    }

    if (!is.null(img$ldir)) {
      raw <- load_image_raster(img$ldir)
      if (!is.null(raw)) ldir_raw_image(raw)
    } else {
      ldir_raw_image(NULL)
    }
  }, ignoreNULL = FALSE)

  # Handle uploaded images
  observeEvent(input$ftir_image_upload, {
    raw <- load_image_raster(input$ftir_image_upload$datapath)
    message("FTIR upload raster dim: ", paste(dim(raw), collapse=" x "))
    if (!is.null(raw)) ftir_raw_image(raw)
  })

  observeEvent(input$raman_image_upload, {
    raw <- load_image_raster(input$raman_image_upload$datapath)
    if (!is.null(raw)) {
      raman_image(raw)
      raman_image_path(input$raman_image_upload$datapath)
    }
  })

  observeEvent(input$overlay_image_upload, {
    raw <- load_image_raster(input$overlay_image_upload$datapath)
    if (!is.null(raw)) {
      raman_image(raw)
      raman_image_path(input$overlay_image_upload$datapath)
    }
  })

  observeEvent(input$ldir_image_upload, {
    raw <- load_image_raster(input$ldir_image_upload$datapath)
    if (!is.null(raw)) ldir_raw_image(raw)
  })

  # ------------------------------------------------------------------
  # Update filter controls from data
  # ------------------------------------------------------------------
  observe({
    ftir_d <- ftir_df_full()
    raman_d <- raman_df_full()
    ldir_d <- ldir_df_full()

    # --- FTIR controls (individual tab + overlay) ---
    if (!is.null(ftir_d) && nrow(ftir_d) > 0) {
      ftir <- ftir_d
      ftir_mats <- sort(unique(ftir$material))
      ftir_ids  <- natural_sort_ids(unique(ftir$particle_id))
      q_range   <- range(ftir$quality, na.rm = TRUE)
      s_max     <- ceiling(max(ftir$feret_max, na.rm = TRUE) / 10) * 10

      # Individual tab
      updateSelectInput(session, "ftir_material_filter",
                        choices = c("All", ftir_mats), selected = "All")
      updateSelectInput(session, "ftir_highlight_particle",
                        choices = c("None", ftir_ids))
      updateSliderInput(session, "ftir_quality_range",
                        min = floor(q_range[1] * 100) / 100,
                        max = ceiling(q_range[2] * 100) / 100,
                        value = c(floor(q_range[1] * 100) / 100,
                                  ceiling(q_range[2] * 100) / 100))
      updateSliderInput(session, "ftir_size_range", min = 0, max = s_max,
                        value = c(0, s_max))

      # Overlay per-instrument
      updateSelectizeInput(session, "overlay_ftir_material",
                           choices = c("All", ftir_mats), selected = "All")
      updateSelectizeInput(session, "overlay_ftir_particles",
                           choices = ftir_ids, selected = character(0))
      updateSliderInput(session, "overlay_ftir_quality",
                        min = floor(q_range[1] * 100) / 100,
                        max = ceiling(q_range[2] * 100) / 100,
                        value = c(floor(q_range[1] * 100) / 100,
                                  ceiling(q_range[2] * 100) / 100))
      updateSliderInput(session, "overlay_ftir_size", min = 0, max = s_max,
                        value = c(0, s_max))
    }

    # --- Raman controls (individual tab + overlay) ---
    if (!is.null(raman_d) && nrow(raman_d) > 0) {
      raman <- raman_d
      raman_mats <- sort(unique(raman$material))
      raman_ids  <- natural_sort_ids(unique(raman$particle_id))
      q_range    <- range(raman$quality, na.rm = TRUE)
      s_max      <- ceiling(max(raman$feret_max, na.rm = TRUE) / 10) * 10

      # Individual tab
      updateSelectInput(session, "raman_material_filter",
                        choices = c("All", raman_mats), selected = "All")
      updateSelectInput(session, "raman_highlight_particle",
                        choices = c("None", raman_ids))
      updateSliderInput(session, "raman_quality_range",
                        min = floor(q_range[1]), max = ceiling(q_range[2]),
                        value = c(floor(q_range[1]), ceiling(q_range[2])))
      updateSliderInput(session, "raman_size_range", min = 0, max = s_max,
                        value = c(0, s_max))

      # Overlay per-instrument
      updateSelectizeInput(session, "overlay_raman_material",
                           choices = c("All", raman_mats), selected = "All")
      updateSelectizeInput(session, "overlay_raman_particles",
                           choices = raman_ids, selected = character(0))
      updateSliderInput(session, "overlay_raman_quality",
                        min = floor(q_range[1]), max = ceiling(q_range[2]),
                        value = c(floor(q_range[1]), ceiling(q_range[2])))
      updateSliderInput(session, "overlay_raman_size", min = 0, max = s_max,
                        value = c(0, s_max))
    }

    # --- LDIR controls (individual tab + overlay) ---
    if (!is.null(ldir_d) && nrow(ldir_d) > 0) {
      ldir <- ldir_d
      ldir_mats <- sort(unique(ldir$material))
      ldir_ids  <- natural_sort_ids(unique(ldir$particle_id))
      q_range   <- range(ldir$quality, na.rm = TRUE)
      s_max     <- ceiling(max(ldir$feret_max, na.rm = TRUE) / 10) * 10

      # Individual tab
      updateSelectInput(session, "ldir_material_filter",
                        choices = c("All", ldir_mats), selected = "All")
      updateSelectInput(session, "ldir_highlight_particle",
                        choices = c("None", ldir_ids))
      if (all(is.finite(q_range))) {
        updateSliderInput(session, "ldir_quality_range",
                          min = floor(q_range[1] * 100) / 100,
                          max = ceiling(q_range[2] * 100) / 100,
                          value = c(floor(q_range[1] * 100) / 100,
                                    ceiling(q_range[2] * 100) / 100))
      }
      if (is.finite(s_max)) {
        updateSliderInput(session, "ldir_size_range", min = 0, max = s_max,
                          value = c(0, s_max))
      }
      score_max <- if ("match_score" %in% names(ldir)) {
        ceiling(max(ldir$match_score, na.rm = TRUE) * 10) / 10
      } else NA_real_
      if (is.finite(score_max) && score_max > 0) {
        updateSliderInput(session, "ldir_score_range",
                          min = 0, max = score_max, value = c(0, score_max))
      }
      cost_max <- if ("coord_match_cost" %in% names(ldir)) {
        ceiling(max(ldir$coord_match_cost, na.rm = TRUE) * 20) / 20
      } else NA_real_
      if (is.finite(cost_max) && cost_max > 0) {
        updateSliderInput(session, "ldir_coord_cost_range",
                          min = 0, max = cost_max, value = c(0, cost_max))
      }

      # Overlay per-instrument
      updateSelectizeInput(session, "overlay_ldir_material",
                           choices = c("All", ldir_mats), selected = "All")
      updateSelectizeInput(session, "overlay_ldir_particles",
                           choices = ldir_ids, selected = character(0))
      if (all(is.finite(q_range))) {
        updateSliderInput(session, "overlay_ldir_quality",
                          min = floor(q_range[1] * 100) / 100,
                          max = ceiling(q_range[2] * 100) / 100,
                          value = c(floor(q_range[1] * 100) / 100,
                                    ceiling(q_range[2] * 100) / 100))
      }
      if (is.finite(s_max)) {
        updateSliderInput(session, "overlay_ldir_size", min = 0, max = s_max,
                          value = c(0, s_max))
      }
    }

    # --- FTIR Bruker controls (individual tab only) ---
    ftir_bruker_d <- ftir_bruker_df_full()
    if (!is.null(ftir_bruker_d) && nrow(ftir_bruker_d) > 0) {
      fb_mats  <- sort(unique(ftir_bruker_d$material))
      fb_ids   <- natural_sort_ids(unique(ftir_bruker_d$particle_id))
      q_range  <- range(ftir_bruker_d$quality, na.rm = TRUE)
      s_max    <- ceiling(max(ftir_bruker_d$feret_max, na.rm = TRUE) / 10) * 10
      updateSelectInput(session, "ftir_bruker_material_filter",
                        choices = c("All", fb_mats), selected = "All")
      updateSelectInput(session, "ftir_bruker_highlight_particle",
                        choices = c("None", fb_ids))
      if (all(is.finite(q_range))) {
        updateSliderInput(session, "ftir_bruker_quality_range",
                          min = floor(q_range[1] * 100) / 100,
                          max = ceiling(q_range[2] * 100) / 100,
                          value = c(floor(q_range[1] * 100) / 100,
                                    ceiling(q_range[2] * 100) / 100))
      }
      if (is.finite(s_max)) {
        updateSliderInput(session, "ftir_bruker_size_range", min = 0, max = s_max,
                          value = c(0, s_max))
      }
    }

    # --- Global overlay controls ---
    if (!is.null(run_data()$matched)) {
      max_dist <- ceiling(max(run_data()$matched$match_distance, na.rm = TRUE))
      updateSliderInput(session, "overlay_dist_range",
                        min = 0, max = max_dist, value = c(0, max_dist))
    }
  })

  # Global Feret Max constrains per-instrument size sliders
  observeEvent(input$overlay_size_range, {
    global <- input$overlay_size_range
    for (slider_id in c("overlay_ftir_size", "overlay_raman_size", "overlay_ldir_size")) {
      current <- input[[slider_id]]
      if (!is.null(current)) {
        new_lo <- max(current[1], global[1])
        new_hi <- min(current[2], global[2])
        updateSliderInput(session, slider_id,
                          min = global[1], max = global[2],
                          value = c(new_lo, new_hi))
      }
    }
  })

  # Pattern Apply buttons for particle selectors
  observeEvent(input$overlay_ftir_apply_pattern, {
    pat <- input$overlay_ftir_pattern
    df <- ftir_df_full()
    if (is.null(df) || nrow(df) == 0 || nchar(trimws(pat)) == 0) return()
    matched_ids <- parse_particle_selection(pat, unique(df$particle_id))
    current <- input$overlay_ftir_particles
    new_sel <- unique(c(current, matched_ids))
    updateSelectizeInput(session, "overlay_ftir_particles", selected = new_sel)
    updateTextInput(session, "overlay_ftir_pattern", value = "")
  })

  observeEvent(input$overlay_raman_apply_pattern, {
    pat <- input$overlay_raman_pattern
    df <- raman_df_full()
    if (is.null(df) || nrow(df) == 0 || nchar(trimws(pat)) == 0) return()
    matched_ids <- parse_particle_selection(pat, unique(df$particle_id))
    current <- input$overlay_raman_particles
    new_sel <- unique(c(current, matched_ids))
    updateSelectizeInput(session, "overlay_raman_particles", selected = new_sel)
    updateTextInput(session, "overlay_raman_pattern", value = "")
  })

  observeEvent(input$overlay_ldir_apply_pattern, {
    pat <- input$overlay_ldir_pattern
    df <- ldir_df_full()
    if (is.null(df) || nrow(df) == 0 || nchar(trimws(pat)) == 0) return()
    matched_ids <- parse_particle_selection(pat, unique(df$particle_id))
    current <- input$overlay_ldir_particles
    new_sel <- unique(c(current, matched_ids))
    updateSelectizeInput(session, "overlay_ldir_particles", selected = new_sel)
    updateTextInput(session, "overlay_ldir_pattern", value = "")
  })

  # Pattern Apply buttons for single-instrument viewer highlight text boxes
  # These update the selectInput highlight_particle to the parsed set.
  single_highlight_ids <- reactiveValues(ftir = NULL, raman = NULL, ldir = NULL,
                                         ftir_bruker = NULL)

  observeEvent(input$ftir_highlight_apply, {
    pat <- input$ftir_highlight_pattern
    df  <- ftir_df_full()
    if (is.null(df) || nrow(df) == 0 || is.null(pat) || nchar(trimws(pat)) == 0) return()
    ids <- parse_particle_selection(trimws(pat), unique(df$particle_id))
    single_highlight_ids$ftir <- if (length(ids) == 0) NULL else ids
    updateTextInput(session, "ftir_highlight_pattern", value = "")
  })

  observeEvent(input$raman_highlight_apply, {
    pat <- input$raman_highlight_pattern
    df  <- raman_df_full()
    if (is.null(df) || nrow(df) == 0 || is.null(pat) || nchar(trimws(pat)) == 0) return()
    ids <- parse_particle_selection(trimws(pat), unique(df$particle_id))
    single_highlight_ids$raman <- if (length(ids) == 0) NULL else ids
    updateTextInput(session, "raman_highlight_pattern", value = "")
  })

  observeEvent(input$ldir_highlight_apply, {
    pat <- input$ldir_highlight_pattern
    df  <- ldir_df_full()
    if (is.null(df) || nrow(df) == 0 || is.null(pat) || nchar(trimws(pat)) == 0) return()
    ids <- parse_particle_selection(trimws(pat), unique(df$particle_id))
    single_highlight_ids$ldir <- if (length(ids) == 0) NULL else ids
    updateTextInput(session, "ldir_highlight_pattern", value = "")
  })

  observeEvent(input$ftir_bruker_highlight_apply, {
    pat <- input$ftir_bruker_highlight_pattern
    df  <- ftir_bruker_df_full()
    if (is.null(df) || nrow(df) == 0 || is.null(pat) || nchar(trimws(pat)) == 0) return()
    ids <- parse_particle_selection(trimws(pat), unique(df$particle_id))
    single_highlight_ids$ftir_bruker <- if (length(ids) == 0) NULL else ids
    updateTextInput(session, "ftir_bruker_highlight_pattern", value = "")
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
  # Helper: generate axis breaks at every 1000 µm within a range
  breaks_1000 <- function(rng) {
    seq(floor(rng[1] / 1000) * 1000, ceiling(rng[2] / 1000) * 1000, by = 1000)
  }

  make_scatter <- function(df, img_info, bounds, title,
                            match_colours = NULL, highlight_id = NULL,
                            full_df = NULL, match_labels = NULL) {

    p <- ggplot(df, aes(x = x, y = y))

    # Background image (with per-image bounds)
    p <- add_image_bg(p, img_info)

    # Points
    p <- p + geom_point(aes(colour = match_status, size = feret_max),
                         alpha = 0.7)

    if (!is.null(match_colours)) {
      if (!is.null(match_labels))
        p <- p + scale_colour_manual(values = match_colours, labels = match_labels)
      else
        p <- p + scale_colour_manual(values = match_colours)
    }

    p <- p +
      scale_size_continuous(name = "Feret Max (\u00b5m)", range = c(2, 12)) +
      scale_x_continuous(breaks = breaks_1000(bounds$x)) +
      scale_y_continuous(breaks = breaks_1000(bounds$y)) +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = title, x = "X (\u00b5m)", y = "Y (\u00b5m)") +
      theme_minimal(base_size = 15) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "grey98", colour = NA),
        panel.grid       = element_line(colour = "grey90"),
        legend.position  = "right",
        legend.title     = element_text(size = 13),
        legend.text      = element_text(size = 11)
      )

    # Highlight selected particle(s) — ALWAYS shown even if filtered out.
    # First try the filtered df, then fall back to full_df (unfiltered).
    # highlight_id can be a character vector (multiple IDs from pattern select)
    # or a single ID (from the selectInput dropdown).
    if (!is.null(highlight_id) && length(highlight_id) > 0 &&
        !identical(highlight_id, "None") && !identical(highlight_id, character(0))) {
      hl <- NULL
      if ("particle_id" %in% names(df))
        hl <- df[df$particle_id %in% highlight_id, ]
      if ((is.null(hl) || nrow(hl) == 0) && !is.null(full_df) &&
          "particle_id" %in% names(full_df))
        hl <- full_df[full_df$particle_id %in% highlight_id, ]
      if (!is.null(hl) && nrow(hl) > 0) {
        # Y-offset scales with plot extent so label doesn't overlap the circle
        y_span <- diff(bounds$y)
        y_nudge <- y_span * 0.03   # 3% of visible y-range
        hl$label_y <- hl$y + y_nudge
        p <- p + geom_point(data = hl, aes(x = x, y = y),
                             shape = 21, size = 10, stroke = 2,
                             fill = NA, colour = "#FFD700") +
                 geom_text(data = hl, aes(x = x, y = label_y, label = particle_id),
                            vjust = 0, size = 4.0, fontface = "bold",
                            colour = "#FFD700")
      }
    }

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
  last_hover <- reactiveValues(ftir = NULL, raman = NULL, ldir = NULL, overlay = NULL,
                               ftir_bruker = NULL)

  # Pinned overlay particle: persists across hover events until cleared.
  # Stores a data row (matched or single-instrument) and its source type.
  pinned_overlay <- reactiveVal(NULL)
  pinned_source  <- reactiveVal(NULL)   # "ftir_raman", "ldir_raman", or "single_ftir"/"single_raman"/"single_ldir"

  # ==================================================================
  # Zoom state: NULL means full view, otherwise list(x=c(lo,hi), y=c(lo,hi))
  # ==================================================================
  zoom <- reactiveValues(ftir = NULL, raman = NULL, ldir = NULL, overlay = NULL,
                         ftir_bruker = NULL)

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

  observeEvent(input$ftir_bruker_brush, {
    b <- input$ftir_bruker_brush
    zoom$ftir_bruker <- list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax))
  })
  observeEvent(input$ftir_bruker_dblclick, { zoom$ftir_bruker <- NULL })

  observeEvent(input$overlay_brush, {
    b <- input$overlay_brush
    zoom$overlay <- list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax))
  })
  observeEvent(input$overlay_dblclick, { zoom$overlay <- NULL })

  # Reset-zoom buttons (same effect as double-click)
  observeEvent(input$ftir_reset_zoom,        { zoom$ftir <- NULL })
  observeEvent(input$raman_reset_zoom,       { zoom$raman <- NULL })
  observeEvent(input$ldir_reset_zoom,        { zoom$ldir <- NULL })
  observeEvent(input$ftir_bruker_reset_zoom, { zoom$ftir_bruker <- NULL })
  observeEvent(input$overlay_reset_zoom,     { zoom$overlay <- NULL })

  # ==================================================================
  # Click-to-select handlers for single-instrument viewers
  # ==================================================================

  # Helper: find nearest particle in a df (native x_orig/y_orig space)
  .find_nearest <- function(click, df, zoom_state) {
    if (is.null(click) || is.null(df) || nrow(df) == 0) return(NULL)
    vis <- if (!is.null(zoom_state)) zoom_state else
             list(x = range(df$x_orig, na.rm = TRUE),
                  y = range(df$y_orig, na.rm = TRUE))
    snap <- max(diff(vis$x), diff(vis$y), 200) * 0.05
    d <- sqrt((df$x_orig - click$x)^2 + (df$y_orig - click$y)^2)
    idx <- which.min(d)
    if (length(idx) > 0 && d[idx] <= snap) df$particle_id[idx] else NULL
  }

  observeEvent(input$ftir_click, {
    pid <- .find_nearest(input$ftir_click, ftir_filtered(), zoom$ftir)
    if (!is.null(pid)) {
      cur <- selected_ids$ftir
      selected_ids$ftir <- if (pid %in% cur) cur else c(cur, pid)
    }
  })
  observeEvent(input$ftir_clear_selection, { selected_ids$ftir <- character(0) })

  observeEvent(input$raman_click, {
    pid <- .find_nearest(input$raman_click, raman_filtered(), zoom$raman)
    if (!is.null(pid)) {
      cur <- selected_ids$raman
      selected_ids$raman <- if (pid %in% cur) cur else c(cur, pid)
    }
  })
  observeEvent(input$raman_clear_selection, { selected_ids$raman <- character(0) })

  observeEvent(input$ldir_click, {
    pid <- .find_nearest(input$ldir_click, ldir_filtered(), zoom$ldir)
    if (!is.null(pid)) {
      cur <- selected_ids$ldir
      selected_ids$ldir <- if (pid %in% cur) cur else c(cur, pid)
    }
  })
  observeEvent(input$ldir_clear_selection, { selected_ids$ldir <- character(0) })

  observeEvent(input$ftir_bruker_click, {
    pid <- .find_nearest(input$ftir_bruker_click, ftir_bruker_filtered(), zoom$ftir_bruker)
    if (!is.null(pid)) {
      cur <- selected_ids$ftir_bruker
      selected_ids$ftir_bruker <- if (pid %in% cur) cur else c(cur, pid)
    }
  })
  observeEvent(input$ftir_bruker_clear_selection, { selected_ids$ftir_bruker <- character(0) })

  # Overlay: select / deselect all layers
  observeEvent(input$overlay_toggle_all, {
    all_choices <- c("matched", "unmatched_ftir", "unmatched_raman",
                     "match_lines", "ldir_matched", "ldir_unmatched",
                     "ldir_lines", "triple_only")
    current <- input$overlay_layers
    if (length(current) == length(all_choices)) {
      updateCheckboxGroupInput(session, "overlay_layers", selected = character(0))
    } else {
      updateCheckboxGroupInput(session, "overlay_layers", selected = all_choices)
    }
  })

  # ==================================================================
  # Shared helpers: selection table + plastics summary HTML
  # ==================================================================

  make_selection_table_ui <- function(ids, full_df) {
    if (length(ids) == 0)
      return(tags$p(class = "text-muted", "Click particles to add them to the selection."))
    if (is.null(full_df) || nrow(full_df) == 0)
      return(tags$p(class = "text-muted", "Data not yet loaded."))
    rows <- full_df[full_df$particle_id %in% ids, ]
    rows <- rows[match(ids, rows$particle_id), ]
    rows <- rows[!is.na(rows$particle_id), ]
    if (nrow(rows) == 0)
      return(tags$p(class = "text-muted", "Selected particles not found in loaded data."))
    tags$table(class = "hover-tbl",
      tags$tr(tags$th("ID"), tags$th("Material"),
              tags$th("Feret (\u00b5m)"), tags$th("Quality")),
      lapply(seq_len(nrow(rows)), function(i) {
        r <- rows[i, ]
        tags$tr(tags$td(r$particle_id),
                tags$td(r$material),
                tags$td(round(r$feret_max, 1)),
                tags$td(round(r$quality, 3)))
      })
    )
  }

  make_plastics_summary_ui <- function(df) {
    tbl <- summarise_plastics(df)
    if (nrow(tbl) == 0)
      return(tags$p(class = "text-muted", "No classified materials detected."))
    # Group by category with subheadings
    ui_rows <- list()
    cur_cat <- ""
    for (i in seq_len(nrow(tbl))) {
      if (tbl$category[i] != cur_cat) {
        cur_cat <- tbl$category[i]
        ui_rows[[length(ui_rows) + 1]] <- tags$tr(
          tags$td(colspan = "2", style = "font-weight: bold; padding-top: 6px;",
                  cur_cat))
      }
      ui_rows[[length(ui_rows) + 1]] <- tags$tr(
        tags$td(style = "padding-left: 12px;", tbl$family[i]),
        tags$td(tbl$n[i]))
    }
    tags$table(class = "hover-tbl",
      tags$tr(tags$th("Family"), tags$th("Count")),
      ui_rows)
  }

  # ==================================================================
  # FTIR TAB
  # ==================================================================

  ftir_filtered <- reactive({
    df <- ftir_df_full()
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    filter_instrument(df, ftir_quality_range_d(), ftir_size_range_d(),
                      input$ftir_material_filter, input$ftir_match_filter)
  })
  
  ftir_points_df <- reactive({
    df <- ftir_df_full()
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    df %>%
      dplyr::mutate(
        x = x_orig,
        y = y_orig
      )
  })
  
  ftir_bg_path <- reactive({
    req(selected_run_manifest())
    manifest_image_path(selected_run_manifest(), "ftir_image", preferred = "canonical")
  })
  
  output$ftir_single_plot <- renderPlot({
    req(ftir_points_df())
    build_single_view_plot(
      points_df = ftir_points_df(),
      bg_png_path = ftir_bg_path()
    )
  })

  output$ftir_plot <- renderPlot({
    df <- ftir_filtered()
    # Display in native FTIR instrument frame (x_orig / y_orig).
    df_disp <- df
    if (nrow(df_disp) > 0) {
      df_disp$x <- df_disp$x_orig
      df_disp$y <- df_disp$y_orig
    }

    img <- ftir_native_image_info()

    # Full (unfiltered) FTIR data for highlight fallback and bounds when filtered to 0
    full_ftir <- ftir_df_full()
    if (!is.null(full_ftir) && nrow(full_ftir) > 0) {
      full_ftir$x <- full_ftir$x_orig
      full_ftir$y <- full_ftir$y_orig
    }

    # Bounds derived from particle coordinates — never from scan-area estimation.
    # When filters remove all particles, use the full dataset so the background
    # image remains correctly positioned.
    bounds <- if (!is.null(zoom$ftir)) zoom$ftir else {
      ref <- if (nrow(df_disp) > 0) df_disp
             else if (!is.null(full_ftir) && nrow(full_ftir) > 0) full_ftir
             else NULL
      if (!is.null(ref) && any(is.finite(ref$x))) {
        pad <- 300
        list(x = c(min(ref$x, na.rm = TRUE) - pad, max(ref$x, na.rm = TRUE) + pad),
             y = c(min(ref$y, na.rm = TRUE) - pad, max(ref$y, na.rm = TRUE) + pad))
      } else list(x = c(0, 10000), y = c(0, 10000))
    }

    if (nrow(df_disp) == 0) {
      # still show background in correct orientation/scale even with 0 points
      # use full FTIR to get bounds (or store bounds separately)
      full_ftir <- ftir_df_full()
      if (is.null(full_ftir) || nrow(full_ftir) == 0) {
        # last resort: show image without points using a default 0..1 bounds
        df0 <- data.frame(x = c(0, 1), y = c(0, 1), match_status = "none", feret_max_um = 1)
        bg <- manifest_image_path(selected_run_manifest(), "ftir_image", preferred = "canonical")
        return(build_single_view_plot(df0, bg))
      }
      
      df0 <- full_ftir %>% dplyr::mutate(
        x = x_orig,
        y = y_orig,
        match_status = "none",
        feret_max_um = 1
      )
      bg <- manifest_image_path(selected_run_manifest(), "ftir_image", preferred = "canonical")
      return(build_single_view_plot(df0, bg))
    }
    # Combine selectInput highlight + pattern-matched highlights
    hl_single <- input$ftir_highlight_particle
    hl_ids <- if (!is.null(hl_single) && hl_single != "None") {
      unique(c(hl_single, single_highlight_ids$ftir))
    } else {
      single_highlight_ids$ftir
    }
    make_scatter(df_disp, img, bounds,
                 paste0("FTIR Particles (", nrow(df_disp), " shown)"),
                 match_colours = c(matched = "#2ca02c", unmatched = "#d62728"),
                 match_labels = c(matched = "matched to Raman", unmatched = "unmatched"),
                 highlight_id = hl_ids,
                 full_df = full_ftir)
  })

  output$ftir_summary_text <- renderText({
    df <- ftir_filtered()
    if (nrow(df) == 0) return("No pipeline data loaded")
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched to Raman | ",
           length(unique(df$material)), " materials")
  })

  observeEvent(input$ftir_hover, {
    hover <- input$ftir_hover
    if (is.null(hover)) return()
    df <- ftir_filtered()
    if (nrow(df) == 0) return()
    # Search in native coordinate space (no Y-flip needed)
    dists <- sqrt((df$x_orig - hover$x)^2 + (df$y_orig - hover$y)^2)
    idx   <- which.min(dists)
    threshold <- max(diff(range(df$x_orig, na.rm = TRUE)),
                     diff(range(df$y_orig, na.rm = TRUE)), 500) * 0.05
    if (dists[idx] <= threshold) last_hover$ftir <- df[idx, , drop = FALSE]
  })

  output$ftir_hover_info <- renderUI({
    row <- last_hover$ftir
    single_detail_html(row, "FTIR (PerkinElmer)", "AAU Quality")
  })

  output$ftir_selection_info <- renderUI({
    make_selection_table_ui(selected_ids$ftir, ftir_df_full())
  })

  output$ftir_plastics_summary <- renderUI({
    make_plastics_summary_ui(ftir_filtered())
  })


  # ==================================================================
  # RAMAN TAB
  # ==================================================================

  raman_filtered <- reactive({
    df <- raman_df_full()
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    filter_instrument(df, raman_quality_range_d(), raman_size_range_d(),
                      input$raman_material_filter, input$raman_match_filter)
  })

  output$raman_plot <- renderPlot({
    df <- raman_filtered()
    # Display in native Raman instrument frame (x_orig, y_orig)
    df_disp <- df
    if (nrow(df_disp) > 0) { df_disp$x <- df_disp$x_orig; df_disp$y <- df_disp$y_orig }

    img <- raman_native_image_info()

    # Use image extent for bounds when available (consistent with overlay)
    bounds <- if (!is.null(zoom$raman)) zoom$raman else {
      if (!is.null(img)) {
        pad <- 200
        list(x = c(img$xmin - pad, img$xmax + pad),
             y = c(img$ymin - pad, img$ymax + pad))
      } else if (nrow(df_disp) > 0) {
        pad <- 300
        list(x = c(min(df_disp$x, na.rm = TRUE) - pad, max(df_disp$x, na.rm = TRUE) + pad),
             y = c(min(df_disp$y, na.rm = TRUE) - pad, max(df_disp$y, na.rm = TRUE) + pad))
      } else list(x = c(-1000, 1000), y = c(-1000, 1000))
    }

    # Full (unfiltered) Raman data for highlight fallback
    full_raman <- raman_df_full()
    if (!is.null(full_raman) && nrow(full_raman) > 0) {
      full_raman$x <- full_raman$x_orig; full_raman$y <- full_raman$y_orig
    }

    if (nrow(df_disp) == 0) {
      p <- ggplot() + coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = "Raman — no particles loaded", x = "X (\u00b5m)", y = "Y (\u00b5m)") +
        theme_minimal(base_size = 15) +
        theme(plot.background = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA))
      return(add_image_bg(p, img))
    }
    hl_single <- input$raman_highlight_particle
    hl_ids <- if (!is.null(hl_single) && hl_single != "None") {
      unique(c(hl_single, single_highlight_ids$raman))
    } else {
      single_highlight_ids$raman
    }
    make_scatter(df_disp, img, bounds,
                 paste0("Raman Particles (", nrow(df_disp), " shown)"),
                 match_colours = c(matched = "#1f77b4", unmatched = "#ff7f0e"),
                 match_labels = c(matched = "matched to FTIR", unmatched = "unmatched"),
                 highlight_id = hl_ids,
                 full_df = full_raman)
  })

  output$raman_summary_text <- renderText({
    df <- raman_filtered()
    if (nrow(df) == 0) return("No pipeline data loaded")
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched to FTIR | ",
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

  output$raman_selection_info <- renderUI({
    make_selection_table_ui(selected_ids$raman, raman_df_full())
  })

  output$raman_plastics_summary <- renderUI({
    make_plastics_summary_ui(raman_filtered())
  })


  # ==================================================================
  # LDIR TAB
  # ==================================================================

  ldir_filtered <- reactive({
    df <- ldir_df_full()
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    df <- filter_instrument(df, ldir_quality_range_d(), ldir_size_range_d(),
                            input$ldir_material_filter, input$ldir_match_filter)
    # Apply match score filter (LDIR↔Raman): keep NA (unmatched) + within range
    if ("match_score" %in% names(df) && !is.null(ldir_score_range_d())) {
      score_r <- ldir_score_range_d()
      df <- df[is.na(df$match_score) |
               (df$match_score >= score_r[1] & df$match_score <= score_r[2]), ]
    }
    # Apply coord match cost filter (image↔Excel): keep NA + within range
    if ("coord_match_cost" %in% names(df) && !is.null(ldir_coord_cost_range_d())) {
      cost_r <- ldir_coord_cost_range_d()
      df <- df[is.na(df$coord_match_cost) |
               (df$coord_match_cost >= cost_r[1] & df$coord_match_cost <= cost_r[2]), ]
    }
    df
  })

  # Processed LDIR image: background-corrected via Python (preferred)
  # or saturation mask (fallback). Shows the image after processing to
  # help diagnose extraction quality.
  ldir_processed_image <- reactive({
    raw <- ldir_raw_image()
    if (is.null(raw)) return(NULL)

    # Try Python background correction for the processed view.
    # Use the canonical LDIR path from the active manifest (no hardcoded filenames).
    ldir_img_path <- tryCatch({
      m <- active_manifest()
      cp <- m$images$ldir$canonical_path
      if (!is.null(cp) && nzchar(cp) && file.exists(cp)) cp else NULL
    }, error = function(e) NULL)
    py_ok <- tryCatch({
      if (!is.null(ldir_img_path) &&
          requireNamespace("reticulate", quietly = TRUE)) {
        py_script <- file.path("..", "inst", "python", "particle_detector.py")
        if (file.exists(py_script)) {
          if (!exists("load_and_prepare", envir = globalenv()))
            reticulate::source_python(py_script)
          data <- load_and_prepare(ldir_img_path)
          bg <- correct_background(data$gray)
          corr <- bg$corrected
          # Normalize to 0-1 for display
          mx_val <- max(corr)
          if (mx_val > 0) corr <- corr / mx_val
          # Convert to 3-channel green-tinted visualization
          h <- nrow(corr); w <- ncol(corr)
          out <- array(0.05, dim = c(h, w, 3))
          out[,,2] <- corr * 0.9       # green channel = intensity
          out[,,1] <- corr * 0.2       # slight red
          out[,,3] <- corr * 0.2       # slight blue
          return(out)
        }
      }
      FALSE
    }, error = function(e) FALSE)

    if (is.logical(py_ok) && !py_ok) {
      # Fallback: saturation mask
      if (length(dim(raw)) < 3 || dim(raw)[3] < 3) return(NULL)
      r <- raw[,,1]; g <- raw[,,2]; b <- raw[,,3]
      mx <- pmax(r, g, b)
      mn <- pmin(r, g, b)
      sat <- ifelse(mx > 0, (mx - mn) / mx, 0)
      binary <- sat > 0.3 & mx > 0.08
      h <- nrow(raw); w <- ncol(raw)
      out <- array(0.12, dim = c(h, w, 3))
      out[,,2][binary] <- 0.8
      out[,,1][binary] <- 0.15
      out[,,3][binary] <- 0.15
      return(out)
    }
    py_ok
  })

  # Processed LDIR image info (same bounds as native)
  ldir_processed_image_info <- reactive({
    proc <- ldir_processed_image()
    if (is.null(proc)) return(NULL)
    native <- ldir_native_image_info()
    if (is.null(native)) return(NULL)
    list(raster = proc,
         xmin = native$xmin, xmax = native$xmax,
         ymin = native$ymin, ymax = native$ymax)
  })

  # LDIR image-extracted particles (pre-join, from pipeline)
  ldir_extracted_pts <- reactive({
    d <- run_data()
    if (is.null(d$ldir_image_extracted) || nrow(d$ldir_image_extracted) == 0)
      return(NULL)
    d$ldir_image_extracted
  })

  output$ldir_plot <- renderPlot({
    df <- ldir_filtered()
    overlay_mode <- input$ldir_overlay_mode

    # Display in native LDIR frame (x_orig, y_orig) — no transform applied here.
    # Transforms are only used for the Overlay tab (Raman frame).
    df_disp <- df
    if (nrow(df_disp) > 0) { df_disp$x <- df_disp$x_orig; df_disp$y <- df_disp$y_orig }

    # Choose image first so we can use its bounds as the viewport.
    img <- NULL
    if ("processed_image" %in% overlay_mode) {
      img <- ldir_processed_image_info()
    }
    if (is.null(img) && "raw_image" %in% overlay_mode) {
      img <- ldir_native_image_info()
    }

    # Viewport priority:
    # 1. User zoom (brush) — always honoured
    # 2. Image bounds — when an image is shown, the viewport must cover the full
    #    scan circle; basing it on sparse particle coords distorts the image.
    # 3. Particle data range — fallback when no image is loaded.
    bounds <- if (!is.null(zoom$ldir)) zoom$ldir else if (!is.null(img)) {
      list(x = c(img$xmin, img$xmax), y = c(img$ymin, img$ymax))
    } else if (nrow(df_disp) > 0 && any(is.finite(df_disp$x))) {
      pad <- 500
      list(x = c(min(df_disp$x, na.rm = TRUE) - pad, max(df_disp$x, na.rm = TRUE) + pad),
           y = c(min(df_disp$y, na.rm = TRUE) - pad, max(df_disp$y, na.rm = TRUE) + pad))
    } else list(x = c(-7000, 7000), y = c(-7000, 7000))

    n_extracted <- 0
    extracted <- ldir_extracted_pts()
    if (!is.null(extracted)) n_extracted <- nrow(extracted)

    title_parts <- paste0("LDIR Particles (", nrow(df_disp), " shown")
    if ("extracted_pts" %in% overlay_mode && n_extracted > 0)
      title_parts <- paste0(title_parts, " + ", n_extracted, " image-extracted")
    title_parts <- paste0(title_parts, ")")

    if (nrow(df_disp) == 0 && !("extracted_pts" %in% overlay_mode && n_extracted > 0)) {
      p <- ggplot() + coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = "LDIR — no particles loaded", x = "X (\u00b5m)", y = "Y (\u00b5m)") +
        theme_minimal(base_size = 15) +
        theme(plot.background = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA))
      return(add_image_bg(p, img))
    }

    p <- ggplot() +
      scale_x_continuous(breaks = breaks_1000(bounds$x)) +
      scale_y_continuous(breaks = breaks_1000(bounds$y)) +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = title_parts, x = "X (\u00b5m)", y = "Y (\u00b5m)") +
      theme_minimal(base_size = 15) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "grey98", colour = NA),
        panel.grid       = element_line(colour = "grey90"),
        legend.position  = "right",
        legend.title     = element_text(size = 13),
        legend.text      = element_text(size = 11)
      )

    # Background image
    p <- add_image_bg(p, img)

    # Image-extracted particles (before join): small open circles
    if ("extracted_pts" %in% overlay_mode && n_extracted > 0) {
      ext_df <- data.frame(x = extracted$x_um, y = extracted$y_um,
                           feret_max = extracted$feret_max_um)
      p <- p + geom_point(data = ext_df,
                            aes(x = x, y = y, size = feret_max),
                            shape = 1, colour = "#e377c2", alpha = 0.5,
                            stroke = 0.5)
    }

    # Excel-joined particles (main layer)
    if (nrow(df_disp) > 0) {
      p <- p + geom_point(data = df_disp,
                            aes(x = x, y = y, colour = match_status,
                                size = feret_max),
                            alpha = 0.7) +
        scale_colour_manual(values = c(matched = "#d62728",
                                        unmatched = "#bcbd22"),
                             labels = c(matched = "matched to Raman",
                                        unmatched = "unmatched"))
    }

    # Size legend (single scale for all layers)
    p <- p + scale_size_continuous(name = "Feret Max (\u00b5m)", range = c(2, 12))

    # Highlight selected particle(s) — always shown even if filtered out
    hl_single <- input$ldir_highlight_particle
    hl_ids <- if (!is.null(hl_single) && hl_single != "None") {
      unique(c(hl_single, single_highlight_ids$ldir))
    } else {
      single_highlight_ids$ldir
    }

    if (!is.null(hl_ids) && length(hl_ids) > 0) {
      hl <- if (nrow(df_disp) > 0) df_disp[df_disp$particle_id %in% hl_ids, ] else data.frame()
      # Fall back to full unfiltered data (native coords) if particle is filtered out
      if (nrow(hl) == 0) {
        full_ldir <- ldir_df_full()
        if (!is.null(full_ldir) && nrow(full_ldir) > 0) {
          full_ldir$x <- full_ldir$x_orig; full_ldir$y <- full_ldir$y_orig
          hl <- full_ldir[full_ldir$particle_id %in% hl_ids, ]
        }
      }
      if (nrow(hl) > 0) {
        bounds_ldir <- if (!is.null(zoom$ldir)) zoom$ldir else {
          if (nrow(df_disp) > 0)
            list(x = range(df_disp$x, na.rm = TRUE),
                 y = range(df_disp$y, na.rm = TRUE))
          else list(x = c(0, 13000), y = c(0, 13000))
        }
        y_span   <- diff(bounds_ldir$y)
        y_nudge  <- y_span * 0.03
        hl$label_y <- hl$y + y_nudge
        p <- p + geom_point(data = hl, aes(x = x, y = y),
                             shape = 21, size = 10, stroke = 2,
                             fill = NA, colour = "#FFD700") +
                 geom_text(data = hl, aes(x = x, y = label_y, label = particle_id),
                            vjust = 0, size = 4.0, fontface = "bold",
                            colour = "#FFD700")
      }
    }

    p
  })

  output$ldir_summary_text <- renderText({
    df <- ldir_filtered()
    extracted <- ldir_extracted_pts()
    n_ext <- if (!is.null(extracted)) nrow(extracted) else 0
    if (nrow(df) == 0 && n_ext == 0) return("No LDIR data loaded")
    paste0(nrow(df), " shown | ",
           sum(df$match_status == "matched"), " matched to Raman | ",
           n_ext, " image-extracted | ",
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

  output$ldir_selection_info <- renderUI({
    make_selection_table_ui(selected_ids$ldir, ldir_df_full())
  })

  output$ldir_plastics_summary <- renderUI({
    make_plastics_summary_ui(ldir_filtered())
  })


  # ==================================================================
  # FTIR BRUKER TAB
  # ==================================================================

  ftir_bruker_filtered <- reactive({
    df <- ftir_bruker_df_full()
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    filter_instrument(df, ftir_bruker_quality_range_d(), ftir_bruker_size_range_d(),
                      input$ftir_bruker_material_filter, input$ftir_bruker_match_filter)
  })

  output$ftir_bruker_plot <- renderPlot({
    df <- ftir_bruker_filtered()
    df_disp <- df
    if (nrow(df_disp) > 0) {
      df_disp$x <- df_disp$x_orig
      df_disp$y <- df_disp$y_orig
    }

    full_fb <- ftir_bruker_df_full()
    if (!is.null(full_fb) && nrow(full_fb) > 0) {
      full_fb$x <- full_fb$x_orig
      full_fb$y <- full_fb$y_orig
    }

    bounds <- if (!is.null(zoom$ftir_bruker)) zoom$ftir_bruker else {
      ref <- if (nrow(df_disp) > 0) df_disp
             else if (!is.null(full_fb) && nrow(full_fb) > 0) full_fb
             else NULL
      if (!is.null(ref) && any(is.finite(ref$x))) {
        pad <- 300
        list(x = c(min(ref$x, na.rm = TRUE) - pad, max(ref$x, na.rm = TRUE) + pad),
             y = c(min(ref$y, na.rm = TRUE) - pad, max(ref$y, na.rm = TRUE) + pad))
      } else list(x = c(0, 10000), y = c(0, 10000))
    }

    if (nrow(df_disp) == 0) {
      return(ggplot() +
        coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = "FTIR (Bruker) \u2014 no data loaded",
             x = "X (\u00b5m)", y = "Y (\u00b5m)") +
        theme_minimal(base_size = 15) +
        theme(plot.background  = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA)))
    }

    hl_single <- input$ftir_bruker_highlight_particle
    hl_ids <- if (!is.null(hl_single) && hl_single != "None") {
      unique(c(hl_single, single_highlight_ids$ftir_bruker))
    } else {
      single_highlight_ids$ftir_bruker
    }
    make_scatter(df_disp, NULL, bounds,
                 paste0("FTIR (Bruker) Particles (", nrow(df_disp), " shown)"),
                 match_colours = c(matched = "#9467bd", unmatched = "#8c564b"),
                 match_labels = c(matched = "matched to Raman", unmatched = "unmatched"),
                 highlight_id = hl_ids,
                 full_df = full_fb)
  })

  output$ftir_bruker_summary_text <- renderText({
    df <- ftir_bruker_filtered()
    if (nrow(df) == 0) return("No FTIR (Bruker) data loaded")
    paste0(nrow(df), " particles | ", length(unique(df$material)), " materials")
  })

  observeEvent(input$ftir_bruker_hover, {
    hover <- input$ftir_bruker_hover
    if (is.null(hover)) return()
    df <- ftir_bruker_filtered()
    if (nrow(df) == 0) return()
    dists <- sqrt((df$x_orig - hover$x)^2 + (df$y_orig - hover$y)^2)
    idx   <- which.min(dists)
    threshold <- max(diff(range(df$x_orig, na.rm = TRUE)),
                     diff(range(df$y_orig, na.rm = TRUE)), 500) * 0.05
    if (dists[idx] <= threshold) last_hover$ftir_bruker <- df[idx, , drop = FALSE]
  })

  output$ftir_bruker_hover_info <- renderUI({
    row <- last_hover$ftir_bruker
    single_detail_html(row, "FTIR (Bruker)", "AAU Quality")
  })

  output$ftir_bruker_selection_info <- renderUI({
    make_selection_table_ui(selected_ids$ftir_bruker, ftir_bruker_df_full())
  })

  output$ftir_bruker_plastics_summary <- renderUI({
    make_plastics_summary_ui(ftir_bruker_filtered())
  })


  # ==================================================================
  # OVERLAY TAB
  # ==================================================================

  overlay_matched <- reactive({
    d <- run_data()
    if (is.null(d$matched) || nrow(d$matched) == 0) return(data.frame())
    df <- d$matched

    # Per-instrument quality filters (debounced)
    raman_q <- overlay_raman_quality_d()
    ftir_q  <- overlay_ftir_quality_d()
    if (!is.null(raman_q)) {
      df <- df[!is.na(df$raman_quality) &
               df$raman_quality >= raman_q[1] &
               df$raman_quality <= raman_q[2], ]
    }
    if (!is.null(ftir_q)) {
      df <- df[!is.na(df$ftir_quality) &
               df$ftir_quality >= ftir_q[1] &
               df$ftir_quality <= ftir_q[2], ]
    }

    # Per-instrument size filters (intersection with global)
    ftir_sz  <- overlay_ftir_size_d()
    raman_sz <- overlay_raman_size_d()
    if (!is.null(ftir_sz)) {
      ftir_ok <- !is.na(df$ftir_feret_max_um) &
                 df$ftir_feret_max_um >= ftir_sz[1] &
                 df$ftir_feret_max_um <= ftir_sz[2]
    } else {
      ftir_ok <- rep(TRUE, nrow(df))
    }
    if (!is.null(raman_sz)) {
      raman_ok <- !is.na(df$raman_feret_max_um) &
                  df$raman_feret_max_um >= raman_sz[1] &
                  df$raman_feret_max_um <= raman_sz[2]
    } else {
      raman_ok <- rep(TRUE, nrow(df))
    }
    df <- df[ftir_ok | raman_ok, ]

    # Match distance (debounced)
    dist_r <- overlay_dist_range_d()
    if (!is.null(dist_r)) {
      df <- df[!is.na(df$match_distance) &
               df$match_distance >= dist_r[1] &
               df$match_distance <= dist_r[2], ]
    }

    # Per-instrument material filters
    ftir_mat  <- input$overlay_ftir_material
    raman_mat <- input$overlay_raman_material
    ftir_mat_ok  <- is.null(ftir_mat) || "All" %in% ftir_mat
    raman_mat_ok <- is.null(raman_mat) || "All" %in% raman_mat
    if (!ftir_mat_ok || !raman_mat_ok) {
      keep <- rep(TRUE, nrow(df))
      if (!ftir_mat_ok)  keep <- keep & (df$ftir_material %in% ftir_mat)
      if (!raman_mat_ok) keep <- keep & (df$raman_material %in% raman_mat)
      df <- df[keep, ]
    }
    df
  })

  # LDIR-Raman matched data for overlay (filtered by per-instrument controls)
  overlay_ldir_matched <- reactive({
    d <- run_data()
    if (is.null(d$ldir_raman_matched) || nrow(d$ldir_raman_matched) == 0)
      return(data.frame())
    df <- d$ldir_raman_matched

    # LDIR quality filter
    ldir_q <- overlay_ldir_quality_d()
    if (!is.null(ldir_q) && "ldir_quality" %in% names(df)) {
      df <- df[!is.na(df$ldir_quality) &
               df$ldir_quality >= ldir_q[1] &
               df$ldir_quality <= ldir_q[2], ]
    }

    # LDIR size filter
    ldir_sz <- overlay_ldir_size_d()
    if (!is.null(ldir_sz) && "ldir_feret_max_um" %in% names(df)) {
      df <- df[!is.na(df$ldir_feret_max_um) &
               df$ldir_feret_max_um >= ldir_sz[1] &
               df$ldir_feret_max_um <= ldir_sz[2], ]
    }

    # LDIR material filter
    ldir_mat <- input$overlay_ldir_material
    if (!is.null(ldir_mat) && !("All" %in% ldir_mat) &&
        "ldir_material" %in% names(df)) {
      df <- df[df$ldir_material %in% ldir_mat, ]
    }

    df
  })

  # Triple-match data: particles detected by all three instruments
  overlay_triplets <- reactive({
    d <- run_data()
    if (is.null(d$triplets) || nrow(d$triplets) == 0) return(data.frame())
    d$triplets
  })

  output$overlay_plot <- renderPlot({
    dfs <- list(ftir = ftir_df_full(), raman = raman_df_full(), ldir = ldir_df_full())
    matched <- overlay_matched()
    ldir_m <- overlay_ldir_matched()
    triplets <- overlay_triplets()
    layers <- input$overlay_layers

    # Overlay bounds: use the image extent when available (consistent framing
    # with the Raman tab), otherwise fall back to FTIR + Raman only.
    # LDIR aligned coords are excluded because LDIR alignment is typically
    # much coarser (ICP RMS > 100 µm) and would expand the view excessively.
    bounds <- if (!is.null(zoom$overlay)) zoom$overlay else {
      img_info <- overlay_image_info()
      if (!is.null(img_info)) {
        pad <- 200
        list(x = c(img_info$xmin - pad, img_info$xmax + pad),
             y = c(img_info$ymin - pad, img_info$ymax + pad))
      } else {
        compute_bounds(dfs$ftir, dfs$raman)
      }
    }

    p <- ggplot() +
      scale_x_continuous(breaks = breaks_1000(bounds$x)) +
      scale_y_continuous(breaks = breaks_1000(bounds$y)) +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = "FTIR + Raman + LDIR Overlay (aligned coordinates)",
           x = "X (\u00b5m)", y = "Y (\u00b5m)") +
      theme_minimal(base_size = 15) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "grey98", colour = NA),
        panel.grid       = element_line(colour = "grey90"),
        legend.position  = "right",
        legend.title     = element_text(size = 13),
        legend.text      = element_text(size = 11)
      )

    # Background image: Raman microscope image placed at Raman-normalized bounds
    p <- add_image_bg(p, overlay_image_info())

    # FTIR-Raman match lines: connect each matched FTIR point to its Raman pair.
    # Visible mainly when zoomed in (good alignment = short lines).
    if ("match_lines" %in% layers && nrow(matched) > 0) {
      seg_df <- data.frame(
        x    = matched$ftir_x_aligned, y    = matched$ftir_y_aligned,
        xend = matched$raman_x_norm,     yend = matched$raman_y_norm
      )
      p <- p + geom_segment(data = seg_df,
                              aes(x = x, y = y, xend = xend, yend = yend),
                              colour = "grey40", alpha = 0.6, linewidth = 0.8)
    }

    # LDIR-Raman match lines: connect each matched LDIR point to its Raman pair.
    if ("ldir_lines" %in% layers && nrow(ldir_m) > 0 &&
        "ldir_x_aligned" %in% names(ldir_m) && "raman_x_norm" %in% names(ldir_m)) {
      ldir_seg <- data.frame(
        x    = ldir_m$ldir_x_aligned, y    = ldir_m$ldir_y_aligned,
        xend = ldir_m$raman_x_norm,    yend = ldir_m$raman_y_norm
      )
      p <- p + geom_segment(data = ldir_seg,
                              aes(x = x, y = y, xend = xend, yend = yend),
                              colour = "#9467bd", alpha = 0.5, linewidth = 0.7)
    }

    # Build combined data frame for all active layers.
    # Each instrument gets a colored circle: matched = filled, unmatched = open.
    # Size driven by feret_max only.
    all_pts <- list()

    # Matched FTIR + Raman
    if ("matched" %in% layers && nrow(matched) > 0) {
      all_pts[[length(all_pts) + 1]] <- data.frame(
        x = matched$ftir_x_aligned, y = matched$ftir_y_aligned,
        feret_max = matched$ftir_feret_max_um,
        instrument = "FTIR", match_status = "matched",
        stringsAsFactors = FALSE
      )
      all_pts[[length(all_pts) + 1]] <- data.frame(
        x = matched$raman_x_norm, y = matched$raman_y_norm,
        feret_max = matched$raman_feret_max_um,
        instrument = "Raman", match_status = "matched",
        stringsAsFactors = FALSE
      )
    }

    # Matched LDIR
    if ("ldir_matched" %in% layers && nrow(ldir_m) > 0 &&
        "ldir_x_aligned" %in% names(ldir_m)) {
      all_pts[[length(all_pts) + 1]] <- data.frame(
        x = ldir_m$ldir_x_aligned, y = ldir_m$ldir_y_aligned,
        feret_max = ldir_m$ldir_feret_max_um,
        instrument = "LDIR", match_status = "matched",
        stringsAsFactors = FALSE
      )
    }

    # Unmatched FTIR
    if ("unmatched_ftir" %in% layers && !is.null(dfs$ftir) && nrow(dfs$ftir) > 0) {
      um_f <- dfs$ftir[dfs$ftir$match_status == "unmatched", ]
      if (nrow(um_f) > 0) {
        all_pts[[length(all_pts) + 1]] <- data.frame(
          x = um_f$x, y = um_f$y, feret_max = um_f$feret_max,
          instrument = "FTIR", match_status = "unmatched",
          stringsAsFactors = FALSE
        )
      }
    }

    # Unmatched Raman
    if ("unmatched_raman" %in% layers && !is.null(dfs$raman) && nrow(dfs$raman) > 0) {
      um_r <- dfs$raman[dfs$raman$match_status == "unmatched", ]
      if (nrow(um_r) > 0) {
        all_pts[[length(all_pts) + 1]] <- data.frame(
          x = um_r$x, y = um_r$y, feret_max = um_r$feret_max,
          instrument = "Raman", match_status = "unmatched",
          stringsAsFactors = FALSE
        )
      }
    }

    # Unmatched LDIR
    if ("ldir_unmatched" %in% layers && !is.null(dfs$ldir) && nrow(dfs$ldir) > 0) {
      um_l <- dfs$ldir[dfs$ldir$match_status == "unmatched", ]
      if (nrow(um_l) > 0) {
        all_pts[[length(all_pts) + 1]] <- data.frame(
          x = um_l$x, y = um_l$y, feret_max = um_l$feret_max,
          instrument = "LDIR", match_status = "unmatched",
          stringsAsFactors = FALSE
        )
      }
    }

    # Draw all particles: filled circles for matched, open circles for unmatched.
    # Both match_status and instrument are mapped so they appear in the legend.
    if (length(all_pts) > 0) {
      both <- do.call(rbind, all_pts)
      p <- p + geom_point(data = both,
                            aes(x = x, y = y, size = feret_max,
                                colour = instrument, shape = match_status),
                            alpha = 0.65)
      p <- p +
        scale_colour_manual(
          name   = "Instrument",
          values = c(FTIR = "#2ca02c", Raman = "#1f77b4", LDIR = "#d62728")
        ) +
        scale_shape_manual(
          name   = "Match Status",
          values = c(matched = 19, unmatched = 1),
          labels = c(matched = "Matched (filled)", unmatched = "Unmatched (open)")
        )
    }

    # Triple matches: gold ring around particles detected by all three instruments
    if ("triple_only" %in% layers && nrow(triplets) > 0 &&
        nrow(matched) > 0 && nrow(ldir_m) > 0) {
      triple_pts <- list()
      m_trip <- matched[matched$raman_particle_id %in% triplets$raman_particle_id, ]
      if (nrow(m_trip) > 0) {
        triple_pts[[1]] <- data.frame(
          x = m_trip$ftir_x_aligned, y = m_trip$ftir_y_aligned)
        triple_pts[[2]] <- data.frame(
          x = m_trip$raman_x_norm, y = m_trip$raman_y_norm)
      }
      l_trip <- ldir_m[ldir_m$raman_particle_id %in% triplets$raman_particle_id, ]
      if (nrow(l_trip) > 0 && "ldir_x_aligned" %in% names(l_trip)) {
        triple_pts[[length(triple_pts) + 1]] <- data.frame(
          x = l_trip$ldir_x_aligned, y = l_trip$ldir_y_aligned)
      }
      if (length(triple_pts) > 0) {
        triple_df <- do.call(rbind, triple_pts)
        p <- p + geom_point(data = triple_df, aes(x = x, y = y),
                              shape = 21, size = 6, stroke = 1.5,
                              fill = NA, colour = "#FFD700", alpha = 0.9)
      }
    }

    # Size legend (single scale for all layers)
    p <- p + scale_size_continuous(name = "Feret Max (\u00b5m)", range = c(2, 12))

    # Highlight selected particles (multi-select per instrument).
    # ALWAYS drawn regardless of layer state — allows single-particle inspection.
    # Uses dfs$<instrument> (full unfiltered data with aligned coordinates).
    hl_specs <- list(
      list(ids = input$overlay_ftir_particles,  df = dfs$ftir,  col = "#2ca02c"),
      list(ids = input$overlay_raman_particles, df = dfs$raman, col = "#1f77b4"),
      list(ids = input$overlay_ldir_particles,  df = dfs$ldir,  col = "#d62728")
    )
    y_span_ov <- diff(bounds$y)
    y_nudge_ov <- y_span_ov * 0.03

    for (spec in hl_specs) {
      sel_ids <- spec$ids
      if (is.null(sel_ids) || length(sel_ids) == 0) next
      inst_df <- spec$df
      if (is.null(inst_df) || nrow(inst_df) == 0) next
      hl <- inst_df[inst_df$particle_id %in% sel_ids, ]
      if (nrow(hl) > 0) {
        hl$label_y <- hl$y + y_nudge_ov
        p <- p + geom_point(data = hl, aes(x = x, y = y),
                             shape = 19, size = 5, colour = spec$col) +
                 geom_point(data = hl, aes(x = x, y = y),
                             shape = 21, size = 10, stroke = 2,
                             fill = NA, colour = "#FFD700") +
                 geom_text(data = hl, aes(x = x, y = label_y, label = particle_id),
                            vjust = 0, size = 4.0, fontface = "bold",
                            colour = "#FFD700")
      }
    }

    # Also highlight pinned particle (from click or dropdown selection)
    pin <- pinned_overlay()
    if (!is.null(pin) && !is.null(pin$x) && !is.null(pin$y)) {
      pin_label_y <- pin$y + y_nudge_ov
      pin_df <- data.frame(x = pin$x, y = pin$y,
                           label_y = pin_label_y,
                           label = if (!is.null(pin$particle_id)) pin$particle_id else "")
      p <- p + geom_point(data = pin_df, aes(x = x, y = y),
                           shape = 8, size = 8, stroke = 2,
                           colour = "#FF6600") +
               geom_text(data = pin_df, aes(x = x, y = label_y, label = label),
                          vjust = 0, size = 4.6, fontface = "bold",
                          colour = "#FF6600")
    }

    p
  })

  output$overlay_summary_text <- renderText({
    m <- overlay_matched()
    dfs <- list(ftir = ftir_df_full(), raman = raman_df_full(), ldir = ldir_df_full())
    triplets <- overlay_triplets()
    if (nrow(m) == 0 && is.null(dfs$ftir)) return("No pipeline data loaded")
    n_um_f <- if (!is.null(dfs$ftir)) sum(dfs$ftir$match_status == "unmatched") else 0
    n_um_r <- if (!is.null(dfs$raman)) sum(dfs$raman$match_status == "unmatched") else 0
    n_ldir <- if (!is.null(dfs$ldir)) nrow(dfs$ldir) else 0
    n_ldir_m <- if (!is.null(dfs$ldir)) sum(dfs$ldir$match_status == "matched") else 0
    n_trip <- nrow(triplets)
    paste0(nrow(m), " FTIR\u2194Raman pairs | ",
           n_um_f, " FTIR unmatched (vs Raman) | ",
           n_um_r, " Raman unmatched (vs FTIR) | ",
           n_ldir_m, "/", n_ldir, " LDIR\u2194Raman matched | ",
           n_trip, " triple matches")
  })

  # ==================================================================
  # Helper: find nearest particle across all instruments (for click/hover)
  # Returns list(row, source, dist) or NULL.
  # Checks: matched pairs, LDIR-Raman pairs, then highlighted/selected
  # single-instrument particles (regardless of layer state).
  # ==================================================================
  find_nearest_overlay_particle <- function(px, py, snap_dist,
                                            layers, matched, ldir_m, dfs) {
    best_dist <- Inf
    best_row <- NULL
    best_source <- NULL

    # Check FTIR-Raman matched (if "matched" layer is active or always for click)
    if (!is.null(matched) && nrow(matched) > 0 &&
        ("matched" %in% layers || is.null(layers))) {
      dist_f <- sqrt((matched$ftir_x_aligned - px)^2 +
                      (matched$ftir_y_aligned - py)^2)
      dist_r <- sqrt((matched$raman_x_norm - px)^2 +
                      (matched$raman_y_norm - py)^2)
      d <- pmin(dist_f, dist_r)
      idx <- which.min(d)
      if (length(idx) > 0 && d[idx] < best_dist) {
        best_dist <- d[idx]
        best_row <- matched[idx, ]
        best_source <- "ftir_raman"
      }
    }

    # Check LDIR-Raman matched (if layer active or always for click)
    if (!is.null(ldir_m) && nrow(ldir_m) > 0 &&
        "ldir_x_aligned" %in% names(ldir_m) &&
        ("ldir_matched" %in% layers || is.null(layers))) {
      dist_l <- sqrt((ldir_m$ldir_x_aligned - px)^2 +
                      (ldir_m$ldir_y_aligned - py)^2)
      idx_l <- which.min(dist_l)
      if (length(idx_l) > 0 && dist_l[idx_l] < best_dist) {
        best_dist <- dist_l[idx_l]
        best_row <- ldir_m[idx_l, ]
        best_source <- "ldir_raman"
      }
    }

    # Check unmatched FTIR (if layer active or always for click)
    if (!is.null(dfs$ftir) && nrow(dfs$ftir) > 0 &&
        ("unmatched_ftir" %in% layers || is.null(layers))) {
      um_f <- dfs$ftir[dfs$ftir$match_status == "unmatched", ]
      if (nrow(um_f) > 0) {
        d_f <- sqrt((um_f$x - px)^2 + (um_f$y - py)^2)
        idx_f <- which.min(d_f)
        if (length(idx_f) > 0 && d_f[idx_f] < best_dist) {
          best_dist <- d_f[idx_f]
          best_row <- um_f[idx_f, , drop = FALSE]
          best_source <- "single_ftir"
        }
      }
    }

    # Check unmatched Raman (if layer active or always for click)
    if (!is.null(dfs$raman) && nrow(dfs$raman) > 0 &&
        ("unmatched_raman" %in% layers || is.null(layers))) {
      um_r <- dfs$raman[dfs$raman$match_status == "unmatched", ]
      if (nrow(um_r) > 0) {
        d_r <- sqrt((um_r$x - px)^2 + (um_r$y - py)^2)
        idx_r <- which.min(d_r)
        if (length(idx_r) > 0 && d_r[idx_r] < best_dist) {
          best_dist <- d_r[idx_r]
          best_row <- um_r[idx_r, , drop = FALSE]
          best_source <- "single_raman"
        }
      }
    }

    # Check unmatched LDIR (if layer active or always for click)
    if (!is.null(dfs$ldir) && nrow(dfs$ldir) > 0 &&
        ("ldir_unmatched" %in% layers || is.null(layers))) {
      um_l <- dfs$ldir[dfs$ldir$match_status == "unmatched", ]
      if (nrow(um_l) > 0) {
        d_l <- sqrt((um_l$x - px)^2 + (um_l$y - py)^2)
        idx_l <- which.min(d_l)
        if (length(idx_l) > 0 && d_l[idx_l] < best_dist) {
          best_dist <- d_l[idx_l]
          best_row <- um_l[idx_l, , drop = FALSE]
          best_source <- "single_ldir"
        }
      }
    }

    # ALWAYS check highlighted/selected particles regardless of layer state
    hl_specs <- list(
      list(ids = input$overlay_ftir_particles,  df = dfs$ftir,  src = "single_ftir"),
      list(ids = input$overlay_raman_particles, df = dfs$raman, src = "single_raman"),
      list(ids = input$overlay_ldir_particles,  df = dfs$ldir,  src = "single_ldir")
    )
    for (spec in hl_specs) {
      if (is.null(spec$ids) || length(spec$ids) == 0) next
      inst_df <- spec$df
      if (is.null(inst_df) || nrow(inst_df) == 0) next
      hl <- inst_df[inst_df$particle_id %in% spec$ids, ]
      if (nrow(hl) > 0) {
        d_hl <- sqrt((hl$x - px)^2 + (hl$y - py)^2)
        idx_hl <- which.min(d_hl)
        if (length(idx_hl) > 0 && d_hl[idx_hl] < best_dist) {
          best_dist <- d_hl[idx_hl]
          best_row <- hl[idx_hl, , drop = FALSE]
          best_source <- spec$src
        }
      }
    }

    if (best_dist <= snap_dist && !is.null(best_row)) {
      list(row = best_row, source = best_source, dist = best_dist)
    } else {
      NULL
    }
  }

  # Overlay: sticky hover — update last_hover$overlay only when a new match is found.
  # Checks active layers AND highlighted/selected particles regardless of layer state.
  observeEvent(input$overlay_hover, {
    hover <- input$overlay_hover
    if (is.null(hover)) return()

    layers <- input$overlay_layers
    matched <- overlay_matched()
    ldir_m <- overlay_ldir_matched()
    dfs <- list(ftir = ftir_df_full(), raman = raman_df_full(), ldir = ldir_df_full())

    vis <- if (!is.null(zoom$overlay)) zoom$overlay
           else compute_bounds(dfs$ftir, dfs$raman, dfs$ldir)
    snap_dist <- max(diff(vis$x), diff(vis$y), 500) * 0.05

    result <- find_nearest_overlay_particle(hover$x, hover$y, snap_dist,
                                            layers, matched, ldir_m, dfs)
    if (!is.null(result)) {
      last_hover$overlay <- result$row
      attr(last_hover$overlay, "source") <- result$source
    }
  })

  # ==================================================================
  # Click-to-select: clicking a particle on the overlay plot adds it to the
  # per-instrument dropdown selection (gold ring) and pins it for details.
  # Clicking an already-selected particle toggles it off.
  # ==================================================================
  observeEvent(input$overlay_click, {
    click <- input$overlay_click
    if (is.null(click)) return()

    matched <- overlay_matched()
    ldir_m <- overlay_ldir_matched()
    dfs <- list(ftir = ftir_df_full(), raman = raman_df_full(), ldir = ldir_df_full())

    vis <- if (!is.null(zoom$overlay)) zoom$overlay
           else compute_bounds(dfs$ftir, dfs$raman, dfs$ldir)
    snap_dist <- max(diff(vis$x), diff(vis$y), 500) * 0.05

    # For click, pass NULL layers to search ALL instruments
    result <- find_nearest_overlay_particle(click$x, click$y, snap_dist,
                                            NULL, matched, ldir_m, dfs)
    if (!is.null(result)) {
      pinned_overlay(result$row)
      pinned_source(result$source)

      # --- Accumulate into per-instrument dropdown selections ---
      .toggle_dropdown <- function(input_id, pid, choices) {
        cur <- input[[input_id]]
        if (pid %in% cur) {
          new_sel <- setdiff(cur, pid)
        } else {
          new_sel <- c(cur, pid)
        }
        updateSelectizeInput(session, input_id, choices = choices,
                             selected = new_sel, server = TRUE)
      }

      row <- result$row
      src <- result$source

      if (src == "ftir_raman") {
        # Matched pair: add both FTIR and Raman particle IDs
        if (!is.null(row$ftir_particle_id) && !is.null(dfs$ftir))
          .toggle_dropdown("overlay_ftir_particles", row$ftir_particle_id,
                           sort(dfs$ftir$particle_id))
        if (!is.null(row$raman_particle_id) && !is.null(dfs$raman))
          .toggle_dropdown("overlay_raman_particles", row$raman_particle_id,
                           sort(dfs$raman$particle_id))
      } else if (src == "ldir_raman") {
        if (!is.null(row$ldir_particle_id) && !is.null(dfs$ldir))
          .toggle_dropdown("overlay_ldir_particles", row$ldir_particle_id,
                           sort(dfs$ldir$particle_id))
        if (!is.null(row$raman_particle_id) && !is.null(dfs$raman))
          .toggle_dropdown("overlay_raman_particles", row$raman_particle_id,
                           sort(dfs$raman$particle_id))
      } else if (src == "single_ftir" && !is.null(row$particle_id)) {
        .toggle_dropdown("overlay_ftir_particles", row$particle_id,
                         sort(dfs$ftir$particle_id))
      } else if (src == "single_raman" && !is.null(row$particle_id)) {
        .toggle_dropdown("overlay_raman_particles", row$particle_id,
                         sort(dfs$raman$particle_id))
      } else if (src == "single_ldir" && !is.null(row$particle_id)) {
        .toggle_dropdown("overlay_ldir_particles", row$particle_id,
                         sort(dfs$ldir$particle_id))
      }
    } else {
      # Click on empty space clears the pin
      pinned_overlay(NULL)
      pinned_source(NULL)
    }
  })

  # Clear pin button
  observeEvent(input$overlay_clear_pin, {
    pinned_overlay(NULL)
    pinned_source(NULL)
  })

  # ==================================================================
  # Dropdown selection triggers details: pin the most recently added particle
  # ==================================================================
  observeEvent(input$overlay_ftir_particles, {
    sel <- input$overlay_ftir_particles
    if (is.null(sel) || length(sel) == 0) return()
    df <- ftir_df_full()
    if (is.null(df) || nrow(df) == 0) return()
    # Pin the last selected particle
    pid <- sel[length(sel)]
    row <- df[df$particle_id == pid, ]
    if (nrow(row) > 0) {
      pinned_overlay(row[1, , drop = FALSE])
      pinned_source("single_ftir")
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$overlay_raman_particles, {
    sel <- input$overlay_raman_particles
    if (is.null(sel) || length(sel) == 0) return()
    df <- raman_df_full()
    if (is.null(df) || nrow(df) == 0) return()
    pid <- sel[length(sel)]
    row <- df[df$particle_id == pid, ]
    if (nrow(row) > 0) {
      pinned_overlay(row[1, , drop = FALSE])
      pinned_source("single_raman")
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$overlay_ldir_particles, {
    sel <- input$overlay_ldir_particles
    if (is.null(sel) || length(sel) == 0) return()
    df <- ldir_df_full()
    if (is.null(df) || nrow(df) == 0) return()
    pid <- sel[length(sel)]
    row <- df[df$particle_id == pid, ]
    if (nrow(row) > 0) {
      pinned_overlay(row[1, , drop = FALSE])
      pinned_source("single_ldir")
    }
  }, ignoreNULL = FALSE)

  # ==================================================================
  # Helper: render detail HTML for a single-instrument particle
  # ==================================================================
  single_overlay_detail <- function(row, instrument, quality_label) {
    tags$table(class = "hover-tbl",
      tags$tr(tags$th("Field"), tags$th("Value")),
      tags$tr(tags$td(tags$b("Instrument")), tags$td(instrument)),
      tags$tr(tags$td(tags$b("Particle ID")), tags$td(row$particle_id)),
      tags$tr(tags$td(tags$b("Material")), tags$td(tags$b(row$material))),
      tags$tr(tags$td(tags$b(quality_label)), tags$td(round(row$quality, 3))),
      tags$tr(tags$td(tags$b("Feret Max")),
              tags$td(paste0(round(row$feret_max, 1), " \u00b5m"))),
      tags$tr(tags$td(tags$b("Match Status")), tags$td(row$match_status)),
      tags$tr(tags$td(tags$b("Position")),
              tags$td(paste0("(", round(row$x, 1), ", ", round(row$y, 1), ")")))
    )
  }

  # ==================================================================
  # Detail panel: pinned > hover. Shows pin source label when pinned.
  # ==================================================================
  output$overlay_hover_info <- renderUI({
    # Priority: pinned particle > hover
    pin <- pinned_overlay()
    pin_src <- pinned_source()

    if (!is.null(pin)) {
      row <- pin
      src <- pin_src
    } else {
      row <- last_hover$overlay
      src <- if (!is.null(row)) attr(row, "source") else NULL
    }

    if (is.null(row)) {
      return(tags$p(class = "text-muted",
                    "Hover, click, or select a particle for details"))
    }

    # Single-instrument particle
    if (!is.null(src) && grepl("^single_", src)) {
      inst <- sub("^single_", "", src)
      inst_label <- switch(inst,
                           ftir = "FTIR", raman = "Raman", ldir = "LDIR", inst)
      q_label <- switch(inst,
                        ftir = "AAU Quality", raman = "HQI", ldir = "Quality", "Quality")
      return(single_overlay_detail(row, inst_label, q_label))
    }

    # LDIR-Raman match pair
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

    # Default: FTIR-Raman match pair
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
