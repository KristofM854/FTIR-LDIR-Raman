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

# Per-material show/hide checkboxes. Choices are filled from the loaded run
# (each material family actually present, "unknown" included), so the list is
# exactly what is on the sample — see update_material_choices() server-side.
# Scrolls rather than growing the sidebar when a run carries many families.
# Badge echoing the Summary tab's scope next to each section heading, so a
# filtered figure can never be mistaken for a whole-run one. Pure UI — it
# reads the radio directly in the browser.
summary_scope_badge <- function() {
  conditionalPanel(
    condition = "input.summary_scope == 'filtered'",
    style = "display:inline;",
    tags$span(class = "label label-info",
              style = "margin-left:8px; font-weight:600; vertical-align:middle;",
              "following each tab's filters"))
}

material_filter_ui <- function(input_id, label = "Materials") {
  div(class = "mat-filter",
    div(class = "mat-filter-head",
        tags$label(label, `for` = input_id),
        span(class = "mat-filter-actions",
             actionLink(paste0(input_id, "_all"),  "all"),
             span(class = "mat-filter-sep", "·"),
             actionLink(paste0(input_id, "_none"), "none"))),
    div(class = "mat-filter-box",
        checkboxGroupInput(input_id, NULL, choices = character(0),
                           selected = character(0)))
  )
}

make_detail_row <- function(label, value) {
  tags$tr(tags$td(tags$b(label)), tags$td(value))
}

instrument_panel_ui <- function(id_prefix, quality_label, quality_min, quality_max,
                                quality_step, size_max = 1200,
                                match_choices = c("matched", "unmatched"),
                                coord_toggle = FALSE, quality_default = NULL) {
  if (is.null(quality_default)) quality_default <- c(quality_min, quality_max)

  sidebarLayout(
    sidebarPanel(width = 3,
      h4(paste0(toupper(id_prefix), " Filters")),
      div(id = paste0(id_prefix, "_tour_quality"),
        sliderInput(paste0(id_prefix, "_quality_range"), quality_label,
                    min = quality_min, max = quality_max,
                    value = quality_default, step = quality_step)),
      div(id = paste0(id_prefix, "_tour_size"),
        sliderInput(paste0(id_prefix, "_size_range"), "Feret Max (\u00b5m)",
                    min = 0, max = size_max,
                    value = c(DEFAULT_MIN_SIZE_UM, size_max), step = 5)),
      div(id = paste0(id_prefix, "_tour_material"),
        material_filter_ui(paste0(id_prefix, "_material_filter"), "Materials")),
      div(id = paste0(id_prefix, "_tour_match"),
        checkboxGroupInput(paste0(id_prefix, "_match_filter"), "Match Status",
                           choices = match_choices,
                           selected = unname(match_choices), inline = TRUE),
        # ON by default: a viewer should open showing every particle the
        # instrument detected. Tick it off to colour the points by match status
        # and honour the Match Status boxes above.
        checkboxInput(paste0(id_prefix, "_show_all_detected"),
                      "Show all detected (ignore match status)",
                      value = TRUE),
        checkboxInput(paste0(id_prefix, "_show_all_labels"),
                      "Number all particles",
                      value = FALSE),
        conditionalPanel(
          condition = sprintf("input.%s_show_all_labels", id_prefix),
          sliderInput(paste0(id_prefix, "_label_size"), "Label size",
                      min = 2, max = 10, value = 3, step = 0.5))),
      # View orientation (native mode only) — rotate/mirror the whole native
      # scene into the Raman orientation, so this tab can be read side by side
      # with the Raman tab and the overlay. Display-only; nothing stored moves.
      if (coord_toggle) div(
        hr(),
        h4("View Orientation"),
        selectInput(paste0(id_prefix, "_view_rotation"), "View rotation (native mode)",
                    choices = c("Auto (match Raman)" = "auto",
                                "None (0°)"     = "0",
                                "90° counter-clockwise" = "90",
                                "90° clockwise" = "-90",
                                "180°"          = "180"),
                    selected = "auto"),
        checkboxInput(paste0(id_prefix, "_view_flip_y"), "Mirror (flip Y)",
                      value = FALSE),
        helpText(style = "font-size:11px;",
                 "Auto measures the rotation/mirror that brings this instrument's ",
                 "particles onto the Raman ones. A mirror cannot be undone by ",
                 "rotation alone, so try it if 180° still looks wrong.")
      ),
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
      if (coord_toggle) div(
        style = "margin-bottom: 6px;",
        radioButtons(paste0(id_prefix, "_coord_mode"), NULL,
                     choices = c("Native coordinates" = "native",
                                 "Aligned (Raman space)" = "aligned"),
                     selected = "native", inline = TRUE)
      ),
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
    /* Material show/hide list: header row (title + all/none) above a boxed,
       scrollable column of checkboxes on one consistent indent. */
    .mat-filter { margin-bottom: 14px; }
    .mat-filter-head { display: flex; align-items: baseline;
                       justify-content: space-between; margin-bottom: 4px;
                       gap: 8px; }
    .mat-filter-head > label { margin: 0; font-weight: 700; }
    .mat-filter-actions { font-size: 11px; white-space: nowrap; }
    .mat-filter-actions a { color: #6c757d; text-decoration: none; }
    .mat-filter-actions a:hover,
    .mat-filter-actions a:focus { color: #337ab7; text-decoration: underline; }
    .mat-filter-sep { color: #ccc; margin: 0 3px; }
    .mat-filter-box { max-height: 170px; overflow-y: auto; background: #fff;
                      border: 1px solid #ced4da; border-radius: 4px;
                      padding: 5px 10px; }
    /* Bootstrap gives .form-group a bottom margin; inside the box that is
       just dead space under the last row. */
    .mat-filter-box .form-group { margin-bottom: 0; }
    .mat-filter-box .shiny-options-group { margin-top: 0; }
    .mat-filter-box .checkbox { margin: 0; }
    .mat-filter-box .checkbox label { display: block; font-weight: 400;
                                      font-size: 13px; line-height: 1.35;
                                      padding: 4px 4px 4px 24px;
                                      border-radius: 3px; cursor: pointer; }
    .mat-filter-box .checkbox label:hover { background: #eef2f6; }
    .mat-filter-box .checkbox input[type=\"checkbox\"] { margin-left: -24px;
                                                       margin-top: 5px; }
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
                            "Unmatched (vs Raman)" = "unmatched"),
          coord_toggle = TRUE, quality_default = c(0.7, 1)))
    ),

    # Tab 2: FTIR (Bruker) — shown only when data present
    tabPanel("FTIR (Bruker)",
      div(id = "ftir_bruker_viewer",
        instrument_panel_ui("ftir_bruker", "AAU Quality", 0, 1, 0.01, 800,
          match_choices = c("Matched \u2194 Raman" = "matched",
                            "Unmatched (vs Raman)" = "unmatched"),
          coord_toggle = TRUE, quality_default = c(0.7, 1)))
    ),

    # Tab 3: Raman
    tabPanel("Raman",
      div(id = "raman_viewer",
        instrument_panel_ui("raman", "HQI", 0, 100, 1, 1200,
          match_choices = c("Matched \u2194 FTIR" = "matched",
                            "Unmatched (vs FTIR)" = "unmatched"),
          quality_default = c(70, 100)))
    ),

    # Tab 3: LDIR
    tabPanel("LDIR",
      div(id = "ldir_viewer",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("LDIR Filters"),
          sliderInput("ldir_quality_range", "Quality",
                      min = 0, max = 1, value = c(0.8, 1), step = 0.01),
          sliderInput("ldir_size_range", "Feret Max (\u00b5m)",
                      min = 0, max = 1200,
                      value = c(DEFAULT_MIN_SIZE_UM, 1200), step = 5),
          material_filter_ui("ldir_material_filter", "Materials"),
          # ON by default, same as the other instrument tabs.
          checkboxInput("ldir_show_all_detected",
                        "Show all detected (ignore match status)",
                        value = TRUE),
          checkboxInput("ldir_show_all_labels",
                        "Number all particles",
                        value = FALSE),
          conditionalPanel(
            condition = "input.ldir_show_all_labels",
            sliderInput("ldir_label_size", "Label size",
                        min = 2, max = 10, value = 3, step = 0.5)),
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
          h4("Coordinate System"),
          radioButtons("ldir_coord_mode", NULL,
                       choices = c("Native coordinates" = "native",
                                   "Aligned (Raman space)" = "aligned"),
                       selected = "native", inline = TRUE),
          selectInput("ldir_view_rotation", "View rotation (native mode)",
                      choices = c("Auto (match Raman)" = "auto",
                                  "None (0°)"     = "0",
                                  "90° counter-clockwise" = "90",
                                  "90° clockwise" = "-90",
                                  "180°"          = "180"),
                      selected = "auto"),
          checkboxInput("ldir_view_flip_y", "Mirror (flip Y)", value = FALSE),
          hr(),
          h4("Image Overlay"),
          checkboxGroupInput("ldir_overlay_mode", "Display",
                             choices = c("Raw image" = "raw_image",
                                         "Processed image" = "processed_image",
                                         "Image-extracted particles" = "extracted_pts"),
                             selected = c("raw_image"),
                             inline = FALSE),
          # Background image source — lets you place the LDIR points over the
          # LDIR image (native), the Raman image (to test whether LDIR points
          # land on the Raman particles), or an uploaded image.
          selectInput("ldir_bg_image", "Background image",
                      choices = c("Auto (LDIR native / Raman aligned)" = "auto",
                                  "LDIR image" = "ldir",
                                  "Raman image" = "raman",
                                  "None" = "none"),
                      selected = "auto"),
          # Readable-overlay aids for LDIR<->Raman matches
          checkboxInput("ldir_show_raman_partners",
                        "Show Raman partners + match lines (aligned)", value = FALSE),
          checkboxInput("ldir_hide_unmatched",
                        "Hide unmatched (single-instrument) particles", value = FALSE),
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
          # --- DISPLAY CONTROLS (top) ---
          fluidRow(
            column(6,
              tags$p(tags$strong("INSTRUMENTS"),
                     tags$br(),
                     tags$span("(which to show)", style = "font-size:11px; color:#888;")),
              checkboxGroupInput("overlay_instruments", NULL,
                                 choices = c("FTIR (PerkinElmer)" = "ftir_pe",
                                             "FTIR (Bruker)"      = "ftir_bruker",
                                             "Raman"              = "raman",
                                             "LDIR"               = "ldir"),
                                 selected = c("ftir_pe", "ftir_bruker", "raman", "ldir"))
            ),
            column(6,
              tags$p(tags$strong("RELATIONSHIPS"),
                     tags$br(),
                     tags$span("(how to show them)", style = "font-size:11px; color:#888;")),
              checkboxGroupInput("overlay_relationships", NULL,
                                 choices = c("Matched particles"           = "matched",
                                             "Unmatched particles"         = "unmatched",
                                             "Match lines"                 = "lines",
                                             "Multi-instrument (3+/4)"     = "multi"),
                                 selected = c("matched", "unmatched"))
            )
          ),
          actionLink("overlay_toggle_all", "Select / Deselect All instruments",
                     style = "font-size:11px; margin-bottom:4px; display:block;"),
          checkboxInput("overlay_show_all_labels", "Number all particles",
                        value = FALSE),
          conditionalPanel(
            condition = "input.overlay_show_all_labels",
            sliderInput("overlay_label_size", "Label size",
                        min = 2, max = 10, value = 3, step = 0.5)),
          hr(),
          div(class = "info-box",
              h5("Match Summary"), textOutput("overlay_summary_text")),
          hr(),
          fileInput("overlay_image_upload", "Background Image",
                    accept = c("image/png", "image/jpeg")),
          fluidRow(
            column(6, numericInput("overlay_img_offset_x",
                                   "Img X offset (µm)", value = 0, step = 25)),
            column(6, numericInput("overlay_img_offset_y",
                                   "Img Y offset (µm)", value = 0, step = 25))
          ),
          hr(),

          # --- GLOBAL FILTERS ---
          h4("Global Filters"),
          sliderInput("overlay_size_range", "Feret Max (\u00b5m)",
                      min = 0, max = 1200,
                      value = c(DEFAULT_MIN_SIZE_UM, 1200), step = 5),
          sliderInput("overlay_dist_range", "Match Distance (\u00b5m)",
                      min = 0, max = 100, value = c(0, 100), step = 1),
          hr(),

          # --- FT-IR SECTION ---
          h4("FT-IR", style = "color: #2ca02c; margin-bottom: 4px;"),
          sliderInput("overlay_ftir_quality", "AAU Quality",
                      min = 0, max = 1, value = c(0.7, 1), step = 0.01),
          sliderInput("overlay_ftir_size", "Feret Max (\u00b5m)",
                      min = 0, max = 800,
                      value = c(DEFAULT_MIN_SIZE_UM, 800), step = 5),
          material_filter_ui("overlay_ftir_material", "Materials"),
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
                      min = 0, max = 100, value = c(70, 100), step = 1),
          sliderInput("overlay_raman_size", "Feret Max (\u00b5m)",
                      min = 0, max = 1200,
                      value = c(DEFAULT_MIN_SIZE_UM, 1200), step = 5),
          material_filter_ui("overlay_raman_material", "Materials"),
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
          # Live acceptance gate: a LDIR<->Raman pair counts as matched only when
          # its aligned-coordinate distance is at/under this value. Seeded from
          # the run's recorded gate; drag to retune matched vs unmatched in real
          # time (summary, overlay and triplet counts all follow).
          sliderInput("overlay_ldir_dist_gate", "Match Gate (µm, LDIR↔Raman)",
                      min = 0, max = 500, value = 250, step = 5),
          sliderInput("overlay_ldir_quality", "Quality",
                      min = 0, max = 1, value = c(0.8, 1), step = 0.01),
          sliderInput("overlay_ldir_size", "Feret Max (\u00b5m)",
                      min = 0, max = 1200,
                      value = c(DEFAULT_MIN_SIZE_UM, 1200), step = 5),
          material_filter_ui("overlay_ldir_material", "Materials"),
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

          # --- FTIR BRUKER SECTION ---
          h4("FT-IR (Bruker)", style = "color: #9467bd; margin-bottom: 4px;"),
          sliderInput("overlay_ftir_bruker_quality", "AAU Quality",
                      min = 0, max = 1, value = c(0.7, 1), step = 0.01),
          sliderInput("overlay_ftir_bruker_size", "Feret Max (\u00b5m)",
                      min = 0, max = 800,
                      value = c(DEFAULT_MIN_SIZE_UM, 800), step = 5),
          material_filter_ui("overlay_ftir_bruker_material", "Materials"),
          fluidRow(
            column(8, textInput("overlay_ftir_bruker_pattern", NULL,
                                placeholder = "Range (1-10) or pattern (MP_*)")),
            column(4, actionButton("overlay_ftir_bruker_apply_pattern", "Apply",
                                   class = "btn-sm", style = "margin-top: 25px;"))
          ),
          selectizeInput("overlay_ftir_bruker_particles", "Highlight Particles",
                         choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Select particles...",
                                        plugins = list("remove_button"))),
          hr()

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
          # One scope control for the WHOLE tab, so the numbers here can never
          # silently disagree with what the instrument tabs are showing.
          div(class = "info-box", style = "margin-top: 20px;",
            fluidRow(
              column(5,
                radioButtons("summary_scope", "Particles included",
                             choices = c("All particles in the run" = "all",
                                         "Only those passing each tab's filters" = "filtered"),
                             selected = "filtered")),
              column(7, p(class = "text-muted", style = "margin-top: 26px;",
                "Applies to every panel on this tab — the comparison bar chart, ",
                "the material table, the pie charts and the size distributions. ",
                "Pick the second option to make this tab follow the quality / ",
                "size / material / match-status filters you set on the FTIR, ",
                "Raman, Bruker and LDIR tabs."))
            )
          ),
          div(class = "info-box", style = "margin-top: 12px;",
            h4("Report", style = "margin-top: 0;"),
            fluidRow(
              column(4,
                textInput("report_analyst_name", "Analyst name (optional)",
                          placeholder = "e.g. J. Smith")),
              column(3, style = "padding-top: 25px;",
                downloadButton("download_report", "Download report",
                               class = "btn-primary")),
              column(5, p(class = "text-muted", style = "margin-top: 26px;",
                "Downloads a ZIP containing a PDF and an interactive HTML report. ",
                "Both reproduce the viewer exactly as it is right now \u2014 the scope ",
                "above and the filters on each tab. The HTML includes the interactive ",
                "material bar chart; the PDF is suitable for printing and archiving."))
            )
          ),
          hr(),
          div(class = "info-box", style = "margin-top: 20px;",
            h4("Material Comparison Across Instruments", summary_scope_badge()),
            fluidRow(
              # Dropdown left; toggles stacked in order on the right
              column(4,
                selectInput("summary_material_select", "Material Family",
                            choices = list(
                              "Overview"     = c("All Plastics (stacked)" = "__all_plastics__"),
                              "Plastics"     = c("PE"),
                              "Non-Plastics" = character(0)
                            ),
                            selected = "__all_plastics__")),
              column(8,
                tags$div(style = "font-size: 15px;",
                  radioButtons("summary_bar_materials_mode", NULL,
                               choices = c("Plastics only" = "plastics",
                                           "All particles" = "all"),
                               selected = "plastics", inline = TRUE),
                  conditionalPanel(
                    condition = "input.summary_material_select === '__all_plastics__' && input.summary_bar_materials_mode === 'plastics'",
                    radioButtons("summary_bar_category_mode", NULL,
                                 choices = c("Synthetic only"       = "synthetic",
                                             "Synth. + Semi-synth." = "both"),
                                 selected = "both", inline = TRUE)),
                  radioButtons("summary_bar_display_mode", NULL,
                               choices = c("Absolute counts" = "abs",
                                           "Relative (%)"    = "rel"),
                               selected = "abs", inline = TRUE)
                )
              )
            ),
            p(class = "text-muted",
              "Hover over bars for exact counts. Click legend entries to show/hide families.",
              "Select a specific family from the dropdown to compare it across instruments."),
            plotly::plotlyOutput("summary_material_barplot", height = "500px")
          ),
          hr(),
          div(class = "info-box",
            h4("Plastics by Instrument", summary_scope_badge()),
            p(class = "text-muted",
              "Material family counts per device."),
            uiOutput("summary_plastics_wide")
          ),
          hr(),
          div(class = "info-box",
            h4("Size Distribution by Instrument", summary_scope_badge()),
            p(class = "text-muted",
              "Histogram of particle Feret Max (µm) per instrument.",
              "Solid bars: matched particles. Outline bars: unmatched."),
            fluidRow(
              column(4, plotOutput("size_hist_ftir",    height = "280px")),
              column(4, plotOutput("size_hist_raman",   height = "280px")),
              column(4, plotOutput("size_hist_ldir",    height = "280px"))
            ),
            hr(),
            h5("Size Statistics"),
            uiOutput("size_stats_table")
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
    ),

    # Tab 8: Multi-Run reproducibility view
    tabPanel("Multi-Run",
      div(id = "multirun_panel",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Reproducibility (Multi-Run)"),
          p(class = "text-muted",
            "Overlay repeat runs of one filter on one instrument, produced by ",
            code("tools/reproducibility.R"), "."),
          textInput("repro_base", "Reproducibility folder",
                    value = "C:/Users/moellerkr/OneDrive - IAEA/My Documents/Automatisations/FTIR-LDIR-Raman/output/reproducibility/"),
          div(style = "margin-bottom: 8px;",
              actionButton("repro_refresh", "Scan", class = "btn-sm btn-primary",
                           icon = icon("refresh")),
              uiOutput("repro_browse_ui", inline = TRUE)),
          uiOutput("repro_run_select_ui"),
          hr(),
          checkboxInput("repro_only_nonrepro",
                        "Show only non-reproducible particles", value = FALSE),
          checkboxInput("repro_show_missing_ring",
                        "Show cyan ring (missing in ≥1 run)", value = TRUE),
          checkboxInput("repro_show_discordant_ring",
                        "Show pink ring (material disagreement)", value = TRUE),
          checkboxGroupInput("repro_run_visibility", "Show Runs",
                             choices = character(0), selected = character(0),
                             inline = TRUE),
          checkboxInput("repro_show_image", "Show background image", value = TRUE),
          fileInput("repro_bg_upload", "Background image (optional)",
                    accept = c("image/png", "image/jpeg", "image/tiff",
                               ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp")),
          selectInput("repro_img_rotation", "Rotate background image",
                      choices = c("None (0°)" = "0",
                                  "90° counter-clockwise" = "90",
                                  "90° clockwise" = "-90",
                                  "180°" = "180"),
                      selected = "0"),
          helpText("The image is placed automatically from the run's recorded ",
                   "metadata (matches the single-instrument tab). Only set ",
                   "width/height below to override that, or when no metadata ",
                   "was recorded."),
          fluidRow(
            column(6, numericInput("repro_img_width_um", "Width (µm)",
                                   value = NA, min = 0, step = 100)),
            column(6, numericInput("repro_img_height_um", "Height (µm)",
                                   value = NA, min = 0, step = 100))
          ),
          fluidRow(
            column(6, numericInput("repro_img_offset_x", "X offset (µm)",
                                   value = 0, step = 100)),
            column(6, numericInput("repro_img_offset_y", "Y offset (µm)",
                                   value = 0, step = 100))
          ),
          uiOutput("repro_quality_ui"),
          uiOutput("repro_coord_cost_ui"),
          checkboxInput("repro_show_lines", "Link instances across runs", value = TRUE),
          selectizeInput("repro_material", "Material",
                         choices = c("All"), selected = "All", multiple = TRUE),
          hr(),
          uiOutput("repro_summary_ui")
        ),
        mainPanel(width = 9,
          plotOutput("repro_plot", height = "720px",
                     hover    = hoverOpts("repro_hover", delay = 100,
                                         delayType = "throttle"),
                     brush    = brushOpts("repro_brush", resetOnNew = TRUE),
                     dblclick = "repro_dblclick"),
          fluidRow(
            column(10, tags$p(class = "text-muted",
                   "Drag to zoom in. Double-click to reset.")),
            column(2, actionButton("repro_reset_zoom", "Reset Zoom",
                                   class = "btn-sm btn-default",
                                   style = "float:right; margin-top:2px;"))
          ),
          hr(),
          div(class = "info-box",
              h5("Particle Details (hover)"),
              detail_table_ui("repro_hover_info")),
          uiOutput("repro_ldir_circle_ui"),
          hr(),
          div(class = "info-box",
              h5("Particle Counts per Run"),
              p(class = "text-muted",
                "Total particles detected in each run (independent of the ",
                "material filter, non-reproducible-only toggle, and Show ",
                "Runs selection above)."),
              plotOutput("repro_counts_barplot", height = "280px")),
          hr(),
          div(class = "info-box",
              h5("Polymer Composition per Run"),
              p(class = "text-muted", "Non-plastic materials excluded."),
              plotOutput("repro_composition_pies", height = "320px"))
        )
      )
      ) # end div#multirun_panel
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
             intro   = "The <strong>Materials</strong> list has one checkbox per material family found on this sample — polymers plus <em>unknown</em>. Untick any to hide those particles; the <em>all</em> / <em>none</em> links toggle the whole list. Every viewer and the Overlay tab has its own copy."),
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
    plastic_fams     <- c(synthetic_families, semi_synthetic_families)
    plastics_present <- all_mats[all_mats %in% plastic_fams]
    nonplastics_present <- all_mats[!all_mats %in% plastic_fams]
    grouped_choices <- list(
      "Overview"     = c("All Plastics (stacked)" = "__all_plastics__"),
      "Plastics"     = plastics_present,
      "Non-Plastics" = nonplastics_present
    )
    # Remove empty groups
    grouped_choices <- grouped_choices[vapply(grouped_choices, length, integer(1)) > 0]
    updateSelectInput(session, "summary_material_select",
                      choices = grouped_choices, selected = "__all_plastics__")
  })

  # ------------------------------------------------------------------
  # Summary tab data source — ONE reactive behind every panel on the tab
  # ------------------------------------------------------------------
  # input$summary_scope switches the whole tab between the full run and each
  # instrument's currently-filtered set, so a number here always matches what
  # the corresponding tab is showing. Everything downstream reads summary_dfs()
  # rather than reaching for *_df_full() or *_filtered() directly — that split
  # is what let the size distributions ignore the toggle.
  SUMMARY_DEVICES <- c("FTIR (PerkinElmer)" = "ftir", "FTIR (Bruker)" = "ftir_bruker",
                       "Raman" = "raman", "LDIR" = "ldir")

  summary_filtered <- reactive(identical(input$summary_scope, "filtered"))

  summary_dfs <- reactive({
    if (summary_filtered())
      list(ftir        = ftir_filtered(),
           ftir_bruker = ftir_bruker_filtered(),
           raman       = raman_filtered(),
           ldir        = ldir_filtered())
    else {
      dfs <- instrument_dfs()
      list(ftir        = dfs$ftir,
           ftir_bruker = dfs$ftir_bruker,
           raman       = dfs$raman,
           ldir        = dfs$ldir)
    }
  })

  # Per-instrument material family counts, named by device label.
  active_material_counts <- reactive({
    d <- summary_dfs()
    lapply(SUMMARY_DEVICES, function(key) {
      x <- d[[key]]
      if (is.null(x) || nrow(x) == 0 || !"material" %in% names(x)) return(NULL)
      table(classify_family_vec(x$material))
    })
  })

  # Interactive plotly barplot: material families across instruments
  summary_material_barplot_obj <- reactive({
    sel_fam   <- input$summary_material_select
    rel_mode  <- identical(input$summary_bar_display_mode, "rel")
    show_all  <- identical(input$summary_bar_materials_mode, "all")
    cat_mode  <- input$summary_bar_category_mode %||% "both"
    if (is.null(sel_fam) || !nzchar(sel_fam)) return(NULL)

    device_colors <- c("FTIR (PerkinElmer)" = "#2ca02c", "FTIR (Bruker)" = "#9467bd",
                       "Raman" = "#1f77b4", "LDIR" = "#d62728")
    cts      <- active_material_counts()
    use_filt <- summary_filtered()
    has_data <- !vapply(names(cts), function(l) is.null(cts[[l]]), logical(1))
    cts_present <- cts[has_data]
    inst_levels <- names(cts_present)

    if (length(cts_present) == 0)
      return(plotly::plot_ly() |>
               plotly::layout(title = "No data available",
                              xaxis = list(visible = FALSE),
                              yaxis = list(visible = FALSE)))

    # Determine which families to include
    plastic_fams <- c(synthetic_families, semi_synthetic_families)
    if (sel_fam == "__all_plastics__") {
      keep_fams <- if (show_all) NULL   # NULL = all families
                   else if (identical(cat_mode, "synthetic")) synthetic_families
                   else plastic_fams
    } else {
      keep_fams <- sel_fam
    }

    # Build long data frame
    rows <- do.call(rbind, lapply(inst_levels, function(dev_label) {
      tbl <- cts_present[[dev_label]]
      fams_here <- names(tbl)
      if (!is.null(keep_fams)) fams_here <- fams_here[fams_here %in% keep_fams]
      if (length(fams_here) == 0) return(NULL)
      data.frame(instrument = dev_label, family = fams_here,
                 count = as.integer(tbl[fams_here]), stringsAsFactors = FALSE)
    }))

    if (is.null(rows) || nrow(rows) == 0)
      return(plotly::plot_ly() |>
               plotly::layout(title = "No matching data",
                              xaxis = list(visible = FALSE),
                              yaxis = list(visible = FALSE)))

    # Compute totals per instrument for relative mode
    inst_totals <- tapply(rows$count, rows$instrument, sum)
    rows$value  <- if (rel_mode)
      round(rows$count / inst_totals[rows$instrument] * 100, 1)
    else
      rows$count
    rows$tooltip <- if (rel_mode)
      paste0(rows$family, ": ", rows$value, "% (n=", rows$count, ")")
    else
      paste0(rows$family, ": ", rows$value)

    y_label <- if (rel_mode) "Share (%)" else "Particle Count"
    title_str <- paste0(
      if (sel_fam == "__all_plastics__") "All Plastics" else sel_fam,
      " across instruments", if (use_filt) " (filtered)" else "")

    # Family stacking order: largest total first (bottom of stack), smallest last (top).
    # Computed from the current value column so it respects abs/rel mode.
    fam_totals <- if (sel_fam == "__all_plastics__") {
      fams_present <- unique(rows$family)
      totals <- vapply(fams_present, function(f) sum(rows$value[rows$family == f]), numeric(1))
      names(totals) <- fams_present
      # largest at bottom = drawn first in plotly stack
      sort(totals, decreasing = TRUE)
    } else {
      setNames(sum(rows$value), sel_fam)
    }
    fam_order <- names(fam_totals)

    # Colour palette: family colours for stacked, device colours for single-family
    .fam_palette <- c(
      PE = "#e41a1c", PP = "#377eb8", PS = "#4daf4a", PET = "#984ea3",
      PVC = "#ff7f00", PA = "#a65628", PU = "#f781bf", PC = "#999999",
      PMMA = "#66c2a5", PTFE = "#fc8d62", ABS = "#e78ac3", Rubber = "#7570b3",
      Cellulose = "#bcbd22", Acrylate = "#17becf", Other = "#e5c494"
    )

    p <- plotly::plot_ly()
    for (fam in fam_order) {
      sub <- rows[rows$family == fam, ]
      # Align to all instruments so every trace spans the same x positions
      sub_aligned <- data.frame(instrument = inst_levels, stringsAsFactors = FALSE)
      sub_aligned <- merge(sub_aligned, sub, by = "instrument", all.x = TRUE)
      sub_aligned$value[is.na(sub_aligned$value)]     <- 0
      sub_aligned$tooltip[is.na(sub_aligned$tooltip)] <- paste0(fam, ": 0")
      # Show label only when segment is large enough to be readable
      bar_label <- ifelse(sub_aligned$value > 0, as.character(sub_aligned$value), "")

      p <- plotly::add_trace(p,
        x             = sub_aligned$instrument,
        y             = sub_aligned$value,
        type          = "bar",
        name          = fam,
        text          = bar_label,
        textposition  = "inside",
        insidetextanchor = "middle",
        textfont      = list(size = 13, color = "white"),
        hovertext     = sub_aligned$tooltip,
        hoverinfo     = "text",
        marker        = list(color = if (sel_fam == "__all_plastics__")
                               (.fam_palette[fam] %||% "#cccccc")
                             else unname(device_colors[sub_aligned$instrument]))
      )
    }

    plotly::layout(p,
      barmode = if (sel_fam == "__all_plastics__") "stack" else "group",
      title   = list(text = title_str, x = 0.5, xanchor = "center",
                     font = list(size = 17)),
      xaxis   = list(title = "", tickfont = list(size = 14)),
      yaxis   = list(title = y_label,
                     titlefont = list(size = 14), tickfont = list(size = 13)),
      legend  = list(title  = list(text = "<b>Family</b>", font = list(size = 14)),
                     font   = list(size = 13)),
      margin  = list(t = 65, r = 20, b = 50, l = 65),
      uniformtext = list(minsize = 10, mode = "hide")
    )
  })

  output$summary_material_barplot <- plotly::renderPlotly(summary_material_barplot_obj()) |>
    bindCache(active_material_counts(), input$summary_material_select,
              input$summary_bar_display_mode, input$summary_bar_materials_mode,
              input$summary_bar_category_mode)

  # Static ggplot version of the barplot used exclusively by the PDF report
  # (report_figure_page expects a ggplot, not a plotly widget).
  summary_material_barplot_gg <- reactive({
    sel_fam  <- input$summary_material_select
    rel_mode <- identical(input$summary_bar_display_mode, "rel")
    show_all <- identical(input$summary_bar_materials_mode, "all")
    cat_mode <- input$summary_bar_category_mode %||% "both"
    if (is.null(sel_fam) || !nzchar(sel_fam)) return(NULL)

    .fam_palette <- c(
      PE = "#e41a1c", PP = "#377eb8", PS = "#4daf4a", PET = "#984ea3",
      PVC = "#ff7f00", PA = "#a65628", PU = "#f781bf", PC = "#999999",
      PMMA = "#66c2a5", PTFE = "#fc8d62", ABS = "#e78ac3", Rubber = "#7570b3",
      Cellulose = "#bcbd22", Acrylate = "#17becf", Other = "#e5c494"
    )
    device_colors <- c("FTIR (PerkinElmer)" = "#2ca02c", "FTIR (Bruker)" = "#9467bd",
                       "Raman" = "#1f77b4", "LDIR" = "#d62728")
    cts      <- active_material_counts()
    use_filt <- summary_filtered()
    has_data <- !vapply(names(cts), function(l) is.null(cts[[l]]), logical(1))
    cts_present <- cts[has_data]
    inst_levels <- names(cts_present)
    if (length(cts_present) == 0) return(NULL)

    plastic_fams <- c(synthetic_families, semi_synthetic_families)
    keep_fams <- if (sel_fam == "__all_plastics__") {
      if (show_all) NULL
      else if (identical(cat_mode, "synthetic")) synthetic_families
      else plastic_fams
    } else sel_fam

    rows <- do.call(rbind, lapply(inst_levels, function(dev_label) {
      tbl <- cts_present[[dev_label]]
      fams_here <- names(tbl)
      if (!is.null(keep_fams)) fams_here <- fams_here[fams_here %in% keep_fams]
      if (length(fams_here) == 0) return(NULL)
      data.frame(instrument = dev_label, family = fams_here,
                 count = as.integer(tbl[fams_here]), stringsAsFactors = FALSE)
    }))
    if (is.null(rows) || nrow(rows) == 0) return(NULL)

    inst_totals <- tapply(rows$count, rows$instrument, sum)
    rows$value  <- if (rel_mode) round(rows$count / inst_totals[rows$instrument] * 100, 1)
                   else rows$count
    fam_order <- if (sel_fam == "__all_plastics__")
      intersect(if (!is.null(keep_fams)) keep_fams else sort(unique(rows$family)), unique(rows$family))
    else sel_fam
    rows$instrument <- factor(rows$instrument, levels = inst_levels)
    rows$family     <- factor(rows$family, levels = fam_order)

    y_label   <- if (rel_mode) "Share (%)" else "Particle Count"
    title_str <- paste0(if (sel_fam == "__all_plastics__") "All Plastics" else sel_fam,
                        " across instruments", if (use_filt) " (filtered)" else "")

    if (sel_fam == "__all_plastics__") {
      fam_colors <- .fam_palette[levels(rows$family)]
      fam_colors[is.na(fam_colors)] <- "#cccccc"
      names(fam_colors) <- levels(rows$family)
      ggplot2::ggplot(rows, ggplot2::aes(x = instrument, y = value, fill = family)) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::scale_fill_manual(values = fam_colors, name = "Family") +
        ggplot2::labs(x = NULL, y = y_label, title = title_str) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    } else {
      ggplot2::ggplot(rows, ggplot2::aes(x = instrument, y = value, fill = instrument)) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::geom_text(ggplot2::aes(label = value), vjust = -0.3, size = 4.5) +
        ggplot2::scale_fill_manual(values = device_colors[levels(rows$instrument)], guide = "none") +
        ggplot2::scale_y_continuous(limits = c(0, max(rows$value, 1) * 1.5)) +
        ggplot2::labs(x = NULL, y = y_label, title = title_str) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    }
  })

  # Family x device counts. The numbers come from report_plastics_table() in
  # global.R, which the PDF report also uses — the HTML below is presentation
  # only, so the tab and the report can never disagree.
  summary_plastics_df <- reactive({
    d <- summary_dfs()
    report_plastics_table(lapply(SUMMARY_DEVICES, function(key) d[[key]]))
  })

  output$summary_plastics_wide <- renderUI({
    tbl <- summary_plastics_df()
    if (nrow(tbl) == 0)
      return(tags$p(class = "text-muted", "No classified materials found."))

    dev_names <- setdiff(names(tbl), c("Family", "Category"))
    header <- tags$tr(tags$th("Family"), tags$th("Category"),
                      lapply(dev_names, tags$th))
    is_total  <- tbl$Family == "Total"
    body_rows <- lapply(which(!is_total), function(i) {
      tags$tr(tags$td(tags$b(tbl$Family[i])), tags$td(tbl$Category[i]),
              lapply(dev_names, function(dv) tags$td(as.character(tbl[[dv]][i]))))
    })
    total_row <- if (any(is_total)) {
      i <- which(is_total)[1]
      tags$tr(style = "border-top: 2px solid #888;",
              tags$td(tags$b("Total")), tags$td(""),
              lapply(dev_names, function(dv) tags$td(tags$b(as.character(tbl[[dv]][i])))))
    } else NULL
    tags$table(class = "hover-tbl", header, body_rows, total_row)
  }) |> bindCache(summary_plastics_df())

  # Helper: plot size distribution for one instrument
  # plot_size_distribution() now lives in global.R -- shared with the pipeline report

  output$size_hist_ftir <- renderPlot({
    df <- summary_dfs()$ftir
    if (is.null(df) || nrow(df) == 0) {
      return(ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No FTIR data"),
                                   size = 4, colour = "grey50") + theme_void())
    }
    plot_size_distribution(df, "FTIR (PerkinElmer)")
  }) |> bindCache(summary_dfs()$ftir)

  output$size_hist_raman <- renderPlot({
    df <- summary_dfs()$raman
    if (is.null(df) || nrow(df) == 0) {
      return(ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No Raman data"),
                                   size = 4, colour = "grey50") + theme_void())
    }
    plot_size_distribution(df, "Raman", color_matched = "#1f77b4")
  }) |> bindCache(summary_dfs()$raman)

  output$size_hist_ldir <- renderPlot({
    df <- summary_dfs()$ldir
    if (is.null(df) || nrow(df) == 0) {
      return(ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No LDIR data"),
                                   size = 4, colour = "grey50") + theme_void())
    }
    plot_size_distribution(df, "LDIR", color_matched = "#ff7f0e")
  }) |> bindCache(summary_dfs()$ldir)

  # Size statistics table
  # Per-instrument size statistics. Values come from report_size_stats_table()
  # in global.R, shared with the PDF report.
  summary_size_stats_df <- reactive(report_size_stats_table(summary_dfs()))

  output$size_stats_table <- renderUI({
    tbl <- summary_size_stats_df()
    if (nrow(tbl) == 0)
      return(tags$p(class = "text-muted", "No instrument data available"))
    cols <- names(tbl)
    rows <- lapply(seq_len(nrow(tbl)), function(i)
      tags$tr(lapply(seq_along(cols), function(j) {
        v <- as.character(tbl[[j]][i])
        if (j == 1L) tags$td(tags$b(v)) else tags$td(v)
      })))
    tags$table(class = "table table-condensed",
      tags$thead(tags$tr(lapply(cols, tags$th))),
      tags$tbody(rows))
  }) |> bindCache(summary_size_stats_df())

  # ==================================================================
  # PDF REPORT
  # ==================================================================
  # Snapshots the Summary tab and every available instrument view into one
  # multi-page PDF. Every figure is the SAME ggplot object the app is showing
  # (the *_plot_obj() reactives), and every table comes from the same builder
  # as the on-screen HTML, so the report cannot drift from the viewer.

  # Human-readable description of a filter state, for the figure captions.
  # An empty checkbox group reads back as character(0), which %||% does NOT
  # catch (it only tests is.null) — the caption then read "shown: ." instead of
  # saying nothing was selected.
  .or_none <- function(x) if (length(x) == 0) "none" else paste(x, collapse = ", ")

  # How the view treated match status. Worth stating explicitly: the viewers
  # default to "show all detected", so a report page shows every particle the
  # instrument found, not only those matched to another instrument.
  .report_match_desc <- function(prefix) {
    if (isTRUE(input[[paste0(prefix, "_show_all_detected")]]))
      return("all detected (match status ignored)")
    paste0("match status: ", .or_none(input[[paste0(prefix, "_match_filter")]]))
  }

  .fmt_range <- function(rng, unit = "", digits = 2) {
    if (is.null(rng) || length(rng) != 2 || any(!is.finite(rng))) return(NA_character_)
    paste0(round(rng[1], digits), "-", round(rng[2], digits), unit)
  }

  .report_inst_caption <- function(prefix, quality_rng, size_rng, mats,
                                   extra = character(0)) {
    bits <- c(
      if (!is.na(.fmt_range(quality_rng))) paste0("quality ", .fmt_range(quality_rng)),
      if (!is.na(.fmt_range(size_rng, " \u00b5m", 0)))
        paste0("size ", .fmt_range(size_rng, " \u00b5m", 0)),
      if (length(mats) && !("All" %in% mats))
        paste0("materials: ", .or_none(mats)) else "materials: all",
      extra)
    paste0(prefix, if (length(bits)) paste0(" | ", paste(bits, collapse = " | ")) else "")
  }

  # The pages, in order. Built inside the download handler's reactive context.
  report_pages <- function() {
    m        <- active_manifest()
    scope    <- if (summary_filtered()) "each tab's active filters"
                else "all particles in the run"
    run_lbl  <- if (!is.null(uploaded_data())) "uploaded data"
                else basename(selected_run_dir() %||% "unknown")
    analyst  <- trimws(input$report_analyst_name %||% "")
    run_id   <- m$run_id %||% "unknown"
    run_id_line <- if (nzchar(analyst))
      paste0(run_id, "  \u2014  Analyst: ", analyst)
    else run_id
    dfs     <- summary_dfs()
    n_of    <- function(k) { d <- dfs[[k]]; if (is.null(d)) 0L else nrow(d) }

    pages <- list()
    # --- 1. Title / provenance -------------------------------------------
    pages <- c(pages, list(report_text_page(
      "Multi-Instrument Particle Matching \u2014 Report",
      c(paste0("Run                 : ", run_lbl),
        paste0("Run ID              : ", run_id_line),
        paste0("Generated           : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        paste0("Particle scope      : ", scope),
        "",
        "Particles included in this report:",
        paste0("  FTIR (PerkinElmer): ", n_of("ftir")),
        paste0("  FTIR (Bruker)     : ", n_of("ftir_bruker")),
        paste0("  Raman             : ", n_of("raman")),
        paste0("  LDIR              : ", n_of("ldir")),
        "",
        # Source files for this run, from the manifest. Same helper the
        # pipeline report uses, so both list them identically. An uploaded
        # session has no run manifest, so this reports that rather than
        # inventing filenames.
        if (!is.null(uploaded_data()))
          c("Input files:", "  (data uploaded in the viewer; no run manifest)")
        else report_input_file_lines(m),
        "",
        "Every figure reproduces the corresponding viewer tab exactly as it was",
        "displayed when this report was generated. Each figure caption records",
        "the filters that produced it."),
      subtitle = "Generated from the Shiny viewer's Summary tab and instrument views")))

    # --- 2. Summary tab content ------------------------------------------
    bar_title <- paste0("Material Comparison Across Instruments \u2014 ",
                        if (identical(input$summary_material_select, "__all_plastics__"))
                          "All Plastics (stacked)" else (input$summary_material_select %||% ""))
    bar_cap   <- paste0(
      "Scope: ", scope, ". Materials: ",
      if (identical(input$summary_bar_materials_mode, "all")) "all particles"
      else if (identical(input$summary_bar_category_mode %||% "both", "synthetic"))
        "synthetic plastics only"
      else "synthetic + semi-synthetic plastics", ".")

    # The report carries BOTH the absolute and the relative view, regardless of
    # which one the Summary tab's toggle is showing -- the tab is a live
    # exploration, the report is a record. Every other setting (material
    # selection, materials/category mode, scope) still comes from the tab.
    # Built with the shared builder rather than the tab's own reactive, so this
    # page is identical to the pipeline report's.
    .bar_counts <- active_material_counts()
    .bar_sel    <- input$summary_material_select %||% "__all_plastics__"
    .bar_all    <- identical(input$summary_bar_materials_mode, "all")
    .bar_cat    <- input$summary_bar_category_mode %||% "both"

    pages <- c(pages, list(report_full_figure_page(
      build_material_barplot_gg(.bar_counts, .bar_sel, rel_mode = FALSE,
                                show_all = .bar_all, cat_mode = .bar_cat),
      paste0(bar_title, " \u2014 absolute"),
      paste0(bar_cap, " Absolute particle counts."))))
    pages <- c(pages, list(report_full_figure_page(
      build_material_barplot_gg(.bar_counts, .bar_sel, rel_mode = TRUE,
                                show_all = .bar_all, cat_mode = .bar_cat),
      paste0(bar_title, " \u2014 relative"),
      paste0(bar_cap, " Share of each instrument's own total (%), so ",
             "composition is comparable across instruments that found ",
             "different numbers of particles."))))

    pages <- c(pages, list(report_table_page(
      summary_plastics_df(), "Plastics by Instrument",
      paste0("Material family counts per device. Scope: ", scope,
             ". 'Unknown' families are excluded; the Total row sums the ",
             "families shown."))))

    size_plots <- list()
    if (n_of("ftir") > 0)
      size_plots <- c(size_plots, list(plot_size_distribution(dfs$ftir, "FTIR (PerkinElmer)")))
    if (n_of("raman") > 0)
      size_plots <- c(size_plots, list(plot_size_distribution(dfs$raman, "Raman",
                                                             color_matched = "#1f77b4")))
    if (n_of("ldir") > 0)
      size_plots <- c(size_plots, list(plot_size_distribution(dfs$ldir, "LDIR",
                                                             color_matched = "#ff7f0e")))
    pages <- c(pages, list(report_grid_page(
      size_plots, "Size Distribution by Instrument",
      paste0("Feret Max (um) per instrument. Solid bars: matched. Outline bars: ",
             "unmatched. Scope: ", scope, "."),
      ncol = min(3, max(1, length(size_plots))))))

    pages <- c(pages, list(report_table_page(
      summary_size_stats_df(), "Size Statistics",
      paste0("Feret Max summary statistics per instrument. Scope: ", scope, "."))))

    # --- 3. Instrument views ---------------------------------------------
    # Only instruments that actually have data get a page; the plot reactives
    # are the very objects the tabs render.
    coord_lbl <- function(id) if (identical(input[[id]], "aligned"))
      "Raman-aligned frame" else "native frame"

    if (!is.null(ftir_df_full()) && nrow(ftir_df_full()) > 0)
      pages <- c(pages, list(report_figure_page(
        ftir_plot_obj(), "FTIR (PerkinElmer) \u2014 particles over instrument image",
        .report_inst_caption(coord_lbl("ftir_coord_mode"),
                             input$ftir_quality_range, input$ftir_size_range,
                             input$ftir_material_filter,
                             extra = .report_match_desc("ftir")))))

    if (!is.null(ftir_bruker_df_full()) && nrow(ftir_bruker_df_full()) > 0)
      pages <- c(pages, list(report_figure_page(
        ftir_bruker_plot_obj(), "FTIR (Bruker) \u2014 particles over instrument image",
        .report_inst_caption(coord_lbl("ftir_bruker_coord_mode"),
                             input$ftir_bruker_quality_range,
                             input$ftir_bruker_size_range,
                             input$ftir_bruker_material_filter,
                             extra = .report_match_desc("ftir_bruker")))))

    if (!is.null(raman_df_full()) && nrow(raman_df_full()) > 0)
      pages <- c(pages, list(report_figure_page(
        raman_plot_obj(), "Raman \u2014 particles over instrument image",
        .report_inst_caption("native frame",
                             input$raman_quality_range, input$raman_size_range,
                             input$raman_material_filter,
                             extra = .report_match_desc("raman")))))

    if (!is.null(ldir_df_full()) && nrow(ldir_df_full()) > 0)
      pages <- c(pages, list(report_figure_page(
        ldir_plot_obj(), "LDIR \u2014 particles over instrument image",
        .report_inst_caption("native frame",
                             input$ldir_quality_range, input$ldir_size_range,
                             input$ldir_material_filter,
                             extra = c(.report_match_desc("ldir"),
                                       paste0("background: ",
                                              input$ldir_bg_image %||% "auto"))))))

    # --- 4. Overlay -------------------------------------------------------
    pages <- c(pages, list(report_figure_page(
      overlay_plot_obj(), "Overlay \u2014 all instruments in the shared Raman frame",
      paste0("Instruments shown: ", .or_none(input$overlay_instruments),
             ". Relationships: ", .or_none(input$overlay_relationships), "."))))

    Filter(Negate(is.null), pages)
  }

  # Helper: derive the reports/ save folder (project root next to shiny_app/).
  # Shiny sets the working directory to the app folder, so the project root is
  # normally its parent; when the app is run from the root itself, use that.
  # normalizePath() needs mustWork = FALSE here - reports/ usually does not
  # exist yet, and on Windows the default warns and returns the input unchanged.
  .reports_dir <- function() {
    parent <- dirname(getwd())
    root <- if (dir.exists(file.path(parent, "shiny_app"))) parent else getwd()
    normalizePath(file.path(root, "reports"), mustWork = FALSE)
  }

  # Archive both report files without depending on an external zip binary.
  # utils::zip() shells out to Sys.getenv("R_ZIPCMD", "zip"); on Windows that
  # is unset and there is no zip.exe unless Rtools happens to be on PATH, so
  # the archive was never created and the browser got an error page instead of
  # a file. Prefer the zip package (pure C, no external process), fall back to
  # utils::zip, and if neither produces a usable archive serve the HTML report
  # on its own - it already carries every page plus the interactive chart, and
  # both files have been written to reports/ by this point regardless.
  # Can this machine produce a ZIP at all? Decided up front, because
  # downloadHandler evaluates filename() before content() and the extension
  # has to match what we are actually going to serve.
  .can_zip <- function() {
    if (requireNamespace("zip", quietly = TRUE)) return(TRUE)
    cmd <- Sys.getenv("R_ZIPCMD", "zip")
    nzchar(Sys.which(cmd))
  }

  .zip_report <- function(dest, files, html_fallback) {
    ok <- FALSE
    if (requireNamespace("zip", quietly = TRUE)) {
      ok <- tryCatch({
        zip::zip(dest, files = basename(files), root = dirname(files[1]))
        file.exists(dest) && file.info(dest)$size > 0
      }, error = function(e) FALSE)
    }
    if (!ok) {
      ok <- tryCatch({
        utils::zip(dest, files = files, flags = "-j")
        file.exists(dest) && file.info(dest)$size > 0
      }, error = function(e) FALSE, warning = function(w) FALSE)
    }
    if (!ok) {
      # Serve the HTML rather than nothing. It carries every page plus the
      # interactive chart, and both files are already in reports/.
      file.copy(html_fallback, dest, overwrite = TRUE)
      showNotification(
        HTML(paste0(
          "<b>Could not build the ZIP archive.</b><br/>",
          "Downloaded the interactive HTML report instead - it contains ",
          "every page. Both the PDF and the HTML were also saved to the ",
          "<code>reports/</code> folder.<br/>",
          "To get the ZIP, install the zip package: ",
          "<code>install.packages(\"zip\")</code>")),
        type = "warning", duration = 15)
    }
    ok
  }

  .save_to_reports <- function(src, filename) {
    rdir <- .reports_dir()
    dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
    dest <- file.path(rdir, filename)
    tryCatch(file.copy(src, dest, overwrite = TRUE),
             error = function(e) message("[report] auto-save failed: ", e$message))
  }

  output$download_report <- downloadHandler(
    filename = function() {
      run <- if (!is.null(uploaded_data())) "uploaded"
             else basename(selected_run_dir() %||% "run")
      paste0("particle_report_", run, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"),
             if (.can_zip()) ".zip" else ".html")
    },
    contentType = if (.can_zip()) "application/zip" else "text/html",
    content = function(file) {
      stem <- isolate({
        run <- if (!is.null(uploaded_data())) "uploaded"
               else basename(selected_run_dir() %||% "run")
        paste0("particle_report_", run, "_",
               format(Sys.time(), "%Y%m%d_%H%M%S"))
      })

      withProgress(message = "Building report", value = 0, {
        incProgress(0.10, detail = "Collecting figures and tables")
        pages <- tryCatch(report_pages(), error = function(e) {
          list(report_text_page("Report failed",
               c("The report could not be assembled.", "",
                 paste0("Error: ", conditionMessage(e)))))
        })

        incProgress(0.10, detail = "Building interactive chart")
        d       <- summary_dfs()
        dcounts <- lapply(
          setNames(nm = names(SUMMARY_DEVICES)),
          function(lbl) {
            x <- d[[SUMMARY_DEVICES[[lbl]]]]
            if (is.null(x) || nrow(x) == 0 || !"material" %in% names(x)) return(NULL)
            table(classify_family_vec(x$material))
          }
        )
        .mk_fig <- function(rel) tryCatch(
          build_plotly_barplot(dcounts, rel_mode = rel),
          error = function(e) {
            message("[report] plotly build failed: ", e$message); NULL }
        )
        # Interactive twins of the two static barplot pages (2 = absolute,
        # 3 = relative), each placed right after its static counterpart.
        plotly_specs <- list(
          list(fig = .mk_fig(FALSE), after = 2L,
               title = "Material Comparison \u2014 Interactive (absolute counts)",
               caption = paste0("Hover over bars for exact counts. ",
                                "Use the legend to show/hide families.")),
          list(fig = .mk_fig(TRUE), after = 3L,
               title = "Material Comparison \u2014 Interactive (relative share)",
               caption = paste0("Share of each instrument's own total (%). ",
                                "Hover for the underlying count.")))

        tmp_dir  <- tempfile()
        dir.create(tmp_dir)
        pdf_name <- paste0(stem, ".pdf")
        html_name <- paste0(stem, ".html")
        pdf_path  <- file.path(tmp_dir, pdf_name)
        html_path <- file.path(tmp_dir, html_name)

        incProgress(0.35, detail = paste0("Writing PDF (", length(pages), " pages)"))
        write_report_pdf(pages, pdf_path)

        incProgress(0.35, detail = "Rendering HTML")
        write_report_html(pages, html_path, plotly_figs = plotly_specs)

        # Save both files to reports/ BEFORE attempting the archive, so the
        # operator still gets them even if the download degrades below.
        .save_to_reports(pdf_path,  pdf_name)
        .save_to_reports(html_path, html_name)

        incProgress(0.05, detail = "Packaging ZIP")
        .zip_report(file, c(pdf_path, html_path), html_path)
        incProgress(0.05, detail = "Done")
      })
    }
  )

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
    !is.null(d$unmatched_ftir_bruker) ||
    !is.null(d$matched_ftir_bruker)
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
  overlay_ldir_gate_d    <- debounce(reactive(input$overlay_ldir_dist_gate), 300)

  # Effective LDIR<->Raman acceptance gate (µm): the live slider value, falling
  # back to the run's recorded default until the slider is initialised.
  eff_ldir_gate <- reactive({
    v <- overlay_ldir_gate_d()
    if (is.null(v) || !is.finite(v) || v <= 0) {
      g <- run_data()$ldir_match_gate_um
      if (!is.null(g) && is.finite(g)) g else 250
    } else v
  })

  # run_data() re-gated against the live slider. Only the LDIR<->Raman
  # classification changes; all other frames pass through untouched. LDIR
  # consumers read this so the gate retunes matched/unmatched in real time.
  run_data_gated <- reactive({ regate_ldir(run_data(), eff_ldir_gate()) })
  overlay_ftir_bruker_quality_d <- debounce(reactive(input$overlay_ftir_bruker_quality), 300)
  overlay_ftir_bruker_size_d    <- debounce(reactive(input$overlay_ftir_bruker_size), 300)

  # --- Overlay per-instrument filtering -------------------------------------
  # The quality/size sliders used to reach only overlay_matched(), which drives
  # the match lines and the tables. The plotted points came straight from
  # *_df_full(), and .add_unmatched() applies only the material filter -- so
  # moving a quality or size slider visibly did nothing to the unmatched
  # points, and the summary text counted particles the plot was not drawing.
  # Both now read overlay_dfs(), so they cannot disagree.
  .ov_filter <- function(df, q, sz) {
    if (is.null(df) || nrow(df) == 0) return(df)
    keep <- rep(TRUE, nrow(df))
    if (!is.null(q) && "quality" %in% names(df)) {
      v <- suppressWarnings(as.numeric(df$quality))
      keep <- keep & !is.na(v) & v >= q[1] & v <= q[2]
    }
    if (!is.null(sz) && "feret_max" %in% names(df)) {
      v <- suppressWarnings(as.numeric(df$feret_max))
      keep <- keep & !is.na(v) & v >= sz[1] & v <= sz[2]
    }
    df[keep, , drop = FALSE]
  }

  # Unfiltered frames, for the VIEWPORT only: like the instrument tabs, the
  # overlay frames on all particles so the view does not jump while a slider
  # is dragged.
  overlay_dfs_full <- reactive({
    list(ftir = ftir_df_full(), raman = raman_df_full(),
         ldir = ldir_df_full(), ftir_bruker = ftir_bruker_df_full())
  })

  overlay_dfs <- reactive({
    d <- overlay_dfs_full()
    list(
      ftir        = .ov_filter(d$ftir,        overlay_ftir_quality_d(),
                               overlay_ftir_size_d()),
      raman       = .ov_filter(d$raman,       overlay_raman_quality_d(),
                               overlay_raman_size_d()),
      ldir        = .ov_filter(d$ldir,        overlay_ldir_quality_d(),
                               overlay_ldir_size_d()),
      ftir_bruker = .ov_filter(d$ftir_bruker, overlay_ftir_bruker_quality_d(),
                               overlay_ftir_bruker_size_d()))
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
  ftir_raw_image        <- reactiveVal(NULL)   # FTIR "Average Abs" image
  ftir_bruker_raw_image <- reactiveVal(NULL)   # FTIR (Bruker) background image
  raman_image           <- reactiveVal(NULL)   # Raman microscope image (Raman tab + Overlay tab)
  raman_image_path      <- reactiveVal(NULL)   # File path of the Raman image (for TIFF metadata extraction)
  ldir_raw_image        <- reactiveVal(NULL)   # LDIR particle map image

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
    ox <- if (isTRUE(is.finite(input$ftir_img_offset_x))) input$ftir_img_offset_x else 0
    oy <- if (isTRUE(is.finite(input$ftir_img_offset_y))) input$ftir_img_offset_y else 0
    list(raster = raw,
         xmin = min(x_vals) + ox, xmax = max(x_vals) + ox,
         ymin = min(y_vals) + oy, ymax = max(y_vals) + oy)
  })

  # FTIR (Bruker) tab: same particle-extent placement as the PerkinElmer tab.
  ftir_bruker_native_image_info <- reactive({
    raw <- ftir_bruker_raw_image()
    if (is.null(raw)) return(NULL)
    fb_d <- ftir_bruker_df_full()
    if (is.null(fb_d) || nrow(fb_d) == 0) return(NULL)
    x_vals <- fb_d$x_orig[is.finite(fb_d$x_orig)]
    y_vals <- fb_d$y_orig[is.finite(fb_d$y_orig)]
    if (length(x_vals) == 0) return(NULL)
    ox <- if (isTRUE(is.finite(input$ftir_bruker_img_offset_x))) input$ftir_bruker_img_offset_x else 0
    oy <- if (isTRUE(is.finite(input$ftir_bruker_img_offset_y))) input$ftir_bruker_img_offset_y else 0
    list(raster = raw,
         xmin = min(x_vals) + ox, xmax = max(x_vals) + ox,
         ymin = min(y_vals) + oy, ymax = max(y_vals) + oy)
  })

  # Raman tab: image placement with 3-tier priority cascade:
  #   1. Physical extent (WITec center + width/height in µm) — resize-invariant
  #   2. Known scale (raman_um_per_px or auto-detected from TIFF) — centroid-centred
  #   3. Fallback: aspect-ratio-preserving bounds via compute_image_bounds()
  # Cross-check the Raman image placement against the data itself.
  #
  # The WITec extent in the config (raman_image_width_um / _height_um /
  # _center_x_um / _center_y_um) is PER-DATASET. Carrying a previous scan's
  # values into a new run draws the micrograph at the wrong size and the
  # particles land nowhere near their blobs. raman_image_extent_from_config()
  # only rejects an extent that fails to CONTAIN the particles, so a stale
  # extent that is merely too large sails through and misplaces silently.
  #
  # The particle areas give an independent estimate of µm/px, so a gross
  # mismatch can be reported instead of drawn without comment.
  raman_scale_warning <- reactive({
    info <- raman_native_image_info(); raw <- raman_image(); df <- raman_df_full()
    if (is.null(info) || is.null(raw) || is.null(df) ||
        !"area_um2" %in% names(df) || nrow(df) == 0) return(NULL)
    placed <- (info$xmax - info$xmin) / ncol(raw)      # µm/px as actually drawn
    est <- image_scale_from_particle_area(raw, df$area_um2)
    if (is.null(est) || !is.finite(placed) || placed <= 0) return(NULL)
    ratio <- placed / est
    if (ratio <= 2 && ratio >= 0.5) return(NULL)
    sprintf(paste0("⚠ Image scale looks wrong: drawn at %.2f µm/px, but the ",
                   "particle areas imply ~%.2f µm/px (%.1f× off). ",
                   "Check the per-dataset raman_image_* values in 00_config.R."),
            placed, est, if (ratio > 1) ratio else 1 / ratio)
  }) |> bindCache(selected_run_dir(), is.null(uploaded_data()),
                  img_key(raman_native_image_info()))

  raman_native_image_info <- reactive({
    raw <- raman_image()
    if (is.null(raw)) return(NULL)
    raman_df <- raman_df_full()
    ox <- if (isTRUE(is.finite(input$raman_img_offset_x))) input$raman_img_offset_x else 0
    oy <- if (isTRUE(is.finite(input$raman_img_offset_y))) input$raman_img_offset_y else 0

    x_vals <- if (!is.null(raman_df)) raman_df$x_orig[is.finite(raman_df$x_orig)] else numeric(0)
    y_vals <- if (!is.null(raman_df)) raman_df$y_orig[is.finite(raman_df$y_orig)] else numeric(0)

    cfg <- tryCatch(active_manifest()$config_snapshot, error = function(e) list())
    h_px <- nrow(raw); w_px <- ncol(raw)

    # --- Priority 1: Physical extent from WITec metadata (resize-invariant) ---
    # raman_image_extent_from_config() (global.R) turns the Particle Scout
    # panel values into stage-frame bounds, resolving WITec's Y-down video
    # frame vs the particle export's Y-up stage frame by scoring both
    # interpretations against the particles.  NULL means the configured
    # extent fits under neither convention (stale per-dataset values) —
    # fall through to Priority 2/3 instead of drawing the image wrong.
    ext <- raman_image_extent_from_config(
      cfg,
      if (!is.null(raman_df)) raman_df$x_orig else numeric(0),
      if (!is.null(raman_df)) raman_df$y_orig else numeric(0))
    if (!is.null(ext)) {
      message("[Particle Viewer] Raman image placed from WITec extent (",
              if (isTRUE(ext$y_negated))
                "Center Y negated: video frame -> stage frame"
              else "Center Y as reported",
              "; ", round(ext$frac_inside * 100), "% of particles inside).")
      return(list(raster = raw,
                  xmin = ext$xmin + ox, xmax = ext$xmax + ox,
                  ymin = ext$ymin + oy, ymax = ext$ymax + oy))
    }
    if (!is.null(cfg$raman_image_width_um) &&
        !is.null(cfg$raman_image_center_x_um)) {
      message("[Particle Viewer] WITec raman_image_* extent does not contain ",
              "this run's particles under either Y convention — the values ",
              "likely belong to a different dataset. Falling back to ",
              "heuristic placement; update Width/Height/Center X/Center Y ",
              "from WITec's Particle Scout for this scan.")
    }

    # --- Priority 2: Known scale (config or auto-detected from TIFF DPI) ---
    um_per_px <- cfg$raman_um_per_px
    if (is.null(um_per_px) || !is.numeric(um_per_px) || um_per_px <= 0) {
      um_per_px <- extract_tiff_um_per_px(raman_image_path())
    }

    if (!is.null(um_per_px) && length(x_vals) > 0) {
      # um_per_px refers to the ORIGINAL file; if the raster was downsized
      # after upload (downsample_raster), rescale it to the reduced raster.
      orig_w <- attr(raw, "orig_width_px")
      if (!is.null(orig_w) && is.numeric(orig_w) && orig_w > 0)
        um_per_px <- um_per_px * orig_w / w_px
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
    ox <- if (isTRUE(is.finite(input$overlay_img_offset_x))) input$overlay_img_offset_x else 0
    oy <- if (isTRUE(is.finite(input$overlay_img_offset_y))) input$overlay_img_offset_y else 0

    # Exact placement when WITec metadata is available: overlay coordinates
    # are Raman stage coordinates minus the Raman centroid (pure translation,
    # normalize_coordinates() applies no scale/rotation to Raman), so the
    # stage-frame extent maps into overlay space by subtracting that same
    # centroid — recovered per-row as (x_orig - x), constant across particles.
    cfg <- tryCatch(active_manifest()$config_snapshot, error = function(e) list())
    ext <- raman_image_extent_from_config(cfg, raman_d$x_orig, raman_d$y_orig)
    fin <- is.finite(raman_d$x) & is.finite(raman_d$x_orig) &
           is.finite(raman_d$y) & is.finite(raman_d$y_orig)
    if (!is.null(ext) && any(fin)) {
      sx <- mean(raman_d$x_orig[fin] - raman_d$x[fin])
      sy <- mean(raman_d$y_orig[fin] - raman_d$y[fin])
      return(list(raster = raw,
                  xmin = ext$xmin - sx + ox, xmax = ext$xmax - sx + ox,
                  ymin = ext$ymin - sy + oy, ymax = ext$ymax - sy + oy))
    }

    # Fallback: aspect-ratio-preserving fit to the Raman particle extent
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
    ox <- if (isTRUE(is.finite(input$ldir_img_offset_x))) input$ldir_img_offset_x else 0
    oy <- if (isTRUE(is.finite(input$ldir_img_offset_y))) input$ldir_img_offset_y else 0

    m <- tryCatch(active_manifest(), error = function(e) NULL)

    # Circle-based bounds: derive exact image placement from scan-circle calibration.
    # cx_px/cy_px are the circle centre in pixel space; scale_um_per_px converts
    # pixels to µm.  The image spans from -cx_px*scale to (w-cx_px)*scale in x
    # and from (cy_px-h)*scale to cy_px*scale in y (y upward, row 0 = ymax).
    ci <- if (!is.null(m)) m$ldir_circle else NULL
    # Require every placement field to be a finite number: a partial/old
    # manifest can carry JSON null (-> R NULL) for cx_px/cy_px/width/height,
    # and `-NULL` errors ("invalid argument to unary operator"). Fall through
    # to the symmetric scan-diameter bounds below when anything is missing.
    .fin1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    if (!is.null(ci) && .fin1(ci$scale_um_per_px) && ci$scale_um_per_px > 0 &&
        .fin1(ci$cx_px) && .fin1(ci$cy_px) &&
        .fin1(ci$image_width_px) && .fin1(ci$image_height_px)) {
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
      ftir_bruker_raw_image(NULL)
      raman_image(NULL)
      raman_image_path(NULL)
      ldir_raw_image(NULL)
      return()
    }
    m   <- load_run_manifest(run_dir)
    img <- get_run_image_paths(m, run_dir)

    # Load run-directory images at full resolution. (An earlier downsample here
    # for render speed block-averaged the raster, washing out crisp instrument
    # images — reverted; correctness of the background wins over the render
    # speed-up. Uploaded images are still downsized in handle_image_upload.)
    load_bg <- function(path) load_image_raster(path)

    if (!is.null(img$ftir)) {
      raw <- load_bg(img$ftir)
      if (!is.null(raw)) ftir_raw_image(raw)
    } else {
      ftir_raw_image(NULL)
    }

    if (!is.null(img$ftir_bruker)) {
      raw <- load_bg(img$ftir_bruker)
      if (!is.null(raw)) ftir_bruker_raw_image(raw)
    } else {
      ftir_bruker_raw_image(NULL)
    }

    if (!is.null(img$raman)) {
      raw <- load_bg(img$raman)
      if (!is.null(raw)) {
        raman_image(raw)
        raman_image_path(img$raman)
      }
    } else {
      raman_image(NULL)
      raman_image_path(NULL)
    }

    if (!is.null(img$ldir)) {
      raw <- load_bg(img$ldir)
      if (!is.null(raw)) ldir_raw_image(raw)
    } else {
      ldir_raw_image(NULL)
    }
  }, ignoreNULL = FALSE)

  # ------------------------------------------------------------------
  # Handle uploaded images
  # ------------------------------------------------------------------
  # Validate + load + downsize one uploaded background image.
  # Returns the raster array, or NULL after notifying the user of the failure
  # (previously a bad read failed silently and multi-file selections crashed
  # the observer before the image was ever stored).
  handle_image_upload <- function(fileinfo, instrument) {
    req(fileinfo)
    if (nrow(fileinfo) > 1) {
      showNotification(
        paste0(nrow(fileinfo), " files selected — using \"", fileinfo$name[1],
               "\" as the ", toupper(instrument),
               " background image (one image per instrument)."),
        type = "warning", duration = 8)
    }
    path    <- fileinfo$datapath[1]
    name    <- fileinfo$name[1]
    size_mb <- round(fileinfo$size[1] / 1024^2, 1)
    message("[Particle Viewer] ", toupper(instrument), " image upload: ", name,
            " (", size_mb, " MB, signature=", sniff_image_type(path), ")")

    raw <- load_image_raster(path)
    if (is.null(raw)) {
      msg <- paste0("Could not read \"", name,
                    "\" as an image (PNG/JPEG/BMP/TIFF/WEBP). ",
                    "TIFF/WEBP and compressed BMPs require the 'magick' ",
                    "package on the server.")
      message("[Particle Viewer] ERROR: ", msg)
      showNotification(msg, type = "error", duration = 10)
      return(NULL)
    }

    dims_in <- dim(raw)
    raw <- downsample_raster(raw, max_dim = BG_IMAGE_MAX_DIM)
    dims_out <- dim(raw)
    if (!identical(dims_in[1:2], dims_out[1:2])) {
      message("[Particle Viewer] ", toupper(instrument), " image downsized: ",
              dims_in[2], "x", dims_in[1], " -> ", dims_out[2], "x", dims_out[1])
      showNotification(
        paste0("\"", name, "\" loaded and downsized from ",
               dims_in[2], "×", dims_in[1], " to ",
               dims_out[2], "×", dims_out[1],
               " px for display (longest edge capped at ",
               BG_IMAGE_MAX_DIM, " px)."),
        type = "message", duration = 8)
    } else {
      showNotification(
        paste0("\"", name, "\" loaded as ", toupper(instrument),
               " background image."),
        type = "message", duration = 5)
    }
    raw
  }

  # Persist an uploaded raster into the active run directory as
  # inputs/<instrument>_image_uploaded.png.  get_run_image_paths() checks this
  # name first, so the upload survives run switches and app restarts instead
  # of living only in this session's memory.  Failures (e.g. read-only deploy
  # dir) are logged to the console but never fatal.
  persist_uploaded_image <- function(raw, instrument) {
    run_dir <- selected_run_dir()
    if (is.null(raw) || is.null(run_dir) || !dir.exists(run_dir))
      return(invisible(NULL))
    target_dir <- file.path(run_dir, "inputs")
    target <- file.path(target_dir, paste0(instrument, "_image_uploaded.png"))
    ok <- tryCatch({
      if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
      png::writePNG(raw, target)
      TRUE
    }, error = function(e) {
      message("[Particle Viewer] WARNING: could not persist ",
              toupper(instrument), " image to ", target, ": ",
              conditionMessage(e))
      FALSE
    })
    if (ok)
      message("[Particle Viewer] Saved uploaded ", toupper(instrument),
              " image to ", target)
    invisible(NULL)
  }

  observeEvent(input$ftir_image_upload, {
    raw <- handle_image_upload(input$ftir_image_upload, "ftir")
    if (!is.null(raw)) {
      ftir_raw_image(raw)
      persist_uploaded_image(raw, "ftir")
    }
  })

  observeEvent(input$ftir_bruker_image_upload, {
    raw <- handle_image_upload(input$ftir_bruker_image_upload, "ftir_bruker")
    if (!is.null(raw)) {
      ftir_bruker_raw_image(raw)
      persist_uploaded_image(raw, "ftir_bruker")
    }
  })

  observeEvent(input$raman_image_upload, {
    raw <- handle_image_upload(input$raman_image_upload, "raman")
    if (!is.null(raw)) {
      raman_image(raw)
      raman_image_path(input$raman_image_upload$datapath[1])
      persist_uploaded_image(raw, "raman")
    }
  })

  observeEvent(input$overlay_image_upload, {
    raw <- handle_image_upload(input$overlay_image_upload, "raman")
    if (!is.null(raw)) {
      raman_image(raw)
      raman_image_path(input$overlay_image_upload$datapath[1])
      persist_uploaded_image(raw, "raman")
    }
  })

  observeEvent(input$ldir_image_upload, {
    raw <- handle_image_upload(input$ldir_image_upload, "ldir")
    if (!is.null(raw)) {
      ldir_raw_image(raw)
      persist_uploaded_image(raw, "ldir")
    }
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
      ftir_mats <- sort(unique(ftir$material_family))
      ftir_ids  <- natural_sort_ids(unique(ftir$particle_id))
      q_range   <- range(ftir$quality, na.rm = TRUE)
      s_max     <- ceiling(max(ftir$feret_max, na.rm = TRUE) / 10) * 10

      # Individual tab
      update_material_choices("ftir_material_filter", ftir_mats)
      updateSelectInput(session, "ftir_highlight_particle",
                        choices = c("None", ftir_ids))
      updateSliderInput(session, "ftir_quality_range",
                        min = 0, max = 1,
                        value = c(0.7, 1))
      updateSliderInput(session, "ftir_size_range", min = 0, max = s_max,
                        value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))

      # Overlay per-instrument
      update_material_choices("overlay_ftir_material", ftir_mats)
      updateSelectizeInput(session, "overlay_ftir_particles",
                           choices = ftir_ids, selected = character(0))
      # Match the single-instrument tabs: fixed scale and fixed default,
      # not the observed data range. Deriving min/max/value from the data
      # reset these to full range on every load, wiping the default.
      updateSliderInput(session, "overlay_ftir_quality",
                        min = 0, max = 1, value = c(0.7, 1))
      updateSliderInput(session, "overlay_ftir_size", min = 0, max = s_max,
                        value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))
    }

    # --- Raman controls (individual tab + overlay) ---
    if (!is.null(raman_d) && nrow(raman_d) > 0) {
      raman <- raman_d
      raman_mats <- sort(unique(raman$material_family))
      raman_ids  <- natural_sort_ids(unique(raman$particle_id))
      q_range    <- range(raman$quality, na.rm = TRUE)
      s_max      <- ceiling(max(raman$feret_max, na.rm = TRUE) / 10) * 10

      # Individual tab
      update_material_choices("raman_material_filter", raman_mats)
      updateSelectInput(session, "raman_highlight_particle",
                        choices = c("None", raman_ids))
      updateSliderInput(session, "raman_quality_range",
                        min = 0, max = 100,
                        value = c(70, 100))
      updateSliderInput(session, "raman_size_range", min = 0, max = s_max,
                        value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))

      # Overlay per-instrument
      update_material_choices("overlay_raman_material", raman_mats)
      updateSelectizeInput(session, "overlay_raman_particles",
                           choices = raman_ids, selected = character(0))
      # Match the single-instrument tabs: fixed scale and fixed default,
      # not the observed data range. Deriving min/max/value from the data
      # reset these to full range on every load, wiping the default.
      updateSliderInput(session, "overlay_raman_quality",
                        min = 0, max = 100, value = c(70, 100))
      updateSliderInput(session, "overlay_raman_size", min = 0, max = s_max,
                        value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))
    }

    # --- LDIR controls (individual tab + overlay) ---
    if (!is.null(ldir_d) && nrow(ldir_d) > 0) {
      ldir <- ldir_d
      ldir_mats <- sort(unique(ldir$material_family))
      ldir_ids  <- natural_sort_ids(unique(ldir$particle_id))
      q_range   <- range(ldir$quality, na.rm = TRUE)
      s_max     <- ceiling(max(ldir$feret_max, na.rm = TRUE) / 10) * 10

      # Individual tab
      update_material_choices("ldir_material_filter", ldir_mats)
      updateSelectInput(session, "ldir_highlight_particle",
                        choices = c("None", ldir_ids))
      if (all(is.finite(q_range))) {
        updateSliderInput(session, "ldir_quality_range",
                          min = 0, max = 1,
                          value = c(0.8, 1))
      }
      if (is.finite(s_max)) {
        updateSliderInput(session, "ldir_size_range", min = 0, max = s_max,
                          value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))
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
      update_material_choices("overlay_ldir_material", ldir_mats)
      updateSelectizeInput(session, "overlay_ldir_particles",
                           choices = ldir_ids, selected = character(0))
      if (all(is.finite(q_range))) {
        # Match the single-instrument tabs: fixed scale and fixed default,
        # not the observed data range. Deriving min/max/value from the data
        # reset these to full range on every load, wiping the default.
        updateSliderInput(session, "overlay_ldir_quality",
                          min = 0, max = 1, value = c(0.8, 1))
      }
      if (is.finite(s_max)) {
        updateSliderInput(session, "overlay_ldir_size", min = 0, max = s_max,
                          value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))
      }
    }

    # --- FTIR Bruker controls (individual tab + overlay) ---
    ftir_bruker_d <- ftir_bruker_df_full()
    if (!is.null(ftir_bruker_d) && nrow(ftir_bruker_d) > 0) {
      fb_mats  <- sort(unique(ftir_bruker_d$material_family))
      fb_ids   <- natural_sort_ids(unique(ftir_bruker_d$particle_id))
      q_range  <- range(ftir_bruker_d$quality, na.rm = TRUE)
      s_max    <- ceiling(max(ftir_bruker_d$feret_max, na.rm = TRUE) / 10) * 10
      update_material_choices("ftir_bruker_material_filter", fb_mats)
      updateSelectInput(session, "ftir_bruker_highlight_particle",
                        choices = c("None", fb_ids))
      if (all(is.finite(q_range))) {
        updateSliderInput(session, "ftir_bruker_quality_range",
                          min = 0, max = 1,
                          value = c(0.7, 1))
      }
      if (is.finite(s_max)) {
        updateSliderInput(session, "ftir_bruker_size_range", min = 0, max = s_max,
                          value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))
      }

      # Overlay per-instrument
      update_material_choices("overlay_ftir_bruker_material", fb_mats)
      updateSelectizeInput(session, "overlay_ftir_bruker_particles",
                           choices = fb_ids, selected = character(0))
      if (all(is.finite(q_range))) {
        # Match the single-instrument tabs: fixed scale and fixed default,
        # not the observed data range. Deriving min/max/value from the data
        # reset these to full range on every load, wiping the default.
        updateSliderInput(session, "overlay_ftir_bruker_quality",
                          min = 0, max = 1, value = c(0.7, 1))
      }
      if (is.finite(s_max)) {
        updateSliderInput(session, "overlay_ftir_bruker_size", min = 0, max = s_max,
                          value = c(min(DEFAULT_MIN_SIZE_UM, s_max), s_max))
      }
    }

    # --- Global overlay controls ---
    if (!is.null(run_data()$matched)) {
      max_dist <- ceiling(max(run_data()$matched$match_distance, na.rm = TRUE))
      updateSliderInput(session, "overlay_dist_range",
                        min = 0, max = max_dist, value = c(0, max_dist))
    }

    # LDIR<->Raman live match gate: seed the slider at the run's recorded gate
    # and open the range up to the largest forced-pair distance so every pair
    # can be admitted if the user drags it all the way up.
    lrm <- run_data()$ldir_raman_matched
    if (!is.null(lrm) && "match_distance" %in% names(lrm) &&
        any(is.finite(lrm$match_distance))) {
      gate0 <- run_data()$ldir_match_gate_um
      if (is.null(gate0) || !is.finite(gate0)) gate0 <- 250
      max_g <- max(ceiling(max(lrm$match_distance, na.rm = TRUE)), gate0)
      updateSliderInput(session, "overlay_ldir_dist_gate",
                        min = 0, max = max_g, value = gate0)
    }
  })

  # Global Feret Max constrains per-instrument size sliders
  observeEvent(input$overlay_size_range, {
    global <- input$overlay_size_range
    for (slider_id in c("overlay_ftir_size", "overlay_raman_size", "overlay_ldir_size",
                        "overlay_ftir_bruker_size")) {
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

  observeEvent(input$overlay_ftir_bruker_apply_pattern, {
    pat <- input$overlay_ftir_bruker_pattern
    df <- ftir_bruker_df_full()
    if (is.null(df) || nrow(df) == 0 || nchar(trimws(pat)) == 0) return()
    matched_ids <- parse_particle_selection(pat, unique(df$particle_id))
    current <- input$overlay_ftir_bruker_particles
    new_sel <- unique(c(current, matched_ids))
    updateSelectizeInput(session, "overlay_ftir_bruker_particles", selected = new_sel)
    updateTextInput(session, "overlay_ftir_bruker_pattern", value = "")
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

  # Effective match-status filter for a viewer: when its "Show all detected"
  # box is ticked, ignore the Match Status checkboxes and keep every status.
  eff_match_filter <- function(prefix) {
    if (isTRUE(input[[paste0(prefix, "_show_all_detected")]]))
      return(c("matched", "unmatched"))
    input[[paste0(prefix, "_match_filter")]]
  }

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
    # mat_filter: NULL = no filter, otherwise the families to keep (possibly
    # none, when every box is unchecked). See mat_keep().
    if (!is.null(mat_filter))
      df <- df[df$material_family %in% mat_filter, ]
    df <- df[df$match_status %in% match_filter, ]
    df
  }

  # ==================================================================
  # Per-material show/hide checkbox groups
  # ==================================================================
  # Every viewer and the overlay lists the material families actually present
  # in its data, all ticked on load, so any of them can be switched off.
  #
  # An empty checkbox group reads back as NULL, which is ambiguous: it is
  # either "not populated yet" (data still loading — must not filter) or "the
  # user unticked everything" (must show nothing). mat_ready distinguishes
  # them: it flips once the group has ever reported a selection.
  MATERIAL_INPUTS <- c(
    "ftir_material_filter", "ftir_bruker_material_filter",
    "raman_material_filter", "ldir_material_filter",
    "overlay_ftir_material", "overlay_raman_material",
    "overlay_ldir_material", "overlay_ftir_bruker_material")

  mat_ready   <- reactiveValues()
  mat_choices <- reactiveValues()   # input_id -> families offered

  # Families to keep for an input: NULL = no filtering at all.
  mat_keep <- function(input_id) {
    sel <- input[[input_id]]
    if (is.null(sel)) {
      if (!isTRUE(mat_ready[[input_id]])) return(NULL)   # not populated yet
      return(character(0))                               # all unticked
    }
    sel
  }

  # Fill a group from the data, keeping the user's ticks where the family still
  # exists; families new to this run arrive ticked, so switching runs never
  # silently hides particles.
  update_material_choices <- function(input_id, families) {
    families <- sort(unique(as.character(families[!is.na(families)])))
    prev_ch  <- mat_choices[[input_id]]
    prev_sel <- isolate(input[[input_id]])
    selected <- if (is.null(prev_ch)) families
                else union(intersect(prev_sel, families), setdiff(families, prev_ch))
    mat_choices[[input_id]] <- families
    updateCheckboxGroupInput(session, input_id, choices = families,
                             selected = selected)
  }

  for (.mid in MATERIAL_INPUTS) local({
    id <- .mid
    # Readiness: ignoreNULL so unticking everything does NOT reset the flag.
    observeEvent(input[[id]], { mat_ready[[id]] <- TRUE }, ignoreNULL = TRUE)
    # Re-supplying choices (rather than selected alone) makes "none" reliable.
    observeEvent(input[[paste0(id, "_all")]], {
      updateCheckboxGroupInput(session, id, choices = mat_choices[[id]],
                               selected = mat_choices[[id]])
    })
    observeEvent(input[[paste0(id, "_none")]], {
      mat_ready[[id]] <- TRUE
      updateCheckboxGroupInput(session, id, choices = mat_choices[[id]],
                               selected = character(0))
    })
  })

  # ==================================================================
  # Helper: add image background to a ggplot
  # ==================================================================
  # add_image_bg() now lives in global.R -- shared with the pipeline report

  # Cheap cache-key identity for a background image: its placement bounds (four
  # numbers) rather than the raster pixels. Used in the plot renderPlot()
  # bindCache() keys below — combined with selected_run_dir() (which changes
  # when the run, and therefore the image pixels, change) this captures image
  # movement (offset sliders, config) without hashing megapixels on every flush.
  img_key <- function(ii) {
    if (is.null(ii)) return("none")
    paste(round(c(ii$xmin, ii$xmax, ii$ymin, ii$ymax), 1), collapse = ",")
  }

  # ==================================================================
  # Helper: ggplot scatter with optional image background
  # ==================================================================
  # Helper: generate axis breaks at every 1000 µm within a range
  # Adaptive axis breaks: choose interval based on zoom level to keep 4-8 breaks visible
  # breaks_adaptive() now lives in global.R -- shared with the pipeline report

  # Positive-width limits for scale_size_continuous. When a view is filtered down
  # to a single particle (or one distinct Feret), the default size domain has
  # zero width, so ggplot rescales from 0/0 -> NaN and the NaN-sized legend key
  # throws "non-finite location/size for viewport", killing the whole plot.
  # Returns NULL when there are no finite values (scale then goes unused).
  # safe_size_limits() now lives in global.R -- shared with the pipeline report

  # Draw the particle_id next to every point (the "Number all particles" toggle).
  # White text with a thin black shadow so it reads over both the light panel and
  # dark instrument images. Plain geom_text (no repel) keeps it fast for the
  # larger instruments; it is an opt-in overlay so occasional overlap is fine.
  # add_particle_labels() now lives in global.R -- shared with the pipeline report

  # make_scatter() now lives in global.R -- shared with the pipeline report

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
                               ftir_bruker = NULL, repro = NULL)

  # Pinned overlay particle: persists across hover events until cleared.
  # Stores a data row (matched or single-instrument) and its source type.
  pinned_overlay <- reactiveVal(NULL)
  pinned_source  <- reactiveVal(NULL)   # "ftir_raman", "ldir_raman", or "single_ftir"/"single_raman"/"single_ldir"

  # ==================================================================
  # Zoom state: NULL means full view, otherwise list(x=c(lo,hi), y=c(lo,hi))
  # ==================================================================
  zoom <- reactiveValues(ftir = NULL, raman = NULL, ldir = NULL, overlay = NULL,
                         ftir_bruker = NULL, repro = NULL)

  # Smallest viewport span we will ever draw (µm). Below this, coord_fixed()
  # with expand = FALSE collapses the panel and grid aborts the whole plot with
  # "non-finite location and/or size for viewport", which reads to the user as
  # the image vanishing.
  .MIN_SPAN <- 1e-6

  # A brush that is really a click (zero/near-zero area) must not become a zoom.
  .brush_bounds <- function(b) {
    if (is.null(b)) return(NULL)
    v <- c(b$xmin, b$xmax, b$ymin, b$ymax)
    if (!all(is.finite(v))) return(NULL)
    if ((b$xmax - b$xmin) < .MIN_SPAN || (b$ymax - b$ymin) < .MIN_SPAN) return(NULL)
    list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax))
  }

  # Last line of defence before coord_fixed(). Limits can still go non-finite
  # without any brush involved — all-NA coordinates make min()/max() return
  # +/-Inf, and unresolved image-placement metadata yields NA extents.
  # sanitize_bounds() now lives in global.R -- shared with the pipeline report

  observeEvent(input$ftir_brush, {
    nb <- .brush_bounds(input$ftir_brush); if (!is.null(nb)) zoom$ftir <- nb
  })
  observeEvent(input$ftir_dblclick, { zoom$ftir <- NULL })

  observeEvent(input$raman_brush, {
    nb <- .brush_bounds(input$raman_brush); if (!is.null(nb)) zoom$raman <- nb
  })
  observeEvent(input$raman_dblclick, { zoom$raman <- NULL })

  observeEvent(input$ldir_brush, {
    nb <- .brush_bounds(input$ldir_brush); if (!is.null(nb)) zoom$ldir <- nb
  })
  observeEvent(input$ldir_dblclick, { zoom$ldir <- NULL })

  observeEvent(input$ftir_bruker_brush, {
    nb <- .brush_bounds(input$ftir_bruker_brush)
    if (!is.null(nb)) zoom$ftir_bruker <- nb
  })
  observeEvent(input$ftir_bruker_dblclick, { zoom$ftir_bruker <- NULL })

  observeEvent(input$overlay_brush, {
    nb <- .brush_bounds(input$overlay_brush); if (!is.null(nb)) zoom$overlay <- nb
  })
  observeEvent(input$overlay_dblclick, { zoom$overlay <- NULL })

  observeEvent(input$repro_brush, {
    nb <- .brush_bounds(input$repro_brush); if (!is.null(nb)) zoom$repro <- nb
  })
  observeEvent(input$repro_dblclick, { zoom$repro <- NULL })

  # Reset-zoom buttons (same effect as double-click)
  observeEvent(input$ftir_reset_zoom,        { zoom$ftir <- NULL })
  observeEvent(input$raman_reset_zoom,       { zoom$raman <- NULL })
  observeEvent(input$ldir_reset_zoom,        { zoom$ldir <- NULL })
  observeEvent(input$ftir_bruker_reset_zoom, { zoom$ftir_bruker <- NULL })
  observeEvent(input$overlay_reset_zoom,     { zoom$overlay <- NULL })
  observeEvent(input$repro_reset_zoom,       { zoom$repro <- NULL })

  # Filter changes invalidate any active zoom (matches other tabs: a coord
  # or filter change can move the data outside the old zoom window).
  observeEvent(input$repro_material,       { zoom$repro <- NULL })
  observeEvent(input$repro_quality_range,  { zoom$repro <- NULL })
  observeEvent(input$repro_only_nonrepro,  { zoom$repro <- NULL })
  observeEvent(input$repro_run_select,     { zoom$repro <- NULL })
  observeEvent(input$repro_run_visibility, { zoom$repro <- NULL })
  observeEvent(input$repro_img_rotation,   { zoom$repro <- NULL })

  observeEvent(input$ftir_coord_mode,        { zoom$ftir        <- NULL })
  observeEvent(input$ftir_view_rotation,     { zoom$ftir        <- NULL })
  observeEvent(input$ftir_view_flip_y,       { zoom$ftir        <- NULL })
  observeEvent(input$ftir_bruker_coord_mode, { zoom$ftir_bruker <- NULL })
  observeEvent(input$ftir_bruker_view_rotation, { zoom$ftir_bruker <- NULL })
  observeEvent(input$ftir_bruker_view_flip_y,   { zoom$ftir_bruker <- NULL })
  observeEvent(input$ldir_coord_mode,        { zoom$ldir        <- NULL })
  observeEvent(input$ldir_view_rotation,     { zoom$ldir        <- NULL })
  observeEvent(input$ldir_view_flip_y,       { zoom$ldir        <- NULL })
  observeEvent(input$ldir_bg_image,          { zoom$ldir        <- NULL })

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
    pid <- .find_nearest(input$ftir_click,
                         .view_click_df(ftir_filtered(), input$ftir_coord_mode,
                                        ftir_view_tf()),
                         zoom$ftir)
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

  # LDIR native display may be view-rotated to match the Raman orientation;
  # the click arrives in rotated plot space, so rotate the lookup coordinates
  # the same way before nearest-particle search.
  .ldir_click_df <- function() {
    df <- ldir_filtered()
    cm <- input$ldir_coord_mode
    if (!is.null(cm) && cm == "native" && !is.null(df) && nrow(df) > 0) {
      vt <- ldir_view_tf()
      if (!.tf_is_identity(vt)) {
        rc <- view_transform_xy(df$x_orig, df$y_orig, vt$deg, vt$flip)
        df$x_orig <- rc$x; df$y_orig <- rc$y
      }
    }
    df
  }

  observeEvent(input$ldir_click, {
    pid <- .find_nearest(input$ldir_click, .ldir_click_df(), zoom$ldir)
    if (!is.null(pid)) {
      cur <- selected_ids$ldir
      selected_ids$ldir <- if (pid %in% cur) cur else c(cur, pid)
    }
  })
  observeEvent(input$ldir_clear_selection, { selected_ids$ldir <- character(0) })

  observeEvent(input$ftir_bruker_click, {
    pid <- .find_nearest(input$ftir_bruker_click,
                         .view_click_df(ftir_bruker_filtered(),
                                        input$ftir_bruker_coord_mode,
                                        ftir_bruker_view_tf()),
                         zoom$ftir_bruker)
    if (!is.null(pid)) {
      cur <- selected_ids$ftir_bruker
      selected_ids$ftir_bruker <- if (pid %in% cur) cur else c(cur, pid)
    }
  })
  observeEvent(input$ftir_bruker_clear_selection, { selected_ids$ftir_bruker <- character(0) })

  # Overlay: select / deselect all instruments
  observeEvent(input$overlay_toggle_all, {
    all_inst <- c("ftir_pe", "ftir_bruker", "raman", "ldir")
    current  <- input$overlay_instruments
    if (length(current) == length(all_inst)) {
      updateCheckboxGroupInput(session, "overlay_instruments", selected = character(0))
    } else {
      updateCheckboxGroupInput(session, "overlay_instruments", selected = all_inst)
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
                      mat_keep("ftir_material_filter"), eff_match_filter("ftir"))
  })

  # ------------------------------------------------------------------
  # Native-view orientation (FTIR tabs)
  # ------------------------------------------------------------------
  # An FTIR export can sit in a different orientation — and a different
  # handedness — from Raman, so the native tab can look nothing like the
  # Raman tab or the overlay. These helpers resolve a display-only
  # rotation + mirror and apply it to the whole native scene (points AND
  # image). "auto" measures it from the two particle clouds; the Mirror
  # checkbox applies on top, so ticking it twice returns to the measured
  # orientation. Aligned mode is already in Raman space and is never
  # transformed.
  # Largest residual (as a fraction of cloud radius) still accepted as a real
  # orientation fit; above this the pairing is not a rigid transform and the
  # cloud-matching fallback is used instead.
  .VIEW_FIT_MAX_RMSE <- 0.35

  .resolve_view_tf <- function(rot_sel, flip_sel, df_full) {
    flip <- isTRUE(flip_sel)
    if (!is.null(rot_sel) && rot_sel != "auto")
      return(list(deg = as.integer(rot_sel), flip = flip))
    if (is.null(df_full) || nrow(df_full) == 0)
      return(list(deg = 0L, flip = flip))

    # Preferred: measure against this instrument's own aligned coordinates,
    # which is the transform the overlay already draws with.
    if (all(c("x", "y", "x_orig", "y_orig") %in% names(df_full))) {
      auto <- auto_view_from_pairs(df_full$x_orig, df_full$y_orig,
                                   df_full$x, df_full$y)
      if (!is.null(auto) && is.finite(auto$rmse) && auto$rmse <= .VIEW_FIT_MAX_RMSE)
        return(list(deg = auto$deg, flip = xor(isTRUE(auto$flip), flip)))
    }
    # Fallback: match the two particle clouds (run never aligned to Raman).
    rd <- raman_df_full()
    if (is.null(rd) || nrow(rd) == 0) return(list(deg = 0L, flip = flip))
    auto <- auto_view_dihedral(df_full$x_orig, df_full$y_orig, rd$x_orig, rd$y_orig)
    list(deg = auto$deg, flip = xor(isTRUE(auto$flip), flip))
  }

  ftir_view_tf <- reactive(
    .resolve_view_tf(input$ftir_view_rotation, input$ftir_view_flip_y,
                     ftir_df_full()))
  ftir_bruker_view_tf <- reactive(
    .resolve_view_tf(input$ftir_bruker_view_rotation, input$ftir_bruker_view_flip_y,
                     ftir_bruker_df_full()))

  .tf_is_identity <- function(vt) vt$deg == 0L && !isTRUE(vt$flip)

  # Transform a points df already carrying display coords in x/y.
  .apply_view_tf <- function(df, vt) {
    if (is.null(df) || nrow(df) == 0 || .tf_is_identity(vt)) return(df)
    rc <- view_transform_xy(df$x, df$y, vt$deg, vt$flip)
    df$x <- rc$x; df$y <- rc$y
    df
  }

  .apply_view_tf_img <- function(img, vt) {
    if (is.null(img) || .tf_is_identity(vt)) return(img)
    ext <- view_transform_extent(img, vt$deg, vt$flip)
    list(raster = view_transform_raster(img$raster, vt$deg, vt$flip),
         xmin = ext$xmin, xmax = ext$xmax, ymin = ext$ymin, ymax = ext$ymax)
  }

  # Spell the direction out rather than printing a signed angle. The sign was
  # ambiguous: +90 is counter-clockwise by the maths convention but clockwise
  # by the screen/image convention (y down), so the title read as inverted.
  # These strings match the dropdown labels exactly.
  .view_tf_suffix <- function(vt) {
    if (.tf_is_identity(vt)) return("")
    rot <- switch(as.character(vt$deg),
                  "90"  = "rotated 90° counter-clockwise",
                  "-90" = "rotated 90° clockwise",
                  "270" = "rotated 90° clockwise",
                  "180" = "rotated 180°",
                  NULL)
    bits <- c(rot, if (isTRUE(vt$flip)) "mirrored")
    paste0(" — view ", paste(bits, collapse = " + "), " to match Raman")
  }

  # Nearest particle under the cursor, in the coordinates actually plotted.
  .hover_pick <- function(hover, df, aligned, vt) {
    if (is.null(hover) || is.null(df) || nrow(df) == 0) return(NULL)
    if (aligned && "x" %in% names(df) && any(!is.na(df$x))) {
      hx <- df$x; hy <- df$y
    } else {
      hx <- df$x_orig; hy <- df$y_orig
      if (!.tf_is_identity(vt)) {
        rc <- view_transform_xy(hx, hy, vt$deg, vt$flip); hx <- rc$x; hy <- rc$y
      }
    }
    dists <- sqrt((hx - hover$x)^2 + (hy - hover$y)^2)
    threshold <- max(diff(range(hx, na.rm = TRUE)),
                     diff(range(hy, na.rm = TRUE)), 500) * 0.05
    idx <- which.min(dists)
    if (length(idx) > 0 && dists[idx] <= threshold) df[idx, , drop = FALSE] else NULL
  }

  # .find_nearest() searches x_orig/y_orig, so restate the plotted coordinates
  # there before a click lookup (aligned frame, or view-transformed native).
  .view_click_df <- function(df, coord_mode, vt) {
    if (is.null(df) || nrow(df) == 0) return(df)
    if (!is.null(coord_mode) && coord_mode == "aligned" &&
        "x" %in% names(df) && any(!is.na(df$x))) {
      df$x_orig <- df$x; df$y_orig <- df$y
      return(df)
    }
    if (!.tf_is_identity(vt)) {
      rc <- view_transform_xy(df$x_orig, df$y_orig, vt$deg, vt$flip)
      df$x_orig <- rc$x; df$y_orig <- rc$y
    }
    df
  }

  ftir_plot_obj <- reactive({
    coord_mode <- input$ftir_coord_mode
    aligned    <- !is.null(coord_mode) && coord_mode == "aligned"

    df <- ftir_filtered()
    df_disp <- df
    if (nrow(df_disp) > 0) {
      if (aligned && "x" %in% names(df_disp) && any(!is.na(df_disp$x))) {
        df_disp$x_orig <- df_disp$x
        df_disp$y_orig <- df_disp$y
      } else {
        df_disp$x <- df_disp$x_orig
        df_disp$y <- df_disp$y_orig
      }
    }

    img <- if (aligned) overlay_image_info() else ftir_native_image_info()

    full_ftir <- ftir_df_full()
    if (!is.null(full_ftir) && nrow(full_ftir) > 0) {
      if (aligned && "x" %in% names(full_ftir) && any(!is.na(full_ftir$x))) {
        full_ftir$x_orig <- full_ftir$x
        full_ftir$y_orig <- full_ftir$y
      } else {
        full_ftir$x <- full_ftir$x_orig
        full_ftir$y <- full_ftir$y_orig
      }
    }

    # Rotate/mirror the whole native scene (image + particles) into the Raman
    # orientation. Aligned mode is already in Raman space, so nothing applies.
    vt <- if (aligned) list(deg = 0L, flip = FALSE) else ftir_view_tf()
    img       <- .apply_view_tf_img(img, vt)
    df_disp   <- .apply_view_tf(df_disp, vt)
    full_ftir <- .apply_view_tf(full_ftir, vt)

    # Default viewport spans ALL this instrument's particles, never just the
    # filtered subset — otherwise unticking a material or moving a slider
    # silently re-fits the axes, which reads as an unrequested zoom. Zoom
    # changes only on brush, double-click and Reset Zoom.
    bounds <- if (!is.null(zoom$ftir)) zoom$ftir else {
      ref <- if (!is.null(full_ftir) && nrow(full_ftir) > 0) full_ftir
             else if (nrow(df_disp) > 0) df_disp
             else NULL
      if (!is.null(ref) && any(is.finite(ref$x))) {
        pad <- 300
        list(x = c(min(ref$x, na.rm=TRUE) - pad, max(ref$x, na.rm=TRUE) + pad),
             y = c(min(ref$y, na.rm=TRUE) - pad, max(ref$y, na.rm=TRUE) + pad))
      } else list(x = c(0, 10000), y = c(0, 10000))
    }
    bounds <- sanitize_bounds(bounds, list(x = c(0, 10000), y = c(0, 10000)))

    if (nrow(df_disp) == 0) {
      # Keep the panel, image and viewport exactly as when particles are shown,
      # so filtering down to nothing does not make the view jump. Mirrors the
      # Raman and Bruker tabs.
      empty_title <- if (is.null(full_ftir) || nrow(full_ftir) == 0)
        "FTIR — no data loaded"
      else
        paste0("FTIR — no particles match the current filters (of ",
               nrow(full_ftir), ")")
      p <- ggplot() + coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = empty_title, x = "X (µm)", y = "Y (µm)") +
        theme_minimal(base_size = 15) +
        theme(plot.background = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA))
      return(add_image_bg(p, img))
    }

    title_suffix <- if (aligned) " (Raman-aligned frame)" else .view_tf_suffix(vt)
    hl_single <- input$ftir_highlight_particle
    hl_ids <- if (!is.null(hl_single) && hl_single != "None") {
      unique(c(hl_single, single_highlight_ids$ftir))
    } else single_highlight_ids$ftir

    make_scatter(df_disp, img, bounds,
                 paste0("FTIR Particles (", nrow(df_disp), " shown)", title_suffix),
                 match_colours = c(matched = "#2ca02c", unmatched = "#d62728"),
                 match_labels  = c(matched = "matched to Raman", unmatched = "unmatched"),
                 highlight_id  = hl_ids,
                 full_df       = full_ftir,
                 plain         = isTRUE(input$ftir_show_all_detected),
                 show_labels   = isTRUE(input$ftir_show_all_labels),
                 label_size    = input$ftir_label_size %||% 3)
  })

  output$ftir_plot <- renderPlot(ftir_plot_obj()) |> bindCache(
    # Cache key must list EVERY input this render reads: an omission both
    # serves a stale plot and stops the render invalidating. ftir_filtered()
    # transitively captures the FTIR quality/size/material filters; run-scoped
    # data (df_full, manifest) and image pixels are captured by selected_run_dir().
    selected_run_dir(), is.null(uploaded_data()),
    ftir_filtered(), input$ftir_coord_mode,
    input$ftir_view_rotation, input$ftir_view_flip_y,
    input$ftir_highlight_particle, single_highlight_ids$ftir,
    input$ftir_show_all_detected, input$ftir_show_all_labels, input$ftir_label_size, zoom$ftir,
    img_key(ftir_native_image_info()), img_key(overlay_image_info())
  )

  output$ftir_summary_text <- renderText({
    df <- ftir_filtered()
    if (nrow(df) == 0) return("No pipeline data loaded")
    paste0(nrow(df), " particles | ",
           sum(df$match_status == "matched"), " matched to Raman | ",
           length(unique(df$material)), " materials")
  })

  observeEvent(input$ftir_hover, {
    aligned <- !is.null(input$ftir_coord_mode) && input$ftir_coord_mode == "aligned"
    hit <- .hover_pick(input$ftir_hover, ftir_filtered(), aligned, ftir_view_tf())
    if (!is.null(hit)) last_hover$ftir <- hit
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
                      mat_keep("raman_material_filter"), eff_match_filter("raman"))
  })

  raman_plot_obj <- reactive({
    df <- raman_filtered()
    # Display in native Raman instrument frame (x_orig, y_orig)
    df_disp <- df
    if (nrow(df_disp) > 0) { df_disp$x <- df_disp$x_orig; df_disp$y <- df_disp$y_orig }

    img <- raman_native_image_info()

    # Use image extent for bounds when available (consistent with overlay)
    # Full (unfiltered) Raman data — highlight fallback, and the reference for
    # the default viewport.
    full_raman <- raman_df_full()
    if (!is.null(full_raman) && nrow(full_raman) > 0) {
      full_raman$x <- full_raman$x_orig; full_raman$y <- full_raman$y_orig
    }

    # Viewport spans all particles, not the filtered subset — see the FTIR tab.
    bounds <- if (!is.null(zoom$raman)) zoom$raman else {
      ref <- if (!is.null(full_raman) && nrow(full_raman) > 0) full_raman
             else if (nrow(df_disp) > 0) df_disp else NULL
      if (!is.null(img)) {
        pad <- 200
        list(x = c(img$xmin - pad, img$xmax + pad),
             y = c(img$ymin - pad, img$ymax + pad))
      } else if (!is.null(ref) && any(is.finite(ref$x))) {
        pad <- 300
        list(x = c(min(ref$x, na.rm = TRUE) - pad, max(ref$x, na.rm = TRUE) + pad),
             y = c(min(ref$y, na.rm = TRUE) - pad, max(ref$y, na.rm = TRUE) + pad))
      } else list(x = c(-1000, 1000), y = c(-1000, 1000))
    }
    bounds <- sanitize_bounds(bounds)

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
                 full_df = full_raman,
                 plain = isTRUE(input$raman_show_all_detected),
                 show_labels = isTRUE(input$raman_show_all_labels),
                 label_size = input$raman_label_size %||% 3,
                 subtitle = raman_scale_warning())
  })

  output$raman_plot <- renderPlot(raman_plot_obj()) |> bindCache(
    selected_run_dir(), is.null(uploaded_data()),
    raman_filtered(), raman_scale_warning(),
    input$raman_highlight_particle, single_highlight_ids$raman,
    input$raman_show_all_detected, input$raman_show_all_labels, input$raman_label_size, zoom$raman,
    img_key(raman_native_image_info())
  )

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
                            mat_keep("ldir_material_filter"), eff_match_filter("ldir"))
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

  # Orientation of the LDIR native display: rotation (multiple of 90 deg) plus
  # an optional Y mirror. "auto" measures it from the per-particle pairing each
  # row already carries — native x_orig/y_orig against aligned x/y — so it
  # reproduces exactly the orientation that Aligned (Raman space) mode uses.
  # Cloud matching is only the fallback, for runs whose LDIR was never aligned
  # to Raman; on its own it could score several rotations alike and leave the
  # view unrotated even when aligned mode was plainly correct.
  ldir_view_tf <- reactive({
    sel  <- input$ldir_view_rotation
    flip <- isTRUE(input$ldir_view_flip_y)
    if (!is.null(sel) && sel != "auto")
      return(list(deg = as.integer(sel), flip = flip))

    ld <- ldir_df_full()
    auto <- if (!is.null(ld) && nrow(ld) > 0 &&
                all(c("x", "y", "x_orig", "y_orig") %in% names(ld)))
              auto_view_from_pairs(ld$x_orig, ld$y_orig, ld$x, ld$y) else NULL
    if (!is.null(auto) && is.finite(auto$rmse) && auto$rmse <= .VIEW_FIT_MAX_RMSE)
      return(list(deg = auto$deg, flip = xor(isTRUE(auto$flip), flip)))

    rd <- raman_df_full()
    if (!is.null(ld) && !is.null(rd) && nrow(ld) > 0 && nrow(rd) > 0) {
      a2 <- auto_view_dihedral(ld$x_orig, ld$y_orig, rd$x_orig, rd$y_orig)
      return(list(deg = a2$deg, flip = xor(isTRUE(a2$flip), flip)))
    }
    list(deg = 0L, flip = flip)
  })

  ldir_plot_obj <- reactive({
    df <- ldir_filtered()
    overlay_mode <- input$ldir_overlay_mode

    coord_mode <- input$ldir_coord_mode
    aligned    <- !is.null(coord_mode) && coord_mode == "aligned"

    df_disp <- df
    if (nrow(df_disp) > 0) {
      if (aligned && "x" %in% names(df_disp) && any(!is.na(df_disp$x))) {
        df_disp$x_orig <- df_disp$x
        df_disp$y_orig <- df_disp$y
      } else {
        df_disp$x <- df_disp$x_orig
        df_disp$y <- df_disp$y_orig
      }
    }
    # Optionally drop unmatched (single-instrument) particles for a clean view
    if (isTRUE(input$ldir_hide_unmatched) && nrow(df_disp) > 0)
      df_disp <- df_disp[df_disp$match_status == "matched", ]

    # Background image. "auto" = LDIR image in native mode, Raman image in
    # aligned mode. The explicit choices let the user test whether LDIR points
    # land on the Raman image: "raman" shows the Raman micrograph (placed in
    # the aligned/normalized frame, or — in native mode — scaled to the LDIR
    # particle extent so the two patterns can be compared).
    bg_sel <- input$ldir_bg_image %||% "auto"
    raman_bg_native <- function() {
      raw_r <- raman_image()
      # Scale to the FULL LDIR extent, not the filtered subset, so the
      # background does not move when materials/filters are toggled.
      ld_all <- ldir_df_full()
      if (is.null(raw_r) || is.null(ld_all) || nrow(ld_all) == 0) return(NULL)
      b <- compute_image_bounds(raw_r, ld_all$x_orig, ld_all$y_orig, padding_um = 300)
      list(raster = raw_r, xmin = b$xmin, xmax = b$xmax, ymin = b$ymin, ymax = b$ymax)
    }
    img <- switch(bg_sel,
      none  = NULL,
      raman = if (aligned) overlay_image_info() else raman_bg_native(),
      ldir  = if ("processed_image" %in% overlay_mode) ldir_processed_image_info()
              else ldir_native_image_info(),
      # auto (default)
      if (aligned) overlay_image_info()
      else if ("processed_image" %in% overlay_mode) ldir_processed_image_info()
      else if ("raw_image" %in% overlay_mode) ldir_native_image_info()
      else NULL)

    # Rotate/mirror the whole native scene (image + particles) into the Raman
    # orientation for side-by-side comparison.  Aligned mode is already in
    # Raman space, so nothing applies there.  Display-only.
    vt <- if (aligned) list(deg = 0L, flip = FALSE) else ldir_view_tf()
    img     <- .apply_view_tf_img(img, vt)
    df_disp <- .apply_view_tf(df_disp, vt)

    # All LDIR particles in the same display frame as df_disp — the reference
    # for the default viewport, so filtering never re-fits the axes.
    ref_full <- ldir_df_full()
    if (!is.null(ref_full) && nrow(ref_full) > 0) {
      if (!(aligned && "x" %in% names(ref_full) && any(!is.na(ref_full$x)))) {
        ref_full$x <- ref_full$x_orig; ref_full$y <- ref_full$y_orig
      }
      ref_full <- .apply_view_tf(ref_full, vt)
    } else ref_full <- NULL

    # Viewport priority:
    # 1. User zoom (brush) — always honoured
    # 2. Image bounds — when an image is shown, the viewport must cover the full
    #    scan circle; basing it on sparse particle coords distorts the image.
    # 3. Full particle data range — fallback when no image is loaded. Uses every
    #    particle, not the filtered subset, so toggling filters does not zoom.
    ref_b <- if (!is.null(ref_full)) ref_full else df_disp
    bounds <- if (!is.null(zoom$ldir)) zoom$ldir else if (!is.null(img)) {
      list(x = c(img$xmin, img$xmax), y = c(img$ymin, img$ymax))
    } else if (nrow(ref_b) > 0 && any(is.finite(ref_b$x))) {
      pad <- 500
      list(x = c(min(ref_b$x, na.rm = TRUE) - pad, max(ref_b$x, na.rm = TRUE) + pad),
           y = c(min(ref_b$y, na.rm = TRUE) - pad, max(ref_b$y, na.rm = TRUE) + pad))
    } else list(x = c(-7000, 7000), y = c(-7000, 7000))
    bounds <- sanitize_bounds(bounds, list(x = c(-7000, 7000), y = c(-7000, 7000)))

    n_extracted <- 0
    extracted <- ldir_extracted_pts()
    if (!is.null(extracted)) n_extracted <- nrow(extracted)

    title_parts <- paste0("LDIR Particles (", nrow(df_disp), " shown")
    if ("extracted_pts" %in% overlay_mode && n_extracted > 0)
      title_parts <- paste0(title_parts, " + ", n_extracted, " image-extracted")
    title_parts <- paste0(title_parts, ")")
    title_parts <- paste0(title_parts, .view_tf_suffix(vt))

    if (nrow(df_disp) == 0 && !("extracted_pts" %in% overlay_mode && n_extracted > 0)) {
      p <- ggplot() + coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
        labs(title = "LDIR — no particles loaded", x = "X (\u00b5m)", y = "Y (\u00b5m)") +
        theme_minimal(base_size = 15) +
        theme(plot.background = element_rect(fill = "white", colour = NA),
              panel.background = element_rect(fill = "grey98", colour = NA))
      return(add_image_bg(p, img))
    }

    p <- ggplot() +
      scale_x_continuous(breaks = breaks_adaptive(bounds$x)) +
      scale_y_continuous(breaks = breaks_adaptive(bounds$y)) +
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

    # Raman partners + match lines (aligned mode): draw each matched LDIR
    # particle's Raman partner as a hollow blue circle and connect the two, so
    # a match reads as "two dots joined by a short line" and the residual
    # LDIR<->Raman centroid scatter (~130 um) is visible rather than mistaken
    # for misalignment.
    if (aligned && isTRUE(input$ldir_show_raman_partners)) {
      rd <- tryCatch(run_data_gated()$ldir_raman_matched, error = function(e) NULL)
      need <- c("ldir_particle_id", "ldir_x_aligned", "ldir_y_aligned",
                "raman_x_norm", "raman_y_norm")
      if (!is.null(rd) && nrow(rd) > 0 && all(need %in% names(rd))) {
        # Genuine matches only: a partner line for a forced over-gate pair would
        # contradict its "unmatched" status elsewhere.
        if ("within_gate" %in% names(rd))
          rd <- rd[!is.na(rd$within_gate) & rd$within_gate, , drop = FALSE]
        if (nrow(rd) > 0 && nrow(df_disp) > 0)
          rd <- rd[rd$ldir_particle_id %in% df_disp$particle_id, , drop = FALSE]
        if (nrow(rd) > 0) {
          seg <- data.frame(x = rd$ldir_x_aligned, y = rd$ldir_y_aligned,
                            xend = rd$raman_x_norm, yend = rd$raman_y_norm)
          p <- p +
            geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend),
                         colour = "#00CED1", linewidth = 0.4, alpha = 0.85) +
            geom_point(data = seg, aes(x = xend, y = yend), shape = 1,
                       colour = "#1f77b4", size = 3, stroke = 1)
        }
      }
    }

    # Image-extracted particles (before join): small open circles
    if ("extracted_pts" %in% overlay_mode && n_extracted > 0) {
      ext_df <- data.frame(x = extracted$x_um, y = extracted$y_um,
                           feret_max = extracted$feret_max_um)
      ext_df <- .apply_view_tf(ext_df, vt)
      p <- p + geom_point(data = ext_df,
                            aes(x = x, y = y, size = feret_max),
                            shape = 1, colour = "#e377c2", alpha = 0.5,
                            stroke = 0.5)
    }

    # Excel-joined particles (main layer). "Show all detected" draws every
    # particle one colour with no matched/unmatched distinction.
    if (nrow(df_disp) > 0) {
      if (isTRUE(input$ldir_show_all_detected)) {
        p <- p + geom_point(data = df_disp,
                              aes(x = x, y = y, size = feret_max),
                              colour = "#1f77b4", alpha = 0.7)
      } else {
        p <- p + geom_point(data = df_disp,
                              aes(x = x, y = y, colour = match_status,
                                  size = feret_max),
                              alpha = 0.7) +
          scale_colour_manual(values = c(matched = "#d62728",
                                          unmatched = "#bcbd22"),
                               labels = c(matched = "matched to Raman",
                                          unmatched = "unmatched"))
      }
    }

    # Size legend (single scale for all layers). Union the Feret values of the
    # size-mapped layers so a single displayed particle can't collapse the
    # domain to zero width (NaN rescale -> non-finite viewport crash).
    size_vals <- c(
      if ("extracted_pts" %in% overlay_mode && n_extracted > 0) extracted$feret_max_um,
      if (nrow(df_disp) > 0) df_disp$feret_max
    )
    p <- p + scale_size_continuous(name = "Feret Max (\u00b5m)", range = c(2, 12),
                                   limits = safe_size_limits(size_vals))

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

    # Number every displayed LDIR particle (opt-in)
    if (isTRUE(input$ldir_show_all_labels) && nrow(df_disp) > 0) {
      bnds <- if (!is.null(zoom$ldir)) zoom$ldir else
        list(x = range(df_disp$x, na.rm = TRUE), y = range(df_disp$y, na.rm = TRUE))
      p <- add_particle_labels(p, df_disp, bnds, size = input$ldir_label_size %||% 3)
    }

    p
  })

  output$ldir_plot <- renderPlot(ldir_plot_obj()) |> bindCache(
    # ldir_filtered() captures the LDIR filters; ldir_view_tf() and
    # ldir_extracted_pts() are reactives whose values fold in their own inputs;
    # the three image sources are folded in cheaply via img_key().
    selected_run_dir(), is.null(uploaded_data()),
    ldir_filtered(), ldir_extracted_pts(), ldir_view_tf(),
    input$ldir_bg_image, input$ldir_coord_mode, input$ldir_hide_unmatched,
    input$ldir_overlay_mode, input$ldir_show_raman_partners,
    input$ldir_show_all_detected, input$ldir_show_all_labels, input$ldir_label_size,
    input$ldir_highlight_particle,
    single_highlight_ids$ldir, zoom$ldir,
    img_key(ldir_native_image_info()), img_key(ldir_processed_image_info()),
    img_key(overlay_image_info())
  )

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
    df <- .ldir_click_df()   # native coords, view-rotated like the display
    if (nrow(df) == 0) return()
    dists <- sqrt((df$x_orig - hover$x)^2 + (df$y_orig - hover$y)^2)
    idx   <- which.min(dists)
    threshold <- max(diff(range(df$x_orig, na.rm = TRUE)),
                     diff(range(df$y_orig, na.rm = TRUE)), 500) * 0.05
    # Store the UNrotated row so the detail table shows true coordinates
    if (dists[idx] <= threshold) last_hover$ldir <- ldir_filtered()[idx, , drop = FALSE]
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
                      mat_keep("ftir_bruker_material_filter"), eff_match_filter("ftir_bruker"))
  })

  ftir_bruker_plot_obj <- reactive({
    coord_mode <- input$ftir_bruker_coord_mode
    aligned    <- !is.null(coord_mode) && coord_mode == "aligned"

    df <- ftir_bruker_filtered()
    df_disp <- df
    if (nrow(df_disp) > 0) {
      if (aligned && "x" %in% names(df_disp) && any(!is.na(df_disp$x))) {
        df_disp$x_orig <- df_disp$x
        df_disp$y_orig <- df_disp$y
      } else {
        df_disp$x <- df_disp$x_orig
        df_disp$y <- df_disp$y_orig
      }
    }

    img <- if (aligned) overlay_image_info() else ftir_bruker_native_image_info()

    full_fb <- ftir_bruker_df_full()
    if (!is.null(full_fb) && nrow(full_fb) > 0) {
      if (aligned && "x" %in% names(full_fb) && any(!is.na(full_fb$x))) {
        full_fb$x_orig <- full_fb$x; full_fb$y_orig <- full_fb$y
      } else {
        full_fb$x <- full_fb$x_orig; full_fb$y <- full_fb$y_orig
      }
    }

    # Rotate/mirror the whole native scene into the Raman orientation.
    vt <- if (aligned) list(deg = 0L, flip = FALSE) else ftir_bruker_view_tf()
    img     <- .apply_view_tf_img(img, vt)
    df_disp <- .apply_view_tf(df_disp, vt)
    full_fb <- .apply_view_tf(full_fb, vt)

    # Viewport spans all particles, not the filtered subset — see the FTIR tab.
    bounds <- if (!is.null(zoom$ftir_bruker)) zoom$ftir_bruker else {
      ref <- if (!is.null(full_fb) && nrow(full_fb) > 0) full_fb
             else if (nrow(df_disp) > 0) df_disp
             else NULL
      if (!is.null(ref) && any(is.finite(ref$x))) {
        pad <- 300
        list(x = c(min(ref$x, na.rm=TRUE) - pad, max(ref$x, na.rm=TRUE) + pad),
             y = c(min(ref$y, na.rm=TRUE) - pad, max(ref$y, na.rm=TRUE) + pad))
      } else list(x = c(0, 10000), y = c(0, 10000))
    }
    bounds <- sanitize_bounds(bounds, list(x = c(0, 10000), y = c(0, 10000)))

    if (nrow(df_disp) == 0) {
      return(ggplot() +
        coord_fixed(xlim=bounds$x, ylim=bounds$y, expand=FALSE) +
        labs(title="FTIR (Bruker) — no data loaded", x="X (µm)", y="Y (µm)") +
        theme_minimal(base_size=15) +
        theme(plot.background=element_rect(fill="white", colour=NA),
              panel.background=element_rect(fill="grey98", colour=NA)))
    }

    title_suffix <- if (aligned) " (Raman-aligned frame)" else .view_tf_suffix(vt)
    hl_single <- input$ftir_bruker_highlight_particle
    hl_ids <- if (!is.null(hl_single) && hl_single != "None") {
      unique(c(hl_single, single_highlight_ids$ftir_bruker))
    } else single_highlight_ids$ftir_bruker

    make_scatter(df_disp, img, bounds,
                 paste0("FTIR (Bruker) Particles (", nrow(df_disp), " shown)", title_suffix),
                 match_colours = c(matched="#9467bd", unmatched="#8c564b"),
                 match_labels  = c(matched="matched to Raman", unmatched="unmatched"),
                 highlight_id  = hl_ids,
                 full_df       = full_fb,
                 plain         = isTRUE(input$ftir_bruker_show_all_detected),
                 show_labels   = isTRUE(input$ftir_bruker_show_all_labels),
                 label_size    = input$ftir_bruker_label_size %||% 3)
  })

  output$ftir_bruker_plot <- renderPlot(ftir_bruker_plot_obj()) |> bindCache(
    selected_run_dir(), is.null(uploaded_data()),
    ftir_bruker_filtered(), input$ftir_bruker_coord_mode,
    input$ftir_bruker_view_rotation, input$ftir_bruker_view_flip_y,
    input$ftir_bruker_highlight_particle, single_highlight_ids$ftir_bruker,
    input$ftir_bruker_show_all_detected, input$ftir_bruker_show_all_labels, input$ftir_bruker_label_size,
    zoom$ftir_bruker,
    img_key(ftir_bruker_native_image_info()), img_key(overlay_image_info())
  )
  output$ftir_bruker_summary_text <- renderText({
    df <- ftir_bruker_filtered()
    if (nrow(df) == 0) return("No FTIR (Bruker) data loaded")
    n_matched <- sum(df$match_status == "matched", na.rm = TRUE)
    if (n_matched > 0) {
      paste0(nrow(df), " particles | ", n_matched, " matched to Raman | ",
             length(unique(df$material)), " materials")
    } else {
      paste0(nrow(df), " particles | ", length(unique(df$material)), " materials")
    }
  })

  observeEvent(input$ftir_bruker_hover, {
    aligned <- !is.null(input$ftir_bruker_coord_mode) &&
               input$ftir_bruker_coord_mode == "aligned"
    hit <- .hover_pick(input$ftir_bruker_hover, ftir_bruker_filtered(),
                       aligned, ftir_bruker_view_tf())
    if (!is.null(hit)) last_hover$ftir_bruker <- hit
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

    # Per-instrument material filters (NULL = no filtering; see mat_keep())
    ftir_mat  <- mat_keep("overlay_ftir_material")
    raman_mat <- mat_keep("overlay_raman_material")
    if (!is.null(ftir_mat) || !is.null(raman_mat)) {
      keep <- rep(TRUE, nrow(df))
      if (!is.null(ftir_mat))  keep <- keep & (df$ftir_material_family %in% ftir_mat)
      if (!is.null(raman_mat)) keep <- keep & (df$raman_material_family %in% raman_mat)
      df <- df[keep, ]
    }
    df
  })

  # LDIR-Raman matched data for overlay (filtered by per-instrument controls)
  overlay_ldir_matched <- reactive({
    d <- run_data_gated()
    if (is.null(d$ldir_raman_matched) || nrow(d$ldir_raman_matched) == 0)
      return(data.frame())
    df <- d$ldir_raman_matched

    # Genuine matches only: drop forced over-gate pairings so the overlay draws
    # the same matched set the summary counts. Over-gate LDIR then fall through
    # to the unmatched layer (lone red) and their Raman partners to unmatched
    # Raman (lone blue) — image and table stay consistent.
    if ("within_gate" %in% names(df)) {
      df <- df[!is.na(df$within_gate) & df$within_gate, , drop = FALSE]
      if (nrow(df) == 0) return(data.frame())
    }

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

    # LDIR material filter (harmonized family names)
    ldir_mat <- mat_keep("overlay_ldir_material")
    if (!is.null(ldir_mat) && "ldir_material_family" %in% names(df)) {
      df <- df[df$ldir_material_family %in% ldir_mat, ]
    }

    df
  })

  # Bruker FTIR-Raman matched data for overlay (filtered by per-instrument controls)
  overlay_bruker_matched <- reactive({
    d <- run_data()
    if (is.null(d$matched_ftir_bruker) || nrow(d$matched_ftir_bruker) == 0)
      return(data.frame())
    df <- d$matched_ftir_bruker

    bruker_q <- overlay_ftir_bruker_quality_d()
    if (!is.null(bruker_q) && "ftir_quality" %in% names(df)) {
      df <- df[!is.na(df$ftir_quality) &
               df$ftir_quality >= bruker_q[1] &
               df$ftir_quality <= bruker_q[2], ]
    }

    bruker_sz <- overlay_ftir_bruker_size_d()
    if (!is.null(bruker_sz) && "ftir_feret_max_um" %in% names(df)) {
      df <- df[!is.na(df$ftir_feret_max_um) &
               df$ftir_feret_max_um >= bruker_sz[1] &
               df$ftir_feret_max_um <= bruker_sz[2], ]
    }

    dist_r <- overlay_dist_range_d()
    if (!is.null(dist_r) && "match_distance" %in% names(df)) {
      df <- df[!is.na(df$match_distance) &
               df$match_distance >= dist_r[1] &
               df$match_distance <= dist_r[2], ]
    }

    bruker_mat <- mat_keep("overlay_ftir_bruker_material")
    if (!is.null(bruker_mat) && "ftir_material_family" %in% names(df)) {
      df <- df[df$ftir_material_family %in% bruker_mat, ]
    }

    df
  })

  # Triple-match data: particles detected by all three instruments
  overlay_triplets <- reactive({
    d <- run_data_gated()
    if (is.null(d$triplets) || nrow(d$triplets) == 0) return(data.frame())
    tr <- d$triplets
    # A triplet is only real if its LDIR↔Raman leg is a genuine (within-gate)
    # match — otherwise the LDIR partner is a forced over-gate assignment and
    # the "triple" is spurious. Keeps the triple count honest alongside the
    # gated LDIR↔Raman count.
    gate <- d$ldir_match_gate_um
    if (!is.null(gate) && is.finite(gate) && "ldir_raman_distance" %in% names(tr)) {
      tr <- tr[!is.na(tr$ldir_raman_distance) & tr$ldir_raman_distance <= gate, , drop = FALSE]
    }
    tr
  })

  # ---- Cross-instrument (non-Raman) pair registry --------------------------
  # The three Raman pairs keep their dedicated reactives (overlay_matched /
  # overlay_ldir_matched / overlay_bruker_matched). These specs cover the other
  # three pairs so the overlay reacts to every instrument combination. Each spec
  # names, per endpoint, the common-frame coordinate/id/Feret/material-family
  # columns of its match table and the overlay material input to honour.
  CROSS_PAIR_SPECS <- list(
    list(table = "matched_bruker_perkin",
         a = list(inst="ftir_bruker", x="ftir_bruker_x_aligned", y="ftir_bruker_y_aligned",
                  id="ftir_bruker_particle_id", feret="ftir_bruker_feret_max_um",
                  fam="ftir_bruker_material_family", matinput="overlay_ftir_bruker_material",
                  label="FTIR (Bruker)"),
         b = list(inst="ftir_pe", x="ftir_perkin_x_aligned", y="ftir_perkin_y_aligned",
                  id="ftir_perkin_particle_id", feret="ftir_perkin_feret_max_um",
                  fam="ftir_perkin_material_family", matinput="overlay_ftir_material",
                  label="FTIR")),
    list(table = "matched_bruker_ldir",
         a = list(inst="ftir_bruker", x="ftir_bruker_x_aligned", y="ftir_bruker_y_aligned",
                  id="ftir_bruker_particle_id", feret="ftir_bruker_feret_max_um",
                  fam="ftir_bruker_material_family", matinput="overlay_ftir_bruker_material",
                  label="FTIR (Bruker)"),
         b = list(inst="ldir", x="ldir_x_aligned", y="ldir_y_aligned",
                  id="ldir_particle_id", feret="ldir_feret_max_um",
                  fam="ldir_material_family", matinput="overlay_ldir_material",
                  label="LDIR")),
    list(table = "matched_perkin_ldir",
         a = list(inst="ldir", x="ldir_x_aligned", y="ldir_y_aligned",
                  id="ldir_particle_id", feret="ldir_feret_max_um",
                  fam="ldir_material_family", matinput="overlay_ldir_material",
                  label="LDIR"),
         b = list(inst="ftir_pe", x="ftir_x_aligned", y="ftir_y_aligned",
                  id="ftir_particle_id", feret="ftir_feret_max_um",
                  fam="ftir_material_family", matinput="overlay_ftir_material",
                  label="FTIR"))
  )

  # Material-filtered match table for a cross-pair spec, or NULL when either
  # endpoint's instrument is not selected / the table is absent / empty.
  cross_pair_df <- function(spec, inst) {
    if (!(spec$a$inst %in% inst && spec$b$inst %in% inst)) return(NULL)
    df <- run_data_gated()[[spec$table]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    for (ep in list(spec$a, spec$b)) {
      mat <- mat_keep(ep$matinput)
      if (!is.null(mat) && ep$fam %in% names(df))
        df <- df[df[[ep$fam]] %in% mat, , drop = FALSE]
    }
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df
  }

  # All cross-pairs whose two instruments are both currently selected, each with
  # its filtered table. Reused by the plot, the summary and hover/click.
  active_cross_pairs <- function(inst) {
    Filter(Negate(is.null), lapply(CROSS_PAIR_SPECS, function(s) {
      df <- cross_pair_df(s, inst)
      if (is.null(df)) NULL else list(spec = s, df = df)
    }))
  }

  overlay_plot_obj <- reactive({
    # Gate the heaviest render on data being present. run_data() is an empty
    # list() (falsy) only when no run is loaded at all — it stays populated when
    # filters yield zero particles — so this suppresses the pre-data render
    # without hiding any "run loaded, nothing matches" state.
    req(run_data())

    dfs_full <- overlay_dfs_full()
    dfs      <- overlay_dfs()
    matched  <- overlay_matched()
    ldir_m   <- overlay_ldir_matched()
    bruker_m <- overlay_bruker_matched()
    triplets <- overlay_triplets()

    inst <- input$overlay_instruments
    rel  <- input$overlay_relationships
    if (is.null(inst)) inst <- character(0)
    if (is.null(rel))  rel  <- character(0)

    show_pe_raman     <- "ftir_pe"     %in% inst && "raman" %in% inst
    show_bruker_raman <- "ftir_bruker" %in% inst && "raman" %in% inst
    show_ldir_raman   <- "ldir"        %in% inst && "raman" %in% inst

    # Non-Raman pairs (Bruker↔Perkin, Bruker↔LDIR, Perkin↔LDIR) whose two
    # instruments are both selected, each with its material-filtered table.
    cross_active <- active_cross_pairs(inst)

    bounds <- if (!is.null(zoom$overlay)) zoom$overlay else {
      img_info <- overlay_image_info()
      if (!is.null(img_info)) {
        pad <- 200
        list(x = c(img_info$xmin - pad, img_info$xmax + pad),
             y = c(img_info$ymin - pad, img_info$ymax + pad))
      } else {
        compute_bounds(dfs_full$ftir, dfs_full$raman)
      }
    }
    bounds <- sanitize_bounds(bounds)

    p <- ggplot() +
      scale_x_continuous(breaks = breaks_adaptive(bounds$x)) +
      scale_y_continuous(breaks = breaks_adaptive(bounds$y)) +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = "Multi-Instrument Overlay (aligned coordinates)",
           x = "X (µm)", y = "Y (µm)") +
      theme_minimal(base_size = 15) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "grey98", colour = NA),
        panel.grid       = element_line(colour = "grey90"),
        legend.position  = "right",
        legend.title     = element_text(size = 13),
        legend.text      = element_text(size = 11)
      )

    p <- add_image_bg(p, overlay_image_info())

    # Match lines
    if ("lines" %in% rel) {
      if (show_pe_raman && nrow(matched) > 0) {
        seg_df <- data.frame(
          x = matched$ftir_x_aligned, y = matched$ftir_y_aligned,
          xend = matched$raman_x_norm, yend = matched$raman_y_norm)
        p <- p + geom_segment(data = seg_df, aes(x=x, y=y, xend=xend, yend=yend),
                               colour = "grey40", alpha = 0.6, linewidth = 0.8)
      }
      if (show_bruker_raman && nrow(bruker_m) > 0 &&
          "ftir_x_aligned" %in% names(bruker_m) && "raman_x_norm" %in% names(bruker_m)) {
        bruker_seg <- data.frame(
          x = bruker_m$ftir_x_aligned, y = bruker_m$ftir_y_aligned,
          xend = bruker_m$raman_x_norm, yend = bruker_m$raman_y_norm)
        p <- p + geom_segment(data = bruker_seg, aes(x=x, y=y, xend=xend, yend=yend),
                               colour = "#9467bd", alpha = 0.6, linewidth = 0.8)
      }
      if (show_ldir_raman && nrow(ldir_m) > 0 &&
          "ldir_x_aligned" %in% names(ldir_m) && "raman_x_norm" %in% names(ldir_m)) {
        ldir_seg <- data.frame(
          x = ldir_m$ldir_x_aligned, y = ldir_m$ldir_y_aligned,
          xend = ldir_m$raman_x_norm, yend = ldir_m$raman_y_norm)
        p <- p + geom_segment(data = ldir_seg, aes(x=x, y=y, xend=xend, yend=yend),
                               colour = "#d62728", alpha = 0.5, linewidth = 0.7)
      }
      # Cross-instrument (non-Raman) pair lines
      for (ca in cross_active) {
        s <- ca$spec; cdf <- ca$df
        if (all(c(s$a$x, s$a$y, s$b$x, s$b$y) %in% names(cdf))) {
          cseg <- data.frame(x = cdf[[s$a$x]], y = cdf[[s$a$y]],
                             xend = cdf[[s$b$x]], yend = cdf[[s$b$y]])
          p <- p + geom_segment(data = cseg, aes(x=x, y=y, xend=xend, yend=yend),
                                 colour = "grey55", alpha = 0.5, linewidth = 0.6)
        }
      }
    }

    all_pts <- list()

    .inst_label <- function(k) switch(k, ftir_pe="FTIR", ftir_bruker="FTIR (Bruker)",
                                         raman="Raman", ldir="LDIR", k)
    .inst_df    <- function(k) switch(k, ftir_pe=dfs$ftir, ftir_bruker=dfs$ftir_bruker,
                                         raman=dfs$raman, ldir=dfs$ldir, NULL)
    .inst_mat   <- function(k) switch(k, ftir_pe="overlay_ftir_material",
                                         ftir_bruker="overlay_ftir_bruker_material",
                                         raman="overlay_raman_material",
                                         ldir="overlay_ldir_material", NULL)

    if (length(inst) == 1) {
      # Single-instrument mode: show ALL particles of that instrument,
      # ignoring the matched/unmatched relationship filter.
      df_s <- .inst_df(inst)
      if (!is.null(df_s) && nrow(df_s) > 0) {
        mat <- mat_keep(.inst_mat(inst))
        if (!is.null(mat)) df_s <- df_s[df_s$material_family %in% mat, ]
        if (nrow(df_s) > 0)
          all_pts[[1]] <- data.frame(x=df_s$x, y=df_s$y, feret_max=df_s$feret_max,
                                     instrument=.inst_label(inst),
                                     match_status=df_s$match_status,
                                     particle_id=as.character(df_s$particle_id),
                                     stringsAsFactors=FALSE)
      }
    } else {
      # Matched status is evaluated relative to the SELECTED pairings, not each
      # instrument's primary match. So in a LDIR+Raman overlay, matched = the
      # LDIR<->Raman pairs (both sides drawn) and "unmatched Raman" = Raman
      # particles not matched to LDIR (not Raman-vs-FTIR). A particle counts as
      # matched if it is paired with ANY selected partner.
      .as_chr <- function(v) if (is.null(v)) character(0) else as.character(v)
      raman_matched_ids <- unique(c(
        if (show_pe_raman)     .as_chr(matched$raman_particle_id),
        if (show_bruker_raman) .as_chr(bruker_m$raman_particle_id),
        if (show_ldir_raman)   .as_chr(ldir_m$raman_particle_id)))
      ftir_matched_ids   <- if (show_pe_raman)     unique(.as_chr(matched$ftir_particle_id))  else character(0)
      bruker_matched_ids <- if (show_bruker_raman) unique(.as_chr(bruker_m$ftir_particle_id)) else character(0)
      ldir_matched_ids   <- if (show_ldir_raman)   unique(.as_chr(ldir_m$ldir_particle_id))   else character(0)

      # Fold in cross-pair matches: a particle counts as matched (and so drops
      # out of the unmatched layer) if it is paired with ANY selected partner,
      # Raman or otherwise.
      for (ca in cross_active) {
        for (ep in list(ca$spec$a, ca$spec$b)) {
          ids <- .as_chr(ca$df[[ep$id]])
          if      (ep$inst == "ftir_pe")     ftir_matched_ids   <- unique(c(ftir_matched_ids, ids))
          else if (ep$inst == "ftir_bruker") bruker_matched_ids <- unique(c(bruker_matched_ids, ids))
          else if (ep$inst == "ldir")        ldir_matched_ids   <- unique(c(ldir_matched_ids, ids))
          else if (ep$inst == "raman")       raman_matched_ids  <- unique(c(raman_matched_ids, ids))
        }
      }

      if ("matched" %in% rel) {
        .add_matched <- function(label, x, y, feret, id = NULL) {
          if (length(x) == 0) return(invisible())
          all_pts[[length(all_pts)+1]] <<- data.frame(
            x=x, y=y, feret_max=feret, instrument=label,
            match_status="matched",
            particle_id = if (!is.null(id)) as.character(id) else NA_character_,
            stringsAsFactors=FALSE)
        }
        if (show_pe_raman && nrow(matched) > 0) {
          .add_matched("FTIR",  matched$ftir_x_aligned, matched$ftir_y_aligned, matched$ftir_feret_max_um, matched$ftir_particle_id)
          .add_matched("Raman", matched$raman_x_norm,   matched$raman_y_norm,   matched$raman_feret_max_um, matched$raman_particle_id)
        }
        if (show_bruker_raman && nrow(bruker_m) > 0) {
          .add_matched("FTIR (Bruker)", bruker_m$ftir_x_aligned, bruker_m$ftir_y_aligned, bruker_m$ftir_feret_max_um, bruker_m$ftir_particle_id)
          .add_matched("Raman",         bruker_m$raman_x_norm,   bruker_m$raman_y_norm,   bruker_m$raman_feret_max_um, bruker_m$raman_particle_id)
        }
        if (show_ldir_raman && nrow(ldir_m) > 0 && "ldir_x_aligned" %in% names(ldir_m)) {
          .add_matched("LDIR",  ldir_m$ldir_x_aligned, ldir_m$ldir_y_aligned, ldir_m$ldir_feret_max_um, ldir_m$ldir_particle_id)
          .add_matched("Raman", ldir_m$raman_x_norm,   ldir_m$raman_y_norm,   ldir_m$raman_feret_max_um, ldir_m$raman_particle_id)
        }
        # Cross-instrument pairs: draw both endpoints filled.
        for (ca in cross_active) {
          s <- ca$spec; cdf <- ca$df
          if (all(c(s$a$x, s$a$y, s$a$feret) %in% names(cdf)))
            .add_matched(s$a$label, cdf[[s$a$x]], cdf[[s$a$y]], cdf[[s$a$feret]],
                         if (!is.null(s$a$id) && s$a$id %in% names(cdf)) cdf[[s$a$id]] else NULL)
          if (all(c(s$b$x, s$b$y, s$b$feret) %in% names(cdf)))
            .add_matched(s$b$label, cdf[[s$b$x]], cdf[[s$b$y]], cdf[[s$b$feret]],
                         if (!is.null(s$b$id) && s$b$id %in% names(cdf)) cdf[[s$b$id]] else NULL)
        }
      }

      if ("unmatched" %in% rel) {
        # Unmatched for an instrument = its particles not paired with any
        # selected partner (uses particle_id against the pairing's matched ids,
        # so every chosen instrument's unmatched particles are shown).
        .add_unmatched <- function(df_full, inst_key, inst_label, mat_input, matched_ids) {
          if (!(inst_key %in% inst) || is.null(df_full) || nrow(df_full) == 0) return(NULL)
          um <- df_full[!(as.character(df_full$particle_id) %in% matched_ids), ]
          mat <- mat_keep(mat_input)
          if (!is.null(mat) && nrow(um) > 0)
            um <- um[um$material_family %in% mat, ]
          if (nrow(um) == 0) return(NULL)
          data.frame(x=um$x, y=um$y, feret_max=um$feret_max,
                     instrument=inst_label, match_status="unmatched",
                     particle_id = as.character(um$particle_id), stringsAsFactors=FALSE)
        }
        all_pts <- c(all_pts, Filter(Negate(is.null), list(
          .add_unmatched(dfs$ftir,        "ftir_pe",     "FTIR",          "overlay_ftir_material",        ftir_matched_ids),
          .add_unmatched(dfs$raman,       "raman",       "Raman",         "overlay_raman_material",       raman_matched_ids),
          .add_unmatched(dfs$ftir_bruker, "ftir_bruker", "FTIR (Bruker)", "overlay_ftir_bruker_material", bruker_matched_ids),
          .add_unmatched(dfs$ldir,        "ldir",        "LDIR",          "overlay_ldir_material",        ldir_matched_ids)
        )))
      }
    }

    # Size-scale domain, kept stable and positive-width. A tight match gate can
    # leave a single point (or one distinct Feret); scale_size_continuous would
    # then rescale from a zero-width range (0/0 -> NaN) and the NaN-sized legend
    # key throws "non-finite location/size for viewport", killing the whole plot
    # (image and points vanish). NULL until points exist -> default (unused) scale.
    size_limits <- NULL
    if (length(all_pts) > 0) {
      both <- do.call(rbind, all_pts)
      # Defensive: never let a non-finite coordinate reach a geom/viewport.
      both <- both[is.finite(both$x) & is.finite(both$y), , drop = FALSE]
      if (nrow(both) > 0) {
        size_limits <- safe_size_limits(both$feret_max)
        p <- p + geom_point(data = both,
                             aes(x=x, y=y, size=feret_max, colour=instrument, shape=match_status),
                             alpha = 0.65) +
          scale_colour_manual(
            name   = "Instrument",
            values = c(FTIR = "#2ca02c", "FTIR (Bruker)" = "#9467bd",
                       Raman = "#1f77b4", LDIR = "#d62728")
          ) +
          scale_shape_manual(
            name   = "Match Status",
            values = c(matched = 19, unmatched = 1),
            labels = c(matched = "Matched (filled)", unmatched = "Unmatched (open)")
          )
        # Number every displayed particle (opt-in). Endpoints keep each
        # instrument's own particle_id, so a matched pair shows both IDs.
        if (isTRUE(input$overlay_show_all_labels))
          p <- add_particle_labels(p, both, bounds, size = input$overlay_label_size %||% 3)
      }
    }

    # Multi-instrument highlights (3+)
    if ("multi" %in% rel && nrow(triplets) > 0 && nrow(matched) > 0 && nrow(ldir_m) > 0) {
      triple_pts <- list()
      m_trip <- matched[matched$raman_particle_id %in% triplets$raman_particle_id, ]
      if (nrow(m_trip) > 0) {
        triple_pts[[1]] <- data.frame(x=m_trip$ftir_x_aligned, y=m_trip$ftir_y_aligned)
        triple_pts[[2]] <- data.frame(x=m_trip$raman_x_norm,   y=m_trip$raman_y_norm)
      }
      l_trip <- ldir_m[ldir_m$raman_particle_id %in% triplets$raman_particle_id, ]
      if (nrow(l_trip) > 0 && "ldir_x_aligned" %in% names(l_trip))
        triple_pts[[length(triple_pts)+1]] <- data.frame(x=l_trip$ldir_x_aligned, y=l_trip$ldir_y_aligned)
      if (length(triple_pts) > 0) {
        triple_df <- do.call(rbind, triple_pts)
        p <- p + geom_point(data = triple_df, aes(x=x, y=y),
                             shape=21, size=6, stroke=1.5, fill=NA, colour="#FFD700", alpha=0.9)
      }
    }

    p <- p + scale_size_continuous(name = "Feret Max (µm)", range = c(2, 12),
                                   limits = size_limits)

    hl_specs <- list(
      list(ids = input$overlay_ftir_particles,        df = dfs$ftir,        col = "#2ca02c"),
      list(ids = input$overlay_raman_particles,       df = dfs$raman,       col = "#1f77b4"),
      list(ids = input$overlay_ldir_particles,        df = dfs$ldir,        col = "#d62728"),
      list(ids = input$overlay_ftir_bruker_particles, df = dfs$ftir_bruker, col = "#9467bd")
    )
    y_span_ov   <- diff(bounds$y)
    y_nudge_ov  <- y_span_ov * 0.03
    for (spec in hl_specs) {
      sel_ids <- spec$ids
      if (is.null(sel_ids) || length(sel_ids) == 0) next
      inst_df <- spec$df
      if (is.null(inst_df) || nrow(inst_df) == 0) next
      hl <- inst_df[inst_df$particle_id %in% sel_ids, ]
      if (nrow(hl) > 0) {
        hl$label_y <- hl$y + y_nudge_ov
        p <- p + geom_point(data=hl, aes(x=x, y=y), shape=19, size=5, colour=spec$col) +
                 geom_point(data=hl, aes(x=x, y=y), shape=21, size=10, stroke=2,
                             fill=NA, colour="#FFD700") +
                 geom_text(data=hl, aes(x=x, y=label_y, label=particle_id),
                            vjust=0, size=4.0, fontface="bold", colour="#FFD700")
      }
    }

    pin <- pinned_overlay()
    if (!is.null(pin) && !is.null(pin$x) && !is.null(pin$y)) {
      pin_label_y <- pin$y + y_nudge_ov
      pin_df <- data.frame(x=pin$x, y=pin$y, label_y=pin_label_y,
                           label=if (!is.null(pin$particle_id)) pin$particle_id else "")
      p <- p + geom_point(data=pin_df, aes(x=x, y=y), shape=8, size=8, stroke=2, colour="#FF6600") +
               geom_text(data=pin_df, aes(x=x, y=label_y, label=label),
                          vjust=0, size=4.6, fontface="bold", colour="#FF6600")
    }

    p
  })

  output$overlay_plot <- renderPlot(overlay_plot_obj()) |> bindCache(
    # overlay_matched()/overlay_ldir_matched()/overlay_bruker_matched() fold in
    # the per-instrument filters as they affect MATCHED particles; the four
    # material inputs are read dynamically via .inst_mat(inst) in the body, so
    # they must be listed explicitly. Key on mat_keep(), NOT the raw input: a
    # raw NULL means both "not populated yet" (show all) and "every box
    # unticked" (show none), which would collide and serve the wrong plot.
    #
    # The quality/size sliders are keyed explicitly too. They now also filter
    # the UNMATCHED points, and a slider move that drops only unmatched
    # particles leaves overlay_matched() byte-identical -- so keying on the
    # matched frames alone served a stale plot and the sliders looked dead.
    selected_run_dir(), is.null(uploaded_data()),
    overlay_ftir_quality_d(), overlay_ftir_size_d(),
    overlay_raman_quality_d(), overlay_raman_size_d(),
    overlay_ldir_quality_d(), overlay_ldir_size_d(),
    overlay_ftir_bruker_quality_d(), overlay_ftir_bruker_size_d(),
    overlay_matched(), overlay_ldir_matched(), overlay_bruker_matched(),
    overlay_triplets(), input$overlay_instruments, input$overlay_relationships,
    mat_keep("overlay_ftir_material"), mat_keep("overlay_raman_material"),
    mat_keep("overlay_ldir_material"), mat_keep("overlay_ftir_bruker_material"),
    input$overlay_ftir_particles, input$overlay_raman_particles,
    input$overlay_ldir_particles, input$overlay_ftir_bruker_particles,
    input$overlay_show_all_labels, input$overlay_label_size,
    zoom$overlay, pinned_overlay(), img_key(overlay_image_info())
  )
  output$overlay_summary_text <- renderText({
    m       <- overlay_matched()
    bm      <- overlay_bruker_matched()
    dfs     <- overlay_dfs()
    triplets <- overlay_triplets()
    if (nrow(m) == 0 && is.null(dfs$ftir) && nrow(bm) == 0) return("No pipeline data loaded")
    n_um_f    <- if (!is.null(dfs$ftir))        sum(dfs$ftir$match_status == "unmatched")        else 0
    n_um_r    <- if (!is.null(dfs$raman))       sum(dfs$raman$match_status == "unmatched")       else 0
    n_um_fb   <- if (!is.null(dfs$ftir_bruker)) sum(dfs$ftir_bruker$match_status == "unmatched") else 0
    n_ldir    <- if (!is.null(dfs$ldir))        nrow(dfs$ldir)                                   else 0
    n_trip    <- nrow(triplets)
    # Genuine matches only, against the LIVE gate slider: over-gate forced
    # pairings are reported as unmatched LDIR. Count from the re-gated frame so
    # the summary tracks the slider in real time; show the gate for clarity.
    gd        <- run_data_gated()
    lrm       <- gd$ldir_raman_matched
    n_ldir_m  <- if (!is.null(lrm) && nrow(lrm) > 0)
                   sum(if ("within_gate" %in% names(lrm)) (!is.na(lrm$within_gate) & lrm$within_gate)
                       else rep(TRUE, nrow(lrm)))
                 else 0
    gate      <- gd$ldir_match_gate_um
    ldir_lbl  <- if (!is.null(gate) && is.finite(gate))
                   paste0(n_ldir_m, "/", n_ldir, " LDIR\u2194Raman matched (\u2264", gate, "\u00b5m)")
                 else
                   paste0(n_ldir_m, "/", n_ldir, " LDIR\u2194Raman matched")
    # Cross-instrument (non-Raman) pair counts \u2014 only shown when present.
    .cross_n <- function(tbl) { d <- gd[[tbl]]; if (is.null(d)) 0L else nrow(d) }
    cross_bits <- c()
    if (!is.null(gd$matched_bruker_perkin))
      cross_bits <- c(cross_bits, paste0(.cross_n("matched_bruker_perkin"), " Bruker\u2194FTIR pairs"))
    if (!is.null(gd$matched_bruker_ldir))
      cross_bits <- c(cross_bits, paste0(.cross_n("matched_bruker_ldir"), " Bruker\u2194LDIR pairs"))
    if (!is.null(gd$matched_perkin_ldir))
      cross_bits <- c(cross_bits, paste0(.cross_n("matched_perkin_ldir"), " FTIR\u2194LDIR pairs"))
    cross_lbl <- if (length(cross_bits) > 0) paste0(" | ", paste(cross_bits, collapse = " | ")) else ""
    paste0(nrow(m),  " FTIR\u2194Raman pairs | ",
           n_um_f,   " FTIR unmatched | ",
           n_um_r,   " Raman unmatched | ",
           nrow(bm), " Bruker\u2194Raman pairs | ",
           n_um_fb,  " Bruker unmatched | ",
           ldir_lbl, " | ",
           n_trip,   " triple matches",
           cross_lbl)
  })

  # ==================================================================
  # Helper: find nearest particle across all instruments (for click/hover)
  # Returns list(row, source, dist) or NULL.
  # Checks: matched pairs, LDIR-Raman pairs, then highlighted/selected
  # single-instrument particles (regardless of layer state).
  # ==================================================================
  # Helper: given a Raman particle id, gather every instrument partner across the
  # three pairwise match tables. They all share raman_particle_id as the join
  # key, so a Raman particle matched by FTIR, Bruker and LDIR yields all of them.
  # matched_pe/matched_ldir/matched_bruker are the run_data() match frames
  # (may be NULL or empty). Returns single-row partner frames keyed by instrument.
  gather_partners_by_raman <- function(raman_pid, matched_pe, matched_ldir, matched_bruker) {
    out <- list(ftir = NULL, ftir_bruker = NULL, raman = NULL, ldir = NULL)
    if (is.null(raman_pid) || length(raman_pid) == 0 || is.na(raman_pid)) return(out)
    rp <- as.character(raman_pid)

    .lookup <- function(tbl) {
      if (is.null(tbl) || nrow(tbl) == 0 || !"raman_particle_id" %in% names(tbl)) return(NULL)
      hit <- tbl[as.character(tbl$raman_particle_id) == rp, , drop = FALSE]
      if (nrow(hit) == 0) NULL else hit[1, , drop = FALSE]
    }

    r_pe <- .lookup(matched_pe)
    if (!is.null(r_pe)) { out$ftir <- r_pe; out$raman <- r_pe }
    r_br <- .lookup(matched_bruker)
    if (!is.null(r_br)) { out$ftir_bruker <- r_br; if (is.null(out$raman)) out$raman <- r_br }
    r_ld <- .lookup(matched_ldir)
    if (!is.null(r_ld)) { out$ldir <- r_ld; if (is.null(out$raman)) out$raman <- r_ld }

    out$count <- sum(!vapply(out, is.null, logical(1)))
    out
  }

  find_nearest_overlay_particle <- function(px, py, snap_dist,
                                            inst, rel, matched, ldir_m, dfs,
                                            bruker_m = data.frame()) {
    best_dist   <- Inf
    best_row    <- NULL
    best_source <- NULL

    # When inst/rel are NULL the caller wants to search everything (click handler)
    search_all <- is.null(inst)

    chk <- function(rel_val, ...) search_all || (rel_val %in% rel && all(c(...) %in% inst))

    # Single-instrument mode: the plot draws ALL particles of the one selected
    # instrument, ignoring the matched/unmatched relationship filter (see the
    # length(inst)==1 branch in output$overlay_plot). The pairing/unmatched
    # branches below can't reach those points — a matched particle's branch
    # needs its partner instrument selected, and unmatched only fires when
    # "unmatched" is in rel — so hovering a Raman-only view found nothing.
    # Mirror the plot here: search every particle of that instrument.
    if (!search_all && length(inst) == 1) {
      one     <- inst[[1]]
      one_key <- switch(one, ftir_pe = "ftir", ftir_bruker = "ftir_bruker",
                              raman = "raman", ldir = "ldir", NA_character_)
      one_src <- switch(one, ftir_pe = "single_ftir", ftir_bruker = "single_ftir_bruker",
                              raman = "single_raman", ldir = "single_ldir", NA_character_)
      one_mat <- switch(one, ftir_pe = "overlay_ftir_material",
                              ftir_bruker = "overlay_ftir_bruker_material",
                              raman = "overlay_raman_material",
                              ldir = "overlay_ldir_material", NA_character_)
      df_s <- if (!is.na(one_key)) dfs[[one_key]] else NULL
      if (!is.null(df_s) && nrow(df_s) > 0) {
        # Respect the same material filter the plot applies.
        mat <- if (!is.na(one_mat)) mat_keep(one_mat) else NULL
        if (!is.null(mat) && "material_family" %in% names(df_s))
          df_s <- df_s[df_s$material_family %in% mat, ]
        if (nrow(df_s) > 0) {
          d_s <- sqrt((df_s$x - px)^2 + (df_s$y - py)^2); idx_s <- which.min(d_s)
          if (length(idx_s) > 0 && d_s[idx_s] < best_dist) {
            best_dist <- d_s[idx_s]; best_row <- df_s[idx_s, , drop = FALSE]
            best_source <- one_src
          }
        }
      }
    }

    # FTIR-PE <-> Raman matched
    if (chk("matched", "ftir_pe", "raman") && !is.null(matched) && nrow(matched) > 0) {
      dist_f <- sqrt((matched$ftir_x_aligned - px)^2 + (matched$ftir_y_aligned - py)^2)
      dist_r <- sqrt((matched$raman_x_norm   - px)^2 + (matched$raman_y_norm   - py)^2)
      d <- pmin(dist_f, dist_r); idx <- which.min(d)
      if (length(idx) > 0 && d[idx] < best_dist) {
        best_dist <- d[idx]; best_row <- matched[idx, ]; best_source <- "ftir_raman"
      }
    }

    # LDIR <-> Raman matched
    if (chk("matched", "ldir", "raman") && !is.null(ldir_m) && nrow(ldir_m) > 0 &&
        "ldir_x_aligned" %in% names(ldir_m)) {
      dist_l <- sqrt((ldir_m$ldir_x_aligned - px)^2 + (ldir_m$ldir_y_aligned - py)^2)
      idx_l <- which.min(dist_l)
      if (length(idx_l) > 0 && dist_l[idx_l] < best_dist) {
        best_dist <- dist_l[idx_l]; best_row <- ldir_m[idx_l, ]; best_source <- "ldir_raman"
      }
    }

    # FTIR (Bruker) <-> Raman matched
    if (chk("matched", "ftir_bruker", "raman") && !is.null(bruker_m) && nrow(bruker_m) > 0 &&
        "ftir_x_aligned" %in% names(bruker_m)) {
      dist_fb <- sqrt((bruker_m$ftir_x_aligned - px)^2 + (bruker_m$ftir_y_aligned - py)^2)
      dist_rb <- sqrt((bruker_m$raman_x_norm   - px)^2 + (bruker_m$raman_y_norm   - py)^2)
      d <- pmin(dist_fb, dist_rb); idx <- which.min(d)
      if (length(idx) > 0 && d[idx] < best_dist) {
        best_dist <- d[idx]; best_row <- bruker_m[idx, ]; best_source <- "bruker_raman"
      }
    }

    # Unmatched FTIR-PE
    if (chk("unmatched", "ftir_pe") && !is.null(dfs$ftir) && nrow(dfs$ftir) > 0) {
      um_f <- dfs$ftir[dfs$ftir$match_status == "unmatched", ]
      if (nrow(um_f) > 0) {
        d_f <- sqrt((um_f$x - px)^2 + (um_f$y - py)^2); idx_f <- which.min(d_f)
        if (length(idx_f) > 0 && d_f[idx_f] < best_dist) {
          best_dist <- d_f[idx_f]; best_row <- um_f[idx_f, , drop = FALSE]
          best_source <- "single_ftir"
        }
      }
    }

    # Unmatched Raman
    if (chk("unmatched", "raman") && !is.null(dfs$raman) && nrow(dfs$raman) > 0) {
      um_r <- dfs$raman[dfs$raman$match_status == "unmatched", ]
      if (nrow(um_r) > 0) {
        d_r <- sqrt((um_r$x - px)^2 + (um_r$y - py)^2); idx_r <- which.min(d_r)
        if (length(idx_r) > 0 && d_r[idx_r] < best_dist) {
          best_dist <- d_r[idx_r]; best_row <- um_r[idx_r, , drop = FALSE]
          best_source <- "single_raman"
        }
      }
    }

    # Unmatched LDIR
    if (chk("unmatched", "ldir") && !is.null(dfs$ldir) && nrow(dfs$ldir) > 0) {
      um_l <- dfs$ldir[dfs$ldir$match_status == "unmatched", ]
      if (nrow(um_l) > 0) {
        d_l <- sqrt((um_l$x - px)^2 + (um_l$y - py)^2); idx_l <- which.min(d_l)
        if (length(idx_l) > 0 && d_l[idx_l] < best_dist) {
          best_dist <- d_l[idx_l]; best_row <- um_l[idx_l, , drop = FALSE]
          best_source <- "single_ldir"
        }
      }
    }

    # Unmatched FTIR Bruker
    if (chk("unmatched", "ftir_bruker") && !is.null(dfs$ftir_bruker) && nrow(dfs$ftir_bruker) > 0) {
      um_fb <- dfs$ftir_bruker[dfs$ftir_bruker$match_status == "unmatched", ]
      if (nrow(um_fb) > 0) {
        d_fb <- sqrt((um_fb$x - px)^2 + (um_fb$y - py)^2); idx_fb <- which.min(d_fb)
        if (length(idx_fb) > 0 && d_fb[idx_fb] < best_dist) {
          best_dist <- d_fb[idx_fb]; best_row <- um_fb[idx_fb, , drop = FALSE]
          best_source <- "single_ftir_bruker"
        }
      }
    }

    # ALWAYS check highlighted/selected particles regardless of layer state
    hl_specs <- list(
      list(ids = input$overlay_ftir_particles,        df = dfs$ftir,        src = "single_ftir"),
      list(ids = input$overlay_raman_particles,       df = dfs$raman,       src = "single_raman"),
      list(ids = input$overlay_ldir_particles,        df = dfs$ldir,        src = "single_ldir"),
      list(ids = input$overlay_ftir_bruker_particles, df = dfs$ftir_bruker, src = "single_ftir_bruker")
    )
    for (spec in hl_specs) {
      if (is.null(spec$ids) || length(spec$ids) == 0) next
      inst_df <- spec$df
      if (is.null(inst_df) || nrow(inst_df) == 0) next
      hl <- inst_df[inst_df$particle_id %in% spec$ids, ]
      if (nrow(hl) > 0) {
        d_hl <- sqrt((hl$x - px)^2 + (hl$y - py)^2); idx_hl <- which.min(d_hl)
        if (length(idx_hl) > 0 && d_hl[idx_hl] < best_dist) {
          best_dist <- d_hl[idx_hl]; best_row <- hl[idx_hl, , drop = FALSE]
          best_source <- spec$src
        }
      }
    }

    if (best_dist <= snap_dist && !is.null(best_row)) {
      list(row = best_row, source = best_source, dist = best_dist)
    } else NULL
  }

  # Overlay: sticky hover — update last_hover$overlay only when a new match is found.
  # Checks active instruments/relationships AND highlighted/selected particles regardless of state.
  observeEvent(input$overlay_hover, {
    hover <- input$overlay_hover
    if (is.null(hover)) return()

    inst <- input$overlay_instruments
    rel  <- input$overlay_relationships
    matched <- overlay_matched()
    ldir_m <- overlay_ldir_matched()
    bruker_m <- overlay_bruker_matched()
    dfs <- list(ftir = ftir_df_full(), raman = raman_df_full(), ldir = ldir_df_full(),
                ftir_bruker = ftir_bruker_df_full())

    vis <- if (!is.null(zoom$overlay)) zoom$overlay
           else compute_bounds(dfs$ftir, dfs$raman, dfs$ldir)
    snap_dist <- max(diff(vis$x), diff(vis$y), 500) * 0.05

    result <- find_nearest_overlay_particle(hover$x, hover$y, snap_dist,
                                            inst, rel, matched, ldir_m, dfs, bruker_m)
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
    bruker_m <- overlay_bruker_matched()
    dfs <- list(ftir = ftir_df_full(), raman = raman_df_full(), ldir = ldir_df_full(),
                ftir_bruker = ftir_bruker_df_full())

    vis <- if (!is.null(zoom$overlay)) zoom$overlay
           else compute_bounds(dfs$ftir, dfs$raman, dfs$ldir)
    snap_dist <- max(diff(vis$x), diff(vis$y), 500) * 0.05

    # For click, pass NULL inst/rel to search ALL instruments
    result <- find_nearest_overlay_particle(click$x, click$y, snap_dist,
                                            NULL, NULL, matched, ldir_m, dfs, bruker_m)
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
      } else if (src == "bruker_raman") {
        if (!is.null(row$ftir_particle_id) && !is.null(dfs$ftir_bruker))
          .toggle_dropdown("overlay_ftir_bruker_particles", row$ftir_particle_id,
                           sort(dfs$ftir_bruker$particle_id))
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
      } else if (src == "single_ftir_bruker" && !is.null(row$particle_id)) {
        .toggle_dropdown("overlay_ftir_bruker_particles", row$particle_id,
                         sort(dfs$ftir_bruker$particle_id))
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

    # Matched particle: gather every instrument partner via the shared Raman id.
    # run_data() keys: matched = FTIR(PE)<->Raman, matched_ftir_bruker =
    # FTIR(Bruker)<->Raman, ldir_raman_matched = LDIR<->Raman. All join on
    # raman_particle_id, so one Raman particle can carry up to three partners.
    rd <- run_data_gated()
    raman_pid <- if (!is.null(row$raman_particle_id)) row$raman_particle_id else NULL
    # Only surface a genuine (within-gate) LDIR partner — a forced over-gate
    # pairing is reported as unmatched everywhere else, so it must not appear
    # here as a partner. Honours the live gate slider.
    partners <- gather_partners_by_raman(raman_pid, rd$matched,
                                         ldir_genuine_pairs(rd$ldir_raman_matched),
                                         rd$matched_ftir_bruker)

    # Safe accessors: return NULL/em-dash for missing columns rather than erroring.
    .g <- function(frame, col) {
      if (is.null(frame) || !col %in% names(frame)) return(NULL)
      frame[[col]][1]
    }
    .fmt_num <- function(v, digits = 1, prefix = "", suffix = "") {
      if (is.null(v) || length(v) == 0 || (is.atomic(v) && all(is.na(v)))) return("\u2014")
      if (is.numeric(v)) return(paste0(prefix, round(v, digits), suffix))
      as.character(v)
    }
    .fmt_pos <- function(frame, xcol, ycol) {
      xv <- .g(frame, xcol); yv <- .g(frame, ycol)
      if (is.null(xv) || is.null(yv) || is.na(xv) || is.na(yv)) return("\u2014")
      paste0("(", round(xv, 1), ", ", round(yv, 1), ")")
    }

    # One column spec per present partner, in a fixed left-to-right order.
    specs <- list()
    if (!is.null(partners$ftir))
      specs[[length(specs) + 1]] <- list(hdr = "FTIR", f = partners$ftir, p = "ftir_",
        qpfx = "AAU ", qd = 3, xcol = "ftir_x_aligned", ycol = "ftir_y_aligned", is_raman = FALSE)
    if (!is.null(partners$ftir_bruker))
      specs[[length(specs) + 1]] <- list(hdr = "FTIR (Bruker)", f = partners$ftir_bruker, p = "ftir_",
        qpfx = "AAU ", qd = 3, xcol = "ftir_x_aligned", ycol = "ftir_y_aligned", is_raman = FALSE)
    if (!is.null(partners$raman))
      specs[[length(specs) + 1]] <- list(hdr = "Raman", f = partners$raman, p = "raman_",
        qpfx = "HQI ", qd = 2, xcol = "raman_x_um", ycol = "raman_y_um", is_raman = TRUE)
    if (!is.null(partners$ldir))
      specs[[length(specs) + 1]] <- list(hdr = "LDIR", f = partners$ldir, p = "ldir_",
        qpfx = "", qd = 3, xcol = "ldir_x_aligned", ycol = "ldir_y_aligned", is_raman = FALSE)

    # Fallback: no partner rows resolved (e.g. stale ids) \u2014 show minimal detail.
    if (length(specs) == 0) {
      return(tags$table(class = "hover-tbl",
        tags$tr(tags$th("Field"), tags$th("Value")),
        make_detail_row("Raman particle", if (!is.null(raman_pid)) raman_pid else "\u2014"),
        make_detail_row("Note", "No match record found for this particle.")))
    }

    .cell <- function(x) tags$td(x)
    mkrow <- function(label, cellfn)
      do.call(tags$tr, c(list(tags$td(tags$b(label))),
                         lapply(specs, function(s) .cell(cellfn(s)))))

    tags$table(class = "hover-tbl",
      do.call(tags$tr, c(list(tags$th("")), lapply(specs, function(s) tags$th(s$hdr)))),
      mkrow("Particle ID", function(s) {
        v <- .g(s$f, paste0(s$p, "particle_id")); if (is.null(v)) "\u2014" else as.character(v) }),
      mkrow("Material", function(s) {
        v <- .g(s$f, paste0(s$p, "material")); if (is.null(v)) "\u2014" else as.character(v) }),
      mkrow("Quality", function(s)
        .fmt_num(.g(s$f, paste0(s$p, "quality")), s$qd, prefix = s$qpfx)),
      mkrow("Feret Max", function(s)
        .fmt_num(.g(s$f, paste0(s$p, "feret_max_um")), 1, suffix = " \u00b5m")),
      mkrow("Area", function(s)
        .fmt_num(.g(s$f, paste0(s$p, "area_um2")), 1, suffix = " \u00b5m\u00b2")),
      mkrow("Position", function(s) .fmt_pos(s$f, s$xcol, s$ycol)),
      mkrow("Match Dist.", function(s)
        if (isTRUE(s$is_raman)) "\u2014" else .fmt_num(.g(s$f, "match_distance"), 1, suffix = " \u00b5m"))
    )
  })

  # ==================================================================
  # MULTI-RUN (reproducibility) TAB
  # ==================================================================
  repro_dir_r <- reactiveVal(NULL)

  # Scan the base folder for reproducibility runs: any subfolder (or the base
  # itself) containing reproducibility_points.csv. Newest first.
  # Returns a list of records: list(dir, instrument, n_runs) — newest first.
  repro_scan <- reactive({
    input$repro_refresh                       # manual rescan trigger
    base <- input$repro_base
    if (is.null(base) || !nzchar(base)) base <- "output/reproducibility"
    has_pts <- function(d) file.exists(file.path(d, "reproducibility_points.csv"))
    dirs <- character(0)
    if (has_pts(base)) dirs <- c(dirs, base)
    subs <- tryCatch(list.dirs(base, recursive = FALSE, full.names = TRUE),
                     error = function(e) character(0))
    subs <- subs[vapply(subs, has_pts, logical(1))]
    dirs <- unique(c(dirs, sort(subs, decreasing = TRUE)))
    lapply(dirs, function(d) {
      m <- tryCatch(
        read.csv(file.path(d, "reproducibility_meta.csv"),
                 stringsAsFactors = FALSE, nrows = 1),
        error = function(e) NULL
      )
      list(
        dir        = d,
        instrument = if (!is.null(m) && "instrument" %in% names(m)) tolower(m$instrument[1]) else NA_character_,
        n_runs     = if (!is.null(m) && "n_runs"     %in% names(m)) as.integer(m$n_runs[1])  else NA_integer_
      )
    })
  })

  .repro_instrument_group <- function(inst) {
    if (is.na(inst)) return("Unknown")
    switch(inst,
      ldir        = "LDIR - Agilent 8700",
      raman       = "Raman - WITec",
      ftir_perkin = "FTIR - PerkinElmer Spotlight",
      ftir_bruker = "FTIR - Bruker OPUS / ALPHA",
      paste0("Other (", inst, ")")
    )
  }

  # Render the Run selector as a grouped (optgroup) selectInput. Using renderUI
  # rather than updateSelectInput because renderUI reliably supports named-list
  # choices (optgroups) across all Shiny versions.
  output$repro_run_select_ui <- renderUI({
    items <- repro_scan()
    if (length(items) == 0)
      return(selectInput("repro_run_select", "Run",
                         choices = character(0), width = "100%"))
    groups <- list()
    for (item in items) {
      grp <- .repro_instrument_group(item$instrument)
      lbl <- basename(item$dir)
      if (!is.na(item$n_runs)) lbl <- paste0(lbl, " (", item$n_runs, " runs)")
      groups[[grp]] <- c(groups[[grp]], setNames(item$dir, lbl))
    }
    # Preserve current selection if it still exists after a rescan.
    all_dirs <- vapply(items, `[[`, character(1), "dir")
    cur <- isolate(input$repro_run_select)
    sel <- if (!is.null(cur) && nzchar(cur) && cur %in% all_dirs) cur else items[[1]]$dir
    selectInput("repro_run_select", "Run",
                choices = groups, selected = sel, width = "100%")
  })

  observeEvent(input$repro_run_select, {
    if (!is.null(input$repro_run_select) && nzchar(input$repro_run_select))
      repro_dir_r(input$repro_run_select)
  })

  # Optional native folder picker — only when shinyFiles is installed; the
  # scan dropdown works regardless.
  output$repro_browse_ui <- renderUI({
    if (!requireNamespace("shinyFiles", quietly = TRUE)) return(NULL)
    shinyFiles::shinyDirButton("repro_dir_btn", "Browse…",
                               "Choose a reproducibility folder", class = "btn-sm")
  })
  if (requireNamespace("shinyFiles", quietly = TRUE)) {
    .repro_roots <- c(project = normalizePath(".", mustWork = FALSE),
                      home = path.expand("~"))
    shinyFiles::shinyDirChoose(input, "repro_dir_btn", roots = .repro_roots,
                               session = session)
    observeEvent(input$repro_dir_btn, {
      p <- tryCatch(shinyFiles::parseDirPath(.repro_roots, input$repro_dir_btn),
                    error = function(e) character(0))
      if (length(p) == 1 && nzchar(p)) updateTextInput(session, "repro_base", value = p)
    })
  }

  # Load a tools/reproducibility.R output folder (points + meta + summary).
  repro_data <- reactive({
    dir <- repro_dir_r()
    if (is.null(dir) || !nzchar(dir)) return(NULL)
    pts_path <- file.path(dir, "reproducibility_points.csv")
    if (!file.exists(pts_path)) return(NULL)
    pts <- tryCatch(read.csv(pts_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(pts) || nrow(pts) == 0) return(NULL)
    rd <- function(f) tryCatch(read.csv(file.path(dir, f), stringsAsFactors = FALSE),
                               error = function(e) NULL)
    list(points = pts, meta = rd("reproducibility_meta.csv"),
         summary = rd("reproducibility_summary.csv"), dir = dir)
  })

  observeEvent(repro_data(), {
    d <- repro_data()
    fams <- if (!is.null(d) && "material_family" %in% names(d$points))
              sort(unique(d$points$material_family)) else character(0)
    updateSelectizeInput(session, "repro_material",
                         choices = c("All", fams), selected = "All")

    run_ids <- if (!is.null(d) && "run" %in% names(d$points))
                 sort(unique(d$points$run)) else integer(0)
    # paste("Run", integer(0)) returns "Run " (length 1, NOT character(0)),
    # so setNames() below would fail its length check unless run_ids is
    # empty is handled explicitly first.
    run_choices <- if (length(run_ids) > 0)
                     setNames(as.character(run_ids), paste("Run", run_ids))
                   else character(0)
    updateCheckboxGroupInput(session, "repro_run_visibility",
                             choices = run_choices, selected = run_choices)
  }, ignoreNULL = FALSE)

  # Filtered/annotated points shared by the plot AND the hover lookup, so
  # hovering always matches what's actually drawn (mirrors ftir_filtered() /
  # raman_filtered() etc. used the same way for the single-instrument tabs).
  repro_filtered_pts <- reactive({
    d <- repro_data()
    if (is.null(d)) return(NULL)
    pts    <- d$points
    n_runs <- if (!is.null(d$meta) && "n_runs" %in% names(d$meta)) d$meta$n_runs[1]
              else max(pts$run, na.rm = TRUE)

    # Identification-quality gate. Unlike run visibility and the material
    # filter below, this is an INCLUSION criterion rather than a display
    # choice, so it is applied first and the per-consensus flags are
    # recomputed from what survives. Without that recompute a particle whose
    # partner in another run fell below the threshold would still be ringed as
    # reproduced in every run — the rings would describe the unfiltered
    # analysis while the plot showed the filtered one. Rows with NA quality are
    # kept (an instrument or run that never reported one).
    qr <- input$repro_quality_range
    if (!is.null(qr) && length(qr) == 2 && "quality" %in% names(pts)) {
      keep <- is.na(pts$quality) | (pts$quality >= qr[1] & pts$quality <= qr[2])
      if (!all(keep)) {
        pts <- pts[keep, , drop = FALSE]
        if (nrow(pts) == 0) { attr(pts, "empty_reason") <- "quality"; return(pts) }
        cid <- as.character(pts$consensus_id)
        det <- tapply(pts$run, pts$consensus_id, function(r) length(unique(r)))
        pts$n_runs_detected <- as.integer(det[cid])
        # Same rule as the upstream consensus builder: concordant when the
        # surviving runs agree on at most one non-missing family.
        con <- tapply(pts$material_family, pts$consensus_id,
                      function(f) length(unique(f[!is.na(f)])) <= 1L)
        pts$material_concordant <- as.logical(con[cid])
      }
    }

    # Run visibility: same "empty selection -> empty result" convention as the
    # instrument tabs' match-status checkboxGroupInput (eff_match_filter());
    # NULL (widget not yet populated on initial load) is treated as "no
    # filter yet" so data isn't hidden before repro_data()'s observer runs.
    vis <- input$repro_run_visibility
    if (!is.null(vis)) {
      visible_runs <- suppressWarnings(as.integer(vis))
      pts <- pts[pts$run %in% visible_runs, , drop = FALSE]
      if (nrow(pts) == 0) { attr(pts, "empty_reason") <- "runs"; return(pts) }
    }

    mats <- input$repro_material
    if (!is.null(mats) && !("All" %in% mats) && "material_family" %in% names(pts))
      pts <- pts[pts$material_family %in% mats, , drop = FALSE]
    if (nrow(pts) == 0) { attr(pts, "empty_reason") <- "material"; return(pts) }

    pts$status <- ifelse(!as.logical(pts$material_concordant), "discordant",
                  ifelse(pts$n_runs_detected < n_runs, "missing", "consensus"))
    pts$run_lbl <- paste("Run", pts$run)

    if (isTRUE(input$repro_only_nonrepro)) {
      keep <- unique(pts$consensus_id[pts$status != "consensus"])
      pts  <- pts[pts$consensus_id %in% keep, , drop = FALSE]
      if (nrow(pts) == 0) attr(pts, "empty_reason") <- "nonrepro"
    }

    # Coord match cost filter (LDIR only — coord_match_cost is NA for FTIR/Raman).
    # Keeps rows with NA cost (non-LDIR instruments) and within-range rows.
    cost_r <- input$repro_coord_cost_range
    if (!is.null(cost_r) && length(cost_r) == 2 &&
        "coord_match_cost" %in% names(pts)) {
      pts <- pts[is.na(pts$coord_match_cost) |
                 (pts$coord_match_cost >= cost_r[1] &
                  pts$coord_match_cost <= cost_r[2]), , drop = FALSE]
    }

    pts
  })

  # Current background-image rotation (degrees; 0 if unset/invalid). Shared by
  # the plot and the hover lookup so both agree on which orientation is shown.
  # Mirrors the single LDIR tab's view-rotation control (rotate_raster_view/
  # rotate_xy_view already normalise negative values via %% 360, so -90 and
  # 270 are equivalent).
  repro_rotation_deg <- reactive({
    rot <- suppressWarnings(as.integer(input$repro_img_rotation))
    if (is.na(rot)) 0L else rot
  })

  # Rotate a points frame's aligned coordinates in place (no-op at 0 deg).
  .repro_rotate_view <- function(pts, rot) {
    if (is.null(pts) || nrow(pts) == 0 || is.na(rot) || rot == 0L) return(pts)
    rc <- rotate_xy_view(pts$x_aligned, pts$y_aligned, rot)
    pts$x_aligned <- rc$x
    pts$y_aligned <- rc$y
    pts
  }

  # Points in the CURRENTLY DISPLAYED (rotated) frame — what the plot actually
  # draws once a background-image rotation is applied. Used by the hover
  # lookup so it matches rotated positions; the plot itself rotates its own
  # copy after computing the (calibration-based) image placement, which must
  # happen in the native/unrotated frame — see output$repro_plot.
  repro_view_pts <- reactive({
    .repro_rotate_view(repro_filtered_pts(), repro_rotation_deg())
  })

  output$repro_plot <- renderPlot({
    d <- repro_data()
    .msg <- function(txt, col = "grey40", size = 5)
      ggplot() + annotate("text", x = 0, y = 0, label = txt, size = size, colour = col) +
        theme_void()
    if (is.null(d))
      return(.msg(paste0("No reproducibility output found.\nRun tools/reproducibility.R, ",
                         "then set the folder and press Load.")))
    pts <- repro_filtered_pts()
    if (is.null(pts) || nrow(pts) == 0) {
      reason <- attr(pts, "empty_reason")
      if (identical(reason, "nonrepro"))
        return(.msg("No non-reproducible particles \u2014 every particle is consistent.",
                    col = "#2ca02c", size = 6))
      if (identical(reason, "runs"))
        return(.msg("No runs selected \u2014 check at least one run above."))
      if (identical(reason, "quality"))
        return(.msg("No particles pass the quality threshold \u2014 widen the slider."))
      return(.msg("No particles match the material filter."))
    }

    # Background image first (it frames the view). Uploaded image takes
    # precedence over the tool-copied backdrop; an optional 90° rotation lets
    # the user orient the image to the particle layout.
    img_info <- NULL
    if (isTRUE(input$repro_show_image)) {
      bg_path <- NULL
      up <- input$repro_bg_upload
      if (!is.null(up) && !is.null(up$datapath) && file.exists(up$datapath)) {
        bg_path <- up$datapath
      } else if (!is.null(d$meta) && "bg_image" %in% names(d$meta)) {
        bgf <- d$meta$bg_image[1]
        if (!is.na(bgf) && nzchar(bgf) && file.exists(file.path(d$dir, bgf)))
          bg_path <- file.path(d$dir, bgf)
      }
      if (!is.null(bg_path)) {
        raw <- tryCatch(load_image_raster(bg_path), error = function(e) NULL)
        if (!is.null(raw)) {
          # NOTE: rotation is applied LAST, to the finished (native-frame)
          # placement — not here. Calibration-based placement (WITec extent,
          # LDIR scan-circle, um-per-px scale) is defined in the image's
          # native/unrotated pixel frame, so rotating the raster before
          # computing `base` would place a rotated (width<->height swapped)
          # raster into a box sized for the unrotated image, and would
          # silently drift out of registration with the (also unrotated)
          # particle points. See the co-rotation step below, which mirrors
          # the single LDIR tab's ldir_view_tf() handling exactly.
          ox <- if (isTRUE(is.finite(input$repro_img_offset_x))) input$repro_img_offset_x else 0
          oy <- if (isTRUE(is.finite(input$repro_img_offset_y))) input$repro_img_offset_y else 0
          w_um <- input$repro_img_width_um
          h_um <- input$repro_img_height_um
          # Placement anchors: the FULL point set, never the filtered view.
          # The backdrop is a property of the dataset, not of what is currently
          # drawn. Anchoring it to `pts` made the image resize and shift every
          # time the quality slider, run visibility, material filter or "only
          # non-reproducible" changed the point cloud — with "only
          # non-reproducible" on, a handful of surviving particles crushed the
          # whole micrograph into their bounding box.
          place_x <- d$points$x_aligned[is.finite(d$points$x_aligned)]
          place_y <- d$points$y_aligned[is.finite(d$points$y_aligned)]
          if (length(place_x) == 0 || length(place_y) == 0) {
            place_x <- pts$x_aligned; place_y <- pts$y_aligned
          }
          instrument <- if (!is.null(d$meta) && "instrument" %in% names(d$meta))
                          d$meta$instrument[1] else NA_character_
          base <- NULL
          if (!is.null(w_um) && is.finite(w_um) && w_um > 0) {
            # Manual override: physical width (+ optional height) centred on the
            # particle mean. Overrides the metadata placement below.
            half_w <- w_um / 2
            half_h <- if (!is.null(h_um) && is.finite(h_um) && h_um > 0)
                        h_um / 2 else (w_um / (ncol(raw) / nrow(raw))) / 2
            cx0 <- mean(place_x, na.rm = TRUE)
            cy0 <- mean(place_y, na.rm = TRUE)
            base <- list(xmin = cx0 - half_w, xmax = cx0 + half_w,
                         ymin = cy0 - half_h, ymax = cy0 + half_h)
          } else {
            # Default: reproduce the single-instrument tab's exact placement from
            # the metadata the tool recorded (Raman = WITec extent, FTIR/Bruker =
            # particle extent, LDIR = scan-circle). NULL when metadata is absent.
            base <- place_image_multirun(instrument, d$meta,
                                         place_x, place_y,
                                         raw = raw, bg_path = bg_path)
          }
          # Last resort: aspect-preserving fit to the particle extent.
          if (is.null(base))
            base <- compute_image_bounds(raw, place_x, place_y,
                                         padding_um = 200)
          b <- list(xmin = base$xmin + ox, xmax = base$xmax + ox,
                    ymin = base$ymin + oy, ymax = base$ymax + oy)
          img_info <- c(list(raster = raw), b)
        }
      }
    }

    # Co-rotate the finished (native-frame) image AND the particle points
    # together — exactly mirroring the single LDIR tab's view-rotation
    # handling (rotate_extent_view/rotate_raster_view/rotate_xy_view), so the
    # two can never drift out of registration the way rotating only the
    # raster (the previous behaviour) did.
    rot <- repro_rotation_deg()
    if (rot != 0L) {
      if (!is.null(img_info)) {
        ext <- rotate_extent_view(img_info, rot)
        img_info <- list(raster = rotate_raster_view(img_info$raster, rot),
                         xmin = ext$xmin, xmax = ext$xmax,
                         ymin = ext$ymin, ymax = ext$ymax)
      }
      pts <- .repro_rotate_view(pts, rot)
    }

    # Frame on the IMAGE extent when a background is placed (matches the single
    # Raman tab: shows the whole image, not just the particle sub-region), else
    # on the particle extent. A brushed zoom overrides both.
    bounds <- if (!is.null(zoom$repro)) zoom$repro else if (!is.null(img_info)) {
      pad <- 200
      list(x = c(img_info$xmin - pad, img_info$xmax + pad),
           y = c(img_info$ymin - pad, img_info$ymax + pad))
    } else {
      compute_bounds(data.frame(x = pts$x_aligned, y = pts$y_aligned))
    }
    bounds <- sanitize_bounds(bounds)

    run_levels <- paste("Run", sort(unique(pts$run)))
    run_pal <- .REPRO_RUN_PALETTE[run_levels]

    show_miss <- isTRUE(input$repro_show_missing_ring)
    show_disc <- isTRUE(input$repro_show_discordant_ring)
    ring_subtitle <- paste(c(
      if (show_miss) "Cyan ring = missing in \u22651 run",
      if (show_disc) "Pink ring = material disagreement"
    ), collapse = "  \u00b7  ")
    if (!nzchar(ring_subtitle))
      ring_subtitle <- "Diagnostic rings hidden \u2014 see checkboxes at left"

    p <- ggplot() +
      scale_x_continuous(breaks = breaks_adaptive(bounds$x)) +
      scale_y_continuous(breaks = breaks_adaptive(bounds$y)) +
      coord_fixed(xlim = bounds$x, ylim = bounds$y, expand = FALSE) +
      labs(title = "Multi-Run reproducibility overlay",
           subtitle = ring_subtitle,
           x = "X (\u00b5m)", y = "Y (\u00b5m)") +
      theme_minimal(base_size = 15) +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "grey98", colour = NA),
            panel.grid = element_line(colour = "grey90"), legend.position = "right")
    if (!is.null(img_info)) p <- add_image_bg(p, img_info)

    if (isTRUE(input$repro_show_lines)) {
      multi <- pts[pts$consensus_id %in% pts$consensus_id[duplicated(pts$consensus_id)], ]
      if (nrow(multi) > 0)
        p <- p + geom_path(data = multi[order(multi$consensus_id, multi$run), ],
                           aes(x = x_aligned, y = y_aligned, group = consensus_id),
                           colour = "grey50", alpha = 0.3, linewidth = 0.4)
    }

    base_alpha <- if (isTRUE(input$repro_only_nonrepro)) 0.85 else 0.5
    p <- p + geom_point(data = pts,
                        aes(x = x_aligned, y = y_aligned, colour = run_lbl, size = feret_max_um),
                        alpha = base_alpha) +
      scale_colour_manual(name = "Run", values = run_pal) +
      scale_size_continuous(name = "Feret (\u00b5m)", range = c(1.5, 7),
                            limits = safe_size_limits(pts$feret_max_um))

    if (show_miss) {
      miss <- pts[pts$status == "missing", ]
      if (nrow(miss) > 0)
        p <- p + geom_point(data = miss, aes(x = x_aligned, y = y_aligned),
                            shape = 21, size = 5, stroke = 1.2, fill = NA, colour = "#22d3ee")
    }
    if (show_disc) {
      disc <- pts[pts$status == "discordant", ]
      if (nrow(disc) > 0)
        p <- p + geom_point(data = disc, aes(x = x_aligned, y = y_aligned),
                            shape = 21, size = 6, stroke = 1.6, fill = NA, colour = "#ff1493")
    }
    p
  })

  # Nearest-particle lookup in the aligned (plot) frame, mirroring
  # raman_hover/ldir_hover etc. Matches against the SAME filtered points the
  # plot actually draws (repro_filtered_pts()), so hovering never highlights
  # a point that isn't visible.
  observeEvent(input$repro_hover, {
    hover <- input$repro_hover
    if (is.null(hover)) return()
    pts <- repro_view_pts()
    if (is.null(pts) || nrow(pts) == 0) return()
    d <- sqrt((pts$x_aligned - hover$x)^2 + (pts$y_aligned - hover$y)^2)
    idx <- which.min(d)
    threshold <- max(diff(range(pts$x_aligned, na.rm = TRUE)),
                     diff(range(pts$y_aligned, na.rm = TRUE)), 500) * 0.05
    if (d[idx] <= threshold) {
      cid <- pts$consensus_id[idx]
      last_hover$repro <- pts[pts$consensus_id == cid, , drop = FALSE]
    }
  })

  # Particle Details (hover): one row per run the hovered physical particle
  # was detected in — the cross-run comparison this tab exists to show.
  output$repro_hover_info <- renderUI({
    rows <- last_hover$repro
    if (is.null(rows) || nrow(rows) == 0)
      return(tags$p(class = "text-muted",
                    "Hover over a particle to see it across runs"))
    rows <- rows[order(rows$run), ]

    # LDIR only: how confidently this particle's material/quality (an Excel
    # row) was matched to its image-detected position, purely by size
    # similarity (no position is involved in that join - see
    # docs/multirun_image_placement_plan.md). High cost => the material shown
    # here may belong to a different blob than the one actually plotted,
    # which can surface as a false "material disagreement" between runs.
    # Hidden for instruments with native coordinates (always NA there).
    show_cost      <- "coord_match_cost" %in% names(rows) && any(is.finite(rows$coord_match_cost))
    has_area       <- "area_um2"        %in% names(rows) && any(is.finite(rows$area_um2))
    has_major      <- "major_um"        %in% names(rows) && any(is.finite(rows$major_um))
    has_img_area   <- "image_area_um2"  %in% names(rows) && any(is.finite(rows$image_area_um2))
    has_img_feret  <- "image_feret_um"  %in% names(rows) && any(is.finite(rows$image_feret_um))
    ncol <- 5L + has_area + has_major + has_img_area + has_img_feret + show_cost

    # Column header helper: label + small muted origin tag
    .th_orig <- function(label, origin) {
      tags$th(label,
              tags$span(class = "text-muted",
                        style = "font-weight:normal; font-size:0.78em; margin-left:3px;",
                        paste0("(", origin, ")")))
    }

    div(style = "overflow-x: auto; width: 100%;",
    tags$table(class = "hover-tbl", style = "white-space: nowrap; min-width: 100%;",
      tags$tr(
        tags$th("Run"), tags$th("Particle ID"), tags$th("Material"),
        .th_orig("Feret Max",      "machine"),
        if (has_img_feret) .th_orig("Feret Max", "image"),
        if (has_major)     .th_orig("Diameter",  "machine"),
        if (has_area)      .th_orig("Area",       "machine"),
        if (has_img_area)  .th_orig("Area",       "image"),
        tags$th("Position (aligned)"),
        if (show_cost) tags$th("Coord Match Cost")),
      lapply(seq_len(nrow(rows)), function(i) {
        r <- rows[i, ]
        tags$tr(
          tags$td(tags$b(paste("Run", r$run))),
          tags$td(r$particle_id),
          tags$td(r$material),
          tags$td(paste0(round(r$feret_max_um, 1), " µm")),
          if (has_img_feret) tags$td(
            if (is.finite(r$image_feret_um)) paste0(round(r$image_feret_um, 1), " µm") else "—"),
          if (has_major) tags$td(
            if (is.finite(r$major_um)) paste0(round(r$major_um, 1), " µm") else "—"),
          if (has_area) tags$td(
            if (is.finite(r$area_um2)) paste0(round(r$area_um2, 0), " µm²") else "—"),
          if (has_img_area) tags$td(
            if (is.finite(r$image_area_um2)) paste0(round(r$image_area_um2, 0), " µm²") else "—"),
          tags$td(paste0("(", round(r$x_aligned, 1), ", ", round(r$y_aligned, 1), ")")),
          if (show_cost) tags$td(
            if (is.finite(r$coord_match_cost)) round(r$coord_match_cost, 3) else "—")
        )
      }),
      tags$tr(tags$td(colspan = as.character(ncol),
                      style = "border-top: 2px solid #adb5bd; padding-top: 6px;",
        tags$span(class = "text-muted",
          paste0(nrow(rows), " of ", rows$n_runs_detected[1], " run(s) detected",
                 if (!isTRUE(rows$material_concordant[1]))
                   " — material mismatch across runs" else ""))
      ))
    ))   # close tags$table and wrapping div
  })

  # LDIR only: coord match cost slider (image↔Excel join confidence).
  # Auto-sized to the actual range in the loaded dataset; returns NULL for
  # FTIR/Raman runs where coord_match_cost is always NA.
  # Identification-quality gate. The scale is instrument-dependent (AAU 0-1 for
  # FTIR/LDIR, HQI 0-100 for Raman), so the range is taken from the data rather
  # than assumed. Rendered only when the loaded run actually carries quality —
  # reproducibility outputs written before that column existed simply do not
  # get the slider (re-run tools/reproducibility.R to gain it).
  output$repro_quality_ui <- renderUI({
    d <- repro_data()
    if (is.null(d)) return(NULL)
    pts <- d$points
    if (is.null(pts) || !"quality" %in% names(pts)) return(NULL)
    q <- pts$quality[is.finite(pts$quality)]
    if (length(q) == 0) return(NULL)
    lo <- floor(min(q) * 100) / 100
    hi <- ceiling(max(q) * 100) / 100
    if (!(hi > lo)) return(NULL)
    lbl <- switch(as.character(d$meta$instrument[1]),
                  raman = "HQI", ldir = "Quality (AAU)",
                  ftir_perkin = , ftir_bruker = "AAU Quality", "Quality")
    step <- if (hi - lo > 10) 1 else 0.01
    tagList(
      sliderInput("repro_quality_range", paste0(lbl, " — applies to every run"),
                  min = lo, max = hi, value = c(lo, hi), step = step),
      helpText(style = "font-size:11px; margin-top:-8px;",
               "One threshold for all runs — a per-run threshold would make the ",
               "runs non-comparable. Particles dropped here stop counting as ",
               "detected, so the rings track the filtered set.")
    )
  })

  output$repro_coord_cost_ui <- renderUI({
    d <- repro_data()
    if (is.null(d)) return(NULL)
    pts <- d$points
    if (is.null(pts) || !"coord_match_cost" %in% names(pts)) return(NULL)
    finite_costs <- pts$coord_match_cost[is.finite(pts$coord_match_cost)]
    if (length(finite_costs) == 0) return(NULL)
    cost_max <- ceiling(max(finite_costs) * 20) / 20
    sliderInput("repro_coord_cost_range",
                "Coord Match Cost (image↔Excel)",
                min = 0, max = cost_max, value = c(0, cost_max), step = 0.05)
  })

  # LDIR only: each run's own detected scan circle, side by side. Lets you
  # check whether a scan-circle detection drift between repeat runs (rather
  # than a genuine reproducibility failure) explains missing/discordant
  # particles - see docs/multirun_image_placement_plan.md. NULL (nothing
  # rendered) for any other instrument or when no per-run circle was
  # recorded (outputs from before this was added).
  output$repro_ldir_circle_ui <- renderUI({
    d <- repro_data()
    if (is.null(d) || is.null(d$meta) ||
        !identical(d$meta$instrument[1], "ldir")) return(NULL)
    meta <- d$meta
    n_runs <- if ("n_runs" %in% names(meta)) meta$n_runs[1] else 0
    run_ids <- seq_len(n_runs)
    has_any <- any(vapply(run_ids, function(i)
      paste0("ldir_run", i, "_cx_px") %in% names(meta), logical(1)))
    if (!has_any) return(NULL)

    g <- function(i, field, default = NA) {
      col <- paste0("ldir_run", i, "_", field)
      if (col %in% names(meta)) meta[[col]][1] else default
    }
    fmt_num <- function(v, digits = 1)
      if (is.null(v) || length(v) == 0 || is.na(v)) "—" else round(as.numeric(v), digits)

    rows <- lapply(run_ids, function(i) {
      tags$tr(
        tags$td(tags$b(paste("Run", i))),
        tags$td(paste0("(", fmt_num(g(i, "cx_px")), ", ", fmt_num(g(i, "cy_px")), ")")),
        tags$td(fmt_num(g(i, "radius_px"))),
        tags$td(fmt_num(g(i, "scale_um_per_px"), 4)),
        tags$td(as.character(g(i, "export_type", "—"))),
        tags$td(as.character(g(i, "method", "—"))),
        tags$td(if (isTRUE(g(i, "detected", FALSE))) "Yes" else "No")
      )
    })

    tagList(
      hr(),
      div(class = "info-box",
        h5("LDIR Scan-Circle Calibration per Run"),
        p(class = "text-muted",
          "Each run's own detected scan circle. A meaningfully different ",
          "center/radius between runs means their coordinate frames were ",
          "calibrated differently before alignment — a plausible cause of ",
          "spurious mismatches independent of true particle reproducibility. ",
          "“mosaic_*”/“full_field” export types use the ",
          "deterministic image bounds (no edge-detection noise possible); ",
          "a real circular-disc detection (other export types) is the ",
          "case to scrutinize."),
        tags$table(class = "hover-tbl",
          tags$tr(tags$th("Run"), tags$th("Center (px)"), tags$th("Radius (px)"),
                  tags$th("Scale (µm/px)"), tags$th("Export Type"),
                  tags$th("Method"), tags$th("Detected")),
          rows
        )
      )
    )
  })

  # Fixed palette keyed by "Run N" — each run number always gets the same
  # colour regardless of which subset is currently visible.
  # Run 1 is golden yellow, not blue: blue reads poorly on the dark LDIR
  # scan image (light-blue particle blobs blend into the dot).
  .REPRO_RUN_PALETTE <- c(
    "Run 1" = "#e8d44d",   # golden yellow
    "Run 2" = "#ff7f0e",   # orange
    "Run 3" = "#2ca02c",   # green
    "Run 4" = "#9467bd",   # purple
    "Run 5" = "#8c564b"    # brown
  )

  # Run colour palette shared by the counts bar chart and the main scatter.
  .repro_run_palette <- function(run_levels) {
    .REPRO_RUN_PALETTE[paste("Run", run_levels)]
  }

  # Total particles detected per run — the raw, unfiltered per-run counts
  # (every particle from every run always ends up as its own consensus
  # particle if it matches nothing, so this equals the true ingest count;
  # see repro_link_runs() in R/reproducibility.R). Deliberately independent
  # of the material filter / non-reproducible-only toggle, mirroring how
  # the "Physical particles" summary panel is always unfiltered.
  output$repro_counts_barplot <- renderPlot({
    d <- repro_data()
    if (is.null(d)) return(NULL)
    pts <- d$points
    if (is.null(pts) || nrow(pts) == 0 || !"run" %in% names(pts)) return(NULL)

    run_levels <- sort(unique(pts$run))
    run_pal    <- .repro_run_palette(run_levels)
    counts     <- as.integer(table(factor(pts$run, levels = run_levels)))
    bar_df <- data.frame(
      run_lbl = factor(paste("Run", run_levels), levels = paste("Run", run_levels)),
      count   = counts, stringsAsFactors = FALSE
    )
    ggplot2::ggplot(bar_df, ggplot2::aes(x = run_lbl, y = count, fill = run_lbl)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.3, size = 5.2) +
      ggplot2::scale_fill_manual(values = run_pal, guide = "none") +
      ggplot2::scale_y_continuous(limits = c(0, max(bar_df$count) * 1.25)) +
      ggplot2::labs(x = NULL, y = "Particle Count", title = "Particles detected per run") +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                           margin = ggplot2::margin(b = 14)),
        plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10),
        axis.text.x = ggplot2::element_text(size = 14),
        panel.grid.major.x = ggplot2::element_blank()
      )
  })

  # Polymer composition per run: one donut per run, faceted so it scales to
  # however many runs the reproducibility tool was given (not fixed at 3).
  # Same family palette (.pie_palette) and Synthetic+Semi-synthetic filter as
  # the Summary tab's per-instrument pies, for a consistent "polymer
  # composition" definition app-wide. Unfiltered, like the counts bar chart.
  output$repro_composition_pies <- renderPlot({
    d <- repro_data()
    if (is.null(d)) return(NULL)
    pts <- d$points
    if (is.null(pts) || nrow(pts) == 0 || !"material" %in% names(pts)) return(NULL)

    fam <- classify_family_vec(pts$material)
    cat <- classify_category_vec(fam)
    keep <- cat %in% c("Synthetic", "Semi-synthetic")
    if (!any(keep)) return(NULL)

    run_levels <- sort(unique(pts$run))
    df <- data.frame(run = pts$run[keep], material = fam[keep], stringsAsFactors = FALSE)
    df$run_lbl <- factor(paste("Run", df$run), levels = paste("Run", run_levels))

    agg <- as.data.frame(table(run_lbl = df$run_lbl, material = df$material))
    names(agg) <- c("run_lbl", "material", "count")
    agg <- agg[agg$count > 0, , drop = FALSE]
    if (nrow(agg) == 0) return(NULL)
    totals <- stats::aggregate(count ~ run_lbl, agg, sum)
    names(totals)[2] <- "run_total"
    agg <- merge(agg, totals, by = "run_lbl")
    agg$pct   <- agg$count / agg$run_total * 100
    agg$label <- ifelse(agg$pct >= 5, paste0(round(agg$pct), "%"), "")

    strip_n <- setNames(paste0(totals$run_lbl, "\n(n = ", totals$run_total, ")"),
                        totals$run_lbl)
    agg$strip_lbl <- factor(strip_n[as.character(agg$run_lbl)],
                            levels = strip_n[paste("Run", run_levels)])

    fam_colors <- .pie_palette[sort(unique(as.character(agg$material)))]
    fam_colors[is.na(fam_colors)] <- "#cccccc"
    names(fam_colors) <- sort(unique(as.character(agg$material)))

    ggplot2::ggplot(agg, ggplot2::aes(x = "", y = count, fill = material)) +
      ggplot2::geom_col(width = 1, colour = "white", linewidth = 0.3,
                        position = ggplot2::position_fill()) +
      ggplot2::geom_text(ggplot2::aes(label = label),
                        position = ggplot2::position_fill(vjust = 0.5),
                        size = 3.4, colour = "white", fontface = "bold",
                        show.legend = FALSE) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::facet_wrap(~ strip_lbl, nrow = 1) +
      ggplot2::scale_fill_manual(values = fam_colors, name = "Material") +
      ggplot2::labs(title = "Polymer composition per run") +
      ggplot2::theme_void(base_size = 13) +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            margin = ggplot2::margin(b = 10)),
        strip.text  = ggplot2::element_text(face = "bold", size = 12),
        legend.position = "right"
      )
  })

  output$repro_summary_ui <- renderUI({
    d <- repro_data()
    if (is.null(d) || is.null(d$summary) || nrow(d$summary) == 0)
      return(tags$p(class = "text-muted", "Load a reproducibility output to see the summary."))
    s <- d$summary[1, ]
    pct <- function(v) if (is.null(v) || is.na(v)) "\u2014" else paste0(round(100 * v, 1), "%")
    val <- function(nm) if (nm %in% names(s)) s[[nm]] else NA
    tags$table(class = "hover-tbl",
      tags$tr(tags$td("Instrument"), tags$td(tags$b(as.character(val("instrument"))))),
      tags$tr(tags$td("Runs"), tags$td(val("n_runs"))),
      tags$tr(tags$td("Physical particles"), tags$td(val("n_consensus"))),
      tags$tr(tags$td("Detected in all"),
              tags$td(paste0(val("detected_in_all"), " (", pct(val("detected_in_all_frac")), ")"))),
      tags$tr(tags$td("Count CV"), tags$td(pct(val("count_cv")))),
      tags$tr(tags$td("Material concordance"), tags$td(pct(val("material_concordance")))),
      tags$tr(tags$td("Accuracy vs reference"), tags$td(pct(val("accuracy_vs_reference")))),
      tags$tr(tags$td("Median Feret CV"), tags$td(pct(val("median_feret_cv")))),
      tags$tr(tags$td("Median jitter"),
              tags$td(if (is.na(val("median_pos_jitter_um"))) "\u2014"
                      else paste0(round(val("median_pos_jitter_um"), 1), " \u00b5m")))
    )
  })
}

# ============================================================================
# Run
# ============================================================================
shinyApp(ui = ui, server = server)
