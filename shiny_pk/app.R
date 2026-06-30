library(shiny)
library(bslib)
library(shinyjs)
library(plotly)
library(rhandsontable)
library(DT)
library(readxl)
library(writexl)
library(dplyr)
library(PKNCA)
library(gridExtra)
library(grid)
library(ggplot2)

source("R/pk_functions.R")
source("R/mod_data.R")
source("R/mod_nca.R")
source("R/mod_modeling.R")
source("R/mod_export.R")

# ── Theme ─────────────────────────────────────────────────────────────────────

pk_theme <- bs_theme(
  bootswatch = "flatly",
  base_font  = font_google("Inter")
)

# Disabled nav-link style (class pk-disabled, toggled via shinyjs)
css_nav_disabled <- "
a.nav-link.pk-disabled {
  opacity        : 0.40;
  pointer-events : none;
  cursor         : not-allowed;
}
"

# Tab values that depend on validated data
CONTROLLED_TABS <- c("tab_nca", "tab_modeling", "tab_export")

# JS snippets to add / remove pk-disabled class on a tab by data-value
.js_disable <- function(tab)
  sprintf("$('a.nav-link[data-value=\"%s\"]').addClass('pk-disabled');", tab)

.js_enable <- function(tab)
  sprintf("$('a.nav-link[data-value=\"%s\"]').removeClass('pk-disabled');", tab)

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_navbar(
  title = tagList(
    tags$span(class = "fw-bold",          "PK"),
    tags$span(class = "fw-normal text-muted ms-1", "Shiny")
  ),
  theme = pk_theme,
  id    = "main_nav",

  # Dependencies + global CSS injected once in the page header
  header = tagList(
    useShinyjs(),
    tags$style(css_nav_disabled)
  ),

  # ── Tab 1 : Données ─────────────────────────────────────────────────────────
  nav_panel(
    title = tagList(tags$i(class = "bi bi-table me-1"), "Données"),
    value = "tab_data",
    mod_data_ui("data")
  ),

  # ── Tab 2 : NCA ──────────────────────────────────────────────────────────────
  nav_panel(
    title = tagList(tags$i(class = "bi bi-bar-chart-line me-1"), "NCA"),
    value = "tab_nca",
    mod_nca_ui("nca")
  ),

  # ── Tab 3 : Modélisation ─────────────────────────────────────────────────────
  nav_panel(
    title = tagList(tags$i(class = "bi bi-graph-up me-1"), "Modélisation"),
    value = "tab_modeling",
    mod_modeling_ui("modeling")
  ),

  # ── Tab 4 : Export ───────────────────────────────────────────────────────────
  nav_panel(
    title = tagList(tags$i(class = "bi bi-download me-1"), "Export"),
    value = "tab_export",
    mod_export_ui("export")
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Module wiring
  data_out  <- mod_data_server("data")
  nca_out   <- mod_nca_server("nca",       data_out)
  model_out <- mod_modeling_server("modeling", data_out)
  mod_export_server("export", nca_out, model_out, data_out)

  # ── Tab gating : disable on startup, enable after data validation ────────────

  # Apply initial disabled state as soon as the DOM is ready
  session$onFlushed(function() {
    shinyjs::runjs(paste(.js_disable(CONTROLLED_TABS), collapse = "\n"))
  }, once = TRUE)

  # Re-evaluate whenever data_out changes (validation button click)
  observe({
    validated <- tryCatch(
      !is.null(data_out()),
      error   = function(e) FALSE,
      warning = function(w) FALSE
    )

    if (validated) {
      shinyjs::runjs(paste(.js_enable(CONTROLLED_TABS), collapse = "\n"))
    } else {
      shinyjs::runjs(paste(.js_disable(CONTROLLED_TABS), collapse = "\n"))
    }
  })
}

shinyApp(ui, server)
