library(shiny)
library(bslib)
library(PKNCA)
library(plotly)
library(DT)
library(dplyr)

# ── UI ────────────────────────────────────────────────────────────────────────

mod_nca_ui <- function(id) {
  ns <- NS(id)
  tagList(

    # Mini récap : n animaux × n points
    uiOutput(ns("data_summary_bar")),

    card(
      card_header("Profils concentration-temps individuels"),
      plotly::plotlyOutput(ns("pk_plot"), height = "440px")
    ),

    card(
      card_header("Résultats NCA par animal (1 ligne = résumé de N mesures)"),
      uiOutput(ns("nca_status")),
      DT::dataTableOutput(ns("nca_table"))
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

mod_nca_server <- function(id, pk_data) {
  moduleServer(id, function(input, output, session) {

    # ── Nombre de points par animal (réutilisé dans plot + tableau) ───────────

    obs_per_animal <- reactive({
      d  <- req(pk_data())
      df <- d$data
      tapply(df[[d$time_col]], df[[d$animal_col]], length)
    })

    # ── Bandeau récapitulatif ─────────────────────────────────────────────────

    output$data_summary_bar <- renderUI({
      d   <- req(pk_data())
      ops <- obs_per_animal()

      n_anim  <- length(ops)
      n_total <- sum(ops)
      n_min   <- min(ops)
      n_max   <- max(ops)
      pts_lbl <- if (n_min == n_max)
        paste0(n_min, " mesures / animal")
      else
        paste0(n_min, "–", n_max, " mesures / animal")

      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(title = "Animaux",          value = n_anim,  theme = "primary"),
        value_box(title = "Observations totales", value = n_total, theme = "primary"),
        value_box(title = "Points temporels", value = pts_lbl, theme = "primary")
      )
    })

    # ── NCA computation ───────────────────────────────────────────────────────

    nca_results <- reactive({
      d <- req(pk_data())

      withProgress(message = "Calcul NCA…", value = 0.5, {
        tryCatch(
          run_nca(
            df         = d$data,
            dose_col   = d$dose_col,
            time_col   = d$time_col,
            conc_col   = d$conc_col,
            animal_col = d$animal_col,
            route      = "iv"
          ),
          error = function(e) {
            showNotification(
              paste0("Erreur NCA : ", conditionMessage(e)),
              type = "error", duration = 15
            )
            NULL
          }
        )
      })
    })

    # ── Concentration-time plot ───────────────────────────────────────────────

    output$pk_plot <- plotly::renderPlotly({
      d       <- req(pk_data())
      df      <- d$data
      animals <- unique(df[[d$animal_col]])
      ops     <- obs_per_animal()

      p <- plotly::plot_ly()

      for (anim in animals) {
        sub   <- df[df[[d$animal_col]] == anim, ]
        sub   <- sub[order(sub[[d$time_col]]), ]
        n_pts <- unname(ops[as.character(anim)])

        # Légende : "Animal_A  (n = 8 pts)"
        trace_name <- paste0(
          "<b>", as.character(anim), "</b>",
          "  <span style='color:#888;font-size:11px'>(n = ", n_pts, " pts)</span>"
        )

        p <- plotly::add_trace(
          p,
          x             = sub[[d$time_col]],
          y             = sub[[d$conc_col]],
          type          = "scatter",
          mode          = "lines+markers",
          name          = trace_name,
          marker        = list(size = 7, symbol = "circle"),
          line          = list(width = 1.8),
          hovertemplate = paste0(
            "<b>", d$animal_col, " : ", as.character(anim), "</b><br>",
            d$time_col, " : %{x}<br>",
            d$conc_col, " : %{y:.4g}<extra></extra>"
          )
        )
      }

      p %>% plotly::layout(
        xaxis = list(
          title = paste0("Temps  (", d$time_col, ")"),
          zeroline = FALSE
        ),
        yaxis = list(
          title = paste0("Concentration  (", d$conc_col, ")"),
          type  = "log"
        ),
        legend = list(
          title = list(
            text = paste0(
              d$animal_col,
              "<br><sup style='color:#888'>chaque trait = 1 profil individuel</sup>"
            )
          )
        ),
        hovermode = "closest"
      )
    })

    # ── Status + NCA table ────────────────────────────────────────────────────

    output$nca_status <- renderUI({
      req(pk_data())
      res <- nca_results()
      if (is.null(res)) {
        tags$p(class = "text-danger mb-1",
               tags$i(class = "bi bi-exclamation-triangle-fill me-1"),
               "Impossible de calculer les paramètres NCA.")
      } else {
        tags$p(class = "text-muted small mb-1",
               tags$i(class = "bi bi-info-circle me-1"),
               "Chaque ligne résume le profil concentration-temps complet d'un animal.")
      }
    })

    output$nca_table <- DT::renderDataTable({
      d   <- req(pk_data())
      res <- req(nca_results())
      ops <- obs_per_animal()

      # Ajoute colonne n_obs pour lever toute ambiguïté
      n_obs_vec <- unname(ops[as.character(res[[d$animal_col]])])
      res_aug   <- cbind(res, n_obs = n_obs_vec)

      # Déplace n_obs en 2e colonne (juste après l'ID animal)
      cols_order <- c(names(res_aug)[1], "n_obs",
                      setdiff(names(res_aug), c(names(res_aug)[1], "n_obs")))
      res_aug <- res_aug[, cols_order, drop = FALSE]

      # Rename n_obs for display
      names(res_aug)[names(res_aug) == "n_obs"] <- "n mesures"
      num_cols <- names(res_aug)[sapply(res_aug, is.numeric) &
                                  names(res_aug) != "n mesures"]

      DT::datatable(
        res_aug,
        rownames = FALSE,
        options  = list(scrollX = TRUE, pageLength = 15, dom = "tip"),
        class    = "compact stripe"
      ) %>%
        DT::formatRound(columns = num_cols, digits = 4) %>%
        DT::formatStyle(
          "n mesures",
          fontWeight = "bold",
          color      = "#0d6efd"
        )
    })

    return(nca_results)
  })
}
