library(shiny)
library(bslib)
library(plotly)
library(DT)
library(dplyr)

# ── UI ────────────────────────────────────────────────────────────────────────

mod_modeling_ui <- function(id) {
  ns <- NS(id)

  layout_columns(
    col_widths = c(4, 8),

    # ── Left panel : controls ─────────────────────────────────────────────────
    card(
      card_header("Paramètres du modèle"),

      tags$b("Valeurs initiales"),

      numericInput(ns("init_CL"), "CL  (L/h)",  value = 0.003,  min = 1e-9, step = 0.001),
      numericInput(ns("init_V1"), "V1  (L)",     value = 0.05,  min = 1e-9, step = 0.01),
      numericInput(ns("init_V2"), "V2  (L)",     value = 0.06,  min = 1e-9, step = 0.01),
      numericInput(ns("init_Q"),  "Q   (L/h)",   value = 0.004, min = 1e-9, step = 0.001),

      tags$small(class = "text-muted",
        "Les deux modèles (1-cmt et 2-cmt) sont ajustés automatiquement."
      ),

      hr(),
      actionButton(
        ns("run_fit"),
        tagList(tags$i(class = "bi bi-play-fill me-1"), "Lancer l'optimisation"),
        class = "btn-success w-100"
      ),

      hr(),
      tags$b("Modèle affiché dans le graphique"),
      radioButtons(
        ns("plot_model"), NULL,
        choices  = c("Meilleur (AIC)" = "best",
                     "1 compartiment" = "1",
                     "2 compartiments" = "2"),
        selected = "best"
      )
    ),

    # ── Right panel : results ─────────────────────────────────────────────────
    tagList(

      # Comparison card
      uiOutput(ns("comparison_card")),

      # Detailed parameters of selected model
      card(
        card_header(uiOutput(ns("param_table_header"))),
        DT::dataTableOutput(ns("param_table"))
      ),

      card(
        card_header("Ajustement : observé vs prédit"),
        plotly::plotlyOutput(ns("fit_plot"), height = "420px")
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

mod_modeling_server <- function(id, pk_data) {
  moduleServer(id, function(input, output, session) {

    fit_result <- reactiveVal(NULL)

    # ── Fitting — triggered by button, fits BOTH models ───────────────────────

    observeEvent(input$run_fit, {
      d <- isolate(req(pk_data()))

      inits <- list(
        CL = max(isolate(input$init_CL), 1e-6),
        V1 = max(isolate(input$init_V1), 1e-6),
        V2 = max(isolate(input$init_V2), 1e-6),
        Q  = max(isolate(input$init_Q),  1e-6)
      )

      withProgress(message = "Optimisation PK en cours…", value = 0, {

        setProgress(0.1, detail = "Ajustement 1 compartiment…")

        fit1 <- tryCatch(
          fit_pk_model(
            df          = d$data,
            dose_col    = d$dose_col,
            time_col    = d$time_col,
            conc_col    = d$conc_col,
            animal_col  = d$animal_col,
            n_comp      = 1L,
            init_params = inits[c("CL", "V1")]
          ),
          error = function(e) list(ok = FALSE, message = conditionMessage(e))
        )

        setProgress(0.4, detail = "Ajustement 2 compartiments…")

        fit2 <- tryCatch(
          fit_pk_model(
            df          = d$data,
            dose_col    = d$dose_col,
            time_col    = d$time_col,
            conc_col    = d$conc_col,
            animal_col  = d$animal_col,
            n_comp      = 2L,
            init_params = inits
          ),
          error = function(e) list(ok = FALSE, message = conditionMessage(e))
        )

        setProgress(0.8, detail = "Simulation des courbes…")

        t_max        <- max(d$data[[d$time_col]], na.rm = TRUE)
        t_sim        <- seq(0, t_max, length.out = 300)
        dose_by_grp  <- tapply(d$data[[d$dose_col]], d$data[[d$animal_col]], mean, na.rm = TRUE)
        unique_doses <- sort(unique(round(unname(dose_by_grp), 8)))

        sim_for <- function(fit, nc) {
          if (is.null(fit$params)) return(NULL)
          do.call(rbind, lapply(unique_doses, function(dv) {
            s <- simulate_pk(fit$params, dv, t_sim, nc)
            data.frame(time = s$time, conc = s$conc, dose = dv)
          }))
        }

        sim1 <- if (!is.null(fit1$params)) sim_for(fit1, 1L) else NULL
        sim2 <- if (!is.null(fit2$params)) sim_for(fit2, 2L) else NULL

        setProgress(1)
      })

      # Notifications
      for (info in list(list(fit1, "1-cmt"), list(fit2, "2-cmt"))) {
        ft <- info[[1]]; label <- info[[2]]
        if (!is.null(ft$convergence) && ft$convergence != 0)
          showNotification(
            paste0(label, " — convergence douteuse (code ", ft$convergence, ")"),
            type = "warning", duration = 10
          )
      }

      showNotification("Optimisation terminée.", type = "message", duration = 4)

      fit_result(list(
        fit1  = fit1,  sim1 = sim1,
        fit2  = fit2,  sim2 = sim2
      ))
    })

    # ── Comparison card ───────────────────────────────────────────────────────

    output$comparison_card <- renderUI({
      req(fit_result())
      r <- fit_result()

      aic1 <- r$fit1$AIC;  bic1 <- r$fit1$BIC
      aic2 <- r$fit2$AIC;  bic2 <- r$fit2$BIC

      fmt <- function(x) if (is.null(x) || !is.finite(x)) "—" else sprintf("%.1f", x)

      winner_aic <- if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2)) {
        if (aic1 < aic2) "1-cmt" else "2-cmt"
      } else "—"

      delta_aic <- if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2))
        sprintf("ΔAIC = %.1f", abs(aic2 - aic1)) else ""

      winner_color <- function(model) {
        if (winner_aic == model) "text-success fw-bold" else "text-muted"
      }

      card(
        card_header(tagList(
          tags$i(class = "bi bi-bar-chart-fill me-1"),
          "Comparaison des modèles",
          tags$span(class = "ms-2 badge bg-info text-dark", delta_aic)
        )),
        layout_columns(
          col_widths = c(6, 6),

          # 1-compartiment
          card(
            class = if (winner_aic == "1-cmt") "border-success" else "",
            card_header(tagList(
              tags$span(class = winner_color("1-cmt"), "1 compartiment"),
              if (winner_aic == "1-cmt")
                tags$span(class = "ms-2 badge bg-success", "Meilleur AIC")
            )),
            tags$table(
              class = "table table-sm mb-0",
              tags$tbody(
                tags$tr(tags$td("AIC"), tags$td(class = winner_color("1-cmt"), fmt(aic1))),
                tags$tr(tags$td("BIC"), tags$td(fmt(bic1))),
                tags$tr(tags$td("n param."), tags$td("2  (CL, V1)"))
              )
            )
          ),

          # 2-compartiments
          card(
            class = if (winner_aic == "2-cmt") "border-success" else "",
            card_header(tagList(
              tags$span(class = winner_color("2-cmt"), "2 compartiments"),
              if (winner_aic == "2-cmt")
                tags$span(class = "ms-2 badge bg-success", "Meilleur AIC")
            )),
            tags$table(
              class = "table table-sm mb-0",
              tags$tbody(
                tags$tr(tags$td("AIC"), tags$td(class = winner_color("2-cmt"), fmt(aic2))),
                tags$tr(tags$td("BIC"), tags$td(fmt(bic2))),
                tags$tr(tags$td("n param."), tags$td("4  (CL, V1, V2, Q)"))
              )
            )
          )
        ),
        if (nchar(delta_aic) > 0) tags$p(
          class = "text-muted small mt-2 mb-0 px-1",
          if (abs(r$fit2$AIC - r$fit1$AIC) < 2)
            "ΔAIC < 2 : les deux modèles sont équivalents — préférer le plus parcimonieux (1-cmt)."
          else if (abs(r$fit2$AIC - r$fit1$AIC) < 10)
            "ΔAIC 2–10 : légère préférence pour le meilleur AIC."
          else
            "ΔAIC > 10 : forte préférence pour le meilleur AIC."
        )
      )
    })

    # ── Active model (for table + plot) ──────────────────────────────────────

    active_fit <- reactive({
      req(fit_result())
      r <- fit_result()

      choice <- input$plot_model
      if (choice == "best") {
        aic1 <- r$fit1$AIC; aic2 <- r$fit2$AIC
        if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2))
          choice <- if (aic1 <= aic2) "1" else "2"
        else
          choice <- "2"
      }

      if (choice == "1")
        list(fit = r$fit1, sim = r$sim1, n_comp = 1L, label = "1 compartiment")
      else
        list(fit = r$fit2, sim = r$sim2, n_comp = 2L, label = "2 compartiments")
    })

    # ── Parameter table header ────────────────────────────────────────────────

    output$param_table_header <- renderUI({
      req(active_fit())
      tagList(
        "Paramètres estimés — ",
        tags$span(class = "text-primary", active_fit()$label)
      )
    })

    # ── Parameter table ───────────────────────────────────────────────────────

    output$param_table <- DT::renderDataTable({
      req(active_fit())
      ft <- active_fit()$fit

      rse_vec <- setNames(
        unname(ft$rse),
        gsub("^%RSE_", "", names(ft$rse))
      )

      primary_df <- data.frame(
        Paramètre = names(ft$params),
        Valeur    = unname(ft$params),
        `%RSE`    = unname(rse_vec[names(ft$params)]),
        check.names = FALSE
      )

      derived_vals <- unlist(ft$derived)
      derived_df <- data.frame(
        Paramètre = names(derived_vals),
        Valeur    = unname(derived_vals),
        `%RSE`    = NA_real_,
        check.names = FALSE
      )

      crit_df <- data.frame(
        Paramètre = c("AIC", "BIC", "RSS", "n_obs"),
        Valeur    = c(ft$AIC, ft$BIC, ft$RSS, as.numeric(ft$n_obs)),
        `%RSE`    = NA_real_,
        check.names = FALSE
      )

      tbl <- rbind(primary_df, derived_df, crit_df)

      DT::datatable(
        tbl,
        rownames = FALSE,
        options  = list(
          dom        = "t",
          pageLength = 25,
          scrollX    = FALSE,
          columnDefs = list(
            list(
              targets = 1,
              render  = DT::JS(
                "function(data, type, row) {",
                "  if (type !== 'display' || data === null || data === undefined) return data;",
                "  var v = parseFloat(data);",
                "  if (isNaN(v)) return data;",
                "  if (Math.abs(v) > 0 && Math.abs(v) < 0.001) return v.toExponential(4);",
                "  return v.toPrecision(5);",
                "}"
              )
            )
          )
        ),
        class = "compact stripe"
      ) %>%
        DT::formatRound("%RSE", digits = 2)
    })

    # ── Fit plot — both curves overlaid ──────────────────────────────────────

    output$fit_plot <- plotly::renderPlotly({
      req(fit_result(), pk_data())
      r  <- fit_result()
      d  <- pk_data()
      df <- d$data

      animals <- unique(df[[d$animal_col]])
      p <- plotly::plot_ly()

      # Observed points
      for (anim in animals) {
        sub <- df[df[[d$animal_col]] == anim, ]
        p <- plotly::add_trace(
          p,
          x      = sub[[d$time_col]],
          y      = sub[[d$conc_col]],
          type   = "scatter",
          mode   = "markers",
          name   = as.character(anim),
          marker = list(size = 9, symbol = "circle-open", line = list(width = 2))
        )
      }

      # 1-cmt curves (dotted)
      if (!is.null(r$sim1)) {
        for (dv in sort(unique(r$sim1$dose))) {
          sub_sim <- r$sim1[r$sim1$dose == dv, ]
          p <- plotly::add_trace(
            p,
            x          = sub_sim$time,
            y          = sub_sim$conc,
            type       = "scatter",
            mode       = "lines",
            name       = paste0("1-cmt ", dv, " mg/kg"),
            line       = list(width = 1.5, dash = "dot", color = "grey"),
            showlegend = TRUE,
            legendgroup = "1cmt"
          )
        }
      }

      # 2-cmt curves (dashed)
      if (!is.null(r$sim2)) {
        for (dv in sort(unique(r$sim2$dose))) {
          sub_sim <- r$sim2[r$sim2$dose == dv, ]
          p <- plotly::add_trace(
            p,
            x          = sub_sim$time,
            y          = sub_sim$conc,
            type       = "scatter",
            mode       = "lines",
            name       = paste0("2-cmt ", dv, " mg/kg"),
            line       = list(width = 2, dash = "dash"),
            showlegend = TRUE,
            legendgroup = "2cmt"
          )
        }
      }

      p %>% plotly::layout(
        xaxis     = list(title = paste0("Temps  (", d$time_col, ")")),
        yaxis     = list(
          title = paste0("Concentration  (", d$conc_col, ")"),
          type  = "log"
        ),
        legend    = list(title = list(text = "Série")),
        hovermode = "x unified"
      )
    })

    return(fit_result)
  })
}
