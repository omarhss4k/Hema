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
      tags$b("Modèle affiché"),
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

      uiOutput(ns("comparison_card")),

      card(
        card_header(uiOutput(ns("param_table_header"))),
        DT::dataTableOutput(ns("param_table"))
      ),

      navset_card_tab(
        nav_panel(
          title = tagList(tags$i(class = "bi bi-graph-up me-1"), "Ajustement"),
          plotly::plotlyOutput(ns("fit_plot"), height = "420px")
        ),
        nav_panel(
          title = tagList(tags$i(class = "bi bi-scatter-chart me-1"), "Résidus"),
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header("Résidus vs Temps"),
              plotly::plotlyOutput(ns("resid_vs_time"), height = "340px")
            ),
            card(
              card_header("Observé vs Prédit"),
              plotly::plotlyOutput(ns("obs_vs_pred"), height = "340px")
            )
          )
        ),
        nav_panel(
          title = tagList(tags$i(class = "bi bi-people me-1"), "Individuel"),
          layout_columns(
            col_widths = c(12),
            card(
              card_header(uiOutput(ns("ind_table_header"))),
              DT::dataTableOutput(ns("ind_table"))
            )
          ),
          card(
            card_header("Courbes individuelles"),
            plotly::plotlyOutput(ns("ind_plot"), height = "420px")
          )
        )
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

mod_modeling_server <- function(id, pk_data) {
  moduleServer(id, function(input, output, session) {

    fit_result  <- reactiveVal(NULL)
    ind_result  <- reactiveVal(NULL)

    # ── Fitting ───────────────────────────────────────────────────────────────

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
          fit_pk_model(d$data, d$dose_col, d$time_col, d$conc_col, d$animal_col,
                       n_comp = 1L, init_params = inits[c("CL","V1")]),
          error = function(e) list(ok = FALSE, message = conditionMessage(e))
        )

        setProgress(0.4, detail = "Ajustement 2 compartiments…")
        fit2 <- tryCatch(
          fit_pk_model(d$data, d$dose_col, d$time_col, d$conc_col, d$animal_col,
                       n_comp = 2L, init_params = inits),
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

        setProgress(0.9, detail = "Estimation individuelle…")

        ind1 <- tryCatch(
          fit_pk_individual(d$data, d$dose_col, d$time_col, d$conc_col, d$animal_col,
                            n_comp = 1L, init_params = inits[c("CL","V1")]),
          error = function(e) NULL
        )
        ind2 <- tryCatch(
          fit_pk_individual(d$data, d$dose_col, d$time_col, d$conc_col, d$animal_col,
                            n_comp = 2L, init_params = inits),
          error = function(e) NULL
        )

        setProgress(1)
      })

      ind_result(list(ind1 = ind1, ind2 = ind2))

      for (info in list(list(fit1,"1-cmt"), list(fit2,"2-cmt"))) {
        ft <- info[[1]]; label <- info[[2]]
        if (!is.null(ft$convergence) && ft$convergence != 0)
          showNotification(paste0(label, " — convergence douteuse (code ", ft$convergence, ")"),
                           type = "warning", duration = 10)
      }
      showNotification("Optimisation terminée.", type = "message", duration = 4)

      fit_result(list(fit1=fit1, sim1=sim1, fit2=fit2, sim2=sim2))
    })

    # ── Comparison card ───────────────────────────────────────────────────────

    output$comparison_card <- renderUI({
      req(fit_result())
      r <- fit_result()

      aic1 <- r$fit1$AIC; bic1 <- r$fit1$BIC
      aic2 <- r$fit2$AIC; bic2 <- r$fit2$BIC
      fmt  <- function(x) if (is.null(x) || !is.finite(x)) "—" else sprintf("%.1f", x)

      winner_aic <- if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2))
        if (aic1 < aic2) "1-cmt" else "2-cmt" else "—"

      delta_aic <- if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2))
        sprintf("ΔAIC = %.1f", abs(aic2 - aic1)) else ""

      winner_color <- function(m) if (winner_aic == m) "text-success fw-bold" else "text-muted"

      card(
        card_header(tagList(
          tags$i(class = "bi bi-bar-chart-fill me-1"), "Comparaison des modèles",
          tags$span(class = "ms-2 badge bg-info text-dark", delta_aic)
        )),
        layout_columns(
          col_widths = c(6, 6),
          card(
            class = if (winner_aic == "1-cmt") "border-success" else "",
            card_header(tagList(
              tags$span(class = winner_color("1-cmt"), "1 compartiment"),
              if (winner_aic == "1-cmt") tags$span(class = "ms-2 badge bg-success", "Meilleur AIC")
            )),
            tags$table(class = "table table-sm mb-0", tags$tbody(
              tags$tr(tags$td("AIC"), tags$td(class = winner_color("1-cmt"), fmt(aic1))),
              tags$tr(tags$td("BIC"), tags$td(fmt(bic1))),
              tags$tr(tags$td("n param."), tags$td("2  (CL, V1)"))
            ))
          ),
          card(
            class = if (winner_aic == "2-cmt") "border-success" else "",
            card_header(tagList(
              tags$span(class = winner_color("2-cmt"), "2 compartiments"),
              if (winner_aic == "2-cmt") tags$span(class = "ms-2 badge bg-success", "Meilleur AIC")
            )),
            tags$table(class = "table table-sm mb-0", tags$tbody(
              tags$tr(tags$td("AIC"), tags$td(class = winner_color("2-cmt"), fmt(aic2))),
              tags$tr(tags$td("BIC"), tags$td(fmt(bic2))),
              tags$tr(tags$td("n param."), tags$td("4  (CL, V1, V2, Q)"))
            ))
          )
        ),
        if (nchar(delta_aic) > 0) tags$p(
          class = "text-muted small mt-2 mb-0 px-1",
          if (abs(aic2 - aic1) < 2)
            "ΔAIC < 2 : modèles équivalents — préférer le plus parcimonieux (1-cmt)."
          else if (abs(aic2 - aic1) < 10)
            "ΔAIC 2–10 : légère préférence pour le meilleur AIC."
          else
            "ΔAIC > 10 : forte préférence pour le meilleur AIC."
        )
      )
    })

    # ── Active model ──────────────────────────────────────────────────────────

    active_fit <- reactive({
      req(fit_result())
      r      <- fit_result()
      choice <- input$plot_model
      if (choice == "best") {
        aic1 <- r$fit1$AIC; aic2 <- r$fit2$AIC
        choice <- if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2))
          if (aic1 <= aic2) "1" else "2" else "2"
      }
      if (choice == "1") list(fit=r$fit1, sim=r$sim1, n_comp=1L, label="1 compartiment")
      else               list(fit=r$fit2, sim=r$sim2, n_comp=2L, label="2 compartiments")
    })

    # ── Residuals helper ──────────────────────────────────────────────────────

    resid_df <- reactive({
      req(active_fit(), pk_data())
      af <- active_fit()
      d  <- pk_data()
      df <- d$data

      dose_by_grp <- tapply(df[[d$dose_col]], df[[d$animal_col]], mean, na.rm = TRUE)

      rows <- lapply(unique(df[[d$animal_col]]), function(id) {
        sub  <- df[df[[d$animal_col]] == id, ]
        dv   <- unname(dose_by_grp[as.character(id)])
        t_obs <- sub[[d$time_col]]
        C_obs <- sub[[d$conc_col]]
        C_pred <- if (af$n_comp == 1L)
          .pred_1comp(dv, af$fit$params["CL"], af$fit$params["V1"], t_obs)
        else
          .pred_2comp(dv, af$fit$params["CL"], af$fit$params["V1"],
                      af$fit$params["V2"], af$fit$params["Q"], t_obs)
        ok <- is.finite(C_pred) & C_pred > 0 & C_obs > 0
        data.frame(
          animal   = as.character(id),
          time     = t_obs[ok],
          obs      = C_obs[ok],
          pred     = C_pred[ok],
          resid    = log(C_obs[ok]) - log(C_pred[ok]),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    })

    # ── Parameter table header ────────────────────────────────────────────────

    output$param_table_header <- renderUI({
      req(active_fit())
      tagList("Paramètres estimés — ",
              tags$span(class = "text-primary", active_fit()$label))
    })

    # ── Parameter table ───────────────────────────────────────────────────────

    output$param_table <- DT::renderDataTable({
      req(active_fit())
      ft <- active_fit()$fit

      rse_vec <- setNames(unname(ft$rse), gsub("^%RSE_","", names(ft$rse)))

      primary_df <- data.frame(
        Parametre = names(ft$params),
        Valeur    = unname(ft$params),
        `%RSE`    = unname(rse_vec[names(ft$params)]),
        check.names = FALSE
      )
      derived_df <- data.frame(
        Parametre = names(unlist(ft$derived)),
        Valeur    = unname(unlist(ft$derived)),
        `%RSE`    = NA_real_,
        check.names = FALSE
      )
      crit_df <- data.frame(
        Parametre = c("AIC","BIC","RSS","n_obs"),
        Valeur    = c(ft$AIC, ft$BIC, ft$RSS, as.numeric(ft$n_obs)),
        `%RSE`    = NA_real_,
        check.names = FALSE
      )
      tbl <- rbind(primary_df, derived_df, crit_df)

      DT::datatable(tbl, rownames = FALSE,
        options = list(dom="t", pageLength=25, scrollX=FALSE,
          columnDefs = list(list(targets=1, render=DT::JS(
            "function(data,type,row){",
            "  if(type!=='display'||data===null||data===undefined) return data;",
            "  var v=parseFloat(data); if(isNaN(v)) return data;",
            "  if(Math.abs(v)>0&&Math.abs(v)<0.001) return v.toExponential(4);",
            "  return v.toPrecision(5);","}"
          )))
        ), class="compact stripe"
      ) %>% DT::formatRound("%RSE", digits=2)
    })

    # ── Fit plot ──────────────────────────────────────────────────────────────

    output$fit_plot <- plotly::renderPlotly({
      req(fit_result(), pk_data())
      r  <- fit_result()
      d  <- pk_data()
      df <- d$data
      p  <- plotly::plot_ly()

      for (anim in unique(df[[d$animal_col]])) {
        sub <- df[df[[d$animal_col]] == anim, ]
        p <- plotly::add_trace(p, x=sub[[d$time_col]], y=sub[[d$conc_col]],
          type="scatter", mode="markers", name=as.character(anim),
          marker=list(size=9, symbol="circle-open", line=list(width=2)))
      }

      if (!is.null(r$sim1)) {
        for (dv in sort(unique(r$sim1$dose))) {
          s <- r$sim1[r$sim1$dose==dv,]
          p <- plotly::add_trace(p, x=s$time, y=s$conc, type="scatter", mode="lines",
            name=paste0("1-cmt ",dv," mg/kg"),
            line=list(width=1.5, dash="dot", color="grey"),
            legendgroup="1cmt", showlegend=TRUE)
        }
      }
      if (!is.null(r$sim2)) {
        for (dv in sort(unique(r$sim2$dose))) {
          s <- r$sim2[r$sim2$dose==dv,]
          p <- plotly::add_trace(p, x=s$time, y=s$conc, type="scatter", mode="lines",
            name=paste0("2-cmt ",dv," mg/kg"),
            line=list(width=2, dash="dash"),
            legendgroup="2cmt", showlegend=TRUE)
        }
      }

      p %>% plotly::layout(
        xaxis=list(title=paste0("Temps (",d$time_col,")")),
        yaxis=list(title=paste0("Concentration (",d$conc_col,")"), type="log"),
        legend=list(title=list(text="Série")), hovermode="x unified"
      )
    })

    # ── Résidus vs Temps ──────────────────────────────────────────────────────

    output$resid_vs_time <- plotly::renderPlotly({
      req(resid_df())
      rd <- resid_df()
      p  <- plotly::plot_ly()

      for (anim in unique(rd$animal)) {
        sub <- rd[rd$animal == anim, ]
        p <- plotly::add_trace(p, x=sub$time, y=sub$resid,
          type="scatter", mode="markers",
          name=anim, marker=list(size=8))
      }

      p %>% plotly::layout(
        xaxis = list(title = "Temps"),
        yaxis = list(title = "Résidu log (obs - prédit)", zeroline = TRUE,
                     zerolinecolor = "#e74c3c", zerolinewidth = 2),
        shapes = list(list(type="line", x0=0, x1=1, xref="paper",
                           y0=0, y1=0, line=list(color="#e74c3c", width=2, dash="dash"))),
        hovermode = "closest",
        legend = list(title = list(text = "Animal"))
      )
    })

    # ── Observé vs Prédit ─────────────────────────────────────────────────────

    output$obs_vs_pred <- plotly::renderPlotly({
      req(resid_df())
      rd  <- resid_df()
      lim <- range(c(rd$obs, rd$pred), na.rm = TRUE)

      p <- plotly::plot_ly()

      for (anim in unique(rd$animal)) {
        sub <- rd[rd$animal == anim, ]
        p <- plotly::add_trace(p, x=sub$pred, y=sub$obs,
          type="scatter", mode="markers",
          name=anim, marker=list(size=8))
      }

      # Ligne d'identité y = x
      p <- plotly::add_trace(p,
        x = lim, y = lim,
        type = "scatter", mode = "lines",
        name = "y = x",
        line = list(color="#e74c3c", width=2, dash="dash"),
        showlegend = TRUE
      )

      p %>% plotly::layout(
        xaxis = list(title = "Concentration prédite", type = "log"),
        yaxis = list(title = "Concentration observée", type = "log"),
        hovermode = "closest",
        legend = list(title = list(text = "Animal"))
      )
    })

    # ── Individual results ────────────────────────────────────────────────────

    active_ind <- reactive({
      req(ind_result())
      r      <- ind_result()
      choice <- input$plot_model
      if (choice == "best") {
        aic1 <- fit_result()$fit1$AIC; aic2 <- fit_result()$fit2$AIC
        choice <- if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2))
          if (aic1 <= aic2) "1" else "2" else "2"
      }
      if (choice == "1") list(ind=r$ind1, n_comp=1L, label="1 compartiment")
      else               list(ind=r$ind2, n_comp=2L, label="2 compartiments")
    })

    output$ind_table_header <- renderUI({
      req(active_ind())
      tagList("Parametres individuels — ",
              tags$span(class="text-primary", active_ind()$label))
    })

    output$ind_table <- DT::renderDataTable({
      req(active_ind())
      df <- active_ind()$ind
      if (is.null(df)) return(NULL)

      num_cols <- names(df)[sapply(df, is.numeric)]
      DT::datatable(df, rownames=FALSE,
        options=list(dom="t", scrollX=TRUE, pageLength=20),
        class="compact stripe"
      ) %>% DT::formatRound(columns=num_cols, digits=4)
    })

    output$ind_plot <- plotly::renderPlotly({
      req(active_ind(), pk_data())
      ai <- active_ind()
      d  <- pk_data()
      df <- d$data

      if (is.null(ai$ind)) return(plotly::plot_ly())

      t_max <- max(df[[d$time_col]], na.rm=TRUE)
      t_sim <- seq(0, t_max, length.out=300)

      p <- plotly::plot_ly()

      for (i in seq_len(nrow(ai$ind))) {
        row   <- ai$ind[i, ]
        id    <- row$animal
        sub   <- df[df[[d$animal_col]] == id, ]
        dv    <- mean(sub[[d$dose_col]], na.rm=TRUE)

        params <- c(CL=row$CL, V1=row$V1)
        if (ai$n_comp == 2L) params <- c(params, V2=row$V2, Q=row$Q)

        s <- simulate_pk(params, dv, t_sim, ai$n_comp)

        p <- plotly::add_trace(p, x=sub[[d$time_col]], y=sub[[d$conc_col]],
          type="scatter", mode="markers", name=id,
          marker=list(size=9, symbol="circle-open", line=list(width=2)),
          legendgroup=id, showlegend=TRUE)

        p <- plotly::add_trace(p, x=s$time, y=s$conc,
          type="scatter", mode="lines", name=paste0(id," (fit)"),
          line=list(width=2), legendgroup=id, showlegend=FALSE)
      }

      p %>% plotly::layout(
        xaxis=list(title=paste0("Temps (",d$time_col,")")),
        yaxis=list(title=paste0("Concentration (",d$conc_col,")"), type="log"),
        hovermode="x unified"
      )
    })

    return(fit_result)
  })
}
