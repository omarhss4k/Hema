library(shiny)
library(bslib)
library(writexl)
library(DT)
library(dplyr)
library(ggplot2)

# ── UI ────────────────────────────────────────────────────────────────────────

mod_export_ui <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      card_header("Recapitulatif des resultats"),
      uiOutput(ns("export_status")),
      DT::dataTableOutput(ns("recap_table"))
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      card(
        class = "text-center p-3",
        downloadButton(
          ns("dl_excel"),
          tagList(tags$i(class = "bi bi-file-earmark-excel me-1"), "Excel (.xlsx)"),
          class = "btn-success"
        ),
        tags$small(class = "text-muted d-block mt-2",
                   "Onglets : NCA + PK_parametres")
      ),
      card(
        class = "text-center p-3",
        downloadButton(
          ns("dl_csv"),
          tagList(tags$i(class = "bi bi-filetype-csv me-1"), "CSV"),
          class = "btn-info text-white"
        ),
        tags$small(class = "text-muted d-block mt-2",
                   "NCA + parametres PK en une table")
      ),
      card(
        class = "text-center p-3",
        downloadButton(
          ns("dl_pdf"),
          tagList(tags$i(class = "bi bi-file-earmark-pdf me-1"), "Rapport PDF"),
          class = "btn-danger"
        ),
        tags$small(class = "text-muted d-block mt-2",
                   "Comparaison, parametres et graphiques")
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

mod_export_server <- function(id, nca_results, model_results, pk_data) {
  moduleServer(id, function(input, output, session) {

    # ── Active fit (best AIC) ─────────────────────────────────────────────────

    active_model <- reactive({
      m <- model_results()
      if (is.null(m)) return(NULL)
      aic1 <- m$fit1$AIC; aic2 <- m$fit2$AIC
      if (!is.null(aic1) && !is.null(aic2) && is.finite(aic1) && is.finite(aic2)) {
        if (aic1 <= aic2) list(fit=m$fit1, sim=m$sim1, n_comp=1L, label="1 compartiment")
        else               list(fit=m$fit2, sim=m$sim2, n_comp=2L, label="2 compartiments")
      } else {
        list(fit=m$fit2, sim=m$sim2, n_comp=2L, label="2 compartiments")
      }
    })

    # ── Status badges ─────────────────────────────────────────────────────────

    output$export_status <- renderUI({
      has_nca   <- !is.null(nca_results())
      has_model <- !is.null(model_results())
      tags$div(class = "mb-2 d-flex gap-3",
        tags$span(class = if (has_nca)   "badge bg-success" else "badge bg-secondary",
                  if (has_nca)   "NCA OK"    else "NCA manquant"),
        tags$span(class = if (has_model) "badge bg-success" else "badge bg-secondary",
                  if (has_model) "Modele PK OK" else "Modele PK manquant")
      )
    })

    # ── Recap data.frame ──────────────────────────────────────────────────────

    recap_df <- reactive({
      nca   <- nca_results()
      am    <- active_model()

      pk_global <- if (!is.null(am)) {
        ft <- am$fit
        cbind(
          as.data.frame(t(ft$params)),
          as.data.frame(t(unlist(ft$derived))),
          data.frame(AIC=ft$AIC, BIC=ft$BIC, RSS=ft$RSS, n_obs=as.numeric(ft$n_obs))
        )
      } else NULL

      if (!is.null(nca) && !is.null(pk_global)) {
        pk_rep <- pk_global[rep(1L, nrow(nca)), , drop=FALSE]
        row.names(pk_rep) <- NULL
        cbind(nca, pk_rep)
      } else if (!is.null(nca))       nca
        else if (!is.null(pk_global)) pk_global
        else                          NULL
    })

    output$recap_table <- DT::renderDataTable({
      df      <- req(recap_df())
      num_cols <- names(df)[sapply(df, is.numeric)]
      DT::datatable(df, rownames=FALSE,
                    options=list(scrollX=TRUE, pageLength=15, dom="tip"),
                    class="compact stripe") %>%
        DT::formatRound(columns=num_cols, digits=4)
    })

    # ── Helper: PK sheet for Excel ────────────────────────────────────────────

    .build_pk_sheet <- function(ft) {
      rse_clean    <- setNames(unname(ft$rse), gsub("^%RSE_","",names(ft$rse)))
      derived_vals <- unlist(ft$derived)
      data.frame(
        Parametre = c(names(ft$params), names(derived_vals), "AIC","BIC","RSS","n_obs"),
        Valeur    = c(unname(ft$params), unname(derived_vals),
                      ft$AIC, ft$BIC, ft$RSS, as.numeric(ft$n_obs)),
        `%RSE`    = c(unname(rse_clean[names(ft$params)]),
                      rep(NA_real_, length(derived_vals)+4)),
        check.names = FALSE
      )
    }

    # ── Helper: residuals data.frame ──────────────────────────────────────────

    .build_resid <- function(am, d) {
      df          <- d$data
      dose_by_grp <- tapply(df[[d$dose_col]], df[[d$animal_col]], mean, na.rm=TRUE)
      rows <- lapply(unique(df[[d$animal_col]]), function(id) {
        sub    <- df[df[[d$animal_col]]==id,]
        dv     <- unname(dose_by_grp[as.character(id)])
        t_obs  <- sub[[d$time_col]]
        C_obs  <- sub[[d$conc_col]]
        C_pred <- if (am$n_comp==1L)
          .pred_1comp(dv, am$fit$params["CL"], am$fit$params["V1"], t_obs)
        else
          .pred_2comp(dv, am$fit$params["CL"], am$fit$params["V1"],
                      am$fit$params["V2"], am$fit$params["Q"], t_obs)
        ok <- is.finite(C_pred) & C_pred>0 & C_obs>0
        data.frame(animal=as.character(id), time=t_obs[ok],
                   obs=C_obs[ok], pred=C_pred[ok],
                   resid=log(C_obs[ok])-log(C_pred[ok]),
                   stringsAsFactors=FALSE)
      })
      do.call(rbind, rows)
    }

    # ── Excel ─────────────────────────────────────────────────────────────────

    output$dl_excel <- downloadHandler(
      filename = function() paste0("pk_results_", format(Sys.time(),"%Y%m%d_%H%M%S"), ".xlsx"),
      content  = function(file) {
        sheets <- list()
        nca <- nca_results()
        if (!is.null(nca)) sheets[["NCA"]] <- as.data.frame(nca)
        am <- active_model()
        if (!is.null(am))  sheets[["PK_parametres"]] <- .build_pk_sheet(am$fit)
        if (length(sheets)==0)
          sheets[["Vide"]] <- data.frame(Message="Aucun resultat.")
        writexl::write_xlsx(sheets, path=file)
      }
    )

    # ── CSV ───────────────────────────────────────────────────────────────────

    output$dl_csv <- downloadHandler(
      filename = function() paste0("pk_resultats_", format(Sys.time(),"%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) {
        df <- recap_df()
        if (is.null(df)) df <- data.frame(Message="Aucun resultat.")
        write.csv(df, file, row.names=FALSE)
      }
    )

    # ── PDF ───────────────────────────────────────────────────────────────────

    output$dl_pdf <- downloadHandler(
      filename = function() paste0("rapport_pk_", format(Sys.time(),"%Y%m%d_%H%M%S"), ".pdf"),
      content  = function(file) {
        am <- req(active_model())
        d  <- req(pk_data())
        m  <- model_results()
        ft <- am$fit

        rd <- .build_resid(am, d)

        # Observed data
        df          <- d$data
        dose_by_grp <- tapply(df[[d$dose_col]], df[[d$animal_col]], mean, na.rm=TRUE)
        t_max       <- max(df[[d$time_col]], na.rm=TRUE)
        t_sim       <- seq(0, t_max, length.out=300)
        unique_doses <- sort(unique(round(unname(dose_by_grp), 8)))

        sim_df <- do.call(rbind, lapply(unique_doses, function(dv) {
          s <- simulate_pk(ft$params, dv, t_sim, am$n_comp)
          data.frame(time=s$time, conc=s$conc, dose=factor(dv))
        }))

        obs_df <- data.frame(
          animal = as.character(df[[d$animal_col]]),
          time   = df[[d$time_col]],
          conc   = df[[d$conc_col]],
          dose   = factor(round(unname(dose_by_grp[as.character(df[[d$animal_col]])]),8))
        )

        # Plot 1 : Obs vs Pred semi-log
        p1 <- ggplot() +
          geom_line(data=sim_df, aes(x=time, y=conc, color=dose), linewidth=1) +
          geom_point(data=obs_df, aes(x=time, y=conc, color=dose), size=2.5, shape=21, fill="white") +
          scale_y_log10() +
          labs(title="Ajustement : observe vs predit (semi-log)",
               x=paste0("Temps (",d$time_col,")"),
               y=paste0("Concentration (",d$conc_col,")"),
               color="Dose", subtitle=paste("Modele :", am$label)) +
          theme_bw(base_size=11)

        # Plot 2 : Residuals vs Time
        p2 <- ggplot(rd, aes(x=time, y=resid, color=animal)) +
          geom_hline(yintercept=0, color="red", linewidth=1, linetype="dashed") +
          geom_point(size=2.5) +
          labs(title="Residus log vs Temps",
               x=paste0("Temps (",d$time_col,")"),
               y="Residu log (obs - predit)", color="Animal") +
          theme_bw(base_size=11)

        # Plot 3 : Obs vs Pred scatter
        lim <- range(c(rd$obs, rd$pred), na.rm=TRUE)
        p3 <- ggplot(rd, aes(x=pred, y=obs, color=animal)) +
          geom_abline(slope=1, intercept=0, color="red", linewidth=1, linetype="dashed") +
          geom_point(size=2.5) +
          scale_x_log10() + scale_y_log10() +
          labs(title="Observe vs Predit", x="Concentration predite",
               y="Concentration observee", color="Animal") +
          theme_bw(base_size=11)

        # Parameters table as data.frame
        pk_sheet <- .build_pk_sheet(ft)

        # AIC comparison
        aic1 <- m$fit1$AIC; bic1 <- m$fit1$BIC
        aic2 <- m$fit2$AIC; bic2 <- m$fit2$BIC
        cmp_df <- data.frame(
          Modele   = c("1 compartiment","2 compartiments"),
          AIC      = c(round(aic1,1), round(aic2,1)),
          BIC      = c(round(bic1,1), round(bic2,1)),
          `n param`= c(2, 4),
          check.names=FALSE
        )

        # Write PDF
        grDevices::pdf(file, width=10, height=7, onefile=TRUE)

        # -- Page 1 : title + comparison + parameters --------------------------
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,1,
          heights=grid::unit(c(0.12,0.22,0.66),"npc"))))

        # Title
        grid::pushViewport(grid::viewport(layout.pos.row=1))
        grid::grid.text(paste0("Rapport PK — ", format(Sys.time(),"%d/%m/%Y %H:%M")),
                        gp=grid::gpar(fontsize=16, fontface="bold"))
        grid::grid.text(paste0("Modele retenu : ", am$label,
                               "   |   N observations : ", ft$n_obs),
                        y=0.3, gp=grid::gpar(fontsize=11, col="grey40"))
        grid::popViewport()

        # Comparison table
        grid::pushViewport(grid::viewport(layout.pos.row=2))
        tbl_cmp <- gridExtra::tableGrob(cmp_df, rows=NULL,
          theme=gridExtra::ttheme_default(base_size=10))
        grid::grid.draw(tbl_cmp)
        grid::popViewport()

        # Parameters table
        grid::pushViewport(grid::viewport(layout.pos.row=3))
        pk_show <- pk_sheet
        pk_show$Valeur <- signif(pk_show$Valeur, 4)
        pk_show$`%RSE` <- round(pk_show$`%RSE`, 1)
        tbl_pk <- gridExtra::tableGrob(pk_show, rows=NULL,
          theme=gridExtra::ttheme_default(base_size=9))
        grid::grid.draw(tbl_pk)
        grid::popViewport()

        # -- Page 2 : fit plot --------------------------------------------------
        print(p1)

        # -- Page 3 : residuals -------------------------------------------------
        gridExtra::grid.arrange(p2, p3, ncol=2)

        grDevices::dev.off()
      }
    )
  })
}
