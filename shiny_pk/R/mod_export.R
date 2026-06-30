library(shiny)
library(bslib)
library(writexl)
library(DT)
library(dplyr)

# ── UI ────────────────────────────────────────────────────────────────────────

mod_export_ui <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      card_header("Récapitulatif des résultats"),
      uiOutput(ns("export_status")),
      DT::dataTableOutput(ns("recap_table"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(
        class = "text-center p-3",
        downloadButton(
          ns("dl_excel"),
          tagList(
            tags$i(class = "bi bi-file-earmark-excel me-1"),
            "Télécharger Excel (.xlsx)"
          ),
          class = "btn-success"
        ),
        tags$small(class = "text-muted d-block mt-2",
                   "Deux onglets : « NCA » et « PK_parametres »")
      ),
      card(
        class = "text-center p-3",
        downloadButton(
          ns("dl_csv"),
          tagList(
            tags$i(class = "bi bi-filetype-csv me-1"),
            "Télécharger CSV"
          ),
          class = "btn-info text-white"
        ),
        tags$small(class = "text-muted d-block mt-2",
                   "Résultats NCA + paramètres PK en une table")
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

mod_export_server <- function(id, nca_results, model_results) {
  moduleServer(id, function(input, output, session) {

    # ── Status badges ─────────────────────────────────────────────────────────

    output$export_status <- renderUI({
      has_nca   <- !is.null(nca_results())
      has_model <- !is.null(model_results()) && isTRUE(model_results()$ok)

      tags$div(
        class = "mb-2 d-flex gap-3",
        tags$span(
          class = if (has_nca) "badge bg-success" else "badge bg-secondary",
          if (has_nca) "NCA ✓" else "NCA manquant"
        ),
        tags$span(
          class = if (has_model) "badge bg-success" else "badge bg-secondary",
          if (has_model) "Modèle PK ✓" else "Modèle PK manquant"
        )
      )
    })

    # ── Combined recap data.frame ─────────────────────────────────────────────

    recap_df <- reactive({
      nca   <- nca_results()
      model <- model_results()

      has_model <- !is.null(model) && isTRUE(model$ok)

      pk_global <- if (has_model) {
        ft          <- model$fit
        params_row  <- as.data.frame(t(ft$params))
        derived_row <- as.data.frame(t(unlist(ft$derived)))
        crit_row    <- data.frame(
          AIC   = ft$AIC,
          BIC   = ft$BIC,
          RSS   = ft$RSS,
          n_obs = as.numeric(ft$n_obs)
        )
        cbind(params_row, derived_row, crit_row)
      } else NULL

      if (!is.null(nca) && !is.null(pk_global)) {
        pk_rep <- pk_global[rep(1L, nrow(nca)), , drop = FALSE]
        row.names(pk_rep) <- NULL
        cbind(nca, pk_rep)
      } else if (!is.null(nca)) {
        nca
      } else if (!is.null(pk_global)) {
        pk_global
      } else {
        NULL
      }
    })

    # ── Recap DT table ────────────────────────────────────────────────────────

    output$recap_table <- DT::renderDataTable({
      df <- req(recap_df())
      num_cols <- names(df)[sapply(df, is.numeric)]

      DT::datatable(
        df,
        rownames = FALSE,
        options  = list(scrollX = TRUE, pageLength = 15, dom = "tip"),
        class    = "compact stripe"
      ) %>%
        DT::formatRound(columns = num_cols, digits = 4)
    })

    # ── Internal helper : build PK params data.frame for Excel ───────────────

    .build_pk_sheet <- function(ft) {
      rse_clean <- setNames(
        unname(ft$rse),
        gsub("^%RSE_", "", names(ft$rse))
      )
      derived_vals <- unlist(ft$derived)

      data.frame(
        Parametre = c(
          names(ft$params),
          names(derived_vals),
          "AIC", "BIC", "RSS", "n_obs"
        ),
        Valeur = c(
          unname(ft$params),
          unname(derived_vals),
          ft$AIC, ft$BIC, ft$RSS, as.numeric(ft$n_obs)
        ),
        `%RSE` = c(
          unname(rse_clean[names(ft$params)]),
          rep(NA_real_, length(derived_vals) + 4)
        ),
        check.names = FALSE
      )
    }

    # ── Excel download ─────────────────────────────────────────────────────────

    output$dl_excel <- downloadHandler(
      filename = function() {
        paste0("pk_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
      },
      content = function(file) {
        sheets <- list()

        nca <- nca_results()
        if (!is.null(nca))
          sheets[["NCA"]] <- as.data.frame(nca)

        model <- model_results()
        if (!is.null(model) && isTRUE(model$ok))
          sheets[["PK_parametres"]] <- .build_pk_sheet(model$fit)

        if (length(sheets) == 0)
          sheets[["Vide"]] <- data.frame(
            Message = "Aucun resultat. Validez les donnees et lancez les analyses."
          )

        writexl::write_xlsx(sheets, path = file)
      }
    )

    # ── CSV download ───────────────────────────────────────────────────────────

    output$dl_csv <- downloadHandler(
      filename = function() {
        paste0("pk_resultats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        df <- recap_df()
        if (is.null(df))
          df <- data.frame(Message = "Aucun resultat disponible.")
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}
