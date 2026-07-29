library(shiny)
library(bslib)
library(readxl)
library(rhandsontable)
library(DT)
library(dplyr)
library(tidyr)

# ── Helpers ───────────────────────────────────────────────────────────────────

.manual_default <- data.frame(
  Animal = c("A","A","A","A","A", "B","B","B","B","B", rep(NA_character_, 6)),
  Time   = c(0, 0.5, 1, 2, 4,    0, 0.5, 1, 2, 4,    rep(NA_real_, 6)),
  Conc   = c(10, 7.5, 5.2, 2.8, 0.9, 12, 8.9, 6.1, 3.2, 1.1, rep(NA_real_, 6)),
  Dose   = c(rep(10, 5), rep(10, 5), rep(NA_real_, 6)),
  stringsAsFactors = FALSE
)

.guess_col <- function(cols, patterns) {
  hit <- grep(paste(patterns, collapse = "|"), cols,
              ignore.case = TRUE, perl = TRUE, value = TRUE)
  if (length(hit)) hit[1] else cols[1]
}

.pivot_wide_to_long <- function(df, time_col, group_cols) {
  df_long <- tidyr::pivot_longer(
    df,
    cols      = tidyr::all_of(group_cols),
    names_to  = "Groupe",
    values_to = "Conc"
  )
  names(df_long)[names(df_long) == time_col] <- "Time"

  df_long$Dose <- sapply(df_long$Groupe, function(g) {
    m <- regmatches(g, regexpr(
      "[0-9]+[.,]?[0-9]*\\s*mg/kg", g, ignore.case = TRUE))
    if (length(m) == 0 || nchar(m) == 0) return(NA_real_)
    as.numeric(gsub(",", ".", gsub("[^0-9.,]", "", m)))
  })

  df_long <- df_long[!is.na(df_long$Conc), ]
  as.data.frame(df_long)
}

.looks_wide <- function(df) {
  if (ncol(df) < 3) return(FALSE)
  n_num <- sum(sapply(df, function(x)
    is.numeric(x) || suppressWarnings(!any(is.na(as.numeric(x[!is.na(x)]))))
  ))
  n_num >= (ncol(df) - 1)
}

# ── UI ────────────────────────────────────────────────────────────────────────

mod_data_ui <- function(id) {
  ns <- NS(id)

  layout_columns(
    col_widths = c(4, 8),

    # ── Colonne gauche : contrôles ────────────────────────────────────────────
    tagList(
      card(
        card_header("Source des données"),
        radioButtons(
          ns("mode"), NULL,
          choices  = c("Importer un fichier" = "import",
                       "Saisie manuelle"      = "manual"),
          selected = "import", inline = TRUE
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] === 'import'", ns("mode")),
          fileInput(
            ns("file"),
            label       = "Fichier Excel (.xlsx / .xls) ou CSV",
            accept      = c(".csv", ".xlsx", ".xls"),
            buttonLabel = "Parcourir…",
            placeholder = "Aucun fichier sélectionné"
          )
        )
      ),

      uiOutput(ns("format_card")),
      uiOutput(ns("mapping_card")),

      actionButton(
        ns("validate"),
        label = tagList(tags$i(class = "bi bi-check2-circle me-1"), "Valider les données"),
        class = "btn-primary w-100 mt-2"
      ),
      uiOutput(ns("validation_panel"))
    ),

    # ── Colonne droite : tableau de données ───────────────────────────────────
    tagList(
      conditionalPanel(
        condition = sprintf("input['%s'] === 'manual'", ns("mode")),
        card(
          card_header(tagList(
            tags$i(class = "bi bi-pencil-square me-1"),
            "Saisie manuelle",
            tags$span(class = "text-muted small ms-3",
              "Double-cliquez pour éditer · Clic droit : ajouter / supprimer des lignes")
          )),
          rhandsontable::rHandsontableOutput(ns("manual_table"), height = "620px")
        )
      ),
      uiOutput(ns("editor_card"))
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

mod_data_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── 1. Raw file data ──────────────────────────────────────────────────────

    raw_data <- reactive({
      req(input$mode == "import", input$file)
      ext <- tolower(tools::file_ext(input$file$name))
      df <- tryCatch({
        if (ext %in% c("xlsx", "xls"))
          readxl::read_excel(input$file$datapath)
        else
          read.csv(input$file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      }, error = function(e)
        shiny::validate(shiny::need(FALSE,
          paste0("Impossible de lire le fichier : ", conditionMessage(e))))
      )
      as.data.frame(df)
    })

    # ── 2. Format card (appears when file loaded) ─────────────────────────────

    output$format_card <- renderUI({
      req(input$mode == "import")
      df   <- req(raw_data())
      cols <- names(df)

      card(
        card_header("Format des données"),
        radioButtons(
          ns("data_format"), NULL,
          choices = c(
            "Format long  — une ligne = une mesure  (colonnes : animal, temps, concentration)" = "long",
            "Format large — une colonne = un groupe  (tableau croisé temps × groupes)" = "wide"
          ),
          selected = if (.looks_wide(df)) "wide" else "long"
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] === 'wide'", ns("data_format")),
          tags$div(
            class = "alert alert-warning py-2 mb-0",
            tags$i(class = "bi bi-shuffle me-1"),
            "Format large détecté. Sélectionnez la colonne Temps ; les autres colonnes ",
            "deviendront la colonne Groupe. La dose sera extraite du nom de colonne si ",
            "elle contient ", tags$code("mg/kg"), "."
          ),
          selectInput(
            ns("pivot_time_col"), "Colonne Temps",
            choices  = cols,
            selected = .guess_col(cols, c("^time", "^t\\s", "temps", "hour", "heure", "^h$"))
          )
        )
      )
    })

    # ── 3. Base data (file + optional pivot, recalculated on source change) ───

    base_data <- reactive({
      req(input$mode == "import")
      df  <- req(raw_data())
      fmt <- req(input$data_format)

      if (fmt == "wide") {
        req(input$pivot_time_col)
        gc <- setdiff(names(df), input$pivot_time_col)
        shiny::validate(
          shiny::need(length(gc) >= 1, "Au moins une colonne groupe est requise."))
        tryCatch(
          .pivot_wide_to_long(df, input$pivot_time_col, gc),
          error = function(e) {
            showNotification(paste0("Erreur pivot : ", conditionMessage(e)),
                             type = "error", duration = 12)
            df
          }
        )
      } else {
        df
      }
    })

    # ── 4. editable_data : reactiveVal, reset on source, updated on edit ──────

    editable_data <- reactiveVal(NULL)

    # Reset whenever file, format, or pivot column changes
    observeEvent(base_data(), {
      editable_data(base_data())
    }, ignoreNULL = TRUE)

    # Update from user edits in the editor rhandsontable
    observeEvent(input$edit_table, {
      req(!is.null(input$edit_table))
      editable_data(rhandsontable::hot_to_r(input$edit_table))
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # ── 5. Editor card ────────────────────────────────────────────────────────

    output$editor_card <- renderUI({
      req(input$mode == "import")
      req(editable_data())

      df      <- editable_data()
      n_rows  <- nrow(df)
      n_cols  <- ncol(df)
      n_grp   <- if ("Groupe" %in% names(df))
        length(unique(df$Groupe)) else "—"

      card(
        card_header(tagList(
          tags$i(class = "bi bi-table me-1"),
          paste0("Données — ", n_rows, " lignes · ", n_cols, " colonnes",
                 if (is.numeric(n_grp)) paste0(" · ", n_grp, " groupe(s)") else ""),
          tags$span(class = "text-muted small ms-3",
            "Double-cliquez pour modifier · Clic droit : insérer / supprimer une ligne")
        )),
        rhandsontable::rHandsontableOutput(ns("edit_table"), height = "620px")
      )
    })

    output$edit_table <- rhandsontable::renderRHandsontable({
      df <- req(editable_data())

      rhandsontable::rhandsontable(
        df,
        rowHeaders  = FALSE,
        useTypes    = TRUE,
        stretchH    = "all",
        contextMenu = TRUE,
        height      = 600
      )
    })

    # ── 6. Manual rhandsontable ───────────────────────────────────────────────

    output$manual_table <- rhandsontable::renderRHandsontable({
      rhandsontable::rhandsontable(
        .manual_default, rowHeaders = TRUE,
        useTypes = TRUE, stretchH = "all", contextMenu = TRUE,
        height = 600
      ) %>%
        rhandsontable::hot_col("Animal", type = "text") %>%
        rhandsontable::hot_col("Time",   type = "numeric", format = "0.000") %>%
        rhandsontable::hot_col("Conc",   type = "numeric", format = "0.0000") %>%
        rhandsontable::hot_col("Dose",   type = "numeric", format = "0.00")
    })

    # ── 7. working_data : editable (import) or manual ─────────────────────────

    working_data <- reactive({
      if (input$mode == "import") {
        req(editable_data())
      } else {
        req(!is.null(input$manual_table))
        rhandsontable::hot_to_r(input$manual_table)
      }
    })

    # ── 8. Mapping card ───────────────────────────────────────────────────────

    output$mapping_card <- renderUI({
      df   <- req(working_data())
      cols <- names(df)

      card(
        card_header("Mapping des colonnes"),
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          selectInput(ns("col_time"),   "Temps",
            choices = cols,
            selected = .guess_col(cols, c("^Time$","^time$","^t$","temps","hour","heure","hr"))
          ),
          selectInput(ns("col_conc"),   "Concentration",
            choices = cols,
            selected = .guess_col(cols, c("^Conc$","conc","cp","plasma","concentration","dv"))
          ),
          selectInput(ns("col_animal"), "Animal / Groupe",
            choices = cols,
            selected = .guess_col(cols, c("^Groupe$","^Group","animal","subject","sujet","^id$"))
          ),
          selectInput(ns("col_dose"),   "Dose",
            choices = cols,
            selected = .guess_col(cols, c("^Dose$","^dose$","^amt$","amount"))
          )
        )
      )
    })

    # ── 9. Validated reactive ─────────────────────────────────────────────────

    pk_data <- eventReactive(input$validate, {
      df         <- req(working_data())
      time_col   <- req(input$col_time)
      conc_col   <- req(input$col_conc)
      animal_col <- req(input$col_animal)
      dose_col   <- req(input$col_dose)

      df[[time_col]] <- suppressWarnings(as.numeric(df[[time_col]]))
      df[[conc_col]] <- suppressWarnings(as.numeric(df[[conc_col]]))
      df[[dose_col]] <- suppressWarnings(as.numeric(df[[dose_col]]))

      shiny::validate(
        shiny::need(!all(is.na(df[[time_col]])),
          paste0("Colonne Temps « ", time_col, " » : aucune valeur numérique.")),
        shiny::need(!all(is.na(df[[conc_col]])),
          paste0("Colonne Concentration « ", conc_col, " » : aucune valeur numérique.")),
        shiny::need(!all(is.na(df[[dose_col]])),
          paste0("Colonne Dose « ", dose_col,
                 " » : aucune valeur numérique. ",
                 "Renseignez la dose directement dans l'éditeur ci-dessus."))
      )

      n_before <- nrow(df)
      df <- df[!is.na(df[[time_col]]) & !is.na(df[[conc_col]]) &
               !is.na(df[[animal_col]]) & !is.na(df[[dose_col]]), ]
      if ((n_before - nrow(df)) > 0)
        showNotification(paste0(n_before - nrow(df),
          " ligne(s) avec valeurs manquantes supprimée(s)."),
          type = "warning", duration = 10)

      n_before2 <- nrow(df)
      df <- df[df[[conc_col]] > 0, ]
      if ((n_before2 - nrow(df)) > 0)
        showNotification(paste0(n_before2 - nrow(df),
          " observation(s) avec concentration ≤ 0 exclue(s)."),
          type = "warning", duration = 10)

      shiny::validate(
        shiny::need(nrow(df) > 0, "Aucune observation valide après nettoyage."))

      pts      <- tapply(df[[time_col]], df[[animal_col]], length)
      bad_anim <- names(pts[pts < 3])
      shiny::validate(
        shiny::need(length(bad_anim) == 0,
          paste0("Groupe(s) avec moins de 3 points : ",
                 paste(bad_anim, collapse = ", "), ".")))

      list(data = df, time_col = time_col, conc_col = conc_col,
           animal_col = animal_col, dose_col = dose_col)
    })

    # ── 10. Validation panel ──────────────────────────────────────────────────

    output$validation_panel <- renderUI({
      req(input$validate > 0)
      result <- tryCatch(pk_data(), error = function(e) e)

      if (inherits(result, "error")) {
        card(
          class = "mt-3 border-danger",
          card_header(class = "bg-danger text-white",
            tagList(tags$i(class = "bi bi-exclamation-triangle-fill me-2"),
                    "Erreur de validation")),
          tags$p(class = "text-danger mb-0 fw-semibold", result$message)
        )
      } else {
        df        <- result$data
        n_animals <- length(unique(df[[result$animal_col]]))
        t_range   <- range(df[[result$time_col]], na.rm = TRUE)
        c_range   <- range(df[[result$conc_col]], na.rm = TRUE)

        card(
          class = "mt-3 border-success",
          card_header(class = "bg-success text-white",
            tagList(tags$i(class = "bi bi-check-circle-fill me-2"),
                    "Données validées")),
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(title = "Observations",      value = nrow(df),    theme = "success"),
            value_box(title = "Groupes / Animaux", value = n_animals,   theme = "success"),
            value_box(title = "Plage de temps",
                      value = paste0(round(t_range[1], 1), " – ", round(t_range[2], 1)),
                      theme = "success"),
            value_box(title = "Cmax globale",
                      value = round(c_range[2], 2), theme = "success")
          )
        )
      }
    })

    return(pk_data)
  })
}
