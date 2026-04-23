# =============================================================================
# Tableau TGI observé — Fc-silent FGFR2-huBPA-LP1
# Adapté depuis tgi_simple_FGFR2.R
# =============================================================================

library(readxl)
library(ggplot2)
library(gt)          # install.packages("gt") si nécessaire

# =============================================================================
# 1. LECTURE DES DONNÉES (deux blocs : moyennes + SEM)
# =============================================================================

raw_all <- suppressMessages(
  read_xlsx("TumorVolume_FGFR2.xlsx",
            col_names = FALSE, .name_repair = "minimal")
)
raw_all <- raw_all[, colSums(!is.na(raw_all)) > 0]

header_idx <- which(raw_all[[1]] == "Group")

times_all <- suppressWarnings(as.numeric(as.character(unlist(raw_all[header_idx[1], -1]))))
valid     <- !is.na(times_all)
times_d   <- times_all[valid]

extract_block <- function(start, end) {
  blk <- raw_all[start:end, ]
  blk <- blk[!is.na(blk[[1]]) & blk[[1]] != "", ]
  dat <- blk[, c(TRUE, valid)]
  for (j in 2:ncol(dat)) storage.mode(dat[[j]]) <- "numeric"
  dat
}

blk_mean <- extract_block(header_idx[1] + 1, header_idx[2] - 1)
blk_sem  <- extract_block(header_idx[2] + 1, nrow(raw_all))

get_row <- function(blk, pat)
  as.numeric(unlist(blk[grep(pat, blk[[1]], ignore.case = TRUE)[1], -1]))

# Moyennes (mm³)
w_ctrl <- get_row(blk_mean, "Group 01")
w_d10  <- get_row(blk_mean, "Group 03")
w_d3   <- get_row(blk_mean, "Group 04")

# SEM (mm³)
sem_ctrl <- get_row(blk_sem, "Group 01")
sem_d10  <- get_row(blk_sem, "Group 03")
sem_d3   <- get_row(blk_sem, "Group 04")

# =============================================================================
# 2. CALCUL TGI OBSERVÉ
# =============================================================================

w0_ctrl <- w_ctrl[times_d == 0]
w0_d10  <- w_d10 [times_d == 0]
w0_d3   <- w_d3  [times_d == 0]

ok <- !is.na(w_ctrl) & !is.na(w_d10) & !is.na(w_d3) & times_d > 0

tgi_pct <- function(wt, wt0, wc, wc0)
  round((1 - (wt - wt0) / (wc - wc0)) * 100, 1)

df_tgi <- data.frame(
  jour   = times_d[ok],
  tgi_3  = tgi_pct(w_d3[ok],  w0_d3,  w_ctrl[ok], w0_ctrl),
  tgi_10 = tgi_pct(w_d10[ok], w0_d10, w_ctrl[ok], w0_ctrl),
  mean_ctrl = w_ctrl[ok],
  mean_3    = w_d3[ok],
  mean_10   = w_d10[ok],
  sem_ctrl  = sem_ctrl[ok],
  sem_3     = sem_d3[ok],
  sem_10    = sem_d10[ok]
)

# =============================================================================
# 3. TABLEAU GT
# =============================================================================

# Interprétation du TGI (seuils standards précliniques)
interprete <- function(tgi) {
  dplyr::case_when(
    tgi >= 100            ~ "Régression",
    tgi >=  60            ~ "Réponse marquée",
    tgi >=  30            ~ "Réponse modérée",
    tgi >=   0            ~ "Réponse faible",
    TRUE                  ~ "Progression"
  )
}

df_table <- data.frame(
  Temps       = paste0("j", df_tgi$jour),
  ctrl_moy    = round(df_tgi$mean_ctrl, 0),
  ctrl_sem    = round(df_tgi$sem_ctrl,  0),
  d3_moy      = round(df_tgi$mean_3,   0),
  d3_sem      = round(df_tgi$sem_3,    0),
  tgi_3       = df_tgi$tgi_3,
  interp_3    = interprete(df_tgi$tgi_3),
  d10_moy     = round(df_tgi$mean_10,  0),
  d10_sem     = round(df_tgi$sem_10,   0),
  tgi_10      = df_tgi$tgi_10,
  interp_10   = interprete(df_tgi$tgi_10)
)

tbl <- gt(df_table) |>
  tab_header(
    title    = md("**TGI observé — Fc-silent FGFR2-huBPA-LP1**"),
    subtitle = "Doses répétées Q2W (j0, j14, j28, j42)"
  ) |>
  tab_spanner(label = "Contrôle",  columns = c(ctrl_moy, ctrl_sem)) |>
  tab_spanner(label = "3 mg/kg",   columns = c(d3_moy,  d3_sem,  tgi_3,  interp_3))  |>
  tab_spanner(label = "10 mg/kg",  columns = c(d10_moy, d10_sem, tgi_10, interp_10)) |>
  cols_label(
    Temps     = "Temps",
    ctrl_moy  = "Moy (mm³)", ctrl_sem  = "SEM",
    d3_moy    = "Moy (mm³)", d3_sem    = "SEM",
    tgi_3     = "TGI (%)",   interp_3  = "Réponse",
    d10_moy   = "Moy (mm³)", d10_sem   = "SEM",
    tgi_10    = "TGI (%)",   interp_10 = "Réponse"
  ) |>
  fmt_number(columns = c(tgi_3, tgi_10), decimals = 1, pattern = "{x} %") |>
  # Couleur TGI par seuil
  tab_style(
    style     = cell_fill(color = "#d4edda"),
    locations = cells_body(columns = c(tgi_3, tgi_10),
                           rows    = tgi_3 >= 60 | tgi_10 >= 60)
  ) |>
  tab_style(
    style     = cell_fill(color = "#fff3cd"),
    locations = cells_body(columns = c(tgi_3, tgi_10),
                           rows    = (tgi_3 >= 30 & tgi_3 < 60) |
                                     (tgi_10 >= 30 & tgi_10 < 60))
  ) |>
  # Mettre en gras le dernier timepoint
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = nrow(df_table))
  ) |>
  tab_source_note("Vert : réponse marquée (TGI ≥ 60 %)  |  Jaune : réponse modérée (30–60 %)")  |>
  tab_options(
    heading.background.color  = "#2c3e50",
    heading.title.font.size   = 14,
    column_labels.font.weight = "bold",
    table.font.size           = 12
  )

print(tbl)
gtsave(tbl, "scripts/tableau_TGI_FGFR2.png")
cat("Tableau → scripts/tableau_TGI_FGFR2.png\n")

# =============================================================================
# 4. EXPORT CSV
# =============================================================================

df_csv <- data.frame(
  Temps        = paste0("j", df_tgi$jour),
  `TGI_3mg_pct`  = df_tgi$tgi_3,
  `TGI_10mg_pct` = df_tgi$tgi_10,
  `Ctrl_mean_mm3`  = round(df_tgi$mean_ctrl, 1),
  `Ctrl_SEM_mm3`   = round(df_tgi$sem_ctrl,  1),
  `D3_mean_mm3`    = round(df_tgi$mean_3,    1),
  `D3_SEM_mm3`     = round(df_tgi$sem_3,     1),
  `D10_mean_mm3`   = round(df_tgi$mean_10,   1),
  `D10_SEM_mm3`    = round(df_tgi$sem_10,    1),
  check.names = FALSE
)
write.csv(df_csv, "scripts/tableau_TGI_FGFR2.csv", row.names = FALSE)
cat("CSV     → scripts/tableau_TGI_FGFR2.csv\n")

# =============================================================================
# 5. RÉSUMÉ CONSOLE
# =============================================================================

j_fin  <- tail(df_tgi$jour, 1)
cat(sprintf("\n=== TGI au dernier timepoint (j%d) ===\n", j_fin))
cat(sprintf("   3 mg/kg  : %5.1f %%\n",   tail(df_tgi$tgi_3,  1)))
cat(sprintf("  10 mg/kg  : %5.1f %%\n\n", tail(df_tgi$tgi_10, 1)))
