# =============================================================================
# Tableau TGI observé — Fc-silent FGFR2-huBPA-LP1
# Sans gt — ggplot2 uniquement
# =============================================================================

library(readxl)
library(ggplot2)

# =============================================================================
# 1. LECTURE (deux blocs : moyennes + SEM)
# =============================================================================

raw_all <- suppressMessages(
  read_xlsx("TumorVolume_FGFR2.xlsx",
            col_names = FALSE, .name_repair = "minimal")
)
raw_all <- raw_all[, colSums(!is.na(raw_all)) > 0]
header_idx <- which(raw_all[[1]] == "Group")

times_all <- suppressWarnings(as.numeric(as.character(unlist(raw_all[header_idx[1], -1]))))
valid  <- !is.na(times_all)
times_d <- times_all[valid]

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

w_ctrl   <- get_row(blk_mean, "Group 01")
w_d10    <- get_row(blk_mean, "Group 03")
w_d3     <- get_row(blk_mean, "Group 04")
sem_ctrl <- get_row(blk_sem,  "Group 01")
sem_d10  <- get_row(blk_sem,  "Group 03")
sem_d3   <- get_row(blk_sem,  "Group 04")

# =============================================================================
# 2. CALCUL TGI
# =============================================================================

w0_ctrl <- w_ctrl[times_d == 0]
w0_d10  <- w_d10 [times_d == 0]
w0_d3   <- w_d3  [times_d == 0]
ok <- !is.na(w_ctrl) & !is.na(w_d10) & !is.na(w_d3) & times_d > 0

tgi_pct <- function(wt, wt0, wc, wc0)
  round((1 - (wt - wt0) / (wc - wc0)) * 100, 1)

df_tgi <- data.frame(
  jour      = times_d[ok],
  tgi_3     = tgi_pct(w_d3[ok],  w0_d3,  w_ctrl[ok], w0_ctrl),
  tgi_10    = tgi_pct(w_d10[ok], w0_d10, w_ctrl[ok], w0_ctrl),
  mean_ctrl = round(w_ctrl[ok],    0),
  sem_ctrl  = round(sem_ctrl[ok],  0),
  mean_3    = round(w_d3[ok],      0),
  sem_3     = round(sem_d3[ok],    0),
  mean_10   = round(w_d10[ok],     0),
  sem_10    = round(sem_d10[ok],   0)
)

# =============================================================================
# 3. RÉSUMÉ CONSOLE
# =============================================================================

cat(sprintf("\n%-6s  %12s  %14s  %12s  %14s\n",
            "Temps",
            "Ctrl (mm³)", "3 mg/kg (mm³)",
            "TGI 3mg%",   "TGI 10mg%"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(df_tgi))) {
  cat(sprintf("j%-5d  %5.0f±%-5.0f  %5.0f±%-6.0f  %8.1f %%  %10.1f %%\n",
              df_tgi$jour[i],
              df_tgi$mean_ctrl[i], df_tgi$sem_ctrl[i],
              df_tgi$mean_3[i],    df_tgi$sem_3[i],
              df_tgi$tgi_3[i],     df_tgi$tgi_10[i]))
}
cat(strrep("-", 65), "\n")
cat(sprintf("Dernier point (j%d)  →  3 mg/kg : %.1f %%  |  10 mg/kg : %.1f %%\n\n",
            tail(df_tgi$jour,1), tail(df_tgi$tgi_3,1), tail(df_tgi$tgi_10,1)))

# =============================================================================
# 4. TABLEAU VISUEL (ggplot2)
# =============================================================================

# Construire les cellules
col_names <- c("Temps",
               "Contrôle\nmoy (mm³)", "Contrôle\nSEM",
               "3 mg/kg\nmoy (mm³)",  "3 mg/kg\nSEM",  "TGI 3mg\n(%)",
               "10 mg/kg\nmoy (mm³)", "10 mg/kg\nSEM", "TGI 10mg\n(%)")
n_col <- length(col_names)
n_row <- nrow(df_tgi)

# Ligne d'en-tête
header <- data.frame(
  row   = 0,
  col   = seq_len(n_col),
  label = col_names,
  fill  = "#2c3e50",
  tcol  = "white",
  bold  = TRUE
)

# Lignes de données
rows_list <- lapply(seq_len(n_row), function(i) {
  tgi3  <- df_tgi$tgi_3[i]
  tgi10 <- df_tgi$tgi_10[i]
  bg    <- if (i %% 2 == 0) "#f2f2f2" else "white"

  fill_tgi3  <- if (tgi3  >= 60) "#c3e6cb" else if (tgi3  >= 30) "#ffeeba" else bg
  fill_tgi10 <- if (tgi10 >= 60) "#c3e6cb" else if (tgi10 >= 30) "#ffeeba" else bg
  is_last    <- i == n_row

  data.frame(
    row   = i,
    col   = seq_len(n_col),
    label = c(paste0("j", df_tgi$jour[i]),
              df_tgi$mean_ctrl[i], df_tgi$sem_ctrl[i],
              df_tgi$mean_3[i],    df_tgi$sem_3[i],
              sprintf("%.1f %%", tgi3),
              df_tgi$mean_10[i],   df_tgi$sem_10[i],
              sprintf("%.1f %%", tgi10)),
    fill  = c(rep(bg, 5), fill_tgi3, rep(bg, 2), fill_tgi10),
    tcol  = "black",
    bold  = is_last
  )
})

cells <- do.call(rbind, c(list(header), rows_list))
cells$y    <- max(cells$row) - cells$row   # inverser : en-tête en haut
cells$face <- ifelse(cells$bold, "bold", "plain")

p_tbl <- ggplot(cells, aes(x = col, y = y)) +
  geom_tile(aes(fill = fill), color = "grey70", linewidth = 0.3) +
  geom_text(aes(label = label, color = tcol, fontface = face),
            size = 3.0, lineheight = 0.95) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.05, 0)) +
  labs(
    title    = "TGI observé — Fc-silent FGFR2-huBPA-LP1",
    subtitle = "Q2W (j0/j14/j28/j42) | Vert : TGI ≥ 60 %  ·  Jaune : 30–60 %  ·  Dernier point en gras"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey40", margin = margin(b = 6)),
    plot.margin   = margin(12, 12, 12, 12)
  )

ggsave("scripts/tableau_TGI_FGFR2.png", p_tbl,
       width = 11, height = 0.55 * (n_row + 2) + 1.5, dpi = 150)
cat("Tableau → scripts/tableau_TGI_FGFR2.png\n")

# =============================================================================
# 5. EXPORT CSV
# =============================================================================

write.csv(
  data.frame(
    Temps         = paste0("j", df_tgi$jour),
    Ctrl_moy_mm3  = df_tgi$mean_ctrl,
    Ctrl_SEM      = df_tgi$sem_ctrl,
    D3_moy_mm3    = df_tgi$mean_3,
    D3_SEM        = df_tgi$sem_3,
    TGI_3mg_pct   = df_tgi$tgi_3,
    D10_moy_mm3   = df_tgi$mean_10,
    D10_SEM       = df_tgi$sem_10,
    TGI_10mg_pct  = df_tgi$tgi_10
  ),
  "scripts/tableau_TGI_FGFR2.csv", row.names = FALSE
)
cat("CSV     → scripts/tableau_TGI_FGFR2.csv\n")
