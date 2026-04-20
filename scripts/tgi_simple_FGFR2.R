# TGI simple — données observées, pas de modèle
# TGI(%) = (1 - ΔW_traité / ΔW_contrôle) × 100
#   ΔW = W_fin - W_j0

library(ggplot2)
library(readxl)

# ── Lecture données ────────────────────────────────────────────────────────────
raw <- read_xlsx("data/TumorVolume_FGFR2.xlsx",
                 n_max     = 3,
                 col_types = c("text", rep("numeric", 15)))

time_d <- suppressWarnings(as.numeric(colnames(raw)[-1]))

get_group <- function(pattern)
  as.numeric(unlist(raw[grep(pattern, raw[[1]], ignore.case=TRUE)[1], -1]))

w_ctrl <- get_group("Group 01")
w_d10  <- get_group("Group 03")   # 10 mg/kg
w_d3   <- get_group("Group 04")   #  3 mg/kg

# ── w0 : valeur initiale commune (j0) ─────────────────────────────────────────
w0_ctrl <- w_ctrl[time_d == 0]
w0_d10  <- w_d10 [time_d == 0]
w0_d3   <- w_d3  [time_d == 0]

# ── TGI à chaque point de mesure ─────────────────────────────────────────────
tgi <- function(w_trt, w_trt0, w_c, w_c0)
  round((1 - (w_trt - w_trt0) / (w_c - w_c0)) * 100, 1)

ok <- !is.na(w_ctrl) & !is.na(w_d10) & !is.na(w_d3) & time_d > 0

df_tgi <- data.frame(
  jour    = time_d[ok],
  TGI_3   = tgi(w_d3[ok],  w0_d3,  w_ctrl[ok], w0_ctrl),
  TGI_10  = tgi(w_d10[ok], w0_d10, w_ctrl[ok], w0_ctrl)
)

cat("=== TGI observé (%) ===\n")
cat(sprintf("%-6s  %10s  %10s\n", "Jour", "3 mg/kg", "10 mg/kg"))
cat(strrep("-", 30), "\n")
for (i in seq_len(nrow(df_tgi)))
  cat(sprintf("j%-5d  %9.1f%%  %9.1f%%\n",
              df_tgi$jour[i], df_tgi$TGI_3[i], df_tgi$TGI_10[i]))

cat("\n--- Dernier point (j", tail(df_tgi$jour, 1), ") ---\n", sep="")
cat("   3 mg/kg  :", tail(df_tgi$TGI_3,  1), "%\n")
cat("  10 mg/kg  :", tail(df_tgi$TGI_10, 1), "%\n")

# ── Graphique TGI au cours du temps ──────────────────────────────────────────
df_long <- rbind(
  data.frame(jour=df_tgi$jour, TGI=df_tgi$TGI_3,  Dose="3 mg/kg"),
  data.frame(jour=df_tgi$jour, TGI=df_tgi$TGI_10, Dose="10 mg/kg")
)
df_long$Dose <- factor(df_long$Dose, levels=c("3 mg/kg","10 mg/kg"))

ggplot(df_long, aes(x=jour, y=TGI, color=Dose, group=Dose)) +
  geom_line(linewidth=1) +
  geom_point(size=2.5) +
  geom_hline(yintercept=c(0, 100), linetype="dashed", color="grey50") +
  labs(title = "TGI observé — Fc-silent FGFR2-huBPA-LP1",
       subtitle = "TGI(%) = (1 − ΔW_traité / ΔW_contrôle) × 100",
       x = "Temps (jours)", y = "TGI (%)") +
  theme_bw(base_size=13)

ggsave("scripts/plot_TGI_simple_FGFR2.png", width=8, height=5, dpi=150)
cat("\nGraphique → scripts/plot_TGI_simple_FGFR2.png\n")
