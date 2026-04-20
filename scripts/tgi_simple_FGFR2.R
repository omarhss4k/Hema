# TGI + profil PK — Fc-silent FGFR2-huBPA-LP1
# TGI calculé sur données observées | PK simulé (2-comp, doses répétées)

library(deSolve)
library(ggplot2)
library(readxl)
library(patchwork)

# ── PK 2-comp ─────────────────────────────────────────────────────────────────
load("scripts/resultats_PK2comp_FGFR2.RData")

CL <- unname(pk2comp$CL);  V1 <- unname(pk2comp$V1)
V2 <- unname(pk2comp$V2);  Q  <- unname(pk2comp$Q)
k10 <- CL/V1;  k12 <- Q/V1;  k21 <- Q/V2

pk_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    list(c(-(k10+k12)*A1 + k21*A2,
             k12*A1 - k21*A2))
  })
}

sim_pk <- function(dose, times, doses_h=0) {
  state0 <- c(A1=dose, A2=0)
  pars   <- c(k10=k10, k12=k12, k21=k21)
  events <- if (length(doses_h) > 1)
    data.frame(var="A1", time=doses_h[-1], value=dose, method="add")
  else NULL
  out <- as.data.frame(lsoda(state0, times, pk_ode, pars,
           events=if (!is.null(events)) list(data=events) else NULL,
           rtol=1e-6, atol=1e-8))
  out$C1 <- out$A1 / V1
  out
}

# ── Données tumorales ─────────────────────────────────────────────────────────
raw <- read_xlsx("data/TumorVolume_FGFR2.xlsx",
                 n_max=3, col_types=c("text", rep("numeric", 15)))

time_d <- suppressWarnings(as.numeric(colnames(raw)[-1]))

get_w <- function(pat)
  as.numeric(unlist(raw[grep(pat, raw[[1]], ignore.case=TRUE)[1], -1])) / 1000

w_ctrl <- get_w("Group 01")
w_d10  <- get_w("Group 03")
w_d3   <- get_w("Group 04")

# ── TGI observé ───────────────────────────────────────────────────────────────
w0_ctrl <- w_ctrl[time_d == 0]
w0_d10  <- w_d10 [time_d == 0]
w0_d3   <- w_d3  [time_d == 0]

ok <- !is.na(w_ctrl) & !is.na(w_d10) & !is.na(w_d3) & time_d > 0

tgi_pct <- function(wt, wt0, wc, wc0)
  round((1 - (wt - wt0) / (wc - wc0)) * 100, 1)

df_tgi <- data.frame(
  jour   = time_d[ok],
  TGI_3  = tgi_pct(w_d3[ok],  w0_d3,  w_ctrl[ok], w0_ctrl),
  TGI_10 = tgi_pct(w_d10[ok], w0_d10, w_ctrl[ok], w0_ctrl)
)

cat("=== TGI observé (%) ===\n")
cat(sprintf("%-6s  %10s  %10s\n", "Jour", "3 mg/kg", "10 mg/kg"))
cat(strrep("-", 30), "\n")
for (i in seq_len(nrow(df_tgi)))
  cat(sprintf("j%-5d  %9.1f%%  %9.1f%%\n",
              df_tgi$jour[i], df_tgi$TGI_3[i], df_tgi$TGI_10[i]))
cat("\n--- Dernier point (j", tail(df_tgi$jour,1), ") ---\n", sep="")
cat("   3 mg/kg  :", tail(df_tgi$TGI_3,  1), "%\n")
cat("  10 mg/kg  :", tail(df_tgi$TGI_10, 1), "%\n")

# ── Simulation PK (doses répétées q14j) ──────────────────────────────────────
doses_h  <- c(0, 14, 28, 42) * 24
times_pk <- seq(0, 49*24, by=1)

pk3  <- sim_pk(3000,  times_pk, doses_h)
pk10 <- sim_pk(10000, times_pk, doses_h)

df_pk <- rbind(
  data.frame(jour=pk3$time/24,  C1=pk3$C1,  Dose="3 mg/kg"),
  data.frame(jour=pk10$time/24, C1=pk10$C1, Dose="10 mg/kg")
)
df_pk$Dose <- factor(df_pk$Dose, c("3 mg/kg","10 mg/kg"))

# ── Graphiques ────────────────────────────────────────────────────────────────
dose_days <- c(0, 14, 28, 42)
cols <- c("3 mg/kg"="#F8766D", "10 mg/kg"="#00BFC4")

p_tgi <- ggplot() +
  geom_line(data=data.frame(
      jour=rep(df_tgi$jour,2),
      TGI =c(df_tgi$TGI_3, df_tgi$TGI_10),
      Dose=rep(c("3 mg/kg","10 mg/kg"), each=nrow(df_tgi))),
    aes(x=jour, y=TGI, color=Dose), linewidth=1) +
  geom_point(data=data.frame(
      jour=rep(df_tgi$jour,2),
      TGI =c(df_tgi$TGI_3, df_tgi$TGI_10),
      Dose=rep(c("3 mg/kg","10 mg/kg"), each=nrow(df_tgi))),
    aes(x=jour, y=TGI, color=Dose), size=2.5) +
  geom_vline(xintercept=dose_days, linetype="dashed",
             color="grey60", linewidth=0.4) +
  geom_hline(yintercept=c(0,100), linetype="dotted", color="grey40") +
  scale_color_manual(values=cols) +
  labs(title="TGI observé — Fc-silent FGFR2-huBPA-LP1",
       x=NULL, y="TGI (%)") +
  theme_bw(base_size=12) +
  theme(legend.position="none")

p_pk <- ggplot(df_pk, aes(x=jour, y=C1, color=Dose)) +
  geom_line(linewidth=1) +
  geom_vline(xintercept=dose_days, linetype="dashed",
             color="grey60", linewidth=0.4) +
  scale_color_manual(values=cols) +
  scale_y_log10() +
  labs(x="Temps (jours)", y="Concentration (µg/L, log)",
       color="Dose") +
  theme_bw(base_size=12)

p_tgi / p_pk + plot_layout(heights=c(2,1))

ggsave("scripts/plot_TGI_simple_FGFR2.png", width=9, height=7, dpi=150)
cat("\nGraphique → scripts/plot_TGI_simple_FGFR2.png\n")
