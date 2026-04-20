# TGI — Fc-silent FGFR2-huBPA-LP1
# PK 2-comp fixé + modèle de croissance Simeoni (deSolve)
# TGI(%) = (1 - ΔW_traité / ΔW_contrôle) × 100  à chaque point de mesure

library(deSolve)
library(ggplot2)
library(readxl)

# ── PK fixé ───────────────────────────────────────────────────────────────────
load("scripts/resultats_PK2comp_FGFR2.RData")   # → pk2comp

CL  <- unname(pk2comp$CL);  V1 <- unname(pk2comp$V1)
V2  <- unname(pk2comp$V2);  Q  <- unname(pk2comp$Q)
k10 <- CL/V1;  k12 <- Q/V1;  k21 <- Q/V2

# ── Paramètres PD (littérature, en heures) ────────────────────────────────────
l0 <- 0.146/24;  l1 <- 0.334/24;  p <- 20

# ── Données tumorales ─────────────────────────────────────────────────────────
raw <- read_xlsx("data/TumorVolume_FGFR2.xlsx",
                 n_max=3, col_types=c("text", rep("numeric", 15)))

time_d <- suppressWarnings(as.numeric(colnames(raw)[-1]))
time_h <- time_d * 24

get_w <- function(pat)
  as.numeric(unlist(raw[grep(pat, raw[[1]], ignore.case=TRUE)[1], -1]))

# mm³ → g
w_ctrl <- get_w("Group 01") / 1000
w_d10  <- get_w("Group 03") / 1000
w_d3   <- get_w("Group 04") / 1000

ok_all <- !is.na(w_ctrl) & !is.na(w_d10) & !is.na(w_d3)
time_d <- time_d[ok_all];  time_h <- time_h[ok_all]
w_ctrl <- w_ctrl[ok_all];  w_d10 <- w_d10[ok_all];  w_d3 <- w_d3[ok_all]

w0 <- mean(c(w_ctrl[time_d==0], w_d10[time_d==0], w_d3[time_d==0]))
dose10 <- 10000;  dose3 <- 3000;  dose0 <- 0

# ── Modèle ODE ────────────────────────────────────────────────────────────────
pkpd_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    dA1 <- -(k10+k12)*A1 + k21*A2
    dA2 <-   k12*A1 - k21*A2
    C1  <- A1/V1
    w   <- x1+x2+x3+x4
    g   <- l0*x1/(1+(l0/l1*w)^p)^(1/p)
    list(c(dA1, dA2, g-k2*C1*x1, k2*C1*x1-k1*x2,
           k1*(x2-x3), k1*(x3-x4)))
  })
}

pars_base <- c(k10=k10, k12=k12, k21=k21, V1=V1, l0=l0, l1=l1, p=p)

simulate <- function(dose, k1, k2, times, doses_h=0) {
  pars   <- c(pars_base, k1=k1, k2=k2)
  state0 <- c(A1=dose, A2=0, x1=w0, x2=0, x3=0, x4=0)
  all_t  <- sort(unique(c(0, times, doses_h)))

  events <- if (length(doses_h) > 1)
    data.frame(var="A1", time=doses_h[-1], value=dose, method="add")
  else NULL

  out <- tryCatch(
    as.data.frame(lsoda(state0, all_t, pkpd_ode, pars,
                        events=if (!is.null(events)) list(data=events) else NULL,
                        rtol=1e-6, atol=1e-8)),
    error=function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(times)))
  w_tot <- pmax(out$x1+out$x2+out$x3+out$x4, 1e-9)
  approx(out$time, w_tot, xout=times, rule=2)$y
}

# ── Optimisation k1, k2 (doses répétées) ─────────────────────────────────────
doses_h <- c(0, 14, 28, 42) * 24
obs_h   <- time_h[time_d > 0]
obs_idx <- time_d > 0

objective <- function(logpar) {
  k1 <- unname(exp(logpar[1]));  k2 <- unname(exp(logpar[2]))
  pc  <- simulate(dose0,  k1, k2, obs_h)
  p3  <- simulate(dose3,  k1, k2, obs_h, doses_h)
  p10 <- simulate(dose10, k1, k2, obs_h, doses_h)
  bad <- function(x) is.null(x) || any(is.na(x) | x <= 0)
  if (bad(pc)||bad(p3)||bad(p10)) return(1e10)
  val <- mean((log(w_ctrl[obs_idx])-log(pc))^2) +
         mean((log(w_d3[obs_idx])  -log(p3))^2) +
         mean((log(w_d10[obs_idx]) -log(p10))^2)
  if (!is.finite(val)) 1e10 else val
}

cat("Optimisation k1/k2...\n")
starts <- list(c(k1=0.5,k2=1e-4), c(k1=0.1,k2=1e-5),
               c(k1=1.0,k2=1e-4), c(k1=0.2,k2=5e-5))
best <- Inf;  fit <- NULL
for (s in starts) {
  f <- tryCatch(nlminb(log(s), objective,
                       control=list(eval.max=3000, iter.max=1500,
                                    rel.tol=1e-12, x.tol=1e-12)),
                error=function(e) NULL)
  if (!is.null(f) && is.finite(f$objective) && f$objective < best) {
    best <- f$objective;  fit <- f
  }
}
k1 <- unname(exp(fit$par[1]));  k2 <- unname(exp(fit$par[2]))
cat("  k1 =", round(k1,5), "/h   k2 =", round(k2,8), "\n")

# ── TGI à chaque point de mesure ─────────────────────────────────────────────
all_h <- time_h
wc  <- simulate(dose0,  k1, k2, all_h)
w3  <- simulate(dose3,  k1, k2, all_h, doses_h)
w10 <- simulate(dose10, k1, k2, all_h, doses_h)

tgi <- function(wt, wc) round((1 - (wt - w0) / (wc - w0)) * 100, 1)

df_tgi <- data.frame(
  jour   = time_d[time_d > 0],
  TGI_3  = tgi(w3[time_d > 0],  wc[time_d > 0]),
  TGI_10 = tgi(w10[time_d > 0], wc[time_d > 0])
)

cat("\n=== TGI simulé (doses répétées toutes les 14j) ===\n")
cat(sprintf("%-6s  %10s  %10s\n", "Jour", "3 mg/kg", "10 mg/kg"))
cat(strrep("-", 30), "\n")
for (i in seq_len(nrow(df_tgi)))
  cat(sprintf("j%-5d  %9.1f%%  %9.1f%%\n",
              df_tgi$jour[i], df_tgi$TGI_3[i], df_tgi$TGI_10[i]))

cat("\n--- Dernier point (j", tail(df_tgi$jour,1), ") ---\n", sep="")
cat("   3 mg/kg  :", tail(df_tgi$TGI_3,  1), "%\n")
cat("  10 mg/kg  :", tail(df_tgi$TGI_10, 1), "%\n")

# ── Graphique TGI ─────────────────────────────────────────────────────────────
df_long <- rbind(
  data.frame(jour=df_tgi$jour, TGI=df_tgi$TGI_3,  Dose="3 mg/kg"),
  data.frame(jour=df_tgi$jour, TGI=df_tgi$TGI_10, Dose="10 mg/kg")
)
df_long$Dose <- factor(df_long$Dose, c("3 mg/kg","10 mg/kg"))

ggplot(df_long, aes(x=jour, y=TGI, color=Dose)) +
  geom_line(linewidth=1) +
  geom_point(size=2.5) +
  geom_vline(xintercept=c(0,14,28,42), linetype="dashed",
             color="grey60", linewidth=0.4) +
  geom_hline(yintercept=c(0,100), linetype="dotted", color="grey40") +
  labs(title   = "TGI simulé — Fc-silent FGFR2-huBPA-LP1",
       subtitle = paste0("Doses répétées j0/j14/j28/j42 | ",
                         "k1=", round(k1,4), " /h  k2=", round(k2,8)),
       x = "Temps (jours)", y = "TGI (%)") +
  theme_bw(base_size=13)

ggsave("scripts/plot_TGI_simple_FGFR2.png", width=8, height=5, dpi=150)
cat("\nGraphique → scripts/plot_TGI_simple_FGFR2.png\n")
