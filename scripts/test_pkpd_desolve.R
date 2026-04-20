# Test PKPD Simeoni — deSolve + données synthétiques
# k2 calibré : k2*C1_max ~ 2×l0 → effet visible sans éradication complète

library(deSolve)
library(ggplot2)

# ── Paramètres PK (estimation 2-comp) ─────────────────────────────────────────
CL  <- 7.02e-4;  V1 <- 0.0711;  V2 <- 0.1326
k10 <- CL / V1                       # 0.00987 /h
k12 <- 0.0089                        # Q/V1 (approx)
k21 <- 0.0048                        # Q/V2 (approx)

# ── Paramètres PD fixés ───────────────────────────────────────────────────────
l0 <- 0.146 / 24   # 0.00608 /h
l1 <- 0.334 / 24   # 0.01392 g/h
p  <- 20

# ── Vrais paramètres à récupérer ──────────────────────────────────────────────
# k2*C1_max(10mg) = k2 * 10000/0.0711 ≈ 140 647 * k2
# Pour effet 2×l0 à pic : k2 ≈ 2*l0/140647 ≈ 8.6e-8
k1_true <- 0.10      # /h  (transit)
k2_true <- 1e-8      # k2*C1max ≈ 0.14*l0 → inhibition partielle
w0      <- 0.15      # g

dose10 <- 10000;  dose3 <- 3000;  dose0 <- 0   # µg/kg

# ── Modèle ODE ────────────────────────────────────────────────────────────────
pkpd_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    dA1 <- -(k10 + k12)*A1 + k21*A2
    dA2 <-   k12*A1 - k21*A2
    C1  <- A1 / V1

    w      <- x1 + x2 + x3 + x4
    growth <- l0*x1 / (1 + (l0/l1*w)^p)^(1/p)

    dx1 <- growth - k2*C1*x1
    dx2 <- k2*C1*x1 - k1*x2
    dx3 <- k1*(x2 - x3)
    dx4 <- k1*(x3 - x4)

    list(c(dA1, dA2, dx1, dx2, dx3, dx4))
  })
}

simulate_group <- function(dose, pars, times) {
  state0 <- c(A1=dose, A2=0, x1=w0, x2=0, x3=0, x4=0)  # A1=quantité
  out <- tryCatch(
    as.data.frame(lsoda(state0, times, pkpd_ode, pars,
                        rtol=1e-6, atol=1e-8)),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(times)))
  pmax(out$x1 + out$x2 + out$x3 + out$x4, 1e-6)
}

# ── Données synthétiques (bruit log-normal 10%) ───────────────────────────────
obs_days <- c(0, 4, 7, 11, 14, 18, 21, 25, 28, 32, 35, 39, 42, 46, 49)
obs_h    <- obs_days * 24

pars_true <- c(k10=k10, k12=k12, k21=k21, V1=V1,
               l0=l0, l1=l1, p=p, k1=k1_true, k2=k2_true)

set.seed(42)
noise <- function(x) pmax(x * exp(rnorm(length(x), 0, 0.10)), 1e-6)

w_ctrl <- noise(simulate_group(dose0,  pars_true, obs_h))
w_d3   <- noise(simulate_group(dose3,  pars_true, obs_h))
w_d10  <- noise(simulate_group(dose10, pars_true, obs_h))

cat("Vérification données synthétiques :\n")
cat("  Contrôle j0=", round(w_ctrl[1],3), " j49=", round(tail(w_ctrl,1),3), "\n")
cat("  3  mg/kg j0=", round(w_d3[1],3),   " j49=", round(tail(w_d3,1),3),   "\n")
cat("  10 mg/kg j0=", round(w_d10[1],3),  " j49=", round(tail(w_d10,1),3),  "\n")

dat_ctrl <- data.frame(t=obs_h, w=w_ctrl)
dat_d3   <- data.frame(t=obs_h, w=w_d3)
dat_d10  <- data.frame(t=obs_h, w=w_d10)

# ── Fonction objective ─────────────────────────────────────────────────────────
objective <- function(logpar) {
  k1  <- unname(exp(logpar[1]));  k2 <- unname(exp(logpar[2]))
  par <- c(k10=k10, k12=k12, k21=k21, V1=V1,
           l0=l0, l1=l1, p=p, k1=k1, k2=k2)

  pc  <- simulate_group(dose0,  par, dat_ctrl$t)
  p3  <- simulate_group(dose3,  par, dat_d3$t)
  p10 <- simulate_group(dose10, par, dat_d10$t)

  bad <- function(x) is.null(x) || any(is.na(x))
  if (bad(pc) || bad(p3) || bad(p10)) return(1e10)

  mean((log(dat_ctrl$w) - log(pc))^2) +
  mean((log(dat_d3$w)   - log(p3))^2) +
  mean((log(dat_d10$w)  - log(p10))^2)
}

# ── Multi-start nlminb ────────────────────────────────────────────────────────
cat("\nOptimisation...\n")
starts <- list(c(k1=0.1, k2=5e-8), c(k1=0.2, k2=1e-7),
               c(k1=0.05, k2=2e-8), c(k1=0.3, k2=1e-7))

best_obj <- Inf;  best_fit <- NULL
for (s in starts) {
  fit_try <- tryCatch(
    nlminb(log(s), objective,
           control=list(eval.max=2000, iter.max=1000,
                        rel.tol=1e-10, x.tol=1e-10)),
    error = function(e) NULL
  )
  if (!is.null(fit_try) && fit_try$objective < best_obj) {
    best_obj <- fit_try$objective;  best_fit <- fit_try
  }
}

k1_est <- unname(exp(best_fit$par[1]));  k2_est <- unname(exp(best_fit$par[2]))
cat("Vrais   : k1 =", k1_true,          "  k2 =", k2_true,          "\n")
cat("Estimés : k1 =", round(k1_est, 5), "  k2 =", round(k2_est, 9), "\n")
cat("Objectif:", round(best_obj, 6), "\n")

# ── Graphique ─────────────────────────────────────────────────────────────────
times_full <- seq(0, 49*24, by=2)
pars_est   <- c(k10=k10, k12=k12, k21=k21, V1=V1,
                l0=l0, l1=l1, p=p, k1=k1_est, k2=k2_est)

lev <- c("Contrôle","3 mg/kg","10 mg/kg")
df_sim <- rbind(
  data.frame(t=times_full/24, w=simulate_group(dose0,  pars_est, times_full), Groupe="Contrôle"),
  data.frame(t=times_full/24, w=simulate_group(dose3,  pars_est, times_full), Groupe="3 mg/kg"),
  data.frame(t=times_full/24, w=simulate_group(dose10, pars_est, times_full), Groupe="10 mg/kg")
)
df_obs <- rbind(
  data.frame(t=obs_days, w=w_ctrl, Groupe="Contrôle"),
  data.frame(t=obs_days, w=w_d3,   Groupe="3 mg/kg"),
  data.frame(t=obs_days, w=w_d10,  Groupe="10 mg/kg")
)
df_sim$Groupe <- factor(df_sim$Groupe, lev)
df_obs$Groupe <- factor(df_obs$Groupe, lev)

ggplot() +
  geom_line(data=df_sim,  aes(x=t, y=w, color=Groupe), linewidth=1) +
  geom_point(data=df_obs, aes(x=t, y=w, color=Groupe), size=2.5) +
  labs(title="Test PKPD Simeoni (deSolve) — données synthétiques",
       subtitle=paste0("k1=", round(k1_est,4), " /h   k2=", round(k2_est,9),
                       "   (vrais: k1=", k1_true, " k2=", k2_true, ")"),
       x="Temps (jours)", y="Masse tumorale (g)") +
  theme_bw(base_size=13)

ggsave("scripts/test_pkpd_desolve.png", width=9, height=5, dpi=150)
cat("Test OK → scripts/test_pkpd_desolve.png\n")
