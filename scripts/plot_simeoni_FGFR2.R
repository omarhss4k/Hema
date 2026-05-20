# =============================================================================
# Plot PK/PD Simeoni — charge les résultats sauvegardés et génère le graphe
# =============================================================================

library(deSolve)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. CHARGEMENT DES RÉSULTATS
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2
load("scripts/resultats_PKPD_simeoni_FGFR2.RData")      # → simeoni_results

pk_fixed <- c(
  CL = unname(pk2comp_rxode2$CL) * 24,
  V1 = unname(pk2comp_rxode2$V1),
  V2 = unname(pk2comp_rxode2$V2),
  Q  = unname(pk2comp_rxode2$Q)  * 24
)

L0  <- simeoni_results$L0
L1  <- simeoni_results$L1
k1  <- simeoni_results$k1
k2  <- simeoni_results$k2
TSC <- simeoni_results$TSC_ugL
MTT <- simeoni_results$MTT_days

# =============================================================================
# 2. DONNÉES TUMORALES
# =============================================================================

raw_all <- suppressMessages(
  read_xlsx("TumorVolume_FGFR2.xlsx", col_names = FALSE, .name_repair = "minimal")
)
raw_all    <- raw_all[, colSums(!is.na(raw_all)) > 0]
header_idx <- which(raw_all[[1]] == "Group")

times_all <- suppressWarnings(
  as.numeric(as.character(unlist(raw_all[header_idx[1], -1])))
)
valid   <- !is.na(times_all)
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

tv_ctrl  <- get_row(blk_mean, "Group 01")
tv_iso   <- get_row(blk_mean, "Group 02")
tv_d3    <- get_row(blk_mean, "Group 04")
tv_d10   <- get_row(blk_mean, "Group 03")

sem_ctrl <- pmax(get_row(blk_sem, "Group 01"), 1)
sem_iso  <- pmax(get_row(blk_sem, "Group 02"), 1)
sem_d3   <- pmax(get_row(blk_sem, "Group 04"), 1)
sem_d10  <- pmax(get_row(blk_sem, "Group 03"), 1)

get_tv0 <- function(tv) {
  v <- tv[times_d == 0]
  if (length(v) == 0 || is.na(v[1])) tv[!is.na(tv)][1] else v[1]
}
tv0_ctrl <- get_tv0(tv_ctrl)
tv0_iso  <- get_tv0(tv_iso)
tv0_d3   <- get_tv0(tv_d3)
tv0_d10  <- get_tv0(tv_d10)

ok_ctrl <- !is.na(tv_ctrl)
ok_iso  <- !is.na(tv_iso)
ok_d3   <- !is.na(tv_d3)
ok_d10  <- !is.na(tv_d10)

# =============================================================================
# 3. MODÈLE ODE
# =============================================================================

PSI <- 20

simeoni_rhs <- function(t, state, parms) {
  A1 <- max(state["A1"], 0); A2 <- max(state["A2"], 0)
  x1 <- max(state["x1"], 0); x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0); x4 <- max(state["x4"], 0)

  CL <- parms["CL"]; V1 <- parms["V1"]
  V2 <- parms["V2"]; Q  <- parms["Q"]
  L0 <- parms["L0"]; L1 <- parms["L1"]
  k1 <- parms["k1"]; k2 <- parms["k2"]

  C1  <- A1 / V1
  dA1 <- -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
  dA2 <-  (Q/V1)*A1 - (Q/V2)*A2

  w   <- x1 + x2 + x3 + x4
  gw  <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else L0

  list(c(A1 = dA1, A2 = dA2,
         x1 = (growth_rate - k2*C1)*x1,
         x2 = k2*C1*x1 - k1*x2,
         x3 = k1*x2 - k1*x3,
         x4 = k1*x3 - k1*x4))
}

simeoni_ctrl_rhs <- function(t, state, parms) {
  x1 <- max(state["x1"], 0)
  w  <- x1 + max(state["x2"],0) + max(state["x3"],0) + max(state["x4"],0)
  gw <- parms["L0"] * w / (1 + (parms["L0"]*w/parms["L1"])^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else parms["L0"]
  list(c(x1 = growth_rate * x1, x2 = 0, x3 = 0, x4 = 0))
}

# =============================================================================
# 4. FONCTIONS DE SIMULATION
# =============================================================================

dose_days <- c(0, 14, 28, 42)

sim_ctrl_fn <- function(tv0, times_out) {
  t_all <- sort(unique(c(0, times_out)))
  out <- tryCatch(
    as.data.frame(lsoda(
      y     = c(x1 = tv0, x2 = 0, x3 = 0, x4 = 0),
      times = t_all,
      func  = simeoni_ctrl_rhs,
      parms = c(L0 = L0, L1 = L1)
    )),
    error = function(e) NULL
  )
  if (is.null(out)) return(rep(NA_real_, length(times_out)))
  w <- pmax(out$x1,0) + pmax(out$x2,0) + pmax(out$x3,0) + pmax(out$x4,0)
  approx(out$time, w, xout = times_out, rule = 2)$y
}

sim_treated <- function(dose_ugkg, tv0, times_out) {
  params <- c(pk_fixed, L0 = L0, L1 = L1, k1 = k1, k2 = k2)
  breaks <- c(dose_days, max(times_out) + 1)
  state  <- c(A1 = dose_ugkg, A2 = 0, x1 = tv0, x2 = 0, x3 = 0, x4 = 0)
  t_all  <- numeric(0); w_all <- numeric(0)

  for (i in seq_along(dose_days)) {
    t_start <- dose_days[i]
    t_stop  <- min(breaks[i + 1], max(times_out))
    if (t_start >= max(times_out)) break
    if (i > 1) state["A1"] <- state["A1"] + dose_ugkg

    t_seg <- sort(unique(c(t_start,
                           times_out[times_out > t_start & times_out <= t_stop],
                           t_stop)))
    if (length(t_seg) < 2) t_seg <- c(t_start, t_stop)

    seg <- tryCatch(
      as.data.frame(lsoda(state, t_seg, simeoni_rhs, params, atol = 1e-6, rtol = 1e-6)),
      error = function(e) NULL
    )
    if (is.null(seg)) return(rep(NA_real_, length(times_out)))

    w_seg <- pmax(seg$x1,0) + pmax(seg$x2,0) + pmax(seg$x3,0) + pmax(seg$x4,0)
    keep  <- if (length(t_all) > 0) seg$time > tail(t_all, 1) else rep(TRUE, nrow(seg))
    t_all <- c(t_all, seg$time[keep]); w_all <- c(w_all, w_seg[keep])

    last  <- seg[nrow(seg), ]
    state <- c(A1=last$A1, A2=last$A2, x1=last$x1, x2=last$x2, x3=last$x3, x4=last$x4)
  }
  if (length(t_all) == 0) return(rep(NA_real_, length(times_out)))
  approx(t_all, w_all, xout = times_out, rule = 2)$y
}

# =============================================================================
# 5. SIMULATION
# =============================================================================

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)

# λ0 re-estimé sur les données contrôle seules (régression log-linéaire j0-j25)
ok_early <- ok_ctrl & times_d <= 25
L0_ctrl_obs <- max(coef(lm(log(tv_ctrl[ok_early]) ~ times_d[ok_early]))[2], 0.005)
cat(sprintf("λ0 contrôle (régression obs) = %.4f /j\n", L0_ctrl_obs))

pred_ctrl <- sim_ctrl_fn(tv0_ctrl, times_sim)
# Remplace L0 par la valeur observée pour la prédiction contrôle
environment(sim_ctrl_fn)  # juste pour rappel — on redéfinit localement
pred_ctrl_obs <- {
  L0_save <- L0; L0 <<- L0_ctrl_obs
  res <- sim_ctrl_fn(tv0_ctrl, times_sim)
  L0 <<- L0_save
  res
}
pred_d3   <- sim_treated(3000,  tv0_d3,  times_sim)
pred_d10  <- sim_treated(10000, tv0_d10, times_sim)

niv <- c("Contrôle", "Isotype 10mg/kg", "3 mg/kg", "10 mg/kg")

# Prédiction uniquement pour véhicule et groupes traités (pas pour isotype)
df_sim <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_obs, Groupe = "Contrôle"),
  data.frame(jour = times_sim, TV = pred_d3,   Groupe = "3 mg/kg"),
  data.frame(jour = times_sim, TV = pred_d10,  Groupe = "10 mg/kg")
)

# Points observés pour tous les groupes y compris isotype
df_obs <- rbind(
  data.frame(jour = times_d[ok_ctrl], TV = tv_ctrl[ok_ctrl], sem = sem_ctrl[ok_ctrl], Groupe = "Contrôle"),
  data.frame(jour = times_d[ok_iso],  TV = tv_iso[ok_iso],   sem = sem_iso[ok_iso],   Groupe = "Isotype 10mg/kg"),
  data.frame(jour = times_d[ok_d3],   TV = tv_d3[ok_d3],     sem = sem_d3[ok_d3],     Groupe = "3 mg/kg"),
  data.frame(jour = times_d[ok_d10],  TV = tv_d10[ok_d10],   sem = sem_d10[ok_d10],   Groupe = "10 mg/kg")
)

df_sim$Groupe <- factor(df_sim$Groupe, levels = niv)
df_obs$Groupe <- factor(df_obs$Groupe, levels = niv)

# =============================================================================
# 6. GRAPHIQUE
# =============================================================================

cols <- c("Contrôle"        = "#888888",
          "Isotype 10mg/kg" = "#CC4444",
          "3 mg/kg"         = "#4393C3",
          "10 mg/kg"        = "#2166AC")

ymax <- max(df_obs$TV + df_obs$sem, na.rm = TRUE)

# Plafonne la prédiction contrôle pour garder l'échelle après j25
df_sim$TV[df_sim$Groupe == "Contrôle" & !is.na(df_sim$TV) & df_sim$TV > ymax * 1.05] <- NA

subtitle_txt <- paste0(
  "λ0 = ", round(L0, 4), " /j  |  ",
  "k2 = ", formatC(k2, digits = 3, format = "e"), " L/µg/j  |  ",
  "k1 = ", round(k1, 3), " /j  |  ",
  "TSC = ", round(TSC, 0), " µg/L  |  ",
  "MTT = ", round(MTT, 1), " j"
)

p <- ggplot() +
  geom_vline(xintercept = dose_days, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_line(data = df_sim[df_sim$Groupe != "Contrôle", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1) +
  geom_line(data = df_sim[df_sim$Groupe == "Contrôle", ],
            aes(x = jour, y = TV, color = Groupe, group = Groupe),
            linewidth = 1, linetype = "dashed") +
  geom_errorbar(data = df_obs,
                aes(x = jour, ymin = TV - sem, ymax = TV + sem, color = Groupe),
                width = 0.8, linewidth = 0.5) +
  geom_point(data = df_obs,
             aes(x = jour, y = TV, color = Groupe),
             size = 2.5) +
  annotate("point", x = dose_days, y = -ymax * 0.06,
           shape = 17, size = 3.5, color = "#CC0000") +
  annotate("text",  x = max(dose_days) + 1.5, y = -ymax * 0.06,
           label = "= Traitement", hjust = 0, size = 3.2, color = "#CC0000") +
  scale_x_continuous(breaks = seq(0, max(times_d, na.rm = TRUE), by = 7)) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(-ymax * 0.12, ymax * 1.1), clip = "off") +
  labs(
    title    = "Modèle PK/PD TGI — Simeoni (2004) — Fc-silent FGFR2-huBPA-LP1",
    subtitle = subtitle_txt,
    x        = "Temps (jours)",
    y        = "Volume tumoral (mm³)",
    color    = NULL,
    caption  = "Lignes = modèle   ● = données observées (mean ± SEM)   -- = Contrôle (prédiction)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.subtitle   = element_text(size = 9, color = "grey50"),
    plot.margin     = margin(t = 5, r = 10, b = 25, l = 5)
  )

print(p)
ggsave("scripts/plot_PKPD_simeoni_FGFR2.png", p,
       width = 9, height = 5.5, dpi = 150)
cat("Graphique → scripts/plot_PKPD_simeoni_FGFR2.png\n")
