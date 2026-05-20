# =============================================================================
# diag_isotype_vs_ctrl.R
# Diagnostic : isotype (Group 02) vs contrôle véhicule (Group 01)
#
# Objectif : vérifier si l'exclusion de l'isotype de la fonction objective
# est justifiée, en quantifiant l'écart entre les deux groupes.
#
# Dépendances : doit être exécuté APRÈS pkpd_TGI_simeoni_FGFR2.R
#   (utilise resultats_PKPD_simeoni_FGFR2.RData et resultats_PK2comp_rxode2_FGFR2.RData)
#
# Sorties :
#   - Console  : taux de croissance, résidus, tests statistiques
#   - Figure   : scripts/diag_isotype_vs_ctrl.png
# =============================================================================

library(deSolve)
library(ggplot2)
library(readxl)

# =============================================================================
# 1. CHARGEMENT DES PARAMÈTRES ESTIMÉS
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")  # → pk2comp_rxode2
load("scripts/resultats_PKPD_simeoni_FGFR2.RData")    # → simeoni_results

L0 <- simeoni_results$L0
L1 <- simeoni_results$L1

# =============================================================================
# 2. DONNÉES TUMORALES (même extraction que le script principal)
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
sem_ctrl <- pmax(get_row(blk_sem, "Group 01"), 1)
sem_iso  <- pmax(get_row(blk_sem, "Group 02"), 1)

ok_ctrl <- !is.na(tv_ctrl)
ok_iso  <- !is.na(tv_iso)

get_tv0 <- function(tv) {
  v <- tv[times_d == 0]
  if (length(v) == 0 || is.na(v[1])) tv[!is.na(tv)][1] else v[1]
}
tv0_ctrl <- get_tv0(tv_ctrl)
tv0_iso  <- get_tv0(tv_iso)

# Timepoints communs aux deux groupes
common_times <- times_d[ok_ctrl & ok_iso]
tv_ctrl_c    <- tv_ctrl[ok_ctrl & ok_iso]
tv_iso_c     <- tv_iso[ok_ctrl & ok_iso]
sem_ctrl_c   <- sem_ctrl[ok_ctrl & ok_iso]
sem_iso_c    <- sem_iso[ok_ctrl & ok_iso]

# =============================================================================
# 3. TAUX DE CROISSANCE — RÉGRESSION LOG-LINÉAIRE
#    Estimé sur la phase précoce (j0–j25) où la croissance est exponentielle
# =============================================================================

early <- times_d <= 25

lm_ctrl_all  <- lm(log(tv_ctrl[ok_ctrl])        ~ times_d[ok_ctrl])
lm_iso_all   <- lm(log(tv_iso[ok_iso])           ~ times_d[ok_iso])
lm_ctrl_early <- lm(log(tv_ctrl[ok_ctrl & early]) ~ times_d[ok_ctrl & early])
lm_iso_early  <- lm(log(tv_iso[ok_iso   & early]) ~ times_d[ok_iso  & early])

L0_ctrl_all   <- coef(lm_ctrl_all)[2]
L0_iso_all    <- coef(lm_iso_all)[2]
L0_ctrl_early <- coef(lm_ctrl_early)[2]
L0_iso_early  <- coef(lm_iso_early)[2]

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DIAGNOSTIC ISOTYPE vs CONTRÔLE                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("── 1. TAUX DE CROISSANCE (régression log-linéaire) ──────────\n")
cat(sprintf("  Contrôle (tous timepoints)  : λ0 = %.4f /j  (t½ = %.1f j)\n",
            L0_ctrl_all, log(2)/L0_ctrl_all))
cat(sprintf("  Isotype  (tous timepoints)  : λ0 = %.4f /j  (t½ = %.1f j)\n",
            L0_iso_all,  log(2)/L0_iso_all))
cat(sprintf("  Contrôle (j0–j25)           : λ0 = %.4f /j\n", L0_ctrl_early))
cat(sprintf("  Isotype  (j0–j25)           : λ0 = %.4f /j\n", L0_iso_early))
cat(sprintf("  Δλ0 (iso − ctrl, j0–j25)   : %.4f /j  (%.1f %%)\n\n",
            L0_iso_early - L0_ctrl_early,
            100 * (L0_iso_early - L0_ctrl_early) / L0_ctrl_early))

# =============================================================================
# 4. PRÉDICTION CONTRÔLE AUX TEMPS ISOTYPE
#    Modèle Simeoni sans drogue, avec tv0_iso comme CI
# =============================================================================

PSI <- 20

simeoni_ctrl_rhs <- function(t, state, parms) {
  x1 <- max(state["x1"], 0)
  x2 <- max(state["x2"], 0)
  x3 <- max(state["x3"], 0)
  x4 <- max(state["x4"], 0)
  w  <- x1 + x2 + x3 + x4
  gw <- parms["L0"] * w / (1 + (parms["L0"] * w / parms["L1"])^PSI)^(1/PSI)
  growth_rate <- if (w > 1e-12) gw / w else parms["L0"]
  list(c(x1 = growth_rate * x1, x2 = 0, x3 = 0, x4 = 0))
}

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
  w <- pmax(out$x1, 0) + pmax(out$x2, 0) + pmax(out$x3, 0) + pmax(out$x4, 0)
  approx(out$time, w, xout = times_out, rule = 2)$y
}

# Prédiction contrôle aux timepoints communs (CI = tv0_ctrl ou tv0_iso)
pred_ctrl_at_ctrl <- sim_ctrl_fn(tv0_ctrl, common_times)
pred_ctrl_at_iso  <- sim_ctrl_fn(tv0_iso,  common_times)   # même lambda, CI isotype

# =============================================================================
# 5. RÉSIDUS : isotype vs prédiction contrôle (log-espace)
# =============================================================================

resid_log_iso_vs_ctrl_pred <- log(tv_iso_c) - log(pred_ctrl_at_iso)
resid_log_ctrl_vs_ctrl_pred <- log(tv_ctrl_c) - log(pred_ctrl_at_ctrl)

# Ratio TV isotype / TV contrôle observé
ratio_iso_ctrl <- tv_iso_c / tv_ctrl_c

cat("── 2. RÉSIDUS LOG : isotype vs prédiction contrôle ─────────\n")
cat(sprintf("  Biais moyen (log)  : %.4f  (0 = pas d'écart)\n",
            mean(resid_log_iso_vs_ctrl_pred, na.rm = TRUE)))
cat(sprintf("  RMSE (log)         : %.4f\n",
            sqrt(mean(resid_log_iso_vs_ctrl_pred^2, na.rm = TRUE))))
cat(sprintf("  Max |résidu| (log) : %.4f  (tj = %.0f)\n",
            max(abs(resid_log_iso_vs_ctrl_pred), na.rm = TRUE),
            common_times[which.max(abs(resid_log_iso_vs_ctrl_pred))]))
cat(sprintf("\n  Pour comparaison — résidus contrôle vs pred contrôle :\n"))
cat(sprintf("  Biais moyen (log)  : %.4f\n",
            mean(resid_log_ctrl_vs_ctrl_pred, na.rm = TRUE)))
cat(sprintf("  RMSE (log)         : %.4f\n\n",
            sqrt(mean(resid_log_ctrl_vs_ctrl_pred^2, na.rm = TRUE))))

cat("── 3. RATIO TV_isotype / TV_contrôle (par timepoint) ───────\n")
for (i in seq_along(common_times)) {
  cat(sprintf("  j%-3.0f : iso=%.0f  ctrl=%.0f  ratio=%.2f  %s\n",
              common_times[i], tv_iso_c[i], tv_ctrl_c[i], ratio_iso_ctrl[i],
              if (ratio_iso_ctrl[i] > 1.3 | ratio_iso_ctrl[i] < 0.77)
                "<-- divergence >30%" else ""))
}

# =============================================================================
# 6. TEST STATISTIQUE : les deux courbes log-linéaires sont-elles parallèles ?
#    Modèle avec interaction Groupe × Temps sur la phase exponentielle
# =============================================================================

df_lm <- data.frame(
  logTV = c(log(tv_ctrl[ok_ctrl & early]), log(tv_iso[ok_iso & early])),
  temps = c(times_d[ok_ctrl & early],      times_d[ok_iso & early]),
  groupe = c(rep("ctrl", sum(ok_ctrl & early)), rep("iso", sum(ok_iso & early)))
)

fit_additive    <- lm(logTV ~ temps + groupe,          data = df_lm)
fit_interaction <- lm(logTV ~ temps * groupe,          data = df_lm)
anova_res       <- anova(fit_additive, fit_interaction)

cat("\n── 4. TEST D'INTERACTION Groupe × Temps (j0–j25) ───────────\n")
cat(sprintf("  Modèle additif     (pentes parallèles)  : AIC = %.2f\n",
            AIC(fit_additive)))
cat(sprintf("  Modèle interaction (pentes différentes) : AIC = %.2f\n",
            AIC(fit_interaction)))
cat(sprintf("  F(1, df) = %.3f   p = %.4f\n",
            anova_res$F[2], anova_res$`Pr(>F)`[2]))

p_val <- anova_res$`Pr(>F)`[2]
conclusion <- if (is.na(p_val) || p_val > 0.05) {
  "NON SIGNIFICATIF (p > 0.05) → pentes parallèles → exclusion isotype JUSTIFIÉE"
} else if (p_val > 0.01) {
  "MARGINAL (0.01 < p ≤ 0.05) → légère divergence → à interpréter avec prudence"
} else {
  "SIGNIFICATIF (p ≤ 0.01) → taux de croissance différents → exclusion DISCUTABLE"
}
cat(sprintf("\n  Conclusion : %s\n\n", conclusion))

# =============================================================================
# 7. VERDICT GLOBAL
# =============================================================================

biais   <- abs(mean(resid_log_iso_vs_ctrl_pred, na.rm = TRUE))
rmse    <- sqrt(mean(resid_log_iso_vs_ctrl_pred^2, na.rm = TRUE))
max_ratio_dev <- max(abs(ratio_iso_ctrl - 1), na.rm = TRUE)

cat("── 5. VERDICT GLOBAL ────────────────────────────────────────\n")
cat(sprintf("  Biais log moyen    : %.3f  (seuil acceptable : < 0.15)\n", biais))
cat(sprintf("  RMSE log           : %.3f  (seuil acceptable : < 0.20)\n", rmse))
cat(sprintf("  Max écart ratio    : %.0f%%  (seuil acceptable : < 30%%)\n",
            100 * max_ratio_dev))

ok_biais <- biais < 0.15
ok_rmse  <- rmse  < 0.20
ok_ratio <- max_ratio_dev < 0.30
ok_pval  <- is.na(p_val) || p_val > 0.05

n_ok <- sum(c(ok_biais, ok_rmse, ok_ratio, ok_pval))
if (n_ok == 4) {
  verdict <- "EXCLUSION JUSTIFIÉE — l'isotype se comporte comme le contrôle véhicule"
} else if (n_ok >= 3) {
  verdict <- "EXCLUSION PROBABLEMENT JUSTIFIÉE — écart mineur, surveiller"
} else {
  verdict <- "EXCLUSION DISCUTABLE — l'isotype diverge du contrôle, envisager son inclusion"
}
cat(sprintf("\n  *** %s ***\n\n", verdict))

# =============================================================================
# 8. FIGURE DIAGNOSTIC
# =============================================================================

times_sim <- seq(0, max(times_d, na.rm = TRUE) * 1.05, by = 0.5)
pred_ctrl_sim <- sim_ctrl_fn(tv0_ctrl, times_sim)
pred_iso_sim  <- sim_ctrl_fn(tv0_iso,  times_sim)

# Panel A : trajectoires observées + prédictions
df_pred <- rbind(
  data.frame(jour = times_sim, TV = pred_ctrl_sim, Groupe = "Contrôle (pred)"),
  data.frame(jour = times_sim, TV = pred_iso_sim,  Groupe = "Isotype — pred ctrl sans drogue")
)
df_pts <- rbind(
  data.frame(jour = common_times, TV = tv_ctrl_c, sem = sem_ctrl_c, Groupe = "Contrôle (obs)"),
  data.frame(jour = common_times, TV = tv_iso_c,  sem = sem_iso_c,  Groupe = "Isotype (obs)")
)

cols_pred <- c("Contrôle (pred)"                 = "#888888",
               "Isotype — pred ctrl sans drogue"  = "#CC4444")
cols_pts  <- c("Contrôle (obs)" = "#888888",
               "Isotype (obs)"  = "#CC4444")

p_traj <- ggplot() +
  geom_line(data = df_pred,
            aes(x = jour, y = TV, color = Groupe, linetype = Groupe),
            linewidth = 0.9) +
  geom_errorbar(data = df_pts,
                aes(x = jour, ymin = TV - sem, ymax = TV + sem, color = Groupe),
                width = 0.8, linewidth = 0.5) +
  geom_point(data = df_pts,
             aes(x = jour, y = TV, color = Groupe),
             size = 2.5) +
  scale_color_manual(values = c(cols_pred, cols_pts)) +
  scale_linetype_manual(values = c("Contrôle (pred)" = "dashed",
                                   "Isotype — pred ctrl sans drogue" = "dashed")) +
  labs(title    = "Panel A — Trajectoires observées vs prédiction contrôle",
       subtitle = sprintf("λ0 estimé = %.4f /j  |  Biais log moyen = %.3f  |  RMSE log = %.3f",
                          L0, biais, rmse),
       x = "Temps (jours)", y = "Volume tumoral (mm³)", color = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 8.5, color = "grey50"))

# Panel B : résidus log isotype vs prédiction contrôle
df_resid <- data.frame(
  jour   = common_times,
  residu = resid_log_iso_vs_ctrl_pred,
  label  = sprintf("%.2f", resid_log_iso_vs_ctrl_pred)
)

p_resid <- ggplot(df_resid, aes(x = jour, y = residu)) +
  geom_hline(yintercept = 0,    color = "grey50", linetype = "dashed") +
  geom_hline(yintercept =  0.15, color = "#E69F00", linetype = "dotted", linewidth = 0.7) +
  geom_hline(yintercept = -0.15, color = "#E69F00", linetype = "dotted", linewidth = 0.7) +
  geom_col(aes(fill = residu > 0), width = 1.2, alpha = 0.75) +
  geom_text(aes(label = label,
                vjust = ifelse(residu >= 0, -0.4, 1.3)),
            size = 3) +
  scale_fill_manual(values = c(`TRUE` = "#4393C3", `FALSE` = "#D6604D"),
                    labels = c("Iso < pred", "Iso > pred"),
                    name = NULL) +
  scale_x_continuous(breaks = common_times) +
  labs(title    = "Panel B — Résidus log(iso_obs) − log(pred_ctrl)",
       subtitle = "Bandes oranges = seuil d'acceptabilité ±0.15",
       x = "Temps (jours)", y = "Résidu log") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 8.5, color = "grey50"))

# Panel C : ratio TV_isotype / TV_contrôle
df_ratio <- data.frame(
  jour  = common_times,
  ratio = ratio_iso_ctrl
)

p_ratio <- ggplot(df_ratio, aes(x = jour, y = ratio)) +
  geom_hline(yintercept = 1,    color = "grey50", linetype = "dashed") +
  geom_hline(yintercept = 1.30, color = "#E69F00", linetype = "dotted", linewidth = 0.7) +
  geom_hline(yintercept = 0.77, color = "#E69F00", linetype = "dotted", linewidth = 0.7) +
  geom_line(color = "#CC4444", linewidth = 0.9) +
  geom_point(color = "#CC4444", size = 2.5) +
  geom_text(aes(label = sprintf("%.2f", ratio),
                vjust = ifelse(ratio >= 1, -0.6, 1.4)),
            size = 3) +
  scale_x_continuous(breaks = common_times) +
  scale_y_continuous(limits = c(
    min(0.5, min(ratio_iso_ctrl, na.rm = TRUE) * 0.9),
    max(1.8, max(ratio_iso_ctrl, na.rm = TRUE) * 1.1)
  )) +
  labs(title    = "Panel C — Ratio TV_isotype / TV_contrôle",
       subtitle = "Bandes oranges = seuil ±30%  |  1.0 = comportement identique",
       x = "Temps (jours)", y = "Ratio") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5, color = "grey50"))

# Assemblage
p_final <- gridExtra::grid.arrange(
  p_traj, p_resid, p_ratio,
  ncol = 1,
  top = grid::textGrob(
    sprintf("Diagnostic isotype vs contrôle — Fc-silent FGFR2-huBPA-LP1\n%s", verdict),
    gp = grid::gpar(fontface = "bold", fontsize = 11)
  )
)

ggsave("scripts/diag_isotype_vs_ctrl.png", p_final,
       width = 9, height = 13, dpi = 150)
cat("Figure → scripts/diag_isotype_vs_ctrl.png\n")
