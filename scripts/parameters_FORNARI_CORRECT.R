############################################################
# parameters_FORNARI_CORRECT.R
# Paramètres DÉRIVÉS — Equation S4 exacte (Fornari 2019)
# Sources parameters_rat.R
############################################################

source("parameters_human.R")

cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  PARAMÈTRES FORNARI — EQUATION S4 EXACTE               ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

ip <- init_pars

# ══════════════════════════════════════════════════════════
# PARAMÈTRES λ (Table 1)
# ══════════════════════════════════════════════════════════
ip$lambda1 <- 2.0   # Ret transit amplification
ip$lambda2 <- 2.0   # Plt transit amplification
ip$lambda3 <- 1.8   # CMP proliferation rate
ip$lambda4 <- 1.8   # MEP proliferation rate
ip$lambda5 <- 1.8   # MPP proliferation rate

cat("── λ (Table 1) ──\n")
cat(sprintf("  λ1=%.0f  λ2=%.0f  λ3=%.1f  λ4=%.1f  λ5=%.1f\n\n",
            ip$lambda1, ip$lambda2, ip$lambda3, ip$lambda4, ip$lambda5))

# ══════════════════════════════════════════════════════════
# 1. k_circ_Ret (Eq S4)
# ══════════════════════════════════════════════════════════
ip$k_circ_Ret <- ip$k_circ_RBC * ip$RBC0 / ip$Ret0

cat("── 1. k_circ_Ret (Eq S4) ──\n")
cat(sprintf("  k_circ_Ret = k_circ_RBC × RBC0/Ret0\n"))
cat(sprintf("             = %.4f × %.0f/%.0f = %.5f /h  (lifespan = %.1fh)\n\n",
            ip$k_circ_RBC, ip$RBC0, ip$Ret0,
            ip$k_circ_Ret, 1/ip$k_circ_Ret))

# ══════════════════════════════════════════════════════════
# 2. Constantes a (maturation transit)
# ══════════════════════════════════════════════════════════
ip$a_Ret <- 3 / ip$MTT_Ret
ip$a_Plt <- 3 / ip$MTT_Plt

# ══════════════════════════════════════════════════════════
# 3. k_tr_Neut, k_tr_Mono (Eq S4)
# ══════════════════════════════════════════════════════════
ip$k_tr_Neut <- ip$k_circ_Neut * ip$Neut0 / ip$CMP0
ip$k_tr_Mono <- ip$k_circ_Mono * ip$Mono0 / ip$CMP0

cat("── 2. k_tr_Neut, k_tr_Mono (Eq S4) ──\n")
cat(sprintf("  k_tr_Neut = %.6f /h\n", ip$k_tr_Neut))
cat(sprintf("  k_tr_Mono = %.6f /h\n\n", ip$k_tr_Mono))

# ══════════════════════════════════════════════════════════
# 4. k_tr_Ret, k_tr_Plt (Eq S4) ← CORRECTION CRITIQUE
#    k_tr = a × (T1_0/MEP0) × (T1_0/T_0)
#    Simplifié : k_tr = k_circ × Cell0 / (λ² × MEP0)
# ══════════════════════════════════════════════════════════
ip$k_tr_Ret <- ip$k_circ_Ret * ip$Ret0 / (ip$lambda1^2 * ip$MEP0)
ip$k_tr_Plt <- ip$k_circ_Plt * ip$Plt0 / (ip$lambda2^2 * ip$MEP0)

cat("── 3. k_tr_Ret, k_tr_Plt (Eq S4) ← CORRECTION: ÷λ² ──\n")
cat(sprintf("  k_tr_Ret = k_circ_Ret × Ret0 / (λ1² × MEP0)\n"))
cat(sprintf("           = %.5f × %.0f / (%.0f² × %.0f) = %.6f /h\n",
            ip$k_circ_Ret, ip$Ret0, ip$lambda1, ip$MEP0, ip$k_tr_Ret))
cat(sprintf("  k_tr_Plt = %.6f /h\n\n", ip$k_tr_Plt))

# ══════════════════════════════════════════════════════════
# 5. k_prol_Ret, k_prol_Plt (Eq S4)
# ══════════════════════════════════════════════════════════
ip$k_prol_Ret <- ip$a_Ret * (1 - 1/ip$lambda1)
ip$k_prol_Plt <- ip$a_Plt * (1 - 1/ip$lambda2)

cat("── 4. k_prol_Ret, k_prol_Plt (Eq S4) ──\n")
cat(sprintf("  k_prol_Ret = a_Ret × (1 - 1/λ1) = %.6f /h\n", ip$k_prol_Ret))
cat(sprintf("  k_prol_Plt = a_Plt × (1 - 1/λ2) = %.6f /h\n\n", ip$k_prol_Plt))

# ══════════════════════════════════════════════════════════
# 6. k_prol_CMP, k_prol_MEP (Eq S4)
# ══════════════════════════════════════════════════════════
ip$k_prol_CMP <- (ip$k_tr_Neut + ip$k_tr_Mono) / ip$lambda3
ip$k_prol_MEP <- (ip$k_tr_Ret  + ip$k_tr_Plt)  / ip$lambda4

cat("── 5. k_prol_CMP, k_prol_MEP (Eq S4) ──\n")
cat(sprintf("  k_prol_CMP = (k_tr_Neut + k_tr_Mono) / λ3 = %.6f /h\n", ip$k_prol_CMP))
cat(sprintf("  k_prol_MEP = (k_tr_Ret  + k_tr_Plt)  / λ4 = %.6f /h\n\n", ip$k_prol_MEP))

# ══════════════════════════════════════════════════════════
# 7. k_tr_CMP, k_tr_MEP (Eq S4)
# ══════════════════════════════════════════════════════════
ip$k_tr_CMP <- (ip$k_tr_Neut + ip$k_tr_Mono - ip$k_prol_CMP) * ip$CMP0 / ip$MPP0
ip$k_tr_MEP <- (ip$k_tr_Ret  + ip$k_tr_Plt  - ip$k_prol_MEP) * ip$MEP0 / ip$MPP0

cat("── 6. k_tr_CMP, k_tr_MEP (Eq S4) ──\n")
cat(sprintf("  k_tr_CMP = %.6f /h\n", ip$k_tr_CMP))
cat(sprintf("  k_tr_MEP = %.6f /h\n\n", ip$k_tr_MEP))

# ══════════════════════════════════════════════════════════
# 8. k_prol_MPP (Eq S4)
# ══════════════════════════════════════════════════════════
ip$k_prol_MPP <- (ip$k_tr_CMP + ip$k_tr_MEP) / ip$lambda5

cat("── 7. k_prol_MPP (Eq S4) ──\n")
cat(sprintf("  k_prol_MPP = (k_tr_CMP + k_tr_MEP) / λ5 = %.6f /h\n\n", ip$k_prol_MPP))

# ══════════════════════════════════════════════════════════
# 9. k_stem (Eq S4)
# ══════════════════════════════════════════════════════════
ip$k_stem <- (ip$k_tr_CMP + ip$k_tr_MEP - ip$k_prol_MPP) * ip$MPP0

cat("── 8. k_stem (Eq S4) ──\n")
cat(sprintf("  k_stem = (k_tr_CMP + k_tr_MEP - k_prol_MPP) × MPP0 = %.6f\n\n", ip$k_stem))

# ══════════════════════════════════════════════════════════
# ÉTATS INITIAUX TRANSIT (Equation S3)
# ══════════════════════════════════════════════════════════
a_Neut <- 3 / ip$MTT_Neut
a_Mono <- 3 / ip$MTT_Mono

T_Neut  <- ip$k_circ_Neut * ip$Neut0 / a_Neut
T_Mono  <- ip$k_circ_Mono * ip$Mono0 / a_Mono
T2_Ret  <- ip$k_circ_Ret  * ip$Ret0  / ip$a_Ret
T1_Ret  <- T2_Ret / ip$lambda1
T2_Plt  <- ip$k_circ_Plt  * ip$Plt0  / ip$a_Plt
T1_Plt  <- T2_Plt / ip$lambda2

cat("── États initiaux transit (Eq S3) ──\n")
cat(sprintf("  T1_Ret = T2_Ret/λ1 = %.4f\n",  T1_Ret))
cat(sprintf("  T2_Ret = T3_Ret    = %.4f\n",  T2_Ret))
cat(sprintf("  T1_Plt = T2_Plt/λ2 = %.4f\n",  T1_Plt))
cat(sprintf("  T2_Plt = T3_Plt    = %.4f\n\n", T2_Plt))

# ── Mise à jour init_state ──
is <- init_state
is["T1_Neut"] <- T_Neut;  is["T2_Neut"] <- T_Neut;  is["T3_Neut"] <- T_Neut
is["T1_Mono"] <- T_Mono;  is["T2_Mono"] <- T_Mono;  is["T3_Mono"] <- T_Mono
is["T1_Ret"]  <- T1_Ret;  is["T2_Ret"]  <- T2_Ret;  is["T3_Ret"]  <- T2_Ret
is["T1_Plt"]  <- T1_Plt;  is["T2_Plt"]  <- T2_Plt;  is["T3_Plt"]  <- T2_Plt

cat("═══════════════════════════════════════════════════════════\n")
cat("  TOUTES LES FORMULES EQ S4 APPLIQUÉES ✓\n")
cat("═══════════════════════════════════════════════════════════\n\n")

init_pars  <- ip
init_state <- is
