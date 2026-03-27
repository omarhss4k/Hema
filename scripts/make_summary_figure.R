#!/usr/bin/env Rscript
# Figure de synthèse — présentation tutrice
library(ggplot2)
library(gridExtra)
library(deSolve)

setwd(dirname(sys.frame(1)$ofile) |> tryCatch(error=function(e) "/home/user/Hema/scripts"))
setwd("/home/user/Hema/scripts")

source("pkpd_model_FORNARI.R")
source("parameters_human.R")
source("parameters_FORNARI_CORRECT.R")
source("vpc.R")
source("plots_grades.R")

AUC_target <- 5
GFR_mLmin  <- 78

# ── Simulation déterministe ──────────────────────────────
times <- seq(0, 63*24, by=1)
init_pars$rate_fun <- make_repeated_infusion(
  dose_mg = AUC_target*(GFR_mLmin+25), Tinfu_h=1,
  interval_h=21*24, n_cycles=2
)

sim <- as.data.frame(lsoda(
  y=init_state, times=times, func=pkpd_fornari,
  parms=init_pars, rtol=1e-4, atol=1e-6, maxsteps=100000
))
sim$days <- sim$time/24

# ── Palette & thème ─────────────────────────────────────
th <- theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="#f0f4ff"),
        plot.title=element_text(face="bold", size=12),
        legend.position="bottom")

dose_lines <- geom_vline(xintercept=c(0,21), linetype="dashed",
                         color="grey50", linewidth=0.4)

# ── Panel A : Neutrophiles ───────────────────────────────
neut_bands <- data.frame(
  ymin=c(0,  0.5, 1.0, 1.5),
  ymax=c(0.5,1.0, 1.5, 2.0),
  grade=c("G4","G3","G2","G1"),
  fill=c("#7b0404","#d73027","#fc8d59","#fee08b")
)
pA <- ggplot() +
  dose_lines +
  geom_rect(data=neut_bands,
            aes(xmin=-Inf,xmax=Inf,ymin=ymin,ymax=ymax,fill=grade),
            alpha=0.15, inherit.aes=FALSE) +
  geom_hline(yintercept=c(0.5,1.0,1.5,2.0), linetype="dotted",
             color="grey60", linewidth=0.3) +
  geom_line(data=sim, aes(days, Neut), color="#d7191c", linewidth=1) +
  scale_fill_manual(values=setNames(neut_bands$fill, neut_bands$grade), guide="none") +
  scale_x_continuous(breaks=seq(0,63,7)) +
  scale_y_continuous(limits=c(0, init_pars$Neut0*1.1)) +
  annotate("text",x=63,y=c(0.25,0.75,1.25,1.75),
           label=c("G4","G3","G2","G1"),size=2.8,hjust=1,color="grey30") +
  labs(title="A — Neutrophil count (ANC)",
       x="Day", y=expression(ANC~(10^9/L))) +
  th

# ── Panel B : Plaquettes ─────────────────────────────────
plt_bands <- data.frame(
  ymin=c(0, 25, 50, 75),
  ymax=c(25,50, 75,150),
  grade=c("G4","G3","G2","G1"),
  fill=c("#7b0404","#d73027","#fc8d59","#fee08b")
)
pB <- ggplot() +
  dose_lines +
  geom_rect(data=plt_bands,
            aes(xmin=-Inf,xmax=Inf,ymin=ymin,ymax=ymax,fill=grade),
            alpha=0.15, inherit.aes=FALSE) +
  geom_hline(yintercept=c(25,50,75,150), linetype="dotted",
             color="grey60", linewidth=0.3) +
  geom_line(data=sim, aes(days, Plt), color="#1b9e77", linewidth=1) +
  scale_fill_manual(values=setNames(plt_bands$fill, plt_bands$grade), guide="none") +
  scale_x_continuous(breaks=seq(0,63,7)) +
  scale_y_continuous(limits=c(0, init_pars$Plt0*1.15)) +
  annotate("text",x=63,y=c(12,37,62,112),
           label=c("G4","G3","G2","G1"),size=2.8,hjust=1,color="grey30") +
  labs(title="B — Platelet count",
       x="Day", y=expression(Platelets~(10^9/L))) +
  th

# ── Panels C/D : VPC bandes population ──────────────────
set.seed(42)
vpc_all  <- run_vpc(sim, n_sim=500)
vpc_plt  <- vpc_all$Plt
vpc_neut <- vpc_all$Neut

pC <- ggplot(vpc_neut, aes(days)) +
  dose_lines +
  geom_rect(data=neut_bands,
            aes(xmin=-Inf,xmax=Inf,ymin=ymin,ymax=ymax,fill=grade),
            alpha=0.12, inherit.aes=FALSE) +
  geom_ribbon(aes(ymin=p05,ymax=p95), fill="#d7191c", alpha=0.25) +
  geom_line(aes(y=p50), color="#d7191c", linewidth=1) +
  geom_hline(yintercept=c(0.5,1.0,1.5,2.0), linetype="dotted",
             color="grey50", linewidth=0.3) +
  annotate("text",x=63,y=c(0.25,0.75,1.25,1.75),
           label=c("G4","G3","G2","G1"),size=2.8,hjust=1,color="grey30") +
  scale_fill_manual(values=setNames(neut_bands$fill, neut_bands$grade), guide="none") +
  scale_x_continuous(breaks=seq(0,63,7)) +
  scale_y_continuous(limits=c(0, init_pars$Neut0*1.3)) +
  labs(title="C — Neutrophil VPC [5th–95th]",
       x="Day", y=expression(ANC~(10^9/L))) +
  th

pD <- ggplot(vpc_plt, aes(days)) +
  dose_lines +
  geom_rect(data=plt_bands,
            aes(xmin=-Inf,xmax=Inf,ymin=ymin,ymax=ymax,fill=grade),
            alpha=0.12, inherit.aes=FALSE) +
  geom_ribbon(aes(ymin=p05,ymax=p95), fill="#1b9e77", alpha=0.25) +
  geom_line(aes(y=p50), color="#1b9e77", linewidth=1) +
  geom_hline(yintercept=c(25,50,75,150), linetype="dotted",
             color="grey50", linewidth=0.3) +
  annotate("text",x=63,y=c(12,37,62,112),
           label=c("G4","G3","G2","G1"),size=2.8,hjust=1,color="grey30") +
  scale_fill_manual(values=setNames(plt_bands$fill, plt_bands$grade), guide="none") +
  scale_x_continuous(breaks=seq(0,63,7)) +
  scale_y_continuous(limits=c(0, 500)) +
  labs(title="D — Platelet VPC [5th–95th]",
       x="Day", y=expression(Platelets~(10^9/L))) +
  th

# ── Panels E/F : Barplot grades population (Figure 4c) ──
set.seed(42)
n_pat     <- 500
omega_CL  <- 0.35; omega_Slope_CMP <- 0.624; omega_Slope_MEP <- 0.547
omega_Slope_MPP <- 0.624; omega_Neut0 <- 0.326; omega_Plt0 <- 0.268
CL_typical <- (125 + 25) * 60 / 1000

eta_CL    <- rnorm(n_pat, 0, omega_CL)
eta_CMP   <- rnorm(n_pat, 0, omega_Slope_CMP)
eta_MEP   <- rnorm(n_pat, 0, omega_Slope_MEP)
eta_MPP   <- rnorm(n_pat, 0, omega_Slope_MPP)
eta_Neut0 <- rnorm(n_pat, 0, omega_Neut0)
eta_Plt0  <- rnorm(n_pat, 0, omega_Plt0)
neut0_i   <- init_pars$Neut0 * exp(eta_Neut0)
plt0_i    <- init_pars$Plt0  * exp(eta_Plt0)

times_g    <- seq(0, 2*21*24 + 21*24, by=4)
grade_neut <- integer(n_pat)
grade_plt  <- integer(n_pat)

cat(sprintf("  → Grades population : %d patients...\n", n_pat))
for (i in seq_len(n_pat)) {
  pars_i <- init_pars
  pars_i$CL        <- CL_typical * exp(eta_CL[i])
  dose_i            <- 5 * (pars_i$CL * 1000/60)
  pars_i$Slope_MPP <- init_pars$Slope_MPP * exp(eta_MPP[i])
  pars_i$Slope_CMP <- init_pars$Slope_CMP * exp(eta_CMP[i])
  pars_i$Slope_MEP <- init_pars$Slope_MEP * exp(eta_MEP[i])
  pars_i$Neut0     <- neut0_i[i]
  pars_i$Plt0      <- plt0_i[i]
  pars_i$rate_fun  <- make_repeated_infusion(dose_i, 1, 21*24, 2)
  a_Neut <- 3/pars_i$MTT_Neut; a_Plt <- 3/pars_i$MTT_Plt
  T_Neut <- pars_i$k_circ_Neut * neut0_i[i] / a_Neut
  T_Plt  <- pars_i$k_circ_Plt  * plt0_i[i]  / a_Plt
  si <- init_state
  si["Neut"] <- neut0_i[i]; si["Plt"] <- plt0_i[i]
  si["T1_Neut"] <- T_Neut; si["T2_Neut"] <- T_Neut; si["T3_Neut"] <- T_Neut
  si["T1_Plt"]  <- T_Plt/pars_i$lambda2
  si["T2_Plt"]  <- T_Plt; si["T3_Plt"] <- T_Plt
  tryCatch({
    out <- as.data.frame(lsoda(si, times_g, pkpd_fornari, pars_i,
                               rtol=1e-3, atol=1e-5, maxsteps=50000))
    grade_neut[i] <- assign_grade_neut(min(out$Neut, na.rm=TRUE))
    grade_plt[i]  <- assign_grade_plt( min(out$Plt,  na.rm=TRUE))
  }, error=function(e) NULL)
  if (i %% 100 == 0) cat(sprintf("    %d/%d\n", i, n_pat))
}

pct_neut <- sapply(1:4, function(g) 100*mean(grade_neut==g, na.rm=TRUE))
pct_plt  <- sapply(1:4, function(g) 100*mean(grade_plt ==g, na.rm=TRUE))
cat(sprintf("  Neut : G1=%.1f%% G2=%.1f%% G3=%.1f%% G4=%.1f%%\n", pct_neut[1],pct_neut[2],pct_neut[3],pct_neut[4]))
cat(sprintf("  Plt  : G1=%.1f%% G2=%.1f%% G3=%.1f%% G4=%.1f%%\n", pct_plt[1], pct_plt[2], pct_plt[3], pct_plt[4]))

GRADE_COLORS <- c(G1="#fee08b", G2="#fc8d59", G3="#d73027", G4="#7b0404")

df_grades <- data.frame(
  Grade   = rep(paste0("G",1:4), 2),
  Pct     = c(pct_neut, pct_plt),
  Lineage = rep(c("Neutropenia","Thrombocytopenia"), each=4)
)
df_grades$Grade   <- factor(df_grades$Grade,   levels=paste0("G",1:4))
df_grades$Lineage <- factor(df_grades$Lineage, levels=c("Neutropenia","Thrombocytopenia"))

pE <- ggplot(df_grades[df_grades$Lineage=="Neutropenia",],
             aes(Grade, Pct, fill=Grade)) +
  geom_col(width=0.65, alpha=0.9) +
  geom_text(aes(label=sprintf("%.1f%%",Pct)), vjust=-0.4, size=3.5) +
  scale_fill_manual(values=GRADE_COLORS, guide="none") +
  scale_y_continuous(limits=c(0, max(pct_neut)*1.25),
                     labels=function(x) paste0(x,"%")) +
  labs(title="E — Neutropenia grade distribution (n=500)",
       subtitle="AUC=5  Q21D×2  IIV: Slope + CL + Neut0",
       x="NCI-CTCAE v5.0 Grade", y="% patients") +
  th

pF <- ggplot(df_grades[df_grades$Lineage=="Thrombocytopenia",],
             aes(Grade, Pct, fill=Grade)) +
  geom_col(width=0.65, alpha=0.9) +
  geom_text(aes(label=sprintf("%.1f%%",Pct)), vjust=-0.4, size=3.5) +
  scale_fill_manual(values=GRADE_COLORS, guide="none") +
  scale_y_continuous(limits=c(0, max(pct_plt)*1.25),
                     labels=function(x) paste0(x,"%")) +
  labs(title="F — Thrombocytopenia grade distribution (n=500)",
       subtitle="AUC=5  Q21D×2  IIV: Slope_MEP + CL + Plt0",
       x="NCI-CTCAE v5.0 Grade", y="% patients") +
  th

# ── Assemblage ───────────────────────────────────────────
pdf("results_HUMAN/Figure_synthese_tutrice.pdf", width=12, height=13)
grid.arrange(
  pA, pB,
  pC, pD,
  pE, pF,
  ncol=2,
  top=grid::textGrob(
    "QSP Model — Carboplatin Myelotoxicity  |  Fornari 2019  |  AUC=5, Q21D×2",
    gp=grid::gpar(fontsize=14, fontface="bold")
  )
)
dev.off()

cat("✓ Figure sauvegardée : results_HUMAN/Figure_synthese_tutrice.pdf\n")
