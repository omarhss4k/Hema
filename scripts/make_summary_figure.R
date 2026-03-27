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

# Unité µM pour le graphe PK
MW_carbo  <- 371.25
sim$C_uM  <- sim$C1 * 1000 / MW_carbo

# ── Palette & thème ─────────────────────────────────────
th <- theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="#f0f4ff"),
        plot.title=element_text(face="bold", size=12),
        legend.position="bottom")

dose_lines <- geom_vline(xintercept=c(0,21), linetype="dashed",
                         color="grey50", linewidth=0.4)

# ── Panel A : PK ────────────────────────────────────────
pA <- ggplot(sim, aes(days, C_uM)) +
  dose_lines +
  geom_line(color="#2166ac", linewidth=0.9) +
  scale_x_continuous(breaks=seq(0,63,7)) +
  labs(title="A — Carboplatin plasma concentration",
       x="Day", y="Concentration (µM)") +
  th

# ── Panel B : Damage ────────────────────────────────────
pB <- ggplot(sim, aes(days, Damage)) +
  dose_lines +
  geom_line(color="#d7301f", linewidth=0.9) +
  geom_hline(yintercept=1, linetype="dotted", color="grey40") +
  annotate("text", x=35, y=1.05, label="Kill threshold", size=3, color="grey40") +
  scale_x_continuous(breaks=seq(0,63,7)) +
  labs(title="B — DNA damage (PD driver)",
       x="Day", y="Damage (a.u.)") +
  th

# ── Panel C : Neutrophiles ───────────────────────────────
neut_bands <- data.frame(
  ymin=c(0,  0.5, 1.0, 1.5),
  ymax=c(0.5,1.0, 1.5, 2.0),
  grade=c("G4","G3","G2","G1"),
  fill=c("#7b0404","#d73027","#fc8d59","#fee08b")
)
pC <- ggplot() +
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
  labs(title="C — Neutrophil count (ANC)",
       x="Day", y=expression(ANC~(10^9/L))) +
  th

# ── Panel D : Plaquettes ─────────────────────────────────
plt_bands <- data.frame(
  ymin=c(0, 25, 50, 75),
  ymax=c(25,50, 75,150),
  grade=c("G4","G3","G2","G1"),
  fill=c("#7b0404","#d73027","#fc8d59","#fee08b")
)
pD <- ggplot() +
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
  labs(title="D — Platelet count",
       x="Day", y=expression(Platelets~(10^9/L))) +
  th

# ── Panel E : VPC bandes population ─────────────────────
set.seed(42)
vpc_all <- run_vpc(sim, n_sim=500)
vpc_plt <- vpc_all$Plt
vpc_plt$days <- vpc_plt$days

pE <- ggplot(vpc_plt, aes(days)) +
  dose_lines +
  geom_ribbon(aes(ymin=p05,ymax=p95), fill="#1b9e77", alpha=0.20) +
  geom_ribbon(aes(ymin=p05,ymax=p50), fill="#1b9e77", alpha=0.15) +
  geom_line(aes(y=p50), color="#1b9e77", linewidth=1) +
  geom_hline(yintercept=c(25,50,75,150), linetype="dotted",
             color="grey50", linewidth=0.3) +
  annotate("text",x=63,y=c(12,37,62,112),
           label=c("G4","G3","G2","G1"),size=2.8,hjust=1,color="grey30") +
  scale_x_continuous(breaks=seq(0,63,7)) +
  scale_y_continuous(limits=c(0,500)) +
  labs(title="E — Platelet VPC (5th–95th percentile)",
       x="Day",y=expression(Platelets~(10^9/L))) +
  th

vpc_neut <- vpc_all$Neut

pF <- ggplot(vpc_neut, aes(days)) +
  dose_lines +
  geom_ribbon(aes(ymin=p05,ymax=p95), fill="#d7191c", alpha=0.20) +
  geom_ribbon(aes(ymin=p05,ymax=p50), fill="#d7191c", alpha=0.15) +
  geom_line(aes(y=p50), color="#d7191c", linewidth=1) +
  geom_hline(yintercept=c(0.5,1.0,1.5,2.0), linetype="dotted",
             color="grey50", linewidth=0.3) +
  annotate("text",x=63,y=c(0.25,0.75,1.25,1.75),
           label=c("G4","G3","G2","G1"),size=2.8,hjust=1,color="grey30") +
  scale_x_continuous(breaks=seq(0,63,7)) +
  scale_y_continuous(limits=c(0, init_pars$Neut0*1.3)) +
  labs(title="F — Neutrophil VPC (5th–95th percentile)",
       x="Day",y=expression(ANC~(10^9/L))) +
  th

# ── Assemblage ───────────────────────────────────────────
pdf("results_HUMAN/Figure_synthese_tutrice.pdf", width=12, height=10)
grid.arrange(
  pA, pB,
  pC, pD,
  pF, pE,
  ncol=2,
  top=grid::textGrob(
    "QSP Model — Carboplatin Myelotoxicity  |  Fornari 2019  |  AUC=5, Q21D×2",
    gp=grid::gpar(fontsize=14, fontface="bold")
  )
)
dev.off()

cat("✓ Figure sauvegardée : results_HUMAN/Figure_synthese_tutrice.pdf\n")
