# =============================================================================
# PRÉSENTATION TUTRICE — Fc-silent FGFR2-huBPA-LP1  (Modèle Simeoni 2004)
# Génère 9 figures PNG prêtes à intégrer dans une présentation :
#   01_donnees_PK_brutes.png
#   02_PK_predit_vs_observe.png
#   03_parametres_PK.png
#   04_schema_modele_PK.png
#   05_croissance_vs_inhibition.png
#   06_schema_modele_Simeoni.png
#   07_PKPD_predit_vs_observe.png
#   08_tableau_TGI.png
#   09_parametres_PKPD.png
# =============================================================================

library(rxode2)
library(deSolve)
library(ggplot2)
library(readxl)
library(gridExtra)

OUT <- "scripts/presentation_simeoni"
dir.create(OUT, showWarnings = FALSE)

SAVE <- function(name, p, w = 10, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".png")), p,
         width = w, height = h, dpi = 180, bg = "white")
  cat("  →", file.path(OUT, paste0(name, ".png")), "\n")
}

COLS <- c("10 mg/kg" = "#1B4F9E",
          "3 mg/kg"  = "#27AE60",
          "Contrôle" = "#888888")

THEME <- theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        legend.background = element_rect(fill = "white", color = NA))

PSI <- 20   # fixé — Simeoni 2004

# =============================================================================
# CHARGEMENT DES RÉSULTATS
# =============================================================================

load("scripts/resultats_PK2comp_rxode2_FGFR2.RData")   # → pk2comp_rxode2
load("scripts/resultats_PKPD_simeoni_FGFR2.RData")      # → simeoni_results

# Paramètres PK (jours)
PK <- c(
  CL = unname(pk2comp_rxode2$CL) * 24,
  V1 = unname(pk2comp_rxode2$V1),
  V2 = unname(pk2comp_rxode2$V2),
  Q  = unname(pk2comp_rxode2$Q)  * 24
)

# Paramètres PD Simeoni
L0  <- simeoni_results$L0
L1  <- simeoni_results$L1
k1  <- simeoni_results$k1
k2  <- simeoni_results$k2
TSC <- simeoni_results$TSC_ugL
MTT <- simeoni_results$MTT_days

DOSE_DAYS <- c(0, 14, 28, 42)

# =============================================================================
# LECTURE DONNÉES TUMORALES (Excel)
# =============================================================================

read_tumor_excel <- function(path) {
  raw   <- suppressMessages(read_xlsx(path, col_names=FALSE, .name_repair="minimal"))
  raw   <- raw[, colSums(!is.na(raw)) > 0]
  hidx  <- which(raw[[1]] == "Group")
  times <- suppressWarnings(as.numeric(as.character(unlist(raw[hidx[1], -1]))))
  valid <- !is.na(times); times <- times[valid]
  extr  <- function(s, e) {
    b   <- raw[s:e, ]; b <- b[!is.na(b[[1]]) & b[[1]] != "", ]
    dat <- b[, c(TRUE, valid)]
    for (j in 2:ncol(dat)) dat[[j]] <- as.numeric(dat[[j]])
    dat
  }
  list(times    = times,
       blk_mean = extr(hidx[1]+1, hidx[2]-1),
       blk_sem  = extr(hidx[2]+1, nrow(raw)))
}

.gr <- function(blk, pat)
  as.numeric(unlist(blk[grep(pat, blk[[1]], ignore.case=TRUE)[1], -1]))

tv      <- read_tumor_excel("TumorVolume_FGFR2.xlsx")
times_d <- tv$times

w_ctrl   <- .gr(tv$blk_mean, "Group 01")
w_d10_tv <- .gr(tv$blk_mean, "Group 03")
w_d3_tv  <- .gr(tv$blk_mean, "Group 04")
sem_ctrl <- pmax(.gr(tv$blk_sem, "Group 01"), 1)
sem_d10  <- pmax(.gr(tv$blk_sem, "Group 03"), 1)
sem_d3   <- pmax(.gr(tv$blk_sem, "Group 04"), 1)

get_tv0 <- function(tv_vec) {
  v <- tv_vec[times_d == 0]
  if (length(v) == 0 || is.na(v[1])) tv_vec[!is.na(tv_vec)][1] else v[1]
}
tv0_ctrl <- get_tv0(w_ctrl)
tv0_d3   <- get_tv0(w_d3_tv)
tv0_d10  <- get_tv0(w_d10_tv)

# Données obs formatées
lev_pd <- c("Contrôle", "3 mg/kg", "10 mg/kg")

ok_ctrl <- !is.na(w_ctrl)
ok_d3   <- !is.na(w_d3_tv)
ok_d10  <- !is.na(w_d10_tv)

df_obs_pd <- rbind(
  data.frame(jour=times_d[ok_ctrl], w=w_ctrl[ok_ctrl],
             sem=sem_ctrl[ok_ctrl], Groupe="Contrôle"),
  data.frame(jour=times_d[ok_d3],   w=w_d3_tv[ok_d3],
             sem=sem_d3[ok_d3],     Groupe="3 mg/kg"),
  data.frame(jour=times_d[ok_d10],  w=w_d10_tv[ok_d10],
             sem=sem_d10[ok_d10],   Groupe="10 mg/kg")
)
df_obs_pd$Groupe <- factor(df_obs_pd$Groupe, levels=lev_pd)

# =============================================================================
# MODÈLES
# =============================================================================

# ---- PK 2-compartiments ----
mod2comp <- rxode2({
  C1       <- A1 / V1
  d/dt(A1) <- -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2
  d/dt(A2) <-  (Q/V1)*A1 - (Q/V2)*A2
})

sim_pk <- function(dose, params, times_h) {
  ev <- eventTable()
  ev$add.dosing(dose=dose, nbr.doses=1, dosing.to=1)
  ev$add.sampling(sort(unique(c(0, times_h))))
  out <- tryCatch(rxSolve(mod2comp, params, ev), error=function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, length(times_h)))
  approx(out$time, out$C1, xout=times_h, rule=2)$y
}

# ---- Simeoni ODE (jours, µg/L) ----
simeoni_rhs_pres <- function(t, state, pars) {
  A1 <- max(state[1], 0); A2 <- max(state[2], 0)
  x1 <- max(state[3], 0); x2 <- max(state[4], 0)
  x3 <- max(state[5], 0); x4 <- max(state[6], 0)

  CL <- pars$CL; V1 <- pars$V1; V2 <- pars$V2; Q <- pars$Q
  L0 <- pars$L0; L1 <- pars$L1; k1 <- pars$k1; k2 <- pars$k2

  C1  <- A1 / V1
  w   <- x1 + x2 + x3 + x4
  gw  <- L0 * w / (1 + (L0 * w / L1)^PSI)^(1/PSI)
  gr  <- if (w > 1e-12) gw/w else L0

  list(c(
    -(CL/V1 + Q/V1)*A1 + (Q/V2)*A2,
     (Q/V1)*A1 - (Q/V2)*A2,
    (gr - k2*C1)*x1,
    k2*C1*x1 - k1*x2,
    k1*x2 - k1*x3,
    k1*x3 - k1*x4
  ))
}

simeoni_ctrl_rhs <- function(t, state, pars) {
  x1 <- max(state[1], 0)
  w  <- x1 + max(state[2],0) + max(state[3],0) + max(state[4],0)
  gw <- pars$L0 * w / (1 + (pars$L0*w/pars$L1)^PSI)^(1/PSI)
  gr <- if (w > 1e-12) gw/w else pars$L0
  list(c(gr*x1, 0, 0, 0))
}

# Simulation contrôle
sim_ctrl_fn <- function(tv0, pars_list, times_out) {
  t_all <- sort(unique(c(0, times_out)))
  out <- tryCatch(
    as.data.frame(lsoda(c(x1=tv0, x2=0, x3=0, x4=0), t_all,
                        simeoni_ctrl_rhs, pars_list, atol=1e-6, rtol=1e-6)),
    error=function(e) NULL)
  if (is.null(out)) return(rep(NA_real_, length(times_out)))
  w <- pmax(out$x1,0)+pmax(out$x2,0)+pmax(out$x3,0)+pmax(out$x4,0)
  approx(out$time, w, xout=times_out, rule=2)$y
}

# Simulation traitée — redémarrage à chaque dose
sim_simeoni <- function(dose_ugkg, tv0, pars_list, times_out) {
  breaks <- c(DOSE_DAYS, max(times_out)+1)
  state  <- c(A1=dose_ugkg, A2=0, x1=tv0, x2=0, x3=0, x4=0)
  t_all <- numeric(0); w_all <- numeric(0)

  for (i in seq_along(DOSE_DAYS)) {
    t_start <- DOSE_DAYS[i]
    t_stop  <- min(breaks[i+1], max(times_out))
    if (t_start >= max(times_out)) break
    if (i > 1) state["A1"] <- state["A1"] + dose_ugkg

    t_seg <- sort(unique(c(t_start,
                            times_out[times_out > t_start & times_out <= t_stop],
                            t_stop)))
    if (length(t_seg) < 2) t_seg <- c(t_start, t_stop)

    seg <- tryCatch(
      as.data.frame(lsoda(state, t_seg, simeoni_rhs_pres, pars_list,
                          atol=1e-6, rtol=1e-6)),
      error=function(e) NULL)
    if (is.null(seg)) return(rep(NA_real_, length(times_out)))

    w_seg <- pmax(seg$x1,0)+pmax(seg$x2,0)+pmax(seg$x3,0)+pmax(seg$x4,0)
    keep  <- if (length(t_all)>0) seg$time > tail(t_all,1) else rep(TRUE,nrow(seg))
    t_all <- c(t_all, seg$time[keep]); w_all <- c(w_all, w_seg[keep])

    last  <- seg[nrow(seg),]
    state <- c(A1=last$A1, A2=last$A2,
               x1=last$x1, x2=last$x2, x3=last$x3, x4=last$x4)
  }
  if (length(t_all)==0) return(rep(NA_real_, length(times_out)))
  approx(t_all, w_all, xout=times_out, rule=2)$y
}

pars_sim <- list(CL=PK["CL"], V1=PK["V1"], V2=PK["V2"], Q=PK["Q"],
                 L0=L0, L1=L1, k1=k1, k2=k2)

# =============================================================================
# FIG 1 — DONNÉES PK BRUTES
# =============================================================================
cat("\n[1/9] Données PK brutes\n")

raw_pk <- read_xlsx("PK souris FGFR2.xlsx")
time_h <- as.numeric(raw_pk[[1]])

.mk <- function(t, c, lbl) {
  ok <- !is.na(c) & c > 0
  data.frame(Temps_h=t[ok], Concentration=c[ok], Dose=lbl)
}
df_pk_raw <- rbind(
  .mk(time_h, as.numeric(raw_pk[[2]]), "10 mg/kg"),
  .mk(time_h, as.numeric(raw_pk[[3]]), "5 mg/kg"),
  .mk(time_h, as.numeric(raw_pk[[4]]), "1 mg/kg")
)
COLS_PK <- c("10 mg/kg"="#1B4F9E", "5 mg/kg"="#2E86C1", "1 mg/kg"="#85C1E9")
df_pk_raw$Dose <- factor(df_pk_raw$Dose, levels=names(COLS_PK))

p1a <- ggplot(df_pk_raw, aes(x=Temps_h, y=Concentration, color=Dose)) +
  geom_point(size=3, alpha=0.9) +
  geom_line(aes(group=Dose), linewidth=0.6, linetype="dotted") +
  scale_color_manual(values=COLS_PK) +
  labs(title="Données PK observées — Échelle linéaire",
       x="Temps (heures)", y="Concentration (µg/L)", color=NULL) + THEME

p1b <- ggplot(df_pk_raw, aes(x=Temps_h, y=Concentration, color=Dose)) +
  geom_point(size=3, alpha=0.9) +
  geom_line(aes(group=Dose), linewidth=0.6, linetype="dotted") +
  scale_y_log10() +
  scale_color_manual(values=COLS_PK) +
  labs(title="Données PK observées — Échelle log",
       x="Temps (heures)", y="Concentration log(µg/L)", color=NULL) + THEME

png(file.path(OUT,"01_donnees_PK_brutes.png"),
    width=13, height=5, units="in", res=180, bg="white")
grid.arrange(p1a+theme(legend.position="none"),
             p1b+theme(legend.position="bottom"), ncol=2)
dev.off()
cat("  →", file.path(OUT,"01_donnees_PK_brutes.png"),"\n")

# =============================================================================
# FIG 2 — AJUSTEMENT PK PRÉDIT VS OBSERVÉ
# =============================================================================
cat("[2/9] Ajustement PK\n")

# PK en heures pour la figure
PK_h <- c(CL=unname(pk2comp_rxode2$CL), V1=unname(pk2comp_rxode2$V1),
           V2=unname(pk2comp_rxode2$V2), Q=unname(pk2comp_rxode2$Q))
t_pk <- seq(0, max(df_pk_raw$Temps_h)*1.05, by=1)

df_pk_sim <- do.call(rbind, list(
  data.frame(Temps_h=t_pk, Concentration=sim_pk(10000, PK_h, t_pk), Dose="10 mg/kg"),
  data.frame(Temps_h=t_pk, Concentration=sim_pk(5000,  PK_h, t_pk), Dose="5 mg/kg"),
  data.frame(Temps_h=t_pk, Concentration=sim_pk(1000,  PK_h, t_pk), Dose="1 mg/kg")
))
df_pk_sim$Dose <- factor(df_pk_sim$Dose, levels=names(COLS_PK))

p2 <- ggplot() +
  geom_line(data=df_pk_sim, aes(x=Temps_h, y=Concentration, color=Dose), linewidth=1.2) +
  geom_point(data=df_pk_raw, aes(x=Temps_h, y=Concentration, color=Dose), size=3) +
  scale_y_log10() +
  scale_color_manual(values=COLS_PK) +
  labs(title="Modèle PK 2-compartiments — Prédit vs Observé",
       subtitle=paste0(
         "V1=", round(PK_h["V1"],3), " L/kg   ",
         "V2=", round(PK_h["V2"],3), " L/kg   ",
         "CL=", round(PK_h["CL"],4), " L/h/kg   ",
         "t½β=", round(pk2comp_rxode2$t_half_beta/24,1), " j"),
       x="Temps (heures)", y="Concentration (µg/L — échelle log)", color=NULL,
       caption="Lignes = modèle   ● = données observées") +
  THEME

SAVE("02_PK_predit_vs_observe", p2, w=9, h=5.5)

# =============================================================================
# FIG 3 — TABLEAU PARAMÈTRES PK
# =============================================================================
cat("[3/9] Tableau paramètres PK\n")

pk_tbl <- data.frame(
  Paramètre   = c("CL","V1","V2","Vss","Q","k10","k12","k21","α","β","t½α","t½β"),
  Description = c("Clairance systémique","Volume central","Volume périphérique",
                  "Volume distribution ss","Clairance intercompart.",
                  "Constante élimination","Transfert C1→C2","Transfert C2→C1",
                  "Constante phase rapide","Constante phase lente",
                  "Demi-vie distribution","Demi-vie élimination"),
  Valeur      = c(
    sprintf("%.4f L/h/kg  (%.3f L/j/kg)", PK_h["CL"], PK["CL"]),
    sprintf("%.4f L/kg", PK_h["V1"]),
    sprintf("%.4f L/kg", PK_h["V2"]),
    sprintf("%.4f L/kg", pk2comp_rxode2$Vss),
    sprintf("%.4f L/h/kg", PK_h["Q"]),
    sprintf("%.5f /h", pk2comp_rxode2$k10),
    sprintf("%.5f /h", pk2comp_rxode2$k12),
    sprintf("%.5f /h", pk2comp_rxode2$k21),
    sprintf("%.5f /h", pk2comp_rxode2$alpha),
    sprintf("%.6f /h", pk2comp_rxode2$beta),
    sprintf("%.1f h", pk2comp_rxode2$t_half_alpha),
    sprintf("%.1f h  (%.1f j)", pk2comp_rxode2$t_half_beta, pk2comp_rxode2$t_half_beta/24)
  )
)

nr <- nrow(pk_tbl); nc <- 3
cells3 <- do.call(rbind, lapply(seq_len(nr), function(i) {
  bg  <- if (i %% 2 == 0) "#EBF5FB" else "white"
  hl  <- i > 9
  data.frame(row=i, col=1:nc,
             lbl=c(pk_tbl$Paramètre[i], pk_tbl$Description[i], pk_tbl$Valeur[i]),
             fill=if(hl) "#D6EAF8" else bg,
             bold=hl)
}))
hdr3 <- data.frame(row=0, col=1:nc, lbl=c("Paramètre","Description","Valeur"),
                   fill="#1A5276", bold=TRUE)
cells3 <- rbind(hdr3, cells3)
cells3$y    <- max(cells3$row) - cells3$row
cells3$tcol <- ifelse(cells3$row==0, "white", "black")
cells3$face <- ifelse(cells3$bold, "bold", "plain")
cells3$xpos <- c(0.45, 1.55, 3.15)[cells3$col]
cells3$twd  <- c(0.75, 1.75, 1.65)[cells3$col]

p3 <- ggplot(cells3, aes(x=xpos, y=y)) +
  geom_tile(aes(fill=fill, width=twd), color="grey80", linewidth=0.3) +
  geom_text(aes(label=lbl, color=tcol, fontface=face), size=3.3, hjust=0.5) +
  scale_fill_identity() + scale_color_identity() +
  labs(title="Paramètres PK — Modèle 2-compartiments IV bolus",
       subtitle="Fc-silent FGFR2-huBPA-LP1 | Souris | Optimisation nlminb (log-espace)") +
  theme_void(base_size=11) +
  theme(plot.title=element_text(face="bold", size=13, margin=margin(b=4)),
        plot.subtitle=element_text(size=9, color="grey40", margin=margin(b=8)),
        plot.margin=margin(12,20,12,20))

SAVE("03_parametres_PK", p3, w=11, h=5.5)

# =============================================================================
# FIG 4 — SCHÉMA MODÈLE PK
# =============================================================================
cat("[4/9] Schéma PK\n")

.seg <- function(x1,y1,x2,y2)
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),
               arrow=grid::arrow(length=unit(0.3,"cm"), type="closed"),
               linewidth=1, color="#2C3E50")
.box <- function(xmin,xmax,ymin,ymax,fill="#D6EAF8",col="#2980B9")
  annotate("rect", xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
           fill=fill, color=col, linewidth=1.2)
.lbl <- function(x,y,txt,sz=4.5,bold=FALSE,col="black")
  annotate("text", x=x,y=y,label=txt,size=sz,
           fontface=if(bold)"bold"else"plain",color=col,hjust=0.5)

p4 <- ggplot() + xlim(0,10) + ylim(0,6) +
  .box(0.2,1.8,3.8,5.2,fill="#FDEBD0",col="#E67E22") +
  .lbl(1,4.8,"Dose IV bolus",sz=4,bold=TRUE,col="#E67E22") +
  .lbl(1,4.3,"D = dose × kg",sz=3.5,col="#7F8C8D") +
  .box(2.8,6.2,3.2,5.8,fill="#D6EAF8",col="#2980B9") +
  .lbl(4.5,5.1,"Compartiment central",sz=4.5,bold=TRUE,col="#1A5276") +
  .lbl(4.5,4.6,"A1  (µg/kg)",sz=3.8,col="#2C3E50") +
  .lbl(4.5,4.0,"C1 = A1/V1  (µg/L)",sz=3.8,col="#2C3E50") +
  .box(2.8,6.2,0.5,2.8,fill="#D5F5E3",col="#27AE60") +
  .lbl(4.5,2.1,"Compartiment périphérique",sz=4.5,bold=TRUE,col="#1E8449") +
  .lbl(4.5,1.6,"A2  (µg/kg)",sz=3.8,col="#2C3E50") +
  .box(7.2,9.8,3.8,5.2,fill="#FADBD8",col="#E74C3C") +
  .lbl(8.5,4.8,"Élimination",sz=4,bold=TRUE,col="#C0392B") +
  .lbl(8.5,4.3,"CL · C1",sz=3.5,col="#7F8C8D") +
  .seg(1.8,4.5,2.8,4.5) + .seg(6.2,4.5,7.2,4.5) +
  .seg(4.5,3.2,4.5,2.8) + .seg(4.9,2.8,4.9,3.2) +
  .lbl(3.9,3.0,"Q/V1",sz=3.5,col="#27AE60") +
  .lbl(5.5,3.0,"Q/V2",sz=3.5,col="#27AE60") +
  annotate("label", x=4.5, y=0.25,
           label="dA1/dt = -(CL/V1+Q/V1)·A1 + (Q/V2)·A2\ndA2/dt =  (Q/V1)·A1 - (Q/V2)·A2",
           size=3.5, hjust=0.5, fill="#F8F9FA", color="#2C3E50",
           label.padding=unit(0.4,"lines")) +
  labs(title="Structure du modèle PK — 2 compartiments (IV bolus)",
       subtitle="Paramètres : CL, V1, V2, Q  |  Optimisation nlminb en log-espace") +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold",size=14,margin=margin(b=4)),
        plot.subtitle=element_text(size=9,color="grey40",margin=margin(b=6)),
        plot.margin=margin(10,10,10,10))

SAVE("04_schema_modele_PK", p4, w=11, h=6)

# =============================================================================
# FIG 5 — CROISSANCE TUMORALE VS INHIBITION
# =============================================================================
cat("[5/9] Croissance tumorale\n")

df_tv <- rbind(
  data.frame(jour=times_d, mean=w_ctrl,   sem=sem_ctrl, Groupe="Contrôle"),
  data.frame(jour=times_d, mean=w_d10_tv, sem=sem_d10,  Groupe="10 mg/kg"),
  data.frame(jour=times_d, mean=w_d3_tv,  sem=sem_d3,   Groupe="3 mg/kg")
)
df_tv <- df_tv[!is.na(df_tv$mean), ]
df_tv$Groupe <- factor(df_tv$Groupe, levels=lev_pd)
ymax <- max(df_tv$mean + df_tv$sem, na.rm=TRUE)

p5 <- ggplot(df_tv, aes(x=jour, y=mean, color=Groupe, group=Groupe)) +
  geom_vline(xintercept=DOSE_DAYS, linetype="dashed", color="grey70", linewidth=0.4) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem, fill=Groupe), alpha=0.15, color=NA) +
  geom_line(linewidth=1.2) + geom_point(size=2.8) +
  annotate("point", x=DOSE_DAYS, y=-ymax*0.07, shape=17, size=4, color="#C0392B") +
  annotate("text", x=max(DOSE_DAYS)+1.5, y=-ymax*0.07,
           label="= Traitement Q2W", hjust=0, size=3.2, color="#C0392B") +
  scale_color_manual(values=COLS) + scale_fill_manual(values=COLS) +
  scale_x_continuous(breaks=seq(0,max(df_tv$jour),by=7)) +
  coord_cartesian(ylim=c(-ymax*0.14, ymax*1.1), clip="off") +
  labs(title="Croissance tumorale vs Inhibition par le traitement",
       subtitle="Fc-silent FGFR2-huBPA-LP1 | mean ± SEM | IV Q2W × 4",
       x="Temps (jours)", y="Volume tumoral (mm³)", color=NULL, fill=NULL) +
  THEME + theme(legend.position="bottom", plot.margin=margin(5,10,30,5))

SAVE("05_croissance_vs_inhibition", p5, w=10, h=6)

# =============================================================================
# FIG 6 — SCHÉMA MODÈLE SIMEONI (effet linéaire k2·C)
# =============================================================================
cat("[6/9] Schéma modèle Simeoni\n")

p6 <- ggplot() + xlim(0,14) + ylim(0,8) +

  # Bloc PK
  annotate("rect",xmin=0.2,xmax=3.8,ymin=5.5,ymax=7.5,
           fill="#D6EAF8",color="#2980B9",linewidth=1.2) +
  .lbl(2,6.9,"PK 2-compartiments",sz=4.2,bold=TRUE,col="#1A5276") +
  .lbl(2,6.4,"C1(t) = A1(t) / V1",sz=3.6,col="#2C3E50") +
  .lbl(2,5.9,"CL, V1, V2, Q  (fixés)",sz=3.4,col="#7F8C8D") +
  .seg(3.8,6.5,5.2,4.7) +
  .lbl(4.8,5.8,"C1(t)",sz=3.5,col="#E74C3C") +

  # x1 Proliférant
  annotate("rect",xmin=5.0,xmax=8.0,ymin=3.5,ymax=5.5,
           fill="#D5F5E3",color="#27AE60",linewidth=1.2) +
  .lbl(6.5,5.1,"x1",sz=5.5,bold=TRUE,col="#1E8449") +
  .lbl(6.5,4.6,"Prolifération",sz=3.8,col="#2C3E50") +
  .lbl(6.5,4.1,"g(w)/w · x1",sz=3.4,col="#1E8449") +
  .lbl(6.5,3.7,"− k2·C1·x1",sz=3.4,col="#E74C3C") +
  .seg(8.0,4.5,9.0,4.5) +
  .lbl(8.5,4.9,"k2·C1",sz=3.4,bold=TRUE,col="#E74C3C") +

  # x2, x3, x4 Transit
  annotate("rect",xmin=9.0,xmax=10.3,ymin=3.5,ymax=5.5,
           fill="#FDEBD0",color="#E67E22",linewidth=1.0) +
  .lbl(9.65,4.5,"x2",sz=4.5,bold=TRUE,col="#784212") +
  .seg(10.3,4.5,10.8,4.5) + .lbl(10.55,4.85,"k1",sz=3,col="#784212") +

  annotate("rect",xmin=10.8,xmax=11.9,ymin=3.5,ymax=5.5,
           fill="#FDEBD0",color="#E67E22",linewidth=1.0) +
  .lbl(11.35,4.5,"x3",sz=4.5,bold=TRUE,col="#784212") +
  .seg(11.9,4.5,12.4,4.5) + .lbl(12.15,4.85,"k1",sz=3,col="#784212") +

  annotate("rect",xmin=12.4,xmax=13.5,ymin=3.5,ymax=5.5,
           fill="#FADBD8",color="#E74C3C",linewidth=1.0) +
  .lbl(12.95,4.5,"x4",sz=4.5,bold=TRUE,col="#922B21") +
  .lbl(12.95,3.7,"k1·x4 → mort",sz=3,col="#922B21") +

  # Volume total
  annotate("rect",xmin=5.5,xmax=12.8,ymin=0.3,ymax=1.8,
           fill="#F9F9F9",color="#7F8C8D",linewidth=0.8,linetype="dashed") +
  .lbl(9.2,1.4,"w(t) = x1+x2+x3+x4",sz=4,bold=TRUE,col="#2C3E50") +
  .lbl(9.2,0.8,"Volume tumoral total (comparé aux données)",sz=3.3,col="#7F8C8D") +

  # Rétroaction croissance w → x1
  geom_curve(aes(x=9.0,y=1.8,xend=6.5,yend=3.5),
             arrow=grid::arrow(length=unit(0.28,"cm"),type="closed"),
             curvature=-0.3,linewidth=0.8,color="#27AE60") +
  .lbl(7.3,2.8,"g(w) =\nλ0·w/(1+(λ0w/λ1)^ψ)^(1/ψ)",sz=2.9,col="#1E8449") +

  # TSC
  annotate("label", x=1.8, y=3.0,
           label=paste0("Paramètres PD :\n",
                        sprintf("  λ0  = %.4f /j  (croissance exponent.)\n", L0),
                        "  λ1  = taux croissance linéaire\n",
                        sprintf("  k1  = %.3f /j  (transit, MTT=%.0f j)\n", k1, MTT),
                        sprintf("  k2  = %.2e L/µg/j  (potence drug)\n", k2),
                        sprintf("  TSC = %.0f µg/L  (= λ0/k2)", TSC)),
           size=3.2, hjust=0, fill="#F4F6F7", color="#2C3E50",
           label.padding=unit(0.4,"lines")) +

  labs(title="Structure du modèle PKPD — Simeoni (2004)  [effet linéaire k2·C]",
       subtitle="Compartiments de transit pour les cellules endommagées — TSC ≈ λ0/k2") +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold",size=14,margin=margin(b=4)),
        plot.subtitle=element_text(size=9,color="grey40",margin=margin(b=6)),
        plot.margin=margin(10,10,10,10))

SAVE("06_schema_modele_Simeoni", p6, w=13, h=7.5)

# =============================================================================
# FIG 7 — PKPD SIMEONI PRÉDIT VS OBSERVÉ + RÉSIDUS
# =============================================================================
cat("[7/9] PKPD prédit vs observé\n")

times_sim <- seq(0, max(times_d)*1.05, by=0.5)

df_sim_pd <- rbind(
  data.frame(jour=times_sim,
             w=sim_ctrl_fn(tv0_ctrl, list(L0=L0,L1=L1), times_sim),
             Groupe="Contrôle"),
  data.frame(jour=times_sim,
             w=sim_simeoni(3000,  tv0_d3,  pars_sim, times_sim),
             Groupe="3 mg/kg"),
  data.frame(jour=times_sim,
             w=sim_simeoni(10000, tv0_d10, pars_sim, times_sim),
             Groupe="10 mg/kg")
)
df_sim_pd$Groupe <- factor(df_sim_pd$Groupe, levels=lev_pd)

# Résidus relatifs aux points observés
df_res <- do.call(rbind, lapply(lev_pd, function(g) {
  obs  <- df_obs_pd[df_obs_pd$Groupe==g, ]
  pred <- sim_simeoni(
    switch(g, "Contrôle"=0, "3 mg/kg"=3000, "10 mg/kg"=10000),
    switch(g, "Contrôle"=tv0_ctrl, "3 mg/kg"=tv0_d3, "10 mg/kg"=tv0_d10),
    pars_sim, obs$jour)
  data.frame(jour=obs$jour, resid=(obs$w-pred)/pmax(pred,1)*100, Groupe=g)
}))
df_res$Groupe <- factor(df_res$Groupe, levels=lev_pd)

ymax7 <- max(df_obs_pd$w + df_obs_pd$sem, na.rm=TRUE)

p7a <- ggplot() +
  geom_vline(xintercept=DOSE_DAYS, linetype="dashed", color="grey70", linewidth=0.4) +
  geom_line(data=df_sim_pd, aes(x=jour, y=w, color=Groupe), linewidth=1.2) +
  geom_errorbar(data=df_obs_pd,
                aes(x=jour, ymin=w-sem, ymax=w+sem, color=Groupe),
                width=0.8, linewidth=0.5) +
  geom_point(data=df_obs_pd, aes(x=jour, y=w, color=Groupe, shape=Groupe), size=3) +
  annotate("point", x=DOSE_DAYS, y=-ymax7*0.07,
           shape=17, size=3.5, color="#C0392B") +
  annotate("text", x=max(DOSE_DAYS)+1.5, y=-ymax7*0.07,
           label="= Traitement", hjust=0, size=3.2, color="#C0392B") +
  scale_color_manual(values=COLS) +
  scale_shape_manual(values=c(16,17,15)) +
  scale_x_continuous(breaks=seq(0,max(times_d),by=7)) +
  coord_cartesian(ylim=c(-ymax7*0.12, ymax7*1.1), clip="off") +
  labs(title="Modèle PK/PD TGI — Simeoni (2004) — Prédit vs Observé",
       subtitle=paste0(
         "λ0=", round(L0,4), " /j  |  ",
         "k2=", formatC(k2,format="e",digits=2), " L/µg/j  |  ",
         "k1=", round(k1,3), " /j  |  ",
         "TSC=", round(TSC,0), " µg/L  |  ",
         "MTT=", round(MTT,1), " j"),
       x="Temps (jours)", y="Volume tumoral (mm³)",
       color=NULL, shape=NULL,
       caption="Lignes = modèle   ●▲■ = données observées (mean ± SEM)") +
  THEME + theme(legend.position="bottom", plot.margin=margin(5,10,30,5))

p7b <- ggplot(df_res, aes(x=jour, y=resid, color=Groupe)) +
  geom_hline(yintercept=0, linewidth=0.8, color="grey40") +
  geom_hline(yintercept=c(-30,30), linetype="dashed", color="grey70") +
  geom_point(size=2.5, alpha=0.9) +
  scale_color_manual(values=COLS) +
  labs(title="Résidus relatifs",
       subtitle="(Obs − Pred)/Pred × 100",
       x="Temps (jours)", y="Résidu (%)", color=NULL) +
  THEME + theme(legend.position="none")

png(file.path(OUT,"07_PKPD_predit_vs_observe.png"),
    width=14, height=6, units="in", res=180, bg="white")
grid.arrange(p7a+theme(legend.position="bottom"),
             p7b, ncol=2, widths=c(2,1))
dev.off()
cat("  →", file.path(OUT,"07_PKPD_predit_vs_observe.png"),"\n")

# =============================================================================
# FIG 8 — TABLEAU TGI
# =============================================================================
cat("[8/9] Tableau TGI\n")

tgi_f <- function(wt, wt0, wc, wc0) {
  dc <- wc - wc0
  if (is.na(dc) || dc <= 0) return(NA_real_)
  round((1 - (wt-wt0)/dc)*100, 1)
}

t_fine   <- seq(0, max(times_d), by=0.2)
wc_sim   <- sim_ctrl_fn(tv0_ctrl, list(L0=L0,L1=L1), t_fine)
w3_sim   <- sim_simeoni(3000,  tv0_d3,  pars_sim, t_fine)
w10_sim  <- sim_simeoni(10000, tv0_d10, pars_sim, t_fine)

ok_tgi   <- !is.na(w_ctrl) & !is.na(w_d10_tv) & !is.na(w_d3_tv) & times_d > 0
t_match  <- times_d[ok_tgi]

wc_s  <- approx(t_fine, wc_sim,  xout=t_match)$y
w3_s  <- approx(t_fine, w3_sim,  xout=t_match)$y
w10_s <- approx(t_fine, w10_sim, xout=t_match)$y

df_tgi <- data.frame(
  Temps       = paste0("j", t_match),
  TGI_3_obs   = mapply(tgi_f, w_d3_tv[ok_tgi],  tv0_d3,  w_ctrl[ok_tgi], tv0_ctrl),
  TGI_3_pred  = mapply(tgi_f, w3_s,              tv0_d3,  wc_s,           tv0_ctrl),
  TGI_10_obs  = mapply(tgi_f, w_d10_tv[ok_tgi], tv0_d10, w_ctrl[ok_tgi], tv0_ctrl),
  TGI_10_pred = mapply(tgi_f, w10_s,             tv0_d10, wc_s,           tv0_ctrl)
)

tgi_hdr <- c("Temps","TGI 3mg obs (%)","TGI 3mg pred (%)","TGI 10mg obs (%)","TGI 10mg pred (%)")
nc8 <- length(tgi_hdr); nr8 <- nrow(df_tgi)

.fill_tgi <- function(v, bg) {
  if (is.na(v)) return(bg)
  if (v >= 60)  return("#A9DFBF")
  if (v >= 30)  return("#FAD7A0")
  bg
}

rows8 <- lapply(seq_len(nr8), function(i) {
  bg  <- if (i%%2==0) "#F2F3F4" else "white"
  r   <- df_tgi[i,]
  lbs <- c(r$Temps,
            ifelse(is.na(r$TGI_3_obs),  "—", sprintf("%.1f %%", r$TGI_3_obs)),
            ifelse(is.na(r$TGI_3_pred), "—", sprintf("%.1f %%", r$TGI_3_pred)),
            ifelse(is.na(r$TGI_10_obs), "—", sprintf("%.1f %%", r$TGI_10_obs)),
            ifelse(is.na(r$TGI_10_pred),"—", sprintf("%.1f %%", r$TGI_10_pred)))
  fls <- c(bg,
            .fill_tgi(r$TGI_3_obs, bg),  .fill_tgi(r$TGI_3_pred, bg),
            .fill_tgi(r$TGI_10_obs, bg), .fill_tgi(r$TGI_10_pred, bg))
  data.frame(row=i, col=1:nc8, lbl=lbs, fill=fls, bold=(i==nr8))
})
hdr8 <- data.frame(row=0, col=1:nc8, lbl=tgi_hdr, fill="#1A5276", bold=TRUE)
cells8 <- rbind(hdr8, do.call(rbind, rows8))
cells8$y    <- max(cells8$row) - cells8$row
cells8$tcol <- ifelse(cells8$row==0, "white", "black")
cells8$face <- ifelse(cells8$bold, "bold", "plain")

p8 <- ggplot(cells8, aes(x=col, y=y)) +
  geom_tile(aes(fill=fill), color="grey70", linewidth=0.3) +
  geom_text(aes(label=lbl, color=tcol, fontface=face), size=3.2) +
  scale_fill_identity() + scale_color_identity() +
  labs(title="Tumour Growth Inhibition — Observé vs Prédit (Simeoni)",
       subtitle="Vert : TGI ≥ 60%  ·  Orange : 30–60%  ·  TGI = (1 − ΔW_traité/ΔW_ctrl) × 100") +
  theme_void(base_size=11) +
  theme(plot.title=element_text(face="bold",size=13,margin=margin(b=4)),
        plot.subtitle=element_text(size=9,color="grey40",margin=margin(b=6)),
        plot.margin=margin(12,12,12,12))

SAVE("08_tableau_TGI", p8, w=11, h=0.5*(nr8+2)+2)

# =============================================================================
# FIG 9 — TABLEAU PARAMÈTRES PKPD
# =============================================================================
cat("[9/9] Tableau paramètres PKPD\n")

pkpd_tbl <- data.frame(
  Param = c("λ0",  "λ1",  "ψ",    "k1",  "k2",
            "TSC", "MTT",
            "V1",  "V2",  "CL",   "Q",  "tv0"),
  Type  = c("PD","PD","PD (fixé)","PD","PD",
            "Dérivé","Dérivé",
            "PK (fixé)","PK (fixé)","PK (fixé)","PK (fixé)","Initial"),
  Desc  = c("Taux de croissance exponentielle",
            "Taux de croissance linéaire",
            "Coeff. Hill (Simeoni 2004)",
            "Constante de transit cellulaire",
            "Potence cytotoxique (effet linéaire)",
            "Conc. statique tumorale  =  λ0/k2",
            "Mean Transit Time  =  4/k1",
            "Volume compartiment central",
            "Volume périphérique",
            "Clairance systémique",
            "Clairance intercompartimentale",
            "Volume tumoral initial"),
  Valeur = c(
    sprintf("%.4f /j",         L0),
    "non identifiable sur ces données",
    "20",
    sprintf("%.3f /j   (MTT = %.0f j)", k1, MTT),
    sprintf("%.3e L/µg/j", k2),
    sprintf("%.0f µg/L", TSC),
    sprintf("%.0f jours", MTT),
    sprintf("%.4f L/kg", PK["V1"]),
    sprintf("%.4f L/kg", PK["V2"]),
    sprintf("%.4f L/j/kg  (%.4f L/h/kg)", PK["CL"], PK_h["CL"]),
    sprintf("%.4f L/j/kg", PK["Q"]),
    sprintf("%.0f mm³", mean(c(tv0_ctrl, tv0_d3, tv0_d10)))
  ),
  Méthode = c("DEoptim + nlminb","DEoptim + nlminb","Littérature",
              "DEoptim + nlminb (borne = 0.2 /j)","DEoptim + nlminb",
              "Calculé","Calculé",
              "nlminb (PK séparé)","nlminb (PK séparé)",
              "nlminb (PK séparé)","nlminb (PK séparé)",
              "Moyenne j0 tous groupes")
)

type_fill <- c("PD"="#D5F5E3","PD (fixé)"="#FDFEFE",
               "Dérivé"="#FEF9E7","PK (fixé)"="#D6EAF8","Initial"="#F4ECF7")

nr9 <- nrow(pkpd_tbl); nc9 <- 5
rows9 <- lapply(seq_len(nr9), function(i) {
  bg <- type_fill[pkpd_tbl$Type[i]]
  data.frame(row=i, col=1:nc9,
             lbl=c(pkpd_tbl$Param[i], pkpd_tbl$Type[i], pkpd_tbl$Desc[i],
                   pkpd_tbl$Valeur[i], pkpd_tbl$Méthode[i]),
             fill=bg, bold=FALSE)
})
hdr9 <- data.frame(row=0, col=1:nc9,
                   lbl=c("Param.","Type","Description","Valeur","Méthode"),
                   fill="#1A5276", bold=TRUE)
cells9 <- rbind(hdr9, do.call(rbind, rows9))
cells9$y    <- max(cells9$row) - cells9$row
cells9$tcol <- ifelse(cells9$row==0, "white", "black")
cells9$face <- ifelse(cells9$bold, "bold", "plain")

p9 <- ggplot(cells9, aes(x=col, y=y)) +
  geom_tile(aes(fill=fill), color="grey70", linewidth=0.3) +
  geom_text(aes(label=lbl, color=tcol, fontface=face), size=3.1) +
  scale_fill_identity() + scale_color_identity() +
  labs(title="Récapitulatif des paramètres — Modèle PK/PD TGI Simeoni (2004)",
       subtitle="Vert=PD estimés  ·  Bleu=PK fixés  ·  Jaune=dérivés  ·  Violet=conditions initiales") +
  theme_void(base_size=11) +
  theme(plot.title=element_text(face="bold",size=13,margin=margin(b=4)),
        plot.subtitle=element_text(size=9,color="grey40",margin=margin(b=6)),
        plot.margin=margin(12,12,12,12))

SAVE("09_parametres_PKPD", p9, w=13, h=5.5)

# =============================================================================
cat("\n=== Toutes les figures générées dans", OUT, "===\n")
print(list.files(OUT, pattern="\\.png$"))
