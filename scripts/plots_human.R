############################################################
# plots_human_neut_plt.R
# Visualisation humaine -- Neutrophils + Platelets uniquement
############################################################
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(readr)
library(dplyr)

# -- Panneau individuel -----------------------------------
plot_cell_panel <- function(sim, yvar, title, color,
                            baseline  = NULL,
                            dose_days = NULL,
                            obs       = NULL) {
  
  y_vals     <- sim[[yvar]]
  y_vals_pos <- y_vals[is.finite(y_vals) & y_vals > 0]
  
  all_pos <- y_vals_pos
  if (!is.null(obs) && length(obs$value) > 0)
    all_pos <- c(all_pos, obs$value[obs$value > 0])
  
  ymin <- 10^floor(log10(min(all_pos) * 0.3))
  ymax <- 10^ceiling(log10(max(all_pos) * 3.0))
  breaks_y <- 10^seq(log10(ymin), log10(ymax), by = 1)
  
  df_sim <- data.frame(
    x = sim$days,
    y = pmax(sim[[yvar]], ymin * 0.5)
  )
  
  p <- ggplot() +
    geom_line(data = df_sim,
              aes(x = x, y = y),
              color = color, linewidth = 1.0, alpha = 0.95)
  
  if (!is.null(obs) && length(obs$time) > 0) {
    df_obs <- data.frame(
      x = obs$time,
      y = pmax(obs$value, ymin * 0.5)
    )
    p <- p + geom_point(data   = df_obs,
                        aes(x = x, y = y),
                        shape  = 21,
                        fill   = "white",
                        color  = "black",
                        size   = 1.6,
                        stroke = 0.6,
                        alpha  = 0.75)
  }
  
  if (!is.null(baseline) && baseline > 0)
    p <- p + geom_hline(yintercept = baseline, linetype = "dashed",
                        color = "grey50", linewidth = 0.5, alpha = 0.6)
  
  if (!is.null(dose_days)) {
    vd <- dose_days[dose_days <= max(sim$days)]
    if (length(vd) > 0)
      p <- p + geom_vline(xintercept = vd, linetype = "dotted",
                          color = "grey55", linewidth = 0.4, alpha = 0.55)
  }
  
  p +
    scale_y_log10(limits = c(ymin, ymax), breaks = breaks_y,
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = title, x = "Time (d)",
         y = expression(10^9~cells~L^{-1})) +
    theme_bw(base_size = 9.5) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92"),
          plot.title  = element_text(face = "bold", size = 9, hjust = 0.5),
          axis.title  = element_text(size = 7.5),
          axis.text   = element_text(size = 7))
}

# -- Grille 1×2 : Neut + Plt uniquement ------------------
plot_human_neut_plt <- function(sim, pars,
                                titre     = "Hematopoiesis - Human",
                                dose_days = NULL,
                                obs_list  = NULL) {
  
  get_obs <- function(name) {
    if (!is.null(obs_list) && name %in% names(obs_list))
      obs_list[[name]]
    else NULL
  }
  
  panels <- list(
    plot_cell_panel(sim, "Neut", "Neutrophils", "#d7191c",
                    pars$Neut0, dose_days, get_obs("Neut")),
    plot_cell_panel(sim, "Plt",  "Platelets",   "#1b9e77",
                    pars$Plt0,  dose_days, get_obs("Plt"))
  )
  
  grid.arrange(grobs = panels, ncol = 2,
               top = textGrob(titre, gp = gpar(fontface = "bold", fontsize = 11)))
}

# -- Lecture des données observées ------------------------
read_wpd <- function(path) {
  df <- read_csv(path, col_names = c("time", "value"),
                 skip = 1, show_col_types = FALSE)
  df <- df %>% arrange(time) %>% distinct()
  list(time = df$time, value = df$value)
}

data_dir <- if (file.exists("Neut_H.csv")) "." else
            if (file.exists("../data/Neut_H.csv")) "../data" else "data"

obs_list_human <- list(
  Neut = read_wpd(file.path(data_dir, "Neut_H.csv")),
  Plt  = read_wpd(file.path(data_dir, "Plt_H.csv"))
)

# -- Appel principal --------------------------------------
# Doses humaines : 2 cycles Q21D (jours 0 et 21)
if (exists("sim_hu")) {
  plot_human_neut_plt(
    sim       = sim_hu,
    pars      = init_pars,
    titre     = "Hematopoiesis - Carboplatin (Human) Q21D × 2",
    dose_days = c(0, 21),
    obs_list  = obs_list_human
  )
}

# -- Sauvegarde PDF ---------------------------------------
save_human_neut_plt <- function(sim, pars, file,
                                titre     = "Hematopoiesis - Human",
                                dose_days = NULL,
                                obs_list  = NULL,
                                width = 10, height = 5) {
  pdf(file, width = width, height = height)
  plot_human_neut_plt(sim, pars, titre = titre,
                      dose_days = dose_days, obs_list = obs_list)
  dev.off()
  message("✓ Sauvegardé : ", file)
}
