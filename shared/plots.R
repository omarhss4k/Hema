############################################################
# plots.R
# Fonctions de visualisation -- axes adaptatifs par panneau
# Avec superposition des données observées (Fornari 2019)
############################################################
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# -- Panneau individuel -----------------------------------
plot_cell_panel <- function(sim, yvar, title, color,
                            baseline  = NULL,
                            dose_days = NULL,
                            obs       = NULL) {

  y_vals     <- sim[[yvar]]
  y_vals_pos <- y_vals[is.finite(y_vals) & y_vals > 0]

  # Limites dynamiques (inclut les obs si présentes)
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

  # Points observés (cercles ouverts noirs)
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

  # Baseline (tiretée grise)
  if (!is.null(baseline) && baseline > 0)
    p <- p + geom_hline(yintercept = baseline, linetype = "dashed",
                        color = "grey50", linewidth = 0.5, alpha = 0.6)

  # Marqueurs doses (pointillés verticaux)
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
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92"),
          plot.title  = element_text(face = "bold", size = 13, hjust = 0.5),
          axis.title  = element_text(size = 11),
          axis.text   = element_text(size = 10))
}

# -- Grille 4×2 -------------------------------------------
plot_all_cells <- function(sim, pars,
                           titre     = "Hematopoiesis",
                           dose_days = NULL,
                           obs_list  = NULL) {

  get_obs <- function(name) {
    if (!is.null(obs_list) && name %in% names(obs_list))
      obs_list[[name]]
    else NULL
  }

  panels <- list(
    plot_cell_panel(sim, "MPP",  "Multi-potent progenitors",   "#2166ac",
                    pars$MPP0,  dose_days, get_obs("MPP")),
    plot_cell_panel(sim, "Neut", "Neutrophils",                "#d7191c",
                    pars$Neut0, dose_days, get_obs("Neut")),
    plot_cell_panel(sim, "CMP",  "Common myeloid progenitors", "#4dac26",
                    pars$CMP0,  dose_days, get_obs("CMP")),
    plot_cell_panel(sim, "Mono", "Monocytes",                  "#a6611a",
                    pars$Mono0, dose_days, get_obs("Mono")),
    plot_cell_panel(sim, "MEP",  "MEP",                        "#7b3294",
                    pars$MEP0,  dose_days, get_obs("MEP")),
    plot_cell_panel(sim, "Plt",  "Platelets",                  "#1b9e77",
                    pars$Plt0,  dose_days, get_obs("Plt")),
    plot_cell_panel(sim, "Ret",  "Reticulocytes",              "#e377c2",
                    pars$Ret0,  dose_days, get_obs("Ret")),
    plot_cell_panel(sim, "RBC",  "Red blood cells",            "#cc3399",
                    pars$RBC0,  dose_days, get_obs("RBC"))
  )

  arrangeGrob(grobs = panels, ncol = 2,
              top = textGrob(titre, gp = gpar(fontface = "bold", fontsize = 15)))
}

# -- Sauvegarder en PDF + PNG -----------------------------
save_all_cells <- function(sim, pars, file, titre,
                           dose_days = NULL, obs_list = NULL,
                           width = 12, height = 11) {
  p <- plot_all_cells(sim, pars, titre = titre,
                      dose_days = dose_days, obs_list = obs_list)

  pdf(file, width = width, height = height)
  grid::grid.draw(p)
  dev.off()
  message("✓ Sauvegardé : ", file)

  png_file <- sub("\\.pdf$", ".png", file)
  png(png_file, width = width, height = height, units = "in", res = 300)
  grid::grid.draw(p)
  dev.off()
  message("✓ Sauvegardé : ", png_file)
}

# -- Résumé nadir/peak ------------------------------------
print_summary <- function(sim, pars) {
  cells <- list(c("Neut","Neut0"), c("Mono","Mono0"),
                c("Plt","Plt0"),   c("Ret","Ret0"),
                c("RBC","RBC0"),   c("MPP","MPP0"),
                c("CMP","CMP0"),   c("MEP","MEP0"))
  cat(sprintf("\n%-6s | %-10s | %-10s | %7s | %7s\n",
              "Cell","Baseline","Nadir","Nadir%","Peak%"))
  cat(paste(rep("-",50), collapse=""), "\n")
  for (cell in cells) {
    name <- cell[1]; base <- pars[[cell[2]]]; y <- sim[[name]]
    cat(sprintf("%-6s | %10.2f | %10.4f | %6.1f%% | %6.1f%%\n",
                name, base, min(y), min(y)/base*100, max(y)/base*100))
  }
  cat("\n")
}
