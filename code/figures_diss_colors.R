install.packages(c('latex2exp', 'scales', 'zoo', 'stable', 'rmutil', 'ggridges', 
                   'ggnewscale', 'cowplot', 'reshape2', 'scico', 'forcats', 
                   'stringr', 'dplyr', 'purrr', 'readr', 'tidyr', 'tibble', 'ggplot2', 'tidyverse'), 
                 version = c('0.9.6', '1.3.0', '1.8-11', '1.1.6', '1.1.10', '0.5.4', 
                             '0.4.8',  '1.1.1', '1.4.4', '1.5.0', '0.5.2', 
                             '1.5.0', '1.0.10', '1.0.1', '2.1.3', '1.2.1', '3.1.8', '3.5.1', '1.3.2'), 
                 repos='https://cran.rstudio.com/')

if (require(devtools)) {
  install_github('shabbychef/ggallin')
}

library(tidyverse)
library(scico)
library(reshape2)
library(stringr)
library(cowplot)
library(ggallin)
library(ggnewscale)
library(ggridges)
library(stable)
library(zoo)
library(scales)
library(latex2exp)

setwd("~/Desktop/PhD/stableEWS_paper")

theme_set(
  theme_classic() + 
    theme(
      axis.text = element_text(color = "black", size = 15),
      axis.title = element_text(color = "black", size = 15),
      plot.title = element_text(color = "black", size = 15),
      plot.subtitle = element_text(color = "black", size = 15),
      plot.caption = element_text(color = "black", size = 15),
      strip.text = element_text(color = "black", size = 15),
      legend.text = element_text(color = "black", size = 15),
      legend.title = element_text(color = "black", size = 15),
      axis.line = element_line(color = "black"),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
      legend.background = element_rect(fill='transparent', color = NA),
      legend.box.background = element_rect(fill='transparent', color = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),  
      plot.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "transparent", color = NA)
    )
)

#overview

colors = c("2" = "black", "1.8" =  "#0072B2", "1.3" = "#009E73", "1" =  "#56B4E9")
highlight = "#D55E00"

alpha_stable_dist = function() {
  x = seq(-5, 5, by = .1)
  
  df = data.frame("x" = x,
                  "2" = dstable(x, tail = 2),
                  "1.8" = dstable(x, tail = 1.8),
                  "1.5" = dstable(x, tail = 1.5),
                  "1.3" = dstable(x, tail = 1.3)) %>%
    melt(id.vars = "x") %>%
    mutate(variable = str_remove(variable, "^X"))
  

  colors = c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9")
  
  df$variable = factor(df$variable, levels = c("2", "1.8",  "1.5", "1.3"))
  
  (p2 = ggplot() + theme_bw() +
      geom_line(data = df, 
                aes(x = x, y = value, color = variable, linewidth = variable),
                linewidth = 0.75, alpha = .9) +
      scale_y_continuous(name = "p(X)", expand = c(0,0), limits = c(-0.001,.43)) +
      scale_x_continuous(name = "X", expand = c(0,0)) +
      theme(text = element_text(size = 12), 
            legend.position = c(.2, .75),
            legend.background = element_rect(fill='white', color = "darkgrey"),
            legend.box.background = element_rect(fill='transparent', color = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),  
            plot.background = element_rect(fill = "transparent", colour = NA),
            strip.background = element_rect(fill = "transparent", color = NA),
            plot.margin = margin(t = 5.5, l = 5.5, b = 5.5, r = 15, unit = "pt")) +
      scale_color_manual(values = colors, name = expression(alpha)))
  
  return(p2)
}

(p1 = alpha_stable_dist())

plot_var_eq = function() {
  
  df = read.csv(paste0("data/var_lin_eq_alpha2.csv")) %>%
    filter(k != 0) %>%
    mutate(system = "Linear")
  
  df_nol = read.csv(paste0("data/var_nol_eq_alpha2.csv")) %>%
    filter(k != 0) %>%
    mutate(system = "Non-linear \n(Fold)")
  
  df = df_nol %>%
    full_join(df) %>%
    mutate(var_true = 0.02/(2*k))
  
  df$system = factor(df$system, levels = c("Linear", "Non-linear \n(Fold)"))
  
  colors = c("#000000", "#56B4E9")
  colors_means = c("#56B4E9", "#E69F00")
  line_color = c("True variance" = "grey20")
  sizes = c(2.5, 2)
  names(colors) <- levels(df$system)
  names(colors_means) <- levels(df$system)
  names(sizes) <- levels(df$system)
  
  df_mean = df %>%
    group_by(k, system) %>%
    summarize(mean = mean(variance))
  
  df_mean$system = factor(df_mean$system, levels = c("Non-linear \n(Fold)", "Linear"))
  
  (p = ggplot() + theme_bw() +
      geom_line(data = df, aes(x = k, y = var_true, linetype = "Theory"), linewidth =.5, alpha = .7) +
      geom_point(data = df, aes(x = k, y = variance, color = system, shape = system), 
                stroke = 0.2,  alpha = .5, position = position_jitter(width = .25), size = 0.75) +
      geom_point(data = df_mean, aes(x = k, y = mean, color = system, shape = system), stroke =  1) +
      scale_linetype_manual(values = c("Theory" = "dashed"), name = "") +
      scale_size_manual(values = sizes, name = "System")  +
      scale_shape_manual(values = c(3, 4), name = "System") +
      scale_color_manual(values = colors_means, name = "System") +
      scale_y_continuous(trans = 'log10', name = "EWS Var X ", expand = c(0,0)) +
      scale_x_continuous(trans = 'log10', name = "Bifurcation parameter k", breaks = c(0.01, 1, 100), expand = c(0,0), labels = c(0.01, 1, 100)) +
      theme(text = element_text(size = 12),
            legend.position = c(.75, .75),
            legend.title = element_blank(),
            legend.background = element_rect(fill='white', color = "darkgrey"),
            legend.box.background = element_rect(fill='transparent', color = NA),
            plot.margin = margin(t = 5.5, l = 5.5, b = 5.5, r = 5.5, unit = "pt")) + 
      guides(size = "none"))
  
  return(p)
  
}

(p2 = plot_var_eq())


plot_grid(p1 , p2, nrow = 1, labels = c("(a)", "(b)"),  hjust = .05)
ggsave("reports/paper/theory_diss.pdf", width = 8.2, height = 4,  bg = "white", dpi = 600)

#gamma equilibrium runs
plot_gamma_eq = function(A =c(1.3, 1.5, 1.8, 2)) {
  
  data = list()
  
  for (a in A) {
    df = read.csv(paste0("data/gamma_lin_eq_alpha", a, ".csv")) %>%
      filter(k != 0) %>%
      mutate(system = "Linear",
             alpha = paste0(a))
    
    
    df_nol = read.csv(paste0("data/gamma_nol_eq_alpha", a, ".csv")) %>%
      filter(k != 0) %>%
      mutate(system = "Non-linear \n(Fold)",
             alpha = paste0(a))
    
    df = df_nol %>%
      full_join(df) %>%
      filter(gamma_sample < 10^40) %>%
      mutate(gamma_true = 0.1*(1/(a*k))^(1/a))
    
    data = append(data, list(df))
  }
  
  df = purrr::reduce(data, bind_rows)
  
  
  df$system = factor(df$system, levels = c("Linear", "Non-linear \n(Fold)"))
  
  colors = c("#000000", "#98C6EA")
  colors_means = c("#56B4E9", "#E69F00")
  line_color = c("Theory" = "#D55E00")
  sizes = c(2.5, 2)
  names(colors) <- levels(df$system)
  names(colors_means) <- levels(df$system)
  names(sizes) <- levels(df$system)
  
  df_mean = df %>%
    group_by(alpha, k, system) %>%
    summarize(mean = mean(gamma_sample))
  
  df_mean$system = factor(df_mean$system, levels = c("Non-linear \n(Fold)", "Linear"))
  
  df$alpha <- factor(df$alpha,
                     levels = rev(c(2, 1.8, 1.5, 1.3)),
                     labels=rev(c('2'=parse(text=TeX('$\\alpha$ = 2')),
                                  '1.8'=parse(text=TeX('$\\alpha$ = 1.8')),
                                  '1.5'=parse(text=TeX('$\\alpha$ = 1.5')),
                                  '1.3'=parse(text=TeX('$\\alpha$ = 1.3')))))
  
  df_mean$alpha <- factor(df_mean$alpha,
                          levels = rev(c(2, 1.8, 1.5, 1.3)),
                          labels=rev(c('2'=parse(text=TeX('$\\alpha$ = 2')),
                                       '1.8'=parse(text=TeX('$\\alpha$ = 1.8')),
                                       '1.5'=parse(text=TeX('$\\alpha$ = 1.5')),
                                       '1.3'=parse(text=TeX('$\\alpha$ = 1.3')))))
  
  
  (p = ggplot() + 
      geom_vline(xintercept = 0, color = "black", linewidth = .75) +
      geom_line(data = df, aes(x = k, y = gamma_true, color = "Theory", linetype = "Theory", shape = "Theory"), linewidth = .75) +
      facet_wrap(~rev(alpha), labeller = label_parsed, nrow = 1) + 
      geom_point(data = df, aes(x = k, y = gamma_sample, color = system, shape = system, linetype = system), 
                 stroke = 0.2, position = position_jitter(width = .25), size = 0.25, alpha = .4) +
      geom_point(data = df_mean, aes(x = k, y = mean, color = system, shape = system, size = system, linetype = system), 
                 stroke = 1) +
      scale_size_manual(values = c("Linear" = 3, "Non-linear \n(Fold)" = 3, "Theory" = 3), name = "System") +
      scale_shape_manual(values = c("Linear" = 3, "Non-linear \n(Fold)" = 4, "Theory" = 21), name = "System") +
      scale_color_manual(values = c("Linear" = "#56B4E9", "Non-linear \n(Fold)" = "#E69F00", "Theory" = "#D55E00"), name = "System") +
      scale_linetype_manual(values = c("Linear" = "solid", "Non-linear \n(Fold)" = "solid", "Theory" = "dashed"), name = "System") +
      scale_y_continuous(trans = 'log10', name = expression(EWS~gamma[X]), expand = c(0,0)) +
      scale_x_continuous(trans = 'log10', name = "Bifurcation parameter k", 
                         breaks = c(0.01, 1, 100), labels = c("0.01", "1", "100"),expand = c(0,0)) +
      theme(axis.text.x = element_text(hjust = 0.75, size = 15),
            legend.position = "bottom",
            legend.direction = "horizontal") + 
      guides(size = "none"))
  
  return(p)
  
}
plot_gamma_eq()
ggsave("reports/paper/gamma_equilibrium_diss.pdf", width = 9.5, height = 4, scale = 1)


#gamma non-equilibrium runs
plot_gamma_neq = function(A =c(1.3, 1.5, 1.8, 2), cutoff_theory = .15) {
  
  data = list()
  
  for (a in A) {
   df = read.csv(paste0("data/gamma_lin_neq_alpha", a, ".csv")) %>%
      filter(k != 0) %>%
      mutate(system = "Linear",
             alpha = paste0(a)) %>%
      filter(gamma_sample != 0)
    
    df_nol = read.csv(paste0("data/gamma_nol_neq_alpha", a, ".csv")) %>%
      filter(k != 0) %>%
      mutate(system = "Non-linear \n(Fold)",
             alpha = paste0(a)) %>%
      filter(gamma_sample != 0) 
    
    df = df_nol %>%
      full_join(df) %>%
      mutate(gamma_true = 0.1*(1/(a*k))^(1/a)) 
    
    data = append(data, list(df))
  }
  
  df = purrr::reduce(data, bind_rows) %>%
    filter(timestep > 1000) %>%
    group_by(alpha, system, sample) %>%
    arrange(alpha, system, sample, k) %>%
    mutate(gamma_smoothed = slider::slide_dbl(gamma_sample, mean, .before =50, .after = 50)) %>%
    filter(gamma_true < cutoff_theory,
           gamma_smoothed < cutoff_theory)
  
  df_mean = df %>%
    group_by(alpha, k, system) %>%
    summarize(mean_gamma = mean(gamma_smoothed, na.rm = TRUE)) %>%
    group_by(system) %>%
    mutate(mean_max = max(mean_gamma),
           mean_min = min(mean_gamma)) %>%
    ungroup()
  
  df$alpha_label <- factor(df$alpha,
                     levels = rev(c(2, 1.8, 1.5, 1.3)),
                     labels=rev(c('2'=parse(text=TeX('$\\alpha$ = 2')),
                              '1.8'=parse(text=TeX('$\\alpha$ = 1.8')),
                              '1.5'=parse(text=TeX('$\\alpha$ = 1.5')),
                              '1.3'=parse(text=TeX('$\\alpha$ = 1.3')))))
  
  df_mean$alpha_label <- factor(df_mean$alpha,
                          levels = rev(c(2, 1.8, 1.5, 1.3)),
                          labels=rev(c('2'=parse(text=TeX('$\\alpha$ = 2')),
                                   '1.8'=parse(text=TeX('$\\alpha$ = 1.8')),
                                   '1.5'=parse(text=TeX('$\\alpha$ = 1.5')),
                                   '1.3'=parse(text=TeX('$\\alpha$ = 1.3')))))
  
  legend_labels = (c(TeX("$\\alpha$ = 2"), TeX("$\\alpha$ = 1.8"), TeX("$\\alpha$ = 1.5"), TeX("$\\alpha$ = 1.3"), "Theory"))
  
  colors = rev(c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9", "Theory" = "#D55E00"))
  
  (p = ggplot() + 
    geom_vline(xintercept = 0, color = "black", linewidth = .75) +
    geom_hline(yintercept = 0.025, color = "black", linewidth = 0.75) + 
    geom_line(data = df, aes(x = k, y = gamma_smoothed, group = interaction(run, alpha), color = as.factor(alpha), linetype = "Simulation"), 
              alpha = .75, linewidth = .001) +
    geom_line(data = df[df$system == "Linear",], aes(x = k, y = gamma_true, group = interaction(system, alpha), 
                                                     color = "Theory", linetype = "Theory"), linewidth = .6) +
    geom_line(data = df_mean, aes(x = k, y = mean_gamma, color = as.factor(alpha), linetype = "Simulation"), 
              alpha = .9, linewidth = .7) +
    scale_color_manual(values = colors,  name = expression(~gamma[X])) + 
    scale_linetype_manual(values = c("Theory" = "dashed", "Simulation" = "solid"), name = expression(~gamma[X])) + 
    scale_x_continuous(limits = c(0.25, 5), expand = c(0,0), breaks = c(1, 3, 5), name = expression(symbol('\254')~~~Bifurcation~parameter~k)) +
    scale_y_continuous(name = expression(gamma[X]), limits = c(0.025, 0.09), expand = c(0,0), breaks = c(0, 0.05, 0.075)) +
    facet_grid(rows = vars(system), cols = vars(alpha_label), scales = "free", labeller = label_parsed) + 
    theme(legend.position = "bottom",
            legend.direction = "horizontal") +
      guides(color = guide_legend(ncol = 4, byrow = TRUE),
             linetype = "none"))
  
  return(p)
}
plot_gamma_neq() 

ggsave("reports/paper/gamma_noneq_diss.pdf",  width = 10, height = 6, scale = 1)


#gamma trajectories
extract_numbers = function(string) {
  pattern <- "i(\\d+)_k([0-9.]+)"
  
  num1 = integer(length(string))
  num2 = numeric(length(string))
  
  for (i in seq_along(string)) {
    matches = str_match(string[i], pattern)
    if (!is.na(matches[1, 2]) && !is.na(matches[1, 3])) {
      num1[i] = as.integer(matches[1, 2])
      num2[i] = as.numeric(matches[1, 3])
    } else {
      num1[i] = NA
      num2[i] = NA
    }
  }
  
  return(data.frame("run" = num1, "k" =  num2))
} #helper function to extract k and run from column identifiers 
plot_traj_lin = function(system.type = "lin", A = c(2, 1.5, 1)) {
  data = list()
  for (a in A) {
    traj = read.csv(paste0("data/trajectories_", system.type, "_neq_", a, ".csv")) %>%
      select(-template) %>%
      mutate(timestep = row_number()) %>%
      melt(id.vars = c("timestep", "K")) %>%
      mutate(run = extract_numbers(variable)$run,
             "alpha" = a, x_star = 0) %>%
      filter(value != 0,
             abs(value) < 1.5,
             K < 1,
             run == 1)
    
    
    
    data = append(data, list(traj))
  }
  
  df = purrr::reduce(data, bind_rows)
  
  df$alpha = factor(df$alpha, levels = c(2, 1.5, 1))
  df$alpha = with(df, factor(alpha, levels = rev(levels(alpha))))
  colors = c("2" = "black", "1.5" = "#009E73", "1" =  "#56B4E9", "Fixed points X*" = "#D55E00")
  line_color = c("Fixed points X*" = "#D55E00")
  
  fixed_points = data.frame(x = c(0, 0, 0),
                            k = c(-Inf, 1, 2),
                            type = c("unstable", "stable", "ghost"))
  
  (p = ggplot() + 
      coord_cartesian(clip = "off") +
      scale_x_continuous(name = expression(symbol('\254')~~~Bifurcation~parameter~k), limits = c(-0.05, 1), expand = c(0, 0)) +
      scale_y_continuous(name = "State X", trans=pseudolog10_trans, expand = c(0,0)) + 
      geom_hline(yintercept = 0, color = "darkgrey") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey20") +
      geom_segment(aes(x = -0.05, y = 0, yend = 0, xend = 0.1), color = "#D55E00", linetype = "dashed", linewidth = .75) + 
      geom_line(data = df, aes(x = K, y = x_star, color = "Fixed points X*"), linewidth = .75)  +
      geom_point(data = fixed_points, aes(x = k, y = x, shape = type), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1) +
      geom_line(data = df, aes (x = K, y = value, color = alpha, group = interaction(alpha, run)), alpha = .9, linewidth = .5) +
      scale_color_manual(values = colors, name = expression(alpha)) +
      scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            plot.margin = unit(c(0, 0.75, 0, 0.5), "cm")) +
      guides(color = guide_legend(override.aes = list(linewidth = 1))))
  
  return(p)
}
plot_traj_nol = function(system.type = "nol", A = c(2, 1.5, 1)) {
  
  data = list()
  for (a in A) {
    traj = read.csv(paste0("data/trajectories_", system.type, "_neq_", a, ".csv")) %>%
      select(-template) %>%
      mutate(timestep = row_number()) %>%
      melt(id.vars = c("timestep", "K")) %>%
      mutate(run = extract_numbers(variable)$run,
             "alpha" = a, 
             x_star = sqrt(K),
             unstable_x_star = -x_star) %>%
      filter(value != 0,
             #value > -1,
             K < 1,
             run == 1)
    
    
    
    data = append(data, list(traj))
  }
  
  df = purrr::reduce(data, bind_rows)
  
  df$alpha = factor(df$alpha, levels = c(2, 1.5, 1, 0.5))
  df$alpha = with(df, factor(alpha, levels = rev(levels(alpha))))
  colors = c("2" = "black", "1.5" =  "#009E73", "1" =  "#56B4E9")
  
  fixed_points = data.frame(x = c(-1, 0, 1),
                            k = c(1, 0,  1),
                            type = c("unstable", "ghost", "stable"))
  
  
  (p = ggplot() +
      coord_cartesian(clip = "off") +
      scale_x_continuous(name = expression(symbol('\254')~~~Bifurcation~parameter~k), limits = c(-0.05, 1), expand = c(0, 0)) +
      scale_y_continuous(name = "State X", trans=pseudolog10_trans, limits = c(-1, 1.25), expand = c(0,0)) + 
      geom_hline(yintercept = 0, color = "darkgrey") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey20") +
      geom_line(data = df, aes(x = K, y = x_star), linewidth = .75, color = "#D55E00")  +
      geom_line(data = df, aes(x = K, y = unstable_x_star), linewidth = .75, linetype = "dashed", color = "#D55E00")  +
      geom_point(data = fixed_points, aes(x = k, y = x, shape = type), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1) +
      geom_line(data = df, aes (x = K, y = value, color = alpha, group = interaction(alpha, run)), alpha = .9, linewidth = .5) +
      scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
      scale_color_manual(values = colors, name = expression(alpha)) +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            plot.margin = unit(c(0, 0.75, 0, 0.5), "cm")))
  return(p)
}
plot_trajectories = function(A = c(2, 1.5, 1)) {
  (p1 = plot_traj_lin(A = A) +
    theme(legend.position = "right",
          legend.direction = "horizontal"))
  
  legend = get_legend(p1)
  
  p2 = plot_traj_nol(A = A)
  
  p = plot_grid(plot_grid(p1 + theme(legend.position = "None"), 
                      p2 + theme(legend.position = "None"), ncol = 1, align = "hv", labels = c("(a)", "(b)")), 
            legend,
            ncol = 1, rel_heights = c(6,1))
  
  return(p)
}
plot_trajectories()


ggsave("reports/paper/trajectories_diss.pdf", width = 9, height = 7, scale = 1.1)

variance_convergance = function() {
  data = list()
  
  A = c(2, 1.5)
  
  for (a in A) {
    df = read.csv(paste0("data/convergance_alpha", a, ".csv")) %>%
      filter(gamma_sample != 0) %>%
      mutate(alpha = a) %>%
      melt(measure.vars = c("var", "gamma_sample"))
    
    data = append(data, list(df))
  }
  
  df = purrr::reduce(data, bind_rows) %>%
    group_by(alpha, variable, run) %>%
    arrange(alpha, variable, run, time_step) %>%
    mutate(value_smoothed = slider::slide_dbl(value, mean, .before = 0, .after = 0),
           time_step = time_step - 40000) 
  df$alpha = factor(df$alpha, levels = c(1.5, 2))
  
  df_mean = df %>%
    group_by(alpha, variable, time_step) %>%
    summarize(value = mean(value_smoothed))
  colors_lines = c("#0072B2", "black")
  colors = c("#0072B2", "black")
  
  levels(df$variable) = c("Var~X", expression(gamma[" "~X]))
  levels(df_mean$variable) = c("Var~X", expression(gamma[" "~X]))
  
  (p = ggplot() + 
      geom_vline(xintercept = 0, color = "black", linewidth = .75) +
      geom_line(data = df, aes(x = time_step, y = value_smoothed, colour = alpha, group = interaction(run, alpha)), linewidth = .1, alpha = .5) +
      scale_color_manual(values = colors_lines, name = expression(alpha~" ")) +
      ggnewscale::new_scale_color() +
      geom_line(data = df_mean, aes(x = time_step, y = value, colour = alpha, group = interaction(alpha)), alpha = .9, linewidth = .6) +
      scale_color_manual(values = colors, name = expression(alpha~" ")) +
      scale_y_continuous(trans = "log10", name = "EWS", limits = c(0.0001, .25), expand = c(0,0)) +
      scale_x_continuous(name = "Simulation timestep", limits = c(0, 10000), expand = c(0,0), breaks = c(0, 2500, 7500)) + 
      facet_wrap(~variable, labeller=label_parsed) +
      theme(legend.position = "bottom",
            legend.direction = "horizontal"))
  
  return(p)
  
}
variance_convergance()

ggsave("reports/paper/variance_convergance_diss.pdf", width = 9, height = 4, scale = 1) 
ggsave("reports/paper/variance_convergance_diss.png", width = 9, height = 4, scale = 1, dpi = 600) 

#benchmark
benchmark_gL = function() {
  df = read.csv("data/benchmark_gammaL.csv") %>%
    filter(gamma !=0)
  
  df$alpha = factor(df$alpha, levels = c(2, 1.8, 1.5, 1.3))
  df$alpha = with(df, factor(alpha, levels = rev(levels(alpha))))
  
  colors_long =  c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9")
  
  (p1 = ggplot() + 
      geom_vline(xintercept = 1, color = "black", linewidth = .5) +
      geom_density_ridges(data = df, 
                          aes(x = gamma, y = alpha, fill = alpha, group = alpha),
                          jittered_points = TRUE,
                          position = position_points_jitter(width = 0.0, height = 0),
                          point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.4, size = 0.6) +
      scale_fill_manual(values = colors_long, name = expression(alpha)) +
      scale_x_continuous(limits = c(0.7, 1.3), name = expression(gamma[~L]), breaks = c(0.75, 1, 1.25), expand = c(0,0)) +
      scale_y_discrete(name = expression(alpha[~~L]), expand = expand_scale(mult = c(0.25, .7))) +
      theme(legend.position=c(.5,.05),
            legend.direction = "horizontal",
            legend.box.background = element_rect(fill='white', color = "darkgrey")))
  
  return(p1)
}
(p1 = benchmark_gL())
benchmark_gX = function() {
  A =c(2, 1.3, 1.5, 1.8)
  
  data = list()
  
  for (a in A){
    df = read.csv(paste0("data/benchmark_gammaX_alpha", a, ".csv")) %>%
      filter(W_in != 0,
             T_in !=1,
             gamma_sample < 2000)
    
    data = append(data, list(df))
    
    df_onetraj = read.csv(paste0("data/benchmark_gammaX_alpha", a, ".csv")) %>%
      filter(W_in != 0,
             gamma_sample < 2000)
    
    data = append(data, list(df_onetraj))
  }
  
  df = purrr::reduce(data, bind_rows) %>%
    group_by(alpha, W_in, T_in)  %>%
    summarise(gamma_sd = sqrt(var(gamma_sample, na.rm = TRUE)),
              gamma_sample = mean(gamma_sample)) %>% 
    mutate(gamma_true = 0.1*(1/(alpha))^(1/alpha),
           Error = ((gamma_sample - gamma_true)^2)/gamma_true,
           gamma_sd = gamma_sd/gamma_true)
  
  box = df %>%
    filter(W_in == 70,
           T_in == 5)
  
  df$alpha = factor(df$alpha,
                     levels = rev(c(2.0, 1.8, 1.5, 1.3)),
                     labels=rev(c('2.0'=parse(text=TeX('$\\alpha$ = 2')),
                                  '1.8'=parse(text=TeX('$\\alpha$ = 1.8')),
                                  '1.5'=parse(text=TeX('$\\alpha$ = 1.5')),
                                  '1.3'=parse(text=TeX('$\\alpha$ = 1.3')))))
  
  box$alpha = factor(box$alpha,
                    levels = rev(c(2.0, 1.8, 1.5, 1.3)),
                    labels=rev(c('2.0'=parse(text=TeX('$\\alpha$ = 2')),
                                 '1.8'=parse(text=TeX('$\\alpha$ = 1.8')),
                                 '1.5'=parse(text=TeX('$\\alpha$ = 1.5')),
                                 '1.3'=parse(text=TeX('$\\alpha$ = 1.3')))))
  
  (p = ggplot() + 
      geom_point(data = df, aes(x = as.factor(W_in), y = as.factor(T_in), fill = (Error), size = (abs(gamma_sd))), 
                 shape = 21, color = "black") +
      geom_point(data = box, aes(x = as.factor(W_in), y = as.factor(T_in), size =(abs(gamma_sd))), 
                 shape = 0, color = "darkgrey",size = 8) +
    scico::scale_fill_scico(palette = "lapaz", direction = -1,  name = "MSE") +
    scale_x_discrete(name = "Window size") +
    scale_y_discrete(name = "Number of trajectories") +
    scale_size_continuous(name = "Standard Error", range = c(1, 7)) +
    facet_wrap(~(factor(alpha)), labeller=label_parsed)  +
    theme(strip.background = element_rect(fill = "transparent")))
  
  ggplot() +
    geom_jitter(data = df[df$alpha == 2 ,], aes(x = as.factor(W_in), y = as.factor(T_in), color = delta)) +
    facet_wrap(~alpha, scales = "free") +
    scico::scale_color_scico(palette = "roma", midpoint = 0)
  
  return(p)
}
p2 = benchmark_gX()
plot_grid(p1, p2, nrow = 1, labels = c("(a)", "(b)"), rel_widths = c(4,5))
ggsave("reports/paper/benchmark_diss.pdf", width = 10) 


