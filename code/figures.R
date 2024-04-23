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

setwd("~/02_Science/alphastableEWS_paper")

#overview


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
  
  colors = c("#000000", "#98C6EA")
  colors_means = c("#E37222", "#98C6EA")
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
      geom_point(data = df, aes(x = k, y = variance, fill = system), 
                 color = "grey20", stroke = 0.001, shape = 21, alpha = .5, position = position_jitter(width = .15), size = 1.5) +
      geom_point(data = df_mean, aes(x = k, y = mean, fill = system, size = system), color = "black",  shape = 21, stroke = .5) +
      scale_color_manual(values = line_color, name = " ", guide = guide_legend(title = NULL)) +
      scale_linetype_manual(values = c("Theory" = "dashed"), name = "") +
      scale_size_manual(values = sizes, name = "System")  +
      scale_fill_manual(values = colors_means, name = "System") +
      scale_y_continuous(trans = 'log10', name = "EWS Var[X]") +
      scale_x_continuous(trans = 'log10', name = "Bifurcation parameter k", labels = scales::comma, breaks = c(0.01, 1, 100)) +
      theme(text = element_text(size = 12),
            legend.position = c(.75, .75),
            legend.title = element_blank(),
            legend.background = element_rect(fill='white', color = "darkgrey"),
            legend.box.background = element_rect(fill='transparent', color = NA)) + guides(size = "none"))
  
  
  return(p)
  
}
p1 = plot_var_eq()
alpha_stable_dist = function() {
  x = seq(-5, 5, by = .1)
  
  df = data.frame("x" = x,
                  "0.5" = dstable(x, tail = 0.5),
                  "1" = dstable(x, tail = 1),
                  "1.5" = dstable(x, tail = 1.5),
                  "2" = dstable(x, tail = 2)) %>%
    melt(id.vars = "x")

  colors = c("#E37222", "#532d5b", "#090979", "#64A0C8")
  
  (p2 = ggplot() + theme_bw() +
    geom_line(data = df, 
              aes(x = x, y = value, color = variable),
              linewidth = 1, alpha = .9) +
    scale_y_continuous(name = "p(X)") +
    theme(text = element_text(size = 12), 
          legend.position = c(.2, .75),
          legend.background = element_rect(fill='white', color = "darkgrey"),
          legend.box.background = element_rect(fill='transparent', color = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),  
          plot.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_rect(fill = "transparent", color = NA)) +
    scale_color_manual(values = colors, name = expression(alpha),
                               labels = c(0.5, 1, 1.5, 2)))
  
  return(p2)
}
p2 = alpha_stable_dist()
plot_grid(p1 , p2, nrow = 1, labels = c("(a)", "(c)"), rel_widths = c(1, .8))
ggsave("reports/paper/theory.png", width = 8.2, height = 4,  bg = "white")
ggsave("reports/paper/theory.pdf", width = 8.2, height = 4,  bg = "white", dpi = 600)

#gamma equilibrium runs
plot_gamma_eq = function(A =c(1.3, 1.5, 1.8, 2)) {
  
  data = list()
  
  for (a in A) {
    df = read.csv(paste0("data/gamma_lin_eq_alpha", a, "_rev.csv")) %>%
      filter(k != 0) %>%
      mutate(system = "Linear",
             alpha = paste0(a))
    
  
    df_nol = read.csv(paste0("data/gamma_nol_eq_alpha", a, "_rev.csv")) %>%
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
  colors_means = c("#E37222", "#98C6EA")
  line_color = c("True \U03B3" = "grey50")
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

  
  (p = ggplot() + theme_bw() +
      geom_line(data = df, aes(x = k, y = gamma_true, color = "Theory", linetype = "Theory"), linewidth =.75) +
      facet_wrap(~rev(alpha), labeller = label_parsed) + 
      geom_point(data = df, aes(x = k, y = gamma_sample, fill = system), 
                 color = "grey20", stroke = 0.001, shape = 21, alpha = .5, position = position_jitter(width = .15), size = 1.5) +
      geom_point(data = df_mean, aes(x = k, y = mean, fill = system, size = system), color = "black",  shape = 21, stroke = .5) +
      scale_color_manual(values = line_color, name = " ") +
      scale_size_manual(values = sizes, name = "System")  +
      scale_fill_manual(values = colors_means, name = "System") +
      scale_linetype_manual(values = c("Theory" = "dashed"), name = "") +
      scale_y_continuous(trans = 'log10', name = expression(EWS~gamma[X])) +
      scale_x_continuous(trans = 'log10', name = "Bifurcation parameter k", labels = scales::comma, breaks = c(0.01, 1, 100)) +
      theme(text = element_text(size = 12)) + 
      guides(size = "none"))
  
  
  return(p)
  
}
plot_gamma_eq()
ggsave("reports/paper/gamma_equilibrium.png", width = 7, height = 5)
ggsave("reports/paper/gamma_equilibrium.pdf", width = 7, height = 5, dpi = 600)


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
  
  legend_labels = rev(c(TeX("$\\alpha$ = 2"), TeX("$\\alpha$ = 1.8"), TeX("$\\alpha$ = 1.5"), TeX("$\\alpha$ = 1.3")))
  
  colors = c("#E37222", "#532d5b",  "#090979",  "#64A0C8")
  (p = ggplot() + theme_bw() +
      geom_line(data = df, aes(x = k, y = gamma_smoothed, group = interaction(run, alpha), color = alpha), alpha = .5, linewidth = .05) +
      geom_line(data = df[df$system == "Linear",], aes(x = k, y = gamma_true, group = interaction(system, alpha), color = alpha, linetype = "Theory"), alpha = .8,  linewidth = .6) +
      geom_line(data = df_mean, aes(x = k, y = mean_gamma, color = alpha), alpha = .9, linewidth = .7) +
      scale_color_manual(values = colors,  name = expression(Estimated~gamma), labels = legend_labels) + 
      scale_linetype_manual(values = c("Theory" = "dashed"), name = "") + 
      facet_grid(rows = vars(system), cols = vars(alpha), scales = "free", labeller = label_parsed) + 
      scale_x_continuous(limits = c(0, 5.25), expression(symbol('\254')~~~Bifurcation~parameter~k)) +
      scale_y_continuous(name = expression(gamma[X])) +
      theme(panel.grid.minor = element_blank()))
  
  return(p)
}
plot_gamma_neq() 
ggsave("reports/paper/gamma_noneq.png",  width = 8, height = 4, dpi = 1000)
ggsave("reports/paper/gamma_noneq.pdf",  width = 8, height = 4, dpi = 1000)

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
  colors = c("#532d5b", "#090979", "#64A0C8")
  line_color = c("Fixed points X*" = "black")
  
  (p = ggplot() + theme_bw() + 
      scale_x_continuous(name = expression(symbol('\254')~~~Bifurcation~parameter~k), limits = c(-0.05, 1), expand = c(0, 0)) +
      scale_y_continuous(name = "State X", trans = pseudolog10_trans) + 
      geom_hline(yintercept = 0, color = "darkgrey") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey20") +
      geom_segment(aes(x = -0.05, y = 0, yend = 0, xend = 0.1), color = "black", linetype = "dashed", linewidth = .75) + 
      geom_line(data = df, aes(x = K, y = x_star, color = "Fixed points X*"), linewidth = .75)  +
      scale_color_manual(values = line_color, name = "") +
      ggnewscale::new_scale_color() +
      geom_line(data = df, aes (x = K, y = value, color = alpha, group = interaction(alpha, run)), alpha = .75, linewidth = .5) +
      scale_color_manual(values = colors, name = expression(alpha)) +
      theme(text = element_text(size = 15),
            legend.position = "bottom",
            legend.direction = "horizontal"))
  
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
             run ==1 )
    
    
    
    data = append(data, list(traj))
  }
  
  df = purrr::reduce(data, bind_rows)
  
  df$alpha = factor(df$alpha, levels = c(2, 1.5, 1, 0.5))
  df$alpha = with(df, factor(alpha, levels = rev(levels(alpha))))
  colors = c("#532d5b", "#090979", "#64A0C8")
  colors_long = c("#64A0C8", "#090979", "#532d5b", "#E37222")
  
  (p = ggplot() + theme_bw() + 
      scale_x_continuous(name = expression(symbol('\254')~~~Bifurcation~parameter~k), limits = c(-0.05, 1), expand = c(0, 0)) +
      scale_y_continuous(name = "State X", trans = pseudolog10_trans, limits = c(-0.4, 1.25)) + 
      geom_hline(yintercept = 0, color = "darkgrey") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey20") +
      geom_line(data = df, aes(x = K, y = x_star), linewidth = .75)  +
      geom_line(data = df, aes(x = K, y = unstable_x_star), linewidth = .75, linetype = "dashed")  +
      geom_line(data = df, aes (x = K, y = value, color = alpha, group = interaction(alpha, run)), alpha = .75, linewidth = .5) +
      scale_color_manual(values = colors, name = expression(alpha)) +
      theme(text = element_text(size = 15)))
  
  return(p)
}
plot_trajectories = function(A = c(2, 1.5, 1)) {
  p1 = plot_traj_lin(A = A)
  legend = get_legend(p1)
  
  
  p2 = plot_traj_nol(A = A)
  
  p = plot_grid(plot_grid(p1 + theme(legend.position = "None"), 
                      p2 + theme(legend.position = "None"), ncol = 1, align = "hv", labels = c("(a)", "(b)")), 
            legend,
            ncol = 1, rel_heights = c(6,1))
  
  return(p)
}
plot_trajectories()
ggsave("reports/paper/trajectories.png", width = 7, height = 5, bg = "white")
ggsave("reports/paper/trajectories.pdf", width = 7, height = 5, bg = "white", dpi = 600)

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
  colors_lines = c("#532d5b", "#eb975c")
  colors = c("#532d5b", "#E37222")
  
  levels(df$variable) = c("Var(X)", expression(gamma[X]))
  levels(df_mean$variable) = c("Var(X)", expression(gamma [X]))
  
  (p = ggplot() + theme_bw() + 
      geom_line(data = df, aes(x = time_step, y = value_smoothed, colour = alpha, group = interaction(run, alpha)), linewidth = .1, alpha = .5) +
      scale_color_manual(values = colors_lines, name = expression(alpha)) +
      ggnewscale::new_scale_color() +
      geom_line(data = df_mean, aes(x = time_step, y = value, colour = alpha, group = interaction(alpha)), alpha = .9, linewidth = .6) +
      scale_color_manual(values = colors, name = expression(alpha)) +
      scale_y_continuous(trans = "log10", name = "EWS", limits = c(0.0001, .2)) +
      scale_x_continuous(name = "Simulation timestep", limits = c(0, 10000)) + 
      facet_wrap(~variable, labeller=label_parsed) +
      theme(text = element_text(size = 12),
            legend.background = element_rect(fill='transparent', color = NA),
            legend.box.background = element_rect(fill='transparent', color = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),  
            plot.background = element_rect(fill = "transparent", colour = NA),
            strip.background = element_rect(fill = "transparent")))
  
  return(p)
  
}
variance_convergance()
ggsave("reports/paper/variance_convergance.png", width = 7, height = 3, bg = "white", dpi = 600) 
ggsave("reports/paper/variance_convergance.pdf", width = 7, height = 3, bg = "white", dpi = 600) 

#benchmark
benchmark_gL = function() {
  df = read.csv("data/benchmark_gammaL.csv") %>%
    filter(gamma !=0)
  
  df$alpha = factor(df$alpha, levels = c(2, 1.8, 1.5, 1.3))
  df$alpha = with(df, factor(alpha, levels = rev(levels(alpha))))
  
  colors_long = c("#E37222", "#532d5b", "#090979", "#64A0C8"  )
  (p1 = ggplot() + theme_bw() + 
      geom_vline(xintercept = 1, color = "black", linewidth = .5) +
      geom_density_ridges(data = df, 
                          aes(x = gamma, y = alpha, fill = alpha, group = alpha),
                          jittered_points = TRUE,
                          position = position_points_jitter(width = 0.0, height = 0),
                          point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.4, size = 0.6) +
      scale_fill_manual(values = colors_long, name = expression(alpha)) +
      scale_x_continuous(limits = c(0.7, 1.3), name = expression(gamma[~L]), breaks = c(0.75, 1, 1.25)) +
      scale_y_discrete(name = expression(alpha[~~L]), expand = expand_scale(mult = c(0.25, .7))) +
      theme(text = element_text(size = 15),
            legend.position=c(.5,.05),
            legend.direction = "horizontal",
            legend.background = element_rect(fill='transparent', color = NA),
            legend.box.background = element_rect(fill='white', color = "darkgrey"),
            panel.background = element_rect(fill = "transparent", colour = NA),  
            plot.background = element_rect(fill = "transparent", colour = NA),
            strip.background = element_rect(fill = "transparent", color = NA)))
  
  return(p1)
}
p1 = benchmark_gL()
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
  
  (p = ggplot() + theme_bw() + 
      geom_point(data = df, aes(x = as.factor(W_in), y = as.factor(T_in), fill = (Error), size = (abs(gamma_sd))), 
                 shape = 21, color = "black") +
      geom_point(data = box, aes(x = as.factor(W_in), y = as.factor(T_in), size =(abs(gamma_sd))), 
                 shape = 0, color = "darkgrey",size = 8) +
    scico::scale_fill_scico(palette = "lapaz", direction = -1,  name = "MSE") +
    scale_x_discrete(name = "Window size") +
    scale_y_discrete(name = "Number of trajectories") +
    scale_size_continuous(name = "Standard Error", range = c(1, 7)) +
    facet_wrap(~(factor(alpha)), labeller=label_parsed)  +
    theme(text = element_text(size = 12),
          strip.background = element_rect(fill = "transparent")))
  
  ggplot() +
    geom_jitter(data = df[df$alpha == 2 ,], aes(x = as.factor(W_in), y = as.factor(T_in), color = delta)) +
    facet_wrap(~alpha, scales = "free") +
    scico::scale_color_scico(palette = "roma", midpoint = 0)
  
  return(p)
}
p2 = benchmark_gX()
plot_grid(p1, p2, nrow = 1, labels = c("(a)", "(c)"), rel_widths = c(4,5))
ggsave("reports/paper/benchmark.png", width = 8, height = 4, bg = "white") 
ggsave("reports/paper/benchmark.pdf", width = 8, height = 4, bg = "white", dpi = 600) 


