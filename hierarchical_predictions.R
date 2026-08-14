#######################################################################
#                                                                     #
# R script for plotting posterior predictive trajectories from        #
# hierarchical model fits (hierarchical_mm<code>.Rsave outputs of     #
# hierarchical_data_and_running.R), in the spirit of predictions1.R   #
# from the non-hierarchical pipeline (github.com/mmorenozam/ecbfm).   #
#                                                                     #
# There, one Stan fit covered one trial, so x_hat's posterior was     #
# pulled directly and rescaled by that trial's own `scl`. Here, one   #
# hierarchical fit covers all 28 trials, so this script pulls out     #
# each trial's [1:T[g]] slice of x_hat_rep (already unscaled to raw   #
# concentration units in hierarchical_mm<code>.stan's generated       #
# quantities block) and plots it the same way: posterior median +     #
# 95% interval ribbon per state variable, against the observed data,  #
# with the t0 initial condition prepended.                            #
#                                                                     #
# Usage: put this in the directory containing your hierarchical_mm*   #
# .Rsave files and Rscript hierarchical_predictions.R. Writes one PDF #
# per (trial, model iteration) as <trial-label>_mm<code>.pdf.         #
#                                                                     #
#######################################################################

library(rstan)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/home")       # directory name, change accordingly if needed
dpath <- "/home/"    # directory where hierarchical_mm*.Rsave files live
outpth <- "/home/"   # directory where plots are written

source("trial_data.R")

G <- 28
trials <- lapply(1:G, get_trial)
T_vec <- sapply(trials, function(tr) tr$T)
lbl_vec <- sub("_$", "", sapply(trials, function(tr) tr$lbl))
Tmax <- max(T_vec)

metb <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)")

# builds the tidy dataframe (median/95% interval ribbon + observed data +
# t0 initial condition) for one trial's posterior predictive trajectory,
# analogous to swf1's data-prep block in predictions1.R
trial_pred_df <- function(x_hat_rep, g, mi) {
  Tg <- T_vec[g]
  tr <- trials[[g]]
  x_g <- x_hat_rep[, g, 1:Tg, ] # [iterations, Tg, 8] (g and any size-1 dims auto-dropped)

  sp <- list()
  for (n in 1:8) {
    q <- apply(x_g[, , n], 2, quantile, probs = c(0.025, 0.5, 0.975))
    df <- data.frame(t(q))
    colnames(df) <- c("lb", "median", "ub")
    df$obs <- tr$x[, n]
    df$time <- tr$ts
    df <- rbind(data.frame(lb = tr$x0[n], median = tr$x0[n], ub = tr$x0[n], obs = tr$x0[n], time = 0), df)
    df$met <- metb[n]
    sp[[n]] <- df
  }
  dfout <- do.call(rbind, sp)
  dfout$data <- lbl_vec[g]
  dfout$model <- mi
  dfout
}

fname <- list.files(dpath, pattern = "^hierarchical_mm.*\\.Rsave$")
if (length(fname) == 0) stop("No hierarchical_mm*.Rsave files found in ", dpath)

for (f in fname) {
  mi <- gsub("hierarchical_mm|\\.Rsave", "", f)
  load(paste0(dpath, f)) # loads `fit`

  x_hat_rep <- rstan::extract(fit, pars = "x_hat_rep")$x_hat_rep # [iterations, G, Tmax, 8]

  for (g in 1:G) {
    dfout <- trial_pred_df(x_hat_rep, g, mi)

    p1 <- ggplot(data = dfout, aes(x = time, y = median)) +
      geom_line(color = "#D43F3AFF", linewidth = 0.7) +
      geom_point(aes(x = time, y = obs)) +
      geom_line(aes(x = time, y = ub), color = "#EEA2367F", linetype = "dashed", linewidth = 1) +
      geom_line(aes(x = time, y = lb), color = "#EEA2367F", linetype = "dashed", linewidth = 1) +
      geom_ribbon(aes(ymin = lb, ymax = ub), fill = "#EEA2367F", alpha = 0.3) +
      facet_wrap(~met, scales = "free_y", ncol = 4) +
      theme_bw() +
      theme(strip.text.x = element_text(hjust = -0.01),
            panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 14)) +
      xlab("time (h)") +
      ylab(expression(paste("mg g(pulp) "^{-1})))

    save_plot(paste0(outpth, lbl_vec[g], "_mm", mi, ".pdf"), p1, base_height = 8.3, base_width = 17.2)
  }
}
