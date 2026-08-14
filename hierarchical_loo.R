#######################################################################
#                                                                     #
# R script for PSIS-LOO model comparison across hierarchical model    #
# iterations (hierarchical_mm<code>.Rsave outputs of                  #
# hierarchical_data_and_running.R), in the spirit of loo1.R/loo2.R    #
# from the non-hierarchical pipeline (github.com/mmorenozam/ecbfm).   #
#                                                                     #
# There, PSIS-LOO was computed once per (dataset, model iteration)    #
# fit and faceted by dataset. Here, one hierarchical fit already      #
# covers all 28 trials at once (log_lik has shape [G,Tmax,8], with    #
# entries beyond T[g] left at the -1 sentinel from                    #
# hierarchical_mm<code>.stan's generated quantities block), so this   #
# script:                                                              #
#   - computes a per-trial PSIS-LOO by subsetting each fit's log_lik  #
#     to that trial's valid entries (facet analogous to loo1.R), and  #
#   - computes one pooled PSIS-LOO per model over all trials, which   #
#     is the single number to compare mechanism variants against.     #
#                                                                     #
# The [G,Tmax,8] -> flat-N column order rstan/loo use for a Stan      #
# `array[G,Tmax,8]` matches R's own column-major array layout (first  #
# declared index fastest), i.e. identical to `array(1:N,              #
# dim=c(G,Tmax,8))`; that identity is what the index masks below rely #
# on to pick out a given trial's columns from the flattened log_lik.  #
#                                                                     #
# Usage: put this in the directory containing your hierarchical_mm*   #
# .Rsave files and Rscript hierarchical_loo.R. Writes loo_results.csv #
# and psis_loo.pdf (per-trial) / psis_loo_pooled.pdf (per-model).     #
#                                                                     #
#######################################################################

library(loo)
library(tidyverse)
library(cowplot)

setwd("/home")       # directory name, change accordingly if needed
dpath <- "/home/"    # directory where hierarchical_mm*.Rsave files live
outpth <- "/home/"   # directory where outputs are written

source("trial_data.R")

G <- 28
trials <- lapply(1:G, get_trial)
T_vec <- sapply(trials, function(tr) tr$T)
lbl_vec <- sub("_$", "", sapply(trials, function(tr) tr$lbl))
Tmax <- max(T_vec)

# flat_idx[[g]] = column indices (into the flattened G*Tmax*8 log_lik
# array) holding trial g's real (non-padding) log-lik values
idx_cube <- array(1:(G * Tmax * 8), dim = c(G, Tmax, 8))
flat_idx <- lapply(1:G, function(g) as.vector(idx_cube[g, 1:T_vec[g], ]))
pooled_idx <- unlist(flat_idx)

fname <- list.files(dpath, pattern = "^hierarchical_mm.*\\.Rsave$")
if (length(fname) == 0) stop("No hierarchical_mm*.Rsave files found in ", dpath)

rows <- list()
for (f in fname) {
  mi <- gsub("hierarchical_mm|\\.Rsave", "", f)
  load(paste0(dpath, f)) # loads `fit`

  ll_arr <- loo::extract_log_lik(fit, "log_lik", merge_chains = FALSE) # [iter, chain, N]
  r_eff <- loo::relative_eff(exp(ll_arr), cores = 1)

  summarize_loo <- function(cols, trial_lbl) {
    lr <- loo::loo(ll_arr[, , cols, drop = FALSE], r_eff = r_eff[cols], cores = 1)
    data.frame(
      model = mi,
      trial = trial_lbl,
      elpd_loo = lr$estimates["elpd_loo", "Estimate"],
      elpd_se = lr$estimates["elpd_loo", "SE"],
      looic = lr$estimates["looic", "Estimate"],
      looic_se = lr$estimates["looic", "SE"],
      max_pareto_k = max(lr$diagnostics$pareto_k),
      n_high_k = sum(lr$diagnostics$pareto_k > 0.7)
    )
  }

  for (g in 1:G) rows[[length(rows) + 1]] <- summarize_loo(flat_idx[[g]], lbl_vec[g])
  rows[[length(rows) + 1]] <- summarize_loo(pooled_idx, "pooled (all trials)")
}

loo_results <- do.call(rbind, rows)
write.csv(loo_results, paste0(outpth, "loo_results.csv"), row.names = FALSE)

if (any(loo_results$n_high_k > 0)) {
  warning(sum(loo_results$n_high_k > 0), " (model, trial) combinations have ",
          "Pareto k > 0.7 observations -- PSIS-LOO is unreliable there; see ",
          "loo_results.csv's n_high_k column before trusting those comparisons.")
}

# per-trial comparison, faceted like loo1.R/loo2.R
per_trial <- loo_results %>% filter(trial != "pooled (all trials)")
l_plt <- ggplot(per_trial, aes(y = looic, x = model)) +
  geom_point(color = "#D43F3AFF") +
  geom_errorbar(aes(ymin = looic - looic_se, ymax = looic + looic_se), color = "#EEA2367F", linewidth = 0.4) +
  facet_wrap(~trial, scales = "free_y", ncol = 4) +
  xlab("MI( )") +
  ylab("PSIS-LOO (looic)") +
  theme_bw() +
  theme(strip.text.x = element_text(hjust = -0.01),
        panel.grid.major = element_line(colour = "grey90", linetype = "dashed", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95))
save_plot(paste0(outpth, "psis_loo.pdf"), l_plt, base_height = 18, base_width = 14)

# pooled, single-panel comparison across model iterations
pooled <- loo_results %>% filter(trial == "pooled (all trials)")
p_plt <- ggplot(pooled, aes(y = looic, x = model)) +
  geom_point(color = "#D43F3AFF", size = 2) +
  geom_errorbar(aes(ymin = looic - looic_se, ymax = looic + looic_se), color = "#EEA2367F", width = 0.2) +
  xlab("MI( )") +
  ylab("Pooled PSIS-LOO (looic), all 28 trials") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90", linetype = "dashed", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95))
save_plot(paste0(outpth, "psis_loo_pooled.pdf"), p_plt, base_height = 6, base_width = 8)
