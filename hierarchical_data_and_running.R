#######################################################################
#                                                                     #
# R script for fitting a hierarchical (partial-pooling) version of    #
# the cocoa bean fermentation kinetic model jointly across all 28     #
# fermentation trials.                                                #
#                                                                     #
# Unlike data_and_running.R (which fits baseline_model.stan / one of  #
# the mm<code>.stan variants to a single trial at a time), this       #
# script pools all 28 trials into one Stan call: every trial gets its #
# own parameter vector, and those vectors share a common              #
# population-level log-normal distribution (see the header comment    #
# in hierarchical_mm0.stan for the modeling rationale). It follows     #
# the same mm<code> naming/model-iteration convention as              #
# data_and_running.R, but for hierarchical_mm<code>.stan files.       #
#                                                                     #
# This is computationally heavy (28 ragged ODE solves per leapfrog    #
# step, hundreds of parameters) and is intended to run on a           #
# server/cluster, not interactively.                                  #
#                                                                     #
# Usage:                                                              #
#   $ Rscript hierarchical_data_and_running.R mi ad st mt s ni        #
# where (same meaning as in data_and_running.R):                      #
#   mi = model iteration code (matches hierarchical_mm<mi>.stan)      #
#   ad = adapt_delta, st = step_size, mt = max_treedepth,              #
#   s  = seed, ni = number of iterations                              #
#                                                                     #
# NOTE: the `ini()` below matches hierarchical_mm0.stan's 24-kinetic- #
# parameter structure only. Every additional mechanism (mm1, mm2,     #
# ... in your mechanism-combination scheme) changes how many mu/ks/yc #
# /ev parameters exist, so ini() must be updated to match whichever   #
# hierarchical_mm<mi>.stan you're running -- the guard below stops    #
# the script rather than silently fitting the wrong number of         #
# parameters.                                                         #
#                                                                     #
#######################################################################

setwd('/home') # directory name, change accordingly if needed
rm(list = ls())
args = commandArgs(TRUE)
mi = args[1]

library(rstan)

rstan_options(auto_write = TRUE)

source("trial_data.R")

if (mi != "0") {
  stop("ini() below only matches hierarchical_mm0.stan (24 kinetic parameters). ",
       "Extend the ini() list (and this check) to match hierarchical_mm", mi, ".stan ",
       "before running model iteration ", mi, ".")
}

G <- 28
trials <- lapply(1:G, get_trial)

Tmax <- max(sapply(trials, function(tr) tr$T))

T_vec <- sapply(trials, function(tr) tr$T)
x0_mat <- t(sapply(trials, function(tr) tr$x0))
scl_mat <- t(sapply(trials, function(tr) tr$scl))

ts_pad <- matrix(0, nrow = G, ncol = Tmax)
x_pad <- array(0, dim = c(G, Tmax, 8))
for (g in 1:G) {
  Tg <- trials[[g]]$T
  ts_pad[g, 1:Tg] <- trials[[g]]$ts
  x_pad[g, 1:Tg, ] <- trials[[g]]$x
  if (Tg < Tmax) {
    # padding times are never read by the Stan model (only [1:T[g]] is
    # used); kept strictly increasing purely for tidiness/debuggability
    ts_pad[g, (Tg+1):Tmax] <- ts_pad[g, Tg] + (1:(Tmax-Tg))
  }
}

stan_data <- list(
  G = G,
  Tmax = Tmax,
  T = T_vec,
  t0 = 0,
  x0 = x0_mat,
  ts = ts_pad,
  x = x_pad,
  scl = scl_mat,
  eps = 1e-3
)

ini <- function() {
  list(
    mu_log_theta = rnorm(24, log(0.5), 0.3),
    sigma_log_theta = abs(rnorm(24, 0.3, 0.1)),
    z_theta = matrix(rnorm(24 * G, 0, 0.3), nrow = 24, ncol = G),
    sigma = abs(rnorm(1, 0.5, 0.2))
  )
}

fit = stan(paste0("hierarchical_mm", mi, ".stan"),
           data = stan_data,
           control = list(adapt_delta = as.numeric(args[2]),
                           stepsize = as.numeric(args[3]),
                           max_treedepth = as.numeric(args[4])),
           warmup = round(as.numeric(args[6])*1/3, 0),
           init = ini,
           refresh = 5,
           cores = min(4, parallel::detectCores()),
           chains = 4, iter = as.numeric(args[6]), seed = as.numeric(args[5]))

save(fit, file = paste0("hierarchical_mm", mi, ".Rsave")) # saving output
