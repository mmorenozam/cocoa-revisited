# cocoa-revisited

Revisiting the PhD-era cocoa bean fermentation kinetic model to fit it
**hierarchically** across all 28 fermentation trials, instead of one trial
at a time as originally done. This file is a running log of what's been
done, so picking the project back up doesn't require reconstructing
context from scratch.

## Background

The original model and data come from:

- Moreno-Zambrano, M., Ullrich, M. S., Hütt, M-T. (2022). *Exploring cocoa
  bean fermentation mechanisms by kinetic modelling*. Royal Society Open
  Science, 9(2), 210274. https://doi.org/10.1098/rsos.210274
  Code/data: https://github.com/mmorenozam/ecbfm
  (a baseline mechanistic ODE model, extended with 5 candidate mechanisms
  M1-M5 added in isolation/combination — 32 model variants "mm0" through
  "mm12345" — each fit *separately per trial* with Stan/rstan.)

The baseline model was originally intended to be fit **hierarchically**
across trials (partial pooling), but that was never carried out during the
PhD. That's what this repo is doing now, following the approach in:

- Rosenbaum, B., Raatz, M., Weithoff, G., Fussmann, G. F., Gaedke, U.
  (2019). *Estimating Parameters From Multiple Time Series of Population
  Dynamics Using Bayesian Inference*. Frontiers in Ecology and Evolution.
  https://doi.org/10.3389/fevo.2018.00234
  (partial-pooling Bayesian ODE fitting across replicate time series in
  Stan: population-level log-scale hyperpriors, per-replicate parameter
  vectors, trajectory matching.)

## Repo contents

**Original (non-hierarchical), per-trial pipeline:**
- `baseline_model.stan` — the 24-kinetic-parameter ODE model (no M1-M5
  mechanisms added; this is model iteration `mm0` in the original scheme).
- `data_and_running.R` — fits one `mm<code>.stan` variant to one trial at a
  time (CLI: `Rscript data_and_running.R d mi ad st mt s ni`). Now sources
  `trial_data.R` instead of inlining the data (behavior unchanged from the
  original; see "Refactor" below).
- `trial_data.R` — `get_trial(d)`, returns the raw data (T, x0, t0, ts, x,
  lbl, x1, scl) for trial `d` (1:28). Extracted verbatim from
  `data_and_running.R`'s original inline `if (args[1]==d)` blocks (checked
  programmatically to be byte-identical apart from the `args[1]` → `d`
  rename) so both pipelines share one source of truth for the data.

**New hierarchical pipeline:**
- `hierarchical_mm0.stan` — hierarchical (partial-pooling) version of
  `baseline_model.stan`. All 28 trials fit in one Stan call; every trial
  gets its own 24-parameter vector, drawn from shared population-level
  log-normal hyperpriors (non-centered parameterization). See the file's
  header comment for the full design writeup, and "Design decisions" below
  for the reasoning.
- `hierarchical_data_and_running.R` — fits `hierarchical_mm<mi>.stan` to
  all 28 trials at once (CLI: `Rscript hierarchical_data_and_running.R mi
  ad st mt s ni`). Currently only `mi=0` is wired up (`ini()` matches
  `hierarchical_mm0.stan`'s 24 parameters); the script stops loudly for
  any other `mi` until its `ini()` is extended to match that variant.
- `hierarchical_loo.R` — PSIS-LOO model comparison across whichever
  `hierarchical_mm*.Rsave` fits are present, both per-trial (faceted,
  like the original `loo1.R`/`loo2.R`) and pooled across all 28 trials
  (the number that actually matters for comparing mechanism variants).
  Flags any Pareto-k > 0.7. Outputs `loo_results.csv`, `psis_loo.pdf`,
  `psis_loo_pooled.pdf`.
- `hierarchical_predictions.R` — posterior predictive plots (median + 95%
  interval ribbon vs. observed data) for every trial in every
  `hierarchical_mm*.Rsave` file present, styled like the original
  `predictions1.R`. One PDF per (trial, model) as `<trial>_mm<code>.pdf`.

**Not ported yet** — `PCAs.R` and `heatmap.R` from the original `ecbfm`
repo compare parameter estimates across mechanism variants using per-trial
metadata (cultivar, country, fermentation method, controlled temperature,
turning) that isn't in this repo, and only make sense once more than one
hierarchical mechanism variant exists. Revisit once both are available.

## Design decisions for the hierarchical model

(Decided in conversation, recorded here so they don't need re-deciding.)

- **Pooling**: all 28 trials as one exchangeable population (not grouped
  by country/method) — mirrors Rosenbaum et al.'s approach directly, and
  keeps the most data per population-level parameter.
- **Likelihood**: log-normal on the raw (not linear-normal on scaled)
  state, per Rosenbaum et al. — better suited to concentrations/microbial
  counts spanning multiple orders of magnitude than the original's normal-
  on-scaled-state likelihood.
- **Numerical conditioning kept from the original**: each trial's state is
  still divided by its own per-variable max (`scl`, computed exactly as in
  `baseline_model.stan`) before ODE integration, since dropping that and
  integrating raw magnitudes (10⁻⁹ to 10²) risked RK45 instability. The
  log-normal likelihood is applied to the *scaled* state; `x_hat_rep` in
  generated quantities is rescaled back to raw units for convenience.
- **Zero handling**: some trials have literal zero measurements (fully
  consumed substrate), which log-normal can't take the log of. A small
  `eps` (data value, default `1e-3` in scaled units) is added before
  logging both data and prediction.
- **Non-centered parameterization**: given many trials have only 4-16
  observations, a centered hierarchical parameterization would likely
  produce funnel pathologies in NUTS; used `z_theta` raw deviates instead.
- **Stan ODE interface modernized**: `integrate_ode_rk45` (deprecated) →
  `ode_rk45_tol`, same tolerances/step cap as the original
  (`rel_tol=abs_tol=1e-6`, `max_num_steps=10000`).
- **`x0[8]` (AAB biomass) is never literally zero** in any of the 28
  trials — confirmed by inspection — which is what keeps `v5`'s
  denominator (`ks5*x8 + x4`) away from a `0/0` singularity even in trials
  where `x4` (lactic acid) starts at exactly zero. This was a deliberate
  choice in the original PhD data (tiny scientific-notation floors instead
  of exact zero), carried through unchanged into `trial_data.R`.

## Known gaps / next steps

1. **No compile check has been possible.** There's no R/Stan toolchain in
   the dev environment this was built in — `hierarchical_mm0.stan` was
   reviewed by hand (types, array shapes, ragged-loop bounds) but never
   run through `stanc`/`rstan::stan_model()`. Do that first.
2. **`renv.lock` only pins `renv` itself.** Add `rstan` (and `loo`,
   `tidyverse`, `cowplot` for the post-processing scripts) and snapshot
   wherever this actually gets run.
3. **Only `mm0` (baseline, no mechanisms) is hierarchical so far.** The
   remaining 31 mechanism-combination variants (`mm1` ... `mm12345`) need
   their own `hierarchical_mm<code>.stan` files, following
   `hierarchical_mm0.stan`'s structure (each mechanism changes how many
   `mu`/`ks`/`yc`/`ev` parameters exist, so the parameter block, `theta`
   indexing, and matching `ini()` in the runner all need updating per
   variant) — `hierarchical_loo.R` and `hierarchical_predictions.R` will
   work unchanged on top of those once they exist, since they only touch
   `log_lik`/`x_hat_rep`, which stay `[G,Tmax,8]`-shaped regardless of
   `theta`'s size.
4. **This is computationally heavy**: 28 ragged ODE solves per leapfrog
   step, ~721 parameters for `mm0` alone. Meant for a server/cluster, not
   an interactive session — the original scripts' `setwd('/home')` reflects
   that.
5. `PCAs.R`/`heatmap.R` equivalents — needs a per-trial metadata table
   (cultivar, country, method, controlled temperature, turning) that
   doesn't exist in this repo yet, plus more than one hierarchical variant
   to compare.
