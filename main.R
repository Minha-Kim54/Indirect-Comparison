# main.R
# -----------------------------------------------------------------------
# Main execution script.
# Loads all source files, builds the scenario grid, runs the simulation,
#
# Run order:
#   1. scenario_setup.R
#   2. data_generation.R
#   3. methods.R
#   4. simulation_run.R
# -----------------------------------------------------------------------

library(dplyr)
library(tibble)
library(mvtnorm)
library(Matrix)
library(rjags)
library(gemtc)

source("scenario_setup.R")
source("data_generation.R")
source("methods.R")
source("simulation_run.R")

# --- 1. Covariate definitions ------------------------------------------
covdefs <- define_covariates()

# --- 2. Scenario grid --------------------------------------------------
# All factor combinations are defined by the defaults in scenario_grid().
# Only network_type is varied explicitly here.

scens <- list()
for (nw in c("basic_3", "star_4", "one_loop_closed")) {
  scens <- c(scens, scenario_grid(network_type = nw))
}

# --- 3. Run simulation -------------------------------------------------

res <- run_simulations_min(
  scenarios      = scens,
  R              = 1000,
  include_bucher = TRUE,
  include_bayes  = TRUE,
  jags_ctrl      = list(n.chains=2, n.adapt=1000, n.iter=3000, thin=3, seed=1000),
  base_seed      = 1000,
  covdefs        = covdefs,
  parallel       = TRUE
)

# --- 4. Performance summary --------------------------------------------

perf <- summarize_performance_min(res)
as.data.frame(perf)
