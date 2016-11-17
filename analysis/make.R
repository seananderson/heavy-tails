# All steps to run the analysis
# =============================
#
# - Alternatively, run `make` from the Terminal
# - Run `make clean` in the Terminal to remove all cached `.rds` files
# - Note that re-running all the model fitting will take a long time (days)
# - If you have cloned the repository from GitHub, then most of the cached
#   files will already be available.

# Start in "heavy-tails" folder, then:
setwd("analysis")

if(!file.exists("gpdd-clean.rds"))
  source("0-make-data.R") # takes a few minutes
# Compile Stan models:
stan_models <- c("stan-gomp.rds", "stan-gomp-bda.rds",
  "stan-gomp-obs.rds", "stan-t.rds", "stan-logistic.rds",
  "stan-gomp2-ar1.rds", "stan-rate.rds", "stan-gomp-uniform.rds",
  "stan-rw.rds", "stan-gomp-gamma.rds")
if(any(!file.exists(stan_models)))
  source("1-compile-models.R")
# Functions used in the fitting:
source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")
# Run main models:
# warning: these take a long time to run (hours to a day each)
if(!file.exists("logistic-hat.rds"))
  source("2-run-model-logistic.R")
if(!file.exists("gomp-ar1-hat.rds"))
  source("2-run-models-ar1.R")
if(!file.exists("gomp-base-stronger-hat.rds"))
  source("2-run-models-base-stronger-prior.R")
if(!file.exists("gomp-base-weaker-hat.rds"))
  source("2-run-models-base-weaker-prior.R")
if(!file.exists("gomp-base-hat.rds"))
  source("2-run-models-base.R")
if(!file.exists("gomp-obs-0.2-hat.rds"))
  source("2-run-models-obs-0.2.R")
if(!file.exists("rate-hat.rds"))
  source("2-run-models-rate.R")
if(!file.exists("gomp-base-mean-sd.rds"))
  source("2.1-get-base-mean-sd.R")
# TODO FINISH THIS
# Simulation testing
if(!file.exists("sample-t-sim-check.rds"))
  source("3.0-test-t-sampling.R") # warning: takes a long time
source("3.1-plot-sample-test.R")
if(!file.exists("nu_effective_seeds.rda"))
  source("3.2-get-effect-nu-seeds.R") # warning: takes a long time
if(!file.exists("check_nu.rds"))
  source("3.3-test-gomp-models.R") # warning: takes a long time
source("3.4-plot-test-gomp-models.R")
source("5-shape-data.R")
# warning: takes a long time; caching implemented
source("5.8-stan-beta-modelling.R")
# warning: takes a long time; caching implemented
source("5.9-order-level-posteriors.R")
source("5.10-extract-skew-samples.R")
source("6-plot-alt-models.R")
source("6-plot-correlates.R") # must run 5.8.. first
source("6-plot-eg-ts-gpdd.R")
source("6-plot-nu-coefs.R")
source("6-plot-order-correlate-posteriors.R") # must run 5.9... first
source("6-plot-prior.R")
source("6-plot-sparks.R")
source("6-plot-t-nu-eg.R")
source("6-causes-table.R")
source("6-plot-skewness.R")
source("8.1-plot-skewness.R")
source("9-values-for-paper.R")

# Created in Vim with :read !ls *.R
