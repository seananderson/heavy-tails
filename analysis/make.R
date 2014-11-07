# All steps to run the analysis
# Created in Vim with :read !ls *.R

source("0-make-data.R") # takes a few minutes
# Begin: on the westgrid server
source("1-compile-models.R") # warning: takes a long time
source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")
source("2-run-model-logistic.R") # warning: takes a long time
source("2-run-models-ar1.R") # warning: takes a long time
source("2-run-models-base-stronger-prior.R") # warning: takes a long time
source("2-run-models-base-weaker-prior.R") # warning: takes a long time
source("2-run-models-base.R") # warning: takes a long time
source("2-run-models-obs-0.2.R") # warning: takes a long time
source("2-run-models-rate.R") # warning: takes a long time
source("2.1-get-base-mean-sd.R")
source("3.0-test-t-sampling.R") # warning: takes a long time
source("3.1-plot-sample-test.R")
source("3.2-get-effect-nu-seeds.R") # warning: takes a long time
source("3.3-test-gomp-models.R") # warning: takes a long time
# End: on the westgrid server
source("3.4-plot-test-gomp-models.R")
source("5-shape-data.R")
# Begin: on the westgrid server
source("5.8-stan-beta-modelling.R") # warning: takes a long time; caching implemented
source("5.9-order-level-posteriors.R") # warning: takes a long time; caching implemented
# End: on the westgrid server
source("6-plot-alt-models.R")
source("6-plot-correlates.R") # must run 5.8.. first
source("6-plot-eg-ts-gpdd.R")
# source("6-plot-gomp-samples.R") # old
source("6-plot-nu-coefs.R")
source("6-plot-order-correlate-posteriors.R") # must run 5.9... first
source("6-plot-prior.R")
source("6-plot-sparks.R")
source("6-plot-t-nu-eg.R")
source("7-values-for-paper.R")
