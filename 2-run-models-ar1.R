# This file runs the Gompertz models with an AR1 coefficient on the residuals

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-gomp2-ar1.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-ar1"

fit_gpdd_model(gpdd_dat = gpdd, model = model, sub_folder = id)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = TRUE, type = "gompertz")

saveRDS(out, file = paste0(id, "-hat.rds"))
