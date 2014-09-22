gpdd <- readRDS("gpdd-clean.rds")
ids <- unique(gpdd$main_id)
main_id_vec <- ids[201:400]
source("2-run-models-ar1.R")
