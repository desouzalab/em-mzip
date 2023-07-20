args <- as.numeric(commandArgs(trailingOnly = TRUE))

ver  <- args[1]
mod  <- args[2]
sval <- args[3]

source("../data_fn.R")
data_fn("Klein", ver, mod, sval, S = 32, seed = 5885598)
