path <- "~/scratch/Klein/"
versions <- list.files(path)
dir.create("RData", showWarnings = FALSE)

for (v in versions) {
    cat("Loading", v, "\n")
    vdir  <- paste0(path, v)
    cases <- list.files(vdir)

    ver <- as.numeric(substr(v, 2, 2))

    N <- 1616
    G <- c(5000, 4514, 500, 500, 100, 100)[ver]

    source("../data_measures.R")
    save(smin, zhat_mat, mat, param_mat, file = sprintf("RData/%s.RData", v))
}
