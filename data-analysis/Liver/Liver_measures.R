path <- "~/scratch/Liver/"
versions <- list.files(path)
dir.create("RData", showWarnings = FALSE)

N.list <- rep(c(4685, 0, 1000), c(7, 13, 3))
G.list <- c(1000, 700, 500, 438, 1001, 700, 500, rep(0, 13), 100, 100, 116)

versions <- c("v21m3s1", "v21m3s2", "v22m3s1", "v22m3s2", "v23m3s1", "v23m3s2")

for (v in versions) {
    cat("Loading", v, "\n")
    vdir  <- paste0(path, v)
    cases <- list.files(vdir)

    ver <- as.numeric(substr(v, 2, 3))

    N <- N.list[ver]
    G <- G.list[ver]

    source("../data_measures.R")
    save(smin, zhat_mat, mat, param_mat, file = sprintf("RData/%s.RData", v))
}
