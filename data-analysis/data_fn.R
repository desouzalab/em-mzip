data_fn <- function(Data, ver, mod, sval, S = 1, K = 2:9, seed = 0) {

    seed <- seed + mod * 1000 + ver * 100 + sval * 10

    if (mod == 1) {
        source("../../mzip/dzip.R")
        source("../../mzip/llmzip.R")
        source("../../mzip/rmzip.R")
        source("../../mzip/EMmzip.R")
    } else if (mod == 2) {
        source("../../mzip/dzip.R")
        source("../../mzip/llmgzip.R")
        source("../../mzip/rmgzip.R")
        source("../../mzip/EMmgzip.R")
        source("../../mzip/NRmgzip.R")
    } else if (mod == 3) {
        source("../../mzinb/ZINB_disp.R")
    } else if (mod == 4) {
        source("../../mzinb/ZINB_forg.R")
    } else if (mod == 5) {
        source("../../mzinb/ZINB_full.R")
    }

    source("../data_estim.R")

    path <- sprintf("%s/v%dm%ds%d", Data, ver, mod, sval)

    dir.create(       "~/scratch",         showWarnings = FALSE)
    dir.create(paste0("~/scratch/", Data), showWarnings = FALSE)
    dir.create(paste0("~/scratch/", path), showWarnings = FALSE)

    load(sprintf("Data/%s v%d.RData", Data, ver))

	init <- c("random", "kmeans")[sval]

	for (k in K) {
		data_estim(mod, Y, Tn, k, S, path, seed + k, init)
	}
}
