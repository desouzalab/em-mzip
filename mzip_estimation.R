mzip_estimation <- function(main_dir = NULL, parallel = TRUE) {

	main_dir <- ifelse(is.null(main_dir), "Data", paste0("Data/", main_dir))

	if (parallel) {
		library(parallel)
		ncl <- detectCores()
		cat(sprintf("Using %d cores\n", ncl))
	}

	em_fun <- function(i, ipath, G, ...) {
		Y <- scan(paste0(ipath, "/Data ", i, " Y.dat"), quiet = TRUE)
		Y <- matrix(Y, ncol = G, byrow = TRUE)

		elapsed <- system.time(
			out <- EMmzip(Y, ...)
		)["elapsed"]

		for (o in names(out)) {
			estim <- out[[o]]
			if (is.matrix(estim)) {
				write(t(estim),
					  file = paste0(ipath, "/EM ", i, " ", o, ".dat"),
					  ncolumns = ncol(estim))
			} else {
				cat(estim, file = paste0(ipath, "/EM ", i, " ", o, ".dat"))
			}
		}

		return(elapsed)
	}

	cases <- list.files(main_dir)
	for (case in cases) {
		case_dir <- paste0(main_dir, "/", case)

		t0 <- Sys.time()
		source(paste0(case_dir, "/_Parameters.txt"))

		n <- list.files(case_dir)
		n <- regmatches(n, regexec("Data ([0-9]+) Y.dat", n))
		n <- as.numeric(sapply(n, `[`, 2))
		n <- sort(n[!is.na(n)])

		if (parallel) {
			# parallel computing
			cl <- makeForkCluster(ncl)
			times <- parSapply(cl, n, em_fun, case_dir, G, K, phi, prob, lambda)
			stopCluster(cl)

			cat(times, file = paste0(case_dir, "/_times.dat"), sep = "\n")

		} else {
			# single-thread computing
			times <- sapply(n, em_fun, case_dir, G, K, phi, prob, lambda)
			cat(times, file = paste0(case_dir, "/_times.dat"), sep = "\n")

		}
		t1 <- Sys.time()

		cat("\t", length(n), "estimations,", format(round(t1 - t0, 2)), "\n")
	}
}
