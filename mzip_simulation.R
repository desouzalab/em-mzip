mzip_simulation <- function(N, G, K, phi, prob, lambda,
							S = 10,
							seed = as.integer(Sys.time()),
							main_dir = NULL,
							parallel = TRUE) {
	set.seed(seed)
	if (parallel) {
		library(parallel)
		ncl <- detectCores()
		cat(sprintf("Using %d cores\n", ncl))
	}

	main_dir <- ifelse(is.null(main_dir), "Data", paste0("Data/", main_dir))

	if (!dir.exists(main_dir)) {
		if (!dir.exists("Data")) {
			cat("Creating directory Data\n")
			dir.create("Data")
		}
		cat("Creating directory", main_dir, "\n")
		dir.create(main_dir)

		n <- 1

	} else {
		n <- list.files(main_dir)

		if (length(n) == 0) {
			n <- 1
		} else {
			n <- list.files(main_dir)
			n <- regmatches(n, regexec("Case ([0-9]+)", n))
			n <- as.numeric(sapply(n, `[`, 2))
			n <- max(n, na.rm = TRUE) + 1
		}
	}

	file_dir <- paste0(main_dir, "/Case ", n)
	cat("Creating directory", file_dir, "\n")
	dir.create(file_dir)


	sink(paste0(file_dir, "/_Parameters.txt"))
	cat("S      <- ", S, "\n", sep = "")
	cat("seed   <- ", seed, "\n", sep = "")
	cat("\n")
	cat("N      <- ", N, "\n", sep = "")
	cat("G      <- ", G, "\n", sep = "")
	cat("K      <- ", K, "\n", sep = "")
	cat("phi    <- c(", paste(phi, collapse = ", "), ")\n", sep = "")
	cat("prob   <- c(", paste(prob, collapse = ", "), ")\n", sep = "")
	cat("lambda <- matrix(c(", paste(t(lambda), collapse = ", "),
		"),\n\tK, G, byrow = TRUE)\n", sep = "")
	sink()

	if (parallel) {
		simfun <- function(j) {
			files <- paste0(file_dir, "/Data ", j, " ", c("Y", "Z", "U"), ".dat")

			Y <- rmzip(N, G, K, phi, prob, lambda)
			write(t(Y$Y), file = files[1], ncolumns = G)
			write(t(Y$Z), file = files[2], ncolumns = 1)
			write(t(Y$U), file = files[3], ncolumns = G)
		}

		t0 <- Sys.time()
		cl <- makeForkCluster(ncl)
		clusterSetRNGStream(cl, seed)
		parSapply(cl, seq_len(S), simfun)
		stopCluster(cl)
		t1 <- Sys.time()

	} else{
		t0 <- Sys.time()
		for (j in 1:S) {
			files <- paste0(file_dir, "/Data ", j, " ", c("Y", "Z", "U"), ".dat")

			Y <- rmzip(N, G, K, phi, prob, lambda)
			write(t(Y$Y), file = files[1], ncolumns = G)
			write(t(Y$Z), file = files[2], ncolumns = 1)
			write(t(Y$U), file = files[3], ncolumns = G)
		}
		t1 <- Sys.time()
	}


	cat("\t", S, "simulations saved,", format(round(t1 - t0, 2)), "\n")
}
