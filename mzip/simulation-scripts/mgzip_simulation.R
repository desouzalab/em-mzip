mgzip_simulation <- function(N, G, K, phi, prob, beta0, rho,
							 mu, sigma, S = 1,
							 seed = as.integer(Sys.time()),
							 main_dir = "Untitled",
							 parallel = TRUE) {


	main_dir <- paste0("~/scratch/Data/", main_dir)
	if (!dir.exists(main_dir)) {
		cat("Creating directory", main_dir, "\n")
		dir.create(main_dir)
	}

	n <- list.files(main_dir)

	if (length(n) == 0) {
		n <- 1
	} else {
		n <- list.files(main_dir)
		n <- regmatches(n, regexec("Case ([0-9]+)", n))
		n <- as.numeric(sapply(n, `[`, 2))
		n <- max(n, na.rm = TRUE) + 1
	}

	file_dir <- paste0(main_dir, "/Case ", n)
	cat("Creating directory", file_dir, "\n")
	dir.create(file_dir)

	{
		sink(paste0(file_dir, "/_Parameters.txt"))
		cat("# Mixture of generalized ZIP model\n")
		cat("S     <- ", S, "\n", sep = "")
		cat("seed  <- ", seed, "\n", sep = "")
		cat("\n")
		cat("N     <- ", N, "\n", sep = "")
		cat("G     <- ", G, "\n", sep = "")
		cat("K     <- ", K, "\n", sep = "")
		cat("phi   <- c(", paste(phi, collapse = ", "), ")\n", sep = "")
		cat("prob  <- c(", paste(prob, collapse = ", "), ")\n", sep = "")
		cat("rho   <- matrix(c(", paste(t(rho), collapse = ", "),
			"),\n\tK, G, byrow = TRUE)\n", sep = "")
		cat("beta0 <- c(", paste(beta0, collapse = ", "), ")\n", sep = "")
		sink()
	}

	simfun <- function(j) {
		files <- paste0(file_dir, "/Data ", j, " ", c("Y", "Z", "U", "Sn"), ".dat")

		data <- rmgzip(N, G, K, phi, prob, beta0, rho, mu, sigma)
		write(t(data$Y),  file = files[1], ncolumns = G)
		write(t(data$Z),  file = files[2], ncolumns = 1)
		write(t(data$U),  file = files[3], ncolumns = G)
		write(t(data$Sn), file = files[4], ncolumns = 1)
	}

	t0 <- Sys.time()

	if (parallel) {
		library(parallel)
		ncl <- detectCores()
		cat(sprintf("\tUsing %d cores\n", ncl))

		cl <- makeForkCluster(ncl)
		clusterSetRNGStream(cl, seed)
		parSapply(cl, seq_len(S), simfun)
		stopCluster(cl)
	} else{
		set.seed(seed)
		sapply(seq_len(S), simfun)
	}
	t1 <- Sys.time()

	cat(sprintf("\t%d simulations saved,%s\n",
		S, format(round(t1 - t0, 2))))
}
