em_fun <- function(i, ipath, G, K, phi, prob, beta0, rho, init) {
	Y <- scan(paste0(ipath, "/Data ", i, " Y.dat"), quiet = TRUE)
	Y <- matrix(Y, ncol = G, byrow = TRUE)

	Sn <- scan(paste0(ipath, "/Data ", i, " Sn.dat"), quiet = TRUE)

	if (!is.na(init) & init == "k-means") {
		Z <- kmeans(Y, centers = K)$cluster

		prob <- as.numeric(table(Z))
		prob <- prob / sum(prob)

		phi <- tapply(1:N, Z, function(x) {mean(Y[x, ] == 0)})

		rho <- tapply(1:N, Z, function(x) {colMeans(Y[x, ])})
		rho <- t(simplify2array(rho))
		beta0 <- colMeans(rho)
		rho <- rho - rep(beta0, each = K)
	}
	elapsed <- system.time({
		out <- EMmgzip(Y, Sn, K, phi, prob, beta0, rho)
	})["elapsed"]

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

mgzip_estimation <- function(main_dir = "Untitled", parallel = TRUE, init = "true") {

	main_dir <- paste0("~/scratch/Data/", main_dir)
	cat(sprintf("Estimating %s\n", main_dir))

	if (parallel) {
		library(parallel)
		ncl <- detectCores()
		cat(sprintf("\tUsing %d cores\n", ncl))
	}


	cases <- list.files(main_dir)
	init <- rep(init, length.out = length(cases))

	for (i in seq_along(cases)) {
		case <- cases[i]
		case_dir <- paste0(main_dir, "/", case)
		if (!is.na(init[i]) & init[i] == "k-means") {
			cat("k-means method of initializing\n")
		}

		t0 <- Sys.time()
		source(paste0(case_dir, "/_Parameters.txt"))

		n <- list.files(case_dir)
		n <- regmatches(n, regexec("Data ([0-9]+) Y.dat", n))
		n <- as.numeric(sapply(n, `[`, 2))
		n <- sort(n[!is.na(n)])

		if (parallel) {
			cl <- makeForkCluster(ncl)
			times <- parSapply(cl, n, em_fun, case_dir, G, K, phi, prob, beta0, rho, init[i])
			stopCluster(cl)
		} else {
			times <- sapply(n, em_fun, case_dir, G, K, phi, prob, beta0, rho, init[i])
		}
		cat(times, file = paste0(case_dir, "/_times.dat"), sep = "\n")
		t1 <- Sys.time()

		cat(sprintf("\t%s, %d estimations, %s\n",
			case, length(n), format(round(t1 - t0, 2))))
	}
}
