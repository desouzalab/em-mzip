mgzip_summary <- function(main_dir = NULL,
	estim = c("times", "iter", "phi", "prob", "rho", "beta0")) {

	if (is.null(main_dir)) {
		main_dir <- "Data"
	}

	times <- "times" %in% estim
	estim <- estim[estim %in% c("iter", "loglik", "phi", "prob", "rho", "beta0")]

	mat <- list()
	cases <- list.files(main_dir)

	for (case in cases) {
		case_dir <- paste0(main_dir, "/", case)
		source(paste0(case_dir, "/_Parameters.txt"))

		n <- list.files(case_dir)
		n <- regmatches(n, regexec("Data ([0-9]+) Y.dat", n))
		n <- as.numeric(sapply(n, `[`, 2))
		n <- sort(n[!is.na(n)])

		if (times) {
			mat$times[[case]] <-
				t(scan(paste0(case_dir, "/_times.dat"), quiet = TRUE))
		}

		for (e in estim) {
			mat[[e]][[case]] <- t(sapply(n, function(x) {
				scan(paste0(case_dir, "/EM ", x, " ", e, ".dat"), quiet = TRUE)
			}))
		}

		if ("rho" %in% estim) {
			vrho <- as.numeric(t(rho))
			k <- rep(1:nrow(rho), each = ncol(rho))
			mse <- rowMeans((t(mat$rho[[case]]) - vrho) ^ 2)
			mse <- tapply(mse, k, mean)

			mad <- apply(abs(t(mat$rho[[case]]) - vrho), 1, median)
			mad <- tapply(mad, k, median)

			mat$rho[[case]] <- mse
			mat$rho_mad[[case]] <- mad
		}

		if ("beta0" %in% estim) {
			vbeta0 <- as.numeric(t(beta0))
			vb <- sort(unique(vbeta0))
			k <- sapply(beta0, function(x)which(vb == x))

			mse <- colMeans((mat$beta0[[case]] - vbeta0) ^ 2)
			mse <- tapply(mse, k, mean)
			names(mse) <- vb

			mad <- apply(abs(mat$beta0[[case]] - vbeta0), 2, median)
			mad <- tapply(mad, k, median)
			names(mad) <- vb

			mat$beta0[[case]] <- mse
			mat$beta0_mad[[case]] <- mad
		}

	}

	if (times) {
		estim <- c(estim, "times")
	}

	for (e in estim) {
		j <- which.max(sapply(mat[[e]], length))
		nj <- length(mat[[e]][[j]])
		aux <- NULL

		for (case in cases) {
			aux <- c(aux,
					 mat[[e]][[case]],
					 rep(NA, nj - length(mat[[e]][[case]])))
		}

		mat[[e]] <- array(aux, c(dim(mat[[e]][[j]]), length(cases)))

		if (length(dim(mat[[e]])) == 3) {
			if (dim(mat[[e]])[1] == 1) {
				mat[[e]] <- aperm(mat[[e]], c(2, 3, 1))[,, 1]
			}  else {
				mat[[e]] <- aperm(mat[[e]], c(1, 3, 2))
			}
		}
	}

	return(mat)
}


mgzip_summary_Z <- function(main_dir = NULL) {
	if (is.null(main_dir)) {
		main_dir <- "Data"
	}

	Vmat <- list()
	cases <- list.files(main_dir)

	for (case in cases) {
		case_dir <- paste0(main_dir, "/", case)

		n <- list.files(case_dir)
		n <- regmatches(n, regexec("Data ([0-9]+) Y.dat", n))
		n <- as.numeric(sapply(n, `[`, 2))
		n <- sort(n[!is.na(n)])

		aux <- function(x) {
			z0 <- scan(paste0(case_dir, "/Data ", x, " Z.dat"), quiet = TRUE)
			zhat <- scan(paste0(case_dir, "/EM ", x, " Z.dat"), quiet = TRUE)
			zhat <- matrix(zhat, nrow = length(z0), byrow = TRUE)
			zhat <- apply(zhat, 1, which.max)
			return(v_measure(z0, zhat))
		}

		Vmat[[case]] <- sapply(n, aux)
	}


	j <- which.max(sapply(Vmat, length))
	nj <- length(Vmat[[j]])
	aux <- NULL

	for (case in cases) {
		aux <- c(aux,
				 Vmat[[case]],
				 rep(NA, nj - length(Vmat[[case]])))
	}

	Vmat <- array(aux, c(length(Vmat[[j]]), length(cases)))

	return(Vmat)
}
