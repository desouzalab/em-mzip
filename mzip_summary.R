mzip_summary <- function(main_dir = NULL, estim = c("iter", "phi",
						 "prob", "times", "lambda")) {

	main_dir <- ifelse(is.null(main_dir), "Data", paste0("Data/", main_dir))

	times <- "times" %in% estim
	estim <- estim[estim %in% c("iter", "loglik", "phi", "prob", "lambda")]

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

		if ("lambda" %in% estim) {
			vlambda <- as.numeric(t(lambda))
			k <- rep(1:nrow(lambda), each = ncol(lambda))
			vl <- sort(unique(vlambda))
			mse <- rowMeans((t(mat$lambda[[case]]) - vlambda) ^ 2)
			mse <- tapply(mse, k, mean)

			mat$lambda[[case]] <- mse
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

mzip_summary_Z <- function(main_dir = NULL) {
	# clevr
	library(FDRSeg)
	# main_dir = "Scenario 1"

	main_dir <- ifelse(is.null(main_dir), "Data", paste0("Data/", main_dir))

	Vmat <- list()
	cases <- list.files(main_dir)

	for (case in cases) {
		# case = cases[1]
		case_dir <- paste0(main_dir, "/", case)
		# source(paste0(case_dir, "/_Parameters.txt"))

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

