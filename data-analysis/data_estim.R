data_estim <- function(mod, Y, Sn, K, S, path, seed,
                       init = "random", parallel = TRUE) {

    path <- sprintf("~/scratch/%s", path)
    if (!dir.exists(path)) {
        cat("Creating directory", path, "\n")
        dir.create(path)
    }

    path <- sprintf("%s/%d clusters", path, K)
    if (dir.exists(path)) {
        files <- list.files(path)
        if ("_times.dat" %in% files) {
            cat("Directory", path, " already exists\n")
            return(NULL)
        } else {
            cat("Cleaning directory", path, "\n")
            file.remove(paste0(path, "/", files))
        }
    } else {
        cat("Creating directory", path, "\n")
        dir.create(path)
    }

    t0 <- Sys.time()

    if (parallel) {
        library(parallel)
        ncl <- min(detectCores(), S)
        cat(sprintf("Using %d cores\n", ncl))

        cl <- makeForkCluster(ncl)
        times <- parSapply(cl, seq_len(S), data_estim_aux,
                           Y, Sn, K, path, init, seed, S)

        stopCluster(cl)
    } else {
		times <- sapply(seq_len(S), data_estim_aux,
                        Y, Sn, K, path, init, seed, S)
    }

    cat(times, file = sprintf("%s/_times.dat", path), sep = "\n")
    t1 <- Sys.time()

    cat("\t", S, "estimations,", format(round(t1 - t0, 2)), "\n")
    return(1)
}


data_estim_aux <- function(s, Y, Sn, K, path, init, seed, S) {

	set.seed(seed + s)
    N <- nrow(Y)
    G <- ncol(Y)

    if (init == "kmeans") {
        Z <- kmeans(Y, centers = K)$cluster
    } else {
        Z <- sample(K, N, replace = TRUE)
    }
    cat(Z, file = paste0(path, "/Z0 ", s, ".dat"))

	beta0 <- rep(0, G)
    z0    <- theta0(Y, Z, K)

    elapsed <- tryCatch({system.time(
		if (mod == 1) {
			out <- EMmzip(Y, K, z0$phi, z0$pi, z0$lambda)
		} else if (mod == 2) {
			out <- EMmgzip(Y, Sn, K, z0$phi, z0$pi, beta0, z0$lambda)
		} else if (mod == 3) {
			out <- ZINB_mix_K_comp(Y, z0$pi, z0$phi, z0$lambda, z0$size,
			    eps = 1e-2)
		} else if (mod %in% c(4, 5)) {
			out <- ZINB_Mix_EM(Y, Sn, beta0, z0$lambda, z0$size, z0$pi, z0$phi,
			    eps = 1e-2)
		}
    	)["elapsed"]},
    error = function(e) {
    	message('An error occurred while executing')
    	return(NA)
    })

    if (!is.na(elapsed)) {
        for (o in names(out)) {
            x <- out[[o]]
            if (length(dim(x)) == 2) {
                write(t(x),
                    file = sprintf("%s/EM %d %s.dat", path, s, o),
                    ncolumns = ifelse(is.matrix(x), ncol(x), length(x)))
            } else {
                write(x,
                    file = sprintf("%s/EM %d %s.dat", path, s, o),
                    ncolumns = ifelse(is.matrix(x), ncol(x), length(x)))
            }
        }
    }

    while (is.na(elapsed)) {
		elapsed <- data_estim_aux(s, Y, Sn, K, path, init, seed + S, S)
	}

    return(elapsed)
}


theta0 <- function(Y, Z, K) {

	# proportions
	prob <- as.numeric(table(Z))
    prob <- prob / sum(prob)

	mat <- t(sapply(1:K, function(k) {theta0.aux(Y[Z == k, ])}))

	return(list(
		pi     = prob,
		phi    = mat[, 2],
		lambda = mat[, -c(1, 2)],
		size   = mat[, 1]
	))
}


theta0.aux <- function(x) {
	if (!is.matrix(x)) {
		x <- t(as.matrix(x))
	}

	# raw proportions of zeros
	phi <- apply(x, 2, function(xg){mean(xg == 0)})
	phi[phi < .01] <- .01
	phi[phi > .99] <- .99
	phi <- mean(phi, na.rm = TRUE)

	# raw mean
	lambda <- colMeans(x)

	# alpha = size = 1/dispersion
	mhat <- mean(x)
	shat <- sd(x)
	size <- mhat^2 / (shat^2 - mhat)

	return(c(size, phi, lambda))
}
