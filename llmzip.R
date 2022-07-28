llmzip_true <- function(Y, K, phi, lambda, prob) {

	bmat <- matrix(nrow = nrow(Y), ncol = K)
	for (k in 1:K) {
		bk <- dzip(Y, phi[k], rep(lambda[k, ], each = nrow(Y)), log = TRUE)
		bk <- rowSums(bk)

		bk <- log(prob[k]) + bk
		bmat[, k] <- exp(bk)
	}

	z  <- bmat / rowSums(bmat)
	ll <- sum(log(rowSums(bmat)))

	return(list(
		loglik = ll,
		z.hat  = z,
		trick  = FALSE
	))
}


llmzip_trick <- function(Y, K, phi, lambda, prob) {
	eps <- log(.Machine$double.xmin)

	bmat <- matrix(nrow = nrow(Y), ncol = K)
	for (k in 1:K) {
		bk <- dzip(Y, phi[k], rep(lambda[k, ], each = nrow(Y)), log = TRUE)
		bk <- rowSums(bk)
		bk[bk == -Inf] <- min(bk[bk > -Inf])

		bk <- log(prob[k]) + bk
		bmat[, k] <- exp(bk / (min(bk) / eps))
	}

	z  <- bmat / rowSums(bmat)
	ll <- sum(log(rowSums(bmat)))

	return(list(
		loglik = ll,
		z.hat  = z,
		trick  = TRUE
	))
}


llmzip <- function(Y, K, phi, lambda, prob, trick = FALSE) {

	if (trick) {
		return(llmzip_trick(Y, K, phi, lambda, prob))

	} else {
		out <- llmzip_true(Y, K, phi, lambda, prob)

		if (any(is.infinite(out$z.hat), is.nan(out$z.hat),
				is.infinite(out$loglik), is.nan(out$loglik))) {
		    warning("Proportial log-likelihhod function")
			return(llmzip(Y, K, phi, lambda, prob, TRUE))
		}

		return(out)
	}
}
