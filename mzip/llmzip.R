llmzip_true <- function(Y, K, phi, lambda, prob, thresh) {

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
		trick  = FALSE,
	    thresh = thresh
	))
}


llmzip_trick <- function(Y, K, phi, lambda, prob, thresh) {
	eps <- log(.Machine$double.xmin)

	bmat <- matrix(nrow = nrow(Y), ncol = K)
	for (k in 1:K) {
		bk <- dzip(Y, phi[k], rep(lambda[k, ], each = nrow(Y)), log = TRUE)
		bk <- rowSums(bk)

		thresh[k] <- max(thresh[k], min(bk[bk > -Inf]))
		bk[bk < thresh[k]] <- thresh[k]

		bk <- log(prob[k]) + bk
		bmat[, k] <- exp(bk / (min(bk) / eps))
	}

	z  <- bmat / rowSums(bmat)
	ll <- sum(log(rowSums(bmat)))

	return(list(
		loglik = ll,
		z.hat  = z,
		trick  = TRUE,
		thresh = thresh
	))
}


llmzip <- function(Y, K, phi, lambda, prob,
                   trick = FALSE, thresh = rep(-Inf, K)) {

	if (trick) {
		out <- llmzip_trick(Y, K, phi, lambda, prob, thresh)

	} else {
		out <- llmzip_true(Y, K, phi, lambda, prob, thresh)

		if (any(is.infinite(out$z.hat), is.nan(out$z.hat),
				is.infinite(out$loglik), is.nan(out$loglik))) {
		    warning("Proportional log-likelihood function")
		    out <- llmzip_trick(Y, K, phi, lambda, prob, thresh)
		}
	}

	return(out)
}
