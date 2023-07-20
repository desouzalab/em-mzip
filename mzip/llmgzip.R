llmgzip_true <- function(Y, K, phi, prob, lambda, thresh) {

	bmat <- matrix(nrow = nrow(Y), ncol = K)

	for (k in 1:K) {
		bk <- dzip(Y, phi[k], lambda[, , k], log = TRUE)
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


llmgzip_trick <- function(Y, K, phi, prob, lambda, thresh) {
	eps <- log(.Machine$double.xmin)

	bmat <- matrix(nrow = nrow(Y), ncol = K)
	for (k in 1:K) {
		bk <- dzip(Y, phi[k], lambda[, , k], log = TRUE)
		bk <- rowSums(bk)

		minbk <- min(bk[bk > -Inf])
		thresh[k] <- ifelse(thresh[k] < minbk & minbk < Inf, minbk, thresh[k])
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


llmgzip <- function(Y, K, phi, prob, lambda,
                    trick = FALSE, thresh = rep(-Inf, K)) {

	if (trick) {
		out <- llmgzip_trick(Y, K, phi, prob, lambda, thresh)

	} else {
		out <- llmgzip_true(Y, K, phi, prob, lambda, thresh)

		if (any(is.infinite(out$z.hat), is.nan(out$z.hat),
				is.infinite(out$loglik), is.nan(out$loglik))) {
			warning("Proportional log-likelihood function")
			out <- llmgzip_trick(Y, K, phi, prob, lambda, thresh)
		}
	}

	return(out)
}
