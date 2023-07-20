# EM algorithm
# Mixture of ZIP
# No covariates

u.mzip <- function(Y, K, phi, prob, lambda) {

	u <- array(dim = c(dim(Y), K))

	for (k in 1:K) {
		num <- phi[k]
		den <- (phi[k] + (1 - phi[k]) * exp(-lambda[k, ]))

		uk <- num / den
		uk <- matrix(uk, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
		uk[Y != 0] <- 0

		u[,, k] <- uk
	}

	return(u)
}


EMmzip <- function(Y, K, phi, prob, lambda, tol = 1e-8, maxiter = 1000) {

    N <- nrow(Y)
	G <- ncol(Y)

	phi.h    <- phi
	lambda.h <- lambda
	pi.h     <- prob

	ll   <- llmzip(Y, K, phi.h, lambda.h, pi.h)
	ll.0 <- ll$loglik

	dl   <- Inf
	iter <- 0

	while (dl > tol & iter < maxiter) {

		# E step
		z.hat <- ll$z
		u.hat <- u.mzip(Y, K, phi.h, pi.h, lambda.h)

		if (any(colSums(z.hat) == 0)) {
			z.hat <- z.hat + (1 / N)
		}

		# M step
		pi.h <- colSums(z.hat)
		pi.h <- pi.h / sum(pi.h)

		phi.num <- numeric(K)
		for (k in 1:K) {
			phi.num[k] <- sum(z.hat[, k] * u.hat[,, k])
		}
		phi.den <- G * colSums(z.hat)
		phi.h   <- phi.num / phi.den

		lambda.num <- matrix(nrow = K, ncol = G)
		lambda.den <- matrix(nrow = K, ncol = G)
		for (k in 1:K) {
			lambda.num[k, ] <- colSums(z.hat[, k] * (1 - u.hat[,, k]) * Y)
			lambda.den[k, ] <- colSums(z.hat[, k] * (1 - u.hat[,, k]))
		}
		lambda.h <- lambda.num / lambda.den

		# updated log-likelihood
		ll    <- llmzip(Y, K, phi.h, lambda.h, pi.h, ll$trick, ll$thresh)
		ll.n  <- ll$loglik

		dl   <- abs(ll.n - ll.0)
		ll.0 <- ll.n
		iter <- iter + 1
	}

	out <- list(
		prob   = pi.h,
		phi    = phi.h,
		lambda = lambda.h,
		U      = u.hat,
		Z      = z.hat,
		loglik = ll.0,
		iter   = iter
	)

	return(out)
}
