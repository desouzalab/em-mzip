u.mzip <- function(Y, K, phi, prob, lambda) {

	num <- sum(prob * phi)

	den <- 0
	for (k in 1:K) {
		den <- den + prob[k] * (phi[k] + (1 - phi[k]) * exp(-lambda[k, ]))
	}

	u <- num / den
	u <- matrix(u, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
	u[Y != 0] <- 0

	return(u)
}


EMmzip <- function(Y, K, phi, prob, lambda, tol = 1e-8, maxiter = 1000) {

	G <- ncol(Y)

	K <- as.integer(K)
	if (K <= 0) {
		stop("K must be an integer greater than or equal to 1")
	}

	if (K != length(phi)) {
		phi <- rep(phi, length.out = K)
	}
	if (K != length(prob)) {
		prob <- rep(prob, length.out = K)
	}

	phi.h    <- phi
	lambda.h <- lambda
	pi.h     <- prob

	ll    <- llmzip(Y, K, phi.h, lambda.h, pi.h)
	ll.0  <- ll$loglik

	dl   <- 1 + tol
	iter <- 0

	while (dl > tol & iter < maxiter) {
	    # cat("Iteration", iter, "\n")
	    # cat("Log-lik", ll.0, "\n")

		iter <- iter + 1

		# E step
		z.hat <- ll$z
		u.hat <- u.mzip(Y, K, phi.h, pi.h, lambda.h)

		# M step
		pi.h <- colSums(z.hat)
		pi.h <- pi.h / sum(pi.h)

		phi.num <- numeric(K)
		for (k in 1:K) {
			phi.num[k] <- sum(z.hat[, k] * u.hat)
		}
		phi.den <- G * colSums(z.hat)
		phi.h   <- phi.num / phi.den

		lambda.num <- matrix(nrow = K, ncol = G)
		lambda.den <- matrix(nrow = K, ncol = G)
		for (k in 1:K) {
			lambda.num[k, ] <- colSums(z.hat[, k] * (1 - u.hat) * Y)
			lambda.den[k, ] <- colSums(z.hat[, k] * (1 - u.hat))
		}
		lambda.h <- lambda.num / lambda.den

		# log-likelihood
		ll    <- llmzip(Y, K, phi.h, lambda.h, pi.h, ll$trick)
		ll.n  <- ll$loglik
		trick <- ll$trick

		dl    <- abs(ll.n - ll.0)
		ll.0  <- ll.n
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
