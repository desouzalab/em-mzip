# EM algorithm
# Mixture of ZIP
# Newton-Raphson with full Hessian
# Size factors

u.mgzip <- function(Y, K, phi, prob, lambda) {

    u <- array(dim = c(dim(Y), K))

	for (k in 1:K) {
	    num <- phi[k]
		den <- phi[k] + (1 - phi[k]) * exp(-lambda[, , k])

        uk <- num / den
        uk <- matrix(uk, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
        uk[Y != 0] <- 0

        u[,, k] <- uk
	}

	return(u)
}


EMmgzip <- function(Y, totals, K, phi, prob, beta0, rho,
                    tol = 1e-8, maxiter = 1000, maxiter_NR = 100) {

	N <- nrow(Y)
	G <- ncol(Y)

	phi.h   <- phi
	pi.h    <- prob
	rho.h   <- rho
	beta0.h <- beta0

	lambda.h <- array(
		exp(rep(log(totals), times = G * K) +
				rep(rep(beta0, each = N), times = K) +
				rep(t(rho.h), each = N)),
		c(N, G, K))

	ll   <- llmgzip(Y, K, phi.h, pi.h, lambda.h)
	ll.0 <- ll$loglik

	dl   <- Inf
	iter <- 1

	while (dl > tol & iter < maxiter) {

		# E step
		z.hat <- ll$z
		u.hat <- u.mgzip(Y, K, phi.h, pi.h, lambda.h)

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

		# NR algorithm
		nr       <- NRmgzip(Y, totals, z.hat, u.hat, rho.h, beta0.h, lambda.h,
							tol = tol, maxiter = maxiter_NR)
		rho.h    <- nr$rho
		beta0.h  <- nr$beta0
		lambda.h <- nr$lambda

		# updated log-likelihood
		ll    <- llmgzip(Y, K, phi.h, pi.h, lambda.h, ll$trick, ll$thresh)
		ll.n  <- ll$loglik
		trick <- ll$trick

		dl   <- abs(ll.n - ll.0)
		ll.0 <- ll.n
		iter <- iter + 1

		if (is.infinite(dl) | is.nan(dl)) {
		    warning("dl is infinite or NaN")
		    dl <- 0
		}
	}

	out <- list(
		prob   = pi.h,
		phi    = phi.h,
		rho    = rho.h,
		beta0  = beta0.h,
		U      = u.hat,
		Z      = z.hat,
		loglik = ll.0,
		iter   = iter
	)

	return(out)
}
