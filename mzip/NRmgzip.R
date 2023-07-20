NRmgzip <- function(Y, Tn, Z, U, rho, beta0, lambda, tol = 1e-8, maxiter = 100) {
	N <- nrow(Y)
	G <- ncol(Y)
	K <- ncol(Z)

	aux <- array(Z[, rep(1:K, each = G)], dim = c(N, G, K)) * (1 - U)

	x  <- c(t(rho) + rep(beta0, times = K))
	d1 <- numeric(K * G)
	d2 <- numeric(K * G)

	dl   <- 1 + tol
	iter <- 0
	while (iter < maxiter & dl > tol) {

		# Gradient and Hessian
		d1 <- c(colSums(aux * (rep(Y, K) - lambda)))
		d2 <- c(colSums(aux * (-lambda)))

		# Update the parameters
		dx <- -(d1 / d2)
		dx[is.na(dx) | is.infinite(dx)] <- 0

		x      <- x + dx
		tau    <- matrix(x, nrow = K, byrow = TRUE)
		lambda <- array(
			exp(rep(log(Tn), times = G * K) +
					rep(t(tau), each = N)),
			c(N, G, K))

		dl <- sum(abs(dx))
		iter <- iter + 1
	}

	beta0 <- colMeans(tau)
	rho   <- tau - rep(beta0, each = K)


	return(list(
		beta0  = beta0,
		rho    = rho,
		lambda = lambda,
		iter   = iter
	))
}
