rmgzip <- function(N, G, K, phi, prob, beta0, rho, mu, sigma) {

	# cluster assignments
	Z <- sample(K, size = N, replace = TRUE, prob = prob)

	# indicators of always-zero
	U <- matrix(rbinom(n = N * G, size = 1, prob = phi[Z]), N, G)

	# known size
	Sn <- pmax(0, round(rnorm(N, mean = mu, sd = sigma)))

	# rate parameters
	loglambda <- log(rep(Sn, G)[U == 0]) +
		rho[Z, ][U == 0] +
		rep(beta0, each = N)[U == 0]

	# generate the data
	Y <- matrix(0, nrow = N, ncol = G)
	Y[U == 0] <- rpois(n = sum(U == 0), lambda = exp(loglambda))

	return(list(
		Y  = Y,
		Z  = Z,
		U  = U,
		Sn = Sn
	))
}
