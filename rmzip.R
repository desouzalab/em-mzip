
# Simulate the Data ----
rmzip <- function(N, G, K, phi, prob, lambda) {

	# 1. Generate Z vector, cluster assignments
	Z <- sample(K, N, replace = TRUE, prob = prob)

	# 2. Generate U matrix, indicators of always-zero
	U <- matrix(rbinom(n = N * G, size = 1, prob = phi[Z]), N, G)

	# 3. Generate Y
	# Y <- matrix(rpois(n = N * G, lambda = lambda[Z, ]), N, G)
	# Y[U == 1] <- 0
	Y <- matrix(0, N, G)
	Y[U == 0] <- rpois(n = sum(U == 0), lambda = lambda[Z, ][U == 0])

	return(list(
	    Y = Y,
	    Z = Z,
	    U = U
	))
}
