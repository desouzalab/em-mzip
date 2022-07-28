source("dzip.R")
source("llmzip.R")
source("rmzip.R")
source("EMmzip.R")
source("mzip_simulation.R")

S <- 128

# Scenario 1 - Vary N ----

N.list  <- c(12, 60, 120, 600, 1200)
G       <- 120
K       <- 3
phi     <- rep(0.1, K)
prob    <- rep(1/K, K)
lambda0 <- c(5, 10, 15)
lambda  <- matrix(rep(lambda0[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / 3),
				 K, G, byrow = TRUE)
seed    <- 100

for (i in seq_along(N.list)) {
	N <- N.list[i]
	mzip_simulation(N, G, K, phi, prob, lambda, S,
					seed = seed + i, main_dir = "Scenario 1")
}


# Scenario 2 - Vary G ----

N       <- 1200
G.list  <- c(12, 60, 120, 600, 1500)
K       <- 3
phi     <- rep(0.1, K)
prob    <- rep(1/K, K)
lambda0 <- c(5, 10, 15)
seed    <- 200

for (i in seq_along(G.list)) {
	G      <- G.list[i]
	lambda <- matrix(rep(lambda0[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / 3),
					 K, G, byrow = TRUE)

	mzip_simulation(N, G, K, phi, prob, lambda, S,
					seed = seed + i, main_dir = "Scenario 2")

}


# Scenario 3 - Vary K ----

N       <- 1200
G       <- 120
K.list  <- c(1, 2, 3, 5)
lambda0 <- c(5, 10, 15, 20, 25)
seed    <- 300

lambda.list <- list(
	2,
	c(2, 3, 3, 2),
	c(1, 2, 3, 2, 3, 1, 3, 1, 2),
	c(1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4)
)

for (i in seq_along(K.list)) {
	K      <- K.list[i]
	phi    <- rep(0.1, K)
	prob   <- rep(1/K, K)
	lambda <- matrix(rep(lambda0[lambda.list[[i]]], each = G / K),
					 K, G, byrow = TRUE)

	mzip_simulation(N, G, K, phi, prob, lambda, S,
					seed = seed + i, main_dir = "Scenario 3")
}


# Scenario 4 - Vary cluster balance ----

N <- 1200
G <- 120
K <- 2
phi <- rep(0.1, K)
pi.list <- c(.50, .25, .10, .05)
lambda0 <- c(10, 15)
lambda <- matrix(rep(lambda0[c(1, 2, 2, 1)], each = G / K),
				 K, G, byrow = TRUE)
seed <- 400

for (i in seq_along(pi.list)) {
	prob <- c(pi.list[i], 1 - pi.list[i])
	mzip_simulation(N, G, K, phi, prob, lambda, S,
					seed = seed + i, main_dir = "Scenario 4")
}


# Scenario 5 - Vary cluster separation ----

N <- 1200
G <- 120
K <- 2
phi <- rep(0.1, K)
prob <- rep(1/K, K)
lambda0 <- c(4, 5)
lambda.sep <- c(0, 1/6, 1/3, 1/2, 2/3)
seed <- 500

for (i in seq_along(lambda.sep)) {
	ni <- round(lambda.sep[i] * G)
	nlambda <- c(c(1/3, 2/3, 1/3) * G, 2/3 * G - ni, ni)
	lambda <- matrix(rep(lambda0[c(1, 2, 2, 1, 2)], nlambda),
					 K, G, byrow = TRUE)
	mzip_simulation(N, G, K, phi, prob, lambda, S,
					seed = seed + i, main_dir = "Scenario 5")
}
