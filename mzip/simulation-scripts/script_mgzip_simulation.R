source("../dzip.R")
source("../llmgzip.R")
source("../rmgzip.R")
source("../EMmgzip.R")

source("mgzip_simulation.R")

S <- 256

aux_dir <- "~/scratch/Data/"
dir.create("~/scratch", showWarnings = FALSE)
dir.create(aux_dir, showWarnings = FALSE)

sim.mean <- 1000
sim.sd   <- 100


# Scenario 21 - Vary N ----

N.list <- c(12, 60, 120, 600, 1200)

G     <- 120
K     <- 3
phi   <- rep(0.1, K)
prob  <- rep(1/K, K)
beta0 <- rep(-4.7, G)
rho_  <- c(-0.6, 0, 0.6)
rho   <- matrix(rep(rho_[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / K),
				 K, G, byrow = TRUE)
seed  <- 2100

for (i in seq_along(N.list)) {
	N <- N.list[i]
	mgzip_simulation(N, G, K, phi, prob, beta0, rho, sim.mean, sim.sd, S,
					seed = seed + i, main_dir = "Scenario 21")
}


# Scenario 22 - Vary G ----

G.list <- c(12, 60, 120, 600, 1200, 6000)

N    <- 1200
K    <- 3
phi  <- rep(0.1, K)
prob <- rep(1/K, K)
rho_ <- c(-0.6, 0, 0.6)
seed <- 2200

for (i in seq_along(G.list)) {
	G     <- G.list[i]
	beta0 <- rep(-4.7, G)
	rho   <- matrix(rep(rho_[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / 3),
					  K, G, byrow = TRUE)

	mgzip_simulation(N, G, K, phi, prob, beta0, rho, sim.mean, sim.sd, S,
					seed = seed + i, main_dir = "Scenario 22")

}


# Scenario 23 - Vary K ----

K.list   <- c(1, 2, 3, 5)

N        <- 1200
G        <- 120
seed     <- 2300

beta0    <- rep(-4.7, G)
rho_     <- c(-0.8, -0.6, 0, 0.6, 0.8)
rho.list <- list(
	3,
	c(2, 4, 4, 2),
	c(2, 3, 4, 3, 4, 2, 4, 2, 3),
	c(1, 2, 3, 4, 5, 2, 3, 4, 5, 1, 3, 4, 5, 1, 2, 4, 5, 1, 2, 3, 5, 1, 2, 3, 4)
)

for (i in seq_along(K.list)) {
	K    <- K.list[i]
	phi  <- rep(0.1, K)
	prob <- rep(1/K, K)
	rho  <- matrix(rep(rho_[rho.list[[i]]], each = G / K),
					 K, G, byrow = TRUE)

	mgzip_simulation(N, G, K, phi, prob, beta0, rho, sim.mean, sim.sd, S,
					seed = seed + i, main_dir = "Scenario 23")
}


# Scenario 25 - Vary initial values ----

N     <- 1200
G     <- 120
K     <- 3
phi   <- rep(0.1, K)
prob  <- rep(1/K, K)
beta0 <- rep(-4.7, G)
rho_  <- c(-0.6, 0, 0.6)
rho   <- matrix(rep(rho_[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / 3),
	K, G, byrow = TRUE)
seed  <- 2500

mgzip_simulation(N, G, K, phi, prob, beta0, rho, sim.mean, sim.sd, S,
				 seed = seed, main_dir = "Scenario 25")

{
	n <- list.files(paste0(aux_dir, "Scenario 25"))
	n <- regmatches(n, regexec("Case ([0-9]+)", n))
	n <- as.numeric(sapply(n, `[`, 2))
	n <- max(n, na.rm = TRUE)

	file_old1 <- paste0(aux_dir, "Scenario 25/Case ", n)
	file_old2 <- paste0(aux_dir, "Case ", n)
	file_new  <- paste0(aux_dir, "Scenario 25/Case ", n + 1)

	cat("Creating directory", file_new, "\n")
	aux1 <- file.copy(file_old1, aux_dir, recursive = TRUE, copy.date = TRUE)
	aux2 <- file.rename(file_old2, file_new)
	if (!aux1 | !aux2) {
		cat("Directory not copied\n")
	} else {
		cat("Simulations copied\n")
	}
}


# Scenario 26 - Vary pi, cluster balance ----

N       <- 1200
G       <- 120
K       <- 2
phi     <- rep(0.1, K)
pi.list <- c(.50, .25, .10, .05)
beta0   <- rep(-4.7, G)
rho_    <- c(-0.6, 0.6)
rho     <- matrix(rep(rho_[c(1, 2, 2, 1)], each = G / K),
	K, G, byrow = TRUE)
seed    <- 2600

for (i in seq_along(pi.list)) {
	prob <- c(pi.list[i], 1 - pi.list[i])
	mgzip_simulation(N, G, K, phi, prob, beta0, rho, sim.mean, sim.sd, S,
		seed = seed + i, main_dir = "Scenario 26")
}


# Scenario 28 - Vary phi, proportion of always-zero ----

N        <- 1200
G        <- 120
K        <- 3
phi.list <- c(.05, .1, .25, .5, .9)
prob     <- rep(1/K, K)
beta0    <- rep(-4.7, G)
rho_     <- c(-0.6, 0, 0.6)
rho      <- matrix(rep(rho_[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / K),
	K, G, byrow = TRUE)
seed     <- 2800

for (i in seq_along(phi.list)) {
	phi <- rep(phi.list[i], K)
	mgzip_simulation(N, G, K, phi, prob, beta0, rho, sim.mean, sim.sd, S,
		seed = seed + i, main_dir = "Scenario 28")
}

i   <- i + 1
phi <- c(.05, .25, .5)
mgzip_simulation(N, G, K, phi, prob, beta0, rho, sim.mean, sim.sd, S,
	seed = seed + i, main_dir = "Scenario 28")
