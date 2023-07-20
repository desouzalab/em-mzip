# Examples of simulating and estimating a data matrix
# using a mixture of zero-inflated Poisson model without covariates

setwd("examples")
library(pheatmap)


# Simulate a single data matrix ----

# Function for simulating the data
source("../mzip/rmzip.R")

# Parameters for the simulations
N       <- 120
G       <- 120
K       <- 3
phi     <- rep(0.1, K)
prob    <- rep(1/K, K)
lambda0 <- c(5, 10, 15)
lambda  <- matrix(rep(lambda0[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / 3),
				  K, G, byrow = TRUE)

# Calling the simulation function
set.seed(1)
simulation <- rmzip(N, G, K, phi, prob, lambda)
Y <- simulation$Y

{
	annotation_row <- data.frame(k = as.factor(simulation$Z))
	rownames(Y) <- paste("Cell", 1:N)
	rownames(annotation_row) <- rownames(Y)
	pheatmap(
		Y[order(simulation$Z), ],
		cellwidth         = 3.1,
		cellheight        = 3.1,
		cluster_rows      = FALSE,
		cluster_cols      = FALSE,
		annotation_row    = annotation_row,
		annotation_legend = FALSE,
		show_rownames     = FALSE,
		gaps_row          = c(45, 80),
		gaps_col          = c(40, 80),
		filename          = "example.pdf",
		width             = 6,
		height            = 5.5,
		asp               = 1
	)
}

# Saving data into a csv file, with row and column names
dimnames(Y) <- list(
	paste0("cell_gr", simulation$Z, "_", 1:N),
    paste0("gene_", 1:G)
)
write.csv(Y, "example-data.csv")


# Reading and estimating for a single data file ----

# Source files for the EM function
source("../mzip/dzip.R")
source("../mzip/llmzip.R")
source("../mzip/EMmzip.R")

# Reading the sample data
data <- read.csv("example-data.csv", row.names = 1)
data_mat <- as.matrix(data)

# Required initial values for the EM algorithm
K       <- 3
phi     <- rep(0.1, K)
prob    <- rep(1/K, K)
lambda0 <- c(5, 10, 15)
lambda  <- matrix(
	rep(lambda0[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = 40),
    nrow = K,
    ncol = 120,
    byrow = TRUE
)

# Calling the EM function
mod <- EMmzip(data_mat, K, phi, prob, lambda)


# Results

str(mod)
# List of 7
#  $ prob  : num [1:3] 0.333 0.275 0.392
#  $ phi   : num [1:3] 0.103 0.1 0.106
#  $ lambda: num [1:3, 1:120] 5.25 9.43 14.81 4.99 9.84 ...
#  $ U     : num [1:120, 1:120, 1:3] 0.957 0 0 0 0 ...
#  $ Z     : num [1:120, 1:3] 5.01e-122 8.26e-146 9.70e-137 1.00 7.82e-111 ...
#  $ loglik: num -37118
#  $ iter  : num 5
