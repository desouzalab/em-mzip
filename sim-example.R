library(pheatmap)
source("rmzip.R")

N       <- 120
G       <- 120
K       <- 3
phi     <- rep(0.1, K)
prob    <- rep(1/K, K)
lambda0 <- c(5, 10, 15)
lambda  <- matrix(rep(lambda0[c(1, 2, 3, 2, 3, 1, 3, 1, 2)], each = G / 3),
				  K, G, byrow = TRUE)

set.seed(19201)
x <- rmzip(N, G, K, phi, prob, lambda)

mat <- x$Y
rownames(mat) = paste("Cell", 1:N)


annotation_row = data.frame(
	k = as.factor(x$Z)
)
rownames(annotation_row) = paste("Cell", 1:N)

mycol <- RColorBrewer::brewer.pal(n = 9, name = "Purples")
mycol <- c(mycol[1:6], rep(mycol[7:9], each = 2))
mycol <- colorRampPalette(mycol)(30)

pheatmap(
	mat[order(x$Z), ],
	color = mycol,
	cluster_rows = FALSE,
	cluster_cols = FALSE,
	annotation_row = annotation_row,
	annotation_legend = FALSE,
	show_rownames = FALSE,
	gaps_row = c(45, 80),
	gaps_col = c(40, 80),
	# filename = "sim1-ex.pdf",
	# width = 6,
	# height = 6,
	asp = 1
)
