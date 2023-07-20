# Klein AM, Mazutis L, Akartuna I, Tallapragada N et al.
# Droplet barcoding for single-cell transcriptomics
# applied to embryonic stem cells.
# Cell 2015 May 21;161(5):1187-1201. PMID: 26000487

dir.create("Data", showWarnings = FALSE)

y0 <- t(read.csv("../../../Raw Data/Klein/d0.csv", header = FALSE)[, -1])
y4 <- t(read.csv("../../../Raw Data/Klein/d4.csv", header = FALSE)[, -1])

y <- rbind(y0, y4)
z <- c(rep(0, nrow(y0)), rep(4, nrow(y4)))
totals <- rowSums(y)

dim(y)  # 1616 cells, 24175 genes

propzero <- apply(y, 2, function(x) {mean(x == 0)})
hist(propzero, xlim = c(0, 1), freq = FALSE)

giqr <- apply(y, 2, IQR)
table(giqr)
y <- y[, giqr > 1]
ncol(y)  # 4514 genes


# Version 6

gsd <- apply(y, 2, sd)
sdstar <- sort(gsd, decreasing = TRUE)[100]
y <- y[, gsd >= sdstar]
ncol(y) # 100

Y <- y
Ztrue <- z
Tn <- totals
names(Tn) <- NULL

save(Y, Ztrue, Tn, file = "Data/Klein v6.RData")
