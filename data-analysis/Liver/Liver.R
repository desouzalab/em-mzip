# Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S,
# Saadatpour A, Zhou Z, Chen H, Ye F, Huang D.
# Mapping the mouse cell atlas by microwell-seq.
# Cell. 2018 Feb 22;172(5):1091-107.

dir.create("Data", showWarnings = FALSE)

d1 <- t(read.csv("../../../Raw Data/Han/Liver1.txt", sep = " "))
d2 <- t(read.csv("../../../Raw Data/Han/Liver2.txt", sep = " "))

genes <- intersect(colnames(d1), colnames(d2))
d1 <- d1[, colnames(d1) %in% genes]
d2 <- d2[, colnames(d2) %in% genes]

d1 <- d1[, order(colnames(d1))]
d2 <- d2[, order(colnames(d2))]
dim(d1) # 6137 15491
dim(d2) #  289 15491

d3 <- rbind(d1, d2)
dim(d3) # 6426 15491


z1 <- read.csv("../../../Raw Data/Han/MCA_CellAssignments.csv", header = TRUE)

cell <- intersect(rownames(d3), z1$Cell.name)

d3 <- d3[rownames(d3) %in% cell, ]
z1 <- z1[z1$Cell.name %in% cell, ]

d3 <- d3[order(rownames(d3)), ]
z1 <- z1[order(z1$Cell.name), ]

dim(d3) # 4685 cells, 15491 genes


# Sample cells

set.seed(7920728)
n <- sample(nrow(d3), 1000)
y <- d3[n, ]
z <- z1[n, ]

totals <- rowSums(y)


# Version 21

gsd <- apply(y, 2, sd)
sdstar <- sort(gsd, decreasing = TRUE)[100]
y1 <- y[, gsd >= sdstar]
ncol(y1) # 100 genes

Y <- y1
Ztrue <- z
Tn <- totals
names(Tn) <- NULL

save(Y, Ztrue, Tn, file = "Data/Liver v21.RData")
