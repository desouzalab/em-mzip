# V-Measure function
# FDRSeg package
# https://cran.r-project.org/web/packages/FDRSeg/
#
# Rosenberg, A., and Hirschberg, J. (2007).
# V-measure: a conditional entropy-based external cluster evaluation measures.
# Proc. Conf. Empirical Methods Natural Lang. Process., (June):410--420.

.entropy_stepF <- function(y, N) {
	n <- length(y)
	p <- diff(c(0, which(diff(y) != 0), n))
	sum(-p/N * log2(p/n))
}


v_measure  <- function(sig, est, beta = 1) {
	n <- length(sig)
	C <- c(0, which(diff(sig) != 0), n)
	K <- c(0, which(diff(est) != 0), n)
	HC <- .entropy_stepF(sig, n)
	if (HC == 0) {
		h <- 1
	}
	else {
		HCK <- 0
		for (i in 1:(length(K) - 1)) {
			HCK <- HCK + .entropy_stepF(sig[(K[i] + 1):K[i +
														 	1]], n)
		}
		h <- 1 - HCK/HC
	}
	HK <- .entropy_stepF(est, n)
	if (HK == 0) {
		c <- 1
	}
	else {
		HKC <- 0
		for (i in 1:(length(C) - 1)) {
			HKC <- HKC + .entropy_stepF(est[(C[i] + 1):C[i +
														 	1]], n)
		}
		c <- 1 - HKC/HK
	}
	(1 + beta) * h * c/(beta * h + c)
}
