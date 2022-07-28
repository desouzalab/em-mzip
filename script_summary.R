source("mzip_summary.R")
library(xtable)
library(FDRSeg)
dir.create("fig", showWarnings = FALSE)
dir.create("tab", showWarnings = FALSE)

options(xtable.include.rownames = FALSE,
		xtable.booktabs = TRUE,
		xtable.caption.placement = "top",
		xtable.sanitize.text.function = function(x) {x})


# scen1 <- mzip_summary("Scenario 1")
# vmat1 <- mzip_summary_Z("Scenario 1")
#
# scen2 <- mzip_summary("Scenario 2")
# vmat2 <- mzip_summary_Z("Scenario 2")
#
# scen3 <- mzip_summary("Scenario 3")
# vmat3 <- mzip_summary_Z("Scenario 3")
#
# scen4 <- mzip_summary("Scenario 4")
# vmat4 <- mzip_summary_Z("Scenario 4")
#
# scen5 <- mzip_summary("Scenario 5")
# vmat5 <- mzip_summary_Z("Scenario 5")

tab_captions <- list(
	phi    = c(paste0("\\scenario{", 1:5, "} Mean and standard error for the estimates of $\\phi_k$ for each $k$ and each ", c("$N$", "$G$", "$K$", "case", "case"), ", obtained using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", 1:5, "}."),
			   paste0("Scenario ", 1:5, ". Summary for $\\hat\\phi_k$")),
	pi    = c(paste0("\\scenario{", 1:5, "} Mean and standard error for the estimates of $\\pi_k$ for each $k$ and each ", c("$N$", "$G$", "$K$", "case", "case"), ", obtained using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", 1:5, "}."),
			  paste0("Scenario ", 1:5, ". Summary for $\\hat\\pi_k$")),
	lambda = c(paste0("\\scenario{", 1:5, "} Mean squared error for the estimates of $\\lambda_k$ for each $k$ and each ", c("$N$", "$G$", "$K$", "case", "case"), ", using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", 1:5, "}."),
			   paste0("Scenario ", 1:5, ". MSE for $\\hat\\lambda_k$")),
	iter   = c(paste0("\\scenario{", 1:5, "} Mean and standard error for the number of iterations until the EM algorithm converged, by ", c("$N$", "$G$", "$K$", "case", "case"), ", for the datasets simulated from the settings described in Table~\\ref{t:sim", 1:5, "}."),
			   paste0("Scenario ", 1:5, ". Summary for the number of iterations")),
	times  = c(paste0("\\scenario{", 1:5, "} Mean and standard error for the EM algorithm computing times, in seconds, by ", c("$N$", "$G$", "$K$", "case", "case"), ", for the datasets simulated from the settings described in Table~\\ref{t:sim", 1:5, "}."),
			   paste0("Scenario ", 1:5, ". Summary of computing times"))
)


# Scenario 1 --------------------------------------------------------------

scenario <- 1

estim <- setdiff(names(scen1), "lambda")
labels   <- c("N = 12", 60, 120, 600, 1200)
cases    <- c(12, 60, 120, 600, 1200)

{
	pdf("fig/sim1.pdf", width = 5.2, height = 2.5)
	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:3) {
		boxplot(scen1$phi[, , k],
				main = captions[k],
				xlab = "N",
				names = cases,
				ylim = range(scen1$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3])
	for (k in 1:3) {
		boxplot(scen1$prob[, , k],
				main = captions[k],
				xlab = "N",
				names = cases,
				ylim = range(scen1$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen1$iter, xlab = "N", names = cases)
	barplot(colSums(scen1$times), xlab = "N", names = cases)
	boxplot(vmat1, xlab = "N", names = cases)

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen1[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen1[[e]]),
			SE   = apply(scen1[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$N$")
	} else if (length(dim(scen1[[e]])) == 3) {
		aux <- NULL
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(2, 1, 2)),
				case = cases,
				Mean = colMeans(scen1[[e]][, , j]),
				SE   = apply(scen1[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1:2] <- c("$k$", "$N$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim1-times",
			 digits = 5),
	  file = "tab/sim1-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim1-iter",
			 digits = 5),
	  file = "tab/sim1-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim1-phi",
			 digits = 5),
	  file = "tab/sim1-phi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim1-pi",
			 digits = 5),
	  file = "tab/sim1-pi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(scen1$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim1-lamdba",
			 digits = 5),
	  file = "tab/sim1-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{5}{c}{$N$}\\\\\n",
	  				  			"\\cmidrule(l){2-6}\n",
	  				  			"$k$ & 12 & 60 & 120 & 600 & 1200\\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 2 --------------------------------------------------------------

scenario <- 2

estim <- setdiff(names(scen2), "lambda")
labels   <- c("G = 12", 60, 120, 600, 1500)
cases    <- c(12, 60, 120, 600, 1500)

{
	pdf("fig/sim2.pdf", width = 5.2, height = 2.5)
	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:3) {
		boxplot(scen2$phi[, , k],
				main = captions[k],
				names = cases,
				xlab = "G",
				ylim = range(scen2$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3])
	for (k in 1:3) {
		boxplot(scen2$prob[, , k],
				main = captions[k],
				names = cases,
				xlab = "G",
				ylim = range(scen2$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen2$iter, names = cases, xlab = "G")
	barplot(colSums(scen2$times), names = cases, xlab = "G")
	boxplot(vmat2, names = cases, xlab = "G")

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen2[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen2[[e]]),
			SE   = apply(scen2[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$G$")
	} else if (length(dim(scen2[[e]])) == 3) {
		aux <- NULL
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(2, 1, 2)),
				case = cases,
				Mean = colMeans(scen2[[e]][, , j]),
				SE   = apply(scen2[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1:2] <- c("$k$", "$G$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim2-times"),
	  file = "tab/sim2-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim2-iter"),
	  file = "tab/sim2-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim2-phi"),
	  file = "tab/sim2-phi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim2-pi"),
	  file = "tab/sim2-pi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(scen2$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim2-lamdba"),
	  file = "tab/sim2-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{5}{c}{$G$}\\\\\n",
	  				  			"\\cmidrule(l){2-6}\n",
	  				  			"$k$ & 12 & 60 & 120 & 600 & 1500\\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 3 --------------------------------------------------------------

scenario <- 3

estim  <- setdiff(names(scen3), "lambda")
labels <- c("K = 1", 2, 3, 5)
cases  <- c(1, 2, 3, 5)
xlab_k <- c(rep("K", 3), rep("K = 5", 2))

{
	pdf("fig/sim3.pdf", width = 5.2, height = 2.5)
	par(mfrow = c(1, 5), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3], phi[4], phi[5])

	for (k in 1:5) {
		j <- k <= cases
		boxplot(scen3$phi[, j, k],
				main = captions[k],
				names = cases[j],
				xlab = xlab_k[j],
				ylim = range(scen3$phi, na.rm = TRUE))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3], pi[4], pi[5])
	for (k in 1:5) {
		j <- k <= cases
		boxplot(scen3$prob[, j, k],
				main = captions[k],
				names = cases[j],
				xlab = xlab_k[j],
				ylim = range(scen3$prob, na.rm = TRUE))
		abline(h = 1 / cases[j], col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen3$iter, names = cases, xlab = "K")
	barplot(colSums(scen3$times), names = cases, xlab = "K")
	boxplot(vmat3, names = cases, xlab = "K")
	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen3[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen3[[e]]),
			SE   = apply(scen3[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$K$")
	} else if (length(dim(scen3[[e]])) == 3) {
		aux <- NULL
		for (j in 1:5) {
			aux <- rbind(aux, data.frame(
				k    = "",
				case = cases,
				Mean = colMeans(scen3[[e]][, , j]),
				SE   = apply(scen3[[e]][, , j], 2, sd)
			))
		}
		aux[c(2, 7, 11, 16, 20), 1] <- 1:5
		colnames(aux)[1:2] <- c("$k$", "$K$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim3-times"),
	  file = "tab/sim3-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim3-iter"),
	  file = "tab/sim3-iter.tex")

print(xtable(tables$phi[complete.cases(tables$phi), ], auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim3-phi"),
	  file = "tab/sim3-phi.tex",
	  hline.after = c(-1, 0, 4, 7, 9, 10, 11))

print(xtable(tables$prob[complete.cases(tables$prob), ], auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim3-pi"),
	  file = "tab/sim3-pi.tex",
	  hline.after = c(-1, 0, 4, 7, 9, 10, 11))

print(xtable(scen3$lambda, auto = TRUE,# align = "rrrrr",
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim3-lamdba"),
	  file = "tab/sim3-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{4}{c}{$K$}\\\\\n",
	  				  			"\\cmidrule(l){2-5}\n",
	  				  			"$k$ & 1 & 2 & 3 & 5 \\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 4 --------------------------------------------------------------

scenario <- 4

estim <- setdiff(names(scen4), "lambda")
labels   <- c("Case = 1", 2, 3, 4)
cases    <- 1:4

{
	pdf("fig/sim4.pdf", width = 5.2, height = 2.5)
	par(mfrow = c(1, 2), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2])
	for (k in 1:2) {
		boxplot(scen4$phi[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen4$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	probs <- matrix(c(.5, .25, .1, .05, .5, .75, .9, .95), ncol = 2)
	captions <- expression(pi[1], pi[2])
	for (k in 1:2) {
		boxplot(scen4$prob[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen4$prob))
		abline(h = probs[, k], col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen4$iter, names = cases, xlab = "Case")
	barplot(colSums(scen4$times), names = cases, xlab = "Case")
	boxplot(vmat4, names = cases, xlab = "Case")
	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen4[[e]])) == 2) {
		tables[[e]] <- data.frame(
			Case = cases,
			Mean = colMeans(scen4[[e]]),
			SE   = apply(scen4[[e]], 2, sd)
		)
	} else if (length(dim(scen4[[e]])) == 3) {
		aux <- NULL
		for (j in 1:2) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(1, 1, 2)),
				Case = cases,
				Mean = colMeans(scen4[[e]][, , j]),
				SE   = apply(scen4[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1] <- c("$k$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim4-times"),
	  file = "tab/sim4-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim4-iter"),
	  file = "tab/sim4-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim4-phi"),
	  file = "tab/sim4-phi.tex",
	  hline.after = c(-1, seq(0, 8, 4)))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim4-pi"),
	  file = "tab/sim4-pi.tex",
	  hline.after = c(-1, seq(0, 8, 4)))

print(xtable(scen4$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim4-lamdba"),
	  file = "tab/sim4-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{4}{c}{Case}\\\\\n",
	  				  			"\\cmidrule(l){2-5}\n",
	  				  			"$k$ & 1 & 2 & 3 & 4 \\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 5 --------------------------------------------------------------

scenario <- 5

estim <- setdiff(names(scen5), "lambda")
labels   <- c("Case = 1", 2:5)
cases    <- 1:5

{
	pdf("fig/sim5.pdf", width = 5.2, height = 2.5)
	par(mfrow = c(1, 2), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2])
	for (k in 1:2) {
		boxplot(scen5$phi[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen5$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	probs <- matrix(c(.5, .25, .1, .05, .5, .75, .9, .95), ncol = 2)
	captions <- expression(pi[1], pi[2])
	for (k in 1:2) {
		boxplot(scen5$prob[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen5$prob))
		abline(h = .5, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen5$iter, names = cases, xlab = "Case")
	barplot(colSums(scen5$times), names = cases, xlab = "Case")
	boxplot(vmat5, names = cases, xlab = "Case")
	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen5[[e]])) == 2) {
		tables[[e]] <- data.frame(
			Case = cases,
			Mean = colMeans(scen5[[e]]),
			SE   = apply(scen5[[e]], 2, sd)
		)
		# colnames(tables[[e]])[1] <- c("$N$")
	} else if (length(dim(scen5[[e]])) == 3) {
		aux <- NULL
		for (j in 1:2) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(2, 1, 2)),
				Case = cases,
				Mean = colMeans(scen5[[e]][, , j]),
				SE   = apply(scen5[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1] <- c("$k$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim5-times"),
	  file = "tab/sim5-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim5-iter"),
	  file = "tab/sim5-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim5-phi"),
	  file = "tab/sim5-phi.tex",
	  hline.after = c(-1, 0, 5, 10))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim5-pi"),
	  file = "tab/sim5-pi.tex",
	  hline.after = c(-1, 0, 5, 10))

print(xtable(scen5$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim5-lamdba"),
	  file = "tab/sim5-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{5}{c}{Case}\\\\\n",
	  				  			"\\cmidrule(l){2-6}\n",
	  				  			"$k$ & 1 & 2 & 3 & 4 & 5 \\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)
