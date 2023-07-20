library(xtable)
options(xtable.include.rownames = FALSE,
		xtable.booktabs = TRUE,
		xtable.caption.placement = "top",
		xtable.sanitize.text.function = function(x) {x})

load("summary_scen_11-16.RData")

dir.create("fig-1", showWarnings = FALSE)
dir.create("tab-1", showWarnings = FALSE)

estim <- names(scen11)

n_digits     <- 5
n_digits_max <- 9

# Captions ----------------------------------------------------------------

nums <- 1:6
vars <- c("$N$", "$G$", "$K$", "case", "case", "case")

tab_captions <- list(
    phi    = c(paste0("\\scenario{", nums, "} Mean and standard error for the estimates of $\\phi_k$ for each $k$ and each ", vars, ", obtained using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}."),
               paste0("Scenario ", nums, ". Summary for $\\hat\\phi_k$")),
    pi    = c(paste0("\\scenario{", nums, "} Mean and standard error for the estimates of $\\pi_k$ for each $k$ and each ", vars, ", obtained using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}."),
              paste0("Scenario ", nums, ". Summary for $\\hat\\pi_k$")),
    lambda = c(paste0("\\scenario{", nums, "} Mean squared error for the estimates of $\\lambda_k$ for each $k$ and each ", vars, ", using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}."),
               paste0("Scenario ", nums, ". Mean squared error for $\\hat\\lambda_k$")),
    iter   = c(paste0("\\scenario{", nums, "} Mean and standard error for the number of iterations until the EM algorithm converged, by ", vars, ", across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}."),
               paste0("Scenario ", nums, ". Summary for the number of iterations")),
    times  = c(paste0("\\scenario{", nums, "} Mean and standard error for the EM algorithm computing times, in seconds, by ", vars, ", across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}."),
               paste0("Scenario ", nums, ". Summary of computing times"))
)

fig_captions <- list(
    phi = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-1/sim1", nums, "-1.pdf}
\\caption[Scenario ", nums, ". Boxplots of \\(\\hat\\phi_k\\)]{\\scenario{", nums, "} Boxplots for the estimates of \\(\\phi_k\\) using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}. Red lines correspond to true values.
See also Table~\\ref{t:sim1", nums, "-phi}.}
\\label{f:sim1", nums, "-phi}
\\end{figure}"),
    pi = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-1/sim1", nums, "-2.pdf}
\\caption[Scenario ", nums, ". Boxplots of \\(\\hat\\pi_k\\)]{\\scenario{", nums, "} Boxplots for the estimates of \\(\\pi_k\\) using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}. Red lines correspond to true values.
See also Table~\\ref{t:sim1", nums, "-pi}.}
\\label{f:sim1", nums, "-pi}
\\end{figure}"),
    iter = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-1/sim1", nums, "-3.pdf}
\\caption[Scenario ", nums, ". Boxplots of the number of iterations]{\\scenario{", nums, "} Boxplots for the number of iterations until the EM algorithm converged across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}.
See also Table~\\ref{t:sim1", nums, "-iter}.}
\\label{f:sim1", nums, "-iter}
\\end{figure}"),
    times = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-1/sim1", nums, "-4.pdf}
\\caption[Scenario ", nums, ". Bar plot of the computing times]{\\scenario{", nums, "} Bar plot for the total computing times, in seconds, taken for the EM algorithm to converge across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}.
See also Table~\\ref{t:sim1", nums, "-times}.}
\\label{f:sim1", nums, "-times}
\\end{figure}"),
    vm = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-1/sim1", nums, "-5.pdf}
\\caption[Scenario ", nums, ". Boxplots of the \\(V\\)-measures]{\\scenario{", nums, "} Boxplots for the \\(V\\)-measures of the clustering obtained by the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim1", nums, "}.}
\\label{f:sim1", nums, "-vm}
\\end{figure}"))

for (x in names(fig_captions)) {
    for (i in nums) {
        cat(fig_captions[[x]][i],
            file = paste0("fig-1/sim1", i, "-", x, ".tex"))
    }
}


# Scenario 11 -------------------------------------------------------------

scenario <- 1

estim <- setdiff(names(scen11), "lambda")
labels   <- c("N = 12", 60, 120, 600, 1200)
cases    <- c(12, 60, 120, 600, 1200)

{
	pdf("fig-1/sim11-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:3) {
		boxplot(scen11$phi[, , k],
				main = captions[k],
				xlab = "N",
				names = cases,
				ylim = range(scen11$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3])
	for (k in 1:3) {
		boxplot(scen11$prob[, , k],
				main = captions[k],
				xlab = "N",
				names = cases,
				ylim = range(scen11$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen11$iter, xlab = "N", names = cases)
	barplot(colSums(scen11$times), xlab = "N", names = cases)
	boxplot(vmat11, xlab = "N", names = cases)

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen11[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen11[[e]]),
			SE   = apply(scen11[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$N$")
	} else if (length(dim(scen11[[e]])) == 3) {
		aux <- NULL
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(2, 1, 2)),
				case = cases,
				Mean = colMeans(scen11[[e]][, , j]),
				SE   = apply(scen11[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1:2] <- c("$k$", "$N$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim11-times",
			 digits = 5),
	  file = "tab-1/sim11-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim11-iter",
			 digits = 5),
	  file = "tab-1/sim11-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim11-phi",
			 digits = 5),
	  file = "tab-1/sim11-phi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim11-pi",
			 digits = 5),
	  file = "tab-1/sim11-pi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(scen11$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim11-lambda",
			 digits = 5),
	  file = "tab-1/sim11-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{5}{c}{$N$}\\\\\n",
	  				  			"\\cmidrule(l){2-6}\n",
	  				  			"$k$ & 12 & 60 & 120 & 600 & 1200\\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 12 -------------------------------------------------------------

scenario <- 2

estim <- setdiff(names(scen12), "lambda")
labels   <- c("G = 12", 60, 120, 600, 1500)
cases    <- c(12, 60, 120, 600, 1500)

{
    pdf("fig-1/sim12-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:3) {
		boxplot(scen12$phi[, , k],
				main = captions[k],
				names = cases,
				xlab = "G",
				ylim = range(scen12$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3])
	for (k in 1:3) {
		boxplot(scen12$prob[, , k],
				main = captions[k],
				names = cases,
				xlab = "G",
				ylim = range(scen12$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen12$iter, names = cases, xlab = "G")
	barplot(colSums(scen12$times), names = cases, xlab = "G")
	boxplot(vmat12, names = cases, xlab = "G")

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen12[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen12[[e]]),
			SE   = apply(scen12[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$G$")
	} else if (length(dim(scen12[[e]])) == 3) {
		aux <- NULL
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(2, 1, 2)),
				case = cases,
				Mean = colMeans(scen12[[e]][, , j]),
				SE   = apply(scen12[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1:2] <- c("$k$", "$G$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim12-times",
			 digits = 5),
	  file = "tab-1/sim12-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim12-iter",
			 digits = 5),
	  file = "tab-1/sim12-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim12-phi",
			 digits = 5),
	  file = "tab-1/sim12-phi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim12-pi",
			 digits = 5),
	  file = "tab-1/sim12-pi.tex",
	  hline.after = c(-1, seq(0, 15, 5)))

print(xtable(scen12$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim12-lambda",
			 digits = 5),
	  file = "tab-1/sim12-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{5}{c}{$G$}\\\\\n",
	  				  			"\\cmidrule(l){2-6}\n",
	  				  			"$k$ & 12 & 60 & 120 & 600 & 1500\\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 13 -------------------------------------------------------------

scenario <- 3

estim  <- setdiff(names(scen13), "lambda")
labels <- c("K = 1", 2, 3, 5)
cases  <- c(1, 2, 3, 5)
xlab_k <- c(rep("K", 3), rep("K = 5", 2))

{
    pdf("fig-1/sim13-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 5), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3], phi[4], phi[5])

	for (k in 1:5) {
		j <- k <= cases
		boxplot(scen13$phi[, j, k],
				main = captions[k],
				names = cases[j],
				xlab = xlab_k[j],
				ylim = range(scen13$phi, na.rm = TRUE))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)
	for (k in 2:4) {
	    boxplot(as.numeric(scen13$prob[, k, ]),
	        main = expression(pi[k]),
	        xlab = paste("K =", cases[k]))
	    abline(h = 1 / cases[k], col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen13$iter, names = cases, xlab = "K")
	barplot(colSums(scen13$times), names = cases, xlab = "K")
	boxplot(vmat13, names = cases, xlab = "K")
	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen13[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen13[[e]]),
			SE   = apply(scen13[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$K$")
	} else if (length(dim(scen13[[e]])) == 3) {
		aux <- NULL
		for (j in 1:5) {
			aux <- rbind(aux, data.frame(
				k    = "",
				case = cases,
				Mean = colMeans(scen13[[e]][, , j]),
				SE   = apply(scen13[[e]][, , j], 2, sd)
			))
		}
		aux[c(2, 7, 11, 16, 20), 1] <- 1:5
		colnames(aux)[1:2] <- c("$k$", "$K$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim13-times",
			 digits = 5),
	  file = "tab-1/sim13-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim13-iter",
			 digits = 5),
	  file = "tab-1/sim13-iter.tex")

print(xtable(tables$phi[complete.cases(tables$phi), ], auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim13-phi",
			 digits = 5),
	  file = "tab-1/sim13-phi.tex",
	  hline.after = c(-1, 0, 4, 7, 9, 10, 11))

print(xtable(tables$prob[complete.cases(tables$prob), ], auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim13-pi",
			 digits = 5),
	  file = "tab-1/sim13-pi.tex",
	  hline.after = c(-1, 0, 4, 7, 9, 10, 11))

print(xtable(scen13$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim13-lambda",
			 digits = 5),
	  file = "tab-1/sim13-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{4}{c}{$K$}\\\\\n",
	  				  			"\\cmidrule(l){2-5}\n",
	  				  			"$k$ & 1 & 2 & 3 & 5 \\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 14 -------------------------------------------------------------

scenario <- 4

estim <- setdiff(names(scen14), "lambda")
labels   <- c("Case = 1", 2, 3, 4)
cases    <- 1:4

{
    pdf("fig-1/sim14-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 2), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2])
	for (k in 1:2) {
		boxplot(scen14$phi[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen14$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	probs <- matrix(c(.5, .25, .1, .05, .5, .75, .9, .95), ncol = 2)
	captions <- expression(pi[1], pi[2])
	for (k in 1:2) {
		boxplot(scen14$prob[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen14$prob))
		abline(h = probs[, k], col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen14$iter, names = cases, xlab = "Case")
	barplot(colSums(scen14$times), names = cases, xlab = "Case")
	boxplot(vmat14, names = cases, xlab = "Case")
	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen14[[e]])) == 2) {
		tables[[e]] <- data.frame(
			Case = cases,
			Mean = colMeans(scen14[[e]]),
			SE   = apply(scen14[[e]], 2, sd)
		)
	} else if (length(dim(scen14[[e]])) == 3) {
		aux <- NULL
		for (j in 1:2) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(1, 1, 2)),
				Case = cases,
				Mean = colMeans(scen14[[e]][, , j]),
				SE   = apply(scen14[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1] <- c("$k$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim14-times",
			 digits = 5),
	  file = "tab-1/sim14-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim14-iter",
			 digits = 5),
	  file = "tab-1/sim14-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim14-phi",
			 digits = 5),
	  file = "tab-1/sim14-phi.tex",
	  hline.after = c(-1, seq(0, 8, 4)))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim14-pi",
			 digits = 5),
	  file = "tab-1/sim14-pi.tex",
	  hline.after = c(-1, seq(0, 8, 4)))

print(xtable(scen14$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim14-lambda",
			 digits = 5),
	  file = "tab-1/sim14-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{4}{c}{Case}\\\\\n",
	  				  			"\\cmidrule(l){2-5}\n",
	  				  			"$k$ & 1 & 2 & 3 & 4 \\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 15 -------------------------------------------------------------

scenario <- 5

estim <- setdiff(names(scen15), "lambda")
labels   <- c("Case = 1", 2:5)
cases    <- 1:5

{
    pdf("fig-1/sim15-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 2), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2])
	for (k in 1:2) {
		boxplot(scen15$phi[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen15$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	probs <- matrix(c(.5, .25, .1, .05, .5, .75, .9, .95), ncol = 2)
	captions <- expression(pi[1], pi[2])
	for (k in 1:2) {
		boxplot(scen15$prob[, , k],
				main = captions[k],
				names = cases,
				xlab = "Case",
				ylim = range(scen15$prob))
		abline(h = .5, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen15$iter, names = cases, xlab = "Case")
	barplot(colSums(scen15$times), names = cases, xlab = "Case")
	boxplot(vmat15, names = cases, xlab = "Case")
	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen15[[e]])) == 2) {
		tables[[e]] <- data.frame(
			Case = cases,
			Mean = colMeans(scen15[[e]]),
			SE   = apply(scen15[[e]], 2, sd)
		)
	} else if (length(dim(scen15[[e]])) == 3) {
		aux <- NULL
		for (j in 1:2) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", j, ""), c(2, 1, 2)),
				Case = cases,
				Mean = colMeans(scen15[[e]][, , j]),
				SE   = apply(scen15[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1] <- c("$k$")
		tables[[e]] <- aux
	}
}

print(xtable(tables$times, auto = TRUE,
			 caption = tab_captions$times[scenario + c(0, 5)],
			 label = "t:sim15-times",
			 digits = 5),
	  file = "tab-1/sim15-times.tex")

print(xtable(tables$iter, auto = TRUE,
			 caption = tab_captions$iter[scenario + c(0, 5)],
			 label = "t:sim15-iter",
			 digits = 5),
	  file = "tab-1/sim15-iter.tex")

print(xtable(tables$phi, auto = TRUE,
			 caption = tab_captions$phi[scenario + c(0, 5)],
			 label = "t:sim15-phi",
			 digits = 5),
	  file = "tab-1/sim15-phi.tex",
	  hline.after = c(-1, 0, 5, 10))

print(xtable(tables$prob, auto = TRUE,
			 caption = tab_captions$pi[scenario + c(0, 5)],
			 label = "t:sim15-pi",
			 digits = 5),
	  file = "tab-1/sim15-pi.tex",
	  hline.after = c(-1, 0, 5, 10))

print(xtable(scen15$lambda, auto = TRUE,
			 caption = tab_captions$lambda[scenario + c(0, 5)],
			 label = "t:sim15-lambda",
			 digits = 5),
	  file = "tab-1/sim15-lambda.tex",
	  add.to.row = list(pos = list(0, 0, 0),
	  				  command = c("& \\multicolumn{5}{c}{Case}\\\\\n",
	  				  			"\\cmidrule(l){2-6}\n",
	  				  			"$k$ & 1 & 2 & 3 & 4 & 5 \\\\\n")),
	  include.rownames = TRUE, include.colnames = FALSE)


# Scenario 16 -------------------------------------------------------------

scenario <- 6
cases    <- 1:6

# Figures
{
    pdf("fig-1/sim16-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)

    par(mfrow = c(1, 6), mar = c(4, 2, 2, 0.1), cex = .6)
    for (k in 1:5) {
        boxplot(as.numeric(scen16$phi[, k, ]),
            main = expression(phi[k]),
            xlab = paste("Case", k))
        abline(h = c(.05, .1, .25, .5, .9)[k], col = 2, lty = 2, lwd = 2)
    }
    boxplot(scen16$phi[, 6, ],
        main = expression(phi[k]),
        xlab = "Case 6")
    abline(h = c(.05, .25, .5), col = 2, lty = 2, lwd = 2)

    par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)
    captions <- expression(pi[1], pi[2], pi[3])
    for (k in 1:3) {
        boxplot(scen16$prob[, , k],
                main = captions[k],
                xlab = "Case",
                names = cases,
                ylim = range(scen16$prob))
        abline(h = 1/3, col = 2, lty = 2, lwd = 2)
    }

    par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
    boxplot(scen16$iter, xlab = "Case", names = cases)
    barplot(colSums(scen16$times), xlab = "Case", names = cases)
    boxplot(vmat16, xlab = "Case", names = cases)

    dev.off()
}

tables <- list()
for (e in estim) {
    if (length(dim(scen16[[e]])) == 2) {
        tables[[e]] <- data.frame(
            case = cases,
            Mean = colMeans(scen16[[e]]),
            SE   = apply(scen16[[e]], 2, sd)
        )
        colnames(tables[[e]])[1] <- c("Case")
    } else if (length(dim(scen16[[e]])) == 3) {
        aux <- NULL
        e_tex <- ifelse(e == "phi", "\\phi_", "\\pi_")
        for (j in 1:3) {
            aux <- rbind(aux, data.frame(
                k    = rep(c("", sprintf("$%s%d$", e_tex, j), ""), c(2, 1, 3)),
                case = cases,
                Mean = colMeans(scen16[[e]][, , j]),
                SE   = apply(scen16[[e]][, , j], 2, sd)
            ))
        }
        colnames(aux)[1:2] <- c(sprintf("$%sk$", e_tex), "Case")
        tables[[e]] <- aux
    }
}

# Tables
{
    print(xtable(tables$times, auto = TRUE,
                 caption = tab_captions$times[scenario + c(0, 5)],
                 label = "t:sim16-times",
                 digits = n_digits),
          file = "tab-1/sim16-times.tex")

    print(xtable(tables$iter, auto = TRUE,
                 caption = tab_captions$iter[scenario + c(0, 5)],
                 label = "t:sim16-iter",
                 digits = n_digits),
          file = "tab-1/sim16-iter.tex")

    print(xtable(tables$phi, auto = TRUE,
                 caption = tab_captions$phi[scenario + c(0, 5)],
                 label = "t:sim16-phi",
                 digits = n_digits),
          file = "tab-1/sim16-phi.tex",
          hline.after = c(-1, seq(0, 18, 6)))

    print(xtable(tables$prob, auto = TRUE,
                 caption = tab_captions$pi[scenario + c(0, 5)],
                 label = "t:sim16-pi",
                 digits = n_digits),
          file = "tab-1/sim16-pi.tex",
          hline.after = c(-1, seq(0, 18, 6)))

    print(xtable(scen16$lambda, auto = TRUE,
                 caption = tab_captions$lambda[scenario + c(0, 5)],
                 label = "t:sim16-lambda",
                 digits = 5),
          file = "tab-1/sim16-lambda.tex",
          add.to.row = list(pos = list(0, 0, 0),
                            command = c("& \\multicolumn{6}{c}{Case}\\\\\n",
                                        "\\cmidrule(l){2-7}\n",
                                        "$k$ & 1 & 2 & 3 & 4 & 5 & 6 \\\\\n")),
          include.rownames = TRUE, include.colnames = FALSE)
}


# Moving files ----

dir.create("Results", showWarnings = FALSE)
file.copy("tab-1", "Results", recursive = TRUE, copy.date = TRUE)
file.copy("fig-1", "Results", recursive = TRUE, copy.date = TRUE)
unlink("tab-1", recursive = TRUE)
unlink("fig-1", recursive = TRUE)
