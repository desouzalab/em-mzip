library(xtable)
options(xtable.include.rownames = FALSE,
		xtable.booktabs = TRUE,
		xtable.caption.placement = "top",
		xtable.sanitize.text.function = function(x) {x})

load("summary_scen_21-28.RData")

dir.create("fig-2", showWarnings = FALSE)
dir.create("tab-2", showWarnings = FALSE)

estim <- setdiff(names(scen21), c("rho", "beta0"))

n_digits     <- 5
n_digits_max <- 9

# Captions ----------------------------------------------------------------

scenarios <- c(1:3, 0, 4:5, 0, 6)
nums <- 21:28
vars <- c("$N$", "$G$", "$K$", "$G$", "case", "case", "case", "case")

tab_captions <- list(
    phi    = c(paste0("\\scenario{", scenarios, "} Mean and standard error for the estimates of $\\phi_k$ for each $k$ and each ", vars, ", obtained using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}."),
               paste0("Scenario ", scenarios, ". Summary for $\\hat\\phi_k$")),
    pi    = c(paste0("\\scenario{", scenarios, "} Mean and standard error for the estimates of $\\pi_k$ for each $k$ and each ", vars, ", obtained using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}."),
              paste0("Scenario ", scenarios, ". Summary for $\\hat\\pi_k$")),
    iter   = c(paste0("\\scenario{", scenarios, "} Mean and standard error for the number of iterations until the EM algorithm converged, by ", vars, ", across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}."),
               paste0("Scenario ", scenarios, ". Summary for the number of iterations")),
    times  = c(paste0("\\scenario{", scenarios, "} Mean and standard error for the EM algorithm computing times, in seconds, by ", vars, ", across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}."),
               paste0("Scenario ", scenarios, ". Summary of computing times")),
    lambda = c(paste0("\\scenario{", scenarios, "} Mean squared error for the estimates of $\\lambda_k$ for each $k$ and each ", vars, ", using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}."),
               paste0("Scenario ", scenarios, ". Mean squared error for $\\hat\\lambda_k$")),
    rho = c(paste0("\\scenario{", scenarios, "} Mean squared error for the estimates of $\\rho_{gk}$ for each $k$ and each ", vars, ", using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}."),
            paste0("Scenario ", scenarios, ". Mean squared error for $\\hat\\rho_{gk}$")),
    beta0 = c(paste0("\\scenario{", scenarios, "} Mean squared error for the estimates of $\\beta_{0g}$ for each ", vars, ", using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}."),
              paste0("Scenario ", scenarios, ". Mean squared error for $\\hat\\beta_{0g}$"))
)

fig_captions <- list(
    phi = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-2/sim", nums, "-1.pdf}
\\caption[Scenario ", scenarios, ". Boxplots of \\(\\hat\\phi_k\\)]{\\scenario{", scenarios, "} Boxplots for the estimates of \\(\\phi_k\\) using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}. Red lines correspond to true values.
See also Table~\\ref{t:sim", nums, "-phi}.}
\\label{f:sim", nums, "-phi}
\\end{figure}"),
    pi = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-2/sim", nums, "-2.pdf}
\\caption[Scenario ", scenarios, ". Boxplots of \\(\\hat\\pi_k\\)]{\\scenario{", scenarios, "} Boxplots for the estimates of \\(\\pi_k\\) using the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}. Red lines correspond to true values.
See also Table~\\ref{t:sim", nums, "-pi}.}
\\label{f:sim", nums, "-pi}
\\end{figure}"),
    iter = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-2/sim", nums, "-3.pdf}
\\caption[Scenario ", scenarios, ". Boxplots of the number of iterations]{\\scenario{", scenarios, "} Boxplots for the number of iterations until the EM algorithm converged across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}.
See also Table~\\ref{t:sim", nums, "-iter}.}
\\label{f:sim", nums, "-iter}
\\end{figure}"),
    times = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-2/sim", nums, "-4.pdf}
\\caption[Scenario ", scenarios, ". Bar plot of the computing times]{\\scenario{", scenarios, "} Bar plot for the total computing times, in seconds, taken for the EM algorithm to converge across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}.
See also Table~\\ref{t:sim", nums, "-times}.}
\\label{f:sim", nums, "-times}
\\end{figure}"),
    vm = paste0("\\begin{figure}
\\centering
\\includegraphics{fig-2/sim", nums, "-5.pdf}
\\caption[Scenario ", scenarios, ". Boxplots of the \\(V\\)-measures]{\\scenario{", scenarios, "} Boxplots for the \\(V\\)-measures of the clustering obtained by the EM algorithm across the datasets simulated from the settings described in Table~\\ref{t:sim", nums, "}.}
\\label{f:sim", nums, "-vm}
\\end{figure}"))

for (x in names(fig_captions)) {
    for (i in seq_along(nums)) {
        cat(fig_captions[[x]][i],
              file = paste0("fig-2/sim", nums[i], "-", x, ".tex"))
    }
}

tab_captions$beta0[c(1, 5, 9, 13)] <- sub(
    "Mean squared error", "Median absolute deviation",
    tab_captions$beta0[c(1, 5, 9, 13)])

tab_captions$rho[c(1, 5, 9, 13)] <- sub(
    "Mean squared error", "Median absolute deviation",
    tab_captions$rho[c(1, 5, 9, 13)])

aux <- function(mat) {
    estim <- c("rho_mad", "beta0_mad")

    for (e in estim) {
        j <- which.max(sapply(mat[[e]], length))
        nj <- length(mat[[e]][[j]])
        aux <- NULL

        for (case in names(mat[[e]])) {
            aux <- c(aux,
                     mat[[e]][[case]],
                     rep(NA, nj - length(mat[[e]][[case]])))
        }

        mat[[e]] <- array(aux, c(dim(mat[[e]][[j]]), length(mat[[e]])))

        if (length(dim(mat[[e]])) == 3) {
            if (dim(mat[[e]])[1] == 1) {
                mat[[e]] <- aperm(mat[[e]], c(2, 3, 1))[,, 1]
            }  else {
                mat[[e]] <- aperm(mat[[e]], c(1, 3, 2))
            }
        }
    }

    return(mat)
}

scen21 <- aux(scen21)
scen25 <- aux(scen25)


# Scenario 21 -------------------------------------------------------------

scenario <- 1
cases    <- c(12, 60, 120, 600, 1200)

scen21$rho   <- scen21$rho_mad
scen21$beta0 <- scen21$beta0_mad

# Figures
{
	pdf("fig-2/sim21-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:3) {
		boxplot(scen21$phi[, , k],
				main = captions[k],
				xlab = "N",
				names = cases,
				ylim = range(scen21$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3])
	for (k in 1:3) {
		boxplot(scen21$prob[, , k],
				main = captions[k],
				xlab = "N",
				names = cases,
				ylim = range(scen21$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen21$iter, xlab = "N", names = cases)
	barplot(colSums(scen21$times), xlab = "N", names = cases)
	boxplot(vmat21, xlab = "N", names = cases)

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen21[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen21[[e]]),
			SE   = apply(scen21[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$N$")
	} else if (length(dim(scen21[[e]])) == 3) {
		aux <- NULL
		e_tex <- ifelse(e == "phi", "\\phi_", "\\pi_")
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", sprintf("$%s%d$", e_tex, j), ""), c(2, 1, 2)),
				case = cases,
				Mean = colMeans(scen21[[e]][, , j]),
				SE   = apply(scen21[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1:2] <- c(sprintf("$%sk$", e_tex), "$N$")
		tables[[e]] <- aux
	}
}

# Tables
{
	print(xtable(tables$times, auto = TRUE,
				 caption = tab_captions$times[scenario + c(0, 5)],
				 label = "t:sim21-times",
				 digits = n_digits),
		  file = "tab-2/sim21-times.tex")

	print(xtable(tables$iter, auto = TRUE,
				 caption = tab_captions$iter[scenario + c(0, 5)],
				 label = "t:sim21-iter",
				 digits = n_digits),
		  file = "tab-2/sim21-iter.tex")

	print(xtable(tables$phi, auto = TRUE,
				 caption = tab_captions$phi[scenario + c(0, 5)],
				 label = "t:sim21-phi",
				 digits = n_digits),
		  file = "tab-2/sim21-phi.tex",
		  hline.after = c(-1, seq(0, 15, 5)))

	print(xtable(tables$prob, auto = TRUE,
				 caption = tab_captions$pi[scenario + c(0, 5)],
				 label = "t:sim21-pi",
				 digits = n_digits),
		  file = "tab-2/sim21-pi.tex",
		  hline.after = c(-1, seq(0, 15, 5)))

	print(xtable(scen21$rho, auto = TRUE,
				 caption = tab_captions$rho[scenario + c(0, 5)],
				 label = "t:sim21-rho",
				 digits = n_digits_max),
		  file = "tab-2/sim21-rho.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("& \\multicolumn{5}{c}{$N$}\\\\\n",
		  				  			"\\cmidrule(l){2-6}\n",
		  				  			"$k$ & 12 & 60 & 120 & 600 & 1200\\\\\n")),
		  include.rownames = TRUE, include.colnames = FALSE)

	print(xtable(scen21$beta0, auto = TRUE,
				 caption = tab_captions$beta0[scenario + c(0, 5)],
				 label = "t:sim21-beta0",
				 digits = n_digits),
		  file = "tab-2/sim21-beta0.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("\\multicolumn{5}{c}{$N$}\\\\\n",
		  				  			"\\midrule\n",
		  				  			"12 & 60 & 120 & 600 & 1200\\\\\n")),
		  include.rownames = FALSE, include.colnames = FALSE)
}


# Scenario 22 -------------------------------------------------------------

scenario <- 2
cases    <- c(12, 60, 120, 600, 1200, 6000)

# Figures
{
	pdf("fig-2/sim22-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .4)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:3) {
		boxplot(scen22$phi[, , k],
				main = captions[k],
				names = cases,
				xlab = "G",
				ylim = range(scen22$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3])
	for (k in 1:3) {
		boxplot(scen22$prob[, , k],
				main = captions[k],
				names = cases,
				xlab = "G",
				ylim = range(scen22$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .4)
	boxplot(scen22$iter, names = cases, xlab = "G")
	barplot(colSums(scen22$times), names = cases, xlab = "G")
	boxplot(vmat22, names = cases, xlab = "G")

	boxplot(b0_22_1, xlab = expression(hat(beta)[0][g]))
	abline(h = -4.7, col = 2, lty = 2, lwd = 2)

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen22[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen22[[e]]),
			SE   = apply(scen22[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$G$")
	} else if (length(dim(scen22[[e]])) == 3) {
		aux <- NULL
		e_tex <- ifelse(e == "phi", "\\phi_", "\\pi_")
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", sprintf("$%s%d$", e_tex, j), ""), c(2, 1, 3)),
				case = cases,
				Mean = colMeans(scen22[[e]][, , j]),
				SE   = apply(scen22[[e]][, , j], 2, sd)
			))
		}
		colnames(aux)[1:2] <- c(sprintf("$%sk$", e_tex), "$G$")
		tables[[e]] <- aux
	}
}

# Tables
{
	print(xtable(tables$times, auto = TRUE,
				 caption = tab_captions$times[scenario + c(0, 5)],
				 label = "t:sim22-times",
				 digits = n_digits),
		  file = "tab-2/sim22-times.tex")

	print(xtable(tables$iter, auto = TRUE,
				 caption = tab_captions$iter[scenario + c(0, 5)],
				 label = "t:sim22-iter",
				 digits = n_digits),
		  file = "tab-2/sim22-iter.tex")

	print(xtable(tables$phi, auto = TRUE,
				 caption = tab_captions$phi[scenario + c(0, 5)],
				 label = "t:sim22-phi",
				 digits = n_digits),
		  file = "tab-2/sim22-phi.tex",
		  hline.after = c(-1, seq(0, 20, 6)))

	print(xtable(tables$prob, auto = TRUE,
				 caption = tab_captions$pi[scenario + c(0, 5)],
				 label = "t:sim22-pi",
				 digits = n_digits),
		  file = "tab-2/sim22-pi.tex",
		  hline.after = c(-1, seq(0, 20, 6)))

	print(xtable(scen22$rho, auto = TRUE,
				 caption = tab_captions$rho[scenario + c(0, 5)],
				 label = "t:sim22-rho",
				 digits = n_digits_max),
		  file = "tab-2/sim22-rho.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("& \\multicolumn{6}{c}{$G$}\\\\\n",
		  				  			"\\cmidrule(l){2-7}\n",
		  				  			"$k$ & 12 & 60 & 120 & 600 & 1200 & 6000\\\\\n")),
		  include.rownames = TRUE, include.colnames = FALSE)

	print(xtable(scen22$beta0, auto = TRUE,
				 caption = tab_captions$beta0[scenario + c(0, 5)],
				 label = "t:sim22-beta0",
				 digits = n_digits_max),
		  file = "tab-2/sim22-beta0.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("\\multicolumn{6}{c}{$G$}\\\\\n",
		  				  			"\\midrule\n",
		  				  			"12 & 60 & 120 & 600 & 1200 & 6000\\\\\n")),
		  include.rownames = FALSE, include.colnames = FALSE)
}


# Scenario 23 -------------------------------------------------------------

scenario <- 3
cases  <- c(1, 2, 3, 5)
xlab_k <- c(rep("K", 3), rep("K = 5", 2))

# Figures
{
	pdf("fig-2/sim23-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 5), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3], phi[4], phi[5])
	for (k in 1:5) {
		j <- k <= cases
		boxplot(scen23$phi[, j, k],
				main = captions[k],
				names = cases[j],
				xlab = xlab_k[j],
				ylim = range(scen23$phi, na.rm = TRUE))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)
	for (k in 2:4) {
	    boxplot(as.numeric(scen23$prob[, k, ]),
	        main = expression(pi[k]),
	        xlab = paste("K =", cases[k]))
	    abline(h = 1 / cases[k], col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen23$iter, names = cases, xlab = "K")
	barplot(colSums(scen23$times), names = cases, xlab = "K")
	boxplot(vmat23, names = cases, xlab = "K")
	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen23[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen23[[e]]),
			SE   = apply(scen23[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("$K$")
	} else if (length(dim(scen23[[e]])) == 3) {
		aux <- NULL
		e_tex <- ifelse(e == "phi", "\\phi_", "\\pi_")
		for (j in 1:5) {
			aux <- rbind(aux, data.frame(
				k    = "",
				case = cases,
				Mean = colMeans(scen23[[e]][, , j]),
				SE   = apply(scen23[[e]][, , j], 2, sd)
			))
		}
		aux[c(2, 7, 11, 16, 20), 1] <- sprintf("$%s%d$", e_tex, 1:5)
		colnames(aux)[1:2] <- c(sprintf("$%sk$", e_tex), "$K$")
		tables[[e]] <- aux
	}
}

# Tables
{
	print(xtable(tables$times, auto = TRUE,
				 caption = tab_captions$times[scenario + c(0, 5)],
				 label = "t:sim23-times",
				 digits = n_digits),
		  file = "tab-2/sim23-times.tex")

	print(xtable(tables$iter, auto = TRUE,
				 caption = tab_captions$iter[scenario + c(0, 5)],
				 label = "t:sim23-iter",
				 digits = n_digits),
		  file = "tab-2/sim23-iter.tex")

	print(xtable(tables$phi[complete.cases(tables$phi), ], auto = TRUE,
				 caption = tab_captions$phi[scenario + c(0, 5)],
				 label = "t:sim23-phi",
				 digits = n_digits),
		  file = "tab-2/sim23-phi.tex",
		  hline.after = c(-1, 0, 4, 7, 9, 10, 11))

	print(xtable(tables$prob[complete.cases(tables$prob), ], auto = TRUE,
				 caption = tab_captions$pi[scenario + c(0, 5)],
				 label = "t:sim23-pi",
				 digits = n_digits),
		  file = "tab-2/sim23-pi.tex",
		  hline.after = c(-1, 0, 4, 7, 9, 10, 11))

	print(xtable(scen23$rho, auto = TRUE,
				 caption = tab_captions$rho[scenario + c(0, 5)],
				 label = "t:sim23-rho",
				 digits = n_digits_max),
		  file = "tab-2/sim23-rho.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("& \\multicolumn{4}{c}{$K$}\\\\\n",
		  				  			"\\cmidrule(l){2-5}\n",
		  				  			"$k$ & 1 & 2 & 3 & 5\\\\\n")),
		  include.rownames = TRUE, include.colnames = FALSE)

	print(xtable(scen23$beta0, auto = TRUE,
				 caption = tab_captions$beta0[scenario + c(0, 5)],
				 label = "t:sim23-beta0",
				 digits = n_digits_max),
		  file = "tab-2/sim23-beta0.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("\\multicolumn{4}{c}{$K$}\\\\\n",
		  				  			"\\midrule\n",
		  				  			"1 & 2 & 3 & 5\\\\\n")),
		  include.rownames = FALSE, include.colnames = FALSE)
}


# Scenario 25 -------------------------------------------------------------

scenario <- 5
cases    <- 1:2

scen25$rho   <- scen25$rho_mad
scen25$beta0 <- scen25$beta0_mad


# Figures
{
	pdf("fig-2/sim25-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 3), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:3) {
		boxplot(scen25$phi[, , k],
				main = captions[k],
				xlab = "Case",
				names = cases,
				ylim = range(scen25$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2], pi[3])
	for (k in 1:3) {
		boxplot(scen25$prob[, , k],
				main = captions[k],
				xlab = "Case",
				names = cases,
				ylim = range(scen25$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen25$iter, xlab = "Case", names = cases)
	barplot(colSums(scen25$times), xlab = "Case", names = cases)
	boxplot(vmat25, xlab = "Case", names = cases)

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen25[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen25[[e]]),
			SE   = apply(scen25[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("Case")
	} else if (length(dim(scen25[[e]])) == 3) {
		aux <- NULL
		e_tex <- ifelse(e == "phi", "\\phi_", "\\pi_")
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = c(sprintf("$%s%d$", e_tex, j), ""),
				case = cases,
				Mean = colMeans(scen25[[e]][, , j]),
				SE   = apply(scen25[[e]][, , j], 2, sd)
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
				 label = "t:sim25-times",
				 digits = n_digits),
		  file = "tab-2/sim25-times.tex")

	print(xtable(tables$iter, auto = TRUE,
				 caption = tab_captions$iter[scenario + c(0, 5)],
				 label = "t:sim25-iter",
				 digits = n_digits),
		  file = "tab-2/sim25-iter.tex")

	print(xtable(tables$phi, auto = TRUE,
				 caption = tab_captions$phi[scenario + c(0, 5)],
				 label = "t:sim25-phi",
				 digits = n_digits),
		  file = "tab-2/sim25-phi.tex",
		  hline.after = c(-1, seq(0, 6, 2)))

	print(xtable(tables$prob, auto = TRUE,
				 caption = tab_captions$pi[scenario + c(0, 5)],
				 label = "t:sim25-pi",
				 digits = n_digits),
		  file = "tab-2/sim25-pi.tex",
		  hline.after = c(-1, seq(0, 6, 2)))

	print(xtable(scen25$rho, auto = TRUE,
				 caption = tab_captions$rho[scenario + c(0, 5)],
				 label = "t:sim25-rho",
				 digits = n_digits),
		  file = "tab-2/sim25-rho.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("& \\multicolumn{2}{c}{Case}\\\\\n",
		  				  			"\\cmidrule(l){2-3}\n",
		  				  			"$k$ & 1 & 2\\\\\n")),
		  include.rownames = TRUE, include.colnames = FALSE)

	print(xtable(scen25$beta0, auto = TRUE,
				 caption = tab_captions$beta0[scenario + c(0, 5)],
				 label = "t:sim25-beta0",
				 digits = n_digits),
		  file = "tab-2/sim25-beta0.tex",
		  add.to.row = list(pos = list(0, 0, 0),
		  				  command = c("\\multicolumn{2}{c}{Case}\\\\\n",
		  				  			"\\midrule\n",
		  				  			"1 & 2\\\\\n")),
		  include.rownames = FALSE, include.colnames = FALSE)
}


# Scenario 26 -------------------------------------------------------------

scenario <- 6
cases    <- 1:4

# Figures
{
	pdf("fig-2/sim26-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)
	par(mfrow = c(1, 2), mar = c(4, 2, 2, 0.1), cex = .6)

	captions <- expression(phi[1], phi[2], phi[3])
	for (k in 1:2) {
		boxplot(scen26$phi[, , k],
			main = captions[k],
			xlab = "Case",
			names = cases,
			ylim = range(scen26$phi))
		abline(h = .1, col = 2, lty = 2, lwd = 2)
	}

	captions <- expression(pi[1], pi[2])
	{
		k <- 1
		boxplot(scen26$prob[, , k],
			main = captions[k],
			xlab = "Case",
			names = cases,
			ylim = range(scen26$prob))
		abline(h = c(.50, .25, .10, .05), col = 2, lty = 2, lwd = 2)

		k <- 2
		boxplot(scen26$prob[, , k],
			main = captions[k],
			xlab = "Case",
			names = cases,
			ylim = range(scen26$prob))
		abline(h = 1 - c(.50, .25, .10, .05), col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen26$iter, xlab = "Case", names = cases)
	barplot(colSums(scen26$times), xlab = "Case", names = cases)
	boxplot(vmat26, xlab = "Case", names = cases)

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen26[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen26[[e]]),
			SE   = apply(scen26[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("Case")
	} else if (length(dim(scen26[[e]])) == 3) {
		aux <- NULL
		e_tex <- ifelse(e == "phi", "\\phi_", "\\pi_")
		for (j in 1:2) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", sprintf("$%s%d$", e_tex, j), ""), c(1, 1, 2)),
				case = cases,
				Mean = colMeans(scen26[[e]][, , j]),
				SE   = apply(scen26[[e]][, , j], 2, sd)
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
		label = "t:sim26-times",
		digits = n_digits),
		file = "tab-2/sim26-times.tex")

	print(xtable(tables$iter, auto = TRUE,
		caption = tab_captions$iter[scenario + c(0, 5)],
		label = "t:sim26-iter",
		digits = n_digits),
		file = "tab-2/sim26-iter.tex")

	print(xtable(tables$phi, auto = TRUE,
		caption = tab_captions$phi[scenario + c(0, 5)],
		label = "t:sim26-phi",
		digits = n_digits),
		file = "tab-2/sim26-phi.tex",
		hline.after = c(-1, seq(0, 8, 4)))

	print(xtable(tables$prob, auto = TRUE,
		caption = tab_captions$pi[scenario + c(0, 5)],
		label = "t:sim26-pi",
		digits = n_digits),
		file = "tab-2/sim26-pi.tex",
		hline.after = c(-1, seq(0, 8, 4)))

	print(xtable(scen26$rho, auto = TRUE,
		caption = tab_captions$rho[scenario + c(0, 5)],
		label = "t:sim26-rho",
		digits = n_digits_max),
		file = "tab-2/sim26-rho.tex",
		add.to.row = list(pos = list(0, 0, 0),
			command = c("& \\multicolumn{4}{c}{Case}\\\\\n",
				"\\cmidrule(l){2-5}\n",
				"$k$ & 1 & 2 & 3 & 4\\\\\n")),
		include.rownames = TRUE, include.colnames = FALSE)

	print(xtable(scen26$beta0, auto = TRUE,
		caption = tab_captions$beta0[scenario + c(0, 5)],
		label = "t:sim26-beta0",
		digits = n_digits),
		file = "tab-2/sim26-beta0.tex",
		add.to.row = list(pos = list(0, 0, 0),
			command = c("\\multicolumn{4}{c}{Case}\\\\\n",
				"\\midrule\n",
				"1 & 2 & 3 & 4\\\\\n")),
		include.rownames = FALSE, include.colnames = FALSE)
}


# Scenario 28 -------------------------------------------------------------

scenario <- 8
cases    <- 1:6

# Figures
{
	pdf("fig-2/sim28-%d.pdf", onefile = FALSE, width = 5.2, height = 2.5)

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
		boxplot(scen28$prob[, , k],
			main = captions[k],
			xlab = "Case",
			names = cases,
			ylim = range(scen28$prob))
		abline(h = 1/3, col = 2, lty = 2, lwd = 2)
	}

	par(mfrow = c(1, 1), mar = c(4, 2, 0.1, 0.1), cex = .6)
	boxplot(scen28$iter, xlab = "Case", names = cases)
	barplot(colSums(scen28$times), xlab = "Case", names = cases)
	boxplot(vmat28, xlab = "Case", names = cases)

	dev.off()
}

tables <- list()
for (e in estim) {
	if (length(dim(scen28[[e]])) == 2) {
		tables[[e]] <- data.frame(
			case = cases,
			Mean = colMeans(scen28[[e]]),
			SE   = apply(scen28[[e]], 2, sd)
		)
		colnames(tables[[e]])[1] <- c("Case")
	} else if (length(dim(scen28[[e]])) == 3) {
		aux <- NULL
		e_tex <- ifelse(e == "phi", "\\phi_", "\\pi_")
		for (j in 1:3) {
			aux <- rbind(aux, data.frame(
				k    = rep(c("", sprintf("$%s%d$", e_tex, j), ""), c(2, 1, 3)),
				case = cases,
				Mean = colMeans(scen28[[e]][, , j]),
				SE   = apply(scen28[[e]][, , j], 2, sd)
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
		label = "t:sim28-times",
		digits = n_digits),
		file = "tab-2/sim28-times.tex")

	print(xtable(tables$iter, auto = TRUE,
		caption = tab_captions$iter[scenario + c(0, 5)],
		label = "t:sim28-iter",
		digits = n_digits),
		file = "tab-2/sim28-iter.tex")

	print(xtable(tables$phi, auto = TRUE,
		caption = tab_captions$phi[scenario + c(0, 5)],
		label = "t:sim28-phi",
		digits = n_digits),
		file = "tab-2/sim28-phi.tex",
		hline.after = c(-1, seq(0, 18, 6)))

	print(xtable(tables$prob, auto = TRUE,
		caption = tab_captions$pi[scenario + c(0, 5)],
		label = "t:sim28-pi",
		digits = n_digits),
		file = "tab-2/sim28-pi.tex",
		hline.after = c(-1, seq(0, 18, 6)))

	print(xtable(scen28$rho, auto = TRUE,
		caption = tab_captions$rho[scenario + c(0, 5)],
		label = "t:sim28-rho",
		digits = n_digits_max),
		file = "tab-2/sim28-rho.tex",
		add.to.row = list(pos = list(0, 0, 0),
			command = c("& \\multicolumn{6}{c}{Case}\\\\\n",
				"\\cmidrule(l){2-7}\n",
				"$k$ & 1 & 2 & 3 & 4 & 5 & 6\\\\\n")),
		include.rownames = TRUE, include.colnames = FALSE)

	print(xtable(scen28$beta0, auto = TRUE,
		caption = tab_captions$beta0[scenario + c(0, 5)],
		label = "t:sim28-beta0",
		digits = n_digits),
		file = "tab-2/sim28-beta0.tex",
		add.to.row = list(pos = list(0, 0, 0),
			command = c("\\multicolumn{6}{c}{Case}\\\\\n",
				"\\midrule\n",
				"1 & 2 & 3 & 4 & 5 & 6\\\\\n")),
		include.rownames = FALSE, include.colnames = FALSE)
}


# Moving files ----

dir.create("Results", showWarnings = FALSE)
file.copy("tab-2", "Results", recursive = TRUE, copy.date = TRUE)
file.copy("fig-2", "Results", recursive = TRUE, copy.date = TRUE)
unlink("tab-2", recursive = TRUE)
unlink("fig-2", recursive = TRUE)
