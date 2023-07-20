mat <- list()


# Calculate AIC for each case ----

for (case in cases) {
    cdir <- paste0(vdir, "/", case)
    cfiles <- list.files(cdir)

    n <- regmatches(cfiles, regexec("([0-9]+) loglik", cfiles))
    n <- as.numeric(sapply(n, `[`, 2))
    n <- sort(n[!is.na(n)])

    mat$ll[[case]] <- sapply(n, function(x) {
        scan(sprintf("%s/%s", cdir,
            grep(paste("", x, "loglik"), cfiles, value = TRUE)),
            quiet = TRUE)
    })

    mat$ll[[case]] <- mat$ll[[case]][!is.na(mat$ll[[case]])]
    if (length(mat$ll[[case]]) == 0) {
        cases <- setdiff(cases, case)
    } else {
        K <- as.numeric(strsplit(case, " ")[[1]][1])

        mat$aic[[case]] <- -2 * mat$ll[[case]] + (K * (G + 2)) * 2
        mat$bic[[case]] <- -2 * mat$ll[[case]] + (K * (G + 2)) * log(N * G)
    }
}

estim <- c("ll", "aic", "bic")
for (e in estim) {
    j <- which.max(sapply(mat[[e]], length))
    nj <- length(mat[[e]][[j]])
    aux <- NULL

    for (case in cases) {
        aux <- c(aux,
            mat[[e]][[case]],
            rep(NA, nj - length(mat[[e]][[case]])))
    }

    mat[[e]] <- array(aux, c(length(mat[[e]][[j]]), length(cases)))
}

s <- apply(mat$aic, 2, which.min)
smin <- numeric(ncol(mat$aic))
for (j in seq_len(ncol(mat$aic))) {
    smin[j] <- mat$aic[s[j], j]
}


# Heatmap, confusion matrix, and estimates ----

zhat_mat  <- list()
param_mat <- list()
for (i in seq_along(cases)) {
    case <- cases[i]
    Kcase <- (2:K)[i]
    cdir <- sprintf("%s/%s", vdir, case)
    cfiles <- list.files(cdir)

    zhat <- scan(sprintf("%s/EM %d Z.dat", cdir, s[i]), quiet = TRUE)
    zhat <- matrix(zhat, nrow = N, byrow = TRUE)
    param_mat$Z[[case]] <- zhat

    zhat <- apply(zhat, 1, which.max)
    zhat_mat[[case]] <- zhat

    zhat <- scan(sprintf("%s/EM %d Z.dat", cdir, s[i]), quiet = TRUE)
    zhat <- matrix(zhat, nrow = N, byrow = TRUE)
    param_mat$Z[[case]] <- zhat


    phi <- scan(sprintf("%s/EM %d phi.dat", cdir, s[i]), quiet = TRUE)
    param_mat$phi[[case]] <- phi
    prob <- scan(sprintf("%s/EM %d prob.dat", cdir, s[i]), quiet = TRUE)
    param_mat$prob[[case]] <- prob

    flambda <- sprintf("%s/EM %d lambda.dat", cdir, s[i])
    if (file.exists(flambda)) {
        lambda <- scan(flambda, quiet = TRUE)
        lambda <- matrix(lambda, nrow = Kcase, byrow = TRUE)
        param_mat$lambda[[case]] <- lambda
    }

    frho <- sprintf("%s/EM %d rho.dat", cdir, s[i])
    if (file.exists(frho)) {
        rho <- scan(frho, quiet = TRUE)
        rho <- matrix(rho, nrow = Kcase, byrow = TRUE)
        param_mat$rho[[case]] <- rho
    }

    fbeta0 <- sprintf("%s/EM %d beta0.dat", cdir, s[i])
    if (file.exists(fbeta0)) {
        beta0 <- scan(fbeta0, quiet = TRUE)
        param_mat$beta0[[case]] <- beta0
    }

    fmu <- sprintf("%s/EM %d mu.dat", cdir, s[i])
    if (file.exists(fmu)) {
        mu <- scan(fmu, quiet = TRUE)
        mu <- matrix(mu, nrow = Kcase, byrow = TRUE)
        param_mat$mu[[case]] <- mu
    }

    fsize <- sprintf("%s/EM %d size.dat", cdir, s[i])
    if (file.exists(fsize)) {
        size <- scan(fsize, quiet = TRUE)
        param_mat$size[[case]] <- size
    }

    phi <- scan(sprintf("%s/EM %d phi.dat", cdir, s[i]), quiet = TRUE)
    param_mat$phi[[case]] <- phi
}
