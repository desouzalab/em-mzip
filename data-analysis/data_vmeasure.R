data_vmeasure <- function(main_dir = NULL, Z) {
    main_dir <- paste0("~/scratch/", main_dir)

    Vmat <- list()
    cases <- list.files(main_dir)

    for (case in cases) {
        case_dir <- paste0(main_dir, "/", case)

        n <- list.files(case_dir)
        n <- regmatches(n, regexec("EM MZIP ([0-9]+) Z.dat", n))
        n <- as.numeric(sapply(n, `[`, 2))
        n <- sort(n[!is.na(n)])

        aux <- function(x, z0) {
            zh <- scan(sprintf("%s/EM MZIP %d Z.dat", case_dir, x), quiet = TRUE)
            zh <- matrix(zh, nrow = length(z0), byrow = TRUE)
            zh <- apply(zh, 1, which.max)
            return(v_measure(z0, zh))
        }

        Vmat[[case]] <- sapply(n, aux, Z)
    }

    j <- which.max(sapply(Vmat, length))
    nj <- length(Vmat[[j]])
    aux <- NULL

    for (case in cases) {
        aux <- c(aux, Vmat[[case]], rep(NA, nj - length(Vmat[[case]])))
    }

    Vmat <- array(aux, c(length(Vmat[[j]]), length(cases)))

    return(Vmat)
}
