library(pheatmap)
library(xtable)
options(
    xtable.include.rownames = FALSE,
    xtable.booktabs = TRUE,
    xtable.caption.placement = "top",
    xtable.sanitize.text.function = function(x) {x}
)
n_digits     <- 5
n_digits_max <- 9

source("../../helper-functions/elbow_finder.R")

dir.create("fig", showWarnings = FALSE)
dir.create("tab", showWarnings = FALSE)

load("Data/Klein v6.RData")
group <- Ztrue

mycol <- c("#FF83FF", "#FF9288", "#81B7FF", "#00D65C")


# Captions ----

tab_index <- 1
tab_captions <- list(
    c("Confusion matrix for the result of the result from fitting the simple ZIP model, and for each cell taking the cluster with the highest estimated probability, compared with the day label from the data set.", "Data 1. Confusion matrix for the simple ZIP"),
    c("Estimates obtained for the $\\pi$ parameter in the simple ZIP model.", "Data 1. $\\pi$ estimates, simple ZIP"),
    c("Estimates obtained for the $\\phi$ parameter in the simple ZIP model.", "Data 1. $\\phi$ estimates, simple ZIP"),
    c("Confusion matrix for the result of the result from fitting the ZIP with size factor model, and for each cell taking the cluster with the highest estimated probability, compared with the day label from the data set.", "Data 1. Confusion matrix for the ZIP with size factor"),
    c("Estimates obtained for the $\\pi$ parameter in the ZIP with size factor model.", "Data 1. $\\pi$ estimates, ZIP with size factor"),
    c("Estimates obtained for the $\\phi$ parameter in the ZIP with size factor model.", "Data 1. $\\phi$ estimates, ZIP with size factor"),
    c("CH index and probability cot-off points obtained for both simple ZIP (Model 1) and ZIP with size factors (Model 2) on the first data set.", "Data 1. CH-index and cut-off probabilities")
)


# First Step ----

Kmax <- 2:9

rdata_files <- list.files("RData")
vfiles <- grep("v6", rdata_files, value = TRUE)

for (m in 1:2) {
    mod <- paste0("m", m)
    vm <- paste0("v6", mod)
    mfiles <- grep(mod, vfiles, value = TRUE)
    minaic <- NULL

    pdf(sprintf("fig/%s-%s.pdf", vm, "%d"), onefile = FALSE)

    for (mf in mfiles) {

        load(sprintf("RData/%s", mf))
        minaic <- rbind(minaic, c(smin, rep(NA, length(Kmax) - length(smin))))
        K <- Kmax[1:length(smin)]

        # AIC, all estimates
        boxplot(mat$aic,
                names = K,
                xlab = "K",
                ylab = "AIC")

        # AIC, best estimates
        plot(K, smin,
             type = "b",
             xlab = "K",
             ylab = "AIC",
             axes = FALSE,
             frame.plot = TRUE)
        axis(1, at = K, labels = formatC(K, format = "fg"))
        axis(2)
        Kstar <- elbow_finder(K, smin)
        points(Kstar[1], Kstar[2], pch = 16, col = "red")
        lines(c(Kstar[1], Kstar[1], 0), c(0, Kstar[2], Kstar[2]),
              lty = 3, col = "red")
        cat(mf, "\t", Kstar[2], "\n")
    }

    minaic <- apply(minaic, 2, min, na.rm = TRUE)
    K <- Kmax[1:length(minaic)]

    dev.off()
}


# Second Step ----

library(Rtsne)
set.seed(8918221)
y_tsne <- Rtsne(Y)$Y

vfiles <- c("v6m1s2", "v6m2s2")
kstar <- list(4, 4)

names(kstar) <- vfiles
ch_index <- numeric()

for (vf in vfiles) {
    K <- kstar[[vf]]
    vm <- substr(vf, 1, 4)
    load(sprintf("RData/%s.RData", vf))
    cluster <- paste(K, "clusters")

    zhat <- zhat_mat[[cluster]]
    zhat <- factor(zhat, levels = 1:K)


    # Confusion matrix ----

    confmat <- table(group, zhat)

    print(
        xtable::xtable(
            confmat,
            auto = TRUE,
            caption = tab_captions[[tab_index]],
            label = sprintf("t:%s-confmat", vm),
            digits = n_digits
        ),
        add.to.row = list(
            pos = list(0, 0, 0),
            command = c(
                paste0("& \\multicolumn{", K, "}{c}{$\\hat{Z}$}\\\\\n"),
                paste0("\\cmidrule(l){2-", K + 1, "}\n"),
                paste0("Day ", paste("&", 1:K, collapse = " "), "\\\\\n"))),
        include.rownames = TRUE,
        include.colnames = FALSE,
        file = sprintf("tab/%s-confmat.tex", vm)
    )
    tab_index <- tab_index + 1


    # Co-clustering heat map ----

    confmat <- round(confmat / rowSums(confmat) * 100, 2)

    pheatmap(
        confmat,
        color           = colorRampPalette(c("white", "forestgreen"))(100),
        cellwidth       = 100,
        cellheight      = 100,
        cluster_rows    = FALSE,
        cluster_cols    = FALSE,
        legend_breaks   = seq(0, 100, 25),
        fontsize        = 15,
        display_numbers = TRUE,
        number_format   = "%.2f",
        fontsize_number = 20,
        labels_row      = paste("Day", c(0, 4)),
        angle_col       = 0,
        filename        = sprintf("fig/%s-6.png", vm)
    )


    # Clusters heat map ----

    Z <- param_mat$Z[[cluster]]
    zmax <- apply(Z, 1, max)
    zorder <- order(zhat, zmax, decreasing = c(FALSE, TRUE))

    annotation_row <- data.frame(
        Cluster = zhat,
        Day     = as.factor(group)
    )
    rownames(Y) <- paste0("V", 1:nrow(Y))
    rownames(annotation_row) <- rownames(Y)

    Yt <- Y
    ythresh <- quantile(Yt, 0.95)
    Yt[Yt > ythresh] <- ythresh
    pheatmap(
        Yt[zorder, ],
        scale          = "none",
        cluster_rows   = FALSE,
        cluster_cols   = FALSE,
        show_rownames  = FALSE,
        annotation_row = annotation_row,
        gaps_row       = cumsum(table(zhat))[-4],
        filename       = sprintf("fig/%s-7.png", vm)
    )


    # Z matrix ----

    rownames(Z) <- paste0("V", 1:nrow(Z))
    pheatmap(
        Z[zorder, ],
        color           = colorRampPalette(c("white", "forestgreen"))(100),
        cellwidth       = 10,
        cluster_rows   = FALSE,
        cluster_cols   = FALSE,
        show_rownames  = FALSE,
        annotation_row = annotation_row,
        filename       = sprintf("fig/%s-8.png", vm)
    )


    # Lambda / Beta0 and Rho ----

    if ("lambda" %in% names(param_mat)) {
        lambda <- param_mat$lambda[[cluster]]

        pheatmap(
            lambda,
            cellheight     = 15,
            cluster_rows   = FALSE,
            cluster_cols   = FALSE,
            show_rownames  = TRUE,
            labels_row     = paste0("k=", 1:K),
            filename       = sprintf("fig/%s-9.png", vm)
        )
    }

    if ("beta0" %in% names(param_mat) & "rho" %in% names(param_mat)) {
        lambda <- rbind(
            param_mat$beta0[[cluster]],
            param_mat$rho[[cluster]]
        )

        pheatmap(
            lambda,
            cellheight     = 15,
            cluster_rows   = FALSE,
            cluster_cols   = FALSE,
            show_rownames  = TRUE,
            labels_row     = c("beta0", paste0("rho", 1:K)),
            gaps_row       = 1,
            filename       = sprintf("fig/%s-9.png", vm)
        )
        pheatmap(
            t(param_mat$beta0[[cluster]]),
            cellheight     = 15,
            cluster_rows   = FALSE,
            cluster_cols   = FALSE,
            legend_breaks  = c(-40, -10),
            show_rownames  = TRUE,
            labels_row     = "beta0",
            filename       = sprintf("fig/%s-9a.png", vm)
        )
        pheatmap(
            param_mat$rho[[cluster]],
            scale          = "column",
            cellheight     = 15,
            cluster_rows   = FALSE,
            cluster_cols   = FALSE,
            show_rownames  = TRUE,
            labels_row     = paste0("rho", 1:K),
            filename       = sprintf("fig/%s-9b.png", vm)
        )
    }


    # t-SNE plots ----

    pdf(sprintf("fig/%s-10.pdf", vm), onefile = FALSE)
    plot(y_tsne,
        col = mycol[zhat],
        pch = ifelse(group == 0, 1, 2),
        xlab = "t-SNE 1",
        ylab = "t-SNE 2"
    )
    legend("topright",
        legend = c("Day", 0, 4, "", "Cluster", 1, 2, 3, 4),
        col = c(0, 1, 1, 0, 0, mycol),
        pch = c(1, 1, 2, 1, 1, rep(15, 4)),
        bty = "n"
    )
    dev.off()


    # Parameters ----

    prob_mat <- as.matrix(param_mat$prob[[cluster]])
    rownames(prob_mat) <- 1:K

    print(
        xtable::xtable(
            prob_mat,
            auto = TRUE,
            caption = tab_captions[[tab_index]],
            label = sprintf("t:%s-pi", vm),
            digits = n_digits),
        add.to.row = list(
            pos = list(0),
            command = c(paste0("$k$ & $\\hat\\pi_k$\\\\\n"))),
        include.rownames = TRUE,
        include.colnames = FALSE,
        file = sprintf("tab/%s-pi.tex", vm)
    )
    tab_index <- tab_index + 1


    phi_mat <- as.matrix(param_mat$phi[[cluster]])
    rownames(phi_mat) <- 1:K

    print(
        xtable::xtable(
            phi_mat,
            auto = TRUE,
            caption = tab_captions[[tab_index]],
            label = sprintf("t:%s-phi", vm),
            digits = n_digits),
        add.to.row = list(
            pos = list(0),
            command = c(paste0("$k$ & $\\hat\\phi_k$\\\\\n"))),
        include.rownames = TRUE,
        include.colnames = FALSE,
        file = sprintf("tab/%s-phi.tex", vm)
    )
    tab_index <- tab_index + 1


    # CH index ----

    ch_index <- c(ch_index, fpc::calinhara(Y, zhat, K))
}

print(ch_index)
