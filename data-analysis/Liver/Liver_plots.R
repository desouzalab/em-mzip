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

Kmax <- 2:9
rdata_files <- list.files("RData")

group_names <- c(
    "B cell",           "Dendritic cell",
    "Endothelial cell", "Epithelial cell",
    "Erythroblast",     "Granulocyte",
    "Hepatocyte",       "Kuppfer cell",
    "Macrophage",       "Neutrophil",
    "T cell"
)

chosen_vf <- c("v21m1s2", "v21m3s2")
kstar <- list(4, 5)
names(kstar) <- chosen_vf

col_list <- list(
    c("#81B7FF"),
    c("#81B7FF", "#FF9288"),
    c("#81B7FF", "#FF9288", "#00D65C"),
    c("#FF83FF", "#FF9288", "#00D4FF", "#F0AD00"),
    c("#FF83FF", "#FF9288", "#00D4FF", "#F0AD00", "#FF80DE")
)

for (v in c(21)) {
    ver <- paste0("v", v)

    load(sprintf("Data/Liver %s.RData", ver))

    group <- Ztrue$Annotation
    group <- gsub("Epithelia ", "Epithelial ", group, fixed = TRUE)
    group <- gsub("hepatocyte", "Hepatocyte", group, fixed = TRUE)

    group_id <- numeric(length(group))
    for (k in seq_along(group)) {
        group_id[grep(group_names[k], group, fixed = TRUE)] <- k
    }

    vfiles <- grep(ver, rdata_files, value = TRUE)


    # First Step ----

    for (m in c(1, 3)) {
        mod <- paste0("m", m)
        vm <- paste0(ver, mod)
        mfiles <- grep(mod, vfiles, value = TRUE)
        minaic <- NULL

        pdf(sprintf("fig/%s-%%d.pdf", vm), onefile = FALSE)

        for (mf in mfiles) {

            load(sprintf("RData/%s", mf))
            K <- Kmax[1:length(smin)]

            # AIC, all estimates
            boxplot(
                mat$aic,
                names = K,
                xlab  = "K",
                ylab  = "AIC"
            )

            # AIC, best estimates
            plot(
                K, smin,
                type = "b",
                xlab = "K",
                ylab = "AIC",
                axes = FALSE,
                frame.plot = TRUE)
            axis(1, at = K, labels = formatC(K, format = "fg"))
            axis(2)
            Kstar <- elbow_finder(K, smin)
            if (length(Kstar) == 0) {
                Kstar <- c(K[1], smin[1])
            }
            points(Kstar[1], Kstar[2], pch = 16, col = "red")
            lines(c(Kstar[1], Kstar[1], 0), c(0, Kstar[2], Kstar[2]),
                  lty = 3, col = "red")
            cat(mf, "\t", Kstar[2], "\n")
        }

        dev.off()
    }


    # Second Step ----

    set.seed(8918221)
    y_tsne <- Rtsne::Rtsne(Y[!duplicated(Y), ])$Y

    vfiles <- grep(ver, chosen_vf, value = TRUE)

    ch_index <- numeric()

    for (vf in vfiles) {
        K <- kstar[[vf]]
        mycol <- col_list[[K]]
        vm <- substr(vf, 1, 5)
        load(sprintf("RData/%s.RData", vf))
        cluster <- paste(K, "clusters")

        zhat <- zhat_mat[[cluster]]
        zhat <- factor(zhat, levels = 1:K)


        # Confusion matrix ----

        confmat <- table(group_id, zhat)

        rownames(confmat) <- group_names
        print(
            xtable::xtable(
                confmat,
                auto    = TRUE,
                label   = sprintf("t:%s-confmat", vm),
                digits  = n_digits
            ),
            add.to.row = list(
                pos     = list(0, 0, 0),
                command = c(
                    paste0("& \\multicolumn{", K, "}{c}{$\\hat{Z}$}\\\\\n"),
                    paste0("\\cmidrule(l){2-", K + 1, "}\n"),
                    paste0("Day ", paste("&", 1:K, collapse = " "), "\\\\\n"))),
            include.rownames = TRUE,
            include.colnames = FALSE,
            file = sprintf("tab/%s-confmat.tex", vm)
        )


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
            labels_row      = group_names,
            angle_col       = 0,
            filename        = sprintf("fig/%s-6.png", vm)
        )


        # Clusters heat map ----

        Z <- param_mat$Z[[cluster]]
        zmax <- apply(Z, 1, max)
        zorder <- order(zhat, zmax, decreasing = c(FALSE, TRUE))

        annotation_row <- data.frame(
            Cluster    = zhat,
            Annotation = as.factor(group_names[group_id])
        )
        rownames(Y) <- paste0("V", 1:nrow(Y))
        rownames(annotation_row) <- rownames(Y)

        breaks <- seq(
            min(Y, na.rm = TRUE),
            quantile(Y, 0.95, na.rm = TRUE),
            length.out = 101
        )
        if (length(unique(breaks)) == 1) {
            breaks <- seq(
                min(Y, na.rm = TRUE),
                max(Y, na.rm = TRUE),
                length.out = 101
            )
        }

        gaps_row <- cumsum(table(zhat))
        gaps_row <- gaps_row[-length(gaps_row)]
        pheatmap(
            Y[zorder, ],
            breaks         = breaks,
            scale          = "none",
            cluster_rows   = FALSE,
            cluster_cols   = FALSE,
            show_rownames  = FALSE,
            show_colnames  = FALSE,
            annotation_row = annotation_row,
            gaps_row       = gaps_row,
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


        # Lambda / Beta0 and Rho / Mu ----

        if ("lambda" %in% names(param_mat)) {
            lambda <- param_mat$lambda[[cluster]]

            pheatmap(
                lambda,
                scale         = "column",
                cellheight    = 15,
                cluster_rows  = FALSE,
                cluster_cols  = FALSE,
                show_rownames = TRUE,
                labels_row    = paste0("k=", 1:K),
                filename      = sprintf("fig/%s-9.png", vm)
            )
        }

        if ("beta0" %in% names(param_mat) & "rho" %in% names(param_mat)) {
            lambda <- rbind(
                param_mat$beta0[[cluster]],
                param_mat$rho[[cluster]]
            )

            pheatmap(
                lambda,
                cellheight    = 15,
                cluster_rows  = FALSE,
                cluster_cols  = FALSE,
                show_rownames = TRUE,
                labels_row    = c("beta0", paste0("rho", 1:K)),
                gaps_row      = 1,
                filename      = sprintf("fig/%s-9.png", vm)
            )
            pheatmap(
                t(param_mat$beta0[[cluster]]),
                cellheight    = 15,
                cluster_rows  = FALSE,
                cluster_cols  = FALSE,
                legend_breaks = c(-40, -10),
                show_rownames = TRUE,
                labels_row    = "beta0",
                filename      = sprintf("fig/%s-9a.png", vm)
            )
            pheatmap(
                param_mat$rho[[cluster]],
                scale         = "column",
                cellheight    = 15,
                cluster_rows  = FALSE,
                cluster_cols  = FALSE,
                show_rownames = TRUE,
                labels_row    = paste0("rho", 1:K),
                filename      = sprintf("fig/%s-9b.png", vm)
            )
        }

        if ("mu" %in% names(param_mat)) {
            mu <- param_mat$mu[[cluster]]

            pheatmap(
                mu,
                scale          = "column",
                cellheight     = 15,
                cluster_rows   = FALSE,
                cluster_cols   = FALSE,
                show_rownames  = TRUE,
                labels_row     = paste0("k=", 1:K),
                filename       = sprintf("fig/%s-9.png", vm)
            )
        }


        # t-SNE plots ----

        pdf(sprintf("fig/%s-10.pdf", vm), onefile = FALSE)
        xlim <- range(y_tsne[, 1])
        xlim[2] <- xlim[2] + 0.3 * diff(xlim)
        plot(y_tsne,
             col  = mycol[zhat],
             pch  = group_id - 1,
             xlab = "t-SNE 1",
             ylab = "t-SNE 2",
             xlim = xlim
        )
        legend("topright",
               legend = c("Annotation", group_names, "", "Cluster", 1:K),
               col    = c(0, rep(c(1, 0), c(11, 2)), mycol),
               pch    = c(1, 0:10, 1, 1, rep(15, K)),
               bty    = "n"
        )
        dev.off()


        # Parameters ----

        prob_mat <- as.matrix(param_mat$prob[[cluster]])
        rownames(prob_mat) <- 1:K

        print(
            xtable::xtable(
                prob_mat,
                auto    = TRUE,
                label   = sprintf("t:%s-pi", vm),
                digits  = n_digits),
            add.to.row  = list(
                pos = list(0),
                command = c(paste0("$k$ & $\\hat\\pi_k$\\\\\n"))),
            include.rownames = TRUE,
            include.colnames = FALSE,
            file = sprintf("tab/%s-pi.tex", vm)
        )

        phi_mat <- as.matrix(param_mat$phi[[cluster]])
        rownames(phi_mat) <- 1:K

        print(
            xtable::xtable(
                phi_mat,
                auto = TRUE,
                label = sprintf("t:%s-phi", vm),
                digits = n_digits),
            add.to.row = list(
                pos = list(0),
                command = c(paste0("$k$ & $\\hat\\phi_k$\\\\\n"))),
            include.rownames = TRUE,
            include.colnames = FALSE,
            file = sprintf("tab/%s-phi.tex", vm)
        )

        if ("size" %in% names(param_mat)) {
            size_mat <- as.matrix(param_mat$size[[cluster]])
            rownames(size_mat) <- 1:K

            print(
                xtable::xtable(
                    size_mat,
                    auto = TRUE,
                    label = sprintf("t:%s-size", vm),
                    digits = n_digits),
                add.to.row = list(
                    pos = list(0),
                    command = c(paste0("$k$ & $\\hat\\alpha_k$\\\\\n"))),
                include.rownames = TRUE,
                include.colnames = FALSE,
                file = sprintf("tab/%s-size.tex", vm)
            )
        }


        # CH index and Probability cut-off ----

        ch_index <- c(ch_index, fpc::calinhara(Y, zhat, K))
    }

    print(ch_index)
}
