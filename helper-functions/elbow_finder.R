# Elbow finder function
#
# de Souza, C. P. E., Andronescu, M., Masud, T., Kabeer, F., Biele, J.,
# Laks, E., Lai, D., Ye, P., Brimhall, J., Wang, B., et al. (2020).
# Epiclomal: probabilistic clustering of sparse single-cell DNA methylation
# data.
# PLoS computational biology, 16(9):e1008270.

elbow_finder <- function(x_values, y_values) {
    # Max values to create line
    max_x_x <- max(x_values)
    max_x_y <- y_values[which.max(x_values)]
    max_y_y <- max(y_values)
    max_y_x <- x_values[which.max(y_values)]
    max_df <- data.frame(
        x = c(max_y_x, max_x_x),
        y = c(max_y_y, max_x_y))

    # Creating straight line between the max values
    fit <- lm(max_df$y ~ max_df$x)

    # Distance from point to line
    distances <- numeric(length(x_values))
    for (i in seq_along(distances)) {
        distances[i] <- abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2)
    }

    # Max distance point
    x_max_dist <- x_values[which.max(distances)]
    y_max_dist <- y_values[which.max(distances)]

    return(c(x_max_dist, y_max_dist))
}
