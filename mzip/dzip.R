dzip <- function(x, phi, lambda, log = FALSE) {
	if (log == FALSE) {
		y <- ifelse(x == 0,
					phi + ((1 - phi) * exp(-lambda)),
					(1 - phi) * dpois(x, lambda))
	} else {
		y <- ifelse(x == 0,
					log(phi + ((1 - phi) * exp(-lambda))),
					log(1 - phi) + dpois(x, lambda, log = TRUE))
	}

	return(y)
}
