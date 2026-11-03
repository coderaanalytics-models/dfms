#' Update the restricted monthly factor loadings
#' 
#' @param X Observation matrix, with missing values as NA.
#' @param smooth Output from SKFS().
#' @param factor_indices Named list mapping each monthly series to its permitted
#' 		contemporaneous factor-state indices.
#' @param idio_indices Named integer vector mapping each monthly series to its
#' 		idiosyncratic state index.
#' 
#' @return Named list of estimated loading vectors.
#' 
#' @keywords internal
update_monthly_loadings <- function(
	X,
	smooth,
	factor_indices,
	idio_indices
) {
	X <- as.matrix(X)

	F <- smooth$F_smooth
	P <- smooth$P_smooth

	if (nrow(X) != nrow(F)) {
		stop("Observation and smoother lengths do not match.")
	}

	loadings <- list()

	for (series in names(factor_indices)) {
		if (!series %in% colnames(X)) {
			stop("Series not found in X: ", series)
		}

		indices <- factor_indices[[series]]
		idio <- unname(idio_indices[[series]])

		if (is.null(idio) || length(idio) != 1L || is.na(idio)) {
			stop("Invalid idiosyncratic state index for: ", series)
		}

		observed <- which(is.finite(X[, series]))

		if (length(observed) == 0L) {
			stop("No observations available for: ", series)
		}

		k <- length(indices)

		if (k == 0L) {
			loadings[[series]] <- numeric(0)
			next
		}

		S_ff <- matrix(0, k, k)
		S_fy <- numeric(k)
		S_fu <- numeric(k)

		for (t in observed) {
			f <- as.numeric(F[t, indices])
			u <- F[t, idio]
			y <- X[t, series]

			P_ff <- P[indices, indices, t, drop = FALSE][, , 1]

			P_fu <- as.numeric(P[indices, idio, t])

			S_ff <- S_ff + tcrossprod(f) + P_ff

			S_fy <- S_fy + f * y

			S_fu <- S_fu + f * u + P_fu
		}

		lambda <- solve(
			S_ff,
			S_fy - S_fu
		)
		loadings[[series]] <- as.numeric(lambda)
	}
	loadings
}