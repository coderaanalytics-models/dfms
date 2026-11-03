#' Update restricted quarterly factor loadings
#' 
#' @param X Observation matrix with missing values as NA.
#' @param smooth Output from SKFS().
#' @param factor_indices Named list mapping each quarterly series to a list of five
#' 		contemporaneous and lagged factor-state index vectors.
#' @param idio_indices Named integer vector mapping each quarterly series to its current
#' 		idiosyncratic state index.
#' 
#' @return Named list of estimated loading factors.
#' 
#' @keywords internal
update_quarterly_loadings <- function(
	X,
	smooth,
	factor_indices,
	idio_indices
) {
	weights <- c(1, 2, 3, 2, 1)

	result <- list()

	for (series in names(factor_indices)) {
			indices <- factor_indices[[series]]
			idio <- idio_indices[[series]]

			# Each element of indices contains the state indices
			# corresponding to one aggregation lag.
			index_matrix <- do.call(cbind, indices)

			k <- nrow(index_matrix)

			S_gg <- matrix(0, k, k)
			S_gy <- numeric(k)
			S_gu <- numeric(k)

			observed <- which(!is.na(X[, series]))

			for (t in observed) {
					mu <- smooth$F_smooth[t, ]
					P <- smooth$P_smooth[, , t]

					g <- as.vector(index_matrix %*% weights)

					# Posterior mean of aggregated factors.
					g_mean <- vapply(
							seq_len(k),
							function(j) {
									sum(weights * mu[index_matrix[j, ]])
							},
							numeric(1)
					)

					# Posterior covariance of aggregated factors.
					g_cov <- matrix(0, k, k)

					for (a in seq_along(weights)) {
							for (b in seq_along(weights)) {
									g_cov <- g_cov +
											weights[a] * weights[b] *
											P[
													index_matrix[, a],
													index_matrix[, b],
													drop = FALSE
											]
							}
					}

					# Covariance between aggregated factors and idiosyncratic state.
					g_u_cov <- numeric(k)

					for (a in seq_along(weights)) {
							g_u_cov <- g_u_cov +
									weights[a] * P[index_matrix[, a], idio]
					}

					S_gg <- S_gg + tcrossprod(g_mean) + g_cov
					S_gy <- S_gy + X[t, series] * g_mean
					S_gu <- S_gu + g_mean * mu[idio] + g_u_cov
			}

			result[[series]] <- as.vector(
					solve(S_gg, S_gy - S_gu)
			)
	}

	return(result)
}