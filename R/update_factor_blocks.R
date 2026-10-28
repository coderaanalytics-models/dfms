#' Update block-specific factor dynamics
#' 
#' Computes the EM M-step for the VAR coefficients and innovation covariances of each
#' factor block.
#' 
#' @param smooth Output from SKFS().
#' @param layout Complete state layout from build_full_state_layout().
#' 
#' @return A list containing coefficients and covariances for each block.
#' 
#' @keywords internal
update_factor_blocks <- function(
	smooth,
	layout
) {
	F <- smooth$F_smooth
	P <- smooth$P_smooth
	PPm <- smooth$PPm_smooth

	F0 <- smooth$F_smooth_0
	P0 <- smooth$P_smooth_0

	factor_layout <- layout$factor

	T_obs <- nrow(F)

	if (T_obs < 2L) {
		stop("At least two observations are required.")
	}

	coefficients <- list()
	covariances <- list()

	for (block in names(factor_layout$lag_indices)) {
		p <- factor_layout$block_orders[[block]]

		lag_indices <- factor_layout$lag_indices[[block]]

		# Current factor states
		current <- lag_indices[[1L]]

		# Predictor states at time t-1:
		# f_(t-1), ..., f_(t-p)
		predictors <- unlist(
			lag_indices[seq_len(p)],
			use.names = FALSE
		)

		r <- length(current)
		k <- length(predictors)

		S_yy <- matrix(0, r, r)
		S_xx <- matrix(0, k, k)
		S_yx <- matrix(0, r, k)

		for (t in seq_len(T_obs)) {
			f_t <- F[t, current]

			if (t == 1L) {
				f_previous <- F0
				P_previous <- P0
			} else {
				f_previous <- F[t - 1L, ]
				P_previous <- matrix(
						P[, , t - 1L],
						nrow = ncol(F),
						ncol = ncol(F)
				)
			}

			y <- matrix(f_t, ncol = 1L)
			x <- matrix(f_previous[predictors], ncol = 1L)

			# E[f_t f_t']
			S_yy <- S_yy +
				tcrossprod(y) +
				matrix(
					P[current, current, t],
					nrow = r,
					ncol = r
				)

			# E[z_(t-1) z_(t-1)']
			S_xx <- S_xx +
				tcrossprod(x) +
				P_previous[predictors, predictors, drop = FALSE]

			# E[f_t z_(t-1)']
			S_yx <- S_yx +
				y %*% t(x) +
				matrix(
					PPm[current, predictors, t],
					nrow = r,
					ncol = k
				)
		}

		# Solve the normal equations without explicitly calculating the inverse
		Phi <- t(solve(S_xx, t(S_yx)))

		Sigma <- (
			S_yy - Phi %*% t(S_yx)
		) / T_obs
		
		Sigma <- (Sigma + t(Sigma)) / 2

		# Convert coefficient matrix into an array with dimensions (r, r, p)
		Phi_array <- array(0, dim = c(r, r, p))

		for (lag in seq_len(p)) {
			cols <- seq.int(
				(lag - 1L) * r + 1L,
				lag * r
			)
			Phi_array[, , lag] <- Phi[, cols, drop = FALSE]
		}

		coefficients[[block]] <- Phi_array
		covariances[[block]] <- Sigma
	}
	list(
		coefficients = coefficients,
		covariances = covariances
	)
}