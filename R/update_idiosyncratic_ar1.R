update_idiosyncratic_ar1 <- function(smooth, layout) {
	F <- smooth$F_smooth
	P <- smooth$P_smooth
	PPm <- smooth$PPm_smooth
	F0 <- as.vector(smooth$F_smooth_0)
	P0 <- smooth$P_smooth_0

	indices <- unname(layout$idio_current)

	T_obs <- nrow(F)
	n_states <- ncol(F)

	if (T_obs < 2L) {
			stop("At least two observations are required.")
	}

	if (
			length(dim(P)) != 3L ||
			!identical(dim(P), c(n_states, n_states, T_obs)) ||
			!identical(dim(PPm), c(n_states, n_states, T_obs)) ||
			length(F0) != n_states ||
			!identical(dim(P0), c(n_states, n_states))
	) {
			stop("Incompatible smoother dimensions.")
	}

	if (anyNA(indices) ||
			any(indices < 1L | indices > n_states) ||
			anyDuplicated(indices)) {
			stop("Invalid idiosyncratic state indices.")
	}

	n_errors <- length(indices)

	if (n_errors == 0L) {
			return(list(
					coefficients = numeric(),
					variances = numeric()
			))
	}

	S_yy <- numeric(n_errors)
	S_xx <- numeric(n_errors)
	S_yx <- numeric(n_errors)

	for (t in seq_len(T_obs)) {
			y <- F[t, indices]

			if (t == 1L) {
					x <- F0[indices]
					P_previous <- P0
			} else {
					x <- F[t - 1L, indices]
					P_previous <- P[, , t - 1L]
			}

			S_yy <- S_yy +
					y^2 +
					diag(P[indices, indices, t, drop = FALSE][, , 1L])

			S_xx <- S_xx +
					x^2 +
					diag(P_previous[indices, indices, drop = FALSE])

			S_yx <- S_yx +
					y * x +
					diag(PPm[indices, indices, t, drop = FALSE][, , 1L])
	}

	if (any(!is.finite(S_xx)) || any(S_xx <= 0)) {
			stop("Invalid idiosyncratic normal equations.")
	}

	rho <- S_yx / S_xx
	sigma2 <- (S_yy - rho * S_yx) / T_obs

	names(rho) <- names(layout$idio_current)
	names(sigma2) <- names(layout$idio_current)

	list(
			coefficients = rho,
			variances = sigma2
	)
}