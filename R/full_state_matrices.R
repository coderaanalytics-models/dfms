#' Construct complete block-factor state space matrices
#' 
#' Extends the factor transition, innovation covariance, and measurement matrices with
#' monthly and quarterly idiosyncratic AR(1) processes.
#' 
#' @param spec Factor specification from parse_factor_spec().
#' @param layout Complete layout from build_full_state_layout().
#' @param factor_matrices Output of build_factor_matrices().
#' @param factor_measurement Output of build_factor_measurement().
#' @param monthly_rho Named numeric vector of monthly AR(1) coefficients.
#' @param monthly_variances Named numeric vector of monthly innovation variances.
#' @param quarterly_rho Named numeric vector of quarterly AR(1) coefficients.
#' @param quarterly_variances Named numeric vector of quarterly innovation variances.
#' @param observation_noise Positive scalar used to construct R.
#' 
#' @return List containing A, C, Q, and R.
#' 
#' @keywords internal
build_full_state_matrices <- function(
	spec,
	layout,
	factor_matrices,
	factor_measurement,
	monthly_rho,
	monthly_variances,
	quarterly_rho,
	quarterly_variances,
	observation_noise = 1e-6
) {
	if (!isTRUE(layout$idio_ar1)) {
		stop("This builder currently requires idio_ar1 = TRUE.")
	}

	series <- rownames(spec$loading_mask)

	n <- layout$n_states
	n_series <- length(series)
	nf <- layout$n_factor_states

	# Validate factor matrices
	if (
		!is.matrix(factor_matrices$T) ||
			!all(dim(factor_matrices$T) == c(nf, nf)) ||
			any(!is.finite(factor_matrices$T))
	) {
		stop("Invalid factor transition matrix.")
	}

	if (
		!is.matrix(factor_matrices$Q) ||
			!all(dim(factor_matrices$Q) == c(nf, nf)) ||
			any(!is.finite(factor_matrices$Q))
	) {
		stop("Invalid factor innovation covariance.")
	}

	if (
		!is.matrix(factor_measurement) ||
		!all(dim(factor_measurement) == c(n_series, nf)) ||
		is.null(rownames(factor_measurement)) ||
		anyDuplicated(rownames(factor_measurement)) ||
		!setequal(rownames(factor_measurement), series) ||
		any(!is.finite(factor_measurement))
	) {
		stop("Invalid factor measurement matrix.")
	}

	factor_measurement <- factor_measurement[series, , drop = FALSE]

	# Validate named idiosyncratic parameters
	validate_parameters <- function(
		x,
		expected,
		name,
		variance = FALSE
	) {
		if (
			!is.numeric(x) ||
				!is.null(dim(x)) ||
				is.null(names(x)) ||
				anyDuplicated(names(x)) ||
				!setequal(names(x), expected) ||
				any(!is.finite(x))
		) {
			stop("Invalid ", name, ".")
		}

		x <- x[expected]

		if (variance && any(x <= 0)) {
			stop(name, " must contain positive variances.")
		}

		if (!variance && any(abs(x) >= 1)) {
			stop(name, " must contain stationary AR(1) coefficients.")
		}
		x
	}

	monthly_rho <- validate_parameters(
		monthly_rho,
		layout$monthly_vars,
		"monthly_rho"
	)

	monthly_variances <- validate_parameters(
		monthly_variances,
		layout$monthly_vars,
		"monthly_variances",
		variance = TRUE
	)

	quarterly_rho <- validate_parameters(
		quarterly_rho,
		layout$quarterly_vars,
		"quarterly_rho"
	)

	quarterly_variances <- validate_parameters(
		quarterly_variances,
		layout$quarterly_vars,
		"quarterly_variances",
		variance = TRUE
	)
	
	if (
		!is.numeric(observation_noise) ||
			length(observation_noise) != 1L ||
			!is.finite(observation_noise) ||
			observation_noise <= 0
	) {
		stop("observation_noise must be positive and finite.")
	}

	# Allocate complete matrices
	A <- matrix(0, n, n)
	Q <- matrix(0, n, n)
	C <- matrix(
		0,
		n_series,
		n,
		dimnames = list(series, NULL)
	)

	R <- diag(observation_noise, n_series)

	# Factor component
	fi <- layout$factor_indices

	A[fi, fi] <- factor_matrices$T
	Q[fi, fi] <- factor_matrices$Q
	C[, fi] <- factor_measurement

	# Monthly idiosyncratic components
	for (series_name in layout$monthly_vars) {
		idx <- layout$monthly_indices[[series_name]]

		A[idx, idx] <- monthly_rho[[series_name]]
		Q[idx, idx] <- monthly_variances[[series_name]]
		C[series_name, idx] <- 1
	}

	# Quarterly idiosyncratic components
	weights <- c(1, 2, 3, 2, 1)

	for (series_name in layout$quarterly_vars) {
		idx <- layout$quarterly_indices[[series_name]]

		# Current idiosyncratic state
		A[idx[1L], idx[1L]] <- quarterly_rho[[series_name]]

		# Deterministic lag shifting
		for (j in 2:5) {
			A[idx[j], idx[j - 1L]] <- 1
		}

		# Innovations enter the current state only
		Q[idx[1L], idx[1L]] <- quarterly_variances[[series_name]]

		# Mariano-Murasawa aggregation
		C[series_name, idx] <- weights
	}

	list(
		A = A,
		C = C,
		Q = Q,
		R = R
	)
}