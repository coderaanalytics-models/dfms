#' Construct block-factor state space matrices
#' 
#' Builds the factor transition matrix and expanded state innovation covariance from
#' a validated factor-state layout.
#' 
#' @param layout Factor component (layout$factor) from build_state_layout().
#' @param coefficients Named list of VAR coefficient arrays.
#' 	Each array has dimensions (r_b, r_b, p_b).
#' @param covariances Named list of factor innovation covariance matrices, each of
#' 	dimension (r_b, r_b).
#' 
#' @return A list containing A and Q.
#' 
#' @keywords internal
build_factor_matrices <- function(
	layout,
	coefficients,
	covariances
) {
	blocks <- names(layout$block_indices)
	
	if (
		!setequal(names(coefficients), blocks) ||
		!setequal(names(covariances), blocks)
	) {
		stop("Coefficient and covariance lists must match the factor blocks.")
	}
	
	n <- layout$n_states
	
	A <- matrix(0, n, n)
	Q <- matrix(0, n, n)
	
	for (block in blocks) {
		r <- layout$block_dims[[block]]
		L <- layout$block_lags[[block]]
		p <- layout$block_orders[[block]]
	
		lag_idx <- layout$lag_indices[[block]]
		current <- lag_idx[[1L]]
	
		coef_block <- coefficients[[block]]
		Sigma <- covariances[[block]]
	
		if (
			!is.matrix(Sigma) ||
			!identical(dim(Sigma), c(r, r)) ||
			any(!is.finite(Sigma))
		) {
			stop("Invalid covariance matrix for block: ", block)
		}
	
		if (!isTRUE(all.equal(Sigma, t(Sigma)))) {
			stop("Covariance matrix must be symmetric: ", block)
		}

		if (
			!is.array(coef_block) ||
			!identical(dim(coef_block), c(r, r, p)) ||
			any(!is.finite(coef_block))
		) {
			stop("Invalid VAR coefficients for block: ", block)
		}

		for (j in seq_len(p)) {
			A[current, lag_idx[[j]]] <- coef_block[, , j]
		}
	
		# Lag-shifting identities
		if (L > 1L) {
			for (j in 2:L) {
				A[lag_idx[[j]], lag_idx[[j - 1L]]] <- diag(r)
			}
		}
	
		# Factor innovations enter contemporaneous states only
		Q[current, current] <- Sigma
	}
	list(
		A = A,
		Q = Q
	)
}




#' Construct the factor measurement matrix
#' 
#' Builds the factor portion of the measurement matrix for numerical factor loadings
#' and precomputed loading-to-state mappings.
#' 
#' @param spec Output of parse_factor_spec().
#' @param layout Factor component (layout$factor) from build_state_layout().
#' @param loading_maps Output from build_loading_maps().
#' @param loadings Numeric matrix with one row per series and one column per factor.
#' 
#' @return Numeric measurement matrix with dimensions n_series by layout$n_states.
#' 
#' @keywords internal
build_factor_measurement <- function(
	spec,
	layout,
	loading_maps,
	loadings
) {
	series <- rownames(spec$loading_mask)
	factors <- spec$factor_names

	if (
		!is.matrix(loadings) ||
			!is.numeric(loadings) ||
			!identical(dim(loadings), dim(spec$loading_mask))
	) {
		stop("loadings must be a numeric matrix matching the dimensions of the loading mask.")
	}

	if (
		!identical(rownames(loadings), series) ||
			!identical(colnames(loadings), factors)
	) {
		stop("Loading matrix row and column names must match the factor specification.")
	}

	if (any(!is.finite(loadings))) {
		stop("Factor loadings must be finite.")
	}

	if (
		!is.list(loading_maps) ||
		is.null(names(loading_maps)) ||
		anyDuplicated(names(loading_maps)) ||
		!setequal(names(loading_maps), series)
	) {
		stop("loading_maps must contain one entry for every series.")
	}
	
	loading_maps <- loading_maps[series]

	factor_measurement <- matrix(
		0,
		nrow = length(series),
		ncol = layout$n_states,
		dimnames = list(series, NULL)
	)

	for (series_name in series) {
		map <- loading_maps[[series_name]]

		allowed <- rownames(map$indices)

		lambda <- loadings[
			series_name,
			allowed,
			drop = TRUE
		]

		lambda <- as.numeric(lambda)

		factor_measurement[
			series_name,
			as.vector(map$indices)
		] <- as.vector(
			outer(lambda, map$weights)
		)
	}
	
	factor_measurement
}





#' Construct complete block-factor state space matrices
#' 
#' Extends the factor transition, innovation covariance, and measurement matrices with
#' monthly and quarterly idiosyncratic AR(1) processes.
#' 
#' @param spec Factor specification from parse_factor_spec().
#' @param layout Complete layout from build_state_layout() with idio_ar1 = TRUE.
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
		!is.matrix(factor_matrices$A) ||
			!all(dim(factor_matrices$A) == c(nf, nf)) ||
			any(!is.finite(factor_matrices$A))
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

	A[fi, fi] <- factor_matrices$A
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
	for (series_name in layout$quarterly_vars) {
		idx <- layout$quarterly_indices[[series_name]]

		if (length(idx) != length(.mariano_murasawa_weights)) {
			stop(
				"Quarterly state dimension does not match Mariano-Murasawa weights ",
				"for series: ",
				series_name
			)
		}

		current <- idx[[1L]]

		A[current, current] <- quarterly_rho[[series_name]]
		Q[current, current] <- quarterly_variances[[series_name]]

		if (length(idx) > 1L) {
			for (j in 2:length(idx)) {
				A[idx[[j]], idx[[j - 1L]]] <- 1
			}
		}

		C[series_name, idx] <- .mariano_murasawa_weights
	}
		
	list(
		A = A,
		C = C,
		Q = Q,
		R = R
	)
}




#' Initialise the complete state space model
#' 
#' Computes the stationary initial state covariance from the discrete Lyapunov equation:
#'			P = A P A' + Q
#' 
#' @param A State transition matrix.
#' @param Q State innovation covariance matrix.
#' @param tolerance Numerical tolerance for stationarity.
#' 
#' @return A list containing F_0 and P_0.
#' 
#' @keywords internal
initialise_full_state <- function(
	A,
	Q,
	tolerance = 1e-8
) {
	if (
		!is.numeric(tolerance) ||
		length(tolerance) != 1L ||
		!is.finite(tolerance) ||
		tolerance <= 0
	) {
		stop("tolerance must be a positive finite scalar.")
	}

	if (
		!is.numeric(A) ||
		!is.matrix(A) ||
		nrow(A) != ncol(A) ||
		any(!is.finite(A))
	) {
		stop("A must be a finite square matrix.")
	}
	
	n <- nrow(A)
	
	if (
		!is.numeric(Q) ||
		!is.matrix(Q) ||
		!identical(dim(Q), c(n, n)) ||
		any(!is.finite(Q))
	) {
		stop("Q must be a finite covariance matrix.")
	}
	
	if (
		!isTRUE(all.equal(
			Q,
			t(Q),
			tolerance = tolerance,
			check.attributes = FALSE
		))
	) {
		stop("Q must be symmetric.")
	}
	
	eigenvalues <- eigen(A, only.values = TRUE)$values
	spectral_radius <- max(Mod(eigenvalues))
	
	if (spectral_radius >= 1 - tolerance) {
		stop(
			"State transition matrix is not sufficiently stationary. Spectral radius: ",
			signif(spectral_radius, 6)
		)
	}
	
	# Solve:
	#
	# vec(P) = (I - A %x% A)^(-1) vec(Q)
	# 
	# Using column-major vectorisation.
	
	# TODO replace with Schur-based discrete Lyapunov solver
	K <- diag(n * n) - kronecker(A, A)
 
	P_vec <- solve(K, as.vector(Q))
	
	P_0 <- matrix(
		P_vec,
		nrow = n,
		ncol = n
	)
	
	# Remove numerical asymmetry
	P_0 <- (P_0 + t(P_0)) / 2
	
	F_0 <- rep(0, n)
	
	list(
		F_0 = F_0,
		P_0 = P_0
	)
}
