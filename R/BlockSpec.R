#' Parse a Dynamic Factor Specification
#' 
#' Constructs a factor loading map and factor block specification.
#' This function does not estimate model parameters.
#' 
#' @keywords internal
parse_factor_spec <- function(
	series,
	factors,
	factor_multiplicities = 1L,
	factor_orders = 1L,
	quarterly_vars = character()
) {
	valid_int <- function(x, minimum = 0L) {
		is.numeric(x) &&
			length(x) == 1L &&
			!is.na(x) &&
			is.finite(x) &&
			x == floor(x) &&
			x >= minimum
	}

	# Validate observed variables
	if (
		!is.character(series) ||
			!length(series) ||
			anyNA(series) ||
			anyDuplicated(series) ||
			any(!nzchar(series))
	) {
		stop("series must contain unique, non-empty names.")
	}

	if (
		!is.character(quarterly_vars) ||
			anyNA(quarterly_vars) ||
			anyDuplicated(quarterly_vars) ||
			!all(quarterly_vars %in% series)
	) {
		stop("quarterly_vars must be a subset of series.")
	}

	monthly_vars <- setdiff(series, quarterly_vars)

	# Normalise factor assignments to a named list
	if (valid_int(factors, 1L)) {
		names_global <- paste0("factor.", seq_len(factors))
		factors <- setNames(rep(list(names_global), length(series)), series)
	} else if (is.character(factors) && length(factors)) {
		factors <- setNames(
			rep(list(factors), length(series)), series)
	}

	if (
		!is.list(factors) ||
			is.null(names(factors)) ||
			anyDuplicated(names(factors)) ||
			!setequal(names(factors), series)
	) {
		stop(
			"factors must be a positive integer, a character vector, or a named list assigning ",
			"factors to every series."
		)
	}

	factors <- factors[series]

	for (name in series) {
		assignment <- factors[[name]]

		if (
			!is.character(assignment) ||
				!length(assignment) ||
				anyNA(assignment) ||
				any(!nzchar(assignment)) ||
				anyDuplicated(assignment)
		) {
			stop("Invalid factor assignment for series: ", name)
		}
	}

	base_factors <- unique(unlist(factors, use.names = FALSE))

	# Normalise factor multiplicities
	if (valid_int(factor_multiplicities, 1L)) {
		factor_multiplicities <- setNames(
			rep(as.integer(factor_multiplicities),
				length(base_factors)),
			base_factors
		)
	} else {
		if (
			!is.numeric(factor_multiplicities) ||
				is.null(names(factor_multiplicities)) ||
				anyDuplicated(names(factor_multiplicities)) ||
				!setequal(names(factor_multiplicities), base_factors) ||
				!all(vapply(
					factor_multiplicities,
					valid_int,
					logical(1),
					minimum = 0L
				))
		) {
			stop("Invalid factor_multiplicities specification.")
		}

		factor_multiplicities <- factor_multiplicities[base_factors]
	}

	# Expand factors according to multiplicity
	expanded = setNames(
		lapply(base_factors, function(name) {
			n <- factor_multiplicities[[name]]

			if (n == 1L) {
				name
			} else {
				paste0(name, ".", seq_len(n))
			}
		}),
		base_factors
	)

	expanded_names <- unlist(expanded, use.names = FALSE)

	if (anyDuplicated(expanded_names)) {
		stop("Expanded factor names must be unique.")
	}

	expanded_assignments <- lapply(factors, function(assignment) {
		unlist(expanded[assignment], use.names = FALSE)
	})

	# Normalise factor blocks and orders
	if (valid_int(factor_orders, 0L)) {
			# Same AR order for every factor block
			block_definitions <- as.list(base_factors)
			names(block_definitions) <- base_factors

			orders <- rep(
					as.integer(factor_orders),
					length(base_factors)
			)

	} else if (
			is.numeric(factor_orders) &&
			!is.null(names(factor_orders))
	) {
			# Individual AR order for each factor
			if (
					anyDuplicated(names(factor_orders)) ||
					!setequal(names(factor_orders), base_factors) ||
					!all(vapply(
							factor_orders,
							valid_int,
							logical(1),
							minimum = 0L
					))
			) {
					stop("Invalid factor_orders specification.")
			}

			block_definitions <- as.list(base_factors)
			names(block_definitions) <- base_factors

			orders <- as.integer(factor_orders[base_factors])

	} else if (is.list(factor_orders)) {
			# Explicit joint factor blocks
			block_definitions <- factor_orders

			if (
					!length(block_definitions) ||
					is.null(names(block_definitions)) ||
					anyDuplicated(names(block_definitions)) ||
					any(!nzchar(names(block_definitions))) ||
					!all(vapply(
							block_definitions,
							function(x) {
									is.list(x) &&
											is.character(x$factors) &&
											length(x$factors) > 0L &&
											!anyNA(x$factors) &&
											!anyDuplicated(x$factors) &&
											valid_int(x$order, 0L)
							},
							logical(1)
					))
			) {
					stop("Invalid factor block specification.")
			}

			assigned <- unlist(
					lapply(block_definitions, `[[`, "factors"),
					use.names = FALSE
			)

			if (
					anyDuplicated(assigned) ||
					!setequal(assigned, base_factors)
			) {
					stop("Every factor must belong to exactly one block.")
			}

			orders <- vapply(
					block_definitions,
					function(x) as.integer(x$order),
					integer(1)
			)

			block_definitions <- lapply(
					block_definitions,
					`[[`,
					"factors"
			)

	} else {
			stop("Invalid factor_orders specification.")
	}

	# Expand block membership
	blocks <- lapply(
		block_definitions,
		function(names) {
			unlist(expanded[names], use.names = FALSE)
		}
	)

	# Preserve explicit block order. Otherwise, order by loading prevalence, then factor name.
	if (!is.list(factor_orders)) {
		loading_counts <- vapply(
			blocks,
			function(block) {
				mean(vapply(block, function(factor) {
					sum(vapply(
						expanded_assignments,
						function(x) factor %in% x,
						logical(1)
					))
				}, numeric(1)))
			},
			numeric(1)
		)

		ordering <- order(-loading_counts, names(blocks))

		blocks <- blocks[ordering]
		orders <- orders[ordering]
	}

	factor_names <- unlist(blocks, use.names = FALSE)

	# Construct loading map
	loading_mask <- matrix(
		FALSE,
		nrow = length(series),
		ncol = length(factor_names),
		dimnames = list(series, factor_names)
	)

	for (name in series) {
		loading_mask[name, expanded_assignments[[name]]] <- TRUE
	}

	# Construct factor block indices
	block_indices <- vector("list", length(blocks))

	offset <- 0L

	for (i in seq_along(blocks)) {
		n <- length(blocks[[i]])

		block_indices[[i]] <- seq.int(offset + 1L, offset + n)

		offset <- offset + n
	}

	names(block_indices) <- names(blocks)

	# Enforce the Statsmodels identification/initialisation dimensionality requirement.
	if (length(factor_names) > length(monthly_vars)) {
		stop("Number of factors cannot exceed the number of monthly observed variables.")
	}

	list(
		factor_names = factor_names,
		loading_mask = loading_mask,
		blocks = blocks,
		block_indices = block_indices,
		block_orders = setNames(
			as.integer(orders),
			names(blocks)
		),
		factor_multiplicities = factor_multiplicities,
		n_factors = length(factor_names),
		n_blocks = length(blocks)
	)
}


#' Construct factor state space layout
#' 
#' Determines the ordering and indexing of factor states for a block-structured dynamic
#' factor model.
#' 
#' The layout contains no estimated parameters.
#' 
#' @param spec A factor specification from parse_factor_spec().
#' @param quarterly Logical. Whether quarterly aggregation is used.
#' 
#' @return A list describing factor-state dimensions and indices.
#' 
#' @keywords internal
build_factor_state_layout <- function(spec, quarterly = FALSE) {
	if (
		!is.list(spec) ||
			is.null(spec$blocks) ||
			is.null(spec$block_orders)
	) {
	  stop("Invalid factor specification.")
	}
  
  if (
		!is.logical(quarterly) ||
			length(quarterly) != 1L ||
			is.na(quarterly)
	) {
		stop("quarterly must be TRUE or FALSE.")
  }
  
  blocks <- spec$blocks
  orders <- spec$block_orders
  
  if (!identical(names(blocks), names(orders))) {
		stop("Block names and orders must agree.")
  }
  
  # Quarterly aggregation requires five monthly factor positions
  min_lags <- if (quarterly) 5L else 1L
  
  block_lags <- pmax(as.integer(orders), min_lags)
  
  names(block_lags) <- names(blocks)
  
  block_dims <- vapply(
		blocks,
		length,
		integer(1)
	)
  
  state_dims <- block_dims * block_lags
  
  # State indices for each block
  block_indices <- vector("list", length(blocks))
  names(block_indices) <- names(blocks)
  
  # State indices for each lag within each block
  lag_indices <- vector("list", length(blocks))
  names(lag_indices) <- names(blocks)
  
  # Contemporaneous state index for each factor
  current_indices <- integer(0)
  
  offset <- 0L
  
  for (block in names(blocks)) {
		r <- block_dims[[block]]
		L <- block_lags[[block]]
	
		block_indices[[block]] <- seq.int(
			offset + 1L,
			offset + r * L
		)
	
		lag_indices[[block]] <- lapply(
			seq_len(L),
			function(lag) {
			  start <- offset + (lag - 1L) * r + 1L
			  seq.int(
					start,
					start + r - 1L
				)
			}
		)
	
		names(lag_indices[[block]]) <- paste0(
			"lag", seq_len(L) - 1L
		)
	
		current_indices <- c(
			current_indices,
			lag_indices[[block]][[1L]]
		)
	
		offset <- offset + r * L
  }
  
  names(current_indices) <- unlist(blocks, use.names = FALSE)
  
  list(
		n_states = offset,
		n_factors = spec$n_factors,
		block_orders = orders,
		block_lags = block_lags,
		block_dims = block_dims,
		state_dims = state_dims,
		block_indices = block_indices,
		lag_indices = lag_indices,
		current_indices = current_indices,
		quarterly = quarterly
	)
}


#' Construct complete state space layout
#' 
#' Extends a block-factor state layout with monthly and quarterly idiosyncratic states.
#' 
#' Monthly idiosyncratic components occupy one state per series.
#' Quarterly components occupy five consecutive states per series to support
#' Mariano-Murasawa aggregation.
#' 
#' @param spec Factor specification from parse_factor_spec().
#' @param factor_layout Output of build_factor_state_layout().
#' @param quarterly_vars Character vector of quarterly series.
#' @param idio_ar1 Logical. Whether idiosyncratic AR(1) dynamics are enabled.
#' 
#' @return A list containing complete state space indices.
#' 
#' @keywords internal
build_full_state_layout <- function(
	spec,
	factor_layout,
	quarterly_vars = character(),
	idio_ar1 = TRUE
) {
	if (
		!is.logical(idio_ar1) ||
			length(idio_ar1) != 1L ||
			is.na(idio_ar1)
	) {
		stop("idio_ar1 must be TRUE or FALSE.")
	}

	series <- rownames(spec$loading_mask)

	if (
		!is.character(quarterly_vars) ||
			anyNA(quarterly_vars) ||
			anyDuplicated(quarterly_vars) ||
			!all(quarterly_vars %in% series)
	) {
		stop("Invalid quarterly variables.")
	}

	if (
		length(quarterly_vars) > 0L &&
			!isTRUE(factor_layout$quarterly)
	) {
		stop("Quarterly variables require a factor layout with quarterly = TRUE.")
	}

	monthly_vars <- series[!series %in% quarterly_vars]

	# Factor states
	n_factor_states <- factor_layout$n_states

	factor_indices <- seq_len(n_factor_states)

	# Monthly idiosyncratic states
	n_monthly <- length(monthly_vars)

	monthly_indices <- seq.int(
		n_factor_states + 1L,
		length.out = n_monthly
	)

	names(monthly_indices) <- monthly_vars

	# Quarterly idiosyncratic states
	quarterly_indices <- vector("list", length(quarterly_vars))
	names(quarterly_indices) <- quarterly_vars

	offset <- n_factor_states + n_monthly

	for (series_name in quarterly_vars) {
		quarterly_indices[[series_name]] <- seq.int(
			offset + 1L,
			offset + 5L
		)

		names(quarterly_indices[[series_name]]) <- paste0("lag", 0:4)

		offset <- offset + 5L
	}

	quarterly_current <- vapply(
		quarterly_indices,
		function(idx) idx[[1L]],
		integer(1)
	)

	quarterly_lagged <- unlist(
		lapply(
			quarterly_indices,
			function(idx) idx[-1L]
		),
		use.names = FALSE
	)

	# All contemporaneous idiosyncratic states
	idio_current <- c(
		monthly_indices,
		quarterly_current
	)

	# All idiosyncratic states
	idio_indices <- c(
		unname(monthly_indices),
		unlist(quarterly_indices, use.names = FALSE)
	)

	list(
		n_states = offset,
		n_factor_states = n_factor_states,
		n_monthly = n_monthly,
		n_quarterly = length(quarterly_vars),
		factor = factor_layout,
		factor_indices = factor_indices,
		monthly_vars = monthly_vars,
		quarterly_vars = quarterly_vars,
		monthly_indices = monthly_indices,
		quarterly_indices = quarterly_indices,
		quarterly_current = quarterly_current,
		quarterly_lagged = quarterly_lagged,
		idio_current = idio_current,
		idio_indices = idio_indices,
		idio_ar1 = idio_ar1
	)
}


#' Construct block-factor state space matrices
#' 
#' Builds the factor transition matrix and expanded state innovation covariance from
#' a validated factor-state layout.
#' 
#' @param layout Output of build_factor_state_layout()
#' @param coefficients Nmaed list of VAR coefficient arrays.
#' 	Each array has dimensions (r_b, r_b, p_b).
#' @param covariances Named list of factor innovation covariance matrices, each of
#' 	dimension (r_b, r_b).
#' 
#' @return A list containing T and Q.
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
  
  T <- matrix(0, n, n)
  Q <- matrix(0, n, n)
  
  for (block in blocks) {
		r <- layout$block_dims[[block]]
		L <- layout$block_lags[[block]]
		p <- layout$block_orders[[block]]
	
		lag_idx <- layout$lag_indices[[block]]
		current <- lag_idx[[1L]]
	
		A <- coefficients[[block]]
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
			p > 0L &&
				(
					!is.array(A) ||
						!identical(dim(A), c(r, r, p)) ||
						any(!is.finite(A))
				)
		) {
		  stop("Invalid VAR coefficients for block: ", block)
		}
	
		if (p == 0L && !is.null(A)) {
		  stop("Order-zero blocks must have NULL coefficients: ", block)
		}
	
		# Top row of the companion matrix
		if (p > 0L) {
		  for (j in seq_len(p)) {
				T[current, lag_idx[[j]]] <- A[, , j]
		  }
		}
	
		# Lag-shifting identities
		if (L > 1L) {
		  for (j in 2:L) {
				T[lag_idx[[j]], lag_idx[[j - 1L]]] <- diag(r)
		  }
		}
	
		# Factor innovations enter contemporaneous states only
		Q[current, current] <- Sigma
  }
  list(
		T = T,
		Q = Q
	)
}


#' Construct the restricted factor measurement matrix
#' 
#' Builds the measurement matrix for a block-factor state space model, supporting
#' monthly observations and Mariano-Murasawa aggregation for quarterly observations.
#' 
#' The function only constructs the factor portion of the measurement equation.
#' Idiosyncratic states are handled separately.
#' 
#' @param spec Output of parse_factor_spec().
#' @param layout Output of build_factor_state_layout().
#' @param loadings Numeric matrix with one row per series and one column per factor.
#' @param quarterly_vars Character vector identifying quarterly variables.
#' 
#' @return Numeric measurement matrix with dimensions n_series by layout$n_states.
#' 
#' @keywords internal
build_factor_measurement <- function(
	spec,
	layout,
	loadings,
	quarterly_vars = character()
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

	if (!all(quarterly_vars %in% series)) {
		stop("Unknown quarterly variables.")
	}

	if (anyDuplicated(quarterly_vars)) {
		stop("quarterly_vars contains duplicates.")
	}

	if (
		length(quarterly_vars) > 0L &&
			!isTRUE(layout$quarterly)
	) {
		stop("Quarterly variables require a state layout with quarterly = TRUE.")
	}

	if (
		length(quarterly_vars) > 0L &&
			any(layout$block_lags < 5L)
	) {
		stop("Quarterly aggregation requires at least five retained factor lags.")
	}

	# Enforce structural loading restrictions
	restricted_loadings <- loadings
	restricted_loadings[!spec$loading_mask] <- 0

	Z <- matrix(
		0,
		nrow = length(series),
		ncol = layout$n_states,
		dimnames = list(series, NULL)
	)

	quarterly_weights <- c(1, 2, 3, 2, 1)

	for (i in seq_along(series)) {
		series_name <- series[[i]]
		is_quarterly <- series_name %in% quarterly_vars

		for (block in names(spec$blocks)) {
			block_factors <- spec$blocks[[block]]

			lambda <- restricted_loadings[
				series_name,
				block_factors,
				drop = TRUE
			]

			# Preserve vector dimensions for single-factor blocks
			lambda <- as.numeric(lambda)

			lag_indices <- layout$lag_indices[[block]]

			if (is_quarterly) {
				for (lag in seq_along(quarterly_weights)) {
					Z[i, lag_indices[[lag]]] <- quarterly_weights[[lag]] * lambda
				}
			} else {
				Z[i, lag_indices[[1L]]] <- lambda
			}
		}
	}
	Z
}


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


#' Initialise the complete state space model
#' 
#' Computes the stationary inital state covariance from the discrete Lyapunov equation:
#'			P - A P A' + Q
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
  if(
		!is.matrix(A) ||
			nrow(A) != ncol(A) ||
			any(!is.finite(A))
	) {
		stop("A must be a finite square matrix.")
  }
  
  n <- nrow(A)
  
  if (
		!is.matrix(Q) ||
			!all(dim(Q) == c(n, n)) ||
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
  
  if (
		length(tolerance) != 1L ||
			!is.finite(tolerance) ||
			tolerance <= 0
	) {
		stop("tolerance must be a positive finite scalar.")
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


#' Construct restricted measurement loading maps
#'
#' Translates factor restrictions and state-space indices into compact representations
#' for the EM loading update.
#'
#' @param spec Output from parse_factor_spec().
#' @param layout Output from build_full_state_layout().
#'
#' @return Named list of loading maps.
#'
#' @keywords internal
build_loading_maps <- function(spec, layout) {
		mask <- spec$loading_mask
		factor_layout <- layout$factor

		series <- rownames(mask)

		if (!setequal(series, names(layout$idio_current))) {
				stop("Series and idiosyncratic state names do not match.")
		}

		# Construct a factor-name -> lagged state indices mapping.
		factor_lags <- list()

		for (block in names(spec$blocks)) {
				factors <- spec$blocks[[block]]
				lags <- factor_layout$lag_indices[[block]]

				for (j in seq_along(factors)) {
						factor_lags[[factors[[j]]]] <- vapply(
								lags,
								function(idx) idx[[j]],
								integer(1)
						)
				}
		}

		maps <- vector("list", length(series))
		names(maps) <- series

		for (series_name in series) {
				allowed <- colnames(mask)[mask[series_name, ]]

				quarterly <- series_name %in% layout$quarterly_vars

				weights <- if (quarterly) {
						c(1, 2, 3, 2, 1)
				} else {
						1
				}

				n_lags <- length(weights)

				indices <- matrix(
						0L,
						nrow = length(allowed),
						ncol = n_lags,
						dimnames = list(
								allowed,
								paste0("lag", seq_len(n_lags) - 1L)
						)
				)

				for (j in seq_along(allowed)) {
						factor <- allowed[[j]]
						available <- factor_lags[[factor]]

						if (length(available) < n_lags) {
								stop(
										"Insufficient factor lags for series: ",
										series_name
								)
						}

						indices[j, ] <- available[seq_len(n_lags)]
				}

				maps[[series_name]] <- list(
						n_loadings = length(allowed),
						indices = indices,
						weights = weights,
						idio_index = unname(
								layout$idio_current[[series_name]]
						)
				)
		}
		maps
}