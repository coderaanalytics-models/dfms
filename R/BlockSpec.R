#' Parse a Dynamic Factor Specification
#' 
#' Constructs a factor loading map and factor block specification.
#' This function does not estimate model parameters.
#'Use spaces to indent, not tabs. (lintr whitespace_linter) 
#' @keywords internal
parse_factor_spec <- function(
	series,
	factors,
	factor_multiplicities = 1L,
	factor_orders = 1L,
	quarterly_vars = character()
) {
	valid_int <- function(x, minimum = 1L) {
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
		global_factor_names <- paste0("factor.", seq_len(factors))
		factors <- setNames(rep(list(global_factor_names), length(series)), series)
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
					minimum = 1L
				))
		) {
			stop("Invalid factor_multiplicities specification.")
		}

		factor_multiplicities <- factor_multiplicities[base_factors]
	}

	# Expand factors according to multiplicity
	expanded_factor_names <- setNames(
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

	all_expanded_factor_names <- unlist(
		expanded_factor_names,
		use.names = FALSE
	)

	if (anyDuplicated(all_expanded_factor_names)) {
		stop("Expanded factor names must be unique.")
	}

	expanded_assignments <- lapply(
		factors,
		function(assignment) {
			unlist(
				expanded_factor_names[assignment],
				use.names = FALSE
			)
		}
	)

	# Normalise factor blocks and orders
	if (valid_int(factor_orders, 1L)) {
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
							minimum = 1L
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
											valid_int(x$order, 1L)
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
			unlist(expanded_factor_names[names], use.names = FALSE)
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

	# Enforce the Statsmodels identification/initialisation dimensionality requirement.
	if (length(factor_names) > length(monthly_vars)) {
		stop("Number of factors cannot exceed the number of monthly observed variables.")
	}

	list(
		factor_names = factor_names,
		loading_mask = loading_mask,
		blocks = blocks,
		block_orders = setNames(
			as.integer(orders),
			names(blocks)
		),
		factor_multiplicities = factor_multiplicities,
		quarterly_vars = quarterly_vars,
		monthly_vars = monthly_vars,
		n_factors = length(factor_names),
		n_blocks = length(blocks)
	)
}






#' Construct Complete Block-Factor State Layout
#' 
#' Determines the ordering and indexing of factor and idiosyncratic states for a
#' block-structured dynamic factor model.
#' 
#' Quarterly observations require five monthly state positions from
#' Mariano-Murasawa aggregation.
#' 
#' @param spec Factor specifications from parse_factor_spec().
#' @param idio_ar1 Bool for whether idiosyncratic AR(1) dynamics are enabled.
#' 
#' @return A list containing the complete state-space dimensions and indices.
#' 
#' @keywords internal
build_state_layout <- function(
	spec,
	idio_ar1 = FALSE
) {
	if (
		!is.list(spec) ||
		is.null(spec$blocks) ||
		is.null(spec$block_orders) ||
		is.null(spec$loading_mask)
	) {
		stop("Invalid factor specification")
	}

	if (
		!is.logical(idio_ar1) ||
		length(idio_ar1) != 1L ||
		is.na(idio_ar1)
	) {
		stop("idio_ar1 must be TRUE or FALSE")
	}

	series <- rownames(spec$loading_mask)

	if (
		!is.character(spec$quarterly_vars) ||
		anyNA(spec$quarterly_vars) ||
		anyDuplicated(spec$quarterly_vars) ||
		!all(spec$quarterly_vars %in% series)
	) {
		stop("Invalid quarterly variables")
	}

	has_quarterly <- length(spec$quarterly_vars) > 0L

	# Factor state layout
	blocks <- spec$blocks
	block_orders <- spec$block_orders

	if (!identical(names(blocks), names(block_orders))) {
		stop("Block names and orders must agree")
	}

	min_lags <- if (has_quarterly) 5L else 1L

	block_lags <- pmax(as.integer(block_orders), min_lags)
	names(block_lags) <- names(blocks)

	block_dims <- vapply(blocks, length, integer(1))
	state_dims <- block_dims * block_lags

	block_indices <- vector("list", length(blocks))
	lag_indices <- vector("list", length(blocks))
	
	names(block_indices) <- names(blocks)
	names(lag_indices) <- names(blocks)

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
				seq.int(start, start + r - 1L)
			}
		)

		names(lag_indices[[block]]) <- paste0("lag", seq_len(L) - 1L)

		current_indices <- c(
			current_indices,
			lag_indices[[block]][[1L]]
		)

		offset <- offset + r * L
	}

	factor_lag_indices <- setNames(
		vector("list", length(spec$factor_names)),
		spec$factor_names
	)

	for (block in names(blocks)) {
		factors <- blocks[[block]]
		lags <- lag_indices[[block]]

		for (j in seq_along(factors)) {
			factor_lag_indices[[factors[[j]]]] <- vapply(
				lags,
				function(idx) idx[[j]],
				integer(1)
			)
		}
	}

	names(current_indices) <- unlist(blocks, use.names = FALSE)

	n_factor_states <- offset

	factor_layout <- list(
		n_states = n_factor_states,
		n_factors = spec$n_factors,
		block_orders = block_orders,
		block_lags = block_lags,
		block_dims = block_dims,
		state_dims = state_dims,
		block_indices = block_indices,
		lag_indices = lag_indices,
		factor_lag_indices = factor_lag_indices,
		current_indices = current_indices,
		quarterly = has_quarterly
	)

	# Monthly idiosyncratic states
	n_monthly <- length(spec$monthly_vars)

	monthly_indices <- seq.int(
		offset + 1L,
		length.out = n_monthly
	)

	names(monthly_indices) <- spec$monthly_vars
	offset <- offset + n_monthly

	# Quarterly idiosyncratic states
	quarterly_indices <- vector("list", length(spec$quarterly_vars))
	names(quarterly_indices) <- spec$quarterly_vars

	for (series_name in spec$quarterly_vars) {
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
		lapply(quarterly_indices, function(idx) idx[-1L]),
		use.names = FALSE
	)

	idio_current <- c(
		monthly_indices,
		quarterly_current
	)

	idio_indices <- c(
		unname(monthly_indices),
		unlist(quarterly_indices, use.names = FALSE)
	)

	list(
		n_states = offset,
		n_factor_states = n_factor_states,
		n_monthly = n_monthly,
		n_quarterly = length(spec$quarterly_vars),
		factor = factor_layout,
		factor_indices = seq_len(n_factor_states),
		monthly_vars = spec$monthly_vars,
		quarterly_vars = spec$quarterly_vars,
		monthly_indices = monthly_indices,
		quarterly_indices = quarterly_indices,
		quarterly_current = quarterly_current,
		quarterly_lagged = quarterly_lagged,
		idio_current = idio_current,
		idio_indices = idio_indices,
		idio_ar1 = idio_ar1
	)
}




#' Build loading maps
#'
#' Translates factor restrictions and state-space indices into compact representations
#' for the EM loading update.
#'
#' @param spec Output from parse_factor_spec().
#' @param layout Output from build_state_layout().
#'
#' @return Named list of loading maps.
#'
#' @keywords internal
build_loading_maps <- function(spec, layout) {
		mask <- spec$loading_mask
		factor_lag_indices <- layout$factor$factor_lag_indices

		series <- rownames(mask)

		if (!setequal(series, names(layout$idio_current))) {
				stop("Series and idiosyncratic state names do not match.")
		}

		maps <- vector("list", length(series))
		names(maps) <- series

		for (series_name in series) {
				allowed <- colnames(mask)[mask[series_name, ]]

				is_quarterly <- series_name %in% layout$quarterly_vars

				weights <- if (is_quarterly) {
						.mariano_murasawa_weights
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
						available <- factor_lag_indices[[factor]]

						if (is.null(available)) {
							stop("Missing state layout for factor: ", factor)
						}

						if (length(available) < n_lags) {
								stop(
										"Insufficient factor lags for series: ",
										series_name
								)
						}

						indices[j, ] <- available[seq_len(n_lags)]
				}

				idio_indices <- if (is_quarterly) {
					unname(
        		layout$quarterly_indices[[series_name]]
					)
				} else {
					unname(
        		layout$monthly_indices[[series_name]]
					)
				}
				
				maps[[series_name]] <- list(
						n_loadings = length(allowed),
						indices = indices,
						weights = weights,
						idio_indices = idio_indices
				)
		}
		maps
}
