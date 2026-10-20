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
