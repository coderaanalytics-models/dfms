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