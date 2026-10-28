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