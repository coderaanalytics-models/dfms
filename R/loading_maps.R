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