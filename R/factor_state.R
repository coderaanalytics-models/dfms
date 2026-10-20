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
		block_lags = block_lags,
		block_dims = block_dims,
		state_dims = state_dims,
		block_indices = block_indices,
		lag_indices = lag_indices,
		current_indices = current_indices,
		quarterly = quarterly
	)
}