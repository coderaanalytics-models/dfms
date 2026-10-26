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