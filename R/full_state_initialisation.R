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