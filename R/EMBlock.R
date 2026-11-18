#' Perform one EM step for a block dynamic factor model
#'
#' Runs the Kalman smoother using the current state-space parameters, updates
#' factor dynamics, idiosyncratic AR(1) parameters, and factor loadings, then
#' reconstructs the state-space matrices for the next iteration.
#'
#' @param X Numeric observation matrix.
#' @param state Current state-space parameters containing A, C, Q, R, F_0,
#'   and P_0.
#' @param spec Output from parse_factor_spec().
#' @param layout Output from build_state_layout().
#' @param loading_maps Output from build_loading_maps().
#'
#' @return A list containing the updated state, log-likelihood, smoother
#'   output, and parameter updates.
#'
#' @keywords internal
EMstepBlock <- function(
	X,
	state,
	spec,
	layout,
	loading_maps
) {
	# E-step
	smooth <- dfms::SKFS(
		X = X,
		A = state$A,
		C = state$C,
		Q = state$Q,
		R = state$R,
		F_0 = state$F_0,
		P_0 = state$P_0,
		loglik = TRUE
	)

	# M-step: update model parameters
	factor_updates <- update_factor_blocks_cpp(
		F_smooth = smooth$F_smooth,
		P_smooth = smooth$P_smooth,
		PPm_smooth = smooth$PPm_smooth,
		F_smooth_0 = as.vector(smooth$F_smooth_0),
		P_smooth_0 = smooth$P_smooth_0,
		factor_layout = layout$factor
	)

	idio_updates <- update_idiosyncratic_ar1_cpp(
		F_smooth = smooth$F_smooth,
		P_smooth = smooth$P_smooth,
		PPm_smooth = smooth$PPm_smooth,
		F_smooth_0 = as.vector(smooth$F_smooth_0),
		P_smooth_0 = smooth$P_smooth_0,
		idio_indices = layout$idio_current
	)

	loading_updates <- update_loadings_cpp(
		X = X,
		F_smooth = smooth$F_smooth,
		P_smooth = smooth$P_smooth,
		loading_maps = loading_maps
	)

	# Reconstruct factor components
	factor_matrices <- build_factor_matrices(
		layout = layout$factor,
		coefficients = factor_updates$coefficients,
		covariances = factor_updates$covariances
	)

	loadings <- assemble_loadings(
		updates = loading_updates,
		loading_mask = spec$loading_mask
	)

	factor_measurement <- build_factor_measurement(
		spec = spec,
		layout = layout$factor,
		loading_maps = loading_maps,
		loadings = loadings
	)

	# Split idiosyncratic parameters by frequency
	monthly_rho <- idio_updates$coefficients[
		layout$monthly_vars
	]

	monthly_variances <- idio_updates$variances[
		layout$monthly_vars
	]

	quarterly_rho <- idio_updates$coefficients[
		layout$quarterly_vars
	]

	quarterly_variances <- idio_updates$variances[
		layout$quarterly_vars
	]

	# Preserve the fixed observation-noise scale
	observation_noise <- state$R[[1L, 1L]]

	# Rebuild A, C, Q, and R
	state_matrices <- build_full_state_matrices(
		spec = spec,
		layout = layout,
		factor_matrices = factor_matrices,
		factor_measurement = factor_measurement,
		monthly_rho = monthly_rho,
		monthly_variances = monthly_variances,
		quarterly_rho = quarterly_rho,
		quarterly_variances = quarterly_variances,
		observation_noise = observation_noise
	)

	# Keep the initial state moments fixed across EM iterations.
	#
	# F_0 and P_0 are constructed once by initialise_block_state() from the
	# stationary initial state representation. They are not re-estimated in the
	# block M-step.
	new_state <- c(
		state_matrices,
		list(
			F_0 = state$F_0,
			P_0 = state$P_0
		)
	)

	list(
		state = new_state,
		loglik = smooth$loglik,
		smooth = smooth,
		factor_updates = factor_updates,
		idio_updates = idio_updates,
		loading_updates = loading_updates,
		loadings = loadings
	)
}


#' Run EM estimation for a block dynamic factor model
#'
#' Repeatedly applies EMstepBlock() until the log-likelihood converges or the
#' maximum number of iterations is reached.
#'
#' @param X Numeric observation matrix.
#' @param state Initial state-space parameters.
#' @param spec Output from parse_factor_spec().
#' @param layout Output from build_state_layout().
#' @param loading_maps Output from build_loading_maps().
#' @param min.iter Minimum number of EM iterations.
#' @param max.iter Maximum number of EM iterations.
#' @param tol Log-likelihood convergence tolerance.
#' @param check.increased Logical indicating whether likelihood decreases
#'   should prevent convergence.
#'
#' @return A list containing the final state, smoother output, parameter
#'   estimates, log-likelihood history, iteration count, and convergence status.
#'
#' @keywords internal
EMBlock <- function(
	X,
	state,
	spec,
	layout,
	loading_maps,
	min.iter = 25L,
	max.iter = 100L,
	tol = 1e-4,
	check.increased = FALSE
) {
	previous_loglik <- -Inf
	loglik_all <- numeric(max.iter)

	num_iter <- 0L
	converged <- FALSE
	step <- NULL

	while (num_iter < max.iter && !converged) {
		step <- EMstepBlock(
			X = X,
			state = state,
			spec = spec,
			layout = layout,
			loading_maps = loading_maps
		)

		loglik <- step$loglik
		state <- step$state

		num_iter <- num_iter + 1L
		loglik_all[[num_iter]] <- loglik

		if (
			num_iter < min.iter ||
			!is.finite(previous_loglik)
		) {
			converged <- FALSE
		} else {
			convergence <- em_converged(
				loglik,
				previous_loglik,
				tol,
				check.increased
			)

			if (check.increased) {
				converged <- convergence[[1L]] &&
					!convergence[[2L]]
			} else {
				converged <- convergence
			}
		}

		previous_loglik <- loglik
	}

	# Run the smoother once more using the final updated parameters.
	smooth <- dfms::SKFS(
		X = X,
		A = state$A,
		C = state$C,
		Q = state$Q,
		R = state$R,
		F_0 = state$F_0,
		P_0 = state$P_0,
		loglik = TRUE
	)

	if (converged) {
		message(
			"Converged after ",
			num_iter,
			" iterations."
		)
	} else {
		warning(
			"Maximum number of iterations reached without convergence."
		)
	}

	list(
		state = state,
		smooth = smooth,
		loadings = step$loadings,
		factor_updates = step$factor_updates,
		idio_updates = step$idio_updates,
		loglik = loglik_all[seq_len(num_iter)],
		iterations = num_iter,
		converged = converged
	)
}




#' Assemble the complete factor loading matrix
#'
#' Expands the per-series loading estimates returned by update_loadings_cpp()
#' into a complete loading matrix according to the structural loading mask.
#'
#' @param updates Named list of loading vectors.
#' @param loading_mask Logical loading restriction matrix.
#'
#' @return Numeric loading matrix with the same dimensions and names as
#'   loading_mask.
#'
#' @keywords internal
assemble_loadings <- function(
	updates,
	loading_mask
) {
	series <- rownames(loading_mask)

	if (
		!is.list(updates) ||
		is.null(names(updates)) ||
		anyDuplicated(names(updates)) ||
		!setequal(names(updates), series)
	) {
		stop(
			"updates must contain one loading vector for every series."
		)
	}

	updates <- updates[series]

	loadings <- matrix(
		0,
		nrow = nrow(loading_mask),
		ncol = ncol(loading_mask),
		dimnames = dimnames(loading_mask)
	)

	for (series_name in series) {
		allowed <- which(
			loading_mask[series_name, ]
		)

		if (length(updates[[series_name]]) != length(allowed)) {
			stop(
				"Loading dimensions do not match the loading mask for series: ",
				series_name
			)
		}

		loadings[
			series_name,
			allowed
		] <- updates[[series_name]]
	}

	loadings
}
