EMstepBlock <- function(X, state, layout, loading_maps) {
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

	factor_updates <- update_factor_blocks_cpp(
			F_smooth = smooth$F_smooth,
			P_smooth = smooth$P_smooth,
			PPm_smooth = smooth$PPm_smooth,
			F_smooth_0 = as.vector(smooth$F_smooth_0),
			P_smooth_0 = smooth$P_smooth_0,
			factor_layout = layout$factor
	)

	idio_updates <- update_idiosyncratic_ar1(
			smooth = smooth,
			layout = layout
	)

	factor_matrices <- build_factor_matrices(
		layout = layout$factor,
		coefficients = factor_updates$coefficients,
		covariances = factor_updates$covariances
	)

	loading_updates <- update_loadings_cpp(
		X = X,
		F_smooth = smooth$F_smooth,
		P_smooth = smooth$P_smooth,
		loading_maps = loading_maps
	)

	loadings <- assemble_loadings(
		updates = loading_updates,
		loading_mask = spec$loading_mask
	)

	list(
		smooth = smooth,
		factor_updates = factor_updates,
		idio_updates = idio_updates,
		factor_matrices = factor_matrices,
		loading_updates = loading_updates,
		loadings = loadings
	)
}


assemble_loadings <- function(updates, loading_mask) {
	loadings <- matrix(
		0,
		nrow = nrow(loading_mask),
		ncol = ncol(loading_mask),
		dimnames = dimnames(loading_mask)
	)

	for (i in seq_len(nrow(loading_mask))) {
		allowed <- which(loading_mask[i, ])

		if (length(updates[[i]]) != length(allowed)) {
			stop("Loading dimensions do not match the loading mask")
		}
		loadings[i, allowed] <- updates[[i]]
	}
	loadings
}