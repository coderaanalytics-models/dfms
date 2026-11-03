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

	list(
			smooth = smooth,
			factor_updates = factor_updates,
			idio_updates = idio_updates
	)
}