test_that("factor block updates have correct dimensions", {
	ss <- make_test_full_state()

	matrices <- ss$matrices

	initial <- initialise_full_state(
			matrices$A,
			matrices$Q
	)

	simulated <- simulate_full_state(
			matrices,
			initial
	)

	X <- make_mixed_frequency_observations(
			simulated$observations
	)

	smooth <- SKFS(
			X = X,
			A = matrices$A,
			C = matrices$C,
			Q = matrices$Q,
			R = matrices$R,
			F_0 = initial$F_0,
			P_0 = initial$P_0,
			loglik = FALSE
	)

	updated <- update_factor_blocks(
			smooth,
			ss$layout
	)

	factor_layout <- ss$layout$factor

	for (block in names(factor_layout$block_dims)) {
			r <- factor_layout$block_dims[[block]]
			p <- factor_layout$block_orders[[block]]

			expect_equal(
					dim(updated$coefficients[[block]]),
					c(r, r, p)
			)

			expect_equal(
					dim(updated$covariances[[block]]),
					c(r, r)
			)
	}
})

test_that("updated factor covariances are symmetric", {
	ss <- make_test_full_state()

	matrices <- ss$matrices

	initial <- initialise_full_state(
			matrices$A,
			matrices$Q
	)

	simulated <- simulate_full_state(
			matrices,
			initial
	)

	X <- make_mixed_frequency_observations(
			simulated$observations
	)

	smooth <- SKFS(
			X = X,
			A = matrices$A,
			C = matrices$C,
			Q = matrices$Q,
			R = matrices$R,
			F_0 = initial$F_0,
			P_0 = initial$P_0,
			loglik = FALSE
	)

	updated <- update_factor_blocks(
			smooth,
			ss$layout
	)

	for (Sigma in updated$covariances) {
			expect_equal(
					Sigma,
					t(Sigma),
					tolerance = 1e-10
			)
	}
})

test_that("updated factor covariances are positive semidefinite", {
	ss <- make_test_full_state()

	matrices <- ss$matrices

	initial <- initialise_full_state(
			matrices$A,
			matrices$Q
	)

	simulated <- simulate_full_state(
			matrices,
			initial
	)

	X <- make_mixed_frequency_observations(
			simulated$observations
	)

	smooth <- SKFS(
			X = X,
			A = matrices$A,
			C = matrices$C,
			Q = matrices$Q,
			R = matrices$R,
			F_0 = initial$F_0,
			P_0 = initial$P_0,
			loglik = FALSE
	)

	updated <- update_factor_blocks(
			smooth,
			ss$layout
	)

	for (Sigma in updated$covariances) {
			eigenvalues <- eigen(
					Sigma,
					symmetric = TRUE,
					only.values = TRUE
			)$values

			expect_true(all(eigenvalues >= -1e-8))
	}
})

test_that("scalar VAR(1) update matches analytical EM solution", {
	# Three observations and one latent factor.
	# Include the initial state at t = 0.
	smooth <- list(
			F_smooth = matrix(c(1, 2, 3), ncol = 1),
			P_smooth = array(c(0.4, 0.5, 0.6), dim = c(1, 1, 3)),
			PPm_smooth = array(c(0.1, 0.2, 0.3), dim = c(1, 1, 3)),
			F_smooth_0 = 0.5,
			P_smooth_0 = matrix(0.2, 1, 1)
	)

	layout <- list(
			factor = list(
					lag_indices = list(Global = list(1L)),
					block_orders = list(Global = 1L)
			)
	)

	result <- update_factor_blocks(smooth, layout)

	# E[f_t^2]
	S_yy <- (1^2 + 0.4) +
			(2^2 + 0.5) +
			(3^2 + 0.6)

	# E[f_(t-1)^2]
	S_xx <- (0.5^2 + 0.2) +
			(1^2 + 0.4) +
			(2^2 + 0.5)

	# E[f_t f_(t-1)]
	S_yx <- (1 * 0.5 + 0.1) +
			(2 * 1 + 0.2) +
			(3 * 2 + 0.3)

	expected_phi <- S_yx / S_xx
	expected_q <- (S_yy - expected_phi * S_yx) / 3

	expect_equal(
			unname(result$coefficients$Global[1, 1, 1]),
			expected_phi,
			tolerance = 1e-12
	)

	expect_equal(
			unname(result$covariances$Global[1, 1]),
			expected_q,
			tolerance = 1e-12
	)
})

test_that("VAR(1) update agrees with direct smoother moment calculation", {
	ss <- make_test_full_state()

	matrices <- ss$matrices

	initial <- initialise_full_state(
			matrices$A,
			matrices$Q
	)

	simulated <- simulate_full_state(
			matrices,
			initial
	)

	X <- make_mixed_frequency_observations(
			simulated$observations
	)

	smooth <- SKFS(
			X = X,
			A = matrices$A,
			C = matrices$C,
			Q = matrices$Q,
			R = matrices$R,
			F_0 = initial$F_0,
			P_0 = initial$P_0,
			loglik = FALSE
	)

	result <- update_factor_blocks(smooth, ss$layout)

	# Select a block with VAR(1) dynamics.
	block_orders <- ss$layout$factor$block_orders
	block <- names(block_orders)[which(block_orders == 1L)[1L]]

	expect_false(is.na(block))

	indices <- ss$layout$factor$lag_indices[[block]][[1L]]

	T_obs <- nrow(smooth$F_smooth)
	r <- length(indices)

	S_yy <- matrix(0, r, r)
	S_xx <- matrix(0, r, r)
	S_yx <- matrix(0, r, r)

	for (t in seq_len(T_obs)) {
			y <- matrix(smooth$F_smooth[t, indices], ncol = 1)

			if (t == 1L) {
					x <- matrix(smooth$F_smooth_0[indices], ncol = 1)
					P_previous <- smooth$P_smooth_0
			} else {
					x <- matrix(smooth$F_smooth[t - 1L, indices], ncol = 1)
					P_previous <- smooth$P_smooth[, , t - 1L]
			}

			S_yy <- S_yy +
					tcrossprod(y) +
					smooth$P_smooth[indices, indices, t, drop = FALSE][, , 1]

			S_xx <- S_xx +
					tcrossprod(x) +
					P_previous[indices, indices, drop = FALSE]

			S_yx <- S_yx +
					y %*% t(x) +
					smooth$PPm_smooth[indices, indices, t, drop = FALSE][, , 1]
	}

	expected_phi <- t(solve(S_xx, t(S_yx)))

	expected_q <- (
			S_yy - expected_phi %*% t(S_yx)
	) / T_obs

	expected_q <- (expected_q + t(expected_q)) / 2

	expect_equal(
			unname(result$coefficients[[block]][, , 1]),
			unname(expected_phi),
			tolerance = 1e-10
	)

	expect_equal(
			unname(result$covariances[[block]]),
			unname(expected_q),
			tolerance = 1e-10
	)
})

