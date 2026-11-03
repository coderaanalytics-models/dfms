make_idio_fixture <- function(T_obs = 30L, n = 9L) {
	set.seed(42)

	F <- matrix(
			rnorm(T_obs * n),
			nrow = T_obs
	)

	P <- array(0, c(n, n, T_obs))
	PPm <- array(0, c(n, n, T_obs))

	for (t in seq_len(T_obs)) {
			A <- matrix(rnorm(n * n), n, n)

			P[, , t] <- crossprod(A) / n + diag(0.1, n)
			PPm[, , t] <- matrix(rnorm(n * n, sd = 0.01), n)
	}

	list(
			F_smooth = F,
			P_smooth = P,
			PPm_smooth = PPm,
			F_smooth_0 = rnorm(n),
			P_smooth_0 = diag(n)
	)
}

make_idio_layout <- function() {
	list(
			idio_current = c(
					M1 = 2L,
					M2 = 3L,
					M3 = 4L,
					Q1 = 5L
			),
			quarterly_indices = list(
					Q1 = 5:9
			)
	)
}

test_that("idiosyncratic AR1 updates match expected moments", {
	smooth <- make_idio_fixture()
	layout <- make_idio_layout()

	actual <- update_idiosyncratic_ar1(smooth, layout)

	idx <- unname(layout$idio_current)
	T_obs <- nrow(smooth$F_smooth)

	Y <- smooth$F_smooth[, idx, drop = FALSE]

	X <- rbind(
			smooth$F_smooth_0[idx],
			smooth$F_smooth[-T_obs, idx, drop = FALSE]
	)

	S_yy <- colSums(Y^2)
	S_xx <- colSums(X^2)
	S_yx <- colSums(Y * X)

	for (t in seq_len(T_obs)) {
			P_prev <- if (t == 1L) {
					smooth$P_smooth_0
			} else {
					smooth$P_smooth[, , t - 1L]
			}

			S_yy <- S_yy +
					diag(smooth$P_smooth[idx, idx, t])

			S_xx <- S_xx +
					diag(P_prev[idx, idx])

			S_yx <- S_yx +
					diag(smooth$PPm_smooth[idx, idx, t])
	}

	expected_rho <- S_yx / S_xx
	expected_variance <- (S_yy - expected_rho * S_yx) / T_obs

	names(expected_rho) <- names(layout$idio_current)
	names(expected_variance) <- names(layout$idio_current)

	expect_equal(
			actual$coefficients,
			expected_rho,
			tolerance = 1e-12
	)

	expect_equal(
			actual$variances,
			expected_variance,
			tolerance = 1e-12
	)
})

test_that("quarterly idiosyncratic lags are not estimated separately", {
	smooth <- make_idio_fixture()
	layout <- make_idio_layout()

	actual <- update_idiosyncratic_ar1(smooth, layout)

	expect_length(actual$coefficients, 4L)
	expect_length(actual$variances, 4L)

	expect_named(
			actual$coefficients,
			c("M1", "M2", "M3", "Q1")
	)
})

