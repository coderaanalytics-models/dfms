make_factor_smoother_fixture <- function(
	T_obs = 20L,
	n = 6L,
	seed = 42L
) {
	set.seed(seed)

	F_smooth <- matrix(
			rnorm(T_obs * n),
			nrow = T_obs,
			ncol = n
	)

	F_smooth_0 <- rnorm(n)

	P_smooth <- array(0, dim = c(n, n, T_obs))
	PPm_smooth <- array(0, dim = c(n, n, T_obs))

	for (t in seq_len(T_obs)) {
			A <- matrix(rnorm(n * n), n, n)

			P_smooth[, , t] <- crossprod(A) / n + diag(0.1, n)

			PPm_smooth[, , t] <- matrix(
					rnorm(n * n, sd = 0.01),
					n,
					n
			)
	}

	A0 <- matrix(rnorm(n * n), n, n)

	P_smooth_0 <- crossprod(A0) / n + diag(0.1, n)

	list(
			F_smooth = F_smooth,
			P_smooth = P_smooth,
			PPm_smooth = PPm_smooth,
			F_smooth_0 = F_smooth_0,
			P_smooth_0 = P_smooth_0
	)
}

make_factor_layout_fixture <- function() {
	list(
			lag_indices = list(
					Global = list(
							c(1L, 2L),
							c(3L, 4L)
					),
					Real = list(
							5L,
							6L
					)
			),
			block_orders = list(
					Global = 2L,
					Real = 2L
			)
	)
}

test_that("C++ factor updates match the R reference", {
	smooth <- make_factor_smoother_fixture()
	factor_layout <- make_factor_layout_fixture()

	reference <- update_factor_blocks(
			smooth = smooth,
			layout = list(factor = factor_layout)
	)

	actual <- update_factor_blocks_cpp(
			F_smooth = smooth$F_smooth,
			P_smooth = smooth$P_smooth,
			PPm_smooth = smooth$PPm_smooth,
			F_smooth_0 = smooth$F_smooth_0,
			P_smooth_0 = smooth$P_smooth_0,
			factor_layout = factor_layout
	)

	expect_named(
			actual,
			c("coefficients", "covariances")
	)

	expect_named(
			actual$coefficients,
			c("Global", "Real")
	)

	expect_named(
			actual$covariances,
			c("Global", "Real")
	)

	for (block in names(factor_layout$lag_indices)) {
			expect_equal(
					actual$coefficients[[block]],
					reference$coefficients[[block]],
					tolerance = 1e-10
			)

			expect_equal(
					actual$covariances[[block]],
					reference$covariances[[block]],
					tolerance = 1e-10
			)
	}
})

test_that("C++ factor updates support VAR(1) blocks", {
	smooth <- make_factor_smoother_fixture(
			T_obs = 30L,
			n = 3L
	)

	factor_layout <- list(
			lag_indices = list(
					Global = list(c(1L, 2L)),
					Real = list(3L)
			),
			block_orders = list(
					Global = 1L,
					Real = 1L
			)
	)

	reference <- update_factor_blocks(
			smooth,
			list(factor = factor_layout)
	)

	actual <- update_factor_blocks_cpp(
			smooth$F_smooth,
			smooth$P_smooth,
			smooth$PPm_smooth,
			smooth$F_smooth_0,
			smooth$P_smooth_0,
			factor_layout
	)

	expect_equal(
			actual$coefficients,
			reference$coefficients,
			tolerance = 1e-10
	)

	expect_equal(
			actual$covariances,
			reference$covariances,
			tolerance = 1e-10
	)

	expect_identical(
			dim(actual$coefficients$Global),
			c(2L, 2L, 1L)
	)

	expect_identical(
			dim(actual$coefficients$Real),
			c(1L, 1L, 1L)
	)
})

test_that("C++ factor updates validate smoother dimensions", {
	smooth <- make_factor_smoother_fixture()
	factor_layout <- make_factor_layout_fixture()

	expect_error(
			update_factor_blocks_cpp(
					smooth$F_smooth[-1, , drop = FALSE],
					smooth$P_smooth,
					smooth$PPm_smooth,
					smooth$F_smooth_0,
					smooth$P_smooth_0,
					factor_layout
			),
			"Incompatible smoother dimensions"
	)

	expect_error(
			update_factor_blocks_cpp(
					smooth$F_smooth,
					smooth$P_smooth,
					smooth$PPm_smooth,
					smooth$F_smooth_0[-1],
					smooth$P_smooth_0,
					factor_layout
			),
			"Incompatible smoother dimensions"
	)
})

test_that("C++ factor updates reject invalid indices", {
	smooth <- make_factor_smoother_fixture()
	factor_layout <- make_factor_layout_fixture()

	factor_layout$lag_indices$Global[[1]][1] <- 100L

	expect_error(
			update_factor_blocks_cpp(
					smooth$F_smooth,
					smooth$P_smooth,
					smooth$PPm_smooth,
					smooth$F_smooth_0,
					smooth$P_smooth_0,
					factor_layout
			),
			"Invalid current-state index"
	)
})

