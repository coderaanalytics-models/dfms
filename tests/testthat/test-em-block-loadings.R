test_that("compiled monthly loadings match R reference", {
	ss <- make_test_full_state()

	initial <- initialise_full_state(
			ss$matrices$A,
			ss$matrices$Q
	)

	sim <- simulate_full_state(
			ss$matrices,
			initial,
			n_periods = 120L
	)

	X <- make_mixed_frequency_observations(
			sim$observations,
			quarterly_vars = "GDP"
	)

	smooth <- SKFS(
			X = X,
			A = ss$matrices$A,
			C = ss$matrices$C,
			Q = ss$matrices$Q,
			R = ss$matrices$R,
			F_0 = initial$F_0,
			P_0 = initial$P_0,
			loglik = FALSE
	)

	maps <- build_loading_maps(
			ss$spec,
			ss$layout
	)

	compiled <- update_loadings_cpp(
			X,
			smooth$F_smooth,
			smooth$P_smooth,
			maps
	)

	# Compare monthly series against the reference estimator.
	for (series in ss$layout$monthly_vars) {
			map <- maps[[series]]
		
			reference <- update_monthly_loadings(
					X = X,
					smooth = smooth,
					factor_indices = setNames(
							list(as.integer(map$indices[, 1L])),
							series
					),
					idio_indices = setNames(
							map$idio_index,
							series
					)
			)

			expect_equal(
					compiled[[series]],
					reference[[series]],
					tolerance = 1e-10
			)
	}
})

test_that("compiled quarterly loadings match R reference", {
	ss <- make_test_full_state()

	initial <- initialise_full_state(
			ss$matrices$A,
			ss$matrices$Q
	)

	sim <- simulate_full_state(
			ss$matrices,
			initial,
			n_periods = 120L
	)

	X <- make_mixed_frequency_observations(
			sim$observations,
			quarterly_vars = "GDP"
	)

	smooth <- SKFS(
			X = X,
			A = ss$matrices$A,
			C = ss$matrices$C,
			Q = ss$matrices$Q,
			R = ss$matrices$R,
			F_0 = initial$F_0,
			P_0 = initial$P_0,
			loglik = FALSE
	)

	maps <- build_loading_maps(
			ss$spec,
			ss$layout
	)

	compiled <- update_loadings_cpp(
			X,
			smooth$F_smooth,
			smooth$P_smooth,
			maps
	)

	map <- maps$GDP

	lag_indices <- lapply(
			seq_len(ncol(map$indices)),
			function(lag) as.integer(map$indices[, lag])
	)

	reference <- update_quarterly_loadings(
			X = X,
			smooth = smooth,
			factor_indices = list(GDP = lag_indices),
			idio_indices = c(GDP = map$idio_index)
	)

	expect_equal(
			compiled$GDP,
			reference$GDP,
			tolerance = 1e-10
	)
})

test_that("loading update rejects invalid factor-state indices", {
  # 5 observations, 1 observed series, 2 latent states
  X <- matrix(rnorm(5), nrow = 5, ncol = 1)
  F_smooth <- matrix(rnorm(10), nrow = 5, ncol = 2)
  P_smooth <- array(0, dim = c(2, 2, 5))

  loading_maps <- list(
    series_1 = list(
      indices = matrix(3L, nrow = 1, ncol = 1),
      weights = 1,
      idio_index = 2L
    )
  )

  expect_error(
    update_loadings_cpp(X, F_smooth, P_smooth, loading_maps),
    "Invalid factor-state index"
  )
  
  loading_maps[[1]]$indices[1, 1] <- 0L

	expect_error(
		update_loadings_cpp(X, F_smooth, P_smooth, loading_maps),
		"Invalid factor-state index"
	)

	loading_maps[[1]]$indices[1, 1] <- NA_integer_

	expect_error(
		update_loadings_cpp(X, F_smooth, P_smooth, loading_maps),
		"Invalid factor-state index"
	)
})