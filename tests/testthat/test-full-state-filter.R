test_that("SKFS accepts complete block-restricted state space", {
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

	result <- SKFS(
			X = X,
			A = matrices$A,
			C = matrices$C,
			Q = matrices$Q,
			R = matrices$R,
			F_0 = initial$F_0,
			P_0 = initial$P_0,
			loglik = FALSE
	)

	expect_type(result, "list")
	expect_true(length(result) > 0L)
})

test_that("SKFS handles mixed-frequency missing observations", {
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

	expect_equal(
			sum(!is.na(X[, "GDP"])),
			40L
	)

	expect_true(
			anyNA(X[, "PMI"])
	)

	expect_no_error(
			SKFS(
					X = X,
					A = matrices$A,
					C = matrices$C,
					Q = matrices$Q,
					R = matrices$R,
					F_0 = initial$F_0,
					P_0 = initial$P_0,
					loglik = FALSE
			)
	)
})

test_that("SKFS handles singular state innovation covariance", {
	ss <- make_test_full_state()

	matrices <- ss$matrices

	eigenvalues <- eigen(
			matrices$Q,
			symmetric = TRUE,
			only.values = TRUE
	)$values

	expect_true(
			any(abs(eigenvalues) < 1e-12)
	)

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

	expect_no_error(
			SKFS(
					X = X,
					A = matrices$A,
					C = matrices$C,
					Q = matrices$Q,
					R = matrices$R,
					F_0 = initial$F_0,
					P_0 = initial$P_0,
					loglik = FALSE
			)
	)
})

