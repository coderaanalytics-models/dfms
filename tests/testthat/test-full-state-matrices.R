test_that("complete state-space matrices have correct dimensions", {
	model <- make_test_full_state()

	ss <- model$matrices

	expect_identical(dim(ss$A), c(29L, 29L))
	expect_identical(dim(ss$Q), c(29L, 29L))
	expect_identical(dim(ss$C), c(5L, 29L))
	expect_identical(dim(ss$R), c(5L, 5L))
})

test_that("factor dynamics are preserved in complete state", {
	model <- make_test_full_state()

	ss <- model$matrices
	fi <- model$layout$factor_indices

	params <- make_block_parameters()

	factor_ss <- build_factor_matrices(
			model$layout$factor,
			params$coefficients,
			params$covariances
	)

	expect_equal(
			unname(ss$A[fi, fi]),
			unname(factor_ss$T)
	)

	expect_equal(
			unname(ss$Q[fi, fi]),
			unname(factor_ss$Q)
	)
})

test_that("monthly idiosyncratic AR1 dynamics are correct", {
	model <- make_test_full_state()

	ss <- model$matrices
	idx <- model$layout$monthly_indices

	expect_equal(ss$A[idx["PMI"], idx["PMI"]], 0.4)
	expect_equal(ss$A[idx["IP"], idx["IP"]], 0.6)
	expect_equal(ss$A[idx["SALES"], idx["SALES"]], 0.3)
	expect_equal(ss$A[idx["RETAIL"], idx["RETAIL"]], 0.5)

	expect_equal(ss$Q[idx["PMI"], idx["PMI"]], 0.2)
	expect_equal(ss$Q[idx["IP"], idx["IP"]], 0.3)

	expect_equal(unname(ss$C["PMI", idx["PMI"]]), 1)
	expect_equal(unname(ss$C["IP", idx["IP"]]), 1)
})

test_that("quarterly idiosyncratic companion matrix is correct", {
	model <- make_test_full_state()

	ss <- model$matrices
	idx <- unname(model$layout$quarterly_indices$GDP)

	expected <- matrix(0, 5, 5)

	expected[1, 1] <- 0.7

	for (j in 2:5) {
			expected[j, j - 1] <- 1
	}

	expect_equal(
			unname(ss$A[idx, idx]),
			expected
	)
})

test_that("quarterly idiosyncratic measurement uses MM weights", {
	model <- make_test_full_state()

	ss <- model$matrices
	idx <- unname(model$layout$quarterly_indices$GDP)

	expect_equal(
			unname(ss$C["GDP", idx]),
			c(1, 2, 3, 2, 1)
	)
})

test_that("quarterly lag states have zero innovation variance", {
	model <- make_test_full_state()

	ss <- model$matrices

	lagged <- model$layout$quarterly_lagged

	expect_true(all(ss$Q[lagged, ] == 0))
	expect_true(all(ss$Q[, lagged] == 0))
})

test_that("idiosyncratic innovations are independent", {
	model <- make_test_full_state()

	ss <- model$matrices

	fi <- model$layout$factor_indices
	ii <- model$layout$idio_indices

	expect_true(all(ss$Q[fi, ii] == 0))
	expect_true(all(ss$Q[ii, fi] == 0))

	expect_true(all(ss$A[fi, ii] == 0))
	expect_true(all(ss$A[ii, fi] == 0))
})

test_that("observation noise is small and diagonal", {
	model <- make_test_full_state()

	ss <- model$matrices

	expect_equal(
			ss$R,
			diag(1e-6, 5)
	)
})

test_that("nonstationary idiosyncratic AR1 is rejected", {
	spec <- make_block_spec()

	factor_layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	layout <- build_full_state_layout(
			spec,
			factor_layout,
			quarterly_vars = "GDP"
	)

	params <- make_block_parameters()

	factor_ss <- build_factor_matrices(
			factor_layout,
			params$coefficients,
			params$covariances
	)

	factor_C <- build_factor_measurement(
			spec,
			factor_layout,
			make_block_loadings(),
			quarterly_vars = "GDP"
	)

	idio <- make_full_state_parameters()

	idio$quarterly_rho["GDP"] <- 1.1

	expect_error(
			build_full_state_matrices(
					spec,
					layout,
					factor_ss,
					factor_C,
					monthly_rho = idio$monthly_rho,
					monthly_variances = idio$monthly_variances,
					quarterly_rho = idio$quarterly_rho,
					quarterly_variances = idio$quarterly_variances
			),
			"stationary"
	)
})
