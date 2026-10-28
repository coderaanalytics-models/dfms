test_that("full state layout has correct dimensions", {
	spec <- make_block_spec()

	factors <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	layout <- build_full_state_layout(
			spec,
			factors,
			quarterly_vars = "GDP"
	)

	expect_equal(layout$n_states, 29L)
	expect_equal(layout$n_factor_states, 20L)
	expect_equal(layout$n_monthly, 4L)
	expect_equal(layout$n_quarterly, 1L)
})

test_that("factor state indices remain unchanged", {
	spec <- make_block_spec()

	factors <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	layout <- build_full_state_layout(
			spec,
			factors,
			quarterly_vars = "GDP"
	)

	expect_equal(
			layout$factor_indices,
			seq_len(factors$n_states)
	)

	expect_identical(
			layout$factor,
			factors
	)
})

test_that("monthly idiosyncratic states are indexed correctly", {
	spec <- make_block_spec()

	factors <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	layout <- build_full_state_layout(
			spec,
			factors,
			quarterly_vars = "GDP"
	)

	expect_equal(
			unname(layout$monthly_indices),
			21:24
	)

	expect_identical(
			names(layout$monthly_indices),
			c("PMI", "IP", "SALES", "RETAIL")
	)
})

test_that("quarterly idiosyncratic states have five lags", {
	spec <- make_block_spec()

	factors <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	layout <- build_full_state_layout(
			spec,
			factors,
			quarterly_vars = "GDP"
	)

	expect_equal(
			unname(layout$quarterly_indices$GDP),
			25:29
	)

	expect_identical(
			names(layout$quarterly_indices$GDP),
			paste0("lag", 0:4)
	)

	expect_equal(
			unname(layout$quarterly_current),
			25L
	)

	expect_equal(
			layout$quarterly_lagged,
			26:29
	)
})

test_that("all states are assigned exactly once", {
	spec <- make_block_spec()

	factors <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	layout <- build_full_state_layout(
			spec,
			factors,
			quarterly_vars = "GDP"
	)

	indices <- c(
			layout$factor_indices,
			unname(layout$monthly_indices),
			unlist(
					layout$quarterly_indices,
					use.names = FALSE
			)
	)

	expect_equal(
			sort(indices),
			seq_len(layout$n_states)
	)

	expect_equal(
			anyDuplicated(indices),
			0L
	)
})

test_that("quarterly variables require quarterly factor layout", {
	spec <- make_block_spec()

	factors <- build_factor_state_layout(
			spec,
			quarterly = FALSE
	)

	expect_error(
			build_full_state_layout(
					spec,
					factors,
					quarterly_vars = "GDP"
			),
			"quarterly = TRUE"
	)
})

test_that("monthly-only state layouts work", {
	spec <- make_block_spec()

	factors <- build_factor_state_layout(
			spec,
			quarterly = FALSE
	)

	layout <- build_full_state_layout(
			spec,
			factors
	)

	expect_equal(layout$n_factor_states, 5L)
	expect_equal(layout$n_monthly, 5L)
	expect_equal(layout$n_quarterly, 0L)
	expect_equal(layout$n_states, 10L)

	expect_length(layout$quarterly_indices, 0L)
	expect_length(layout$quarterly_current, 0L)
	expect_length(layout$quarterly_lagged, 0L)
})