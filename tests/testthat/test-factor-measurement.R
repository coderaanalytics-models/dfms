test_that("measurement matrix has correct dimensions", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	Z <- build_factor_measurement(
			spec,
			layout,
			make_block_loadings(),
			quarterly_vars = "GDP"
	)

	expect_identical(dim(Z), c(5L, 20L))
	expect_identical(rownames(Z), rownames(spec$loading_mask))
})

test_that("monthly variables load on current factors only", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	Z <- build_factor_measurement(
			spec,
			layout,
			make_block_loadings(),
			quarterly_vars = "GDP"
	)

	current <- unname(layout$current_indices)

	lagged <- setdiff(
			seq_len(layout$n_states),
			current
	)

	expect_equal(
			unname(Z["PMI", current]),
			c(0.8, 0.3, 0.0, 0.5)
	)

	expect_true(all(Z["PMI", lagged] == 0))
})

test_that("quarterly variables use Mariano-Murasawa weights", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	Z <- build_factor_measurement(
			spec,
			layout,
			make_block_loadings(),
			quarterly_vars = "GDP"
	)

	global_states <- layout$block_indices$global

	expect_equal(
			unname(Z["GDP", global_states]),
			c(
					0.9, 0.4,
					1.8, 0.8,
					2.7, 1.2,
					1.8, 0.8,
					0.9, 0.4
			)
	)

	other_states <- setdiff(
			seq_len(layout$n_states),
			global_states
	)

	expect_true(all(Z["GDP", other_states] == 0))
})

test_that("restricted factor loadings remain zero", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	loadings <- make_block_loadings()

	# Deliberately violate the restrictions.
	loadings["PMI", "consumption"] <- 100
	loadings["GDP", "manufacturing"] <- 100

	Z <- build_factor_measurement(
			spec,
			layout,
			loadings,
			quarterly_vars = "GDP"
	)

	expect_true(
			all(
					Z[
							"PMI",
							layout$block_indices$consumption
					] == 0
			)
	)

	expect_true(
			all(
					Z[
							"GDP",
							layout$block_indices$manufacturing
					] == 0
			)
	)
})

test_that("quarterly aggregation requires quarterly state layout", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = FALSE
	)

	expect_error(
			build_factor_measurement(
					spec,
					layout,
					make_block_loadings(),
					quarterly_vars = "GDP"
			),
			"quarterly = TRUE"
	)
})

test_that("monthly-only measurement matrices work", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = FALSE
	)

	Z <- build_factor_measurement(
			spec,
			layout,
			make_block_loadings()
	)

	expect_identical(dim(Z), c(5L, 5L))

	current <- unname(layout$current_indices)

	expect_equal(
			unname(Z["PMI", current]),
			c(0.8, 0.3, 0.0, 0.5)
	)
})
