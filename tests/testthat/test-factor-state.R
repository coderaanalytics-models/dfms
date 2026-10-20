make_block_spec <- function() {
	parse_factor_spec(
			series = c(
					"PMI",
					"IP",
					"SALES",
					"RETAIL",
					"GDP"
			),
			factors = list(
					PMI = c("global", "manufacturing"),
					IP = c("global", "manufacturing"),
					SALES = c("global", "consumption"),
					RETAIL = c("global", "consumption"),
					GDP = "global"
			),
			factor_multiplicities = c(
					global = 2L,
					manufacturing = 1L,
					consumption = 1L
			),
			factor_orders = c(
					global = 1L,
					manufacturing = 2L,
					consumption = 1L
			),
			quarterly_vars = "GDP"
	)
}

test_that("factor state dimensions follow block VAR orders", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = FALSE
	)

	expect_equal(layout$n_states, 5L)

	expect_identical(
			layout$block_lags,
			c(
					global = 1L,
					consumption = 1L,
					manufacturing = 2L
			)
	)

	expect_identical(
			layout$current_indices,
			c(
					global.1 = 1L,
					global.2 = 2L,
					consumption = 3L,
					manufacturing = 4L
			)
	)
})

test_that("quarterly aggregation retains five factor lags", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	expect_equal(layout$n_states, 20L)

	expect_true(all(layout$block_lags >= 5L))

	expect_identical(
			layout$current_indices,
			c(
					global.1 = 1L,
					global.2 = 2L,
					consumption = 11L,
					manufacturing = 16L
			)
	)
})

test_that("lag indices cover the state without overlap", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	all_indices <- unlist(
			layout$lag_indices,
			use.names = FALSE
	)

	expect_setequal(
			all_indices,
			seq_len(layout$n_states)
	)

	expect_equal(
			length(unique(all_indices)),
			layout$n_states
	)
})

test_that("factor blocks occupy contiguous state regions", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	expect_identical(
			layout$block_indices$global,
			1:10
	)

	expect_identical(
			layout$block_indices$consumption,
			11:15
	)

	expect_identical(
			layout$block_indices$manufacturing,
			16:20
	)
})

test_that("higher VAR orders are preserved", {
	spec <- parse_factor_spec(
			series = paste0("X", 1:8),
			factors = "global",
			factor_orders = 6L
	)

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	expect_equal(layout$n_states, 6L)
	expect_identical(layout$block_lags, c(global = 6L))
})
