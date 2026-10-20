test_that("global factors load on all series", {
	spec <- parse_factor_spec(
			series = c("PMI", "IP", "GDP"),
			factors = "global",
			quarterly_vars = "GDP"
	)

	expect_true(all(spec$loading_mask))
	expect_equal(spec$n_factors, 1L)
})

test_that("sector restrictions produce structural zeros", {
	spec <- parse_factor_spec(
			series = c("PMI", "IP", "SALES", "GDP"),
			factors = list(
					PMI = c("global", "manufacturing"),
					IP = c("global", "manufacturing"),
					SALES = c("global", "consumption"),
					GDP = "global"
			),
			quarterly_vars = "GDP"
	)

	expect_true(
			spec$loading_mask["PMI", "manufacturing"]
	)

	expect_false(
			spec$loading_mask["SALES", "manufacturing"]
	)

	expect_false(
			spec$loading_mask["GDP", "consumption"]
	)
})

test_that("factor multiplicities expand correctly", {
	spec <- parse_factor_spec(
			series = c("PMI", "IP", "GDP"),
			factors = "global",
			factor_multiplicities = 2L,
			quarterly_vars = "GDP"
	)

	expect_identical(
			spec$factor_names,
			c("global.1", "global.2")
	)

	expect_equal(spec$n_factors, 2L)
})

test_that("factor blocks are disjoint", {
	spec <- parse_factor_spec(
			series = c("PMI", "IP", "SALES", "GDP"),
			factors = list(
					PMI = c("global", "manufacturing"),
					IP = c("global", "manufacturing"),
					SALES = c("global", "consumption"),
					GDP = "global"
			),
			factor_orders = list(
					global = list(
							factors = "global",
							order = 1L
					),
					sectors = list(
							factors = c("manufacturing", "consumption"),
							order = 2L
					)
			),
			quarterly_vars = "GDP"
	)

	expect_equal(spec$n_blocks, 2L)

	expect_setequal(
			unlist(spec$blocks, use.names = FALSE),
			spec$factor_names
	)

	expect_identical(
			spec$block_orders,
			c(global = 1L, sectors = 2L)
	)
})

test_that("loading map follows observed series ordering", {
	series <- c("PMI", "IP", "SALES", "GDP")

	spec <- parse_factor_spec(
			series = series,
			factors = list(
					GDP = "global",
					SALES = c("global", "consumption"),
					PMI = c("global", "manufacturing"),
					IP = c("global", "manufacturing")
			),
			quarterly_vars = "GDP"
	)

	expect_identical(
			rownames(spec$loading_mask),
			series
	)
})

test_that("invalid specifications are rejected", {
	expect_error(
			parse_factor_spec(
					series = c("PMI", "PMI"),
					factors = "global"
			),
			"unique"
	)

	expect_error(
			parse_factor_spec(
					series = c("PMI", "IP"),
					factors = "global",
					factor_multiplicities = 3L
			),
			"monthly"
	)
})

test_that("named factor orders are parsed correctly", {
	spec <- parse_factor_spec(
			series = c("PMI", "IP", "SALES", "GDP"),
			factors = list(
					PMI = c("global", "manufacturing"),
					IP = c("global", "manufacturing"),
					SALES = c("global", "consumption"),
					GDP = "global"
			),
			factor_orders = c(
					global = 1L,
					manufacturing = 2L,
					consumption = 1L
			),
			quarterly_vars = "GDP"
	)

	expect_identical(
			unname(spec$block_orders["manufacturing"]),
			2L
	)

	expect_equal(spec$n_blocks, 3L)
})