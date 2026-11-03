test_that("monthly loading maps select contemporaneous factors", {
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

	maps <- build_loading_maps(spec, layout)

	for (series in layout$monthly_vars) {
			map <- maps[[series]]

			allowed <- colnames(spec$loading_mask)[
					spec$loading_mask[series, ]
			]

			expect_equal(
					rownames(map$indices),
					allowed
			)

			expect_equal(
					ncol(map$indices),
					1L
			)

			expect_equal(
					as.integer(map$indices[, 1L]),
					unname(factor_layout$current_indices[allowed])
			)

			expect_equal(map$weights, 1)
			expect_equal(
					map$idio_index,
					unname(layout$monthly_indices[[series]])
			)
	}
})

test_that("quarterly loading maps use Mariano-Murasawa weights", {
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

	maps <- build_loading_maps(spec, layout)

	map <- maps$GDP

	expect_equal(
			map$weights,
			c(1, 2, 3, 2, 1)
	)

	expect_equal(ncol(map$indices), 5L)

	expect_equal(
			map$idio_index,
			unname(layout$quarterly_current[["GDP"]])
	)

	expect_equal(
			rownames(map$indices),
			colnames(spec$loading_mask)[
					spec$loading_mask["GDP", ]
			]
	)
})

test_that("loading maps preserve block-major factor ordering", {
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

	maps <- build_loading_maps(spec, layout)

	for (series in names(maps)) {
			map <- maps[[series]]

			for (factor in rownames(map$indices)) {
					block <- names(spec$blocks)[
							vapply(
									spec$blocks,
									function(x) factor %in% x,
									logical(1)
							)
					]

					position <- match(
							factor,
							spec$blocks[[block]]
					)

					for (lag in seq_len(ncol(map$indices))) {
						block_lags <- factor_layout$lag_indices[[block]]
						lag_indices <- block_lags[[lag]]
						expected <- lag_indices[[position]]

							expect_equal(
									unname(map$indices[factor, lag]),
									unname(expected)
							)
					}
			}
	}
})

