params = make_block_parameters()

test_that("block matrices have correct dimensions", {
	spec <- make_block_spec()
	layout <- build_factor_state_layout(spec, quarterly = TRUE)

	ss <- build_factor_matrices(
			layout,
			params$coefficients,
			params$covariances
	)

	expect_identical(dim(ss$T), c(20L, 20L))
	expect_identical(dim(ss$Q), c(20L, 20L))
})

test_that("different factor blocks do not interact", {
	spec <- make_block_spec()
	layout <- build_factor_state_layout(spec, quarterly = TRUE)

	ss <- build_factor_matrices(
			layout,
			params$coefficients,
			params$covariances
	)

	blocks <- layout$block_indices

	for (i in seq_along(blocks)) {
			for (j in seq_along(blocks)) {
					if (i == j) next

					expect_true(
							all(ss$T[blocks[[i]], blocks[[j]]] == 0)
					)

					expect_true(
							all(ss$Q[blocks[[i]], blocks[[j]]] == 0)
					)
			}
	}
})

test_that("lagged factor states shift correctly", {
	spec <- make_block_spec()
	layout <- build_factor_state_layout(spec, quarterly = TRUE)

	ss <- build_factor_matrices(
			layout,
			params$coefficients,
			params$covariances
	)

	lag_idx <- layout$lag_indices$global

	for (j in 2:length(lag_idx)) {
			expect_equal(
					ss$T[
							lag_idx[[j]],
							lag_idx[[j - 1L]]
					],
					diag(2)
			)
	}
})

test_that("manufacturing follows the specified AR(2)", {
	spec <- make_block_spec()
	layout <- build_factor_state_layout(spec, quarterly = TRUE)

	ss <- build_factor_matrices(
			layout,
			params$coefficients,
			params$covariances
	)

	idx <- layout$lag_indices$manufacturing

	expect_equal(ss$T[idx[[1]], idx[[1]]], 0.6)
	expect_equal(ss$T[idx[[1]], idx[[2]]], -0.2)

	expect_true(
			all(ss$T[idx[[1]], unlist(idx[3:5])] == 0)
	)
})

test_that("innovations enter current factors only", {
	spec <- make_block_spec()
	layout <- build_factor_state_layout(spec, quarterly = TRUE)

	ss <- build_factor_matrices(
			layout,
			params$coefficients,
			params$covariances
	)

	lagged <- setdiff(
			seq_len(layout$n_states),
			unname(layout$current_indices)
	)

	expect_true(all(ss$Q[lagged, ] == 0))
	expect_true(all(ss$Q[, lagged] == 0))
})