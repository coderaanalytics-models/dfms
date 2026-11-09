# Block DFM fixtures

make_block_spec <- function() {
	parse_factor_spec(
			series = c("PMI", "IP", "SALES", "RETAIL", "GDP"),
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

make_block_parameters <- function() {
	list(
			coefficients = list(
					global = array(
							c(0.7, 0.1, 0.2, 0.6),
							dim = c(2, 2, 1)
					),
					consumption = array(
							0.5,
							dim = c(1, 1, 1)
					),
					manufacturing = array(
							c(0.6, -0.2),
							dim = c(1, 1, 2)
					)
			),
			covariances = list(
					global = matrix(
							c(1.0, 0.2, 0.2, 0.8),
							nrow = 2
					),
					consumption = matrix(0.5),
					manufacturing = matrix(0.7)
			)
	)
}

make_block_loadings <- function() {
	matrix(
			c(
					# global.1, global.2, consumption, manufacturing
					0.8, 0.3, 0.0, 0.5,  # PMI
					0.7, 0.2, 0.0, 0.6,  # IP
					0.6, 0.1, 0.4, 0.0,  # SALES
					0.5, 0.2, 0.3, 0.0,  # RETAIL
					0.9, 0.4, 0.0, 0.0   # GDP
			),
			nrow = 5,
			byrow = TRUE,
			dimnames = list(
					c("PMI", "IP", "SALES", "RETAIL", "GDP"),
					c(
							"global.1",
							"global.2",
							"consumption",
							"manufacturing"
					)
			)
	)
}

make_full_state_parameters <- function() {
	list(
			monthly_rho = c(
					PMI = 0.4,
					IP = 0.6,
					SALES = 0.3,
					RETAIL = 0.5
			),

			monthly_variances = c(
					PMI = 0.2,
					IP = 0.3,
					SALES = 0.4,
					RETAIL = 0.5
			),

			quarterly_rho = c(
					GDP = 0.7
			),

			quarterly_variances = c(
					GDP = 0.6
			)
	)
}

make_test_full_state <- function() {
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

	factor_matrices <- build_factor_matrices(
			factor_layout,
			params$coefficients,
			params$covariances
	)

	factor_measurement <- build_factor_measurement(
			spec,
			factor_layout,
			make_block_loadings(),
			quarterly_vars = "GDP"
	)

	idio <- make_full_state_parameters()

	matrices <- build_full_state_matrices(
			spec,
			layout,
			factor_matrices,
			factor_measurement,
			monthly_rho = idio$monthly_rho,
			monthly_variances = idio$monthly_variances,
			quarterly_rho = idio$quarterly_rho,
			quarterly_variances = idio$quarterly_variances
	)

	list(
			spec = spec,
			layout = layout,
			matrices = matrices
	)
}

simulate_full_state <- function(
	matrices,
	initial,
	n_periods = 120L,
	seed = 123L
) {
	set.seed(seed)

	A <- matrices$A
	C <- matrices$C
	Q <- matrices$Q
	R <- matrices$R

	n_states <- nrow(A)
	n_series <- nrow(C)

	draw_normal <- function(covariance) {
			eig <- eigen(
					covariance,
					symmetric = TRUE
			)

			values <- pmax(eig$values, 0)

			as.vector(
					eig$vectors %*%
							(sqrt(values) * rnorm(length(values)))
			)
	}

	states <- matrix(
			0,
			nrow = n_periods,
			ncol = n_states
	)

	observations <- matrix(
			NA_real_,
			nrow = n_periods,
			ncol = n_series
	)

	alpha <- initial$F_0 + draw_normal(initial$P_0)

	for (t in seq_len(n_periods)) {
			alpha <- as.vector(
					A %*% alpha + draw_normal(Q)
			)

			y <- as.vector(
					C %*% alpha + draw_normal(R)
			)

			states[t, ] <- alpha
			observations[t, ] <- y
	}

	rownames(observations) <- NULL
	colnames(observations) <- rownames(C)

	list(
			states = states,
			observations = observations
	)
}

make_mixed_frequency_observations <- function(
	observations,
	quarterly_vars = "GDP",
	missing_rate = 0.05,
	seed = 456L
) {
	set.seed(seed)

	X <- observations

	n <- nrow(X)

	# Assume the first observation corresponds to the
	# first month of a quarter.
	quarter_end <- seq.int(3L, n, by = 3L)

	for (variable in quarterly_vars) {
			X[-quarter_end, variable] <- NA_real_
	}

	monthly_vars <- setdiff(
			colnames(X),
			quarterly_vars
	)

	if (missing_rate > 0 && length(monthly_vars) > 0L) {
			missing <- matrix(
					runif(n * length(monthly_vars)) < missing_rate,
					nrow = n,
					ncol = length(monthly_vars)
			)

			monthly_data <- X[, monthly_vars, drop = FALSE]
			monthly_data[missing] <- NA_real_

			X[, monthly_vars] <- monthly_data
	}

	X
}

# Specifications, state construction, and loading maps

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

test_that("VAR orders are distinct from retained state lags", {
	spec <- make_block_spec()

	layout <- build_factor_state_layout(
			spec,
			quarterly = TRUE
	)

	expect_equal(
			layout$block_orders,
			c(
					global = 1L,
					consumption = 1L,
					manufacturing = 2L
			)
	)

	expect_equal(
			layout$block_lags,
			c(
					global = 5L,
					consumption = 5L,
					manufacturing = 5L
			)
	)
})

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



test_that("initialisation has correct dimensions", {
	ss <- make_test_full_state()


	init <- initialise_full_state(
			ss$matrices$A,
			ss$matrices$Q
	)

	n <- nrow(ss$matrices$A)

	expect_length(init$F_0, n)
	expect_equal(dim(init$P_0), c(n, n))

	expect_equal(init$F_0, rep(0, n))
})

test_that("initial covariance solves the Lyapunov equation", {
	ss <- make_test_full_state()

	init <- initialise_full_state(
			ss$matrices$A,
			ss$matrices$Q
	)

	P <- init$P_0

	reconstructed <- ss$matrices$A %*% P %*% t(ss$matrices$A) + ss$matrices$Q

	expect_equal(
			P,
			reconstructed,
			tolerance = 1e-8
	)
})

test_that("initial covariance is symmetric", {
	ss <- make_test_full_state()

	init <- initialise_full_state(
			ss$matrices$A,
			ss$matrices$Q
	)

	expect_equal(
			init$P_0,
			t(init$P_0),
			tolerance = 1e-10
	)
})

test_that("initial covariance is positive semidefinite", {
	ss <- make_test_full_state()

	init <- initialise_full_state(
			ss$matrices$A,
			ss$matrices$Q
	)

	eigenvalues <- eigen(
			init$P_0,
			symmetric = TRUE,
			only.values = TRUE
	)$values

	expect_true(all(eigenvalues >= -1e-10))
})

test_that("nonstationary transition matrices are rejected", {
	A <- matrix(1.05, 1, 1)
	Q <- matrix(1, 1, 1)

	expect_error(
			initialise_full_state(A, Q),
			"not sufficiently stationary"
	)
})

test_that("scalar AR1 initialisation matches analytical variance", {
	rho <- 0.7
	sigma2 <- 0.4

	A <- matrix(rho, 1, 1)
	Q <- matrix(sigma2, 1, 1)

	init <- initialise_full_state(A, Q)

	expected_variance <- sigma2 / (1 - rho^2)

	expect_equal(
			unname(init$P_0[1, 1]),
			expected_variance,
			tolerance = 1e-10
	)
})

test_that("quarterly idiosyncratic covariance has AR1 structure", {
	ss <- make_test_full_state()

	init <- initialise_full_state(
			ss$matrices$A,
			ss$matrices$Q
	)

	idx <- ss$layout$quarterly_indices[["GDP"]]

	rho <- 0.7
	sigma2 <- 0.6

	expected <- outer(
			0:4,
			0:4,
			function(i, j) {
					sigma2 / (1 - rho^2) * rho^abs(i - j)
			}
	)

	actual <- init$P_0[idx, idx]

	expect_equal(
			unname(actual),
			expected,
			tolerance = 1e-8
	)
})

# Compiled M-step updates

make_block_moments <- function() {
  # Deterministic, asymmetric lag-one covariances distinguish current and
  # previous states. The initial moments differ from the observed periods.
  n <- 6L
  periods <- 8L
  list(
    F_smooth = outer(seq_len(periods), seq_len(n), function(t, i) {
      sin(t * i / 3) + cos(t + i / 4)
    }),
    P_smooth = array(rep(diag(n) + 0.04, periods), c(n, n, periods)),
    PPm_smooth = array(
      rep(outer(seq_len(n), seq_len(n), function(i, j) (2 * i - j) / 200),
          periods),
      c(n, n, periods)
    ),
    F_smooth_0 = seq_len(n) / 5,
    P_smooth_0 = diag(seq_len(n) / 3) + 0.02
  )
}

make_block_var_layout <- function() {
  list(
    lag_indices = list(pair = list(1:2, 3:4), single = list(5L, 6L)),
    block_orders = list(pair = 2L, single = 1L)
  )
}

test_that("factor VAR updates satisfy the expected normal equations", {
  moments <- make_block_moments()
  layout <- make_block_var_layout()
  actual <- do.call(update_factor_blocks_cpp, c(moments, list(layout)))
  expect_named(actual, c("coefficients", "covariances"))
  expect_named(actual$coefficients, names(layout$lag_indices))
  expect_named(actual$covariances, names(layout$lag_indices))

  periods <- nrow(moments$F_smooth)
  previous <- rbind(moments$F_smooth_0, moments$F_smooth[-periods, ])
  for (block in names(layout$lag_indices)) {
    lags <- layout$lag_indices[[block]]
    order <- layout$block_orders[[block]]
    current <- lags[[1L]]
    predictors <- unlist(lags[seq_len(order)], use.names = FALSE)
    y <- moments$F_smooth[, current, drop = FALSE]
    x <- previous[, predictors, drop = FALSE]

    xx <- crossprod(x) + moments$P_smooth_0[predictors, predictors]
    for (t in seq_len(periods - 1L)) {
      xx <- xx + moments$P_smooth[predictors, predictors, t]
    }
    yx <- crossprod(y, x)
    yy <- crossprod(y)
    for (t in seq_len(periods)) {
      yx <- yx + moments$PPm_smooth[current, predictors, t]
      yy <- yy + moments$P_smooth[current, current, t]
    }

    coefficients <- actual$coefficients[[block]]
    expect_identical(dim(coefficients), c(rep(length(current), 2L), order))
    phi <- matrix(coefficients, nrow = length(current))
    q <- actual$covariances[[block]]
    expect_equal(phi %*% xx, yx, tolerance = 1e-10)
    expect_equal(periods * q + phi %*% t(yx), yy, tolerance = 1e-10)
    expect_equal(q, t(q), tolerance = 1e-12)
    expect_true(min(eigen(q, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10)
  }
})

test_that("idiosyncratic AR(1) uses initial moments and only requested states", {
  moments <- make_block_moments()
  indices <- c(GDP = 5L, PMI = 2L)
  actual <- do.call(update_idiosyncratic_ar1_cpp, c(moments, list(indices)))
  expect_named(actual, c("coefficients", "variances"))
  expect_named(actual$coefficients, names(indices))
  expect_named(actual$variances, names(indices))

  for (series in names(indices)) {
    i <- indices[[series]]
    y <- moments$F_smooth[, i]
    periods <- length(y)
    x <- c(moments$F_smooth_0[i], y[-periods])
    v <- moments$P_smooth[i, i, ]
    xx <- sum(x^2) + moments$P_smooth_0[i, i] + sum(v[-periods])
    yx <- sum(y * x + moments$PPm_smooth[i, i, ])
    yy <- sum(y^2 + v)
    expect_equal(actual$coefficients[[series]], yx / xx, tolerance = 1e-12)
    expect_equal(actual$variances[[series]], (yy - yx^2 / xx) / periods,
                 tolerance = 1e-12)
  }

  empty <- do.call(update_idiosyncratic_ar1_cpp, c(moments, list(integer())))
  expect_identical(empty, list(coefficients = numeric(), variances = numeric()))
})

test_that("transition updates reject invalid dimensions and state indices", {
  moments <- make_block_moments()
  layout <- make_block_var_layout()
  for (field in names(moments)) {
    broken <- moments
    broken[[field]] <- switch(field,
      F_smooth = moments$F_smooth[-1L, , drop = FALSE],
      P_smooth = moments$P_smooth[, , -1L, drop = FALSE],
      PPm_smooth = moments$PPm_smooth[, , -1L, drop = FALSE],
      F_smooth_0 = moments$F_smooth_0[-1L],
      P_smooth_0 = moments$P_smooth_0[-1L, -1L]
    )
    expect_error(do.call(update_factor_blocks_cpp, c(broken, list(layout))),
                 "Incompatible smoother dimensions")
    expect_error(do.call(update_idiosyncratic_ar1_cpp, c(broken, list(5L))),
                 "Incompatible smoother dimensions")
  }

  short <- moments
  short$F_smooth <- moments$F_smooth[1L, , drop = FALSE]
  expect_error(do.call(update_factor_blocks_cpp, c(short, list(layout))),
               "At least two observations")
  expect_error(do.call(update_idiosyncratic_ar1_cpp, c(short, list(5L))),
               "At least two observations")

  for (indices in list(0L, 7L, NA_integer_, c(2L, 2L))) {
    expect_error(do.call(update_idiosyncratic_ar1_cpp, c(moments, list(indices))),
                 "Invalid idiosyncratic state indices")
  }

  broken <- layout
  broken$lag_indices$pair[[1L]][1L] <- 7L
  expect_error(do.call(update_factor_blocks_cpp, c(moments, list(broken))),
               "Invalid current-state index")
  broken <- layout
  broken$lag_indices$pair[[2L]][1L] <- NA_integer_
  expect_error(do.call(update_factor_blocks_cpp, c(moments, list(broken))),
               "Invalid predictor index")
  broken$lag_indices$pair[[2L]] <- 3L
  expect_error(do.call(update_factor_blocks_cpp, c(moments, list(broken))),
               "Inconsistent factor dimensions")
  broken$block_orders$pair <- 3L
  expect_error(do.call(update_factor_blocks_cpp, c(moments, list(broken))),
               "Invalid VAR order")
})

test_that("singular transition normal equations report estimation errors", {
  moments <- lapply(make_block_moments(), function(x) x * 0)
  expect_error(
    do.call(update_factor_blocks_cpp, c(moments, list(make_block_var_layout()))),
    "Failed to estimate VAR coefficients"
  )
  expect_error(do.call(update_idiosyncratic_ar1_cpp, c(moments, list(5L))),
               "Invalid idiosyncratic normal equations")
})

# Check a loading normal equation in measurement coordinates. Each row of
# H projects one allowed factor and d projects the idiosyncratic component.
expect_loading_equation <- function(lambda, y, moments, H, d) {
  gram <- matrix(0, nrow(H), nrow(H))
  rhs <- numeric(nrow(H))
  for (t in which(is.finite(y))) {
    f <- moments$F_smooth[t, ]
    E_ff <- tcrossprod(f) + moments$P_smooth[, , t]
    gram <- gram + H %*% E_ff %*% t(H)
    rhs <- rhs + drop(H %*% (f * y[t] - E_ff %*% d))
  }
  expect_equal(drop(gram %*% lambda), rhs, tolerance = 1e-10)
}

test_that("monthly loadings use observed rows and respect structural zeros", {
  moments <- make_block_moments()
  X <- cbind(PMI = c(3, NA, 5, NaN, 1, Inf, 8, 9),
             IP = c(NA, 1, 2, 3, 4, 5, 6, 7),
             Restricted = NA_real_)
  maps <- list(
    PMI = list(indices = matrix(c(1L, 3L), 2L), weights = 1, idio_index = 5L),
    IP = list(indices = matrix(2L, 1L), weights = 1, idio_index = 6L),
    Restricted = list(indices = matrix(integer(), 0L, 1L), weights = 1, idio_index = 4L)
  )
  actual <- update_loadings_cpp(X, moments$F_smooth, moments$P_smooth, maps)
  expect_named(actual, colnames(X))
  expect_identical(lengths(actual), c(PMI = 2L, IP = 1L, Restricted = 0L))

  basis <- diag(6)
  expect_loading_equation(actual$PMI, X[, "PMI"], moments,
                          basis[c(1, 3), ], basis[5, ])
  expect_loading_equation(actual$IP, X[, "IP"], moments,
                          basis[2, , drop = FALSE], basis[6, ])

  mask <- rbind(PMI = c(TRUE, FALSE, TRUE), IP = c(FALSE, TRUE, FALSE),
                Restricted = c(FALSE, FALSE, FALSE))
  colnames(mask) <- c("first", "second", "third")
  loadings <- assemble_loadings(actual, mask)
  expect_identical(dimnames(loadings), dimnames(mask))
  expect_true(all(loadings[!mask] == 0))
  expect_equal(unname(loadings["PMI", mask["PMI", ]]), actual$PMI)
  expect_equal(unname(loadings["IP", mask["IP", ]]), actual$IP)
  expect_error(assemble_loadings(list(1, 2, numeric()), mask),
               "Loading dimensions do not match")
})

test_that("loading updates reject malformed maps and singular normal equations", {
  moments <- make_block_moments()
  X <- matrix(seq_len(8))
  run <- function(maps, P = moments$P_smooth) {
    update_loadings_cpp(X, moments$F_smooth, P, maps)
  }
  map <- list(indices = matrix(1L, 1L), weights = 1, idio_index = 2L)
  expect_error(run(list()), "Number of loading maps")
  expect_error(run(list(map), moments$P_smooth[, , -1L, drop = FALSE]),
               "Incompatible smoother dimensions")
  expect_error(run(list(modifyList(map, list(weights = c(1, 2))))),
               "Incompatible loading-map dimensions")
  expect_error(run(list(modifyList(map, list(idio_index = 0L)))),
               "Invalid idiosyncratic state index")
  for (index in c(0L, 7L, NA_integer_)) {
    expect_error(run(list(modifyList(map, list(indices = matrix(index, 1L))))),
                 "Invalid factor-state index")
  }
  X[] <- NA_real_
  expect_error(run(list(map)), "Failed to solve loading normal equations")
})

# Mixed-frequency integration

make_block_em_case <- function() {
  model <- make_test_full_state()
  initial <- initialise_full_state(model$matrices$A, model$matrices$Q)
  state <- c(model$matrices, initial)
  simulated <- simulate_full_state(model$matrices, initial)
  list(
    X = make_mixed_frequency_observations(simulated$observations),
    state = state,
    spec = model$spec,
    layout = model$layout,
    maps = build_loading_maps(model$spec, model$layout)
  )
}

smooth_block_em_case <- function(model) {
  with(model$state, dfms::SKFS(model$X, A, C, Q, R, F_0, P_0, loglik = TRUE))
}

test_that("SKFS supplies finite EM moments for a ragged mixed-frequency panel", {
  model <- make_block_em_case()
  expect_equal(which(!is.na(model$X[, "GDP"])), seq.int(3L, 120L, by = 3L))
  expect_true(anyNA(model$X[, "PMI"]))
  expect_true(any(diag(model$state$Q) == 0))

  smooth <- smooth_block_em_case(model)
  n <- nrow(model$state$A)
  periods <- nrow(model$X)
  expect_identical(dim(smooth$F_smooth), c(periods, n))
  expect_identical(dim(smooth$P_smooth), c(n, n, periods))
  expect_identical(dim(smooth$PPm_smooth), c(n, n, periods))
  expect_length(smooth$F_smooth_0, n)
  expect_identical(dim(smooth$P_smooth_0), c(n, n))
  fields <- c("F_smooth", "P_smooth", "PPm_smooth", "F_smooth_0", "P_smooth_0", "loglik")
  expect_true(all(is.finite(unlist(smooth[fields]))))
})

# Known regressions. These expectations describe the intended behavior and
# currently fail against the supplied implementation.

test_that("quarterly loadings account for all weighted idiosyncratic lags", {
  model <- make_block_em_case()
  moments <- smooth_block_em_case(model)
  actual <- update_loadings_cpp(
    model$X, moments$F_smooth, moments$P_smooth, model$maps
  )

  # Derive the projections from the fitted state-space measurement row and
  # the factor layout, independently of the idiosyncratic map representation.
  mask <- model$spec$loading_mask["GDP", ]
  allowed <- names(mask)[mask]
  n <- nrow(model$state$A)
  H <- matrix(0, length(allowed), n)
  weights <- c(1, 2, 3, 2, 1)
  for (j in seq_along(allowed)) {
    block <- names(Filter(function(factors) allowed[j] %in% factors, model$spec$blocks))
    position <- match(allowed[j], model$spec$blocks[[block]])
    lags <- model$layout$factor$lag_indices[[block]][seq_along(weights)]
    H[j, vapply(lags, function(lag) lag[[position]], integer(1))] <- weights
  }
  d <- numeric(n)
  idio <- model$layout$quarterly_indices[["GDP"]]
  d[idio] <- model$state$C["GDP", idio]
  expect_loading_equation(actual$GDP, model$X[, "GDP"], moments, H, d)
})

test_that("EMstepBlock assembles a complete update from its supplied inputs", {
  model <- make_block_em_case()
  actual <- EMstepBlock(model$X, model$state, model$layout, model$maps)
  expect_named(actual, c("smooth", "factor_updates", "idio_updates",
                         "factor_matrices", "loading_updates", "loadings"))
  expect_true(is.finite(actual$smooth$loglik))
  expect_named(actual$idio_updates$coefficients, names(model$layout$idio_current))
  expect_named(actual$factor_updates$coefficients, names(model$spec$blocks))
  expect_identical(dimnames(actual$loadings), dimnames(model$spec$loading_mask))
  expect_true(all(is.finite(actual$loadings)))
  expect_true(all(actual$loadings[!model$spec$loading_mask] == 0))
})
