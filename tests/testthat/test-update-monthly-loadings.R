test_that("monthly loading update matches analytical EM solution", {
	X <- matrix(
			c(2, 4, 6),
			ncol = 1,
			dimnames = list(NULL, "PMI")
	)

	# State ordering: factor, idiosyncratic component.
	F <- matrix(
			c(
					1, 0.5,
					2, 1.0,
					3, 1.5
			),
			nrow = 3,
			byrow = TRUE
	)

	P <- array(0, dim = c(2, 2, 3))

	for (t in 1:3) {
			P[, , t] <- matrix(
					c(
							0.4, 0.1,
							0.1, 0.2
					),
					nrow = 2
			)
	}

	smooth <- list(
			F_smooth = F,
			P_smooth = P
	)

	result <- update_monthly_loadings(
			X = X,
			smooth = smooth,
			factor_indices = list(PMI = 1L),
			idio_indices = c(PMI = 2L)
	)

	S_ff <- sum(c(1, 2, 3)^2) + 3 * 0.4

	S_fy <- sum(c(1, 2, 3) * c(2, 4, 6))

	S_fu <- sum(
			c(1, 2, 3) * c(0.5, 1, 1.5)
	) + 3 * 0.1

	expected <- (S_fy - S_fu) / S_ff

	expect_equal(
			result$PMI,
			expected,
			tolerance = 1e-12
	)
})

test_that("monthly loading update respects factor restrictions", {
	X <- matrix(
			c(2, 3, 4),
			ncol = 1,
			dimnames = list(NULL, "PMI")
	)

	# States: Global factor, Real factor,
	# Financial factor, PMI idiosyncratic component.
	F <- matrix(
			c(
					1, 2, 10, 0.1,
					2, 1, 20, 0.2,
					3, 4, 30, 0.3
			),
			nrow = 3,
			byrow = TRUE
	)

	P <- array(0, dim = c(4, 4, 3))

	for (t in 1:3) {
			P[, , t] <- diag(0.1, 4)
	}

	smooth <- list(
			F_smooth = F,
			P_smooth = P
	)

	result <- update_monthly_loadings(
			X = X,
			smooth = smooth,
			factor_indices = list(PMI = c(1L, 2L)),
			idio_indices = c(PMI = 4L)
	)

	expect_length(result$PMI, 2L)

	# The Financial factor at index 3 is excluded.
	expect_true(all(is.finite(result$PMI)))
})

test_that("monthly loading update excludes missing observations", {
	X <- matrix(
			c(2, NA, 6),
			ncol = 1,
			dimnames = list(NULL, "PMI")
	)

	F <- matrix(
			c(
					1, 0.5,
					2, 1.0,
					3, 1.5
			),
			nrow = 3,
			byrow = TRUE
	)

	P <- array(0, dim = c(2, 2, 3))

	smooth <- list(
			F_smooth = F,
			P_smooth = P
	)

	result <- update_monthly_loadings(
			X = X,
			smooth = smooth,
			factor_indices = list(PMI = 1L),
			idio_indices = c(PMI = 2L)
	)

	expected <- (
			1 * (2 - 0.5) +
			3 * (6 - 1.5)
	) / (1^2 + 3^2)

	expect_equal(
			result$PMI,
			expected,
			tolerance = 1e-12
	)
})

