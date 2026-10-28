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