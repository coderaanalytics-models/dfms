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