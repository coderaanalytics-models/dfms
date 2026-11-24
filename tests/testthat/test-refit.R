test_that("update.dfm reuses fitted parameters for block models", {
	set.seed(123)

	X_old <- matrix(
		rnorm(400),
		nrow = 100,
		ncol = 4
	)

	colnames(X_old) <- c("x1", "x2", "x3", "x4")

	# Series -> factor mapping
	factors <- list(
		x1 = c("Global", "Real"),
		x2 = c("Global", "Real"),
		x3 = "Global",
		x4 = "Global"
	)

	X_new <- X_old
	X_new[100, 1] <- X_new[100, 1] + 1
	X_new[99,  2] <- X_new[99,  2] - 0.5

	fit <- suppressWarnings(
		DFM(
			X_old,
			factors = factors,
			factor_multiplicities = c(
				Global = 1L,
				Real = 1L
			),
			factor_orders = c(
				Global = 1L,
				Real = 1L
			),
			idio.ar1 = TRUE,
			em.method = "BM",
			min.iter = 1L,
			max.iter = 2L,
			save.full.state = TRUE
		)
	)

	updated <- update(
		fit,
		newdata = X_new,
		refit = FALSE
	)

	expect_s3_class(updated, "dfm")

	# Parameter matrices must be unchanged
	expect_equal(updated$A, fit$A)
	expect_equal(updated$C, fit$C)
	expect_equal(updated$Q, fit$Q)
	expect_equal(updated$R, fit$R)

	# Block-model metadata must be unchanged
	expect_equal(updated$factor.spec, fit$factor.spec)
	expect_equal(updated$state.layout, fit$state.layout)

	# Original scaling stats must be reused
	expect_equal(
		attr(updated$X_imp, "stats"),
		attr(fit$X_imp, "stats")
	)

	# Full state-space system matrices must remain unchanged
	expect_false(is.null(updated$ss_full))
	expect_equal(updated$ss_full$A, fit$ss_full$A)
	expect_equal(updated$ss_full$C, fit$ss_full$C)
	expect_equal(updated$ss_full$Q, fit$ss_full$Q)
	expect_equal(updated$ss_full$R, fit$ss_full$R)
	expect_equal(updated$ss_full$F_0, fit$ss_full$F_0)
	expect_equal(updated$ss_full$P_0, fit$ss_full$P_0)

	# Smoothed state means should update when observed values change
	expect_false(
		isTRUE(all.equal(
			updated$F_qml,
			fit$F_qml
		))
	)

	expect_false(
		isTRUE(all.equal(
			updated$ss_full$F_smooth,
			fit$ss_full$F_smooth
		))
	)

	# With fixed parameters and the same missingness pattern,
	# smoothed state covariances should be unchanged
	expect_equal(
		updated$P_qml,
		fit$P_qml
	)

	expect_equal(
		updated$ss_full$P_smooth,
		fit$ss_full$P_smooth
	)

	# Idiosyncratic state estimates should also update
	expect_false(is.null(updated$e))
	expect_equal(colnames(updated$e), colnames(X_old))

	expect_false(
		isTRUE(all.equal(
			updated$e,
			fit$e
		))
	)

	# AR(1) parameters are fitted parameters and must remain fixed
	expect_equal(updated$rho, fit$rho)

	# Bookkeeping
	expect_false(updated$refit)
	expect_equal(updated$em.method, fit$em.method)
	expect_equal(updated$tol, fit$tol)
	expect_equal(updated$converged, fit$converged)
})


test_that("update.dfm refits model parameters when refit is TRUE", {
	set.seed(123)

	X_old <- matrix(
		rnorm(400),
		nrow = 100,
		ncol = 4
	)

	colnames(X_old) <- paste0("x", 1:4)

	fit <- DFM(
		X_old,
		r = 2,
		p = 1,
		em.method = "none"
	)

	# Make the new sample materially different so re-estimation
	# should produce different parameter estimates.
	X_new <- X_old
	X_new[, 1] <- X_new[, 1] * 2
	X_new[, 2] <- X_new[, 2] + 0.5 * X_new[, 1]

	updated <- update(
		fit,
		newdata = X_new,
		refit = TRUE
	)

	expect_s3_class(updated, "dfm")

	# New fit should use scaling estimated from the new data.
	expect_false(
		isTRUE(all.equal(
			attr(updated$X_imp, "stats"),
			attr(fit$X_imp, "stats")
		))
	)

	# At least some fitted model parameters should change.
	expect_false(
		isTRUE(all.equal(
			updated$C,
			fit$C
		))
	)

	# Original model specification should be retained.
	expect_equal(updated$em.method, fit$em.method)
	expect_equal(nrow(updated$A), nrow(fit$A))
	expect_equal(ncol(updated$A), ncol(fit$A))
})


test_that("update.dfm can override fitting arguments when refitting", {
	set.seed(123)

	X <- matrix(
		rnorm(400),
		nrow = 100,
		ncol = 4
	)

	colnames(X) <- paste0("x", 1:4)

	fit <- DFM(
		X,
		r = 2,
		p = 1,
		em.method = "none"
	)

	updated <- update(
		fit,
		newdata = X,
		refit = TRUE,
		p = 2
	)

	expect_equal(
		ncol(updated$A),
		2L * nrow(updated$A)
	)
})


test_that("update.dfm errors on mismatched number of variables", {
	set.seed(123)

	X <- matrix(
		rnorm(300),
		nrow = 100,
		ncol = 3
	)
	colnames(X) <- c("x1", "x2", "x3")

	fit <- DFM(
		X,
		r = 1,
		p = 1,
		em.method = "none"
	)

	X_new <- cbind(
		X,
		x4 = rnorm(nrow(X))
	)

	expect_error(
		update(
			fit,
			newdata = X_new,
			refit = FALSE
		),
		"different number of variables"
	)
})


test_that("update.dfm errors on mismatched variable names or ordering", {
	set.seed(123)

	X <- matrix(
		rnorm(300),
		nrow = 100,
		ncol = 3
	)
	colnames(X) <- c("x1", "x2", "x3")

	fit <- DFM(
		X,
		r = 1,
		p = 1,
		em.method = "none"
	)

	X_new <- X[, c("x2", "x1", "x3")]

	expect_error(
		update(
			fit,
			newdata = X_new,
			refit = FALSE
		),
		"different variable names or ordering"
	)
})


test_that("update.dfm changes state covariance when missingness changes", {
	set.seed(123)

	X_old <- matrix(
		rnorm(400),
		nrow = 100,
		ncol = 4
	)
	colnames(X_old) <- paste0("x", 1:4)

	# Same underlying data, but the old vintage is missing one observation.
	X_new <- X_old
	X_old[100, 1] <- NA

	fit <- suppressWarnings(
		DFM(
			X_old,
			r = 2,
			p = 1,
			em.method = "BM",
			min.iter = 1L,
			max.iter = 2L,
			save.full.state = TRUE
		)
	)

	updated <- update(
		fit,
		newdata = X_new,
		refit = FALSE
	)

	expect_false(
		isTRUE(all.equal(
			updated$P_qml,
			fit$P_qml
		))
	)

	expect_false(
		isTRUE(all.equal(
			updated$ss_full$P_smooth,
			fit$ss_full$P_smooth
		))
	)
})


test_that("update.dfm works without ss_full for classical models", {
	set.seed(123)

	X_old <- matrix(
		rnorm(300),
		nrow = 100,
		ncol = 3
	)
	colnames(X_old) <- c("x1", "x2", "x3")

	X_new <- X_old
	X_new[100, 1] <- X_new[100, 1] + 1

	fit <- DFM(
		X_old,
		r = 1,
		p = 1,
		em.method = "none",
		save.full.state = FALSE
	)

	expect_null(fit$ss_full)

	updated <- update(
		fit,
		newdata = X_new,
		refit = FALSE
	)

	expect_s3_class(updated, "dfm")
	expect_null(updated$ss_full)

	expect_equal(updated$A, fit$A)
	expect_equal(updated$C, fit$C)
	expect_equal(updated$Q, fit$Q)
	expect_equal(updated$R, fit$R)

	expect_false(
		isTRUE(all.equal(
			updated$F_2s,
			fit$F_2s
		))
	)
})


test_that("update.dfm requires ss_full for block models", {
	set.seed(123)

	X <- matrix(
		rnorm(400),
		nrow = 100,
		ncol = 4
	)
	colnames(X) <- c("x1", "x2", "x3", "x4")

	factors <- list(
		x1 = c("Global", "Real"),
		x2 = c("Global", "Real"),
		x3 = "Global",
		x4 = "Global"
	)

	fit <- suppressWarnings(
		DFM(
			X,
			factors = factors,
			factor_multiplicities = c(
				Global = 1L,
				Real = 1L
			),
			factor_orders = c(
				Global = 1L,
				Real = 1L
			),
			idio.ar1 = TRUE,
			em.method = "BM",
			min.iter = 1L,
			max.iter = 2L,
			save.full.state = FALSE
		)
	)

	expect_null(fit$ss_full)

	expect_error(
		update(
			fit,
			newdata = X,
			refit = FALSE
		),
		"Full state matrices are required"
	)
})


test_that("update.dfm preserves block specification when refitting", {
	set.seed(123)

	X_old <- matrix(
		rnorm(400),
		nrow = 100,
		ncol = 4
	)
	colnames(X_old) <- c("x1", "x2", "x3", "x4")

	factors <- list(
		x1 = c("Global", "Real"),
		x2 = c("Global", "Real"),
		x3 = "Global",
		x4 = "Global"
	)

	fit <- suppressWarnings(
		DFM(
			X_old,
			factors = factors,
			factor_multiplicities = c(
				Global = 1L,
				Real = 1L
			),
			factor_orders = c(
				Global = 1L,
				Real = 1L
			),
			idio.ar1 = TRUE,
			em.method = "BM",
			min.iter = 1L,
			max.iter = 2L
		)
	)

	X_new <- X_old
	X_new[, 1] <- X_new[, 1] * 1.5

	updated <- suppressWarnings(
		update(
			fit,
			newdata = X_new,
			refit = TRUE
		)
	)

	expect_s3_class(updated, "dfm")

	expect_equal(
		updated$factor.spec$blocks,
		fit$factor.spec$blocks
	)

	expect_equal(
		updated$factor.spec$block_orders,
		fit$factor.spec$block_orders
	)

	expect_equal(
		updated$factor.spec$factor_names,
		fit$factor.spec$factor_names
	)

	expect_equal(
		updated$state.layout$factor$current_indices,
		fit$state.layout$factor$current_indices
	)
})


test_that("update.dfm does not mutate the original object", {
	set.seed(123)

	X_old <- matrix(
		rnorm(300),
		nrow = 100,
		ncol = 3
	)
	colnames(X_old) <- c("x1", "x2", "x3")

	fit <- DFM(
		X_old,
		r = 1,
		p = 1,
		em.method = "none"
	)

	original <- unserialize(serialize(fit, NULL))

	X_new <- X_old
	X_new[100, 1] <- X_new[100, 1] + 1

	invisible(
		update(
			fit,
			newdata = X_new,
			refit = FALSE
		)
	)

	expect_equal(fit, original)
})
