test_that("dfm_apply_state reuses fitted parameters", {
	X_old <- matrix(
		rnorm(200),
		nrow = 50,
		ncol = 4
	)

	colnames(X_old) <- paste0("x", 1:4)

	X_new <- X_old
	X_new[50, 1] <- X_new[50, 1] + 0.5

	fit <- DFM(
		X_old,
		r = 2,
		p = 1,
		em.method = "none"
	)

	res <- dfm_apply_state(fit, X_new)

	state <- dfm_state(fit)

	expect_equal(res$state$A, state$A)
	expect_equal(res$state$C, state$C)
	expect_equal(res$state$Q, state$Q)
	expect_equal(res$state$R, state$R)

	expect_equal(
		res$prep$stats,
		attr(fit$X_imp, "stats")
	)

	expect_equal(
		nrow(res$kfs$F_smooth),
		nrow(res$prep$X)
	)

	expect_false(
		isTRUE(all.equal(
			res$kfs$F_smooth,
			state$F_smooth
		))
	)
})


test_that("dfm_apply_state reuses full state for block models", {
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

	X_new <- X_old
	X_new[100, 1] <- X_new[100, 1] + 0.5

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

	res <- dfm_apply_state(fit, X_new)

	expect_equal(res$state$A, fit$ss_full$A)
	expect_equal(res$state$C, fit$ss_full$C)
	expect_equal(res$state$Q, fit$ss_full$Q)
	expect_equal(res$state$R, fit$ss_full$R)

	expect_equal(
		res$prep$stats,
		attr(fit$X_imp, "stats")
	)

	expect_false(
		isTRUE(all.equal(
			res$kfs$F_smooth,
			fit$ss_full$F_smooth
		))
	)
})
