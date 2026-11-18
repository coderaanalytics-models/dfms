library(collapse)
library(magrittr)
library(xts)

#' @srrstats {G5.0} *Where applicable or practicable, tests should use standard data sets with known properties (for example, the [NIST Standard Reference Datasets](https://www.itl.nist.gov/div898/strd/), or data sets provided by other widely-used R packages).*
#' @srrstats {G5.1} *Data sets created within, and used to test, a package should be exported (or otherwise made generally available) so that users can confirm tests and run examples.*
#' @srrstats {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstats {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.3} *For functions which are expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values, the absence of any such values in return objects should be explicitly tested.*
#' @srrstats {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as [RStata](https://github.com/lbraglia/RStata)).*
#' -> I have translated the authors original matlab code into R and run tests with that code (see misc/ directory in the repo). These comparisons yielded that my implementation is equivalent to the original Matlab code.
#'  testing can be improved upon in general, an idea would be to run tests against the r-translations of those Matlab codes in the testthat framework, but I have not come round to do that yet and I would first like to gather some more substantive feedback on the software.


# BM14 Replication Data. Constructing the database:
BM14 = merge(BM14_M, BM14_Q)
BM14[, BM14_Models$log_trans] %<>% log()
BM14[, BM14_Models$freq == "M"] %<>% diff()
BM14[, BM14_Models$freq == "Q"] %<>% diff(3)

### Missing value removal
for(narm in c("LE", "all")) {
  for(m in c("median", "rnorm", "median.ma", "median.ma.spline"))
    expect_false(anyNA(tsnarmimp(BM14, na.rm.method = narm, na.impute = m)))
}


### Small Model ---------------------------------------

ics = ICr(BM14[, BM14_Models$small], max.r = 5)
expect_output(print(ics))

# Estimating the model with 2 factors and 3 lags
dfm_small = DFM(BM14[, BM14_Models$small], 2, 3, em.method = "BM", check.increased = TRUE)

expect_equal(dfm_small$A[1, 1], 1.21068973, tolerance = 1e-5)
expect_equal(dfm_small$A[2, 6], 0.391438990, tolerance = 1e-5)

# Rudimentry testing various methods
expect_output(print(dfm_small))
expect_visible(summary(dfm_small))
expect_true(is.matrix(fitted(dfm_small)))
expect_true(is.matrix(residuals(dfm_small)))
expect_true(ncol(as.data.frame(dfm_small)) == 4L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "wide.method")) == 5L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "wide.factor")) == 4L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "wide")) == 7L)
expect_true(ncol(as.data.frame(dfm_small, method = "qml", pivot = "wide")) == 3L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "t.wide")) == 7L)
expect_true(is.list(predict(dfm_small)))
expect_true(is.list(predict(dfm_small, standardized = FALSE)))
expect_true(is.list(predict(dfm_small, resFUN = function(x, h) predict(ar(na_rm(x)), n.ahead = h)$pred)))

expect_output(print(predict(dfm_small)))
expect_true(ncol(as.data.frame(predict(dfm_small))) == 4L)
expect_true(ncol(as.data.frame(predict(dfm_small), pivot = "wide")) == 4L)
expect_true(ncol(as.data.frame(predict(dfm_small), use = "data", pivot = "wide")) == 16L)
expect_true(ncol(as.data.frame(predict(dfm_small), use = "both", pivot = "wide")) == 18L)

# Other missing value options
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.rm.method = "all"))
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.impute = "median.ma"))
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.impute = "median"))
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.impute = "rnorm"))

### Medium-Sized Model ---------------------------------

# Estimating the model with 3 factors and 3 lags
dfm_medium = DFM(BM14[, BM14_Models$medium], 3, 3, em.method = "BM")

expect_equal(dfm_medium$A[1, 1], 0.74619087, tolerance = 1e-5)
expect_equal(dfm_medium$A[3, 9], 0.15380781, tolerance = 1e-5)

### Large Model ---------------------------------

# Estimating the model with 6 factors and 3 lags
dfm_large = DFM(BM14, 6, 3, em.method = "BM")

expect_equal(dfm_large$A[1, 1], 0.48915420, tolerance = 1e-5)
expect_equal(dfm_large$A[6, 18], 0.10027110, tolerance = 1e-5)

# Twostep estimates
dfm_large_2s = DFM(BM14, 6, 3, em.method = "none")
expect_equal(dfm_large_2s$A[1, 1], 0.6778009, tolerance = 1e-5)
expect_equal(dfm_large_2s$A[6, 18], 0.1314988, tolerance = 1e-5)

# Now testing DGR 2012 method
Xcols = colSums(is.na(BM14_M[70:350, ])) == 0
X_cc = BM14_M[70:350, Xcols]

expect_identical(ICr(diff(X_cc), max.r = 10)$r.star, c(IC1 = 7L, IC2 = 4L, IC3 = 7L))
mod = DFM(diff(X_cc), 4, 3)

expect_equal(mod$F_qml[10,], c(f1 = 0.4147068, f2 = 3.2364928, f3 = 1.6050372, f4 = -3.1349201), tolerance = 1e-5)
expect_equal(mod$A[1, 1], 0.5770928, tolerance = 1e-5)
expect_equal(mod$A[4, 12], -0.02980421, tolerance = 1e-5)

# BM should give almost the same...
mod_BM = DFM(diff(X_cc), 4, 3, em.method = "BM")

expect_equal(mod$F_twostep, mod_BM$F_twostep, tolerance = 1e-2)
expect_equal(mod$F_qml, mod_BM$F_qml, tolerance = 1e-2)
expect_equal(mod$A, mod_BM$A, tolerance = 1e-1)

# testing #73
y<-as.xts(rnorm(100),order.by = seq(from=as.Date("1980-04-01"),length.out=100,by="quarter")-1)
x1<-as.xts(rnorm(100),order.by = seq(from=as.Date("1980-02-01"),length.out=100,by="month")-1)
x2<-as.xts(rnorm(100),order.by = seq(from=as.Date("1980-02-01"),length.out=100,by="month")-1)

data<-cbind(y,x1,x2)
data["1988-03-31",1]<-NA
data1<-data[time(data)<="1987-12-31"]

expect_visible(DFM(data1[,c("x1","x2","y")], r = 1, p = 1, quarterly.vars = c("y")))


### Block Factor Model ---------------------------------

make_public_block_data <- function(
	n = 120L,
	seed = 123L
) {
	set.seed(seed)

	global <- numeric(n)
	sector <- numeric(n)

	global_shock <- rnorm(n, sd = 0.4)
	sector_shock <- rnorm(n, sd = 0.3)

	for (t in 2:n) {
		global[t] <- 0.55 * global[t - 1L] + global_shock[t]
		sector[t] <- 0.35 * sector[t - 1L] + sector_shock[t]
	}

	aggregate <- rep(NA_real_, n)

	for (t in 5:n) {
		aggregate[t] <- sum(
			.mariano_murasawa_weights *
				global[t:(t - 4L)]
		)
	}

	X <- cbind(
		M1 = 0.9 * global + 0.7 * sector + rnorm(n, sd = 0.15),
		M2 = 0.8 * global - 0.5 * sector + rnorm(n, sd = 0.15),
		M3 = 0.6 * global + rnorm(n, sd = 0.15),
		M4 = -0.7 * global + rnorm(n, sd = 0.15),
		Q1 = 1.1 * aggregate + rnorm(n, sd = 0.15)
	)

	quarter_end <- seq.int(6L, n, by = 3L)
	X[-quarter_end, "Q1"] <- NA_real_

	X
}

make_public_block_factors <- function() {
	list(
		M1 = c("global", "sector"),
		M2 = c("global", "sector"),
		M3 = "global",
		M4 = "global",
		Q1 = "global"
	)
}

fit_public_block_model <- function(
	save.full.state = TRUE,
	min.iter = 2L,
	max.iter = 2L
) {
	X <- make_public_block_data()

	expect_warning(
		fit <- DFM(
			X = X,
			factors = make_public_block_factors(),
			factor_multiplicities = c(
				global = 1L,
				sector = 1L
			),
			factor_orders = c(
				global = 1L,
				sector = 1L
			),
			quarterly.vars = "Q1",
			idio.ar1 = TRUE,
			min.iter = min.iter,
			max.iter = max.iter,
			save.full.state = save.full.state
		),
		"Maximum number of iterations"
	)

	fit
}

test_that("DFM dispatches block specifications to the block estimator", {
	fit <- fit_public_block_model()

	expect_s3_class(fit, "dfm")
	expect_identical(fit$em.method, "Block")
	expect_true(is.list(fit$factor.spec))
	expect_true(is.list(fit$state.layout))

	expect_identical(
		fit$factor.spec$factor_names,
		c("global", "sector")
	)
})

test_that("public block DFM matrices retain the classical dfm contract", {
	fit <- fit_public_block_model()

	expect_identical(dim(fit$A), c(2L, 2L))
	expect_identical(dim(fit$C), c(5L, 2L))
	expect_identical(dim(fit$Q), c(2L, 2L))
	expect_identical(dim(fit$R), c(5L, 5L))

	expect_identical(
		rownames(fit$C),
		c("M1", "M2", "M3", "M4", "Q1")
	)
	expect_identical(
		colnames(fit$C),
		c("global", "sector")
	)
})

test_that("public block DFM preserves structural loading restrictions", {
	fit <- fit_public_block_model()

	expect_equal(fit$C["M3", "sector"], 0)
	expect_equal(fit$C["M4", "sector"], 0)
	expect_equal(fit$C["Q1", "sector"], 0)
})

test_that("public block DFM exposes block factor estimates", {
	fit <- fit_public_block_model()
	n <- nrow(fit$X_imp)

	expect_identical(dim(fit$F_pca), c(n, 2L))
	expect_identical(dim(fit$F_2s), c(n, 2L))
	expect_identical(dim(fit$F_qml), c(n, 2L))

	expect_identical(
		colnames(fit$F_qml),
		c("global", "sector")
	)
	expect_true(all(is.finite(fit$F_qml)))
})

test_that("public block DFM exposes idiosyncratic AR1 estimates", {
	fit <- fit_public_block_model()

	expect_length(fit$rho, 5L)
	expect_identical(
		names(fit$rho),
		c("M1", "M2", "M3", "M4", "Q1")
	)
	expect_true(all(is.finite(fit$rho)))
	expect_true(all(abs(fit$rho) < 1))

	expect_identical(
		dim(fit$e),
		c(nrow(fit$X_imp), 5L)
	)
})

test_that("public block DFM retains the complete block state when requested", {
	fit <- fit_public_block_model()

	expect_named(
		fit$ss_full,
		c(
			"A",
			"C",
			"Q",
			"R",
			"F_0",
			"P_0",
			"F_smooth",
			"P_smooth"
		)
	)

	expect_identical(
		nrow(fit$ss_full$A),
		fit$state.layout$n_states
	)
	expect_identical(
		ncol(fit$ss_full$C),
		fit$state.layout$n_states
	)
})

test_that("public block DFM can omit the complete state", {
	fit <- fit_public_block_model(
		save.full.state = FALSE
	)

	expect_null(fit$ss_full)
})

test_that("block DFM final likelihood corresponds to the reported state", {
	fit <- fit_public_block_model()

	X <- unclass(fit$X_imp)
	missing <- attr(fit$X_imp, "missing")

	if (!is.null(missing)) {
		X[missing] <- NA_real_
	}

	expect_equal(
		tail(fit$loglik, 1L),
		with(
			fit$ss_full,
			dfms::SKFS(
				X,
				A,
				C,
				Q,
				R,
				F_0,
				P_0,
				loglik = TRUE
			)$loglik
		),
		tolerance = 1e-8
	)
})

test_that("block DFM rejects conflicting legacy factor count", {
	X <- make_public_block_data()

	expect_error(
		DFM(
			X = X,
			r = 2L,
			factors = make_public_block_factors(),
			quarterly.vars = "Q1",
			idio.ar1 = TRUE
		),
		"r cannot be supplied"
	)
})

test_that("block DFM currently requires idiosyncratic AR1 states", {
	X <- make_public_block_data()

	expect_error(
		DFM(
			X = X,
			factors = make_public_block_factors(),
			quarterly.vars = "Q1",
			idio.ar1 = FALSE
		),
		"require idio.ar1 = TRUE"
	)
})

test_that("block DFM currently rejects legacy DGR and none EM modes", {
	X <- make_public_block_data()

	expect_error(
		DFM(
			X = X,
			factors = make_public_block_factors(),
			quarterly.vars = "Q1",
			idio.ar1 = TRUE,
			em.method = "DGR"
		),
		"support em.method"
	)

	expect_error(
		DFM(
			X = X,
			factors = make_public_block_factors(),
			quarterly.vars = "Q1",
			idio.ar1 = TRUE,
			em.method = "none"
		),
		"support em.method"
	)
})
