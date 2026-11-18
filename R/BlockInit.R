#' Initialise factors for a block factor model
#'
#' Constructs starting factor estimates using sequential principal components
#' while respecting the zero restrictions in the factor loading matrix.
#'
#' Factors sharing an identical loading pattern are initialised jointly. Loading
#' groups are processed from the most pervasive to the least pervasive, with
#' each successive group extracted from the residual variation left by earlier
#' groups.
#'
#' This is an initialisation procedure rather than the final estimator. It adapts
#' the principal-components initialisation of approximate dynamic factor models
#' to a multilevel loading structure with block zero restrictions.
#'
#' @param X_imp Complete numeric observation matrix with rows representing time
#'   and columns corresponding to the series in spec$loading_mask.
#' @param spec Output from parse_factor_spec().
#'
#' @return A numeric matrix containing the initial factor estimates, with one
#'   column per factor in spec$factor_names.
#'
#' @references
#' Bai, J. (2003). Inferential Theory for Factor Models of Large Dimensions.
#' Econometrica, 71(1), 135-171.
#'
#' Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large
#' approximate dynamic factor models based on Kalman filtering.
#' Journal of Econometrics, 164(1), 188-205.
#'
#' Breitung, J., & Eickmeier, S. (2016). Analyzing International Business and
#' Financial Cycles using Multi-Level Factor Models: A Comparison of Alternative
#' Approaches. Advances in Econometrics, 35, 177-214.
#'
#' Choi, I., Kim, D., Kim, Y. J., & Kwark, N.-S. (2018). A multilevel factor
#' model: Identification, asymptotic theory and applications.
#' Journal of Applied Econometrics, 33(3), 355-377.
#'
#' @keywords internal
initialise_block_factors <- function(
	X_imp,
	spec
) {
	if (
		!is.matrix(X_imp) ||
		!is.numeric(X_imp) ||
		any(!is.finite(X_imp))
	) {
		stop("X_imp must be a finite numeric matrix.")
	}

	if (
		!is.list(spec) ||
		is.null(spec$loading_mask) ||
		is.null(spec$factor_names)
	) {
		stop("Invalid factor specification.")
	}

	loading_mask <- spec$loading_mask
	series <- rownames(loading_mask)
	factor_names <- spec$factor_names

	if (
		ncol(X_imp) != length(series) ||
		is.null(colnames(X_imp)) ||
		!setequal(colnames(X_imp), series)
	) {
		stop(
			"X_imp columns must match the series in the factor specification."
		)
	}

	X_imp <- X_imp[, series, drop = FALSE]

	n_obs <- nrow(X_imp)
	n_factors <- length(factor_names)

	factors <- matrix(
		0,
		nrow = n_obs,
		ncol = n_factors,
		dimnames = list(
			rownames(X_imp),
			factor_names
		)
	)

	residuals <- X_imp

	# Factors with identical loading support are extracted jointly.
	support_keys <- vapply(
		seq_along(factor_names),
		function(j) {
			paste(
				which(loading_mask[, j]),
				collapse = ","
			)
		},
		character(1)
	)

	loading_groups <- split(
		factor_names,
		support_keys
	)

	# Extract the most pervasive factor groups first.
	group_prevalence <- vapply(
		loading_groups,
		function(group) {
			sum(loading_mask[, group[[1L]]])
		},
		integer(1)
	)

	first_position <- vapply(
		loading_groups,
		function(group) {
			match(group[[1L]], factor_names)
		},
		integer(1)
	)

	group_order <- order(
		-group_prevalence,
		first_position
	)

	loading_groups <- loading_groups[group_order]

	for (group in loading_groups) {
		support <- loading_mask[, group[[1L]]]
		group_series <- series[support]
		n_group_factors <- length(group)

		if (length(group_series) < n_group_factors) {
			stop(
				"Loading group has fewer observed series than factors: ",
				paste(group, collapse = ", ")
			)
		}

		group_data <- residuals[
			,
			group_series,
			drop = FALSE
		]

		decomposition <- eigen(
			cov(group_data),
			symmetric = TRUE
		)

		vectors <- decomposition$vectors[
			,
			seq_len(n_group_factors),
			drop = FALSE
		]

		group_factors <- group_data %*% vectors
		colnames(group_factors) <- group
		factors[, group] <- group_factors

		# Temporary loadings are used only to residualise the panel before
		# extracting the next loading group.
		group_loadings <- solve(
			crossprod(group_factors),
			crossprod(
				group_factors,
				group_data
			)
		)

		fitted <- group_factors %*% group_loadings

		residuals[, group_series] <-
			residuals[, group_series, drop = FALSE] - fitted
	}

	factors
}



#' Initialise mixed-frequency factor loadings and residuals
#'
#' Estimates initial factor loadings under the block loading restrictions and
#' the mixed-frequency measurement equation.
#'
#' Monthly series are regressed on their permitted contemporaneous factors.
#' Quarterly series are regressed on Mariano-Murasawa weighted aggregates of
#' their permitted factors,
#'
#'   g_t = f_t + 2 f_{t-1} + 3 f_{t-2} + 2 f_{t-3} + f_{t-4}.
#'
#' This is equivalent to estimating the five-lag quarterly measurement equation
#' subject to the fixed Mariano-Murasawa loading restriction.
#'
#' @param X Numeric observation matrix containing the original missing-value
#'   pattern.
#' @param X_imp Complete numeric observation matrix used for initialization.
#' @param factors Numeric matrix of initial factor estimates, with columns
#'   corresponding to spec$factor_names.
#' @param spec Output from parse_factor_spec().
#'
#' @return A list containing:
#'   \itemize{
#'     \item loadings: Initial restricted factor loading matrix.
#'     \item residuals: Measurement residual matrix preserving the missing-value
#'       pattern of X.
#'     \item fitted: Initial factor contribution to the measurement equation.
#'   }
#'
#' @references
#' Mariano, R. S., & Murasawa, Y. (2003). A new coincident index of business
#' cycles based on monthly and quarterly series.
#' Journal of Applied Econometrics, 18(4), 427-443.
#'
#' Bańbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor
#' models on datasets with arbitrary pattern of missing data.
#' Journal of Applied Econometrics, 29(1), 133-160.
#'
#' @keywords internal
initialise_block_measurement <- function(
	X,
	X_imp,
	factors,
	spec
) {
	if (
		!is.matrix(X) ||
		!is.numeric(X)
	) {
		stop("X must be a numeric matrix.")
	}

	if (
		!is.matrix(X_imp) ||
		!is.numeric(X_imp) ||
		any(!is.finite(X_imp))
	) {
		stop("X_imp must be a finite numeric matrix.")
	}

	if (
		!identical(dim(X), dim(X_imp)) ||
		!identical(colnames(X), colnames(X_imp))
	) {
		stop("X and X_imp must have matching dimensions and columns.")
	}

	if (
		!is.list(spec) ||
		is.null(spec$loading_mask) ||
		is.null(spec$factor_names) ||
		is.null(spec$monthly_vars) ||
		is.null(spec$quarterly_vars)
	) {
		stop("Invalid factor specification.")
	}

	series <- rownames(spec$loading_mask)
	factor_names <- spec$factor_names

	if (
		is.null(colnames(X)) ||
		!setequal(colnames(X), series)
	) {
		stop(
			"Observation columns must match the factor specification."
		)
	}

	if (
		!is.matrix(factors) ||
		!is.numeric(factors) ||
		any(!is.finite(factors)) ||
		nrow(factors) != nrow(X) ||
		is.null(colnames(factors)) ||
		!setequal(colnames(factors), factor_names)
	) {
		stop(
			"factors must be a finite matrix matching the factor specification."
		)
	}

	X <- X[, series, drop = FALSE]
	X_imp <- X_imp[, series, drop = FALSE]

	factors <- factors[
		,
		factor_names,
		drop = FALSE
	]

	n_obs <- nrow(X)
	n_factors <- length(factor_names)

	loadings <- matrix(
		0,
		nrow = length(series),
		ncol = n_factors,
		dimnames = list(
			series,
			factor_names
		)
	)

	fitted <- matrix(
		NA_real_,
		nrow = n_obs,
		ncol = length(series),
		dimnames = list(
			rownames(X),
			series
		)
	)

	# Monthly measurement equations
	for (series_name in spec$monthly_vars) {
		allowed <- factor_names[
			spec$loading_mask[series_name, ]
		]

		F_allowed <- factors[
			,
			allowed,
			drop = FALSE
		]

		# Initial loading estimation uses the complete initialization panel.
		y <- X_imp[, series_name]

		lambda <- solve(
			crossprod(F_allowed),
			crossprod(F_allowed, y)
		)

		loadings[
			series_name,
			allowed
		] <- as.vector(lambda)

		fitted[, series_name] <- as.vector(
			F_allowed %*% lambda
		)
	}

	# Quarterly measurement equations
	if (length(spec$quarterly_vars) > 0L) {
		weights <- .mariano_murasawa_weights
		n_lags <- length(weights)

		if (n_obs < n_lags) {
			stop(
				"Insufficient observations for quarterly aggregation."
			)
		}

		# Weighted factor aggregates aligned with the current month.
		#
		# Rows before the fifth observation cannot have a complete
		# Mariano-Murasawa factor history and remain NA.
		aggregated_factors <- matrix(
			NA_real_,
			nrow = n_obs,
			ncol = n_factors,
			dimnames = list(
				rownames(factors),
				factor_names
			)
		)

		for (t in n_lags:n_obs) {
			window <- factors[
				t:(t - n_lags + 1L),
				,
				drop = FALSE
			]

			aggregated_factors[t, ] <- as.vector(
				crossprod(
					weights,
					window
				)
			)
		}

		for (series_name in spec$quarterly_vars) {
			allowed <- factor_names[
				spec$loading_mask[series_name, ]
			]

			G_allowed <- aggregated_factors[
				,
				allowed,
				drop = FALSE
			]

			# Prefer actual quarterly observations.
			observed <- is.finite(X[, series_name]) &
				apply(
					G_allowed,
					1L,
					function(x) all(is.finite(x))
				)

			# As in the existing dfms initialization, fall back to the
			# completed panel when too few actual observations are available.
			if (sum(observed) < length(allowed) + 2L) {
				observed <- apply(
					G_allowed,
					1L,
					function(x) all(is.finite(x))
				)

				y <- X_imp[
					observed,
					series_name
				]
			} else {
				y <- X[
					observed,
					series_name
				]
			}

			G_fit <- G_allowed[
				observed,
				,
				drop = FALSE
			]

			if (nrow(G_fit) <= length(allowed)) {
				stop(
					"Insufficient quarterly observations to initialise ",
					"loadings for series: ",
					series_name
				)
			}

			lambda <- solve(
				crossprod(G_fit),
				crossprod(G_fit, y)
			)

			loadings[
				series_name,
				allowed
			] <- as.vector(lambda)

			valid <- apply(
				G_allowed,
				1L,
				function(x) all(is.finite(x))
			)

			fitted[
				valid,
				series_name
			] <- as.vector(
				G_allowed[
					valid,
					,
					drop = FALSE
				] %*% lambda
			)
		}
	}

	# Preserve the actual missing-data pattern when constructing residuals.
	residuals <- X - fitted

	list(
		loadings = loadings,
		residuals = residuals,
		fitted = fitted
	)
}




#' Initialise factor dynamics for a block factor model
#'
#' Estimates an independent VAR for each dynamic factor block using the initial
#' factor estimates. Each block is estimated at its specified autoregressive
#' order, and the resulting coefficients and innovation covariance matrices are
#' returned in the format required by build_factor_matrices().
#'
#' The factor VARs are estimated by ordinary least squares without an intercept,
#' following the factor-dynamics initialisation used by dfms. Dynamic blocks are
#' estimated independently, imposing zero cross-block transition coefficients.
#'
#' @param factors Numeric matrix containing initial factor estimates, with one
#'   row per observation and columns corresponding to spec$factor_names.
#' @param spec Output from parse_factor_spec().
#'
#' @return A list containing:
#'     coefficients: Named list of VAR coefficient arrays. For block b,
#'       the array has dimensions (r_b, r_b, p_b), with one matrix per lag.
#'     covariances: Named list of factor innovation covariance matrices.
#'     residuals: Named list of VAR residual matrices.
#'
#' @references
#' Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large
#' approximate dynamic factor models based on Kalman filtering.
#' Journal of Econometrics, 164(1), 188-205.
#'
#' Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood
#' approach for large, approximate dynamic factor models.
#' Review of Economics and Statistics, 94(4), 1014-1024.
#'
#' Lütkepohl, H. (2005). New Introduction to Multiple Time Series Analysis.
#' Springer.
#'
#' @keywords internal
initialise_block_dynamics <- function(
	factors,
	spec
) {
	if (
		!is.matrix(factors) ||
		!is.numeric(factors) ||
		any(!is.finite(factors))
	) {
		stop("factors must be a finite numeric matrix.")
	}

	if (
		!is.list(spec) ||
		is.null(spec$blocks) ||
		is.null(spec$block_orders) ||
		is.null(spec$factor_names)
	) {
		stop("Invalid factor specification.")
	}

	if (
		is.null(colnames(factors)) ||
		anyDuplicated(colnames(factors)) ||
		!setequal(colnames(factors), spec$factor_names)
	) {
		stop(
			"Factor columns must match the factor specification."
		)
	}

	factors <- factors[
		,
		spec$factor_names,
		drop = FALSE
	]

	blocks <- spec$blocks
	block_orders <- spec$block_orders

	if (!identical(names(blocks), names(block_orders))) {
		stop("Block names and orders must agree.")
	}

	coefficients <- vector(
		"list",
		length(blocks)
	)

	covariances <- vector(
		"list",
		length(blocks)
	)

	residuals <- vector(
		"list",
		length(blocks)
	)

	names(coefficients) <- names(blocks)
	names(covariances) <- names(blocks)
	names(residuals) <- names(blocks)

	for (block in names(blocks)) {
		block_factors <- blocks[[block]]
		p <- block_orders[[block]]
		r <- length(block_factors)

		if (nrow(factors) <= p) {
			stop(
				"Insufficient observations to estimate VAR(",
				p,
				") for block: ",
				block
			)
		}

		F_block <- factors[
			,
			block_factors,
			drop = FALSE
		]

		var_fit <- .VAR(
			F_block,
			p = p
		)

		# .VAR() returns an (r * p) x r coefficient matrix:
		#
		#        response factors
		# lag 1
		# lag 2
		# ...
		#
		# build_factor_matrices() expects an r x r x p array where
		# slice j contains the coefficient matrix for lag j.
		coef_matrix <- t(var_fit$A)

		coef_array <- array(
			0,
			dim = c(r, r, p),
			dimnames = list(
				block_factors,
				block_factors,
				paste0("lag", seq_len(p))
			)
		)

		for (lag in seq_len(p)) {
			cols <- ((lag - 1L) * r + 1L):(lag * r)

			coef_array[, , lag] <- coef_matrix[
				,
				cols,
				drop = FALSE
			]
		}

		Sigma <- cov(var_fit$res)

		if (r == 1L) {
			Sigma <- matrix(
				Sigma,
				nrow = 1L,
				ncol = 1L,
				dimnames = list(
					block_factors,
					block_factors
				)
			)
		}

		coefficients[[block]] <- coef_array
		covariances[[block]] <- Sigma
		residuals[[block]] <- var_fit$res
	}

	list(
		coefficients = coefficients,
		covariances = covariances,
		residuals = residuals
	)
}





#' Initialise idiosyncratic dynamics
#'
#' Constructs starting parameters for the monthly latent idiosyncratic
#' processes used by the mixed-frequency state-space model.
#'
#' Monthly idiosyncratic components are initialised as AR(1) processes estimated
#' directly from consecutive monthly measurement residuals.
#'
#' Quarterly idiosyncratic components are latent monthly processes but are only
#' observed through Mariano-Murasawa temporal aggregation. Their persistence is
#' therefore not estimated directly from the sparse quarterly residuals.
#' Quarterly persistence is initialised at zero, and the latent innovation
#' variance is obtained from the variance of the aggregated residual divided by
#' the sum of squared Mariano-Murasawa weights.
#'
#' These are starting values only. Quarterly persistence and innovation
#' variances are subsequently estimated by the EM algorithm using the smoothed
#' latent monthly states.
#'
#' @param residuals Numeric matrix of measurement-equation residuals.
#' @param layout Output from build_state_layout().
#' @param min_variance Positive lower bound for initial innovation variances.
#' @param stationarity_bound Absolute bound applied to estimated monthly AR(1)
#'   coefficients to ensure stationary starting values.
#'
#' @return A list containing named vectors:
#'   \itemize{
#'     \item monthly_rho
#'     \item monthly_variances
#'     \item quarterly_rho
#'     \item quarterly_variances
#'   }
#'
#' @references
#' Mariano, R. S., & Murasawa, Y. (2003). A new coincident index of business
#' cycles based on monthly and quarterly series.
#' Journal of Applied Econometrics, 18(4), 427-443.
#'
#' Bańbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor
#' models on datasets with arbitrary pattern of missing data.
#' Journal of Applied Econometrics, 29(1), 133-160.
#'
#' @keywords internal
initialise_idiosyncratic_dynamics <- function(
	residuals,
	layout,
	min_variance = 1e-6,
	stationarity_bound = 0.99
) {
	if (
		!is.matrix(residuals) ||
		!is.numeric(residuals)
	) {
		stop("residuals must be a numeric matrix.")
	}

	if (
		!is.numeric(min_variance) ||
		length(min_variance) != 1L ||
		!is.finite(min_variance) ||
		min_variance <= 0
	) {
		stop("min_variance must be a positive finite scalar.")
	}

	if (
		!is.numeric(stationarity_bound) ||
		length(stationarity_bound) != 1L ||
		!is.finite(stationarity_bound) ||
		stationarity_bound <= 0 ||
		stationarity_bound >= 1
	) {
		stop(
			"stationarity_bound must lie strictly between zero and one."
		)
	}

	series <- c(
		layout$monthly_vars,
		layout$quarterly_vars
	)

	if (
		is.null(colnames(residuals)) ||
		anyDuplicated(colnames(residuals)) ||
		!setequal(colnames(residuals), series)
	) {
		stop(
			"Residual columns must match the observed series in the layout."
		)
	}

	residuals <- residuals[
		,
		series,
		drop = FALSE
	]

	monthly_rho <- setNames(
		numeric(length(layout$monthly_vars)),
		layout$monthly_vars
	)

	monthly_variances <- setNames(
		numeric(length(layout$monthly_vars)),
		layout$monthly_vars
	)

	quarterly_rho <- setNames(
		numeric(length(layout$quarterly_vars)),
		layout$quarterly_vars
	)

	quarterly_variances <- setNames(
		numeric(length(layout$quarterly_vars)),
		layout$quarterly_vars
	)

	for (series_name in layout$monthly_vars) {
		e <- residuals[, series_name]

		current <- e[-1L]
		lagged <- e[-length(e)]

		observed <- is.finite(current) &
			is.finite(lagged)

		current <- current[observed]
		lagged <- lagged[observed]

		if (length(current) < 2L) {
			stop(
				"Insufficient residual observations for series: ",
				series_name
			)
		}

		denominator <- sum(lagged^2)

		if (
			!is.finite(denominator) ||
			denominator <= 0
		) {
			stop(
				"Cannot estimate idiosyncratic AR(1) for series: ",
				series_name
			)
		}

		rho <- sum(current * lagged) / denominator

		if (!is.finite(rho)) {
			stop(
				"Invalid idiosyncratic AR(1) estimate for series: ",
				series_name
			)
		}

		rho <- max(
			-stationarity_bound,
			min(stationarity_bound, rho)
		)

		innovations <- current - rho * lagged
		variance <- var(innovations)

		if (!is.finite(variance)) {
			stop(
				"Invalid idiosyncratic innovation variance for series: ",
				series_name
			)
		}

		monthly_rho[[series_name]] <- rho
		monthly_variances[[series_name]] <- max(
			variance,
			min_variance
		)
	}

	weight_scale <- sum(.mariano_murasawa_weights^2)

	for (series_name in layout$quarterly_vars) {
		e <- residuals[, series_name]
		e <- e[is.finite(e)]

		if (length(e) < 2L) {
			stop(
				"Insufficient quarterly residual observations for series: ",
				series_name
			)
		}

		aggregate_variance <- var(e)

		if (!is.finite(aggregate_variance)) {
			stop(
				"Invalid quarterly residual variance for series: ",
				series_name
			)
		}

		quarterly_rho[[series_name]] <- 0
		quarterly_variances[[series_name]] <- max(
			aggregate_variance / weight_scale,
			min_variance
		)
	}

	list(
		monthly_rho = monthly_rho,
		monthly_variances = monthly_variances,
		quarterly_rho = quarterly_rho,
		quarterly_variances = quarterly_variances
	)
}



#' Initialise a block dynamic factor model state
#'
#' Constructs the complete initial state-space representation for a block
#' dynamic factor model.
#'
#' Initial factors are extracted using the block-factor initialisation
#' procedure, after which mixed-frequency factor loadings, block-specific factor
#' dynamics, and idiosyncratic dynamics are estimated. These numerical
#' parameters are then assembled into the complete state-space matrices.
#'
#' The initial state mean is set to zero and the initial state covariance is
#' obtained from the stationary covariance implied by the initial transition
#' and innovation covariance matrices.
#'
#' @param X Numeric observation matrix containing the original missing-data
#'   pattern.
#' @param X_imp Complete numeric observation matrix used for factor and
#'   parameter initialisation.
#' @param spec Output from parse_factor_spec().
#' @param layout Output from build_state_layout().
#' @param loading_maps Output from build_loading_maps().
#' @param observation_noise Positive scalar used for the fixed measurement
#'   error covariance.
#' @param state_tolerance Positive numerical tolerance passed to
#'   initialise_full_state().
#' @param min_idio_variance Positive lower bound for initial idiosyncratic
#'   innovation variances.
#' @param stationarity_bound Absolute bound applied to estimated monthly
#'   idiosyncratic AR(1) coefficients.
#'
#' @return A list containing the initial state-space parameters A, C, Q, R,
#'   F_0, and P_0.
#'
#' @references
#' Bai, J. (2003). Inferential Theory for Factor Models of Large Dimensions.
#' Econometrica, 71(1), 135-171.
#'
#' Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large
#' approximate dynamic factor models based on Kalman filtering.
#' Journal of Econometrics, 164(1), 188-205.
#'
#' Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood
#' approach for large, approximate dynamic factor models.
#' Review of Economics and Statistics, 94(4), 1014-1024.
#'
#' Mariano, R. S., & Murasawa, Y. (2003). A new coincident index of business
#' cycles based on monthly and quarterly series.
#' Journal of Applied Econometrics, 18(4), 427-443.
#'
#' Bańbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor
#' models on datasets with arbitrary pattern of missing data.
#' Journal of Applied Econometrics, 29(1), 133-160.
#'
#' Breitung, J., & Eickmeier, S. (2016). Analyzing International Business and
#' Financial Cycles using Multi-Level Factor Models: A Comparison of Alternative
#' Approaches. Advances in Econometrics, 35, 177-214.
#'
#' Choi, I., Kim, D., Kim, Y. J., & Kwark, N.-S. (2018). A multilevel factor
#' model: Identification, asymptotic theory and applications.
#' Journal of Applied Econometrics, 33(3), 355-377.
#'
#' @keywords internal
initialise_block_state <- function(
	X,
	X_imp,
	spec,
	layout,
	loading_maps,
	observation_noise = 1e-6,
	state_tolerance = 1e-8,
	min_idio_variance = 1e-6,
	stationarity_bound = 0.99
) {
	if (!isTRUE(layout$idio_ar1)) {
		stop(
			"Block state initialisation currently requires idio_ar1 = TRUE."
		)
	}

	# Initial factor paths
	factors <- initialise_block_factors(
		X_imp = X_imp,
		spec = spec
	)

	# Mixed-frequency measurement equations
	measurement <- initialise_block_measurement(
		X = X,
		X_imp = X_imp,
		factors = factors,
		spec = spec
	)

	# Block-specific factor dynamics
	factor_dynamics <- initialise_block_dynamics(
		factors = factors,
		spec = spec
	)

	# Monthly and quarterly idiosyncratic dynamics
	idio_dynamics <- initialise_idiosyncratic_dynamics(
		residuals = measurement$residuals,
		layout = layout,
		min_variance = min_idio_variance,
		stationarity_bound = stationarity_bound
	)

	# Factor transition and innovation covariance matrices
	factor_matrices <- build_factor_matrices(
		layout = layout$factor,
		coefficients = factor_dynamics$coefficients,
		covariances = factor_dynamics$covariances
	)

	# Factor portion of the measurement matrix
	factor_measurement <- build_factor_measurement(
		spec = spec,
		layout = layout$factor,
		loading_maps = loading_maps,
		loadings = measurement$loadings
	)

	# Complete A, C, Q, and R matrices
	state_matrices <- build_full_state_matrices(
		spec = spec,
		layout = layout,
		factor_matrices = factor_matrices,
		factor_measurement = factor_measurement,
		monthly_rho = idio_dynamics$monthly_rho,
		monthly_variances = idio_dynamics$monthly_variances,
		quarterly_rho = idio_dynamics$quarterly_rho,
		quarterly_variances = idio_dynamics$quarterly_variances,
		observation_noise = observation_noise
	)

	# Stationary initial state moments
	initial_state <- initialise_full_state(
		A = state_matrices$A,
		Q = state_matrices$Q,
		tolerance = state_tolerance
	)

	c(
		state_matrices,
		initial_state
	)
}
