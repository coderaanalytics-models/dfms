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
    PMI = list(indices = matrix(c(1L, 3L), 2L), weights = 1, idio_indices = 5L),
    IP = list(indices = matrix(2L, 1L), weights = 1, idio_indices = 6L),
    Restricted = list(indices = matrix(integer(), 0L, 1L), weights = 1, idio_indices = 4L)
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
  expect_error(
      assemble_loadings(
          setNames(
              list(1, 2, numeric()),
              rownames(mask)
          ),
          mask
      ),
      "Loading dimensions do not match"
  )
})

test_that("loading updates reject malformed maps and singular normal equations", {
  moments <- make_block_moments()
  X <- matrix(seq_len(8))
  run <- function(maps, P = moments$P_smooth) {
    update_loadings_cpp(X, moments$F_smooth, P, maps)
  }
  map <- list(indices = matrix(1L, 1L), weights = 1, idio_indices = 2L)
  expect_error(run(list()), "Number of loading maps")
  expect_error(run(list(map), moments$P_smooth[, , -1L, drop = FALSE]),
               "Incompatible smoother dimensions")
  expect_error(run(list(modifyList(map, list(weights = c(1, 2))))),
               "Incompatible loading-map dimensions")
  expect_error(run(list(modifyList(map, list(idio_indices = 0L)))),
               "Invalid idiosyncratic state index")
  for (index in c(0L, 7L, NA_integer_)) {
    expect_error(run(list(modifyList(map, list(indices = matrix(index, 1L))))),
                 "Invalid factor-state index")
  }
  X[] <- NA_real_
  expect_error(run(list(map)), "Failed to solve loading normal equations")
  expect_error(
      run(
          list(
              modifyList(
                  map,
                  list(
                      weights = c(1, 2),
                      indices = matrix(
                          c(1L, 2L),
                          nrow = 1L
                      ),
                      idio_indices = 2L
                  )
              )
          )
      ),
      "idiosyncratic loading-map dimensions"
  )
})

# Mixed-frequency integration

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

    quarter_end <- seq.int(3L, n, by = 3L)

    for (variable in quarterly_vars) {
        X[-quarter_end, variable] <- NA_real_
    }

    monthly_vars <- setdiff(
        colnames(X),
        quarterly_vars
    )

    if (
        missing_rate > 0 &&
        length(monthly_vars) > 0L
    ) {
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

make_block_em_case <- function() {
    model <- make_test_full_state()

    initial <- initialise_full_state(
        model$matrices$A,
        model$matrices$Q
    )

    state <- c(
        model$matrices,
        initial
    )

    simulated <- simulate_full_state(
        model$matrices,
        initial
    )

    list(
        X = make_mixed_frequency_observations(
            simulated$observations
        ),
        state = state,
        spec = model$spec,
        layout = model$layout,
        maps = model$loading_maps
    )
}

smooth_block_em_case <- function(model) {
    with(
        model$state,
        dfms::SKFS(
            model$X,
            A,
            C,
            Q,
            R,
            F_0,
            P_0,
            loglik = TRUE
        )
    )
}

test_that("SKFS supplies finite EM moments for a ragged mixed-frequency panel", {
    model <- make_block_em_case()

    expect_equal(
        which(!is.na(model$X[, "GDP"])),
        seq.int(3L, 120L, by = 3L)
    )
    expect_true(anyNA(model$X[, "PMI"]))
    expect_true(any(diag(model$state$Q) == 0))

    smooth <- smooth_block_em_case(model)

    n <- nrow(model$state$A)
    periods <- nrow(model$X)

    expect_identical(
        dim(smooth$F_smooth),
        c(periods, n)
    )
    expect_identical(
        dim(smooth$P_smooth),
        c(n, n, periods)
    )
    expect_identical(
        dim(smooth$PPm_smooth),
        c(n, n, periods)
    )
    expect_length(smooth$F_smooth_0, n)
    expect_identical(
        dim(smooth$P_smooth_0),
        c(n, n)
    )

    fields <- c(
        "F_smooth",
        "P_smooth",
        "PPm_smooth",
        "F_smooth_0",
        "P_smooth_0",
        "loglik"
    )

    expect_true(
        all(is.finite(unlist(smooth[fields])))
    )
})

test_that("quarterly loadings account for all weighted idiosyncratic lags", {
    model <- make_block_em_case()
    moments <- smooth_block_em_case(model)

    actual <- update_loadings_cpp(
        model$X,
        moments$F_smooth,
        moments$P_smooth,
        model$maps
    )

    mask <- model$spec$loading_mask["GDP", ]
    allowed <- names(mask)[mask]

    n <- nrow(model$state$A)
    H <- matrix(0, length(allowed), n)
    weights <- .mariano_murasawa_weights

    for (j in seq_along(allowed)) {
        block <- names(
            Filter(
                function(factors) {
                    allowed[[j]] %in% factors
                },
                model$spec$blocks
            )
        )

        position <- match(
            allowed[[j]],
            model$spec$blocks[[block]]
        )

        lags <- model$layout$factor$lag_indices[[block]][
            seq_along(weights)
        ]

        H[
            j,
            vapply(
                lags,
                function(lag) lag[[position]],
                integer(1)
            )
        ] <- weights
    }

    d <- numeric(n)
    idio <- model$layout$quarterly_indices[["GDP"]]
    d[idio] <- model$state$C["GDP", idio]

    expect_loading_equation(
        actual$GDP,
        model$X[, "GDP"],
        moments,
        H,
        d
    )
})

test_that("EMstepBlock assembles a complete state update", {
    model <- make_block_em_case()

    actual <- EMstepBlock(
        X = model$X,
        state = model$state,
        spec = model$spec,
        layout = model$layout,
        loading_maps = model$maps
    )

    expect_named(
        actual,
        c(
            "state",
            "loglik",
            "smooth",
            "factor_updates",
            "idio_updates",
            "loading_updates",
            "loadings"
        )
    )

    expect_named(
        actual$state,
        c("A", "C", "Q", "R", "F_0", "P_0")
    )

    expect_true(is.finite(actual$loglik))
    expect_true(is.finite(actual$smooth$loglik))

    expect_named(
        actual$idio_updates$coefficients,
        names(model$layout$idio_current)
    )
    expect_named(
        actual$factor_updates$coefficients,
        names(model$spec$blocks)
    )

    expect_identical(
        dimnames(actual$loadings),
        dimnames(model$spec$loading_mask)
    )
    expect_true(all(is.finite(actual$loadings)))
    expect_true(
        all(
            actual$loadings[
                !model$spec$loading_mask
            ] == 0
        )
    )

    expect_identical(
        dim(actual$state$A),
        dim(model$state$A)
    )
    expect_identical(
        dim(actual$state$C),
        dim(model$state$C)
    )
    expect_identical(
        dim(actual$state$Q),
        dim(model$state$Q)
    )
})

test_that("EMstepBlock preserves fixed initial state moments", {
    model <- make_block_em_case()

    actual <- EMstepBlock(
        X = model$X,
        state = model$state,
        spec = model$spec,
        layout = model$layout,
        loading_maps = model$maps
    )

    expect_identical(
        actual$state$F_0,
        model$state$F_0
    )
    expect_identical(
        actual$state$P_0,
        model$state$P_0
    )
})

test_that("EMstepBlock preserves fixed observation noise", {
    model <- make_block_em_case()

    actual <- EMstepBlock(
        X = model$X,
        state = model$state,
        spec = model$spec,
        layout = model$layout,
        loading_maps = model$maps
    )

    expect_equal(
        actual$state$R,
        model$state$R
    )
})

test_that("EMstepBlock preserves loading restrictions in the measurement matrix", {
    model <- make_block_em_case()

    actual <- EMstepBlock(
        X = model$X,
        state = model$state,
        spec = model$spec,
        layout = model$layout,
        loading_maps = model$maps
    )

    for (series_name in rownames(model$spec$loading_mask)) {
        forbidden <- model$spec$factor_names[
            !model$spec$loading_mask[series_name, ]
        ]

        for (factor in forbidden) {
            factor_states <- model$layout$factor$factor_lag_indices[[factor]]

            expect_true(
                all(
                    actual$state$C[
                        series_name,
                        factor_states
                    ] == 0
                )
            )
        }
    }
})

test_that("EMBlock runs repeated updates and returns a final smoother", {
    model <- make_block_em_case()

    expect_warning(
        fit <- EMBlock(
            X = model$X,
            state = model$state,
            spec = model$spec,
            layout = model$layout,
            loading_maps = model$maps,
            min.iter = 2L,
            max.iter = 2L,
            tol = 1e-8
        ),
        "Maximum number of iterations"
    )

    expect_identical(fit$iterations, 2L)
    expect_false(fit$converged)
    expect_length(fit$loglik, 2L)
    expect_true(all(is.finite(fit$loglik)))

    expect_named(
        fit$state,
        c("A", "C", "Q", "R", "F_0", "P_0")
    )
    expect_true(is.finite(fit$smooth$loglik))
    expect_true(all(is.finite(fit$smooth$F_smooth)))
})

test_that("EMBlock preserves initial moments across repeated updates", {
    model <- make_block_em_case()

    expect_warning(
        fit <- EMBlock(
            X = model$X,
            state = model$state,
            spec = model$spec,
            layout = model$layout,
            loading_maps = model$maps,
            min.iter = 2L,
            max.iter = 2L
        ),
        "Maximum number of iterations"
    )

    expect_identical(
        fit$state$F_0,
        model$state$F_0
    )
    expect_identical(
        fit$state$P_0,
        model$state$P_0
    )
})

test_that("EMBlock final smoother uses the final updated state", {
    model <- make_block_em_case()

    expect_warning(
        fit <- EMBlock(
            X = model$X,
            state = model$state,
            spec = model$spec,
            layout = model$layout,
            loading_maps = model$maps,
            min.iter = 2L,
            max.iter = 2L
        ),
        "Maximum number of iterations"
    )

    expected <- with(
        fit$state,
        dfms::SKFS(
            model$X,
            A,
            C,
            Q,
            R,
            F_0,
            P_0,
            loglik = TRUE
        )
    )

    expect_equal(
        fit$smooth$loglik,
        expected$loglik,
        tolerance = 1e-10
    )
    expect_equal(
        fit$smooth$F_smooth,
        expected$F_smooth,
        tolerance = 1e-10
    )
})

test_that("EMBlock does not converge before min.iter", {
    model <- make_block_em_case()

    expect_warning(
        fit <- EMBlock(
            X = model$X,
            state = model$state,
            spec = model$spec,
            layout = model$layout,
            loading_maps = model$maps,
            min.iter = 3L,
            max.iter = 2L,
            tol = 1e10
        ),
        "Maximum number of iterations"
    )

    expect_identical(fit$iterations, 2L)
    expect_false(fit$converged)
})
