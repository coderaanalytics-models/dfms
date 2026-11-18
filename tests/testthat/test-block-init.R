make_simple_block_init_fixture <- function(
    n = 120L,
    seed = 123L
) {
    set.seed(seed)

    spec <- parse_factor_spec(
        series = c("M1", "M2", "M3", "Q1"),
        factors = "global",
        factor_orders = 1L,
        quarterly_vars = "Q1"
    )

    factor <- numeric(n)
    shocks <- rnorm(n, sd = 0.5)

    for (t in 2:n) {
        factor[t] <- 0.5 * factor[t - 1L] + shocks[t]
    }

    aggregated_factor <- rep(NA_real_, n)
    weights <- .mariano_murasawa_weights

    for (t in length(weights):n) {
        aggregated_factor[t] <- sum(
            weights * factor[t:(t - length(weights) + 1L)]
        )
    }

    X_imp <- cbind(
        M1 = 1.0 * factor + rnorm(n, sd = 0.15),
        M2 = -0.8 * factor + rnorm(n, sd = 0.15),
        M3 = 0.6 * factor + rnorm(n, sd = 0.15),
        Q1 = ifelse(
            is.finite(aggregated_factor),
            1.2 * aggregated_factor + rnorm(n, sd = 0.15),
            rnorm(n, sd = 0.15)
        )
    )

    X <- X_imp
    quarter_end <- seq.int(6L, n, by = 3L)
    X[-quarter_end, "Q1"] <- NA_real_

    layout <- build_state_layout(
        spec,
        idio_ar1 = TRUE
    )

    loading_maps <- build_loading_maps(
        spec,
        layout
    )

    list(
        X = X,
        X_imp = X_imp,
        spec = spec,
        layout = layout,
        loading_maps = loading_maps,
        factor = factor
    )
}

test_that("block factor initialisation returns finite factors in specification order", {
    fixture <- make_simple_block_init_fixture()

    factors <- initialise_block_factors(
        fixture$X_imp,
        fixture$spec
    )

    expect_true(is.matrix(factors))
    expect_identical(
        dim(factors),
        c(nrow(fixture$X_imp), fixture$spec$n_factors)
    )
    expect_identical(
        colnames(factors),
        fixture$spec$factor_names
    )
    expect_true(all(is.finite(factors)))
})

test_that("block factor initialisation preserves row names", {
    fixture <- make_simple_block_init_fixture()
    rownames(fixture$X_imp) <- paste0("t", seq_len(nrow(fixture$X_imp)))

    factors <- initialise_block_factors(
        fixture$X_imp,
        fixture$spec
    )

    expect_identical(
        rownames(factors),
        rownames(fixture$X_imp)
    )
})

test_that("block factor initialisation rejects incomplete data", {
    fixture <- make_simple_block_init_fixture()
    fixture$X_imp[1, 1] <- NA_real_

    expect_error(
        initialise_block_factors(
            fixture$X_imp,
            fixture$spec
        ),
        "finite numeric matrix"
    )
})

test_that("mixed-frequency measurement recovers exact monthly and quarterly loadings", {
    n <- 15L

    spec <- parse_factor_spec(
        series = c("M1", "M2", "Q1"),
        factors = "global",
        quarterly_vars = "Q1"
    )

    factors <- matrix(
        seq_len(n),
        ncol = 1L,
        dimnames = list(NULL, "global")
    )

    weights <- .mariano_murasawa_weights
    aggregate <- rep(NA_real_, n)

    for (t in length(weights):n) {
        aggregate[t] <- sum(
            weights * factors[
                t:(t - length(weights) + 1L),
                "global"
            ]
        )
    }

    X_imp <- cbind(
        M1 = 1.5 * factors[, "global"],
        M2 = -0.5 * factors[, "global"],
        Q1 = ifelse(is.finite(aggregate), 2 * aggregate, 0)
    )

    X <- X_imp
    X[, "Q1"] <- NA_real_
    X[c(6L, 9L, 12L, 15L), "Q1"] <- X_imp[
        c(6L, 9L, 12L, 15L),
        "Q1"
    ]

    result <- initialise_block_measurement(
        X,
        X_imp,
        factors,
        spec
    )

    expect_equal(result$loadings["M1", "global"], 1.5)
    expect_equal(result$loadings["M2", "global"], -0.5)
    expect_equal(result$loadings["Q1", "global"], 2)

    expect_equal(
        result$residuals[
            is.finite(X[, "Q1"]),
            "Q1"
        ],
        rep(0, sum(is.finite(X[, "Q1"]))),
        tolerance = 1e-10
    )
})

test_that("mixed-frequency measurement preserves structural loading zeros", {
    spec <- parse_factor_spec(
        series = c("M1", "M2", "M3", "Q1"),
        factors = list(
            M1 = c("global", "sector"),
            M2 = c("global", "sector"),
            M3 = "global",
            Q1 = "global"
        ),
        quarterly_vars = "Q1"
    )

    n <- 30L
    factors <- cbind(
        global = sin(seq_len(n) / 4),
        sector = cos(seq_len(n) / 5)
    )

    X_imp <- cbind(
        M1 = factors[, "global"] + factors[, "sector"],
        M2 = factors[, "global"] - factors[, "sector"],
        M3 = factors[, "global"],
        Q1 = rep(0, n)
    )

    aggregate <- rep(NA_real_, n)

    for (t in 5:n) {
        aggregate[t] <- sum(
            .mariano_murasawa_weights *
                factors[t:(t - 4L), "global"]
        )
    }

    X_imp[is.finite(aggregate), "Q1"] <- aggregate[is.finite(aggregate)]

    X <- X_imp
    X[, "Q1"] <- NA_real_
    quarter_end <- seq.int(6L, n, by = 3L)
    X[quarter_end, "Q1"] <- X_imp[quarter_end, "Q1"]

    result <- initialise_block_measurement(
        X,
        X_imp,
        factors,
        spec
    )

    expect_equal(result$loadings["M3", "sector"], 0)
    expect_equal(result$loadings["Q1", "sector"], 0)
})

test_that("block dynamics return one VAR parameter set per dynamic block", {
    set.seed(42)

    spec <- make_block_spec()
    factors <- matrix(
        rnorm(200L * spec$n_factors),
        nrow = 200L,
        ncol = spec$n_factors,
        dimnames = list(NULL, spec$factor_names)
    )

    result <- initialise_block_dynamics(
        factors,
        spec
    )

    expect_identical(
        names(result$coefficients),
        names(spec$blocks)
    )
    expect_identical(
        names(result$covariances),
        names(spec$blocks)
    )

    for (block in names(spec$blocks)) {
        r <- length(spec$blocks[[block]])
        p <- spec$block_orders[[block]]

        expect_identical(
            dim(result$coefficients[[block]]),
            c(r, r, p)
        )
        expect_identical(
            dim(result$covariances[[block]]),
            c(r, r)
        )
        expect_true(
            isTRUE(all.equal(
                result$covariances[[block]],
                t(result$covariances[[block]])
            ))
        )
    }
})

test_that("block dynamics reject factor columns inconsistent with the specification", {
    spec <- make_block_spec()
    factors <- matrix(
        rnorm(100L * spec$n_factors),
        nrow = 100L,
        dimnames = list(
            NULL,
            paste0("wrong", seq_len(spec$n_factors))
        )
    )

    expect_error(
        initialise_block_dynamics(
            factors,
            spec
        ),
        "must match"
    )
})

test_that("monthly idiosyncratic persistence is estimated from consecutive residuals", {
    layout <- build_state_layout(
        parse_factor_spec(
            series = c("M1", "Q1"),
            factors = "global",
            quarterly_vars = "Q1"
        ),
        idio_ar1 = TRUE
    )

    residuals <- matrix(
        NA_real_,
        nrow = 8L,
        ncol = 2L,
        dimnames = list(NULL, c("M1", "Q1"))
    )

    residuals[, "M1"] <- 2 * 0.5^(0:7)
    residuals[c(3L, 6L), "Q1"] <- c(1, 2)

    result <- initialise_idiosyncratic_dynamics(
        residuals,
        layout,
        min_variance = 1e-5
    )

    expect_equal(
        unname(result$monthly_rho["M1"]),
        0.5,
        tolerance = 1e-12
    )
    expect_equal(
        unname(result$monthly_variances["M1"]),
        1e-5
    )
})

test_that("quarterly idiosyncratic persistence starts at zero and variance is MM-scaled", {
    layout <- build_state_layout(
        parse_factor_spec(
            series = c("M1", "Q1"),
            factors = "global",
            quarterly_vars = "Q1"
        ),
        idio_ar1 = TRUE
    )

    residuals <- matrix(
        NA_real_,
        nrow = 12L,
        ncol = 2L,
        dimnames = list(NULL, c("M1", "Q1"))
    )

    residuals[, "M1"] <- rnorm(12L)
    residuals[c(3L, 6L, 9L, 12L), "Q1"] <- c(1, 2, 4, 7)

    result <- initialise_idiosyncratic_dynamics(
        residuals,
        layout
    )

    expected <- var(c(1, 2, 4, 7)) /
        sum(.mariano_murasawa_weights^2)

    expect_equal(
        unname(result$quarterly_rho["Q1"]),
        0
    )
    expect_equal(
        unname(result$quarterly_variances["Q1"]),
        expected
    )
})

test_that("monthly idiosyncratic persistence is clipped to the stationarity bound", {
    layout <- build_state_layout(
        parse_factor_spec(
            series = c("M1", "Q1"),
            factors = "global",
            quarterly_vars = "Q1"
        ),
        idio_ar1 = TRUE
    )

    residuals <- matrix(
        NA_real_,
        nrow = 8L,
        ncol = 2L,
        dimnames = list(NULL, c("M1", "Q1"))
    )

    residuals[, "M1"] <- 2^(0:7)
    residuals[c(3L, 6L), "Q1"] <- c(1, 2)

    result <- initialise_idiosyncratic_dynamics(
        residuals,
        layout,
        stationarity_bound = 0.8
    )

    expect_equal(
        unname(result$monthly_rho["M1"]),
        0.8
    )
})

test_that("complete block state initialisation returns the full state-space model", {
    fixture <- make_simple_block_init_fixture()

    state <- initialise_block_state(
        X = fixture$X,
        X_imp = fixture$X_imp,
        spec = fixture$spec,
        layout = fixture$layout,
        loading_maps = fixture$loading_maps
    )

    n_states <- fixture$layout$n_states
    n_series <- ncol(fixture$X)

    expect_named(
        state,
        c("A", "C", "Q", "R", "F_0", "P_0")
    )
    expect_identical(dim(state$A), c(n_states, n_states))
    expect_identical(dim(state$Q), c(n_states, n_states))
    expect_identical(dim(state$C), c(n_series, n_states))
    expect_identical(dim(state$R), c(n_series, n_series))
    expect_length(state$F_0, n_states)
    expect_identical(dim(state$P_0), c(n_states, n_states))

    expect_true(all(is.finite(state$A)))
    expect_true(all(is.finite(state$C)))
    expect_true(all(is.finite(state$Q)))
    expect_true(all(is.finite(state$R)))
    expect_true(all(is.finite(state$P_0)))
})

test_that("complete block state initialisation produces stationary initial moments", {
    fixture <- make_simple_block_init_fixture()

    state <- initialise_block_state(
        X = fixture$X,
        X_imp = fixture$X_imp,
        spec = fixture$spec,
        layout = fixture$layout,
        loading_maps = fixture$loading_maps
    )

    reconstructed <- (
        state$A %*%
            state$P_0 %*%
            t(state$A) +
            state$Q
    )

    expect_equal(
        state$P_0,
        reconstructed,
        tolerance = 1e-8
    )
    expect_equal(
        state$P_0,
        t(state$P_0),
        tolerance = 1e-10
    )
})

test_that("complete block state initialisation requires idiosyncratic AR1 layout", {
    fixture <- make_simple_block_init_fixture()
    layout <- build_state_layout(
        fixture$spec,
        idio_ar1 = FALSE
    )
    maps <- build_loading_maps(
        fixture$spec,
        layout
    )

    expect_error(
        initialise_block_state(
            X = fixture$X,
            X_imp = fixture$X_imp,
            spec = fixture$spec,
            layout = layout,
            loading_maps = maps
        ),
        "requires idio_ar1 = TRUE"
    )
})
