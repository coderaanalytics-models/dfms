test_that("block matrices have correct dimensions", {
    spec <- make_block_spec()
    params <- make_block_parameters()
    layout <- build_state_layout(spec)$factor

    ss <- build_factor_matrices(
        layout,
        params$coefficients,
        params$covariances
    )

    expect_identical(dim(ss$A), c(20L, 20L))
    expect_identical(dim(ss$Q), c(20L, 20L))
})

test_that("different factor blocks do not interact", {
    spec <- make_block_spec()
    params <- make_block_parameters()
    layout <- build_state_layout(spec)$factor

    ss <- build_factor_matrices(
        layout,
        params$coefficients,
        params$covariances
    )

    blocks <- layout$block_indices

    for (i in seq_along(blocks)) {
        for (j in seq_along(blocks)) {
            if (i == j) {
                next
            }

            expect_true(
                all(ss$A[blocks[[i]], blocks[[j]]] == 0)
            )
            expect_true(
                all(ss$Q[blocks[[i]], blocks[[j]]] == 0)
            )
        }
    }
})

test_that("lagged factor states shift correctly", {
    spec <- make_block_spec()
    params <- make_block_parameters()
    layout <- build_state_layout(spec)$factor

    ss <- build_factor_matrices(
        layout,
        params$coefficients,
        params$covariances
    )

    lag_idx <- layout$lag_indices$global

    for (j in 2:length(lag_idx)) {
        expect_equal(
            ss$A[
                lag_idx[[j]],
                lag_idx[[j - 1L]]
            ],
            diag(2)
        )
    }
})

test_that("manufacturing follows the specified AR(2)", {
    spec <- make_block_spec()
    params <- make_block_parameters()
    layout <- build_state_layout(spec)$factor

    ss <- build_factor_matrices(
        layout,
        params$coefficients,
        params$covariances
    )

    idx <- layout$lag_indices$manufacturing

    expect_equal(ss$A[idx[[1]], idx[[1]]], 0.6)
    expect_equal(ss$A[idx[[1]], idx[[2]]], -0.2)
    expect_true(
        all(ss$A[idx[[1]], unlist(idx[3:5])] == 0)
    )
})

test_that("factor innovations enter current factors only", {
    spec <- make_block_spec()
    params <- make_block_parameters()
    layout <- build_state_layout(spec)$factor

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
    full_layout <- build_state_layout(spec, idio_ar1 = TRUE)
    loading_maps <- build_loading_maps(spec, full_layout)

    C_factor <- build_factor_measurement(
        spec,
        full_layout$factor,
        loading_maps,
        make_block_loadings()
    )

    expect_identical(dim(C_factor), c(5L, 20L))
    expect_identical(
        rownames(C_factor),
        rownames(spec$loading_mask)
    )
})

test_that("monthly variables load on current factors only", {
    spec <- make_block_spec()
    full_layout <- build_state_layout(spec, idio_ar1 = TRUE)
    layout <- full_layout$factor
    loading_maps <- build_loading_maps(spec, full_layout)

    C_factor <- build_factor_measurement(
        spec,
        layout,
        loading_maps,
        make_block_loadings()
    )

    current <- unname(layout$current_indices)
    lagged <- setdiff(
        seq_len(layout$n_states),
        current
    )

    expect_equal(
        unname(C_factor["PMI", current]),
        c(0.8, 0.3, 0.0, 0.5)
    )
    expect_true(all(C_factor["PMI", lagged] == 0))
})

test_that("quarterly variables use Mariano-Murasawa weights", {
    spec <- make_block_spec()
    full_layout <- build_state_layout(spec, idio_ar1 = TRUE)
    layout <- full_layout$factor
    loading_maps <- build_loading_maps(spec, full_layout)

    C_factor <- build_factor_measurement(
        spec,
        layout,
        loading_maps,
        make_block_loadings()
    )

    global_states <- layout$block_indices$global

    expect_equal(
        unname(C_factor["GDP", global_states]),
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

    expect_true(all(C_factor["GDP", other_states] == 0))
})

test_that("restricted factor loadings remain zero", {
    spec <- make_block_spec()
    full_layout <- build_state_layout(spec, idio_ar1 = TRUE)
    layout <- full_layout$factor
    loading_maps <- build_loading_maps(spec, full_layout)

    loadings <- make_block_loadings()
    loadings["PMI", "consumption"] <- 100
    loadings["GDP", "manufacturing"] <- 100

    C_factor <- build_factor_measurement(
        spec,
        layout,
        loading_maps,
        loadings
    )

    expect_true(
        all(
            C_factor[
                "PMI",
                layout$block_indices$consumption
            ] == 0
        )
    )
    expect_true(
        all(
            C_factor[
                "GDP",
                layout$block_indices$manufacturing
            ] == 0
        )
    )
})

test_that("monthly-only measurement matrices work", {
    spec <- make_monthly_block_spec()
    full_layout <- build_state_layout(spec, idio_ar1 = TRUE)
    loading_maps <- build_loading_maps(spec, full_layout)

    C_factor <- build_factor_measurement(
        spec,
        full_layout$factor,
        loading_maps,
        make_block_loadings()
    )

    expect_identical(dim(C_factor), c(5L, 5L))

    current <- unname(full_layout$factor$current_indices)

    expect_equal(
        unname(C_factor["PMI", current]),
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
        unname(factor_ss$A)
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
        .mariano_murasawa_weights
    )
})

test_that("quarterly lag states have zero innovation variance", {
    model <- make_test_full_state()
    ss <- model$matrices
    lagged <- model$layout$quarterly_lagged

    expect_true(all(ss$Q[lagged, ] == 0))
    expect_true(all(ss$Q[, lagged] == 0))
})

test_that("idiosyncratic innovations are independent of factors", {
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

    expect_equal(
        model$matrices$R,
        diag(1e-6, 5)
    )
})

test_that("nonstationary idiosyncratic AR1 is rejected", {
    spec <- make_block_spec()
    layout <- build_state_layout(spec, idio_ar1 = TRUE)
    loading_maps <- build_loading_maps(spec, layout)
    params <- make_block_parameters()

    factor_ss <- build_factor_matrices(
        layout$factor,
        params$coefficients,
        params$covariances
    )

    factor_C <- build_factor_measurement(
        spec,
        layout$factor,
        loading_maps,
        make_block_loadings()
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

test_that("full-state initialisation has correct dimensions", {
    model <- make_test_full_state()

    init <- initialise_full_state(
        model$matrices$A,
        model$matrices$Q
    )

    n <- nrow(model$matrices$A)

    expect_length(init$F_0, n)
    expect_equal(dim(init$P_0), c(n, n))
    expect_equal(init$F_0, rep(0, n))
})

test_that("initial covariance solves the Lyapunov equation", {
    model <- make_test_full_state()

    init <- initialise_full_state(
        model$matrices$A,
        model$matrices$Q
    )

    P <- init$P_0
    reconstructed <- (
        model$matrices$A %*%
            P %*%
            t(model$matrices$A) +
            model$matrices$Q
    )

    expect_equal(
        P,
        reconstructed,
        tolerance = 1e-8
    )
})

test_that("initial covariance is symmetric", {
    model <- make_test_full_state()

    init <- initialise_full_state(
        model$matrices$A,
        model$matrices$Q
    )

    expect_equal(
        init$P_0,
        t(init$P_0),
        tolerance = 1e-10
    )
})

test_that("initial covariance is positive semidefinite", {
    model <- make_test_full_state()

    init <- initialise_full_state(
        model$matrices$A,
        model$matrices$Q
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
    model <- make_test_full_state()

    init <- initialise_full_state(
        model$matrices$A,
        model$matrices$Q
    )

    idx <- model$layout$quarterly_indices[["GDP"]]
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
