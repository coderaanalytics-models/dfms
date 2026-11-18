test_that("global factors load on all series", {
    spec <- parse_factor_spec(
        series = c("PMI", "IP", "GDP"),
        factors = "global",
        quarterly_vars = "GDP"
    )

    expect_true(all(spec$loading_mask))
    expect_equal(spec$n_factors, 1L)
})

test_that("sector restrictions produce structural zeros", {
    spec <- parse_factor_spec(
        series = c("PMI", "IP", "SALES", "GDP"),
        factors = list(
            PMI = c("global", "manufacturing"),
            IP = c("global", "manufacturing"),
            SALES = c("global", "consumption"),
            GDP = "global"
        ),
        quarterly_vars = "GDP"
    )

    expect_true(spec$loading_mask["PMI", "manufacturing"])
    expect_false(spec$loading_mask["SALES", "manufacturing"])
    expect_false(spec$loading_mask["GDP", "consumption"])
})

test_that("factor multiplicities expand correctly", {
    spec <- parse_factor_spec(
        series = c("PMI", "IP", "GDP"),
        factors = "global",
        factor_multiplicities = 2L,
        quarterly_vars = "GDP"
    )

    expect_identical(
        spec$factor_names,
        c("global.1", "global.2")
    )
    expect_equal(spec$n_factors, 2L)
})

test_that("factor blocks are disjoint", {
    spec <- parse_factor_spec(
        series = c("PMI", "IP", "SALES", "GDP"),
        factors = list(
            PMI = c("global", "manufacturing"),
            IP = c("global", "manufacturing"),
            SALES = c("global", "consumption"),
            GDP = "global"
        ),
        factor_orders = list(
            global = list(
                factors = "global",
                order = 1L
            ),
            sectors = list(
                factors = c("manufacturing", "consumption"),
                order = 2L
            )
        ),
        quarterly_vars = "GDP"
    )

    expect_equal(spec$n_blocks, 2L)
    expect_setequal(
        unlist(spec$blocks, use.names = FALSE),
        spec$factor_names
    )
    expect_identical(
        spec$block_orders,
        c(global = 1L, sectors = 2L)
    )
})

test_that("loading mask follows observed series ordering", {
    series <- c("PMI", "IP", "SALES", "GDP")

    spec <- parse_factor_spec(
        series = series,
        factors = list(
            GDP = "global",
            SALES = c("global", "consumption"),
            PMI = c("global", "manufacturing"),
            IP = c("global", "manufacturing")
        ),
        quarterly_vars = "GDP"
    )

    expect_identical(rownames(spec$loading_mask), series)
})

test_that("monthly and quarterly variables are retained in the specification", {
    spec <- make_block_spec()

    expect_identical(spec$quarterly_vars, "GDP")
    expect_identical(
        spec$monthly_vars,
        c("PMI", "IP", "SALES", "RETAIL")
    )
})

test_that("invalid specifications are rejected", {
    expect_error(
        parse_factor_spec(
            series = c("PMI", "PMI"),
            factors = "global"
        ),
        "unique"
    )

    expect_error(
        parse_factor_spec(
            series = c("PMI", "IP"),
            factors = "global",
            factor_multiplicities = 3L
        ),
        "monthly"
    )

    expect_error(
        parse_factor_spec(
            series = c("PMI", "IP"),
            factors = "global",
            quarterly_vars = "GDP"
        ),
        "subset"
    )
})

test_that("named factor orders are parsed correctly", {
    spec <- parse_factor_spec(
        series = c("PMI", "IP", "SALES", "GDP"),
        factors = list(
            PMI = c("global", "manufacturing"),
            IP = c("global", "manufacturing"),
            SALES = c("global", "consumption"),
            GDP = "global"
        ),
        factor_orders = c(
            global = 1L,
            manufacturing = 2L,
            consumption = 1L
        ),
        quarterly_vars = "GDP"
    )

    expect_identical(
        unname(spec$block_orders["manufacturing"]),
        2L
    )
    expect_equal(spec$n_blocks, 3L)
})

test_that("monthly factor state dimensions follow block VAR orders", {
    spec <- make_monthly_block_spec()
    layout <- build_state_layout(spec)$factor

    expect_equal(layout$n_states, 5L)
    expect_identical(
        layout$block_lags,
        c(
            global = 1L,
            consumption = 1L,
            manufacturing = 2L
        )
    )
    expect_identical(
        layout$current_indices,
        c(
            global.1 = 1L,
            global.2 = 2L,
            consumption = 3L,
            manufacturing = 4L
        )
    )
})

test_that("quarterly specification retains five factor lags", {
    spec <- make_block_spec()
    layout <- build_state_layout(spec)$factor

    expect_equal(layout$n_states, 20L)
    expect_true(all(layout$block_lags >= 5L))
    expect_identical(
        layout$current_indices,
        c(
            global.1 = 1L,
            global.2 = 2L,
            consumption = 11L,
            manufacturing = 16L
        )
    )
})

test_that("factor lag indices cover the factor state without overlap", {
    layout <- build_state_layout(make_block_spec())$factor

    all_indices <- unlist(
        layout$lag_indices,
        use.names = FALSE
    )

    expect_setequal(all_indices, seq_len(layout$n_states))
    expect_equal(length(unique(all_indices)), layout$n_states)
})

test_that("factor lag indices map individual factors across retained lags", {
    layout <- build_state_layout(make_block_spec())$factor

    expect_identical(
        unname(layout$factor_lag_indices$global.1),
        c(1L, 3L, 5L, 7L, 9L)
    )
    expect_identical(
        unname(layout$factor_lag_indices$manufacturing),
        16:20
    )
})

test_that("factor blocks occupy contiguous state regions", {
    layout <- build_state_layout(make_block_spec())$factor

    expect_identical(layout$block_indices$global, 1:10)
    expect_identical(layout$block_indices$consumption, 11:15)
    expect_identical(layout$block_indices$manufacturing, 16:20)
})

test_that("higher VAR orders are preserved", {
    spec <- parse_factor_spec(
        series = paste0("X", 1:8),
        factors = "global",
        factor_orders = 6L,
        quarterly_vars = "X8"
    )

    layout <- build_state_layout(spec)$factor

    expect_equal(layout$n_states, 6L)
    expect_identical(layout$block_lags, c(global = 6L))
})

test_that("VAR orders are distinct from retained quarterly state lags", {
    layout <- build_state_layout(make_block_spec())$factor

    expect_equal(
        layout$block_orders,
        c(
            global = 1L,
            consumption = 1L,
            manufacturing = 2L
        )
    )
    expect_equal(
        layout$block_lags,
        c(
            global = 5L,
            consumption = 5L,
            manufacturing = 5L
        )
    )
})

test_that("full mixed-frequency state layout has correct dimensions", {
    layout <- build_state_layout(make_block_spec())

    expect_equal(layout$n_states, 29L)
    expect_equal(layout$n_factor_states, 20L)
    expect_equal(layout$n_monthly, 4L)
    expect_equal(layout$n_quarterly, 1L)
})

test_that("factor indices occupy the leading region of the complete state", {
    layout <- build_state_layout(make_block_spec())

    expect_identical(
        layout$factor_indices,
        seq_len(layout$factor$n_states)
    )
    expect_identical(
        unlist(layout$factor$block_indices, use.names = FALSE),
        layout$factor_indices
    )
    expect_identical(layout$n_factor_states, layout$factor$n_states)
})

test_that("monthly idiosyncratic states are indexed correctly", {
    layout <- build_state_layout(make_block_spec())

    expect_equal(
        unname(layout$monthly_indices),
        21:24
    )
    expect_identical(
        names(layout$monthly_indices),
        c("PMI", "IP", "SALES", "RETAIL")
    )
})

test_that("quarterly idiosyncratic states have five lags", {
    layout <- build_state_layout(make_block_spec())

    expect_equal(
        unname(layout$quarterly_indices$GDP),
        25:29
    )
    expect_identical(
        names(layout$quarterly_indices$GDP),
        paste0("lag", 0:4)
    )
    expect_equal(unname(layout$quarterly_current), 25L)
    expect_equal(layout$quarterly_lagged, 26:29)
})

test_that("all complete-state indices are assigned exactly once", {
    layout <- build_state_layout(make_block_spec())

    indices <- c(
        layout$factor_indices,
        unname(layout$monthly_indices),
        unlist(
            layout$quarterly_indices,
            use.names = FALSE
        )
    )

    expect_equal(sort(indices), seq_len(layout$n_states))
    expect_equal(anyDuplicated(indices), 0L)
})

test_that("quarterly status in the specification determines state lag retention", {
    monthly <- build_state_layout(make_monthly_block_spec())
    mixed <- build_state_layout(make_block_spec())

    expect_false(monthly$factor$quarterly)
    expect_true(mixed$factor$quarterly)
    expect_identical(
        monthly$factor$block_orders,
        mixed$factor$block_orders
    )
    expect_identical(
        monthly$factor$block_lags,
        c(global = 1L, consumption = 1L, manufacturing = 2L)
    )
    expect_identical(
        mixed$factor$block_lags,
        c(global = 5L, consumption = 5L, manufacturing = 5L)
    )
    expect_identical(
        mixed$idio_current,
        c(mixed$monthly_indices, GDP = 25L)
    )
})

test_that("monthly-only state layouts work", {
    layout <- build_state_layout(make_monthly_block_spec())

    expect_equal(layout$n_factor_states, 5L)
    expect_equal(layout$n_monthly, 5L)
    expect_equal(layout$n_quarterly, 0L)
    expect_equal(layout$n_states, 10L)

    expect_length(layout$quarterly_indices, 0L)
    expect_length(layout$quarterly_current, 0L)
    expect_length(layout$quarterly_lagged, 0L)
})

test_that("monthly loading maps select contemporaneous factors", {
    spec <- make_block_spec()
    layout <- build_state_layout(spec, idio_ar1 = TRUE)
    maps <- build_loading_maps(spec, layout)

    for (series_name in layout$monthly_vars) {
        map <- maps[[series_name]]

        allowed <- colnames(spec$loading_mask)[
            spec$loading_mask[series_name, ]
        ]

        expect_equal(rownames(map$indices), allowed)
        expect_equal(ncol(map$indices), 1L)
        expect_equal(
            as.integer(map$indices[, 1L]),
            unname(layout$factor$current_indices[allowed])
        )
        expect_equal(map$weights, 1)
        expect_equal(
            map$idio_indices,
            unname(layout$monthly_indices[[series_name]])
        )
    }
})

test_that("quarterly loading maps use Mariano-Murasawa weights", {
    spec <- make_block_spec()
    layout <- build_state_layout(spec, idio_ar1 = TRUE)
    maps <- build_loading_maps(spec, layout)
    map <- maps$GDP

    expect_equal(map$weights, .mariano_murasawa_weights)
    expect_equal(ncol(map$indices), 5L)
    expect_equal(
        rownames(map$indices),
        colnames(spec$loading_mask)[
            spec$loading_mask["GDP", ]
        ]
    )
})

test_that("loading maps preserve block-major factor ordering", {
    spec <- make_block_spec()
    layout <- build_state_layout(spec, idio_ar1 = TRUE)
    maps <- build_loading_maps(spec, layout)

    for (series_name in names(maps)) {
        map <- maps[[series_name]]

        for (factor in rownames(map$indices)) {
            block <- names(spec$blocks)[
                vapply(
                    spec$blocks,
                    function(x) factor %in% x,
                    logical(1)
                )
            ]

            position <- match(
                factor,
                spec$blocks[[block]]
            )

            for (lag in seq_len(ncol(map$indices))) {
                expected <- layout$factor$lag_indices[[block]][[lag]][[position]]

                expect_equal(
                    unname(map$indices[factor, lag]),
                    unname(expected)
                )
            }
        }
    }
})
