make_block_spec <- function() {
	parse_factor_spec(
			series = c("PMI", "IP", "SALES", "RETAIL", "GDP"),
			factors = list(
					PMI = c("global", "manufacturing"),
					IP = c("global", "manufacturing"),
					SALES = c("global", "consumption"),
					RETAIL = c("global", "consumption"),
					GDP = "global"
			),
			factor_multiplicities = c(
					global = 2L,
					manufacturing = 1L,
					consumption = 1L
			),
			factor_orders = c(
					global = 1L,
					manufacturing = 2L,
					consumption = 1L
			),
			quarterly_vars = "GDP"
	)
}

make_block_parameters <- function() {
	list(
			coefficients = list(
					global = array(
							c(0.7, 0.1, 0.2, 0.6),
							dim = c(2, 2, 1)
					),
					consumption = array(
							0.5,
							dim = c(1, 1, 1)
					),
					manufacturing = array(
							c(0.6, -0.2),
							dim = c(1, 1, 2)
					)
			),
			covariances = list(
					global = matrix(
							c(1.0, 0.2, 0.2, 0.8),
							nrow = 2
					),
					consumption = matrix(0.5),
					manufacturing = matrix(0.7)
			)
	)
}

make_block_loadings <- function() {
	matrix(
			c(
					# global.1, global.2, consumption, manufacturing
					0.8, 0.3, 0.0, 0.5,  # PMI
					0.7, 0.2, 0.0, 0.6,  # IP
					0.6, 0.1, 0.4, 0.0,  # SALES
					0.5, 0.2, 0.3, 0.0,  # RETAIL
					0.9, 0.4, 0.0, 0.0   # GDP
			),
			nrow = 5,
			byrow = TRUE,
			dimnames = list(
					c("PMI", "IP", "SALES", "RETAIL", "GDP"),
					c(
							"global.1",
							"global.2",
							"consumption",
							"manufacturing"
					)
			)
	)
}

