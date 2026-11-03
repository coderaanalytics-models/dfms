#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List update_loadings_cpp(
	const arma::mat& X,
	const arma::mat& F_smooth,
	const arma::cube& P_smooth,
	const List& loading_maps
) {
	const arma::uword T = X.n_rows;
	const arma::uword n_states = F_smooth.n_cols;

	if (F_smooth.n_rows != T ||
		P_smooth.n_rows != n_states ||
		P_smooth.n_cols != n_states ||
		P_smooth.n_slices != T) {
			Rcpp::stop("Incompatible smoother dimensions.");
		}

	if (loading_maps.size() != X.n_cols) {
		Rcpp::stop("Number of loading maps must match X columns.");
	}

	List results(X.n_cols);

	for (arma::uword i = 0; i < X.n_cols; ++i) {
		List map = loading_maps[i];

		IntegerMatrix indices = map["indices"];
		NumericVector weights = map["weights"];
		int idio = as<int>(map["idio_index"]) - 1;

		const arma::uword k = indices.nrow();
		const arma::uword L = indices.ncol();

		if (weights.size() != static_cast<int>(L)) {
			Rcpp::stop("Incompatible loading-map dimensions.");
		}

		if (idio < 0 || static_cast<arma::uword>(idio) >= n_states) {
			Rcpp::stop("Invalid idiosyncratic state index.");
		}

		arma::mat S_gg(k, k, arma::fill::zeros);
		arma::vec S_gy(k, arma::fill::zeros);
		arma::vec S_gu(k, arma::fill::zeros);

		arma::vec g(k);

		for (arma::uword t = 0; t < T; ++t) {
			const double y = X(t, i);

			if (!std::isfinite(y)) {
				continue;
			}

			g.zeros();

			for (arma::uword j = 0; j < k; ++j) {
				for (arma::uword lag = 0; lag < L; ++lag) {
					const int idx = indices(j, lag) - 1;

					g(j) += weights[lag] * F_smooth(t, idx);
				}
			}

			const double u = F_smooth(t, idio);

			S_gg += g * g.t();
			S_gy += g * y;
			S_gu += g * u;

			const arma::mat& P = P_smooth.slice(t);

			for (arma::uword j = 0; j < k; ++j) {
				for (arma::uword a = 0; a < L; ++a) {
					const int idx_j = indices(j, a) - 1;
					const double w_a = weights[a];

					S_gu(j) += w_a * P(idx_j, idio);

					for (arma::uword h = 0; h < k; ++h) {
						for (arma::uword b = 0; b < L; ++b) {
							const int idx_h = indices(h, b) - 1;

							S_gg(j, h) += w_a * weights[b] * P(idx_j, idx_h);
						}
					}
				}
			}
		}

		if (k == 0) {
			results[i] = NumericVector(0);
			continue;
		}

		arma::vec lambda;

		bool success = arma::solve(
			lambda,
			S_gg,
			S_gy - S_gu,
			arma::solve_opts::no_approx
		);

		if (!success) {
			Rcpp::stop("Failed to solve loading normal equations.");
		}

		results[i] = Rcpp::NumericVector(lambda.begin(), lambda.end());
	}
	results.attr("names") = loading_maps.names();
	return(results);
}



// [[Rcpp::export]]
Rcpp::List update_factor_blocks_cpp(
	const arma::mat& F_smooth,
	const arma::cube& P_smooth,
	const arma::cube& PPm_smooth,
	const arma::vec& F_smooth_0,
	const arma::mat& P_smooth_0,
	const Rcpp::List& factor_layout
) {
	const arma::uword T = F_smooth.n_rows;
	const arma::uword n = F_smooth.n_cols;

	if (T < 2) {
		Rcpp::stop("At least two observations are required.");
	}

	if (
		P_smooth.n_rows != n ||
		P_smooth.n_cols != n ||
		P_smooth.n_slices != T ||
		PPm_smooth.n_rows != n ||
		PPm_smooth.n_cols != n ||
		PPm_smooth.n_slices != T ||
		F_smooth_0.n_elem != n ||
		P_smooth_0.n_rows != n ||
		P_smooth_0.n_cols != n
	) {
		Rcpp::stop("Incompatible smoother dimensions.");
	}

	const Rcpp::List lag_indices = factor_layout["lag_indices"];
	const Rcpp::List block_orders = factor_layout["block_orders"];
	const Rcpp::CharacterVector block_names = lag_indices.names();

	Rcpp::List coefficients(lag_indices.size());
	Rcpp::List covariances(lag_indices.size());

	for (R_xlen_t b = 0; b < lag_indices.size(); ++b) {
		const std::string name = Rcpp::as<std::string>(block_names[b]);

		const Rcpp::List block_lags = lag_indices[b];
		const int p = Rcpp::as<int>(block_orders[name]);

		if (p < 1 || block_lags.size() < p) {
			Rcpp::stop("Invalid VAR order for block: " + name);
		}

		const Rcpp::IntegerVector current_r = block_lags[0];

		const arma::uword r = current_r.size();
		const arma::uword k = r * static_cast<arma::uword>(p);

		if (r == 0) {
			Rcpp::stop("Empty factor block: " + name);
		}

		arma::uvec current(r);
		arma::uvec predictors(k);

		for (arma::uword j = 0; j < r; ++j) {
			const int idx = current_r[j];

			if (idx < 1 || static_cast<arma::uword>(idx) > n) {
				Rcpp::stop("Invalid current-state index for block: " + name);
			}

			current[j] = static_cast<arma::uword>(idx - 1);
		}

		for (int lag = 0; lag < p; ++lag) {
			const Rcpp::IntegerVector lag_r = block_lags[lag];

			if (static_cast<arma::uword>(lag_r.size()) != r) {
				Rcpp::stop("Inconsistent factor dimensions for block: " + name);
			}

			for (arma::uword j = 0; j < r; ++j) {
				const int idx = lag_r[j];

				if (idx < 1 || static_cast<arma::uword>(idx) > n) {
					Rcpp::stop("Invalid predictor index for block: " + name);
				}

				predictors[static_cast<arma::uword>(lag) * r + j] =
					static_cast<arma::uword>(idx - 1);
			}
		}

		// Sufficient statistics accumulation
		arma::mat S_yy(r, r, arma::fill::zeros);
		arma::mat S_xx(k, k, arma::fill::zeros);
		arma::mat S_yx(r, k, arma::fill::zeros);

		for (arma::uword t = 0; t < T; ++t) {
			const arma::vec f_t = F_smooth.row(t).t();

			const arma::vec f_prev = (t == 0) ? F_smooth_0: arma::vec(F_smooth.row(t - 1).t());

			const arma::vec y = f_t.elem(current);
			const arma::vec x = f_prev.elem(predictors);

			const arma::mat P_prev = (t == 0) ? P_smooth_0: arma::mat(P_smooth.slice(t - 1));

			S_yy += y * y.t() + P_smooth.slice(t).submat(current, current);
			S_xx += x * x.t() + P_prev.submat(predictors, predictors);
			S_yx += y * x.t() + PPm_smooth.slice(t).submat(current, predictors);
		}

		// Estimate the VAR coefficient and innovation covariance matrices
		arma::mat Phi_t;

		bool success = arma::solve(
			Phi_t,
			S_xx,
			S_yx.t(),
			arma::solve_opts::no_approx
		);

		if (!success) {
			Rcpp::stop("Failed to estimate VAR coefficients for block: " + name);
		}

		const arma::mat Phi = Phi_t.t();

		arma::mat Q = (S_yy - Phi * S_yx.t()) / static_cast<double>(T);

		Q = 0.5 * (Q + Q.t());

		arma::cube Phi_blocks(r, r, p, arma::fill::zeros);

		for (int lag = 0; lag < p; ++lag) {
			const arma::uword start = static_cast<arma::uword>(lag) * r;
			const arma::uword end = start + r - 1;

			Phi_blocks.slice(lag) = Phi.cols(start, end);
		}

		coefficients[b] = Rcpp::wrap(Phi_blocks);
		covariances[b] = Rcpp::wrap(Q);
	}
	
	coefficients.attr("names") = block_names;
	covariances.attr("names") = block_names;
	
	return Rcpp::List::create(
		Rcpp::Named("coefficients") = coefficients,
		Rcpp::Named("covariances") = covariances
	);
}