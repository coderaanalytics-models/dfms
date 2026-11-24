

dfm_stationary_cov <- function(A, Q) {
	rp <- nrow(A)
	M <- diag(rp^2) - kronecker(A, A)
	P_0 <- ainv(M) %*% unattrib(Q)
	dim(P_0) <- c(rp, rp)
	P_0
}


dfm_state <- function(dfm, require_full = FALSE) {
	ss_full <- dfm[["ss_full"]]

	if (!is.null(ss_full)) {
		needed = c("A", "C", "Q", "R", "F_0", "P_0")

		if (any(!needed %in% names(ss_full))) {
			stop("ss_full is missing required state matrices")
		}
		return(ss_full)
	}

	if (isTRUE(require_full)) {
		stop(
			"Full state matrices are required for MQ or idio.ar1 models: ",
			"re-estimate with the updated package"
		)
	}

	# Build full state from compact matrices for legacy baseline objects
	A <- dfm$A
	C <- dfm$C
	Q <- dfm$Q
	R <- dfm$R

	r <- nrow(A)

	if (ncol(A) %% r != 0L) {
		stop("Invalid transition matrix dimensions in dfm object")
	}

	p <- ncol(A) / r
	rp <- r * p

	A_state <- matrix(0, rp, rp)
	A_state[1:r, 1:rp] <- A

	if (rp > r) {
		A_state[(r + 1L):rp, 1:(rp - r)] <- diag(rp - r)
	}

	C_state <- cbind(
		C,
		matrix(0, nrow(C), rp - r)
	)

	Q_state <- matrix(0, rp, rp)
	Q_state[1:r, 1:r] <- Q

	F_0 <- rep(0, rp)
	P_0 <- dfm_stationary_cov(A_state, Q_state)

	list(
		A = A_state,
		C = C_state,
		Q = Q_state,
		R = R,
		F_0 = F_0,
		P_0 = P_0
	)
}


dfm_scale <- function(X, stats) {
	if (is.null(stats)) {
		stop("stats are required to scale X")
	}

	class(stats) <- NULL

	Mx <- stats[, "Mean"]
	Wx <- stats[, "SD"]

	TRA.matrix(
		TRA.matrix(X, Mx, "-"),
		Wx,
		"/"
	)
}


dfm_prepare_data <- function(X, stats = NULL, ...) {
	ax <- attributes(X)
	ilX <- is.list(X)

	if (is.null(stats)) {
		stats <- qsu(X)
		X <- fscale(qM(X), na.rm = TRUE)
	} else {
		X <- dfm_scale(qM(X), stats)
	}

	Xnam <- colnames(X)
	dimnames(X) <- NULL

	anymiss <- anyNA(X)

	if (anymiss) {
		X_imp <- tsnarmimp(X, ...)
		W <- attr(X_imp, "missing")
		rm.rows <- attr(X_imp, "rm.rows")
		attributes(X_imp) <- list(dim = dim(X_imp))

		if (length(rm.rows)) {
			X <- X[-rm.rows, , drop = FALSE]
		}
	} else {
		X_imp <- X
		W <- NULL
		rm.rows <- NULL
	}

	if (anymiss && length(rm.rows)) {
		anymiss <- any(W)
	}

	list(
		X = X,
		X_imp = X_imp,
		stats = stats,
		series = Xnam,
		W = W, # Returing for now due to downstream BM code path expecting it
		missing = if (anymiss) W else NULL,
		rm.rows = rm.rows,
		attributes = ax,
		is.list = ilX,
		anyNA = anymiss
	)
}


dfm_apply_state <- function(object, X, ...) {
	stats <- attr(object$X_imp, "stats")

	if (is.null(stats)) {
		stop("dfm object does not contain scaling statistics")
	}

	X_check <- qM(X)

	# New data must describe the same observed variables
	old_series <- colnames(object$X_imp)
	new_series <- colnames(X_check)
	
	if (ncol(X_check) != ncol(object$X_imp)) {
		stop("New data has a different number of variables")
	}

	if (
		!is.null(old_series) &&
		!is.null(new_series) &&
		!identical(old_series, new_series)
	) {
		stop("New data has different variable names or ordering")
	}

	prep <- dfm_prepare_data(X, stats = stats, ...)

	require_full <- (
		length(object$quarterly.vars) > 0L ||
		!is.null(object$rho) ||
		!is.null(object$factor.spec)
	)

	state <- dfm_state(object, require_full = require_full)

	kfs <- SKFS(
		X = prep$X,
		A = state$A,
		C = state$C,
		Q = state$Q,
		R = state$R,
		F_0 = state$F_0,
		P_0 = state$P_0,
		loglik = FALSE
	)

	list(
		prep = prep,
		state = state,
		kfs = kfs
	)
}


dfm_update_object <- function(object, applied) {
	prep <- applied$prep
	kfs <- applied$kfs

	r <- nrow(object$A)
	p <- ncol(object$A) / r

	# Identify the contemporaneous factor states
	if (!is.null(object$factor.spec)) {
		factor_names <- object$factor.spec$factor_names

		factor_indices <- unname(
			object$state.layout$factor$current_indices[factor_names]
		)
	} else {
		factor_names <- paste0("f", seq_len(r))
		factor_indices <- seq_len(r)
	}

	F_new <- kfs$F_smooth[
		,
		factor_indices,
		drop = FALSE
	]
	colnames(F_new) <- factor_names

	P_new <- kfs$P_smooth[
		factor_indices,
		factor_indices,
		,
		drop = FALSE
	]
	dimnames(P_new) <- list(
		factor_names,
		factor_names,
		NULL
	)

	# Start from the fitted model so all estimated parameters and specification
	# metadata are retained
	updated <- object

	# Replace the data attached to the model
	updated$X_imp <- structure(
		prep$X_imp,
		dimnames = list(NULL, prep$series),
		stats = prep$stats,
		missing = prep$missing,
		attributes = prep$attributes,
		is.list = prep$is.list
	)

	updated$anyNA <- prep$anyNA
	updated$rm.rows <- prep$rm.rows

	# Replace the active state estimates, but do not rerun the estimation stages that
	# produced PCA/two-step estimates.
	if (identical(object$em.method, "none")) {
		updated$F_2s <- F_new
		updated$P_2s <- P_new
	} else {
		updated$F_qml <- F_new
		updated$P_qml <- P_new
	}

	# If the original object retained its full state representaiton, preserve the
	# fitted matrices and replace only the smoothing results.
	if (!is.null(object$ss_full)) {
		updated$ss_full <- applied$state
		updated$ss_full$F_smooth <- kfs$F_smooth
		updated$ss_full$P_smooth <- kfs$P_smooth
	}

	# Update idiosyncratic state estimates where they exist
	if (!is.null(object$rho)) {
		if (!is.null(object$factor.spec)) {
			idio_indices <- unname(
				object$state.layout$idio_current[prep$series]
			)
		} else if (length(object$quarterly.vars)) {
			n <- ncol(prep$X)
			nq <- length(object$quarterly.vars)
			nm <- n - nq
			rpC <- r * max(p, 5L)

			idio_indices <- c(
				(rpC + 1L):(rpC + nm),
				rpC + nm + (seq_len(nq) - 1L) * 5L + 1L
			)
		} else {
			rp <- r * p
			n <- ncol(prep$X)
			idio_indices <- rp + seq_len(n)
		}

		updated$e <- kfs$F_smooth[
			,
			idio_indices,
			drop = FALSE
		]
		colnames(updated$e) <- prep$series
	}
	
updated
}
