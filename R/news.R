dfm_news_kfs <- function(X, state, k) {
  A <- state$A
  C <- state$C
  Q <- state$Q
  R <- state$R
  F_0 <- state$F_0
  P_0 <- state$P_0

  if(any(vapply(list(A, C, Q, R, F_0, P_0), is.null, logical(1L)))) {
    stop("state is missing required matrices")
  }

  state_dim <- ncol(C)
  if(nrow(A) != state_dim || ncol(A) != state_dim) stop("state transition matrix has incompatible dimensions")
  if(nrow(Q) != state_dim || ncol(Q) != state_dim) stop("state covariance matrix has incompatible dimensions")
  if(length(F_0) != state_dim) stop("state initial mean has incompatible dimensions")
  if(nrow(P_0) != state_dim || ncol(P_0) != state_dim) stop("state initial covariance has incompatible dimensions")
  if(nrow(C) != ncol(X)) stop("state observation matrix has incompatible dimensions")

  if(k > 0L) {
    aug_dim <- state_dim * (k + 1L)
    A_aug <- matrix(0, aug_dim, aug_dim)
    A_aug[1:state_dim, 1:state_dim] <- A
    A_aug[(state_dim + 1L):aug_dim, 1:(k * state_dim)] <- diag(k * state_dim)

    C_aug <- cbind(C, matrix(0, nrow(C), aug_dim - state_dim))
    Q_aug <- matrix(0, aug_dim, aug_dim)
    Q_aug[1:state_dim, 1:state_dim] <- Q

    F_0 <- c(F_0, rep(0, aug_dim - state_dim))
    P_0 <- rbind(cbind(P_0, matrix(0, state_dim, aug_dim - state_dim)),
                 matrix(0, aug_dim - state_dim, aug_dim))
    diag(P_0)[(state_dim + 1L):aug_dim] <- 1e-8
  } else {
    A_aug <- A
    C_aug <- C
    Q_aug <- Q
  }
  # Stabilize Q if singular (zero variance states)
  if(any(diag(Q_aug) == 0)) {
    diag(Q_aug)[diag(Q_aug) == 0] <- 1e-8
  }

  kfs_res <- SKFS(X, A_aug, C_aug, Q_aug, R, F_0, P_0, FALSE)

  F_sm <- kfs_res$F_smooth[, 1:state_dim, drop = FALSE]
  X_sm <- tcrossprod(F_sm, C)
  list(X_sm = X_sm, P = kfs_res$P_smooth, C = C, R = R, state_dim = state_dim)
}

dfm_news_restore_missing <- function(X) {
  W <- attr(X, "missing")
  if(!is.null(W)) X[W] <- NA
  qM(X)
}

dfm_news_stats <- function(X) {
  stats <- unclass(attr(X, "stats"))
  if(is.null(stats)) {
    n <- ncol(X)
    return(list(Mx = rep(0, n), Wx = rep(1, n)))
  }
  list(Mx = stats[, "Mean"], Wx = stats[, "SD"])
}

resolve_vars <- function(target.vars, n, names) {
  if(is.null(target.vars)) return(seq_len(n))
  if(is.character(target.vars)) {
    if(is.null(names)) stop("target.vars is a name but data have no column names")
    return(as.integer(ckmatch(target.vars, names, e = "Unknown target.vars:")))
  }
  if(!is.numeric(target.vars)) stop("target.vars must be NULL, numeric indices, or names")
  target.vars <- as.integer(target.vars)
  if(any(target.vars < 1L | target.vars > n)) stop("target.vars is out of bounds")
  unique(target.vars)
}
