#' Cross-Fitted Nuisance Estimation
#'
#' Estimates all nuisance functions required for DML-Wald and DML-TC using
#' K-fold cross-fitting. Supports Lasso, Random Forest, neural networks,
#' OLS, or a user-supplied function.
#'
#' @param Y Numeric vector of outcomes.
#' @param D Numeric vector of treatment.
#' @param G Binary vector of group assignment (0/1).
#' @param Ti Binary vector of time period (0/1).
#' @param X Numeric matrix of covariates.
#' @param method Character string (`"lasso"`, `"rf"`, `"nn"`, `"ols"`) or a
#'   function `f(Y_train, X_train, X_predict)` returning predicted values.
#' @param K Integer, number of cross-fitting folds (default 5).
#' @param estimand Which estimand to prepare nuisance for: `"wald"`, `"tc"`,
#'   or `"both"` (default).
#' @param seed Random seed for fold assignment (default `NULL`).
#'
#' @return A list with cross-fitted predictions for every observation:
#'   \describe{
#'     \item{m_Y_10, m_Y_01, m_Y_00}{Conditional outcome means E\[Y|G=g,T=t,X\].}
#'     \item{m_D_10, m_D_01, m_D_00}{Conditional treatment means E\[D|G=g,T=t,X\].}
#'     \item{pG_raw}{Propensity score Pr(G=1|X).}
#'     \item{pi_11, pi_10, pi_01, pi_00}{Joint cell probabilities.}
#'     \item{mu_Y_101, ..., mu_Y_000}{Treatment-conditional means (TC only).}
#'     \item{pi_101, ..., pi_000}{Treatment-conditional propensities (TC only).}
#'     \item{folds}{Integer vector of fold assignments.}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' N <- 500
#' X <- matrix(rnorm(N * 5), ncol = 5)
#' G <- rbinom(N, 1, 0.5)
#' Ti <- rbinom(N, 1, 0.5)
#' D <- rbinom(N, 1, plogis(X[,1] + G))
#' Y <- 1 + X[,1] + D + rnorm(N)
#' nuis <- fddml_nuisance(Y, D, G, Ti, X, method = "ols", K = 2)
#' }
#'
#' @export
fddml_nuisance <- function(Y, D, G, Ti, X,
                            method = "lasso",
                            K = 5L,
                            estimand = "both",
                            seed = NULL) {
  X <- as.matrix(X)
  N <- length(Y)
  .validate_inputs(Y, D, G, Ti, X)

  folds <- .make_folds(N, K, seed = seed)
  pT <- mean(Ti)

  # Cell indices
  cells <- list(
    "10" = which(G == 1 & Ti == 0),
    "01" = which(G == 0 & Ti == 1),
    "00" = which(G == 0 & Ti == 0)
  )

  # Initialize prediction vectors
  m_Y <- m_D <- list()
  for (key in names(cells)) {
    m_Y[[key]] <- rep(NA_real_, N)
    m_D[[key]] <- rep(NA_real_, N)
  }

  # Cross-fitted conditional expectations within each (g,t) cell
  for (key in names(cells)) {
    idx <- cells[[key]]
    for (k in seq_len(K)) {
      train_k <- idx[folds[idx] != k]
      pred_k <- which(folds == k)

      if (length(train_k) < 5) {
        m_Y[[key]][pred_k] <- mean(Y[idx])
        m_D[[key]][pred_k] <- mean(D[idx])
        next
      }

      m_Y[[key]][pred_k] <- .fit_nuisance_cell(Y[train_k], X[train_k, , drop = FALSE],
                                                 X[pred_k, , drop = FALSE], method)
      m_D[[key]][pred_k] <- .fit_nuisance_cell(D[train_k], X[train_k, , drop = FALSE],
                                                 X[pred_k, , drop = FALSE], method)
    }
  }

  # Cross-fitted propensity score P(G=1|X)
  pG_raw <- rep(NA_real_, N)
  for (k in seq_len(K)) {
    train_k <- which(folds != k)
    pred_k <- which(folds == k)
    pG_raw[pred_k] <- .fit_nuisance_cell(G[train_k], X[train_k, , drop = FALSE],
                                           X[pred_k, , drop = FALSE], method,
                                           family = "binomial")
  }
  pG_raw <- pmax(pmin(pG_raw, 0.999), 0.001)

  # Joint cell probabilities
  pi_11 <- pG_raw * pT
  pi_10 <- pG_raw * (1 - pT)
  pi_01 <- (1 - pG_raw) * pT
  pi_00 <- (1 - pG_raw) * (1 - pT)

  result <- list(
    m_Y_10 = m_Y[["10"]], m_Y_01 = m_Y[["01"]], m_Y_00 = m_Y[["00"]],
    m_D_10 = m_D[["10"]], m_D_01 = m_D[["01"]], m_D_00 = m_D[["00"]],
    pG_raw = pG_raw,
    pi_11 = pi_11, pi_10 = pi_10, pi_01 = pi_01, pi_00 = pi_00,
    folds = folds
  )

  # TC-specific: treatment-conditional expectations in control group
  if (estimand %in% c("tc", "both")) {
    tc_cells <- list(
      "101" = which(D == 1 & G == 0 & Ti == 1),
      "100" = which(D == 1 & G == 0 & Ti == 0),
      "001" = which(D == 0 & G == 0 & Ti == 1),
      "000" = which(D == 0 & G == 0 & Ti == 0)
    )

    for (key in names(tc_cells)) {
      idx <- tc_cells[[key]]
      mu_key <- paste0("mu_Y_", key)
      pi_key <- paste0("pi_", key)
      result[[mu_key]] <- rep(NA_real_, N)
      result[[pi_key]] <- rep(NA_real_, N)

      for (k in seq_len(K)) {
        train_k <- idx[folds[idx] != k]
        pred_k <- which(folds == k)

        if (length(train_k) < 5) {
          result[[mu_key]][pred_k] <- if (length(idx) > 0) mean(Y[idx]) else 0
          result[[pi_key]][pred_k] <- max(length(idx) / N, 0.001)
          next
        }

        result[[mu_key]][pred_k] <- .fit_nuisance_cell(
          Y[train_k], X[train_k, , drop = FALSE],
          X[pred_k, , drop = FALSE], method)

        # Propensity for this (d,g,t) cell
        L_ind <- as.numeric(seq_len(N) %in% idx)
        result[[pi_key]][pred_k] <- .fit_nuisance_cell(
          L_ind[folds != k], X[folds != k, , drop = FALSE],
          X[pred_k, , drop = FALSE], method, family = "binomial")
      }
      result[[pi_key]] <- pmax(result[[pi_key]], 0.001)
    }
  }

  result
}
