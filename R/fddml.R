#' Fuzzy Difference-in-Differences with Machine Learning
#'
#' Main estimation function. Estimates the LATE for fuzzy DID designs using
#' doubly robust GMM with cross-fitted ML nuisance estimation. Implements
#' both the DML-Wald and DML-TC estimators of Mena (2026).
#'
#' @param Y Numeric vector of outcomes.
#' @param D Numeric vector of treatment (binary or continuous).
#' @param G Binary vector of group assignment (0/1).
#' @param Ti Binary vector of time period (0/1).
#' @param X Numeric matrix or data.frame of covariates.
#' @param estimand Character: `"wald"`, `"tc"`, or `"both"` (default).
#' @param method Character (`"lasso"`, `"rf"`, `"nn"`, `"ols"`) or a custom
#'   function `f(Y_train, X_train, X_predict)`.
#' @param K Integer, number of cross-fitting folds (default 5).
#' @param trim Character: `"auto"` (data-driven, default), `"fixed"`,
#'   or `"none"`.
#' @param trim_alpha Numeric, fixed trimming level when `trim = "fixed"`.
#' @param cluster Optional vector of cluster identifiers for cluster-robust SEs.
#' @param se_type Character: `"analytical"` (default), `"cluster"`, or
#'   `"bootstrap"`.
#' @param B Integer, bootstrap replications (default 1000).
#' @param seed Random seed for cross-fitting folds (default `NULL`).
#' @param verbose Logical, print progress messages (default `FALSE`).
#'
#' @return An object of class `"fddml"` containing:
#'   \describe{
#'     \item{estimates}{Data frame with columns: estimand, estimate, se,
#'       ci_lower, ci_upper, p_value.}
#'     \item{wald}{Full DML-Wald results (if requested).}
#'     \item{tc}{Full DML-TC results (if requested).}
#'     \item{nuisance}{Cross-fitted nuisance predictions.}
#'     \item{trim_info}{Trimming diagnostics.}
#'     \item{call}{The matched call.}
#'     \item{N}{Sample size.}
#'     \item{p}{Number of covariates.}
#'     \item{settings}{List of estimation settings.}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' N <- 500
#' X <- matrix(rnorm(N * 5), ncol = 5)
#' G <- rbinom(N, 1, plogis(0.5 * X[,1]))
#' Ti <- rbinom(N, 1, 0.5)
#' D <- rbinom(N, 1, plogis(X[,1] + G + Ti * G))
#' Y <- 1 + X[,1] + 0.3 * G * Ti + D * (G * Ti) + rnorm(N)
#' fit <- fddml(Y, D, G, Ti, X, estimand = "wald", method = "ols", K = 2)
#' print(fit)
#' }
#'
#' @references
#' Mena, A. (2026). Double Debiased Machine Learning for Difference-in-Differences
#' under Imperfect Compliance. Working Paper, Brown University.
#'
#' De Chaisemartin, C. and D'Haultfoeuille, X. (2018). Fuzzy
#' Differences-in-Differences. *Review of Economic Studies*, 85(2), 999-1028.
#'
#' @export
fddml <- function(Y, D, G, Ti, X,
                   estimand = "both",
                   method = "lasso",
                   K = 5L,
                   trim = "auto",
                   trim_alpha = 0.10,
                   cluster = NULL,
                   se_type = "analytical",
                   B = 1000L,
                   seed = NULL,
                   verbose = FALSE) {
  cl <- match.call()
  X <- as.matrix(X)
  N <- length(Y)
  p <- ncol(X)
  estimand <- match.arg(estimand, c("wald", "tc", "both"))

  .validate_inputs(Y, D, G, Ti, X)

  # Step 1: Cross-fitted nuisance estimation
  if (verbose) message("Step 1/4: Nuisance estimation (K=", K, ", method=", method, ")...")
  nuis <- fddml_nuisance(Y, D, G, Ti, X, method = method, K = K,
                          estimand = estimand, seed = seed)

  # Step 2: Trimming
  if (verbose) message("Step 2/4: Propensity trimming (", trim, ")...")
  W <- data.frame(Y = Y, D = D, G = G, Ti = Ti)

  DID_D_hat <- D - nuis$m_D_10 - nuis$m_D_01 + nuis$m_D_00
  trim_info <- fddml_trim(nuis$pG_raw, DID_D = DID_D_hat,
                           method = trim, alpha_fixed = trim_alpha)

  # Apply trimming to propensity
  pT <- mean(Ti)
  nuis$pG_raw <- trim_info$pG_trimmed
  nuis$pi_11 <- nuis$pG_raw * pT
  nuis$pi_10 <- nuis$pG_raw * (1 - pT)
  nuis$pi_01 <- (1 - nuis$pG_raw) * pT
  nuis$pi_00 <- (1 - nuis$pG_raw) * (1 - pT)

  # Step 3: Estimation
  if (verbose) message("Step 3/4: Computing estimators...")
  wald_res <- tc_res <- NULL
  estimates <- data.frame()

  if (estimand %in% c("wald", "both")) {
    wald_res <- fddml_wald(W, nuis)
    wald_inf <- fddml_inference(wald_res, cluster = cluster, se_type = se_type, B = B)
    estimates <- rbind(estimates, data.frame(
      estimand = "DML-Wald",
      estimate = wald_inf$estimate,
      se = wald_inf$se,
      ci_lower = wald_inf$ci_lower,
      ci_upper = wald_inf$ci_upper,
      p_value = wald_inf$p_value,
      stringsAsFactors = FALSE
    ))
    wald_res <- c(wald_res, wald_inf[c("se", "ci_lower", "ci_upper", "p_value", "se_type")])
  }

  if (estimand %in% c("tc", "both")) {
    tc_res <- fddml_tc(W, nuis)
    tc_inf <- fddml_inference(tc_res, cluster = cluster, se_type = se_type, B = B)
    estimates <- rbind(estimates, data.frame(
      estimand = "DML-TC",
      estimate = tc_inf$estimate,
      se = tc_inf$se,
      ci_lower = tc_inf$ci_lower,
      ci_upper = tc_inf$ci_upper,
      p_value = tc_inf$p_value,
      stringsAsFactors = FALSE
    ))
    tc_res <- c(tc_res, tc_inf[c("se", "ci_lower", "ci_upper", "p_value", "se_type")])
  }

  if (verbose) message("Step 4/4: Done.")

  result <- list(
    estimates = estimates,
    wald = wald_res,
    tc = tc_res,
    nuisance = nuis,
    trim_info = trim_info,
    call = cl,
    N = N,
    p = p,
    settings = list(
      estimand = estimand, method = method, K = K,
      trim = trim, se_type = se_type, B = B
    )
  )
  class(result) <- "fddml"
  result
}
