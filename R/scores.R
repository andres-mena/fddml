#' Compute DML-Wald Scores
#'
#' Computes the locally robust augmented score for the Wald-DID estimand.
#' The score is psi = g + phi where g is the identifying moment and phi
#' is the FSIF correction (Ichimura & Newey 2022). Sign convention: (-, -, +)
#' for correction terms w_10, w_01, w_00, matching DID = Y - m_10 - m_01 + m_00.
#'
#' @param W Data frame with columns Y, D, G, Ti.
#' @param pred List of cross-fitted nuisance predictions from \code{fddml_nuisance()}.
#'
#' @return List with `estimate`, `se`, `ci_lower`, `ci_upper`, `variance`,
#'   `scores` (individual psi_i), `wald_naive`.
#'
#' @export
fddml_wald <- function(W, pred) {
  N <- nrow(W)
  Y <- W$Y; D <- W$D; G <- W$G; Ti <- W$Ti
  GT <- G * Ti
  p_11 <- mean(GT)

  # DID of conditional expectations
  DID_Y <- Y - pred$m_Y_10 - pred$m_Y_01 + pred$m_Y_00
  DID_D <- D - pred$m_D_10 - pred$m_D_01 + pred$m_D_00

  wald_naive <- sum(GT * DID_Y) / sum(GT * DID_D)

  # Riesz representer weights (CD'H Eq 68)
  w_10 <- G * (1 - Ti) * pred$pi_11 / (pred$pi_10 * p_11)
  w_01 <- (1 - G) * Ti * pred$pi_11 / (pred$pi_01 * p_11)
  w_00 <- (1 - G) * (1 - Ti) * pred$pi_11 / (pred$pi_00 * p_11)

  # Residuals
  r_Y_10 <- Y - pred$m_Y_10; r_Y_01 <- Y - pred$m_Y_01; r_Y_00 <- Y - pred$m_Y_00
  r_D_10 <- D - pred$m_D_10; r_D_01 <- D - pred$m_D_01; r_D_00 <- D - pred$m_D_00

  # Closed-form GMM: numerator and denominator
  num_i <- (GT / p_11) * DID_Y - w_10 * r_Y_10 - w_01 * r_Y_01 + w_00 * r_Y_00
  den_i <- (GT / p_11) * DID_D - w_10 * r_D_10 - w_01 * r_D_01 + w_00 * r_D_00

  theta_hat <- sum(num_i) / sum(den_i)

  # Scores at theta_hat
  psi_i <- (GT / p_11) * (DID_Y - DID_D * theta_hat) -
    w_10 * (r_Y_10 - r_D_10 * theta_hat) -
    w_01 * (r_Y_01 - r_D_01 * theta_hat) +
    w_00 * (r_Y_00 - r_D_00 * theta_hat)

  # Variance (Theorem 1)
  G_hat <- -mean(den_i)
  V_hat <- mean(psi_i^2) / (G_hat^2)
  se_hat <- sqrt(V_hat / N)

  list(
    estimate = theta_hat,
    se = se_hat,
    ci_lower = theta_hat - 1.96 * se_hat,
    ci_upper = theta_hat + 1.96 * se_hat,
    variance = V_hat,
    scores = psi_i,
    num_i = num_i,
    den_i = den_i,
    wald_naive = wald_naive
  )
}


#' Compute DML-TC Scores
#'
#' Computes the locally robust augmented score for the Time-Corrected estimand.
#' Equivalent to CD'H (2018) supplement p.45, psi*_TC(g), up to normalization.
#' Requires treatment-conditional nuisance functions in `pred`.
#'
#' @param W Data frame with columns Y, D, G, Ti.
#' @param pred List of cross-fitted nuisance predictions from [fddml_nuisance()]
#'   with `estimand = "tc"` or `"both"`.
#'
#' @return List with `estimate`, `se`, `ci_lower`, `ci_upper`, `variance`,
#'   `scores` (individual psi_i).
#'
#' @export
fddml_tc <- function(W, pred) {
  N <- nrow(W)
  Y <- W$Y; D <- W$D; G <- W$G; Ti <- W$Ti
  GT <- G * Ti
  p_11 <- mean(GT)

  m_Y_10 <- pred$m_Y_10; m_D_10 <- pred$m_D_10

  # Time trends by treatment status in control group
  delta_1 <- pred$mu_Y_101 - pred$mu_Y_100
  delta_0 <- pred$mu_Y_001 - pred$mu_Y_000

  # TC moment components
  g_Y <- Y - m_Y_10 - delta_1 * m_D_10 - (1 - m_D_10) * delta_0
  g_D <- D - m_D_10

  # Propensities and weights
  pi_11 <- pred$pi_11; pi_10 <- pred$pi_10
  w_10 <- G * (1 - Ti) * pi_11 / (pi_10 * p_11)
  w_101 <- D * (1 - G) * Ti * pi_11 / (pred$pi_101 * p_11)
  w_100 <- D * (1 - G) * (1 - Ti) * pi_11 / (pred$pi_100 * p_11)
  w_001 <- (1 - D) * (1 - G) * Ti * pi_11 / (pred$pi_001 * p_11)
  w_000 <- (1 - D) * (1 - G) * (1 - Ti) * pi_11 / (pred$pi_000 * p_11)

  # Residuals
  r_Y_10 <- Y - m_Y_10; r_D_10 <- D - m_D_10
  r_101 <- Y - pred$mu_Y_101; r_100 <- Y - pred$mu_Y_100
  r_001 <- Y - pred$mu_Y_001; r_000 <- Y - pred$mu_Y_000

  # Closed-form GMM
  num_i <- (GT / p_11) * g_Y -
    w_10 * (r_Y_10 - r_D_10 * (delta_0 - delta_1)) -
    m_D_10 * (w_101 * r_101 - w_100 * r_100) -
    (1 - m_D_10) * (w_001 * r_001 - w_000 * r_000)

  den_i <- (GT / p_11) * g_D - w_10 * r_D_10

  theta_hat <- sum(num_i) / sum(den_i)

  # Scores at theta_hat
  psi_i <- (GT / p_11) * (g_Y - g_D * theta_hat) -
    w_10 * (r_Y_10 - r_D_10 * (delta_0 - delta_1 + theta_hat)) -
    m_D_10 * (w_101 * r_101 - w_100 * r_100) -
    (1 - m_D_10) * (w_001 * r_001 - w_000 * r_000)

  G_hat <- -mean(den_i)
  V_hat <- mean(psi_i^2) / (G_hat^2)
  se_hat <- sqrt(V_hat / N)

  list(
    estimate = theta_hat,
    se = se_hat,
    ci_lower = theta_hat - 1.96 * se_hat,
    ci_upper = theta_hat + 1.96 * se_hat,
    variance = V_hat,
    scores = psi_i,
    num_i = num_i,
    den_i = den_i
  )
}
