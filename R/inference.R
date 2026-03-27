#' Inference for fddml Results
#'
#' Computes standard errors, confidence intervals, and p-values for
#' DML-Wald and DML-TC estimates. Supports analytical (i.i.d.), clustered,
#' and multiplier bootstrap variance estimation.
#'
#' @param result An `fddml` result object or a list from \code{fddml_wald()} or \code{fddml_tc()}.
#' @param cluster Optional numeric or factor vector of cluster identifiers
#'   (same length as the data). When provided, computes cluster-robust SEs.
#' @param se_type Character: `"analytical"` (default), `"cluster"` (requires
#'   `cluster`), or `"bootstrap"`.
#' @param B Integer, number of bootstrap replications (default 1000).
#' @param alpha Significance level for confidence intervals (default 0.05).
#'
#' @return A list with updated `se`, `ci_lower`, `ci_upper`, `p_value`, and
#'   `se_type`.
#'
#' @export
fddml_inference <- function(result, cluster = NULL, se_type = "analytical",
                             B = 1000L, alpha = 0.05) {
  psi_i <- result$scores
  N <- length(psi_i)
  theta_hat <- result$estimate
  G_hat <- -mean(result$den_i)
  z_alpha <- stats::qnorm(1 - alpha / 2)

  if (se_type == "cluster" || !is.null(cluster)) {
    if (is.null(cluster)) stop("`cluster` must be provided for se_type = 'cluster'.", call. = FALSE)
    if (length(cluster) != N) stop("`cluster` must have the same length as the data.", call. = FALSE)

    # Cluster-robust variance: V = G^{-2} * (1/N^2) * sum_s (sum_{i in s} psi_i)^2
    cluster_sums <- tapply(psi_i, cluster, sum)
    V_hat <- sum(cluster_sums^2) / (N^2 * G_hat^2)
    se_hat <- sqrt(V_hat)
    se_type_out <- "cluster"

  } else if (se_type == "bootstrap") {
    # Multiplier bootstrap: resample scores with Rademacher weights
    boot_theta <- numeric(B)
    for (b in seq_len(B)) {
      w_b <- sample(c(-1, 1), N, replace = TRUE)
      boot_theta[b] <- theta_hat + sum(w_b * psi_i) / (N * G_hat)
    }
    se_hat <- stats::sd(boot_theta)
    V_hat <- se_hat^2
    se_type_out <- "bootstrap"

  } else {
    # Analytical (i.i.d.) — already computed in scores
    V_hat <- result$variance
    se_hat <- sqrt(V_hat / N)
    se_type_out <- "analytical"
  }

  ci_lower <- theta_hat - z_alpha * se_hat
  ci_upper <- theta_hat + z_alpha * se_hat
  p_value <- 2 * stats::pnorm(-abs(theta_hat / se_hat))

  list(
    estimate = theta_hat,
    se = se_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    variance = V_hat,
    se_type = se_type_out,
    alpha = alpha
  )
}
