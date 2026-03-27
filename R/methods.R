#' @export
print.fddml <- function(x, digits = 4, ...) {
  cat("\n  Fuzzy DID with Machine Learning\n\n")
  cat("  N =", x$N, " p =", x$p, " K =", x$settings$K,
      " method =", x$settings$method, "\n")
  if (x$trim_info$alpha_hat > 0)
    cat("  Trimming: alpha =", round(x$trim_info$alpha_hat, 3), "\n")
  cat("\n")
  est <- x$estimates
  for (i in seq_len(nrow(est))) {
    cat(sprintf("  %-10s  %s (SE = %s)  95%% CI [%s, %s]\n",
                est$estimand[i],
                formatC(est$estimate[i], digits = digits, format = "f"),
                formatC(est$se[i], digits = digits, format = "f"),
                formatC(est$ci_lower[i], digits = digits, format = "f"),
                formatC(est$ci_upper[i], digits = digits, format = "f")))
  }
  cat("\n")
  invisible(x)
}

#' @export
summary.fddml <- function(object, ...) {
  cat("\n  Fuzzy DID with Machine Learning — Summary\n")
  cat("  ", paste(rep("-", 50), collapse = ""), "\n")
  cat("  Call: "); print(object$call)
  cat("\n  Sample: N =", object$N, ", p =", object$p, "\n")
  cat("  Method:", object$settings$method, "| Folds:", object$settings$K, "\n")
  cat("  SE type:", object$settings$se_type, "\n")
  cat("  Trimming:", object$settings$trim,
      if (object$trim_info$alpha_hat > 0) paste0("(alpha = ", round(object$trim_info$alpha_hat, 3), ")") else "",
      "\n\n")
  cat("  Estimates:\n")
  print(object$estimates, row.names = FALSE, digits = 4)
  cat("\n")
  if (!is.null(object$wald))
    cat("  Wald naive (no DML):", round(object$wald$wald_naive, 4), "\n")
  cat("\n")
  invisible(object)
}

#' @export
confint.fddml <- function(object, parm, level = 0.95, ...) {
  z <- stats::qnorm((1 + level) / 2)
  est <- object$estimates
  ci <- data.frame(
    lower = est$estimate - z * est$se,
    upper = est$estimate + z * est$se,
    row.names = est$estimand
  )
  names(ci) <- paste0(formatC(100 * c((1 - level) / 2, (1 + level) / 2),
                               format = "f", digits = 1), " %")
  ci
}

#' Tidy an fddml object (broom-style)
#'
#' @param x An `fddml` object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with columns: term, estimate, std.error, statistic,
#'   p.value, conf.low, conf.high.
#' @export
tidy.fddml <- function(x, ...) {
  est <- x$estimates
  data.frame(
    term = est$estimand,
    estimate = est$estimate,
    std.error = est$se,
    statistic = est$estimate / est$se,
    p.value = est$p_value,
    conf.low = est$ci_lower,
    conf.high = est$ci_upper,
    stringsAsFactors = FALSE
  )
}
