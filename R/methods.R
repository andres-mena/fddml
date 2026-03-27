#' Print an fddml object
#'
#' @param x An `fddml` object.
#' @param digits Number of significant digits (default 4).
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @method print fddml
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
    if (is.finite(est$estimate[i])) {
      cat(sprintf("  %-10s  %s (SE = %s)  95%% CI [%s, %s]\n",
                  est$estimand[i],
                  formatC(est$estimate[i], digits = digits, format = "f"),
                  formatC(est$se[i], digits = digits, format = "f"),
                  formatC(est$ci_lower[i], digits = digits, format = "f"),
                  formatC(est$ci_upper[i], digits = digits, format = "f")))
    } else {
      cat(sprintf("  %-10s  NA (estimation failed)\n", est$estimand[i]))
    }
  }
  cat("\n")
  invisible(x)
}

#' Summary of an fddml object
#'
#' @param object An `fddml` object.
#' @param ... Additional arguments (ignored).
#' @return An object of class `summary.fddml`, printed via [print.summary.fddml()].
#' @method summary fddml
#' @export
summary.fddml <- function(object, ...) {
  out <- list(
    call = object$call,
    N = object$N,
    p = object$p,
    estimates = object$estimates,
    settings = object$settings,
    trim_info = object$trim_info,
    wald_naive = if (!is.null(object$wald) && !is.null(object$wald$wald_naive))
      object$wald$wald_naive else NULL
  )
  class(out) <- "summary.fddml"
  out
}

#' Print a summary.fddml object
#'
#' @param x A `summary.fddml` object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @method print summary.fddml
#' @export
print.summary.fddml <- function(x, ...) {
  cat("\n  Fuzzy DID with Machine Learning - Summary\n")
  cat("  ", paste(rep("-", 50), collapse = ""), "\n")
  cat("  Call:", deparse(x$call, width.cutoff = 60), "\n")
  cat("\n  Sample: N =", x$N, ", p =", x$p, "\n")
  cat("  Method:", x$settings$method, "| Folds:", x$settings$K, "\n")
  cat("  SE type:", x$settings$se_type, "\n")
  cat("  Trimming:", x$settings$trim,
      if (x$trim_info$alpha_hat > 0) paste0("(alpha = ", round(x$trim_info$alpha_hat, 3), ")") else "",
      "\n\n")
  cat("  Estimates:\n")
  print(x$estimates, row.names = FALSE, digits = 4)
  cat("\n")
  if (!is.null(x$wald_naive) && is.finite(x$wald_naive))
    cat("  Wald naive (no DML):", round(x$wald_naive, 4), "\n")
  cat("\n")
  invisible(x)
}

#' @method confint fddml
#' @export
confint.fddml <- function(object, parm = NULL, level = 0.95, ...) {
  z <- stats::qnorm((1 + level) / 2)
  est <- object$estimates
  if (!is.null(parm)) {
    est <- est[est$estimand %in% parm, , drop = FALSE]
    if (nrow(est) == 0)
      stop("No matching estimand for parm = '", paste(parm, collapse = "', '"), "'.",
           call. = FALSE)
  }
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
