#' Plot Diagnostics for fddml Results
#'
#' Produces a propensity score distribution plot showing the estimated
#' Pr(G=1|X) with the trimming threshold marked.
#'
#' @param x An `fddml` object.
#' @param type Character: `"propensity"` (default) or `"scores"`.
#' @param ... Additional arguments passed to ggplot.
#'
#' @return A `ggplot` object.
#'
#' @export
plot.fddml <- function(x, type = "propensity", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Install 'ggplot2' for plotting.", call. = FALSE)

  if (type == "propensity") {
    pG_display <- if (!is.null(x$nuisance$pG_raw_original))
      x$nuisance$pG_raw_original else x$nuisance$pG_raw
    df <- data.frame(pG = pG_display)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = pG)) +
      ggplot2::geom_histogram(fill = "#012169", color = "white",
                               bins = 40, alpha = 0.8) +
      ggplot2::labs(x = expression(hat(e)(X) == Pr(G == 1 ~ "|" ~ X)),
                    y = "Count",
                    title = "Estimated propensity score distribution") +
      ggplot2::theme_minimal(base_size = 13)

    if (x$trim_info$alpha_hat > 0) {
      p <- p + ggplot2::geom_vline(
        xintercept = 1 - x$trim_info$alpha_hat,
        linetype = "dashed", color = "#b91c1c", linewidth = 0.8
      ) +
        ggplot2::annotate("text",
                          x = 1 - x$trim_info$alpha_hat - 0.03,
                          y = Inf, vjust = 2,
                          label = paste0("alpha == ", x$trim_info$alpha_hat),
                          parse = TRUE, color = "#b91c1c", size = 4)
    }
    return(p)

  } else if (type == "scores") {
    scores_list <- list()
    if (!is.null(x$wald))
      scores_list[["DML-Wald"]] <- data.frame(score = x$wald$scores, estimand = "DML-Wald")
    if (!is.null(x$tc))
      scores_list[["DML-TC"]] <- data.frame(score = x$tc$scores, estimand = "DML-TC")
    df <- do.call(rbind, scores_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = score)) +
      ggplot2::geom_histogram(fill = "#012169", color = "white",
                               bins = 50, alpha = 0.8) +
      ggplot2::facet_wrap(~estimand, scales = "free") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#b91c1c") +
      ggplot2::labs(x = expression(psi[i](hat(theta))), y = "Count",
                    title = "Score distribution at estimated theta") +
      ggplot2::theme_minimal(base_size = 13)
    return(p)
  }

  stop("Unknown plot type: '", type, "'. Use 'propensity' or 'scores'.", call. = FALSE)
}
