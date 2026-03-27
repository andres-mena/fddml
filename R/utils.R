# Internal helper functions (not exported)

#' Validate inputs for fddml estimation
#' @noRd
.validate_inputs <- function(Y, D, G, Ti, X) {
  N <- length(Y)
  if (length(D) != N) stop("`D` must have the same length as `Y`.", call. = FALSE)
  if (length(G) != N) stop("`G` must have the same length as `Y`.", call. = FALSE)
  if (length(Ti) != N) stop("`Ti` must have the same length as `Y`.", call. = FALSE)
  if (nrow(X) != N) stop("`X` must have the same number of rows as `Y`.", call. = FALSE)
  if (!all(G %in% 0:1)) stop("`G` must be binary (0/1).", call. = FALSE)
  if (!all(Ti %in% 0:1)) stop("`Ti` must be binary (0/1).", call. = FALSE)
  if (any(is.na(Y)) || any(is.na(D))) stop("Missing values in `Y` or `D`.", call. = FALSE)
  invisible(TRUE)
}

#' Create cross-fitting fold assignments
#' @noRd
.make_folds <- function(N, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep(seq_len(K), length.out = N))
}

#' Fit nuisance model within a (g,t) cell using specified method
#' @noRd
.fit_nuisance_cell <- function(Y_train, X_train, X_predict, method, family = "gaussian") {
  N_pred <- nrow(X_predict)
  fallback <- rep(mean(Y_train), N_pred)

  if (length(Y_train) < 10) return(fallback)

  tryCatch({
    if (method == "ols") {
      fit <- stats::lm.fit(x = cbind(1, X_train), y = Y_train)
      cbind(1, X_predict) %*% fit$coefficients

    } else if (method == "lasso") {
      if (family == "binomial") {
        mod <- glmnet::cv.glmnet(X_train, Y_train, alpha = 1, family = "binomial", nfolds = 5)
        as.vector(stats::predict(mod, newx = X_predict, s = "lambda.min", type = "response"))
      } else {
        mod <- glmnet::cv.glmnet(X_train, Y_train, alpha = 1, nfolds = 5)
        as.vector(stats::predict(mod, newx = X_predict, s = "lambda.min"))
      }

    } else if (method == "rf") {
      if (!requireNamespace("randomForest", quietly = TRUE))
        stop("Install 'randomForest' to use method = 'rf'.", call. = FALSE)
      df_train <- data.frame(y = Y_train, X_train)
      mod <- randomForest::randomForest(y ~ ., data = df_train, ntree = 500,
                                         mtry = max(1L, floor(ncol(X_train) / 3)),
                                         nodesize = 5L)
      as.vector(stats::predict(mod, newdata = data.frame(X_predict)))

    } else if (method == "nn") {
      if (!requireNamespace("nnet", quietly = TRUE))
        stop("Install 'nnet' to use method = 'nn'.", call. = FALSE)
      X_scaled <- scale(X_train)
      center <- attr(X_scaled, "scaled:center")
      sc <- attr(X_scaled, "scaled:scale")
      sc[sc == 0] <- 1
      mod <- nnet::nnet(X_scaled, Y_train, size = 5L, decay = 0.01,
                        linout = TRUE, trace = FALSE, maxit = 500)
      X_pred_scaled <- scale(X_predict, center = center, scale = sc)
      as.vector(stats::predict(mod, newdata = X_pred_scaled))

    } else if (is.function(method)) {
      method(Y_train, X_train, X_predict)

    } else {
      stop("Unknown method: '", method, "'. Use 'lasso', 'rf', 'nn', 'ols', or a function.",
           call. = FALSE)
    }
  }, error = function(e) {
    warning("Nuisance estimation failed in cell: ", e$message, ". Using cell mean.", call. = FALSE)
    fallback
  })
}
