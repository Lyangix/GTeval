# Learn column means and SDs from training rows; omit constant/all-missing columns.
fit_numeric_preprocessor <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  for (j in seq_len(ncol(x))) x[is.na(x[, j]), j] <- center[j]
  spread <- apply(x, 2, sd)
  list(center = center, scale = spread, keep = which(is.finite(spread) & spread > 1e-10))
}

apply_numeric_preprocessor <- function(object, x) {
  x <- x[, object$keep, drop = FALSE]
  for (j in seq_along(object$keep)) {
    k <- object$keep[j]
    x[is.na(x[, j]), j] <- object$center[k]
    x[, j] <- (x[, j] - object$center[k]) / object$scale[k]
  }
  x
}

# Covariate adjustment and rank inverse-normal transformation within each group.
fit_outcome_transform <- function(data, outcome, adjustment = ~ 1, group = NULL) {
  formula <- update(adjustment, paste(outcome, "~ ."))
  groups <- if (is.null(group)) rep("all", nrow(data)) else as.character(data[[group]])
  models <- list()
  y <- numeric(nrow(data))
  for (g in unique(groups)) {
    index <- which(groups == g)
    models[[g]] <- lm(formula, data[index, , drop = FALSE], na.action = na.fail)
    residual <- residuals(models[[g]])
    y[index] <- qnorm((rank(residual) - 0.5) / length(index))
  }
  list(train_y = y, models = models, outcome = outcome, group = group,
       epsilon = 1 / (2 * (nrow(data) + 1)))
}

# Evaluate held-out outcomes using the training regressions and residual ECDFs.
transform_outcome <- function(object, data) {
  groups <- if (is.null(object$group)) rep("all", nrow(data)) else as.character(data[[object$group]])
  y <- rep(NA_real_, nrow(data))
  for (g in unique(groups)) {
    index <- which(groups == g)
    fit <- object$models[[g]]
    residual <- data[[object$outcome]][index] - predict(fit, data[index, , drop = FALSE])
    p <- ecdf(residuals(fit))(residual)
    y[index] <- qnorm(pmax(object$epsilon, pmin(1 - object$epsilon, p)))
  }
  y
}

prediction_r2 <- function(y, predicted) {
  total <- sum((y - mean(y))^2)
  if (total <= 1e-12) return(NA_real_)
  1 - sum((y - predicted)^2) / total
}
