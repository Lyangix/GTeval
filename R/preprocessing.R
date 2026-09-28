# Small, in-memory helpers. No function in this release reads study data.

.require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the '", package, "' package first.", call. = FALSE)
  }
}

.numeric_matrix <- function(x, label = "x", allow_na = FALSE) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || nrow(x) < 1L ||
      any(is.infinite(x)) || (!allow_na && anyNA(x))) {
    stop(label, " must be a numeric matrix with valid entries.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

# Explicit identifiers prevent silent row mismatches after phenotype merges.
align_genotypes <- function(genotypes, genotype_ids, sample_ids) {
  genotypes <- .numeric_matrix(genotypes, "genotypes", allow_na = TRUE)
  valid_ids <- function(ids) {
    length(ids) > 0L && !anyNA(ids) && !anyDuplicated(ids) &&
      all(nzchar(as.character(ids)))
  }
  if (length(genotype_ids) != nrow(genotypes) ||
      !valid_ids(genotype_ids) || !valid_ids(sample_ids)) {
    stop("Provide unique, nonmissing genotype_ids and sample_ids.")
  }
  index <- match(sample_ids, genotype_ids)
  if (anyNA(index)) stop("Some sample_ids have no genotype row.")
  result <- genotypes[index, , drop = FALSE]
  rownames(result) <- as.character(sample_ids)
  result
}

# Fit on training rows only. All-missing and constant columns are excluded.
fit_numeric_preprocessor <- function(x, impute = FALSE) {
  x <- .numeric_matrix(x, allow_na = impute)
  center <- colMeans(x, na.rm = impute)
  filled <- x
  for (j in seq_len(ncol(x))) filled[is.na(filled[, j]), j] <- center[j]
  scale <- apply(filled, 2L, stats::sd)
  keep <- which(is.finite(center) & is.finite(scale) & scale > 1e-10)
  list(center = center, scale = scale, keep = keep, n_columns = ncol(x),
       column_names = colnames(x), impute = impute)
}

apply_numeric_preprocessor <- function(object, x) {
  x <- .numeric_matrix(x, allow_na = object$impute)
  if (ncol(x) != object$n_columns ||
      !identical(colnames(x), object$column_names)) {
    stop("Columns and their order must match the training matrix.")
  }
  x <- x[, object$keep, drop = FALSE]
  for (j in seq_along(object$keep)) {
    original <- object$keep[j]
    x[is.na(x[, j]), j] <- object$center[original]
    x[, j] <- (x[, j] - object$center[original]) / object$scale[original]
  }
  x
}

.outcome_groups <- function(data, group) {
  if (is.null(group)) return(rep("all", nrow(data)))
  if (length(group) != 1L || !group %in% names(data) || anyNA(data[[group]])) {
    stop("group must name a complete grouping column, or be NULL.")
  }
  as.character(data[[group]])
}

# Regress the phenotype on adjustment covariates separately in each group,
# then rank inverse-normal transform the training residuals.
fit_outcome_transform <- function(data, outcome, adjustment = ~ 1, group = NULL) {
  if (!inherits(adjustment, "formula") || length(adjustment) != 2L) {
    stop("adjustment must be a one-sided formula, e.g. ~ age + I(age^2).")
  }
  if (length(outcome) != 1L || !outcome %in% names(data) ||
      !is.numeric(data[[outcome]]) || any(!is.finite(data[[outcome]]))) {
    stop("outcome must name a complete, finite numeric column.")
  }
  frame <- stats::model.frame(adjustment, data, na.action = stats::na.fail)
  design <- .numeric_matrix(stats::model.matrix(adjustment, frame))
  groups <- .outcome_groups(data, group)
  if (nrow(design) != nrow(data)) stop("Adjustment rows do not match data.")
  models <- list()
  transformed <- numeric(nrow(data))
  for (g in unique(groups)) {
    index <- which(groups == g)
    x <- design[index, , drop = FALSE]
    if (length(index) <= ncol(x) + 1L || qr(x)$rank != ncol(x)) {
      stop("Insufficient rows or rank-deficient adjustment design in group '", g, "'.")
    }
    fit <- stats::lm.fit(x, data[[outcome]][index])
    residuals <- as.numeric(fit$residuals)
    if (stats::sd(residuals) <= 1e-10) stop("No residual variation in group '", g, "'.")
    transformed[index] <- stats::qnorm((rank(residuals, ties.method = "average") - 0.5) /
                                       length(residuals))
    models[[g]] <- list(beta = fit$coefficients, residuals = residuals)
  }
  factors <- vapply(frame, is.factor, logical(1))
  list(train_y = transformed, models = models, outcome = outcome, group = group,
       adjustment = adjustment, columns = colnames(design),
       xlevels = lapply(frame[factors], levels), contrasts = attr(design, "contrasts"),
       epsilon = 1 / (2 * (nrow(data) + 1)))
}

# Held-out outcomes are evaluated using training coefficients and training ECDFs.
# Use object$train_y for the training sample (midranks, as in the source scripts).
transform_outcome <- function(object, data) {
  y <- data[[object$outcome]]
  if (!is.numeric(y) || length(y) != nrow(data) || any(!is.finite(y))) {
    stop("Held-out outcomes must be complete and finite.")
  }
  frame <- stats::model.frame(object$adjustment, data, xlev = object$xlevels,
                              na.action = stats::na.fail)
  x <- .numeric_matrix(stats::model.matrix(object$adjustment, frame,
                                          contrasts.arg = object$contrasts))
  if (!identical(colnames(x), object$columns)) stop("Adjustment design changed.")
  groups <- .outcome_groups(data, object$group)
  if (any(!groups %in% names(object$models))) stop("A held-out group is absent from training.")
  result <- numeric(nrow(data))
  for (g in unique(groups)) {
    index <- which(groups == g)
    fit <- object$models[[g]]
    residuals <- y[index] - as.numeric(x[index, , drop = FALSE] %*% fit$beta)
    p <- stats::ecdf(fit$residuals)(residuals)
    result[index] <- stats::qnorm(pmax(object$epsilon, pmin(1 - object$epsilon, p)))
  }
  result
}

# Out-of-sample R2 may be negative. Constant outcomes have undefined R2.
prediction_r2 <- function(y, predicted) {
  if (length(y) != length(predicted) || any(!is.finite(c(y, predicted)))) {
    stop("y and predicted must be finite vectors of equal length.")
  }
  denominator <- sum((y - mean(y))^2)
  if (denominator <= 1e-12) return(NA_real_)
  1 - sum((y - predicted)^2) / denominator
}
