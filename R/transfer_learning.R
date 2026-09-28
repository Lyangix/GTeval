# Transfer model: transformed Y = intercept + a * external PRS + C * gamma + Z * delta.
# Only delta is penalized. The external PRS must not use target evaluation outcomes.

stratified_folds <- function(strata, k = 5L, seed = 42L) {
  if (anyNA(strata) || length(k) != 1L || !is.finite(k) ||
      k != as.integer(k) || k < 2L || k > length(strata)) stop("Invalid strata or fold count.")
  # Leave the caller's random-number state unchanged.
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  folds <- integer(length(strata))
  offset <- 0L
  for (s in unique(as.character(strata))) {
    index <- which(as.character(strata) == s)
    index <- index[sample.int(length(index))]
    folds[index] <- (offset + seq_along(index) - 1L) %% k + 1L
    offset <- (offset + length(index)) %% k
  }
  sample.int(k)[folds]
}

.fit_transfer_design <- function(prs, treatments, genotypes) {
  if (!is.numeric(prs) || any(!is.finite(prs)) || stats::sd(prs) <= 1e-10) {
    stop("PRS must be finite and vary in every training fold.")
  }
  treatments <- .numeric_matrix(treatments, "treatments")
  genotypes <- .numeric_matrix(genotypes, "genotypes", allow_na = TRUE)
  if (nrow(treatments) != length(prs) || nrow(genotypes) != length(prs)) {
    stop("PRS, treatments, and genotypes must have the same rows.")
  }
  fixed <- cbind(PRS = prs, treatments)
  fixed_prep <- fit_numeric_preprocessor(fixed)
  genotype_prep <- fit_numeric_preprocessor(genotypes, impute = TRUE)
  if (!length(genotype_prep$keep)) stop("No variable SNPs remain in this training fold.")
  f <- apply_numeric_preprocessor(fixed_prep, fixed)
  if (qr(cbind(1, f))$rank != ncol(f) + 1L) stop("PRS and treatment columns are collinear.")
  z <- apply_numeric_preprocessor(genotype_prep, genotypes)
  list(x = cbind(f, z), fixed = f, fixed_prep = fixed_prep,
       genotype_prep = genotype_prep,
       treatment_positions = which(fixed_prep$keep != 1L),
       penalty_factor = c(rep(0, ncol(f)), rep(1, ncol(z))))
}

.apply_transfer_design <- function(object, prs, treatments, genotypes) {
  f <- apply_numeric_preprocessor(object$fixed_prep, cbind(PRS = prs, treatments))
  z <- apply_numeric_preprocessor(object$genotype_prep, genotypes)
  if (nrow(f) != nrow(z)) stop("Predictor rows do not match.")
  list(x = cbind(f, z), fixed = f)
}

.validate_lambda <- function(lambda) {
  if (!is.numeric(lambda) || !length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("lambda must contain finite positive values.")
  }
  sort(unique(lambda), decreasing = TRUE)
}

.fit_transfer_path <- function(design, y, penalty, lambda) {
  .require_package("glmnet")
  fit <- glmnet::glmnet(design$x, y, alpha = if (penalty == "ridge") 0 else 1,
                        lambda = lambda, penalty.factor = design$penalty_factor,
                        standardize = FALSE, intercept = TRUE, thresh = 1e-9,
                        maxit = 100000L)
  if (fit$jerr != 0L || length(fit$lambda) != length(lambda)) {
    stop("glmnet did not fit the complete requested path; inspect convergence or lambda range.")
  }
  fit
}

# Fit a final model after choosing lambda in training-only CV. y is already on
# the desired analysis scale. Pass raw aligned SNP dosage, not globally imputed SNPs.
fit_transfer_model <- function(y, prs, genotypes, treatments = NULL,
                               penalty = c("ridge", "lasso"), lambda = 0.1) {
  penalty <- match.arg(penalty)
  lambda <- .validate_lambda(lambda)
  if (length(lambda) != 1L) stop("Select one lambda for the final model.")
  if (is.null(treatments)) treatments <- matrix(numeric(0), length(prs), 0L)
  if (!is.numeric(y) || length(y) != length(prs) || any(!is.finite(y)) ||
      stats::sd(y) <= 1e-10) stop("y must be finite, varying, and aligned with PRS.")
  design <- .fit_transfer_design(prs, treatments, genotypes)
  fit <- .fit_transfer_path(design, y, penalty, lambda)
  # Preprocessors and fitted coefficients suffice; do not retain training genotypes.
  design$x <- design$fixed <- NULL
  list(fit = fit, design = design, lambda = lambda, penalty = penalty)
}

predict_transfer_model <- function(model, prs, genotypes, treatments = NULL) {
  if (is.null(treatments)) treatments <- matrix(numeric(0), length(prs), 0L)
  design <- .apply_transfer_design(model$design, prs, treatments, genotypes)
  as.numeric(stats::predict(model$fit, newx = design$x, s = model$lambda))
}

# Nested CV refits outcome adjustment, ECDFs, dosage imputation, and predictor
# scaling in every inner and outer training set. It returns held-out metrics only.
cross_validate_transfer <- function(data, genotypes, outcome, prs,
                                     treatment_cols = character(), adjustment = ~ 1,
                                     group = NULL, strata = NULL,
                                     outer_folds = 5L, inner_folds = 5L,
                                     penalty = c("ridge", "lasso"),
                                     lambda = 10^seq(2, -4, length.out = 50L),
                                     seed = 42L, outer_fold_id = NULL) {
  .require_package("glmnet")
  penalty <- match.arg(penalty)
  lambda <- .validate_lambda(lambda)
  if (length(lambda) < 2L) stop("Supply at least two candidate lambdas for CV.")
  n <- nrow(data)
  genotypes <- .numeric_matrix(genotypes, "genotypes", allow_na = TRUE)
  if (nrow(genotypes) != n || ncol(genotypes) < 1L || length(prs) != 1L ||
      !all(c(prs, treatment_cols) %in% names(data)) ||
      anyDuplicated(c(prs, treatment_cols))) stop("Check predictor columns and row alignment.")
  fixed <- .numeric_matrix(data[, c(prs, treatment_cols), drop = FALSE], "fixed predictors")
  score <- fixed[, 1L]
  treatments <- fixed[, -1L, drop = FALSE]
  groups <- .outcome_groups(data, group)
  if (is.null(strata)) strata <- rep("all", n)
  if (length(strata) != n || anyNA(strata)) stop("strata must have one complete entry per row.")
  strata <- interaction(groups, strata, drop = TRUE)
  if (is.null(outer_fold_id)) {
    outer_fold_id <- stratified_folds(strata, outer_folds, seed)
  } else {
    if (!is.numeric(outer_fold_id) || length(outer_fold_id) != n || anyNA(outer_fold_id) ||
        !identical(sort(unique(as.integer(outer_fold_id))), seq_len(outer_folds)) ||
        any(outer_fold_id != as.integer(outer_fold_id))) stop("Invalid outer_fold_id.")
  }
  prepare <- function(train, test) {
    transform <- fit_outcome_transform(data[train, , drop = FALSE], outcome, adjustment, group)
    design <- .fit_transfer_design(score[train], treatments[train, , drop = FALSE],
                                   genotypes[train, , drop = FALSE])
    list(y_train = transform$train_y,
         y_test = transform_outcome(transform, data[test, , drop = FALSE]),
         design = design,
         test_design = .apply_transfer_design(design, score[test],
           treatments[test, , drop = FALSE], genotypes[test, , drop = FALSE]))
  }
  rows <- tuning <- predictions <- vector("list", outer_folds)
  for (k in seq_len(outer_folds)) {
    train <- which(outer_fold_id != k)
    test <- which(outer_fold_id == k)
    inner <- stratified_folds(strata[train], inner_folds, seed + k)
    mse <- matrix(NA_real_, inner_folds, length(lambda))
    for (j in seq_len(inner_folds)) {
      prepared <- prepare(train[inner != j], train[inner == j])
      fit <- .fit_transfer_path(prepared$design, prepared$y_train, penalty, lambda)
      pred <- stats::predict(fit, newx = prepared$test_design$x, s = lambda)
      mse[j, ] <- colMeans((pred - prepared$y_test)^2)
    }
    mean_mse <- colMeans(mse)
    se_mse <- apply(mse, 2L, stats::sd) / sqrt(inner_folds)
    best <- which.min(mean_mse)
    conservative <- which(mean_mse <= mean_mse[best] + se_mse[best])[1L]
    prepared <- prepare(train, test)
    full <- .fit_transfer_path(prepared$design, prepared$y_train, penalty, lambda)
    pred_full <- stats::predict(full, newx = prepared$test_design$x,
                                s = lambda[c(best, conservative)])
    base_x <- cbind(1, prepared$design$fixed)
    base_test <- cbind(1, prepared$test_design$fixed)
    treat_index <- c(1L, 1L + prepared$design$treatment_positions)
    base_fit <- stats::lm.fit(base_x, prepared$y_train)
    treat_fit <- stats::lm.fit(base_x[, treat_index, drop = FALSE], prepared$y_train)
    pred_base <- as.numeric(base_test %*% base_fit$coefficients)
    pred_treat <- as.numeric(base_test[, treat_index, drop = FALSE] %*% treat_fit$coefficients)
    y <- prepared$y_test
    r_treat <- prediction_r2(y, pred_treat)
    r_base <- prediction_r2(y, pred_base)
    r_min <- prediction_r2(y, pred_full[, 1L])
    r_1se <- prediction_r2(y, pred_full[, 2L])
    rows[[k]] <- data.frame(fold = k, n_train = length(train), n_test = length(test),
      n_snps = length(prepared$design$genotype_prep$keep),
      r2_treatment = r_treat, r2_prs = r_base,
      r2_transfer_min = r_min, r2_transfer_1se = r_1se,
      delta_r2_prs = r_base - r_treat,
      delta_r2_transfer_min = r_min - r_base, delta_r2_transfer_1se = r_1se - r_base,
      lambda_min = lambda[best], lambda_1se = lambda[conservative],
      lambda_min_at_edge = best %in% c(1L, length(lambda)))
    tuning[[k]] <- list(lambda = lambda, mean_mse = mean_mse, se_mse = se_mse,
                        fold_mse = mse, inner_fold_id = inner, outer_train_rows = train)
    predictions[[k]] <- data.frame(row = test, fold = k, transformed_y = y,
      treatment = pred_treat, prs = pred_base, transfer_min = pred_full[, 1L],
      transfer_1se = pred_full[, 2L])
  }
  per_fold <- do.call(rbind, rows)
  metrics <- c("r2_treatment", "r2_prs", "r2_transfer_min", "r2_transfer_1se",
               "delta_r2_prs", "delta_r2_transfer_min", "delta_r2_transfer_1se")
  summary <- data.frame(metric = metrics,
    mean = vapply(per_fold[metrics], mean, numeric(1)),
    sd = vapply(per_fold[metrics], stats::sd, numeric(1)), row.names = NULL)
  predictions <- do.call(rbind, predictions)
  predictions <- predictions[order(predictions$row), ]
  rownames(predictions) <- NULL
  list(per_fold = per_fold, summary = summary, predictions = predictions,
       tuning = tuning, outer_fold_id = outer_fold_id,
       options = list(penalty = penalty, lambda = lambda, seed = seed,
                      outer_folds = outer_folds, inner_folds = inner_folds))
}
