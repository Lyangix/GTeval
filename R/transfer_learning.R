# Balanced folds within strata. Rotating the offset also balances small strata.
stratified_folds <- function(strata, k, seed) {
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

# Nested CV: leave PRS and treatments unpenalized; penalize only SNP effects.
cross_validate_transfer <- function(data, genotypes, outcome, prs,
                                     treatment_cols = character(), adjustment = ~ 1,
                                     group = NULL, strata = NULL,
                                     outer_folds = 5L, inner_folds = 5L,
                                     penalty = c("ridge", "lasso"),
                                     lambda = 10^seq(2, -4, length.out = 50L),
                                     seed = 42L, outer_fold_id = NULL) {
  penalty <- match.arg(penalty)
  lambda <- sort(unique(lambda), decreasing = TRUE)
  fixed <- as.matrix(data[, c(prs, treatment_cols), drop = FALSE])
  groups <- if (is.null(group)) rep("all", nrow(data)) else data[[group]]
  if (is.null(strata)) strata <- rep("all", nrow(data))
  strata <- interaction(groups, strata, drop = TRUE)
  if (is.null(outer_fold_id)) outer_fold_id <- stratified_folds(strata, outer_folds, seed)

  # All learned preprocessing uses only the current training rows.
  prepare <- function(train, test) {
    outcome_fit <- fit_outcome_transform(data[train, ], outcome, adjustment, group)
    fixed_prep <- fit_numeric_preprocessor(fixed[train, , drop = FALSE])
    snp_prep <- fit_numeric_preprocessor(genotypes[train, , drop = FALSE])
    f_train <- apply_numeric_preprocessor(fixed_prep, fixed[train, , drop = FALSE])
    f_test <- apply_numeric_preprocessor(fixed_prep, fixed[test, , drop = FALSE])
    z_train <- apply_numeric_preprocessor(snp_prep, genotypes[train, , drop = FALSE])
    z_test <- apply_numeric_preprocessor(snp_prep, genotypes[test, , drop = FALSE])
    list(y_train = outcome_fit$train_y, y_test = transform_outcome(outcome_fit, data[test, ]),
         x_train = cbind(f_train, z_train), x_test = cbind(f_test, z_test),
         f_train = cbind(Intercept = 1, f_train), f_test = cbind(Intercept = 1, f_test),
         treatment_index = c(1L, 1L + which(fixed_prep$keep != 1L)),
         penalty_factor = c(rep(0, ncol(f_train)), rep(1, ncol(z_train))))
  }
  fit_path <- function(d) {
    glmnet::glmnet(d$x_train, d$y_train, alpha = if (penalty == "ridge") 0 else 1,
      lambda = lambda, penalty.factor = d$penalty_factor, standardize = FALSE,
      intercept = TRUE, thresh = 1e-9, maxit = 100000L)
  }

  results <- predictions <- vector("list", outer_folds)
  for (k in seq_len(outer_folds)) {
    train <- which(outer_fold_id != k)
    test <- which(outer_fold_id == k)
    inner <- stratified_folds(strata[train], inner_folds, seed + k)
    mse <- matrix(NA_real_, inner_folds, length(lambda))
    for (j in seq_len(inner_folds)) {
      d <- prepare(train[inner != j], train[inner == j])
      predicted <- predict(fit_path(d), newx = d$x_test, s = lambda)
      mse[j, ] <- colMeans((predicted - d$y_test)^2)
    }
    cv_mse <- colMeans(mse)
    cv_se <- apply(mse, 2, sd) / sqrt(inner_folds)
    best <- which.min(cv_mse)
    one_se <- which(cv_mse <= cv_mse[best] + cv_se[best])[1L]

    # Refit on the outer training sample; evaluate once on its held-out sample.
    d <- prepare(train, test)
    pred_full <- predict(fit_path(d), newx = d$x_test, s = lambda[c(best, one_se)])
    baseline <- lm.fit(d$f_train, d$y_train)
    treatment <- lm.fit(d$f_train[, d$treatment_index, drop = FALSE], d$y_train)
    pred_prs <- as.numeric(d$f_test %*% baseline$coefficients)
    pred_treatment <- as.numeric(d$f_test[, d$treatment_index, drop = FALSE] %*% treatment$coefficients)
    r_treatment <- prediction_r2(d$y_test, pred_treatment)
    r_prs <- prediction_r2(d$y_test, pred_prs)
    r_min <- prediction_r2(d$y_test, pred_full[, 1L])
    r_1se <- prediction_r2(d$y_test, pred_full[, 2L])
    results[[k]] <- data.frame(fold = k, r2_treatment = r_treatment, r2_prs = r_prs,
      r2_transfer_min = r_min, r2_transfer_1se = r_1se, delta_r2_prs = r_prs - r_treatment,
      delta_r2_transfer_min = r_min - r_prs, delta_r2_transfer_1se = r_1se - r_prs,
      lambda_min = lambda[best], lambda_1se = lambda[one_se])
    predictions[[k]] <- data.frame(row = test, fold = k, transformed_y = d$y_test,
      treatment = pred_treatment, prs = pred_prs, transfer_min = pred_full[, 1L],
      transfer_1se = pred_full[, 2L])
  }
  per_fold <- do.call(rbind, results)
  metrics <- grep("r2", names(per_fold), value = TRUE)
  summary <- data.frame(metric = metrics, mean = sapply(per_fold[metrics], mean),
                         sd = sapply(per_fold[metrics], sd), row.names = NULL)
  predictions <- do.call(rbind, predictions)
  list(per_fold = per_fold, summary = summary,
       predictions = predictions[order(predictions$row), ], outer_fold_id = outer_fold_id)
}
