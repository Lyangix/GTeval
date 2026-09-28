#!/usr/bin/env Rscript
# Integration and statistical invariants; uses only simulated data.
for (file in c("preprocessing", "mixed_model", "transfer_learning", "simulation")) {
  source(file.path("R", paste0(file, ".R")))
}
.require_package("gaston")
.require_package("glmnet")
assert_equal <- function(a, b, tolerance = 1e-7) {
  stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)))
}
assert_error <- function(expr) {
  caught <- tryCatch({force(expr); FALSE}, error = function(e) TRUE)
  stopifnot(caught)
}

sim <- simulate_example_data(n = 150L, p = 25L, seed = 3L, missing_rate = 0.03)
dat <- sim$data
assert_equal(sim$genotypes, simulate_example_data(n = 150L, p = 25L, seed = 3L,
                                                  missing_rate = 0.03)$genotypes)
order <- rev(seq_len(nrow(dat)))
assert_equal(align_genotypes(sim$genotypes[order, ], dat$id[order], dat$id), sim$genotypes)
assert_error(align_genotypes(sim$genotypes, rep("duplicate", nrow(dat)), dat$id))
assert_error(align_genotypes(sim$genotypes, dat$id, c(dat$id[-1L], "unmatched")))

train <- matrix(c(0, 1, NA, 2, 1, 1, 1, 1, NA, NA, NA, NA), 4L, 3L)
prep <- fit_numeric_preprocessor(train, impute = TRUE)
stopifnot(identical(prep$keep, 1L))
assert_equal(apply_numeric_preprocessor(prep, matrix(c(NA, 9, 8), 1L)), matrix(0, 1L))
assert_error(fit_numeric_preprocessor(matrix(c(1, Inf), 2L)))

# Numerical equivalence to the source kernel definitions when no SNP is excluded.
kernels <- build_kernels(sim$complete_genotypes, dat$treatment)
assert_equal(kernels$G, tcrossprod(scale(sim$complete_genotypes)) / 25)
assert_equal(kernels$GxT, tcrossprod(scale(sim$complete_genotypes * dat$treatment)) / 25)
stopifnot(min(eigen(kernels$GxT, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
augmented <- cbind(sim$genotypes, constant = 1, missing = NA_real_)
augmented_k <- build_kernels(augmented, dat$treatment)
stopifnot(attr(augmented_k, "snp_counts")["G"] == 25L)
assert_error(build_kernels(sim$genotypes, rep(0, nrow(dat))))

# Group-specific adjustment preserves row order and training midranks.
transform <- fit_outcome_transform(dat, "outcome", ~ age + I(age^2) + PC1 + PC2, "sex")
for (g in c(0, 1)) {
  i <- which(dat$sex == g)
  residual <- residuals(lm(outcome ~ age + I(age^2) + PC1 + PC2, dat[i, ]))
  assert_equal(transform$train_y[i], qnorm((rank(residual) - 0.5) / length(i)))
}
assert_error(fit_outcome_transform(transform(dat, age = NA_real_), "outcome", ~ age, "sex"))

X <- cbind(Intercept = 1, PRS = dat$prs, Treatment = dat$treatment)
model <- fit_variance_components(sim$y_model, X, kernels)
reference <- gaston::lmm.aireml(sim$y_model, X, unname(kernels), min_tau = 0,
                               max_iter = 200L, eps = 1e-7, verbose = FALSE)
assert_equal(model$components$variance, c(reference$tau, reference$sigma2))
assert_equal(sum(model$components$conditional_share), 1)
assert_equal(model$logLik, gaston::lmm.restricted.likelihood(sim$y_model, X,
  unname(kernels), tau = reference$tau, s2 = reference$sigma2))
null_test <- test_variance_component(model, "GxT")
stopifnot(null_test$statistic >= 0, null_test$p_value >= 0, null_test$p_value <= 1)
assert_equal(.lrt_result(-10, -10, 0)["p_value"], 1)
profile <- profile_variance_component(model, "GxT", c(0, reference$tau[1L], 0.8))
stopifnot(length(profile$accepted) >= 1L,
          profile$table$p_value[which.min(abs(profile$table$null_value - reference$tau[1L]))] > 0.99)
assert_error(fit_variance_components(sim$y_model, cbind(X, X[, 2L]), kernels))
assert_error(fit_variance_components(sim$y_model, X, list(bad = -diag(nrow(dat)))))
comparisons <- compare_prs_models(transform$train_y, dat$prs, dat$treatment)
stopifnot(all(diff(comparisons$r_squared) >= -1e-10))

# Folds are balanced even when each stratum has only one observation.
folds <- stratified_folds(seq_len(17), 5L)
stopifnot(length(unique(folds)) == 5L, diff(range(table(folds))) <= 1L)
assert_equal(prediction_r2(c(0, 1), c(2, 2)), -9)
stopifnot(is.na(prediction_r2(c(1, 1), c(0, 0))))

cv_args <- list(data = dat, genotypes = sim$genotypes, outcome = "outcome", prs = "prs",
  treatment_cols = c("treatment", "chemo"), adjustment = ~ age + I(age^2) + PC1 + PC2,
  group = "sex", strata = interaction(dat$treatment, dat$chemo),
  outer_folds = 3L, inner_folds = 2L, lambda = c(10, 1, 0.1, 0.01), seed = 19L)
cv <- do.call(cross_validate_transfer, cv_args)
stopifnot(nrow(cv$predictions) == nrow(dat), all(is.finite(as.matrix(cv$predictions))),
          all(cv$per_fold$lambda_1se >= cv$per_fold$lambda_min))

# Changing outer-test outcomes cannot change that fold's tuning or predictions.
altered_args <- cv_args
heldout <- which(cv$outer_fold_id == 1L)
altered_args$data$outcome[heldout] <- altered_args$data$outcome[heldout] + 4
altered_args$outer_fold_id <- cv$outer_fold_id
altered <- do.call(cross_validate_transfer, altered_args)
assert_equal(cv$tuning[[1L]]$mean_mse, altered$tuning[[1L]]$mean_mse)
assert_equal(cv$predictions[heldout, c("prs", "transfer_min", "transfer_1se")],
             altered$predictions[heldout, c("prs", "transfer_min", "transfer_1se")])
stopifnot(any(cv$predictions$transformed_y[heldout] != altered$predictions$transformed_y[heldout]))

# A held-out dosage change also cannot change that fold's training-only tuning.
altered_args <- cv_args
altered_args$genotypes[heldout, ] <- 2
altered_args$outer_fold_id <- cv$outer_fold_id
altered <- do.call(cross_validate_transfer, altered_args)
assert_equal(cv$tuning[[1L]]$mean_mse, altered$tuning[[1L]]$mean_mse)

# Lasso and a no-treatment subgroup (constant treatments) remain usable.
cv_args$penalty <- "lasso"
cv_args$data$treatment <- cv_args$data$chemo <- 0
lasso <- do.call(cross_validate_transfer, cv_args)
stopifnot(all(is.finite(lasso$per_fold$r2_transfer_min)))
fitted <- fit_transfer_model(transform$train_y, dat$prs, sim$genotypes,
                             treatments = as.matrix(dat[c("treatment", "chemo")]), lambda = 0.1)
pred <- predict_transfer_model(fitted, dat$prs, sim$genotypes,
                                treatments = as.matrix(dat[c("treatment", "chemo")]))
stopifnot(length(pred) == nrow(dat), all(is.finite(pred)))
cat("All synthetic-data checks passed.\n")
