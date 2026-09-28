#!/usr/bin/env Rscript
# Run from the repository root: Rscript --vanilla examples/run_simulation.R
# Optional first argument is the output directory. Everything generated is synthetic.
for (file in c("preprocessing", "mixed_model", "transfer_learning", "simulation")) {
  source(file.path("R", paste0(file, ".R")))
}
.require_package("gaston")
.require_package("glmnet")
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1L] else "outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sim <- simulate_example_data()
dat <- sim$data
cat("Synthetic example:", nrow(dat), "subjects and", ncol(sim$genotypes), "SNPs\n")

# A. Known-covariance demonstration: use the generated Gaussian outcome before
# nuisance effects or rank transformation, and the complete simulated genotypes.
# The truth is meaningful on THIS scale.
kernels <- build_kernels(sim$complete_genotypes, dat$treatment)
X <- cbind(Intercept = 1, PRS = dat$prs, Treatment = dat$treatment)
known_fit <- fit_variance_components(sim$y_model, X, kernels)
known_results <- known_fit$components
known_results$generating_coefficient <- unname(sim$truth[known_results$component])
cat("\nGenerating variance coefficients and estimates (one replicate):\n")
print(known_results)

# B. Analysis workflow matching the original phenotype preparation.
transform <- fit_outcome_transform(dat, "outcome", ~ age + I(age^2) + PC1 + PC2, "sex")
fixed <- compare_prs_models(transform$train_y, dat$prs, dat$treatment)
analysis_fit <- fit_variance_components(transform$train_y, X,
                                       build_kernels(sim$genotypes, dat$treatment))
tests <- do.call(rbind, lapply(c("GxT", "G"), function(component) {
  test_variance_component(analysis_fit, component)
}))
cat("\nTransformed-outcome variance decomposition:\n")
print(analysis_fit$components)
print(tests)

# C. Prediction: raw outcomes and unimputed genotypes enter nested CV.
# The default example uses 3 x 3 folds for speed; functions default to 5 x 5.
cv <- cross_validate_transfer(dat, sim$genotypes, outcome = "outcome", prs = "prs",
  treatment_cols = c("treatment", "chemo"), adjustment = ~ age + I(age^2) + PC1 + PC2,
  group = "sex", strata = interaction(dat$treatment, dat$chemo),
  outer_folds = 3L, inner_folds = 3L, penalty = "ridge", seed = 42L)
cat("\nHeld-out prediction results (mean and SD across outer folds):\n")
print(cv$summary)
if (any(cv$per_fold$lambda_min_at_edge)) {
  message("A selected lambda is at the grid edge; a larger study should examine its tuning range.")
}

# Optional, slower profile example (uncomment to run). Expand the upper endpoint
# if upper_grid_edge_accepted is TRUE; this is a grid approximation.
# profile <- profile_variance_component(analysis_fit, "GxT", seq(0, 0.8, length.out = 21))
# print(profile$interval)

write.csv(known_results, file.path(output_dir, "known_variance_example.csv"), row.names = FALSE)
write.csv(fixed, file.path(output_dir, "fixed_effect_comparisons.csv"), row.names = FALSE)
write.csv(analysis_fit$components, file.path(output_dir, "variance_components.csv"), row.names = FALSE)
write.csv(tests, file.path(output_dir, "variance_tests.csv"), row.names = FALSE)
write.csv(cv$per_fold, file.path(output_dir, "transfer_per_fold.csv"), row.names = FALSE)
write.csv(cv$summary, file.path(output_dir, "transfer_summary.csv"), row.names = FALSE)

grDevices::pdf(file.path(output_dir, "simulation_diagnostics.pdf"), width = 10, height = 4.5)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
barplot(rbind(known_results$generating_coefficient, known_results$variance),
        beside = TRUE, names.arg = known_results$component, col = c("grey75", "steelblue"),
        ylab = "Variance coefficient", main = "One synthetic replicate")
legend("topright", c("Generating", "Estimated"), fill = c("grey75", "steelblue"), bty = "n")
values <- as.matrix(cv$per_fold[, c("r2_treatment", "r2_prs", "r2_transfer_min")])
matplot(seq_len(nrow(values)), values, type = "b", pch = 1:3, lty = 1,
        col = c("grey40", "steelblue", "darkorange"), xaxt = "n",
        xlab = "Outer test fold", ylab = expression(R^2), main = "Held-out prediction")
axis(1, seq_len(nrow(values)))
abline(h = 0, lty = 3, col = "grey70")
legend("bottomleft", c("Treatment", "+ PRS", "+ SNP adaptation"),
       col = c("grey40", "steelblue", "darkorange"), lty = 1, pch = 1:3, bty = "n")
invisible(grDevices::dev.off())
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("\nSynthetic summaries and plot written to", output_dir, "\n")
