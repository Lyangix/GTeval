# Run from the repository root: Rscript --vanilla examples/run_simulation.R
for (name in c("preprocessing", "mixed_model", "transfer_learning", "simulation")) {
  source(file.path("R", paste0(name, ".R")))
}
sim <- simulate_example_data()
dat <- sim$data

# 1. Adjust the phenotype within sex and inverse-normal transform residuals.
adjusted <- fit_outcome_transform(dat, "outcome", ~ age + I(age^2) + PC1 + PC2, "sex")
fixed <- compare_prs_models(adjusted$train_y, dat$prs, dat$treatment)

# 2. Estimate genetic, genetic-by-treatment, and residual variance.
X <- cbind(Intercept = 1, PRS = dat$prs, Treatment = dat$treatment)
kernels <- build_kernels(sim$genotypes, dat$treatment)
mixed <- fit_variance_components(adjusted$train_y, X, kernels)
print(fixed)
print(mixed$components)
print(mixed$tests)

# 3. Adapt the PRS with penalized SNP effects using nested cross-validation.
# Raw outcomes and unimputed genotypes enter CV. Use penalty = "lasso" for L1.
cv <- cross_validate_transfer(dat, sim$genotypes, outcome = "outcome", prs = "prs",
  treatment_cols = c("treatment", "chemo"), adjustment = ~ age + I(age^2) + PC1 + PC2,
  group = "sex", strata = interaction(dat$treatment, dat$chemo),
  outer_folds = 3L, inner_folds = 3L, penalty = "ridge", seed = 42L)
print(cv$summary)

dir.create("outputs", showWarnings = FALSE)
write.csv(fixed, "outputs/fixed_effect_comparisons.csv", row.names = FALSE)
write.csv(mixed$components, "outputs/variance_components.csv", row.names = FALSE)
write.csv(mixed$tests, "outputs/variance_tests.csv", row.names = FALSE)
write.csv(cv$per_fold, "outputs/transfer_per_fold.csv", row.names = FALSE)
write.csv(cv$summary, "outputs/transfer_summary.csv", row.names = FALSE)
