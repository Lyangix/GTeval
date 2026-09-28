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

dir.create("outputs", showWarnings = FALSE)
write.csv(fixed, "outputs/fixed_effect_comparisons.csv", row.names = FALSE)
write.csv(mixed$components, "outputs/variance_components.csv", row.names = FALSE)
write.csv(mixed$tests, "outputs/variance_tests.csv", row.names = FALSE)
