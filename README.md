# PRS variance decomposition and transfer learning

A small R workflow based on `mixedModel.R` and the transfer learning scripts:
adjust a phenotype, estimate genetic variance components, and adapt an existing
PRS with ridge or lasso SNP effects.

The [SNP lists](snp_lists/README.md) include variant IDs, alleles, and source-GWAS
weights for Height, BMI, and WHR. No participant-level data are included.

## Run the simulation

Install R packages once, then run from this folder:

```sh
Rscript --vanilla -e 'install.packages(c("gaston", "glmnet"), repos="https://cloud.r-project.org")'
Rscript --vanilla examples/run_simulation.R
```

The example generates 240 subjects and 80 independent SNPs with invented PRS
weights. It runs the full workflow, prints results, and writes five summary CSVs
to `outputs/`. It uses 3 outer and 3 inner folds for speed; the CV function
defaults to 5 of each. Tested with R 4.4.0, gaston 1.6, and glmnet 4.1-8.

| File | Purpose |
| --- | --- |
| `R/preprocessing.R` | Covariate adjustment, inverse-normal transformation, imputation, scaling |
| `R/mixed_model.R` | PRS regressions, genetic kernels, REML estimates, variance-component p-values |
| `R/transfer_learning.R` | Ridge/lasso adaptation and nested cross-validation |
| `R/simulation.R` | Synthetic data generation |
| `examples/run_simulation.R` | Complete example |
| `snp_lists/` | Nine SNP tables and matching ID-only lists |

## Use your own data

Prepare a data frame `dat` and a numeric subjects-by-SNPs matrix `G`, with subjects
in exactly the same order. For example, if genotype row names are sample IDs:

```r
G <- G[match(dat$id, rownames(G)), , drop = FALSE]
```

IDs must be unique and matched. Outcome, PRS, treatment, group, and adjustment
columns must be complete. Genotypes may contain `NA`; missing dosages are
mean-imputed and constant/all-missing SNPs are dropped. Each training group needs
enough subjects for the adjustment regression. Use independent subjects and
non-collinear covariates; related subjects or repeated measures need a different
fold design. These functions assume prepared inputs rather than providing a
general data-validation layer.

Load the functions and fit the mixed model:

```r
source("R/preprocessing.R")
source("R/mixed_model.R")
source("R/transfer_learning.R")

adjusted <- fit_outcome_transform(
  dat, outcome = "outcome", adjustment = ~ age + I(age^2) + PC1 + PC2, group = "sex"
)
compare_prs_models(adjusted$train_y, dat$prs, dat$treatment)

X <- cbind(Intercept = 1, PRS = dat$prs, Treatment = dat$treatment)
mixed <- fit_variance_components(
  adjusted$train_y, X, build_kernels(G, dat$treatment)
)
mixed$components
mixed$tests
```

Replace the column names and adjustment formula for your trait. Use `group = NULL`
for pooled adjustment. The mixed model is `y = X beta + u_G + u_GxT + error`.
The genetic kernel uses standardized dosages. The interaction kernel multiplies
raw mean-imputed dosages by treatment, then standardizes the products, preserving
the order used in the research scripts. Treatment must vary for this two-kernel fit.

`components` reports each variance coefficient and its share of
`tau_G + tau_GxT + sigma2`; these shares exclude fixed-effect variation.
`tests` gives approximate zero-variance likelihood-ratio p-values using the
half-point-mass/half-chi-square(1) reference. This approximation can be unreliable
with small samples, similar kernels, or other components on the boundary.

Run transfer learning from the raw phenotype and unimputed genotypes:

```r
cv <- cross_validate_transfer(
  dat, G, outcome = "outcome", prs = "prs",
  treatment_cols = c("treatment", "chemo"),
  adjustment = ~ age + I(age^2) + PC1 + PC2, group = "sex",
  strata = interaction(dat$treatment, dat$chemo), penalty = "ridge", seed = 42
)
cv$per_fold
cv$summary
```

Use `penalty = "lasso"` for L1 adaptation. PRS and treatment coefficients are
unpenalized; SNP coefficients are penalized. Phenotype adjustment, imputation,
and scaling are learned separately within every training split. The source PRS
must be constructed without using target test outcomes.

The inner folds choose minimum-MSE and one-standard-error lambdas; the outer
folds measure treatment-only, PRS, and adapted prediction R-squared. The one-SE
rule chooses a more regularized model and is not a confidence interval.
`summary` gives mean and SD across outer folds; `predictions` holds held-out
predictions. R-squared and its improvements can be negative.
The candidate grid can be set with `lambda`; `outer_fold_id` can reuse splits.

The simulation's `truth` applies to `y_model` with complete genotypes, before
covariate adjustment and rank transformation. A single run illustrates the
process and does not establish performance or guarantee improvement.

## Upload to GitHub

Upload this release folder, including `snp_lists/`. `MANIFEST.txt` lists the
intended files, and `.gitignore` excludes generated outputs and other unlisted
files. Choose your repository name, author/citation details, and license.
