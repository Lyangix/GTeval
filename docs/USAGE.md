# Using your own data locally

Run all examples from the repository root after sourcing the four R files in
the README. Input files stay outside the public repository. The functions accept
R objects; they do not load any institution-specific files or require a scheduler.

## Input contract

| Input | Required structure |
| --- | --- |
| Phenotype table `dat` | One row per subject; finite numeric outcome, PRS, adjustment covariates, and numeric treatment columns |
| Grouping column | Optional complete column such as sex; each training group must have enough observations for its adjustment regression |
| Genotype matrix `G` | Numeric subjects-by-SNPs dosage matrix, normally 0/1/2 or imputed dosages; `NA` allowed |
| IDs | Unique, nonmissing sample IDs; phenotype and genotype tables must identify the same analysis subjects |
| External PRS | Computed using source weights selected without the target test outcomes; consistent allele coding is the caller's responsibility |

The functions reject incomplete phenotypes and covariates. Define the analysis
sample first, removing incomplete rows consistently from the phenotype and
genotype data. A missing treatment value is not automatically an untreated value.
Do not impute or standardize genotypes before calling nested CV: it learns those
operations separately in each training set. Constant or all-missing SNPs in a
training set are excluded; missing retained dosages use that training set's means.

Here is an executable input-format example using simulated objects:

```r
sim <- simulate_example_data()
dat <- sim$data
G <- align_genotypes(sim$genotypes, genotype_ids = rownames(sim$genotypes),
                     sample_ids = dat$id)
stopifnot(identical(rownames(G), dat$id))
```

For real inputs, create `dat` and `G` from files in a private location, then use
the same calls below. `align_genotypes()` reorders by IDs and stops for duplicates
or unmatched subjects. It does not silently drop unmatched rows. The modeling
functions consume matrices in the supplied order; use this alignment step first.
When subsetting, subset `dat` and the aligned `G` with the same row indices.

## Fixed-effect and mixed-model analyses

```r
adjustment <- ~ age + I(age^2) + PC1 + PC2
transform <- fit_outcome_transform(dat, outcome = "outcome",
                                   adjustment = adjustment, group = "sex")
y <- transform$train_y

# Choose the PRS and treatment explicitly for this analysis.
X <- cbind(Intercept = 1, PRS = dat$prs, Treatment = dat$treatment)
fixed_comparisons <- compare_prs_models(y, dat$prs, dat$treatment)
kernels <- build_kernels(G, dat$treatment)
model <- fit_variance_components(y, X, kernels)
model$components
model$beta
test_variance_component(model, "GxT", value = 0)
test_variance_component(model, "G", value = 0)

# Optional slower profile; this tests variance COEFFICIENTS, not proportions.
profile <- profile_variance_component(model, "GxT", seq(0, 1, length.out = 21))
profile$interval
profile$upper_grid_edge_accepted
```

Replace the outcome name, adjustment formula, treatment vector, and PRS column
to analyze height, BMI, WHR, or another continuous phenotype. Add the desired PCs
explicitly to the formula. Pass `group = NULL` for pooled adjustment. The design
matrix must include an intercept if one is wanted; `X` is never silently changed.
The same sample and fixed effects are used in full and restricted REML fits.

`model$components` contains variance coefficients, their normalized shares, the
average diagonal variance contributed by each kernel, and shares based on those
contributions. `beta` contains estimated fixed effects. Iteration count and final
gradient norm are retained in `model`; an iteration-limit failure stops the fit.
Custom named positive semidefinite kernels can also be supplied.

`profile$table` contains each tested value and p-value; `accepted` is the accepted
set of grid values. `interval` is their range, not an interpolated exact interval.
If the upper grid edge is accepted, expand it. If a positive lower grid edge is
accepted, extend downward. If `disconnected` is true, inspect the full confidence
set and optimization rather than reporting its range as a single interval.

## Nested cross-validation for transfer learning

```r
cv <- cross_validate_transfer(
  data = dat, genotypes = G, outcome = "outcome", prs = "prs",
  treatment_cols = c("treatment", "chemo"),
  adjustment = ~ age + I(age^2) + PC1 + PC2, group = "sex",
  strata = interaction(dat$treatment, dat$chemo),
  outer_folds = 5, inner_folds = 5, penalty = "ridge", seed = 42
)
cv$summary
cv$per_fold
```

Use `penalty = "lasso"` for L1 adaptation. Ridge versus lasso is a prespecified
choice here: selecting whichever has better outer-test performance requires an
additional independent evaluation. To compare the methods on identical splits,
reuse `outer_fold_id = cv$outer_fold_id` and the same seed.

The default grid is `10^seq(2, -4, length.out = 50)`. You can supply a positive
numeric `lambda` vector; it is sorted from strongest to weakest penalty. The
grid is fixed independently of target outcomes. Examine the saved inner tuning
curves, particularly `lambda_min_at_edge`. Any grid adaptation must use training
data rather than outer test results.

`cv$per_fold` reports held-out R-squared for treatment-only, PRS plus treatment,
and adaptation under minimum-error and one-SE choices. The PRS delta compares
PRS plus treatment to treatment-only. Adaptation deltas compare adaptation to
PRS plus treatment. `summary` gives the unweighted mean and SD of these metrics
across outer folds. The SD describes fold variation, not a confidence interval.
There is no average that mixes training and test performance.

`cv$predictions` stores row indices, fold numbers, transformed evaluation outcomes,
and predictions. `cv$tuning` stores inner-fold IDs and errors. With real data these
objects contain participant-level or study-derived information; keep them private.
The example writes only summaries of its simulated data.

For an untreated subgroup, subset rows first and pass `treatment_cols = character()`.
Constant treatment columns are otherwise dropped within each training fold.
Collinear remaining fixed predictors raise an error. Every inner training group
must have enough rows to fit the adjustment formula. Reduce fold counts or adjust
the prespecified design when samples are too small. The automatic splitter assumes
independent subjects. Related subjects, sites, or repeated measures require
group-aware inner and outer splitting; this helper does not implement that design.

## Fitting a selected model for new samples

The lower-level functions support a final model once a penalty has been selected
using only its training data:

```r
# For this mechanics example, lambda is prespecified, not claimed to be optimal.
outcome_map <- fit_outcome_transform(dat, "outcome", adjustment, "sex")
final_model <- fit_transfer_model(
  y = outcome_map$train_y, prs = dat$prs, genotypes = G,
  treatments = as.matrix(dat[c("treatment", "chemo")]),
  penalty = "ridge", lambda = 0.1
)

# In actual use these predictors come from independent samples, in the same
# SNP and treatment column order. No outcomes are needed to make predictions.
predicted <- predict_transfer_model(
  final_model, prs = dat$prs, genotypes = G,
  treatments = as.matrix(dat[c("treatment", "chemo")])
)
```

The last call demonstrates the API on the training inputs and is not a validation
result. New-sample predictions are on the adjusted, inverse-normal scale. If new
outcomes are available for evaluation, transform them with
`transform_outcome(outcome_map, new_dat)`. The current functions do not invert this
mapping to original clinical units. Final-model tuning on the entire development
sample is a separate training step from reporting nested CV performance.
