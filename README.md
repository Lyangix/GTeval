# PRS variance decomposition and transfer learning

Reusable R functions for studying genetic variation beyond an existing polygenic
risk score (PRS), and adapting that score to a target cohort with ridge or lasso
regression. This is a portable extraction of the `mixedModel.R` and
`transfer_learning_*_CV.R` workflows, with explicit inputs and synthetic examples.

**No participant data, real genotypes, real PRS weights, or study results are included.**
The example generates all inputs locally from a random seed. Study data are not
distributed through this repository.

## Run the example

Use R with the `gaston` and `glmnet` packages. The release was tested with R 4.4.0,
gaston 1.6, and glmnet 4.1-8. From this repository's root:

```sh
Rscript --vanilla -e 'install.packages(c("gaston", "glmnet"), repos="https://cloud.r-project.org")'
Rscript --vanilla examples/run_simulation.R
Rscript --vanilla tests/run_tests.R
```

Installation needs internet access; running the example needs no downloaded data.
The default example has 240 simulated subjects and 80 independent SNPs, and uses
3 outer and 3 inner folds. It writes synthetic summary tables, a PDF plot, and
`sessionInfo.txt` to `outputs/`, which is excluded from version control. It also
prints the estimates and prediction results. To choose a different output folder:

```sh
Rscript --vanilla examples/run_simulation.R /tmp/prs-demo
```

These are sourceable R functions, not an installable R package. Load them with:

```r
source("R/preprocessing.R")
source("R/mixed_model.R")
source("R/transfer_learning.R")
source("R/simulation.R")
```

## What the models do

The mixed model fits

\[
y = X\beta + u_G + u_{G\times T} + \epsilon, \qquad
\operatorname{Var}(y\mid X) = \tau_G K_G + \tau_{G\times T}K_{G\times T} + \sigma_e^2 I.
\]

`X` contains the intercept, selected PRS, and treatment. The functions construct
the genetic and genetic-by-treatment kernels, fit variance components by REML,
compare full and restricted models, and optionally profile a variance coefficient
over a user-specified grid. Preliminary linear regressions also compare treatment,
treatment plus PRS, and treatment plus PRS plus PRS-by-treatment.

Transfer learning fits

\[
y = \alpha_0 + \alpha_1\mathrm{PRS}^{(G)} + C\gamma + Z\delta + \epsilon.
\]

The external PRS and treatment covariates are unpenalized. SNP adaptation
coefficients `delta` receive a ridge or lasso penalty. Nested cross-validation
selects the penalty and compares treatment-only, PRS, and adapted predictions on
the same held-out subjects. Every fitted preprocessing step uses training rows.

## Functions and files

| File | Main functions | Purpose |
| --- | --- | --- |
| `R/preprocessing.R` | `align_genotypes()`, `fit_outcome_transform()`, `transform_outcome()` | Align samples; fit group-specific covariate adjustment and rank transformation |
| `R/mixed_model.R` | `compare_prs_models()`, `build_kernels()`, `fit_variance_components()` | Fixed-effect comparisons and genetic variance decomposition |
| `R/mixed_model.R` | `test_variance_component()`, `profile_variance_component()` | Named-component likelihood tests and approximate grid confidence sets |
| `R/transfer_learning.R` | `cross_validate_transfer()` | Nested CV, minimum-error and one-SE penalty choices, held-out R-squared |
| `R/transfer_learning.R` | `fit_transfer_model()`, `predict_transfer_model()` | Fit a selected model and predict in new samples |
| `R/simulation.R` | `simulate_example_data()` | Generate an entirely synthetic demonstration |

See [the data and usage guide](docs/USAGE.md) for input formats and examples,
[the statistical notes](docs/METHODS.md) for definitions and differences from the
research scripts, and [the publishing guide](docs/PUBLISHING.md) for creating a
code-only GitHub repository.

## Interpreting the example

The example first fits a Gaussian phenotype whose generating variance coefficients
are known, then demonstrates the original residual-adjustment and inverse-normal
workflow. The generating coefficients apply to the first phenotype only; the
rank transformation changes the scale. A single replicate demonstrates execution,
not estimator bias, coverage, power, or guaranteed prediction improvement.

Variance shares describe random variation conditional on the fixed effects.
They are different from the PRS contribution to R-squared and from held-out
prediction R-squared. Held-out R-squared and improvement can be negative.

The implementation uses dense matrices and is intended as a readable reference
and a small-data example. Kernel storage grows quadratically with sample size;
large studies need an appropriate memory budget or a separately validated
blockwise implementation.

The model backends are documented in the
[gaston reference](https://search.r-project.org/CRAN/refmans/gaston/html/lmm.aireml.html)
and the [glmnet guide](https://glmnet.stanford.edu/articles/glmnet.html).
