# Statistical definitions and extraction notes

## Phenotype preparation

The original height analysis separately regresses phenotype on age, age squared,
and 20 PCs within each sex. Within a training group of size `n_g`, residual `r_i`
is transformed to `qnorm((rank(r_i) - 0.5) / n_g)`, using average ranks for ties.
The public functions let the caller supply the outcome, adjustment formula,
and grouping variable. The example uses two synthetic PCs to keep it small.

In cross-validation, each inner and outer training sample fits its own adjustment
regression and residual empirical CDF. A held-out outcome is adjusted with those
training coefficients and mapped through the training CDF. Probabilities are
clipped to `[epsilon, 1 - epsilon]`, where
`epsilon = 1 / (2 * (total_training_n + 1))`, matching the source CV scripts.
Training values use midranks; evaluation values use the training ECDF. They are
intentionally different mappings at observed residual values. Held-out outcomes
are used to define evaluation targets, never to train the predictor or choose lambda.

## Kernel construction and interpretation

Let `D` be the mean-imputed dosage matrix. Standardize each retained SNP column
with its sample mean and sample standard deviation to obtain `Z`. Separately
form `D * T` by multiplying each subject's dosage row by the raw treatment value,
and standardize those product columns to obtain `W`. Then

```
K_G   = Z Z' / number_of_retained_G_columns
K_GxT = W W' / number_of_retained_GxT_columns
```

This preserves the product-before-standardization convention of `mixedModel.R`.
It is generally different from multiplying the already standardized genetic
kernel by a treatment outer product. After centering the products, untreated
rows of `W` need not be zero. Treatment coding therefore matters: specify and
report whether the variable is binary, continuous, centered, or rescaled.

All-missing and constant dosage columns are excluded. Zero-variance interaction
columns are also excluded. Each kernel is divided by its own retained count, and
the retained indices and counts are stored as attributes. With sample-standardized
columns, each kernel has mean diagonal `(n - 1) / n`.

The model assumes independent random effects with covariances `tau_G * K_G`,
`tau_GxT * K_GxT`, and `sigma2 * I`. It does not estimate covariance between the
two random effects. The PRS and treatment are fixed effects. The two random
components describe residual covariance after accounting for those fixed effects;
they are not orthogonal observed scores or causal treatment effects.

`coefficient_share` divides each fitted coefficient by
`tau_G + tau_GxT + sigma2`. `average_variance` multiplies each coefficient by the
mean diagonal of its kernel (one for the residual), and `conditional_share`
normalizes those contributions. The latter also handles differently scaled
user-supplied kernels. Neither includes fixed-effect variation in its denominator.
Do not add these shares to an ordinary-regression PRS R-squared and interpret
the sum as a total phenotype decomposition.

## Fitting and inference

The wrapper uses constrained AI-REML in
[gaston](https://search.r-project.org/CRAN/refmans/gaston/html/lmm.aireml.html).
Variance coefficients may equal zero. Fits that reach the iteration limit or
return nonfinite estimates raise errors. The returned iteration count and
gradient norm support further review; these checks cannot guarantee identifiability
or a global optimum. Highly similar kernels can make individual components hard
to distinguish. Custom kernels must be symmetric and positive semidefinite.

For a zero null, the named kernel is removed and the remaining model is fitted
with identical `y` and `X`. For a positive null, the named coefficient is fixed
and the remaining nonnegative coefficients are optimized against the
[restricted likelihood](https://search.r-project.org/CRAN/refmans/gaston/html/lmm.restricted.likelihood.html).
The full fit is reused across a profile grid. A materially higher restricted
likelihood raises an error instead of producing a misleading p-value.

For an interior coefficient, the likelihood ratio is compared with chi-square(1).
For a single variance tested at zero, the reference is the approximate mixture
of a point mass at zero and chi-square(1), with weights one half each. A zero
statistic gets p = 1, using the inclusive upper tail at the point mass. This
mixture assumes an identifiable component and suitable regularity, including
interior nuisance parameters. Other variance components on the boundary, strong
kernel similarity, and small samples can invalidate the approximation. It is not
an exact finite-sample test, and the helper does not offer a free-covariance model.
Use a separately designed null simulation or bootstrap when calibration is needed.
For background on boundary asymptotics and nuisance-parameter qualifications, see
[Self and Liang (1987)](https://doi.org/10.1080/01621459.1987.10478472).

The grid profile inverts these pointwise tests. Its accepted range is approximate,
and an accepted grid endpoint does not establish the true confidence limit. No
multiple-testing correction is built in for analyses over many outcomes, treatments,
or PRS thresholds. The preliminary PRS and PRS-by-treatment comparisons instead
use nested linear-model F tests and incremental in-sample R-squared.

## Transfer learning and evaluation

The source-domain information enters as an externally computed PRS. The target
cohort estimates an intercept, PRS coefficient, treatment effects, and SNP-specific
adaptation coefficients. Only the SNP coefficients receive a penalty. Ridge and
lasso use the [glmnet penalty-factor interface](https://glmnet.stanford.edu/articles/glmnet.html).
This implements PRS-based adaptation, not joint retraining on source and target
individual-level datasets. It contains no SNP-by-treatment prediction terms.

Every training split learns genotype imputation means and column scales. Fixed
predictors and retained genotypes are standardized with training sample standard
deviations; `glmnet` receives `standardize = FALSE`. Constant treatments are dropped;
collinear remaining fixed columns and a constant PRS cause a clear error. A
prespecified lambda grid is shared across folds. As in glmnet generally, penalty
factors are internally rescaled; interpret lambda with this exact design and scaling.

Each outer training sample runs a full inner CV, including re-estimation of outcome
transformation and genotype preprocessing. Mean validation MSE is averaged equally
across inner folds. `lambda_min` minimizes it. `lambda_1se` is the largest lambda
within one estimated standard error of that minimum. Outer-test subjects contribute
only to evaluation after selection and refitting on the outer training sample.
Fold-specific adjustment means transformed targets can differ between outer folds;
the reported summary averages fold metrics rather than claiming a single global
transformed target or pooling training and test results.

For each evaluation set, `R2 = 1 - SSE / sum((y - mean(y))^2)`; it is undefined
for a constant outcome and can be negative out of sample. An intercept-only
prediction uses the training mean, so its test R-squared need not be zero. The
reported differences are incremental R-squared, sometimes labeled "partial R2"
in the source scripts. They are not the conventional partial coefficient
`(SSE_reduced - SSE_full) / SSE_reduced` and are not relative percentage improvements.

## Deliberate differences from the research scripts

| Research-script detail | Public extraction |
| --- | --- |
| Institution-specific paths, merges, and SLURM choices | Explicit data, matrices, formulas, PRS, and treatment inputs |
| Hard-coded `1:3220` kernel slice in `mixedModel.R` | Dimensions come from input rows |
| Final `SCORE_random` assignment overrides the selected PRS in `mixedModel.R` | Caller chooses `X`; no implicit override |
| `tau / (tau + sigma2)` in `mixedModel.R` omits the other component from each denominator | Both coefficient shares and mean-diagonal conditional shares are returned |
| Global genotype imputation in the transfer CV scripts | Imputation and scaling fitted within every training split |
| Lambda grid derived from the first inner training fold | Fixed, outcome-independent candidate grid |
| Repeated single-lambda fits | One path per split, using the requested grid |
| Implicit glmnet predictor standardization | Explicit sample-SD training standardization and `standardize = FALSE` |
| Constant/missing SNPs and constant subgroup treatments can cause numerical problems | Drop unusable columns and validate the remaining design |
| Some scripts divide a filtered interaction kernel by the original SNP count | Each kernel uses its retained SNP count, recorded in metadata |
| Boundary p-value expression returns 0.5 at exactly zero LRT | Inclusive tail p = 1 at zero; negative numerical noise is handled explicitly |
| Full fit repeated at every confidence-grid value | Full fit cached; restricted fits use component names |
| Separate trait/subgroup script copies and mixed train/test averages | Parameterized functions; held-out summaries only |

These differences are intentional. This release is a reproducible method
demonstration, not a claim of byte-for-byte reproduction of previous study tables.
The source scripts in the parent working directory are not modified.

## Simulation scope

Independent binomial dosages, invented source weights, Gaussian SNP adaptation,
and Gaussian SNP-by-treatment coefficients generate the example. Genotypes are
drawn independently across SNPs: there is no LD, ancestry structure, or relatedness.
For the known-covariance fit, the generating coefficients correspond to the exact
complete-genotype kernels used in fitting. Missingness is then introduced into
a separate observed dosage matrix for demonstrating preprocessing and prediction.
The rank-transformed analysis does not retain those exact generating coefficients.
Assessing operating characteristics requires repeated, appropriately designed
simulations beyond this one-run example.
