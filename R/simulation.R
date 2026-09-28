# All draws are synthetic. No cohort data, study estimates, or real SNP weights
# are used. This small example assumes independent SNPs and unrelated subjects.
simulate_example_data <- function(n = 240L, p = 80L, seed = 2026L,
                                   tau_g = 0.2, tau_gxt = 0.15, sigma2 = 0.65,
                                   missing_rate = 0.01) {
  if (length(n) != 1L || length(p) != 1L || !is.finite(n) || !is.finite(p) ||
      n < 40L || p < 5L || n != as.integer(n) || p != as.integer(p) ||
      any(!is.finite(c(tau_g, tau_gxt, sigma2, missing_rate))) ||
      tau_g < 0 || tau_gxt < 0 || sigma2 <= 0 || missing_rate < 0 || missing_rate >= 1) {
    stop("Invalid simulation dimensions, variances, or missingness rate.")
  }
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  maf <- stats::runif(p, 0.1, 0.5)
  genotypes <- vapply(maf, function(a) stats::rbinom(n, 2, a), numeric(n))
  colnames(genotypes) <- paste0("synthetic_snp_", seq_len(p))
  # Resample any rare constant columns so the known covariance is well-defined.
  for (j in seq_len(p)) {
    while (stats::sd(genotypes[, j]) == 0) genotypes[, j] <- stats::rbinom(n, 2, maf[j])
  }
  sex <- sample(rep(c(0, 1), length.out = n))
  treatment <- sample(rep(c(0, 0, 1), length.out = n))
  chemo <- stats::rbinom(n, 1, 0.4)
  age <- stats::runif(n, 20, 60)
  pc1 <- stats::rnorm(n)
  pc2 <- stats::rnorm(n)
  z <- scale(genotypes)
  product <- genotypes * treatment
  keep <- which(apply(product, 2L, stats::sd) > 0)
  w <- scale(product[, keep, drop = FALSE])
  # The source weights are invented without using target phenotypes.
  source_weights <- stats::rnorm(p) / sqrt(p)
  prs <- as.numeric(scale(z %*% source_weights))
  delta <- stats::rnorm(p, sd = sqrt(tau_g / p))
  interaction_effect <- stats::rnorm(ncol(w), sd = sqrt(tau_gxt / ncol(w)))
  y_model <- 0.5 * prs + 0.25 * treatment +
    as.numeric(z %*% delta + w %*% interaction_effect) + stats::rnorm(n, sd = sqrt(sigma2))
  outcome <- 10 + 0.8 * sex + 0.03 * (age - 40) + 0.001 * (age - 40)^2 +
    0.12 * pc1 - 0.1 * pc2 + y_model
  observed <- genotypes
  observed[matrix(stats::runif(n * p) < missing_rate, n, p)] <- NA_real_
  data <- data.frame(id = paste0("synthetic_", seq_len(n)), outcome = outcome,
                      sex = sex, age = age, PC1 = pc1, PC2 = pc2,
                      prs = prs, treatment = treatment, chemo = chemo)
  rownames(genotypes) <- rownames(observed) <- data$id
  list(data = data, genotypes = observed, complete_genotypes = genotypes,
       y_model = y_model, source_weights = source_weights,
       truth = c(GxT = tau_gxt, G = tau_g, Residual = sigma2), seed = seed)
}
