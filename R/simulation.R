# Independent synthetic SNPs and invented external PRS weights; no study data.
simulate_example_data <- function(n = 240L, p = 80L, seed = 2026L,
                                   tau_g = 0.2, tau_gxt = 0.15, sigma2 = 0.65,
                                   missing_rate = 0.01) {
  set.seed(seed)
  maf <- runif(p, 0.1, 0.5)
  genotypes <- sapply(maf, function(a) rbinom(n, 2, a))
  colnames(genotypes) <- paste0("synthetic_snp_", seq_len(p))
  sex <- sample(rep(c(0, 1), length.out = n))
  treatment <- sample(rep(c(0, 0, 1), length.out = n))
  chemo <- rbinom(n, 1, 0.4)
  age <- runif(n, 20, 60)
  pc1 <- rnorm(n)
  pc2 <- rnorm(n)
  z <- scale(genotypes)
  w <- scale(genotypes * treatment)
  prs <- as.numeric(scale(z %*% (rnorm(p) / sqrt(p))))
  genetic <- z %*% rnorm(p, sd = sqrt(tau_g / p))
  interaction <- w %*% rnorm(p, sd = sqrt(tau_gxt / p))
  y_model <- as.numeric(0.5 * prs + 0.25 * treatment + genetic + interaction +
                          rnorm(n, sd = sqrt(sigma2)))
  outcome <- 10 + 0.8 * sex + 0.03 * (age - 40) + 0.001 * (age - 40)^2 +
    0.12 * pc1 - 0.1 * pc2 + y_model
  data <- data.frame(id = paste0("synthetic_", seq_len(n)), outcome, sex, age,
                     PC1 = pc1, PC2 = pc2, prs, treatment, chemo)
  rownames(genotypes) <- data$id
  observed <- genotypes
  observed[matrix(runif(n * p) < missing_rate, n, p)] <- NA_real_
  list(data = data, genotypes = observed, complete_genotypes = genotypes,
       y_model = y_model, truth = c(GxT = tau_gxt, G = tau_g, Residual = sigma2))
}
