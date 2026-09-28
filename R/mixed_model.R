# Preliminary fixed-effect comparisons: treatment, PRS, and PRS-by-treatment.
compare_prs_models <- function(y, prs, treatment) {
  data <- data.frame(y, prs, treatment)
  base <- lm(y ~ treatment, data)
  score <- lm(y ~ prs + treatment, data)
  interaction <- lm(y ~ prs * treatment, data)
  r2 <- sapply(list(base, score, interaction), function(fit) summary(fit)$r.squared)
  data.frame(model = c("Treatment", "Treatment + PRS", "Treatment + PRS + PRS:T"),
             r_squared = r2, delta_r_squared = c(NA, diff(r2)),
             p_value = c(NA, anova(base, score)[["Pr(>F)"]][2],
                         anova(score, interaction)[["Pr(>F)"]][2]))
}

# Multiply raw mean-imputed dosages by treatment BEFORE standardizing products.
build_kernels <- function(genotypes, treatment) {
  prep <- fit_numeric_preprocessor(genotypes)
  dosage <- genotypes[, prep$keep, drop = FALSE]
  for (j in seq_along(prep$keep)) {
    dosage[is.na(dosage[, j]), j] <- prep$center[prep$keep[j]]
  }
  z <- apply_numeric_preprocessor(prep, genotypes)
  product <- dosage * treatment
  w <- apply_numeric_preprocessor(fit_numeric_preprocessor(product), product)
  list(GxT = tcrossprod(w) / ncol(w), G = tcrossprod(z) / ncol(z))
}

# Fit both kernels, then omit each in turn for an approximate zero-variance LRT.
fit_variance_components <- function(y, X, kernels, max_iter = 200L) {
  fit_model <- function(K) {
    fit <- gaston::lmm.aireml(y, X, K = unname(K), min_tau = 0,
                             max_iter = max_iter, eps = 1e-7, verbose = FALSE)
    if (fit$niter >= max_iter) stop("REML reached its iteration limit.")
    fit
  }
  full <- fit_model(kernels)
  variance <- c(full$tau, full$sigma2)
  lrt <- sapply(seq_along(kernels), function(j) {
    max(0, 2 * (full$logL - fit_model(kernels[-j])$logL))
  })
  p <- ifelse(lrt <= 1e-8, 1, 0.5 * pchisq(lrt, df = 1, lower.tail = FALSE))
  list(components = data.frame(component = c(names(kernels), "Residual"),
         variance = variance, proportion = variance / sum(variance), row.names = NULL),
       tests = data.frame(component = names(kernels), statistic = lrt, p_value = p),
       beta = setNames(as.numeric(full$BLUP_beta), colnames(X)), fit = full)
}
