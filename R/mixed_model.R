# Fixed-effect comparisons corresponding to the preliminary regressions.
compare_prs_models <- function(y, prs, treatment) {
  dat <- data.frame(y = y, prs = prs, treatment = treatment)
  if (any(!vapply(dat, is.numeric, logical(1))) || any(!is.finite(as.matrix(dat))) ||
      length(y) != length(prs) || length(y) != length(treatment)) {
    stop("Provide aligned finite numeric y, prs, and treatment vectors.")
  }
  design <- stats::model.matrix(~ prs * treatment, dat)
  if (nrow(dat) <= ncol(design) || qr(design)$rank != ncol(design)) {
    stop("Fixed-effect comparison requires a full-rank design and residual degrees of freedom.")
  }
  base <- stats::lm(y ~ treatment, dat)
  score <- stats::lm(y ~ prs + treatment, dat)
  interaction <- stats::lm(y ~ prs * treatment, dat)
  r2 <- vapply(list(base, score, interaction), function(m) summary(m)$r.squared, numeric(1))
  data.frame(model = c("Treatment", "Treatment + PRS", "Treatment + PRS + PRS:T"),
             r_squared = r2, delta_r_squared = c(NA_real_, diff(r2)),
             p_value = c(NA_real_, stats::anova(base, score)$`Pr(>F)`[2L],
                         stats::anova(score, interaction)$`Pr(>F)`[2L]))
}

# Kernel convention follows mixedModel.R: standardize dosage and, separately,
# dosage multiplied by raw treatment. This is not diag(T) %*% K_G %*% diag(T).
build_kernels <- function(genotypes, treatment) {
  genotypes <- .numeric_matrix(genotypes, "genotypes", allow_na = TRUE)
  if (ncol(genotypes) < 1L || nrow(genotypes) < 3L ||
      !is.numeric(treatment) || length(treatment) != nrow(genotypes) ||
      any(!is.finite(treatment)) || stats::sd(treatment) <= 1e-10) {
    stop("Provide genotypes and a finite, varying numeric treatment vector.")
  }
  prep_g <- fit_numeric_preprocessor(genotypes, impute = TRUE)
  if (!length(prep_g$keep)) stop("No variable genotype columns remain.")
  dosage <- genotypes[, prep_g$keep, drop = FALSE]
  for (j in seq_along(prep_g$keep)) {
    dosage[is.na(dosage[, j]), j] <- prep_g$center[prep_g$keep[j]]
  }
  z <- apply_numeric_preprocessor(prep_g, genotypes)
  product <- dosage * treatment
  prep_gt <- fit_numeric_preprocessor(product)
  if (!length(prep_gt$keep)) stop("No variable genotype-by-treatment columns remain.")
  w <- apply_numeric_preprocessor(prep_gt, product)
  # Each kernel uses its own retained column count; metadata records exclusions.
  kernels <- list(GxT = tcrossprod(w) / ncol(w), G = tcrossprod(z) / ncol(z))
  attr(kernels, "snp_counts") <- c(input = ncol(genotypes), G = ncol(z), GxT = ncol(w))
  attr(kernels, "kept_G") <- prep_g$keep
  attr(kernels, "kept_GxT") <- prep_g$keep[prep_gt$keep]
  kernels
}

.validate_lmm <- function(y, X, kernels) {
  X <- .numeric_matrix(X, "X")
  if (!is.numeric(y) || any(!is.finite(y)) || length(y) != nrow(X) ||
      length(y) <= ncol(X) || qr(X)$rank != ncol(X) || stats::var(y) <= 1e-12) {
    stop("y must vary, and X must have matching rows and full column rank.")
  }
  if (!is.list(kernels) || !length(kernels) || is.null(names(kernels)) ||
      any(!nzchar(names(kernels))) || anyDuplicated(names(kernels)) ||
      "Residual" %in% names(kernels)) stop("Supply uniquely named kernels.")
  for (k in kernels) {
    if (!is.matrix(k) || !is.numeric(k) || any(!is.finite(k)) ||
        !identical(dim(k), c(length(y), length(y))) ||
        max(abs(k - t(k))) > 1e-8 * max(1, max(abs(k)))) {
      stop("Each kernel must be a finite symmetric n-by-n matrix.")
    }
    eig <- eigen(k, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) < -1e-8 * max(1, max(eig)) || mean(diag(k)) <= 0) {
      stop("Kernels must be nonzero positive semidefinite matrices.")
    }
  }
  X
}

.fit_aireml <- function(y, X, kernels, max_iter) {
  fit <- gaston::lmm.aireml(y, X, K = unname(kernels), min_tau = 0,
                            max_iter = max_iter, eps = 1e-7, verbose = FALSE)
  if (any(!is.finite(c(fit$tau, fit$sigma2, fit$logL))) ||
      fit$niter >= max_iter) stop("REML fit failed or reached max_iter; inspect the design.")
  fit
}

# Same fixed-effect design must be used in full and restricted REML models.
fit_variance_components <- function(y, X, kernels, max_iter = 200L) {
  .require_package("gaston")
  X <- .validate_lmm(y, X, kernels)
  fit <- .fit_aireml(y, X, kernels, max_iter)
  variance <- c(as.numeric(fit$tau), fit$sigma2)
  diagonal <- c(vapply(kernels, function(k) mean(diag(k)), numeric(1)), 1)
  contribution <- variance * diagonal
  components <- data.frame(component = c(names(kernels), "Residual"),
                           variance = variance,
                           coefficient_share = variance / sum(variance),
                           average_variance = contribution,
                           conditional_share = contribution / sum(contribution),
                           row.names = NULL)
  list(components = components, beta = setNames(as.numeric(fit$BLUP_beta), colnames(X)),
       logLik = fit$logL, iterations = fit$niter, gradient_norm = fit$norm_grad,
       raw_fit = fit, y = y, X = X, kernels = kernels, max_iter = max_iter)
}

.lrt_result <- function(full_logLik, null_logLik, value) {
  statistic <- 2 * (full_logLik - null_logLik)
  if (statistic < -1e-4 * max(1, abs(full_logLik))) {
    stop("The restricted likelihood exceeds the full fit; resolve optimization before inference.")
  }
  statistic <- max(0, statistic)
  # At the atom at zero, use the inclusive upper tail, P(T >= 0) = 1.
  p <- if (value == 0) {
    if (statistic <= 1e-8) 1 else 0.5 * stats::pchisq(statistic, 1, lower.tail = FALSE)
  } else stats::pchisq(statistic, 1, lower.tail = FALSE)
  c(statistic = statistic, p_value = p)
}

# Profiles a named variance coefficient. The zero-boundary mixture is approximate;
# its regularity conditions are described in docs/METHODS.md.
test_variance_component <- function(model, component, value = 0) {
  .require_package("gaston")
  index <- match(component, names(model$kernels))
  if (length(index) != 1L || is.na(index) || length(value) != 1L ||
      !is.finite(value) || value < 0) stop("Specify a kernel name and nonnegative value.")
  m <- length(model$kernels)
  if (value == 0 && m > 1L) {
    null <- .fit_aireml(model$y, model$X, model$kernels[-index], model$max_iter)
    null_logLik <- null$logL
  } else {
    free <- setdiff(seq_len(m + 1L), index)
    start <- c(as.numeric(model$raw_fit$tau), model$raw_fit$sigma2)
    lower <- c(rep(0, m), 1e-6)[free]
    objective <- function(par) {
      theta <- numeric(m + 1L)
      theta[index] <- value
      theta[free] <- par
      -as.numeric(gaston::lmm.restricted.likelihood(model$y, model$X,
        unname(model$kernels), tau = theta[seq_len(m)], s2 = theta[m + 1L]))
    }
    opt <- stats::nlminb(pmax(start[free], lower), objective, lower = lower,
                         control = list(iter.max = 500L, eval.max = 1500L, rel.tol = 1e-9))
    if (opt$convergence != 0L || !is.finite(opt$objective)) {
      stop("Restricted optimization did not converge: ", opt$message)
    }
    null_logLik <- -opt$objective
  }
  result <- .lrt_result(model$logLik, null_logLik, value)
  data.frame(component = component, null_value = value,
             statistic = unname(result["statistic"]), p_value = unname(result["p_value"]),
             reference = if (value == 0) "0.5 point mass at 0 + 0.5 chi-square(1)" else "chi-square(1)")
}

# Returns accepted grid values and their range, not an interpolated exact interval.
profile_variance_component <- function(model, component, grid, level = 0.95) {
  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid)) ||
      any(grid < 0) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("Provide a nonnegative grid and a confidence level between 0 and 1.")
  }
  grid <- sort(unique(grid))
  table <- do.call(rbind, lapply(grid, function(v) test_variance_component(model, component, v)))
  accepted <- table$p_value >= 1 - level
  list(table = table, accepted = grid[accepted], level = level,
       interval = if (any(accepted)) range(grid[accepted]) else c(NA_real_, NA_real_),
       upper_grid_edge_accepted = tail(accepted, 1L),
       lower_grid_edge_accepted = head(accepted, 1L),
       disconnected = sum(diff(c(FALSE, accepted)) == 1L) > 1L)
}
