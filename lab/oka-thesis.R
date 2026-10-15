################################################################################
#
# Models from Motonori Oka's PhD thesis, fitted with inla_vanilla() (October
# 2026). Chapters 2 and 3 of the thesis (Oka, 2026) estimate large latent
# variable models by marginal maximum likelihood, with stochastic
# optimisation and Langevin samplers for the latent variables, and Chapter 4
# by a Gibbs sampler. inla_vanilla() has no latent field, so the latent
# variables must be integrated out inside logpost, and a model fits when that
# integral is cheap enough to evaluate thousands of times. Four of the
# thesis's models qualify, at or scaled down from the sizes of its own
# simulation studies. The designs, simulators and log-posteriors are in
# lab/oka-thesis-models.R, and lab/oka-thesis.qmd writes up one data set per
# model.
#
#   1. Chapter 2, Appendix 6.6: the confirmatory M2PL with K = 2 correlated
#      factors, J = 15 binary items and N = 2,000, by product Gauss-Hermite
#      quadrature (m = 34 parameters, analytic gradient).
#   2. Chapter 2, Section 2.4.1.1: the multilevel logistic regression with a
#      random intercept and two random slopes (K = 3, where the thesis has 5
#      or 10), 10,000 clusters of 10, by product Gauss-Hermite quadrature
#      (m = 9, no gradient).
#   3. Chapter 3, Section 3.3.1.1: the grade of membership model with K = 3
#      profiles, J = 20 items and N = 2,000. Its marginal likelihood is exact,
#      by a Gauss-Jacobi rule on the simplex (m = 63, no gradient).
#   4. Chapter 4, Section 4.4.1: the multilevel latent regression scaling
#      model, with fixed 2PL, 3PL and GPCM items, 5,000 students in 50
#      schools and the predictors chosen. Each student's trait integrates out
#      by a precomputed rule, each school's intercept by adaptive quadrature
#      (m = 3 plus a coefficient per predictor column, 19 here). The thesis's
#      Bayesian cross-validation on held-out items then compares two
#      candidate models.
#
# Out of reach, because the integral is too big or m too large:
#   - The M2PL with K = 5 or 10 (m = 150 or 650) and the IPIP-NEO application
#     (K = 30, m = 1,035), which need a 5- to 30-D integral per respondent.
#   - The multilevel logistic regression with K = 5 or 10 random coefficients
#     (5-D or 10-D integrals over 10,000 clusters in every call).
#   - The grade of membership model with K = 6 (J = 40 or 80, m = 246 or 486)
#     and its NLTCS application (K = 9, J = 16, an 8-D simplex whose exact rule
#     needs about 9^8 nodes per response pattern).
#   - The latent surface model for aggregated relational data (Section
#     3.2.1.2). With K = 3 each respondent's integral is over the sphere and a
#     log degree, by a product rule on the sphere and adaptive Gauss-Hermite
#     quadrature, but at the thesis's N = 2,000 a trial had not finished its
#     mode and Hessian after 12,000 calls. A Karnataka village (N of about 64,
#     m = 26) would be quick, but the data need an openICPSR sign-in and the
#     thesis does not print its anchor vectors. With K = 6 the sphere is S^5.
#   - The scaling model with K = 2 and all 80 predictor columns (m = 168), and
#     the imputation models, which serve a Gibbs sampler.
#
# Reference: Oka, M. (2026). *Towards scalable psychometrics: Estimation,
# selection and evaluation of large-scale latent variable models* [Doctoral
# thesis, London School of Economics and Political Science].
#
# Each model's data are simulated from the thesis's design, with its own seed.
# Needs lme4 only for an optional check. Takes about half an hour, on one
# core. From the package root:
#   Rscript lab/oka-thesis.R
#
# Results (October 2026, one run of 26 minutes on one core):
#
#   1. M2PL: 90 s, and every check of summary() passes. 33 of the 34 true
#      values lie in their 95% intervals. The mean absolute errors are 0.078
#      for the loadings and 0.031 for the intercepts (Table 6.11, K = 5: 0.082
#      and 0.058). Against Stan (NUTS on the full joint, factor scores
#      sampled, 4 chains of 2,000 draws, R-hat 1.003, no divergences), every
#      posterior mean is within 0.014 posterior SDs, every SD within 4%, and
#      every 2.5% and 97.5% quantile within 0.11 SDs. The log-likelihood at
#      the truth is -19,096.8, against Table 6.12's -18,286 (SD 82) from the
#      thesis's own indicator matrix and true values.
#   2. Multilevel logistic regression: 10 minutes, and every check passes but
#      a mean gap of 0.10 SDs on one partial correlation. 8 of the 9 true
#      values lie in their intervals, and the thesis's averaged MAE is 0.0098
#      (Table 2.2, K = 5: 0.0113). glmer()'s Laplace estimates of the three
#      variances are 16% to 25% below these posterior means (0.087, 0.082 and
#      0.077, against 0.1). The 6-node rule's error changes by under 0.01
#      nats across the posterior, so the shortfall is this data set's.
#   3. Grade of membership: 10 minutes. The point estimates recover the truth:
#      alpha has posterior means 1.43, 1.02 and 0.58 (true 1.5, 1 and 0.5),
#      and the mean absolute errors are 0.031 for theta and 0.056 for alpha
#      (Figure 3.3: about 0.03 and 0.09). But only 42 of the 63 true values
#      lie in their intervals. The VB step does not settle, and summary()
#      flags the marginals of the theta_jk near 0 or 1: heavy tails, mass
#      beyond the scan, ridge bends and mean gaps. Stan on the full joint
#      (memberships sampled, each chain relabelled, R-hat 1.03 after that, 80
#      divergences in 4,000 draws) covered 62 of the 63, so these skew-normal
#      marginals are too narrow. With theta_jk ~ U(0, 1), Stan's exact
#      posterior itself missed all three alpha (means 0.91, 0.65 and 0.35) and
#      covered 33 of the 63. On the logit scale, with the uniform prior, over
#      half of some marginals' mass lay beyond the scanned window.
#   4. Scaling model: 2 minutes for both candidates. The selected model has 16
#      coefficients (m = 19), every check passes, and 18 of the 19 true values
#      lie in their intervals. Each student's rule is within 1.4e-4 of a fine
#      grid, and 5 adaptive nodes per school within 6e-8 of 15. On the held-out
#      items the expected log-likelihood is -14,647 with the true predictors
#      and -14,784 without, a difference of 137 (SE 13).
#
################################################################################

## ----- Configuration ---------------------------------------------------------

smoke <- FALSE
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("lab/oka-thesis-models.R")
parts <- c("m2pl", "mlr", "gom", "scaling") # the models to fit
nsamp <- if (smoke) 1000 else 4000 # posterior draws, Models 1 to 3
n_cv <- if (smoke) 50 else 200 # draws for Model 4's cross-validation
m2pl_n <- if (smoke) 500 else 2000 # respondents
mlr_n <- if (smoke) 2000 else 10000 # clusters
gom_n <- if (smoke) 500 else 2000 # respondents
scal_g <- if (smoke) 20 else 50 # schools of 100 students

## =============================================================================
## Model 1 --- Confirmatory M2PL (Chapter 2, Appendix 6.6)
## =============================================================================

## ----- Fit -------------------------------------------------------------------

if ("m2pl" %in% parts) {
  cat("\n===== Model 1: confirmatory M2PL, K = 2, J = 15 =====\n")
  set.seed(2026)
  m2pl_a <- m2pl_q * matrix(runif(length(m2pl_q), 0.5, 1.5), nrow(m2pl_q))
  m2pl_d <- runif(nrow(m2pl_q), -1, 1)
  m2pl_rho <- 0.5
  m2pl_truth <- setNames(
    c(m2pl_a[m2pl_free], m2pl_d, atanh(m2pl_rho)),
    m2pl_names
  )
  m2pl_pat <- collapse_patterns(sim_m2pl(m2pl_n, m2pl_a, m2pl_d, m2pl_rho))
  m2pl_nodes <- gh_grid(m2pl_n_gh, 2)

  # The marginal log-likelihood at the true values, by the thesis's benchmark
  # rule (41 nodes per factor) and by the rule of the fit. Table 6.12 reports
  # -18,286.1 (SD 82.2) over 100 data sets from its own indicator matrix.
  m2pl_ll <- function(n_gh) {
    sum(
      m2pl_pat$count *
        m2pl_parts(m2pl_truth, m2pl_pat$y, gh_grid(n_gh, 2))$ll
    )
  }
  m2pl_ll_ref <- m2pl_ll(41)
  cat(sprintf(
    paste(
      "%d patterns | loglik at the truth %.3f (41 nodes per factor),",
      "%d nodes per factor off by %.1e\n"
    ),
    nrow(m2pl_pat$y),
    m2pl_ll_ref,
    m2pl_n_gh,
    m2pl_ll(m2pl_n_gh) - m2pl_ll_ref
  ))
  cat(sprintf(
    "grad vs central differences at the truth: %.1e\n",
    check_grad(
      m2pl_logpost,
      m2pl_grad,
      m2pl_truth,
      y = m2pl_pat$y,
      count = m2pl_pat$count,
      nodes = m2pl_nodes
    )
  ))

  # Positive loadings at the start, as in Section 2.4.1.3, put the fit on the
  # mode with the true signs of the factors. Flipping a factor's sign flips
  # its loadings and rho, so the posterior has four mirror-image modes, and
  # fit$logml is log(4) short of the whole posterior's.
  m2pl_init <- setNames(
    c(rep(1, length(m2pl_free)), numeric(nrow(m2pl_q)), 0),
    m2pl_names
  )
  m2pl_fit <- inla_vanilla(
    m2pl_logpost,
    m2pl_init,
    m2pl_grad,
    y = m2pl_pat$y,
    count = m2pl_pat$count,
    nodes = m2pl_nodes
  )
  print(summary(m2pl_fit))
  m2pl_rec <- recovery(m2pl_fit, m2pl_truth)
  print_rounded(m2pl_rec)

  # The correlation on its own scale, and the mean absolute errors of the
  # posterior means. Table 6.11 gives 0.082 for the loadings, 0.058 for the
  # intercepts and 0.025 for the correlations (QN-SOMALA, N = 2,000, K = 5,
  # J = 50, 100 data sets), for scale only: here the correlation's error is a
  # single |rho_hat - rho|.
  m2pl_rho_draws <- sampling(
    m2pl_fit,
    nsamp = nsamp,
    fun = function(theta) c(rho = tanh(theta[["atanh_rho"]]))
  )
  m2pl_group <- sub("[0-9_]+$", "", m2pl_names)
  m2pl_mae <- tapply(abs(m2pl_rec$mean - m2pl_rec$true), m2pl_group, mean)
  m2pl_mae[["atanh_rho"]] <- abs(mean(m2pl_rho_draws) - m2pl_rho)
  cat(sprintf(
    "rho: mean %.3f, 95%% interval (%.3f, %.3f), true %.1f\n",
    mean(m2pl_rho_draws),
    quantile(m2pl_rho_draws, 0.025),
    quantile(m2pl_rho_draws, 0.975),
    m2pl_rho
  ))
  cat(sprintf(
    "MAE: a %.3f, d %.3f, rho %.3f | coverage of the 95%% intervals %d/%d\n",
    m2pl_mae[["a"]],
    m2pl_mae[["d"]],
    m2pl_mae[["atanh_rho"]],
    sum(m2pl_rec$covers),
    nrow(m2pl_rec)
  ))
}

## =============================================================================
## Model 2 --- Multilevel logistic regression (Chapter 2, Section 2.4.1.1)
## =============================================================================

## ----- Fit -------------------------------------------------------------------

if ("mlr" %in% parts) {
  cat("\n===== Model 2: multilevel logistic regression, K = 3 =====\n")
  set.seed(2027)
  mlr_mu <- c(0.3, runif(mlr_k - 1, 0.1, 1.1))
  mlr_sigma <- matrix(0.05, mlr_k, mlr_k)
  diag(mlr_sigma) <- 0.1
  mlr_truth <- setNames(
    c(
      mlr_mu,
      log(sqrt(diag(mlr_sigma))),
      corr_to_w(cov2cor(mlr_sigma))
    ),
    mlr_names
  )
  mlr_dat <- sim_mlr(mlr_n, mlr_mu, mlr_sigma)
  mlr_nodes <- gh_grid(mlr_n_gh, mlr_k)
  mlr_ll_ref <- mlr_ll(mlr_truth, mlr_dat$y, mlr_dat$x, gh_grid(10, mlr_k))
  cat(sprintf(
    paste(
      "loglik at the truth %.3f (10 nodes per coefficient),",
      "%d nodes per coefficient off by %.1e\n"
    ),
    mlr_ll_ref,
    mlr_n_gh,
    mlr_ll(mlr_truth, mlr_dat$y, mlr_dat$x, mlr_nodes) - mlr_ll_ref
  ))

  # The start of Section 2.4.1.3: mu from U(0, 1.5) and Sigma = I
  mlr_init <- setNames(
    c(runif(mlr_k, 0, 1.5), numeric(mlr_k + choose(mlr_k, 2))),
    mlr_names
  )
  mlr_fit <- inla_vanilla(
    mlr_logpost,
    mlr_init,
    y = mlr_dat$y,
    x = mlr_dat$x,
    nodes = mlr_nodes
  )
  print(summary(mlr_fit))
  print_rounded(recovery(mlr_fit, mlr_truth))

  # On the natural scale, with the thesis's error measure: the mean of
  # MAE(mu) and MAE(Sigma), Sigma's diagonal included. Table 2.2 gives 0.0113
  # (D-SOMALA with minibatches of 250, K = 5, 100 data sets, after 2,000 s),
  # for scale only.
  mlr_draws <- sampling(mlr_fit, nsamp = nsamp, fun = mlr_natural)
  mlr_true_nat <- mlr_natural(mlr_truth)
  print_rounded(
    data.frame(
      true = mlr_true_nat,
      mean = colMeans(mlr_draws),
      q025 = apply(mlr_draws, 2, quantile, 0.025),
      q975 = apply(mlr_draws, 2, quantile, 0.975)
    ),
    4
  )
  mlr_err <- abs(colMeans(mlr_draws) - mlr_true_nat)
  is_mu <- startsWith(names(mlr_err), "mu")
  cat(sprintf(
    "averaged MAE (mu, Sigma): %.4f\n",
    (mean(mlr_err[is_mu]) + mean(mlr_err[!is_mu])) / 2
  ))

  # lme4's glmer() maximises the Laplace approximation of the same marginal
  # likelihood (its adaptive quadrature takes one random effect only).
  if (requireNamespace("lme4", quietly = TRUE)) {
    long <- data.frame(
      y = as.vector(t(mlr_dat$y)),
      x2 = as.vector(t(mlr_dat$x[,, 1])),
      x3 = as.vector(t(mlr_dat$x[,, 2])),
      cluster = rep(seq_len(mlr_n), each = mlr_j)
    )
    t_glmer <- system.time(
      glmer_fit <- lme4::glmer(
        y ~ x2 + x3 + (x2 + x3 | cluster),
        data = long,
        family = binomial
      )
    )[["elapsed"]]
    vc <- lme4::VarCorr(glmer_fit)$cluster
    cat(sprintf(
      "glmer (%.0f s): mu %s | Sigma %s\n",
      t_glmer,
      paste(sprintf("%.4f", lme4::fixef(glmer_fit)), collapse = " "),
      paste(sprintf("%.4f", vc[lower.tri(vc, diag = TRUE)]), collapse = " ")
    ))
  }
}

## =============================================================================
## Model 3 --- Grade of membership model (Chapter 3, Section 3.3.1.1)
## =============================================================================

## ----- Fit -------------------------------------------------------------------

if ("gom" %in% parts) {
  cat("\n===== Model 3: grade of membership, K = 3, J = 20 =====\n")
  set.seed(2028)
  gom_theta <- matrix(rbeta(gom_j * gom_k, 0.2, 0.2), gom_j, gom_k)
  gom_pat <- collapse_patterns(sim_gom(gom_n, gom_theta, gom_alpha))
  gom_truth <- setNames(c(qnorm(gom_theta), log(gom_alpha)), gom_names)
  gb <- gauss_beta(4, 0.5, 2)
  cat(sprintf(
    "E[t^7] under Beta(0.5, 2), 4-node rule vs exact: %.1e\n",
    sum(gb$w * gb$t^7) - beta(7.5, 2) / beta(0.5, 2)
  ))
  cat(sprintf(
    "%d patterns | loglik at the truth %.4f, 30 nodes per stick %.1e off\n",
    nrow(gom_pat$y),
    gom_ll(gom_truth, gom_pat$y, gom_pat$count),
    gom_ll(gom_truth, gom_pat$y, gom_pat$count, 30) -
      gom_ll(gom_truth, gom_pat$y, gom_pat$count)
  ))

  # The thesis's start, theta_j = (0.4, 0.5, 0.6) for every item and
  # alpha_k = 2, breaks the symmetry between the profiles. The posterior has
  # 3! = 6 copies of the mode, one per labelling, so fit$logml is that of one
  # copy, log(6) short of the whole posterior's.
  gom_init <- setNames(
    c(qnorm(rep(c(0.4, 0.5, 0.6), each = gom_j)), rep(log(2), 3)),
    gom_names
  )
  gom_fit <- inla_vanilla(
    gom_logpost,
    gom_init,
    y = gom_pat$y,
    count = gom_pat$count
  )
  print(summary(gom_fit))

  # Profiles carry no labels, so match them to the true ones first.
  gom_perms <- as.matrix(expand.grid(1:3, 1:3, 1:3))
  gom_perms <- gom_perms[apply(gom_perms, 1, function(p) all(sort(p) == 1:3)), ]
  gom_est <- pnorm(matrix(gom_fit$summary[gom_theta_ix, "Mean"], gom_j))
  gom_gap <- apply(gom_perms, 1, function(p) sum(abs(gom_est - gom_theta[, p])))
  gom_perm <- gom_perms[which.min(gom_gap), ]
  gom_truth_aligned <- setNames(
    c(qnorm(gom_theta[, gom_perm]), log(gom_alpha[gom_perm])),
    gom_names
  )
  gom_rec <- recovery(gom_fit, gom_truth_aligned)
  print_rounded(gom_rec[c(seq_len(6), gom_alpha_ix), ])

  # Errors on the natural scale, against the thesis's Figure 3.3 for this
  # condition after 500 s of stochastic optimisation (read off the plot:
  # MAE(theta) about 0.03 and MAE(alpha) about 0.09 for the manifold MALA,
  # 0.075 and 0.52 for data-augmentation Gibbs).
  gom_draws <- sampling(gom_fit, nsamp = nsamp, fun = function(theta) {
    c(pnorm(theta[gom_theta_ix]), alpha = exp(theta[gom_alpha_ix]))
  })
  gom_post_mean <- colMeans(gom_draws)
  cat(sprintf(
    "MAE: theta %.3f, alpha %.3f | alpha: mean %s, true %s\n",
    mean(abs(gom_post_mean[gom_theta_ix] - gom_theta[, gom_perm])),
    mean(abs(gom_post_mean[gom_alpha_ix] - gom_alpha[gom_perm])),
    paste(round(gom_post_mean[gom_alpha_ix], 3), collapse = " "),
    paste(gom_alpha[gom_perm], collapse = " ")
  ))
  # Coverage falls short for the theta_jk near 0 or 1, whose marginals here
  # are too narrow. In a Stan check of this data set, the exact posterior
  # covered 62 of the 63 true values.
  cat(sprintf(
    "coverage of the 95%% intervals %d/%d\n",
    sum(gom_rec$covers),
    nrow(gom_rec)
  ))
}

## =============================================================================
## Model 4 --- Multilevel latent regression scaling model (Chapter 4)
## =============================================================================

## ----- Fit and cross-validate ------------------------------------------------

if ("scaling" %in% parts) {
  cat("\n===== Model 4: multilevel latent regression, K = 1 =====\n")
  set.seed(2029)
  scal_dat <- sim_scaling(scal_g, scal_per)
  scal_school <- scal_dat$school
  scal_rule <- scal_inner_rule(scal_dat, scal_dat$train, scal_r)
  scal_test <- items_loglik(scal_rule$nodes, scal_dat, scal_dat$test)

  # Two candidate scaling models, as steps 3 to 5 would produce: the true
  # predictors (what a perfect selection keeps) and none at all. The second
  # keeps the school intercepts, unlike the thesis's fixed-intercept null
  # model.
  scal_cand <- list(selected = which(scal_dat$beta != 0), none = integer(0))
  scal_param_names <- function(sel) {
    c(sprintf("b_%s", colnames(scal_dat$x)[sel]), "mu", "log_sigma", "log_tau")
  }
  scal_truth <- setNames(
    c(scal_dat$beta[scal_cand$selected], 0, 0, log(sqrt(0.1))),
    scal_param_names(scal_cand$selected)
  )
  scal_xs <- scal_dat$x[, scal_cand$selected, drop = FALSE]

  # At the true values: the inner rule against a grid of 0.01 on (-12, 12)
  # for 500 students, the adaptive rule against 15 nodes, and the Newton
  # steps against 20.
  check_ix <- seq_len(500)
  dat_check <- scal_dat
  dat_check$resp <- scal_dat$resp[check_ix, ]
  grid <- seq(-12, 12, by = 0.01)
  centre <- scal_dat$eta[scal_school] + drop(scal_dat$x %*% scal_dat$beta)
  centre <- centre[check_ix]
  on_grid <- lse_rows(
    items_loglik(
      matrix(grid, length(check_ix), length(grid), byrow = TRUE),
      dat_check,
      scal_dat$train[check_ix, ]
    ) +
      dnorm(rep(grid, each = length(check_ix)), centre, 1, log = TRUE) +
      log(0.01)
  )
  by_rule <- lse_rows(
    scal_rule$log_omega[check_ix, ] +
      dnorm(scal_rule$nodes[check_ix, ], centre, 1, log = TRUE)
  )
  cat(sprintf(
    "inner rule (%d nodes) vs grid, %d students: largest gap %.1e\n",
    scal_r,
    length(check_ix),
    max(abs(by_rule - on_grid))
  ))
  scal_ll <- function(...) {
    sum(scal_parts(scal_truth, scal_rule, scal_xs, scal_school, ...))
  }
  scal_ll_ref <- scal_ll(n_q = 15)
  scal_ll_fit <- scal_ll()
  scal_ll_newton <- scal_ll(n_newton = 20)
  cat(sprintf(
    paste(
      "loglik at the truth %.4f | %d adaptive nodes per school off by %.1e,",
      "%d Newton steps off by %.1e\n"
    ),
    scal_ll_ref,
    scal_q,
    scal_ll_fit - scal_ll_ref,
    scal_newton,
    scal_ll_fit - scal_ll_newton
  ))

  scal_fits <- list()
  scal_elpd <- list()
  for (nm in names(scal_cand)) {
    sel <- scal_cand[[nm]]
    xs <- scal_dat$x[, sel, drop = FALSE]
    init <- setNames(
      c(numeric(length(sel)), 0, 0, log(0.5)),
      scal_param_names(sel)
    )
    cat(sprintf("\n--- %s: %d coefficients ---\n", nm, length(sel)))
    fit <- inla_vanilla(
      scal_logpost,
      init,
      rule = scal_rule,
      xs = xs,
      school = scal_school
    )
    print(summary(fit))
    scal_fits[[nm]] <- fit

    # Bayesian cross-validation (step 5, Section 4.2.3.6): each student's
    # log-likelihood of their held-out items, averaged over the posterior of
    # their trait and school and over draws of the parameters, then summed
    # over students. The SE below is across students, where the thesis's SD
    # is across its 50 imputed data sets.
    pred <- sampling(
      fit,
      nsamp = n_cv,
      fun = function(theta) {
        scal_parts(theta, scal_rule, xs, scal_school, test = scal_test)
      }
    )
    scal_elpd[[nm]] <- colMeans(pred)
  }
  scal_rec <- recovery(scal_fits$selected, scal_truth)
  print_rounded(scal_rec)
  cat(sprintf(
    "coverage of the 95%% intervals %d/%d\n",
    sum(scal_rec$covers),
    nrow(scal_rec)
  ))
  scal_diff <- scal_elpd$selected - scal_elpd$none
  cat(sprintf(
    paste(
      "held-out items, expected log-likelihood: selected %.2f, none %.2f,",
      "difference %.2f (SE %.2f)\n"
    ),
    sum(scal_elpd$selected),
    sum(scal_elpd$none),
    sum(scal_diff),
    sqrt(length(scal_diff)) * sd(scal_diff)
  ))
}
