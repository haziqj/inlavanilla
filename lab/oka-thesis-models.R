################################################################################
#
# Models from Motonori Oka's PhD thesis (Oka, 2026), written for
# inla_vanilla(): the designs, simulators and log-posteriors that
# lab/oka-thesis.R (the study, with its checks) and lab/oka-thesis.qmd (the
# report) share. Every model's latent variables are integrated out inside
# logpost. Source from the package root after devtools::load_all().
#
################################################################################

## =============================================================================
## Shared helpers
## =============================================================================

## ----- Quadrature and sums on the log scale ----------------------------------

# Nodes and weights of the n-point Gauss-Hermite rule for N(0, 1), by
# Golub-Welsch, and their product grid in k dimensions, copied from
# lab/binary-cfa-model.R.
gh_normal <- function(n) {
  b <- sqrt(seq_len(n - 1))
  jac <- matrix(0, n, n)
  jac[cbind(seq_len(n - 1), 2:n)] <- b
  jac[cbind(2:n, seq_len(n - 1))] <- b
  e <- eigen(jac, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}

gh_grid <- function(n, k) {
  g <- gh_normal(n)
  z <- as.matrix(expand.grid(rep(list(g$z), k)))
  w <- as.matrix(expand.grid(rep(list(g$w), k)))
  list(z = unname(z), log_w = rowSums(log(w)))
}

# Row-wise log-sum-exp
lse_rows <- function(m) {
  m_max <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
  m_max + log(rowSums(exp(m - m_max)))
}

# Unique response patterns and their counts, copied from lab/binary-cfa-model.R
collapse_patterns <- function(y) {
  key <- apply(y, 1, paste, collapse = "")
  first <- !duplicated(key)
  list(
    y = y[first, , drop = FALSE],
    count = as.vector(table(key)[key[first]])
  )
}

## ----- Checks ----------------------------------------------------------------

# The largest gap between grad() and central differences of logpost()
check_grad <- function(logpost, grad, theta, ..., h = 1e-5) {
  fd <- vapply(
    seq_along(theta),
    function(j) {
      e <- replace(numeric(length(theta)), j, h)
      (logpost(theta + e, ...) - logpost(theta - e, ...)) / (2 * h)
    },
    numeric(1)
  )
  max(abs(fd - grad(theta, ...)))
}

# Each parameter against its true value, on the scale inla_vanilla() sees
recovery <- function(fit, truth) {
  s <- fit$summary
  data.frame(
    true = truth,
    mean = s[, "Mean"],
    sd = s[, "SD"],
    q025 = s[, "2.5%"],
    q975 = s[, "97.5%"],
    covers = s[, "2.5%"] <= truth & truth <= s[, "97.5%"]
  )
}

print_rounded <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits)
  print(df)
}

## =============================================================================
## Model 1 --- Confirmatory M2PL (Chapter 2, Appendix 6.6)
## =============================================================================

## ----- Design ----------------------------------------------------------------

# Appendix 6.6 of the thesis: K = 2 correlated factors, J = 15 binary items and
# N = 2,000 respondents, with d_j ~ U(-1, 1), the nonzero loadings
# a_jk ~ U(0.5, 1.5) and a factor correlation of 0.5, as in Section 2.4.1.2.
# The thesis does not give the 15 x 2 indicator matrix, and the recipe of
# Section 2.4.1.2 (three stacked identities, then three copies of every
# two-factor pattern) gives only 9 items when K = 2. Six stacked identities
# and three items on both factors make up the 15. The thesis fits by maximum
# likelihood, so the priors are ours and weak: a_jk, d_j ~ N(0, 3^2) and
# rho ~ U(-1, 1).
m2pl_q <- rbind(do.call(rbind, rep(list(diag(2)), 6)), matrix(1, 3, 2))
m2pl_free <- which(m2pl_q == 1) # column-major positions of the loadings
m2pl_prior_sd <- 3
m2pl_n_gh <- 15 # nodes per factor, 225 in all

m2pl_names <- c(
  paste0("a", row(m2pl_q)[m2pl_free], "_", col(m2pl_q)[m2pl_free]),
  paste0("d", seq_len(nrow(m2pl_q))),
  "atanh_rho"
)

sim_m2pl <- function(n, a, d, rho) {
  xi <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, rho, rho, 1), 2))
  prob <- plogis(xi %*% t(a) + rep(d, each = n))
  y <- 1 * (matrix(runif(length(prob)), n) < prob)
  colnames(y) <- paste0("y", seq_len(ncol(y)))
  y
}

## ----- Log-posterior and gradient --------------------------------------------

# P(y_ij = 1 | xi_i) = logistic(d_j + a_j' xi_i), with xi_i ~ N(0, R) and R
# the 2 x 2 correlation matrix. At node q the factors are xi_q = L z_q, with L
# the Cholesky factor of R, so each response pattern's log-likelihood at every
# node is a matrix product, and its marginal is a log-sum-exp over the nodes.
m2pl_parts <- function(theta, y, nodes) {
  j <- ncol(y)
  n_a <- length(m2pl_free)
  a <- matrix(0, j, 2)
  a[m2pl_free] <- theta[seq_len(n_a)]
  d <- theta[n_a + seq_len(j)]
  rho <- tanh(theta[[n_a + j + 1]])
  z <- nodes$z
  xi <- cbind(z[, 1], rho * z[, 1] + sqrt(1 - rho^2) * z[, 2])
  eta <- tcrossprod(xi, a) + rep(d, each = nrow(z)) # nodes x items
  ell <- tcrossprod(y, plogis(eta, log.p = TRUE)) +
    tcrossprod(1 - y, plogis(-eta, log.p = TRUE)) +
    rep(nodes$log_w, each = nrow(y)) # patterns x nodes
  list(
    a = a,
    d = d,
    rho = rho,
    xi = xi,
    eta = eta,
    ell = ell,
    ll = lse_rows(ell)
  )
}

m2pl_logpost <- function(theta, y, count, nodes) {
  p <- m2pl_parts(theta, y, nodes)
  log_prior_rho <- log1p(-p$rho^2) - log(2) # U(-1, 1), Jacobian of tanh
  sum(count * p$ll) +
    sum(dnorm(p$a[m2pl_free], 0, m2pl_prior_sd, log = TRUE)) +
    sum(dnorm(p$d, 0, m2pl_prior_sd, log = TRUE)) +
    log_prior_rho
}

# Fisher's identity: the score of each pattern is its complete-data score
# averaged over the posterior weights of the nodes, and the complete-data
# score in eta_qj is y_ij - P(y_ij = 1 | xi_q).
m2pl_grad <- function(theta, y, count, nodes) {
  p <- m2pl_parts(theta, y, nodes)
  post <- count * exp(p$ell - p$ll) # patterns x nodes
  d_eta <- crossprod(post, y) - colSums(post) * plogis(p$eta) # nodes x items
  z <- nodes$z
  d_xi2 <- z[, 1] - p$rho / sqrt(1 - p$rho^2) * z[, 2] # d xi_q2 / d rho
  d_rho <- sum(drop(d_eta %*% p$a[, 2]) * d_xi2)
  c(
    t(crossprod(p$xi, d_eta))[m2pl_free] - p$a[m2pl_free] / m2pl_prior_sd^2,
    colSums(d_eta) - p$d / m2pl_prior_sd^2,
    d_rho * (1 - p$rho^2) - 2 * p$rho
  )
}

## =============================================================================
## Model 2 --- Multilevel logistic regression (Chapter 2, Section 2.4.1.1)
## =============================================================================

## ----- Design ----------------------------------------------------------------

# Section 2.4.1.1: N = 10,000 clusters of J_i = 10 binary responses, every
# covariate with a random coefficient, xi_i ~ N(mu, Sigma), with mu_1 = 0.3,
# the slopes mu_k ~ U(0.1, 1.1), Sigma with variances 0.1 and covariances
# 0.05, and covariates x_ijk ~ N(0, R) with correlations 0.25. The thesis has
# K = 5 and K = 10 random coefficients. The integral over xi_i runs over all
# 10,000 clusters in every logpost call, and a product rule of n nodes per
# coefficient costs n^K passes over the 100,000 responses. With n = 6, K = 3
# needs 216 nodes and K = 4 would need 1,296. So this keeps the thesis's N,
# J_i and true values and takes K = 3, a random intercept and two random
# slopes.
mlr_k <- 3
mlr_j <- 10
mlr_n_gh <- 6 # nodes per coefficient, 216 in all

# Builds the Cholesky factor L of a d x d correlation matrix from d(d - 1)/2
# unconstrained values, as Stan does: each w gives a canonical partial
# correlation z = tanh(w), and row i of L takes, from left to right, the share
# z of the squared length still left. It also returns the log of L's diagonal
# and the log-Jacobian of the map, as in lab/lkj-correlation.R.
cholesky_corr <- function(w, d) {
  z <- tanh(w)
  log_1mz2 <- 2 * (log(2) - abs(w) - log1p(exp(-2 * abs(w)))) # log(1 - z^2)
  L <- diag(d)
  log_diag <- numeric(d)
  log_jac <- sum(log_1mz2)
  k <- 0
  for (i in seq_len(d)[-1]) {
    log_left <- 0
    for (j in seq_len(i - 1)) {
      k <- k + 1
      L[i, j] <- z[k] * exp(log_left / 2)
      log_jac <- log_jac + log_left / 2
      log_left <- log_left + log_1mz2[k]
    }
    L[i, i] <- exp(log_left / 2)
    log_diag[i] <- log_left / 2
  }
  list(L = L, log_diag = log_diag, log_jac = log_jac)
}

# The inverse map, from a correlation matrix to the w above
corr_to_w <- function(omega) {
  L <- t(chol(omega))
  d <- nrow(L)
  w <- numeric(0)
  for (i in seq_len(d)[-1]) {
    left <- 1
    for (j in seq_len(i - 1)) {
      z <- L[i, j] / sqrt(left)
      w <- c(w, atanh(z))
      left <- left - L[i, j]^2
    }
  }
  w
}

# The LKJ(eta) log density of L, with its normalising constant
# (Lewandowski et al., 2009), as in lab/lkj-correlation.R.
lkj_cholesky_lpdf <- function(log_diag, eta) {
  d <- length(log_diag)
  col <- sequence(seq_len(d - 1))
  b <- eta + (d - 1 - col) / 2
  log_const <- -sum(lbeta(b, b) + (2 * b - 1) * log(2))
  sum((d - seq_len(d) + 2 * eta - 2) * log_diag) + log_const
}

mlr_names <- c(
  paste0("mu", seq_len(mlr_k)),
  paste0("log_sigma", seq_len(mlr_k)),
  paste0(
    "atanh_z",
    rep(seq_len(mlr_k)[-1], seq_len(mlr_k - 1)),
    sequence(seq_len(mlr_k - 1))
  )
)

# Cluster i in row i of y, and its covariates in x[i, , ], one slice per
# random slope.
sim_mlr <- function(n, mu, sigma) {
  r <- matrix(0.25, mlr_k - 1, mlr_k - 1)
  diag(r) <- 1
  x <- matrix(rnorm(n * mlr_j * (mlr_k - 1)), ncol = mlr_k - 1) %*% chol(r)
  x <- array(x, c(n, mlr_j, mlr_k - 1))
  xi <- matrix(rnorm(n * mlr_k), n) %*% chol(sigma) + rep(mu, each = n)
  eta <- xi[, 1]
  for (k in seq_len(mlr_k)[-1]) {
    eta <- eta + x[,, k - 1] * xi[, k]
  }
  y <- 1 * (matrix(runif(n * mlr_j), n) < plogis(eta))
  list(y = y, x = x)
}

## ----- Log-posterior ---------------------------------------------------------

# Sigma = diag(sigma) R diag(sigma), with R a correlation matrix. At node q,
# xi_q = mu + diag(sigma) L z_q, and cluster i's log-likelihood is a sum over
# its J_i responses, one N x J_i matrix per node. Priors: mu_k ~ N(0, 5^2),
# sigma_k ~ half-N(0, 1) and R ~ LKJ(1), uniform over correlation matrices,
# all weak beside 100,000 responses. No gradient: through the map to L it is
# long to write.
mlr_ll <- function(theta, y, x, nodes) {
  k <- mlr_k
  chol_r <- cholesky_corr(theta[-seq_len(2 * k)], k)
  xi <- theta[seq_len(k)] +
    exp(theta[k + seq_len(k)]) * (chol_r$L %*% t(nodes$z)) # k x nodes
  # sum_j y_ij eta_ij is linear in xi_q, so only log(1 + exp(eta)) needs
  # every response at every node.
  y_x <- cbind(rowSums(y), apply(x * as.vector(y), 3, rowSums))
  ell <- y_x %*% xi
  for (q in seq_len(ncol(xi))) {
    eta <- xi[1, q]
    for (s in seq_len(k)[-1]) {
      eta <- eta + x[,, s - 1] * xi[s, q]
    }
    ell[, q] <- ell[, q] - rowSums(log1p(exp(eta))) # eta stays small here
  }
  sum(lse_rows(ell + rep(nodes$log_w, each = nrow(y))))
}

mlr_logpost <- function(theta, y, x, nodes) {
  k <- mlr_k
  sigma <- exp(theta[k + seq_len(k)])
  chol_r <- cholesky_corr(theta[-seq_len(2 * k)], k)
  mlr_ll(theta, y, x, nodes) +
    sum(dnorm(theta[seq_len(k)], 0, 5, log = TRUE)) +
    sum(dnorm(sigma, 0, 1, log = TRUE) + log(2) + log(sigma)) +
    lkj_cholesky_lpdf(chol_r$log_diag, 1) +
    chol_r$log_jac
}

# The mean and the entries of Sigma on or below the diagonal, from theta
mlr_natural <- function(theta) {
  k <- mlr_k
  sigma <- exp(unname(theta[k + seq_len(k)]))
  L <- sigma * cholesky_corr(theta[-seq_len(2 * k)], k)$L
  big <- tcrossprod(L)
  low <- lower.tri(big, diag = TRUE)
  c(
    setNames(unname(theta[seq_len(k)]), paste0("mu", seq_len(k))),
    setNames(big[low], paste0("sigma", row(big)[low], col(big)[low]))
  )
}

## =============================================================================
## Model 3 --- Grade of membership model (Chapter 3, Section 3.3.1.1)
## =============================================================================

## ----- Design ----------------------------------------------------------------

# Section 3.3.1.1: N = 2,000 respondents answer J = 20 binary items, with
# K = 3 latent profiles. Respondent i has a mixed-membership vector xi_i on the
# simplex, xi_i ~ Dirichlet(alpha), and P(y_ij = 1 | xi_i) = sum_k xi_ik
# theta_jk. The true theta_jk ~ Beta(0.2, 0.2), so about 40% of them are
# within 0.01 of 0 or 1, and the non-uniform condition has
# alpha = (0.5, 1, 1.5). The code below is for K = 3.
gom_j <- 20
gom_k <- 3
gom_alpha <- c(0.5, 1, 1.5)
gom_n_gj <- ceiling((gom_j + 1) / 2) # Gauss-Jacobi nodes per stick, see below
gom_theta_ix <- seq_len(gom_j * gom_k) # positions in theta
gom_alpha_ix <- gom_j * gom_k + seq_len(gom_k)
gom_names <- c(
  paste0(
    "probit_theta",
    seq_len(gom_j),
    "_",
    rep(seq_len(gom_k), each = gom_j)
  ),
  paste0("log_alpha", seq_len(gom_k))
)

sim_gom <- function(n, theta, alpha) {
  g <- matrix(rgamma(n * length(alpha), alpha), n, byrow = TRUE)
  xi <- g / rowSums(g)
  1 * (matrix(runif(n * nrow(theta)), n) < tcrossprod(xi, theta))
}

## ----- The marginal likelihood is exact --------------------------------------

# The thesis integrates xi_i by MCMC on the sphere, and notes that earlier
# maximum likelihood fits with analytical integrals (Potthoff et al., 2000)
# were limited to a few items. But each factor sum_k xi_ik
# f(y_ij | theta_jk) is linear in xi_i, so the integrand is a polynomial of
# degree J in xi_i, and a Gauss rule for the Dirichlet with enough nodes is
# exact. Stick-breaking makes the Dirichlet a product of independent betas,
# xi_1 = u, xi_2 = (1 - u) v and xi_3 = (1 - u)(1 - v), with
# u ~ Beta(alpha_1, alpha_2 + alpha_3) and v ~ Beta(alpha_2, alpha_3), and the
# integrand is a polynomial of degree at most J in u and in v. A Gauss-Jacobi
# rule with n nodes per stick is exact to degree 2n - 1, so n = 11 nodes per
# stick, 121 in all, integrate J = 20 items exactly. The nodes move with
# alpha, so the rule is rebuilt in every call (two 11 x 11 eigenproblems).

# The n-point Gauss rule for Beta(a, b) on (0, 1), by Golub-Welsch on the
# recurrence of the Jacobi polynomials with weight (1 - x)^(b - 1)
# (1 + x)^(a - 1) on (-1, 1), with the first diagonal and off-diagonal entries
# written out to avoid 0 / 0 when a + b is 1 or 2.
gauss_beta <- function(n, a, b) {
  al <- b - 1
  be <- a - 1
  s <- al + be
  i <- seq_len(n) - 1
  centre <- (be^2 - al^2) / ((2 * i + s) * (2 * i + s + 2))
  centre[1] <- (be - al) / (s + 2)
  i <- seq_len(n - 1)
  num <- 4 * i * (i + al) * (i + be) * (i + s)
  off2 <- num / ((2 * i + s)^2 * (2 * i + s + 1) * (2 * i + s - 1))
  off2[1] <- 4 * (1 + al) * (1 + be) / ((2 + s)^2 * (3 + s))
  jac <- diag(centre, n)
  jac[cbind(i, i + 1)] <- jac[cbind(i + 1, i)] <- sqrt(off2)
  e <- eigen(jac, symmetric = TRUE)
  list(t = (1 + e$values) / 2, w = e$vectors[1, ]^2)
}

gom_rule <- function(alpha, n) {
  u <- gauss_beta(n, alpha[1], alpha[2] + alpha[3])
  v <- gauss_beta(n, alpha[2], alpha[3])
  iu <- rep(seq_len(n), times = n)
  iv <- rep(seq_len(n), each = n)
  uu <- u$t[iu]
  list(
    xi = cbind(uu, (1 - uu) * v$t[iv], (1 - uu) * (1 - v$t[iv])),
    log_w = log(u$w[iu]) + log(v$w[iv])
  )
}

## ----- Log-posterior ---------------------------------------------------------

# The thesis fits by maximum likelihood, with alpha kept above 1e-16. The
# priors here are theta_jk ~ Beta(0.2, 0.2), the distribution the thesis draws
# them from, and log alpha_k ~ N(0, 1). With theta_jk ~ U(0, 1) instead, the
# probabilities the data put near 0 or 1 are pulled inwards and alpha drops
# to make up for it: in a Stan check, even the exact posterior missed all
# three true alpha. Near 0 or 1 the likelihood flattens out, so a marginal
# takes the shape of the prior's tail there. On the logit scale a beta tail
# is exponential, which no skew-normal follows. On the probit scale it is
# normal: U(0, 1) becomes N(0, 1), and Beta(0.2, 0.2) has a tail like
# N(0, 5). No gradient: logpost costs a few milliseconds, and inla_vanilla()
# gets the VB step from logpost values alone.
gom_prior_a <- 0.2
gom_ll <- function(theta, y, count, n_gj = gom_n_gj) {
  probit <- matrix(theta[gom_theta_ix], gom_j, gom_k)
  rule <- gom_rule(exp(theta[gom_alpha_ix]), n_gj)
  p1 <- tcrossprod(rule$xi, pnorm(probit)) # nodes x items
  p0 <- tcrossprod(rule$xi, pnorm(-probit)) # 1 - p1, without cancellation
  ell <- tcrossprod(y, log(p1)) +
    tcrossprod(1 - y, log(p0)) +
    rep(rule$log_w, each = nrow(y))
  sum(count * lse_rows(ell))
}

# The Beta(a, a) log density of pnorm(z), with the Jacobian dnorm(z)
log_beta_probit <- function(z, a) {
  (a - 1) *
    (pnorm(z, log.p = TRUE) + pnorm(-z, log.p = TRUE)) -
    lbeta(a, a) +
    dnorm(z, log = TRUE)
}

gom_logpost <- function(theta, y, count) {
  gom_ll(theta, y, count) +
    sum(log_beta_probit(theta[gom_theta_ix], gom_prior_a)) +
    sum(dnorm(theta[gom_alpha_ix], 0, 1, log = TRUE))
}

## =============================================================================
## Model 4 --- Multilevel latent regression scaling model (Chapter 4)
## =============================================================================

## ----- Design ----------------------------------------------------------------

# Section 4.4.1, unidimensional condition: N = 5,000 students in G = 50
# schools (100 each, the thesis does not give the sizes), 150 items in 5
# booklets of 30, one booklet per student at random. Items are 30 GPCM (scores
# 0 to 2, steps -0.5 and 0.5), 60 3PL and 60 2PL, with fixed, known
# parameters. The thesis draws them uniformly between the smallest and largest
# calibrated TIMSS 2023 values, which it does not print: here a ~ U(0.6, 2),
# b ~ U(-2, 2) and c ~ U(0.05, 0.3). There are 60 predictors, 20 continuous,
# 20 binary and 20 with three categories (80 columns once dummy coded), cut
# from N(0, R) with correlations 0.1. In the thesis each predictor has a
# nonzero effect with probability 0.2, and its draw gave 13, so here 13 of the
# 60 variables at random (both dummies of a chosen three-category variable),
# each nonzero coefficient of size U(0.2, 0.8) with a random sign. The scaling
# model is
#   theta_gi = eta_g + x_gi' beta + e_gi,  e_gi ~ N(0, 1),
#   eta_g ~ N(0, 0.1), a variance.
# Predictors are complete here: the thesis's imputation of missing predictors
# and its e-value selection are Gibbs and multiple-imputation machinery
# outside inla_vanilla(). What it fits is the scaling model once the
# predictors are chosen (steps 4 and 6 of the thesis's procedure), and the
# Bayesian cross-validation of step 5 compares two candidate models through
# held-out items: 6 of each booklet's 30 items, at random.
scal_per <- 100 # students per school
scal_r <- 30 # nodes of each student's inner rule
scal_q <- 5 # adaptive Gauss-Hermite nodes per school
scal_newton <- 4 # Newton steps to each school's mode

# Log-probability of the scores y (a vector over students, or one score) of
# one item at the trait values th (a matrix with a row per student, or a
# vector).
item_loglik <- function(th, y, a, b, guess, type) {
  if (type == "gpcm") {
    z1 <- a * (th - b - 0.5) # a (theta - b + d_1), d_1 = -0.5
    z2 <- 2 * a * (th - b) # a (2 (theta - b) + d_1 + d_2), d_2 = 0.5
    z_max <- pmax(0, z1, z2)
    lse <- z_max + log(exp(-z_max) + exp(z1 - z_max) + exp(z2 - z_max))
    (y == 1) * z1 + (y == 2) * z2 - lse
  } else {
    p1 <- guess + (1 - guess) * plogis(a * (th - b))
    y * log(p1) + (1 - y) * log1p(-p1)
  }
}

# Each student's log-likelihood of the items in `use` (a logical matrix, the
# items each student took and that count) at the trait values in th.
items_loglik <- function(th, dat, use) {
  out <- matrix(0, nrow(th), ncol(th))
  for (j in seq_len(ncol(dat$resp))) {
    who <- which(use[, j])
    if (length(who) == 0) {
      next
    }
    it <- dat$items[j, ]
    out[who, ] <- out[who, ] +
      item_loglik(
        th[who, , drop = FALSE],
        dat$resp[who, j],
        it$a,
        it$b,
        it$guess,
        it$type
      )
  }
  out
}

sim_scaling <- function(n_school, n_per) {
  n <- n_school * n_per
  school <- rep(seq_len(n_school), each = n_per)
  items <- data.frame(
    type = sample(rep(c("gpcm", "3pl", "2pl"), c(30, 60, 60))),
    a = runif(150, 0.6, 2),
    b = runif(150, -2, 2)
  )
  items$guess <- ifelse(items$type == "3pl", runif(150, 0.05, 0.3), 0)
  in_booklet <- split(sample(150), rep(1:5, each = 30))
  test_items <- unlist(lapply(in_booklet, sample, 6))
  booklet <- sample(5, n, replace = TRUE)
  took <- matrix(FALSE, n, 150)
  for (k in 1:5) {
    took[booklet == k, in_booklet[[k]]] <- TRUE
  }
  r <- matrix(0.1, 60, 60)
  diag(r) <- 1
  u <- matrix(rnorm(n * 60), n) %*% chol(r)
  cat3 <- (u[, 41:60] > -0.5) + (u[, 41:60] >= 0.5)
  dummies <- 1 * (cat3[, rep(1:20, each = 2)] == rep(1:2, each = n))
  x <- cbind(u[, 1:20], 1 * (u[, 21:40] > 0), dummies)
  colnames(x) <- c(
    paste0("x", 1:40),
    paste0("x", rep(41:60, each = 2), "_", 1:2)
  )
  var_of <- c(1:40, rep(41:60, each = 2))
  beta <- (seq_len(60) %in% sample(60, 13))[var_of] *
    runif(80, 0.2, 0.8) *
    sample(c(-1, 1), 80, replace = TRUE)
  eta <- rnorm(n_school, 0, sqrt(0.1))
  trait <- eta[school] + drop(x %*% beta) + rnorm(n)
  # Scores by inversion: the number of cumulative probabilities below a
  # uniform draw (a binary item has no score 2, so its second cumulative
  # probability is 1).
  resp <- matrix(NA_real_, n, 150)
  for (j in seq_len(150)) {
    who <- which(took[, j])
    it <- items[j, ]
    prob <- vapply(
      0:2,
      function(s) {
        exp(item_loglik(trait[who], s, it$a, it$b, it$guess, it$type))
      },
      numeric(length(who))
    )
    if (it$type != "gpcm") {
      prob[, 3] <- 0
    }
    u_j <- runif(length(who))
    resp[who, j] <- (u_j > prob[, 1]) + (u_j > prob[, 1] + prob[, 2])
  }
  list(
    resp = resp,
    items = items,
    train = took & !(col(took) %in% test_items),
    test = took & col(took) %in% test_items,
    x = x,
    beta = beta,
    eta = eta,
    school = school
  )
}

## ----- Each student's inner rule ---------------------------------------------

# With the item parameters fixed, student i's likelihood L_i(theta) never
# changes, so its rule is built once. Let m_i and v_i be the mean and SD of
# L_i(theta) dnorm(theta, 0, 2), from a grid, and (z_r, w_r) the Gauss-Hermite
# rule for N(0, 1). Then with theta_ir = m_i + v_i z_r,
#   int L_i(theta) dnorm(theta, c, sigma) dtheta
#     = sum_r omega_ir dnorm(theta_ir, c, sigma),
#   omega_ir = w_r L_i(theta_ir) / dnorm(theta_ir, m_i, v_i).
# The dnorm(theta, 0, 2) only places the nodes, and cancels from omega.
scal_inner_rule <- function(dat, use, n_r) {
  grid <- seq(-10, 10, by = 0.1)
  n <- nrow(dat$resp)
  lg <- items_loglik(matrix(grid, n, length(grid), byrow = TRUE), dat, use) +
    rep(dnorm(grid, 0, 2, log = TRUE), each = n)
  w <- exp(lg - lse_rows(lg))
  m <- drop(w %*% grid)
  v <- sqrt(drop(w %*% grid^2) - m^2)
  gh <- gh_normal(n_r)
  nodes <- m + outer(v, gh$z)
  list(
    m = m,
    v = v,
    nodes = nodes,
    log_omega = items_loglik(nodes, dat, use) +
      rep(log(gh$w), each = n) -
      dnorm(nodes, m, v, log = TRUE)
  )
}

## ----- Log-posterior ---------------------------------------------------------

# School g contributes the log of
#   int dnorm(eta, mu, tau) prod_i I_i(eta + x_i' beta) deta,
# with I_i the inner integral above, by adaptive Gauss-Hermite quadrature:
# Newton steps to the mode in eta (the first two derivatives of log I_i come
# from the inner rule's weights), then scal_q nodes scaled by the curvature
# there. It returns one value per school. Given `test`, the log-likelihood of
# each student's held-out items at the inner nodes, it returns instead each
# student's expected log-likelihood of those items, averaged over the
# posterior of their trait and school: the thesis's measure for the Bayesian
# cross-validation (Section 4.2.3.6).
scal_parts <- function(
  theta,
  rule,
  xs,
  school,
  test = NULL,
  n_q = scal_q,
  n_newton = scal_newton
) {
  p <- ncol(xs)
  mu <- theta[[p + 1]]
  sigma <- exp(theta[[p + 2]])
  tau <- exp(theta[[p + 3]])
  xb <- if (p > 0) drop(xs %*% theta[seq_len(p)]) else 0
  inner <- function(eta, moments = TRUE) {
    centre <- eta[school] + xb
    lw <- rule$log_omega - (rule$nodes - centre)^2 / (2 * sigma^2)
    li <- lse_rows(lw)
    out <- list(li = li - log(sigma) - 0.5 * log(2 * pi))
    if (moments) {
      out$w <- exp(lw - li)
      m1 <- rowSums(out$w * rule$nodes)
      m2 <- rowSums(out$w * rule$nodes^2)
      out$d1 <- (m1 - centre) / sigma^2
      out$d2 <- (m2 - m1^2) / sigma^4 - 1 / sigma^2
    }
    out
  }
  # Newton steps from the mode the schools would have if each student's
  # likelihood were normal, with the mean and variance of the inner rule
  # (whose dnorm(theta, 0, 2) is taken out).
  prec <- pmax(1 / rule$v^2 - 1 / 4, 1e-2)
  wt <- 1 / (sigma^2 + 1 / prec)
  gap <- rule$m / (rule$v^2 * prec) - xb # likelihood mean less x' beta
  eta <- (mu / tau^2 + as.vector(rowsum(wt * gap, school))) /
    (1 / tau^2 + as.vector(rowsum(wt, school)))
  for (it in seq_len(n_newton + 1)) {
    f <- inner(eta)
    h2 <- -1 / tau^2 + as.vector(rowsum(f$d2, school))
    if (it > n_newton) {
      break
    }
    h1 <- -(eta - mu) / tau^2 + as.vector(rowsum(f$d1, school))
    eta <- eta - h1 / h2
  }
  s <- 1 / sqrt(-h2)
  gh <- gh_normal(n_q)
  h <- matrix(0, length(eta), n_q)
  if (!is.null(test)) {
    pred <- matrix(0, length(school), n_q)
  }
  for (q in seq_len(n_q)) {
    eta_q <- eta + s * gh$z[q]
    f <- inner(eta_q, moments = !is.null(test))
    h[, q] <- dnorm(eta_q, mu, tau, log = TRUE) +
      as.vector(rowsum(f$li, school)) +
      log(gh$w[q]) -
      dnorm(gh$z[q], log = TRUE)
    if (!is.null(test)) {
      pred[, q] <- rowSums(f$w * test)
    }
  }
  ll_school <- log(s) + lse_rows(h)
  if (is.null(test)) {
    return(ll_school)
  }
  post <- exp(h - ll_school + log(s)) # each school's weights on its nodes
  rowSums(post[school, , drop = FALSE] * pred)
}

# The thesis's prior families, with values of our own (it does not give
# them): beta | s^2 ~ N(0, s^2 I) with s^2 ~ IG(1, 1), which integrates to a
# multivariate t, sigma^2 and tau^2 ~ IG(1, 0.5) (the 1 x 1 inverse Wishart
# with 2 degrees of freedom and scale 1), and mu ~ N(0, 10^2).
scal_logprior <- function(theta, p) {
  beta <- theta[seq_len(p)]
  log_var <- 2 * theta[p + 2:3]
  ig <- log(0.5) - 2 * log_var - 0.5 * exp(-log_var) # IG(1, 0.5) of each
  lgamma(1 + p / 2) -
    p / 2 * log(2 * pi) -
    (1 + p / 2) * log(1 + sum(beta^2) / 2) +
    dnorm(theta[[p + 1]], 0, 10, log = TRUE) +
    sum(ig + log(2) + log_var) # with the Jacobian of log SD -> variance
}

scal_logpost <- function(theta, rule, xs, school) {
  sum(scal_parts(theta, rule, xs, school)) + scal_logprior(theta, ncol(xs))
}
