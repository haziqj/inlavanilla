################################################################################
#
# Binary CFA models shared by lab/binary-cfa.R (theta parameterisation) and
# lab/binary-cfa-delta.R (delta parameterisation). Probit (normal-ogive)
# factor models for binary items, y_ij = 1 when lambda_j eta_i + e_ij > tau_j,
# with standard normal factors and errors (the theta parameterisation, as in
# lavaan). The log-posterior integrates the factors out by Gauss-Hermite
# quadrature, over the counts of the distinct response patterns, so
# inla_vanilla() sees only the loadings, the thresholds and the factor
# correlation, and the inner integral runs inside every logpost call.
#
# Models: one factor with 6 items (m = 12), and two correlated factors with 4
# items each (m = 17, a 2-D inner integral). Priors: lambda_j, tau_j ~ N(0, 2^2)
# and rho ~ U(-1, 1), on atanh(rho). Also the Stan model of the full joint,
# factor scores sampled, which serves as the reference. Source from the
# package root after devtools::load_all().
#
################################################################################

## ----- Priors and true values ------------------------------------------------

prior_sd <- 2

models <- list(
  "one factor" = list(
    item_factor = rep(1, 6),
    lambda = c(0.6, 0.8, 1.0, 1.2, 1.5, 0.7),
    tau = c(-1, -0.5, 0, 0.3, 0.8, 1.2),
    rho = NULL,
    n_gh = 21
  ),
  "two factors" = list(
    item_factor = rep(1:2, each = 4),
    lambda = c(0.6, 0.8, 1.0, 1.2, 1.5, 0.7, 0.9, 1.1),
    tau = c(-1, -0.5, 0, 0.3, 0.8, 1.2, -0.2, 0.5),
    rho = 0.4,
    n_gh = 15
  )
)

## ----- Gauss-Hermite quadrature ----------------------------------------------

# Nodes and weights of the n-point Gauss-Hermite rule for N(0, 1), by
# Golub-Welsch, and their product grid in k dimensions.
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

## ----- Log-posterior and gradient --------------------------------------------

theta_names <- function(model) {
  p <- length(model$lambda)
  c(
    paste0("lambda", seq_len(p)),
    paste0("tau", seq_len(p)),
    if (!is.null(model$rho)) "atanh_rho"
  )
}

true_theta <- function(model) {
  setNames(
    c(model$lambda, model$tau, if (!is.null(model$rho)) atanh(model$rho)),
    theta_names(model)
  )
}

# The start: unit loadings, zero thresholds and no factor correlation. The
# loadings are left unconstrained, and the fit stays on the positive mode, far
# from the mirror image with every sign flipped.
theta_init <- function(model) {
  init <- true_theta(model)
  init[] <- 0
  init[startsWith(names(init), "lambda")] <- 1
  init
}

# Shared by the log-posterior and its gradient. At quadrature node q the
# factors are eta_q = L z_q, with L the Cholesky factor of the factor
# correlation matrix, and item j has linear predictor a_qj = lambda_j eta_q,f(j)
# - tau_j. The log-likelihood of response pattern i at node q is then a matrix
# product, y log Phi(a) + (1 - y) log Phi(-a), and the pattern's marginal
# log-likelihood is a log-sum-exp over the nodes.
cfa_parts <- function(theta, y, item_factor, nodes) {
  p <- ncol(y)
  lambda <- theta[seq_len(p)]
  tau <- theta[p + seq_len(p)]
  z <- nodes$z
  rho <- if (ncol(z) == 2) tanh(theta[[2 * p + 1]])
  eta <- z
  if (!is.null(rho)) {
    eta[, 2] <- rho * z[, 1] + sqrt(1 - rho^2) * z[, 2]
  }
  n_q <- nrow(z)
  a <- eta[, item_factor, drop = FALSE] *
    matrix(lambda, n_q, p, byrow = TRUE) -
    matrix(tau, n_q, p, byrow = TRUE)
  log_p1 <- pnorm(a, log.p = TRUE)
  log_p0 <- pnorm(-a, log.p = TRUE)
  ell <- tcrossprod(y, log_p1) +
    tcrossprod(1 - y, log_p0) +
    matrix(nodes$log_w, nrow(y), n_q, byrow = TRUE)
  ell_max <- apply(ell, 1, max)
  list(
    lambda = lambda,
    tau = tau,
    rho = rho,
    eta = eta,
    a = a,
    log_p1 = log_p1,
    log_p0 = log_p0,
    ell = ell,
    ll = ell_max + log(rowSums(exp(ell - ell_max)))
  )
}

cfa_logpost <- function(theta, y, count, item_factor, nodes) {
  cp <- cfa_parts(theta, y, item_factor, nodes)
  lp <- sum(count * cp$ll) +
    sum(dnorm(cp$lambda, 0, prior_sd, log = TRUE)) +
    sum(dnorm(cp$tau, 0, prior_sd, log = TRUE))
  # rho ~ U(-1, 1), with the Jacobian of rho = tanh(atanh_rho)
  if (!is.null(cp$rho)) {
    lp <- lp + log1p(-cp$rho^2) - log(2)
  }
  lp
}

# The score of each pattern is the average over its posterior node weights of
# the complete-data score, whose derivative in a_qj is the signed Mills ratio
# of the probit.
cfa_grad <- function(theta, y, count, item_factor, nodes) {
  cp <- cfa_parts(theta, y, item_factor, nodes)
  post <- count * exp(cp$ell - cp$ll)
  log_phi <- dnorm(cp$a, log = TRUE)
  d_a <- crossprod(post, y) *
    exp(log_phi - cp$log_p1) -
    crossprod(post, 1 - y) * exp(log_phi - cp$log_p0)
  eta_j <- cp$eta[, item_factor, drop = FALSE]
  out <- c(
    colSums(d_a * eta_j) - cp$lambda / prior_sd^2,
    -colSums(d_a) - cp$tau / prior_sd^2
  )
  if (!is.null(cp$rho)) {
    z <- nodes$z
    d_eta2 <- z[, 1] - cp$rho / sqrt(1 - cp$rho^2) * z[, 2]
    on_2 <- item_factor == 2
    d_rho <- sum(d_a[, on_2, drop = FALSE] * outer(d_eta2, cp$lambda[on_2]))
    out <- c(out, d_rho * (1 - cp$rho^2) - 2 * cp$rho)
  }
  out
}

## ----- Simulation and fitting ------------------------------------------------

sim_cfa <- function(n, model) {
  p <- length(model$lambda)
  eta <- matrix(rnorm(n * max(model$item_factor)), n)
  if (!is.null(model$rho)) {
    eta[, 2] <- model$rho * eta[, 1] + sqrt(1 - model$rho^2) * eta[, 2]
  }
  y_star <- eta[, model$item_factor, drop = FALSE] *
    matrix(model$lambda, n, p, byrow = TRUE) +
    matrix(rnorm(n * p), n)
  y <- 1 * (y_star > matrix(model$tau, n, p, byrow = TRUE))
  colnames(y) <- paste0("y", seq_len(p))
  y
}

# Unique response patterns and their counts: at most 2^p rows whatever n
collapse_patterns <- function(y) {
  key <- apply(y, 1, paste, collapse = "")
  first <- !duplicated(key)
  list(
    y = y[first, , drop = FALSE],
    count = as.vector(table(key)[key[first]])
  )
}

fit_cfa <- function(
  y,
  model,
  logpost = cfa_logpost,
  grad = cfa_grad,
  init = theta_init(model),
  ...
) {
  pat <- collapse_patterns(y)
  inla_vanilla(
    logpost,
    init,
    grad,
    y = pat$y,
    count = pat$count,
    item_factor = model$item_factor,
    nodes = gh_grid(model$n_gh, max(model$item_factor)),
    ...
  )
}

## ----- Stan reference --------------------------------------------------------

# The full joint, factor scores sampled, with the same priors, so the reference
# depends on neither the quadrature nor the approximation. The loadings are
# constrained positive, which picks the same mode as inla_vanilla()'s start:
# the true loadings are all over 4 posterior SDs from zero. For one factor, rho
# is sampled from its prior and ignored.
cfa_stan <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1, upper=2> K;
  array[N, J] int<lower=0, upper=1> y;
  array[J] int<lower=1, upper=K> f;
  real<lower=0> prior_sd;
}
parameters {
  vector<lower=0>[J] lambda;
  vector[J] tau;
  real<lower=-1, upper=1> rho;
  matrix[N, K] z;
}
model {
  matrix[N, K] eta = z;
  if (K == 2) {
    eta[, 2] = rho * z[, 1] + sqrt(1 - square(rho)) * z[, 2];
  }
  to_vector(z) ~ std_normal();
  lambda ~ normal(0, prior_sd);
  tau ~ normal(0, prior_sd);
  for (j in 1:J) {
    y[, j] ~ bernoulli(Phi(lambda[j] * eta[, f[j]] - tau[j]));
  }
}
"

# Draws on the theta scale, with atanh(rho) in place of rho, as inla_vanilla()
# sees them. n_draws per chain, 4 chains.
run_stan <- function(y, model, stan_model, n_draws) {
  fit <- rstan::sampling(
    stan_model,
    data = list(
      N = nrow(y),
      J = ncol(y),
      K = max(model$item_factor),
      y = matrix(as.integer(y), nrow(y)),
      f = model$item_factor,
      prior_sd = prior_sd
    ),
    chains = 4,
    cores = 4,
    iter = 2 * n_draws,
    warmup = n_draws,
    seed = 1,
    refresh = 0,
    pars = c("lambda", "tau", "rho"), # not the factor scores: they swamp RAM
    init = function() {
      list(lambda = rep(1, ncol(y)), tau = rep(0, ncol(y)), rho = 0)
    }
  )
  pars <- c("lambda", "tau", if (!is.null(model$rho)) "rho")
  draws <- as.matrix(fit, pars = pars)
  if (!is.null(model$rho)) {
    draws[, "rho"] <- atanh(draws[, "rho"])
  }
  colnames(draws) <- theta_names(model)
  summ <- rstan::summary(fit, pars = pars)$summary
  elapsed <- rstan::get_elapsed_time(fit)
  list(
    draws = draws,
    time = max(rowSums(elapsed)),
    rhat = max(summ[, "Rhat"]),
    n_div = sum(rstan::get_divergent_iterations(fit))
  )
}
