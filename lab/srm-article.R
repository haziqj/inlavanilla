################################################################################
#
# The multivariate social relations model (SRM) of Bhangale and Jorgensen
# (2025, Multivariate Behavioral Research, 60(5), 930--953), fitted with
# inla_vanilla(). In a round-robin design each of G groups has n people, and
# every person i rates every other person j on p = 3 variables:
# y_k,ij = E_k,i + A_k,j + R_k,ij. The ego and alter effects of a person are
# N(0, Sigma_EA) (6 x 6), and the relationship effects of a dyad,
# (R_k,ij, R_k,ji), are N(0, Sigma_R) (6 x 6, with the symmetry of
# indistinguishable dyads): 9 SDs and 24 correlations in all. The random
# effects are Gaussian, so they integrate out exactly, and in a complete
# round-robin design the marginal likelihood needs only 6 x 6 matrices. The
# data are simulated from the population values of the paper's Tables 1 and 2,
# with zero group means (the paper centres each group instead), in its
# smallest and largest designs, and fitted with its default priors. If rstan is
# installed, the same posterior is sampled with NUTS as a check. Run from the
# package root:
#   Rscript lab/srm-article.R
#
################################################################################

## ----- Configuration ---------------------------------------------------------

library(inlavanilla)
set.seed(2025)
run_stan <- requireNamespace("rstan", quietly = TRUE)
nsamp <- 4000 # draws from each inla_vanilla() fit
stan_iter <- 6000 # NUTS iterations per chain, the first 1000 warmup
p <- 3 # round-robin variables
# fmt: skip
designs <- data.frame(
  G = c(10, 25), # groups
  n = c(6, 20)   # people per group
)

## =============================================================================
## Population values
## =============================================================================

# Fills a d x d correlation matrix from the rows of its lower triangle
corr_from_rows <- function(r, d) {
  m <- diag(d)
  m[upper.tri(m)] <- r
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

# Case level (Table 2), in the paper's order E1, A1, E2, A2, E3, A3
case_names <- c("E1", "A1", "E2", "A2", "E3", "A3")
sd_ea_true <- c(0.775, 0.548, 0.881, 0.415, 0.629, 0.415)
# fmt: skip
cor_ea_true <- corr_from_rows(c(
  0.236,
  0.703, 0.124,
  0.093, 0.528, 0.099,
  0.574, 0.102, 0.606, 0.080,
  0.093, 0.528, 0.099, 0.419, 0.195
), 6)
sigma_ea_true <- cor_ea_true * tcrossprod(sd_ea_true)

# Dyad level (Table 1): the relationship SDs, and the intrapersonal and
# interpersonal correlations between the variables, with the dyadic
# reciprocities on the diagonal of the interpersonal ones
sd_r_true <- c(0.949, 0.940, 1.255)
rho_true <- c(0.167, 0.222, 0.060)
intra_true <- corr_from_rows(c(0.538, 0.705, 0.569), 3)
inter_true <- corr_from_rows(c(0.135, 0.176, 0.142), 3)
diag(inter_true) <- rho_true

# Sigma_R in the order R1,ij, R1,ji, R2,ij, ..., R3,ji. The swap matrix
# exchanges each ij entry with its ji entry, and each E with its A.
swap <- kronecker(diag(p), matrix(c(0, 1, 1, 0), 2))
dyad_cov <- function(sd_r, intra, inter) {
  cor_r <- kronecker(intra, diag(2)) + kronecker(inter, swap[1:2, 1:2])
  cor_r * tcrossprod(rep(sd_r, each = 2))
}
sigma_r_true <- dyad_cov(sd_r_true, intra_true, inter_true)

## =============================================================================
## Data
## =============================================================================

## ----- Simulate a round-robin design -----------------------------------------

# A G x n x n x p array, y[g, i, j, k] being person i's rating of person j on
# variable k in group g, with NA on the diagonals
simulate_srm <- function(G, n) {
  y <- array(NA_real_, c(G, n, n, p))
  e <- 2 * seq_len(p) - 1 # E_k and R_k,ij in the 6-vectors
  a <- 2 * seq_len(p) # A_k and R_k,ji
  for (g in seq_len(G)) {
    u <- matrix(rnorm(n * 2 * p), n) %*% chol(sigma_ea_true)
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        r <- drop(rnorm(2 * p) %*% chol(sigma_r_true))
        y[g, i, j, ] <- u[i, e] + u[j, a] + r[e]
        y[g, j, i, ] <- u[j, e] + u[i, a] + r[a]
      }
    }
  }
  y
}

## ----- Sums of squares -------------------------------------------------------

# Write x_ij = (y1,ij, y1,ji, ..., y3,ji) for the ratings in dyad {ij}, and
# s_i = sum_j x_ij for person i's totals, given and received. The likelihood
# below sees the data only through W, the sum of x_ij x_ij' over the dyads, C,
# the sum of s_i s_i' over the people, and TT, the sum over the groups of the
# outer product of the group total sum_i s_i. All groups must have n people.
srm_stats <- function(y) {
  G <- dim(y)[1]
  n <- dim(y)[2]
  pairs <- which(upper.tri(diag(n)), arr.ind = TRUE) # i < j
  W <- C <- TT <- matrix(0, 2 * p, 2 * p)
  for (g in seq_len(G)) {
    x <- s <- NULL
    for (k in seq_len(p)) {
      yk <- y[g, , , k]
      x <- cbind(x, yk[pairs], yk[pairs[, 2:1]])
      s <- cbind(s, rowSums(yk, na.rm = TRUE), colSums(yk, na.rm = TRUE))
    }
    W <- W + crossprod(x)
    C <- C + crossprod(s)
    TT <- TT + tcrossprod(colSums(s))
  }
  list(G = G, n = n, W = W, C = C, TT = TT)
}

## =============================================================================
## Model
## =============================================================================

## ----- Marginal likelihood ---------------------------------------------------

# With u_i = (E1,i, A1,i, ..., A3,i) for person i's effects, each dyad has
# x_ij = u_i + S u_j + R_ij, where S is the swap matrix, which commutes with
# Q = Sigma_R^-1. Integrating out the u_i by Woodbury's identity needs the
# inverse of the 6n x 6n matrix I_n x K0 + J_n x QS of a group, with
# K0 = Sigma_EA^-1 + (n - 1) Q - QS. Its eigenvectors split it into one block
# K1 = K0 + n QS and n - 1 copies of K0, so a call works only with 6 x 6
# matrices, and its cost grows with neither G nor n.
srm_loglik <- function(sigma_ea, sigma_r, stats) {
  n <- stats$n
  chol_ea <- chol(sigma_ea)
  chol_r <- chol(sigma_r)
  Q <- chol2inv(chol_r)
  QS <- Q %*% swap
  K0 <- chol2inv(chol_ea) + (n - 1) * Q - QS
  K1 <- K0 + n * QS
  chol_k0 <- chol(K0)
  chol_k1 <- chol(K1)
  K0_inv <- chol2inv(chol_k0)
  K1_inv <- chol2inv(chol_k1)
  chols <- list(chol_r, chol_ea, chol_k1, chol_k0)
  log_dets <- vapply(chols, function(m) 2 * sum(log(diag(m))), numeric(1))
  log_det <- sum(c(n * (n - 1) / 2, n, 1, n - 1) * log_dets)
  quad <- sum(Q * stats$W) -
    sum(K0_inv * (Q %*% stats$C %*% Q)) -
    sum((K1_inv - K0_inv) * (Q %*% stats$TT %*% Q)) / n
  -0.5 * (stats$G * (p * n * (n - 1) * log(2 * pi) + log_det) + quad)
}

## ----- From unconstrained values to correlation matrices ---------------------

# log(1 - tanh(w)^2), exact when |w| is large and tanh(w) rounds to 1
log1m_tanh2 <- function(w) 2 * (log(2) - abs(w) - log1p(exp(-2 * abs(w))))

# The Cholesky factor L of a d x d correlation matrix from d(d - 1)/2
# unconstrained values, as Stan builds it (see lab/lkj-correlation.R): each w
# gives a CPC z = tanh(w), and row i of L has unit length, its entries taking,
# from left to right, the share z of the squared length still left
cholesky_corr <- function(w, d) {
  z <- tanh(w)
  log_1mz2 <- log1m_tanh2(w)
  L <- diag(d)
  log_diag <- numeric(d) # log L[i, i]
  log_jac <- sum(log_1mz2) # dz / dw
  k <- 0
  for (i in seq_len(d)[-1]) {
    log_left <- 0 # log of the squared length left in row i
    for (j in seq_len(i - 1)) {
      k <- k + 1
      L[i, j] <- z[k] * exp(log_left / 2)
      log_jac <- log_jac + log_left / 2 # dL[i, j] / dz
      log_left <- log_left + log_1mz2[k]
    }
    L[i, i] <- exp(log_left / 2)
    log_diag[i] <- log_left / 2
  }
  list(L = L, log_diag = log_diag, log_jac = log_jac)
}

# The row and column of L that each CPC fills, in the order of w
cpc_index <- function(d) {
  data.frame(
    row = rep(seq_len(d)[-1], seq_len(d - 1)),
    col = sequence(seq_len(d - 1))
  )
}

## ----- Parameters ------------------------------------------------------------

# Sigma_R is positive definite exactly when the covariance matrices of the sums
# R_k,ij + R_k,ji and of the differences R_k,ij - R_k,ji both are. Their
# correlation matrices, Omega_P and Omega_M, are free, so they are built from
# CPCs, as Omega_EA is, and the intrapersonal and interpersonal correlations
# follow: intra = (a * Omega_P + b * Omega_M) / 2 and
# inter = (a * Omega_P - b * Omega_M) / 2, with
# a_kl = sqrt((1 + rho_k) (1 + rho_l)) and b_kl = sqrt((1 - rho_k) (1 - rho_l)).
# Every theta then gives a valid Sigma_R, where Stan, as lavaan.srm sets the
# model up, rejects the proposals that do not.
cpc_ea <- cpc_index(2 * p)
cpc_r <- cpc_index(p)
theta_names <- c(
  paste0("log_sd_", case_names),
  paste0("cpc_", case_names[cpc_ea$row], "_", case_names[cpc_ea$col]),
  paste0("log_sd_R", seq_len(p)),
  paste0("atanh_rho_R", seq_len(p)),
  paste0("cpc_P", cpc_r$row, cpc_r$col),
  paste0("cpc_M", cpc_r$row, cpc_r$col)
)
idx <- split(seq_along(theta_names), rep(1:6, c(6, 15, 3, 3, 3, 3)))
names(idx) <- c("sd_ea", "ea", "sd_r", "rho", "P", "M")

# Splits theta into the SRM's SDs and correlations, with the log-Jacobian of
# the map from theta to the 9 SDs and the 24 free correlations. Each L -> Omega
# adds sum_i (d - i) log L[i, i] (lab/lkj-correlation.R), and each pair (k, l)
# of (Omega_P, Omega_M) -> (intra, inter) adds log(a_kl b_kl / 2).
unpack <- function(theta) {
  ea <- cholesky_corr(theta[idx$ea], 2 * p)
  P <- cholesky_corr(theta[idx$P], p)
  M <- cholesky_corr(theta[idx$M], p)
  rho <- tanh(theta[idx$rho])
  log1m_rho2 <- log1m_tanh2(theta[idx$rho])
  a <- sqrt(outer(1 + rho, 1 + rho))
  b <- sqrt(outer(1 - rho, 1 - rho))
  intra <- (a * tcrossprod(P$L) + b * tcrossprod(M$L)) / 2
  inter <- (a * tcrossprod(P$L) - b * tcrossprod(M$L)) / 2
  low <- lower.tri(intra)
  log_jac <- sum(theta[c(idx$sd_ea, idx$sd_r)]) +
    ea$log_jac +
    sum((2 * p - seq_len(2 * p)) * ea$log_diag) +
    sum(log1m_rho2) +
    P$log_jac +
    M$log_jac +
    sum((p - seq_len(p)) * (P$log_diag + M$log_diag)) +
    sum(outer(log1m_rho2, log1m_rho2, "+")[low] / 2 - log(2))
  sd_ea <- exp(theta[idx$sd_ea])
  sd_r <- exp(theta[idx$sd_r])
  cor_ea <- tcrossprod(ea$L)
  list(
    sd = c(sd_ea, sd_r),
    cor = c(cor_ea[lower.tri(cor_ea)], rho, intra[low], inter[low]),
    sigma_ea = cor_ea * tcrossprod(sd_ea),
    sigma_r = dyad_cov(sd_r, intra, inter),
    log_jac = log_jac
  )
}

# The 9 SDs and 24 correlations, in the order of unpack()
low_ea <- which(lower.tri(diag(2 * p)), arr.ind = TRUE)
low_r <- which(lower.tri(diag(p)), arr.ind = TRUE)
natural_names <- c(
  paste0("sd_", case_names),
  paste0("sd_R", seq_len(p)),
  paste0("r_", case_names[low_ea[, 1]], "_", case_names[low_ea[, 2]]),
  paste0("rho_R", seq_len(p)),
  paste0("intra_R", low_r[, 1], low_r[, 2]),
  paste0("inter_R", low_r[, 1], low_r[, 2])
)
natural <- function(theta) {
  q <- unpack(theta)
  setNames(c(q$sd, q$cor), natural_names)
}

## ----- Log-posterior ---------------------------------------------------------

# The paper's default priors, as it describes them: each SD ~ t(4, s, s)
# truncated at zero, with s the sample SD of its variable, and each correlation
# r with (r + 1) / 2 ~ Beta(1.5, 1.5). On Omega_EA and on Sigma_R the Beta
# densities join into one, truncated to the positive definite matrices.
log_prior <- function(q, sd_y) {
  loc <- c(rep(sd_y, each = 2), sd_y)
  sum(dt((q$sd - loc) / loc, df = 4, log = TRUE) - log(loc)) -
    length(loc) * pt(1, df = 4, log.p = TRUE) +
    sum(dbeta((q$cor + 1) / 2, 1.5, 1.5, log = TRUE) - log(2))
}

logpost <- function(theta, stats, sd_y) {
  q <- unpack(theta)
  srm_loglik(q$sigma_ea, q$sigma_r, stats) + log_prior(q, sd_y) + q$log_jac
}

## ----- Checks ----------------------------------------------------------------

# The log-likelihood of one group straight from the SRM, with the covariance
# of every pair of ratings built from its definition
dense_loglik <- function(y, sigma_ea, sigma_r) {
  n <- dim(y)[1]
  obs <- expand.grid(k = seq_len(p), i = seq_len(n), j = seq_len(n))
  obs <- obs[obs$i != obs$j, ]
  e <- 2 * obs$k - 1
  a <- 2 * obs$k
  N <- nrow(obs)
  V <- matrix(0, N, N)
  for (r in seq_len(N)) {
    for (s in seq_len(N)) {
      i <- obs$i[r]
      j <- obs$j[r]
      i2 <- obs$i[s]
      j2 <- obs$j[s]
      case <- c(i == i2, j == j2, i == j2, j == i2) *
        c(
          sigma_ea[e[r], e[s]],
          sigma_ea[a[r], a[s]],
          sigma_ea[e[r], a[s]],
          sigma_ea[a[r], e[s]]
        )
      dyad <- c(i == i2 && j == j2, i == j2 && j == i2) *
        c(sigma_r[e[r], e[s]], sigma_r[e[r], a[s]])
      V[r, s] <- sum(case) + sum(dyad)
    }
  }
  x <- y[cbind(obs$i, obs$j, obs$k)]
  chol_v <- chol(V)
  z <- backsolve(chol_v, x, transpose = TRUE)
  -0.5 * (N * log(2 * pi) + 2 * sum(log(diag(chol_v))) + sum(z^2))
}

y_check <- simulate_srm(1, 4)
stopifnot(all.equal(
  srm_loglik(sigma_ea_true, sigma_r_true, srm_stats(y_check)),
  dense_loglik(y_check[1, , , ], sigma_ea_true, sigma_r_true)
))

# The log-Jacobian against central differences of natural()
theta_check <- rnorm(length(theta_names), sd = 0.5)
jac <- sapply(seq_along(theta_check), function(j) {
  h <- replace(numeric(length(theta_check)), j, 1e-5)
  (natural(theta_check + h) - natural(theta_check - h)) / 2e-5
})
stopifnot(all.equal(
  unpack(theta_check)$log_jac,
  determinant(jac)$modulus[[1]],
  tolerance = 1e-6
))

## =============================================================================
## Fit
## =============================================================================

# The SDs start at a third of each variable's variance, and the correlations at
# zero
fit_design <- function(G, n) {
  y <- simulate_srm(G, n)
  stats <- srm_stats(y)
  sd_y <- apply(y, 4, sd, na.rm = TRUE)
  init <- setNames(numeric(length(theta_names)), theta_names)
  init[c(idx$sd_ea, idx$sd_r)] <- log(c(rep(sd_y, each = 2), sd_y) / sqrt(3))
  fit <- inla_vanilla(logpost, init, stats = stats, sd_y = sd_y)
  draws <- sampling(fit, nsamp = nsamp, fun = natural)
  list(G = G, n = n, stats = stats, sd_y = sd_y, fit = fit, draws = draws)
}

results <- Map(fit_design, designs$G, designs$n)
for (res in results) {
  cat(sprintf("\nG = %d groups of n = %d\n", res$G, res$n))
  print(res$fit)
}

## =============================================================================
## Check against Stan
## =============================================================================

# The same posterior: Stan's own transforms add the Jacobians of the SDs, rho
# and the CPCs, and the model block adds the rest of unpack()'s
srm_stan <- "
data {
  int<lower=1> P;
  int<lower=2> n;
  int<lower=1> G;
  matrix[2 * P, 2 * P] W;
  matrix[2 * P, 2 * P] C;
  matrix[2 * P, 2 * P] TT;
  vector<lower=0>[P] sd_y;
}
transformed data {
  int D = 2 * P;
  matrix[D, D] S = rep_matrix(0, D, D);
  vector[D] loc_ea;
  for (k in 1:P) {
    S[2 * k - 1, 2 * k] = 1;
    S[2 * k, 2 * k - 1] = 1;
    loc_ea[2 * k - 1] = sd_y[k];
    loc_ea[2 * k] = sd_y[k];
  }
}
parameters {
  vector<lower=0>[D] sd_ea;
  cholesky_factor_corr[D] L_ea;
  vector<lower=0>[P] sd_r;
  vector<lower=-1, upper=1>[P] rho;
  cholesky_factor_corr[P] L_p;
  cholesky_factor_corr[P] L_m;
}
transformed parameters {
  matrix[D, D] cor_ea = multiply_lower_tri_self_transpose(L_ea);
  matrix[P, P] intra;
  matrix[P, P] inter;
  {
    matrix[P, P] om_p = multiply_lower_tri_self_transpose(L_p);
    matrix[P, P] om_m = multiply_lower_tri_self_transpose(L_m);
    for (k in 1:P) {
      for (l in 1:P) {
        real a = sqrt((1 + rho[k]) * (1 + rho[l]));
        real b = sqrt((1 - rho[k]) * (1 - rho[l]));
        intra[k, l] = (a * om_p[k, l] + b * om_m[k, l]) / 2;
        inter[k, l] = (a * om_p[k, l] - b * om_m[k, l]) / 2;
      }
    }
  }
}
model {
  matrix[D, D] sigma_ea = quad_form_diag(cor_ea, sd_ea);
  matrix[D, D] sigma_r;
  for (k in 1:P) {
    for (l in 1:P) {
      real v = sd_r[k] * sd_r[l];
      sigma_r[2 * k - 1, 2 * l - 1] = v * intra[k, l];
      sigma_r[2 * k, 2 * l] = v * intra[k, l];
      sigma_r[2 * k - 1, 2 * l] = v * inter[k, l];
      sigma_r[2 * k, 2 * l - 1] = v * inter[k, l];
    }
  }
  for (i in 1:D) target += (D - i) * log(L_ea[i, i]);
  for (i in 1:P) target += (P - i) * (log(L_p[i, i]) + log(L_m[i, i]));
  for (k in 2:P) {
    for (l in 1:(k - 1)) {
      target += 0.5 * (log1m(square(rho[k])) + log1m(square(rho[l]))) - log(2);
    }
  }
  sd_ea ~ student_t(4, loc_ea, loc_ea);
  sd_r ~ student_t(4, sd_y, sd_y);
  for (k in 2:D) {
    for (l in 1:(k - 1)) target += beta_lpdf((cor_ea[k, l] + 1) / 2 | 1.5, 1.5);
  }
  target += beta_lpdf((rho + 1) / 2 | 1.5, 1.5);
  for (k in 2:P) {
    for (l in 1:(k - 1)) {
      target += beta_lpdf((intra[k, l] + 1) / 2 | 1.5, 1.5);
      target += beta_lpdf((inter[k, l] + 1) / 2 | 1.5, 1.5);
    }
  }
  {
    matrix[D, D] Q = inverse_spd(sigma_r);
    matrix[D, D] QS = 0.5 * (Q * S + S * Q);
    matrix[D, D] K0 = inverse_spd(sigma_ea) + (n - 1) * Q - QS;
    matrix[D, D] K1 = K0 + n * QS;
    matrix[D, D] K0_inv = inverse_spd(K0);
    matrix[D, D] K1_inv = inverse_spd(K1);
    real log_det = 0.5 * n * (n - 1) * log_determinant(sigma_r) +
      n * log_determinant(sigma_ea) + log_determinant(K1) +
      (n - 1) * log_determinant(K0);
    real quad = sum(Q .* W) - sum(K0_inv .* (Q * C * Q)) -
      sum((K1_inv - K0_inv) .* (Q * TT * Q)) / n;
    target += -0.5 * (G * (P * n * (n - 1) * log(2 * pi()) + log_det) + quad);
  }
}
"

# Stan's draws of the 9 SDs and 24 correlations, in the order of natural()
stan_natural <- function(stan_fit) {
  ea <- sprintf("cor_ea[%d,%d]", low_ea[, 1], low_ea[, 2])
  r <- sprintf("[%d,%d]", low_r[, 1], low_r[, 2])
  pars <- c(
    sprintf("sd_ea[%d]", seq_len(2 * p)),
    sprintf("sd_r[%d]", seq_len(p)),
    ea,
    sprintf("rho[%d]", seq_len(p)),
    paste0("intra", r),
    paste0("inter", r)
  )
  draws <- as.matrix(stan_fit)[, pars]
  colnames(draws) <- natural_names
  draws
}

if (run_stan) {
  stan_model <- rstan::stan_model(model_code = srm_stan)
  for (i in seq_along(results)) {
    res <- results[[i]]
    stan_fit <- rstan::sampling(
      stan_model,
      data = c(list(P = p, sd_y = res$sd_y), res$stats),
      chains = 4,
      iter = stan_iter,
      warmup = 1000,
      seed = 1,
      refresh = 0
    )
    stan_draws <- stan_natural(stan_fit)
    draws <- res$draws

    # Gaps in posterior SDs of the Stan draws: under 0.1 is a close match
    compare <- data.frame(
      truth = c(
        sd_ea_true,
        sd_r_true,
        cor_ea_true[low_ea],
        rho_true,
        intra_true[low_r],
        inter_true[low_r]
      ),
      mean = colMeans(draws),
      stan_mean = colMeans(stan_draws),
      sd = apply(draws, 2, sd),
      stan_sd = apply(stan_draws, 2, sd),
      row.names = natural_names
    )
    compare$mean_gap <- (compare$mean - compare$stan_mean) / compare$stan_sd
    compare$sd_ratio <- compare$sd / compare$stan_sd
    q_gap <- function(prob) {
      (apply(draws, 2, quantile, prob) -
        apply(stan_draws, 2, quantile, prob)) /
        compare$stan_sd
    }
    compare$q025_gap <- q_gap(0.025)
    compare$q975_gap <- q_gap(0.975)
    results[[i]]$stan_fit <- stan_fit
    results[[i]]$stan_draws <- stan_draws
    results[[i]]$compare <- compare
    results[[i]]$stan_time <- sum(rstan::get_elapsed_time(stan_fit))

    cat(sprintf("\nG = %d groups of n = %d\n", res$G, res$n))
    print(round(compare, 3))
    cat(sprintf(
      "inla_vanilla %.1f s, Stan %.1f s (4 chains, without compiling)\n",
      res$fit$timing[["total"]],
      results[[i]]$stan_time
    ))
  }
}

## =============================================================================
## Plot the smallest design
## =============================================================================

# Each SD and reciprocity in the smallest design, with Stan's dashed when it
# ran and the true value in grey
res <- results[[1]]
show <- c(natural_names[1:9], paste0("rho_R", seq_len(p)))
truth <- c(sd_ea_true, sd_r_true, rho_true)
op <- par(mfrow = c(3, 4), mar = c(4, 3, 2, 1))
for (j in seq_along(show)) {
  dens <- density(res$draws[, show[j]])
  stan_dens <- if (run_stan) density(res$stan_draws[, show[j]])
  plot(
    dens,
    main = show[j],
    xlab = "",
    ylim = c(0, max(dens$y, stan_dens$y))
  )
  if (run_stan) {
    lines(stan_dens, lty = 2)
  }
  abline(v = truth[j], col = "grey")
}
par(op)
