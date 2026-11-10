# Gradients and Hessians with RTMB

[RTMB](https://github.com/kaskr/RTMB) (Kristensen, 2026) differentiates
a function written in R. `RTMB::MakeADFun()` runs the function once on
RTMB’s own type of number and records every operation. `obj$fn`,
`obj$gr` and `obj$he` then replay that record, the *tape*, to give the
value, the exact gradient and the exact Hessian at any parameter value.
RTMB is the R interface to TMB (Kristensen et al., 2016), whose models
are written as C++ templates. RTMB needs neither C++ nor a compiler.

This article writes the structural equation model of [Speeding up with
Rcpp](http://haziqj.ml/inlavanilla/articles/rcpp.md), with 39
parameters, for RTMB, and compares the fit with the R and C++ fits of
that article, for speed and for accuracy, with and without RTMB’s
Hessian. The model, its parameterisation and its priors are described
there.

``` r

library(inlavanilla)
data(PoliticalDemocracy, package = "lavaan")
```

The data, the layout of the parameters, the starting values and the R
log-posterior `logpost_r()` are those of the Rcpp article, and so are
the C++ log-posterior `logpost_cpp()` and its gradient `grad_cpp()`,
written by hand from the closed-form gradient of a structural equation
model (Neudecker & Satorra, 1991). Expand the code to see them.

R code from the Rcpp article

``` r

ov <- c("x1", "x2", "x3", paste0("y", 1:8))
Y <- as.matrix(PoliticalDemocracy[, ov])

cov_pairs <- rbind(
  c("y1", "y5"),
  c("y2", "y4"),
  c("y2", "y6"),
  c("y3", "y7"),
  c("y4", "y8"),
  c("y6", "y8")
)

lambda_fixed <- matrix(0, 11, 3)
lambda_fixed[cbind(c(1, 4, 8), 1:3)] <- 1
lambda_at <- cbind(c(2, 3, 5, 6, 7, 9, 10, 11), rep(1:3, c(2, 3, 3)))
lambda_par <- c(12, 13, 14, 15, 16, 14, 15, 16)
beta_at <- cbind(c(2, 3, 3), c(1, 1, 2))
cov_at <- matrix(match(cov_pairs, ov), ncol = 2)

implied_cov <- function(theta) {
  Lambda <- lambda_fixed
  Lambda[lambda_at] <- theta[lambda_par]
  B <- matrix(0, 3, 3)
  B[beta_at] <- theta[17:19]
  psi <- exp(theta[20:22])
  v <- exp(theta[23:33])
  Theta <- diag(v)
  Theta[cov_at] <- tanh(theta[34:39]) * sqrt(v[cov_at[, 1]] * v[cov_at[, 2]])
  Theta[cov_at[, 2:1]] <- Theta[cov_at]
  LA <- Lambda %*% solve(diag(3) - B)
  LA %*% (psi * t(LA)) + Theta # LA Psi LA' + Theta
}

loglik_r <- function(theta, Y) {
  Sigma <- implied_cov(theta)
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) {
    return(rep(-Inf, nrow(Y))) # Sigma is not positive definite
  }
  Z <- backsolve(R, t(Y) - theta[1:11], transpose = TRUE)
  -0.5 * ncol(Y) * log(2 * pi) - sum(log(diag(R))) - 0.5 * colSums(Z^2)
}

logprior_r <- function(theta) {
  logv <- theta[20:33]
  sd <- exp(logv / 2)
  rho <- tanh(theta[34:39])
  sum(dnorm(theta[1:11], 0, 32, log = TRUE)) +
    sum(dnorm(theta[12:19], 0, 10, log = TRUE)) +
    sum(dgamma(sd, shape = 1, rate = 0.5, log = TRUE) + log(sd / 2)) +
    sum(dbeta((rho + 1) / 2, 1, 1, log = TRUE) - log(2) + log1p(-rho^2))
}

logpost_r <- function(theta, Y) {
  sum(loglik_r(theta, Y)) + logprior_r(theta)
}

init <- c(
  colMeans(Y),
  rep(1, 5),
  rep(0, 3),
  rep(0, 3),
  log(apply(Y, 2, var) / 2),
  rep(0, 6)
)
names(init) <- c(
  paste0("nu_", ov),
  c("lam_x2", "lam_x3", "a", "b", "c"),
  c("beta_dem60_ind60", "beta_dem65_ind60", "beta_dem65_dem60"),
  paste0("logv_", c("ind60", "dem60", "dem65")),
  paste0("logv_", ov),
  paste0("z_", cov_pairs[, 1], "_", cov_pairs[, 2])
)
```

C++ code from the Rcpp article

``` r

Rcpp::sourceCpp(
  code = r"(
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// theta: intercepts 0-10, loadings 11-15, regressions 16-18, latent log
// variances 19-21, residual log variances 22-32, residual correlations
// (atanh scale) 33-38. Indicators in the order x1, x2, x3, y1, ..., y8.

// Free loadings: row and column of Lambda, and the element of theta
const int lambda_row[8] = {1, 2, 4, 5, 6, 8, 9, 10};
const int lambda_col[8] = {0, 0, 1, 1, 1, 2, 2, 2};
const int lambda_par[8] = {11, 12, 13, 14, 15, 13, 14, 15};
// Regressions: dem60 ~ ind60, dem65 ~ ind60, dem65 ~ dem60
const int beta_row[3] = {1, 2, 2};
const int beta_col[3] = {0, 0, 1};
// Residual covariances: y1-y5, y2-y4, y2-y6, y3-y7, y4-y8, y6-y8
const int cov_i[6] = {3, 4, 4, 5, 6, 8};
const int cov_j[6] = {7, 6, 8, 9, 10, 10};

struct Sem {
  VectorXd nu, psi, v, rho;
  MatrixXd Lambda, A, M, Theta, Sigma;  // A = (I - B)^-1, M = A Psi A'
};

Sem sem_matrices(const VectorXd& theta) {
  Sem s;
  s.nu = theta.head(11);
  s.Lambda = MatrixXd::Zero(11, 3);
  s.Lambda(0, 0) = s.Lambda(3, 1) = s.Lambda(7, 2) = 1;
  for (int k = 0; k < 8; k++) {
    s.Lambda(lambda_row[k], lambda_col[k]) = theta[lambda_par[k]];
  }
  MatrixXd B = MatrixXd::Zero(3, 3);
  for (int k = 0; k < 3; k++) B(beta_row[k], beta_col[k]) = theta[16 + k];
  s.A = (MatrixXd::Identity(3, 3) - B).inverse();
  s.psi = theta.segment(19, 3).array().exp();
  s.v = theta.segment(22, 11).array().exp();
  s.rho = theta.segment(33, 6).array().tanh();
  s.Theta = s.v.asDiagonal();
  for (int k = 0; k < 6; k++) {
    int i = cov_i[k], j = cov_j[k];
    s.Theta(i, j) = s.Theta(j, i) = s.rho[k] * std::sqrt(s.v[i] * s.v[j]);
  }
  s.M = s.A * s.psi.asDiagonal() * s.A.transpose();
  s.Sigma = s.Lambda * s.M * s.Lambda.transpose() + s.Theta;
  return s;
}

// [[Rcpp::export]]
VectorXd loglik_cpp(const Map<VectorXd> theta, const Map<MatrixXd> Y) {
  Sem s = sem_matrices(theta);
  Eigen::LLT<MatrixXd> llt(s.Sigma);
  if (llt.info() != Eigen::Success) {  // Sigma is not positive definite
    return VectorXd::Constant(Y.rows(), R_NegInf);
  }
  MatrixXd D = Y.rowwise() - s.nu.transpose();
  MatrixXd Z = llt.matrixL().solve(D.transpose());
  // matrixLLT() holds L, whose diagonal gives half the log-determinant
  double half_logdet = llt.matrixLLT().diagonal().array().log().sum();
  double c = -0.5 * Y.cols() * std::log(2 * M_PI) - half_logdet;
  return (c - 0.5 * Z.colwise().squaredNorm().array()).matrix().transpose();
}

// [[Rcpp::export]]
double logprior_cpp(const Map<VectorXd> theta) {
  double lp = 0;
  for (int i = 0; i < 11; i++) lp += R::dnorm(theta[i], 0, 32, true);
  for (int i = 11; i < 19; i++) lp += R::dnorm(theta[i], 0, 10, true);
  for (int i = 19; i < 33; i++) {
    double sd = std::exp(theta[i] / 2);
    // R::dgamma() takes the scale, 1 / rate
    lp += R::dgamma(sd, 1, 1 / 0.5, true) + std::log(sd / 2);
  }
  for (int i = 33; i < 39; i++) {
    double rho = std::tanh(theta[i]);
    lp += R::dbeta((rho + 1) / 2, 1, 1, true) - std::log(2.0) +
      std::log1p(-rho * rho);
  }
  return lp;
}

// [[Rcpp::export]]
double logpost_cpp(const Map<VectorXd> theta, const Map<MatrixXd> Y) {
  return loglik_cpp(theta, Y).sum() + logprior_cpp(theta);
}
// Gradient of logpost_cpp() with respect to theta
// [[Rcpp::export]]
VectorXd grad_cpp(const Map<VectorXd> theta, const Map<MatrixXd> Y) {
  Sem s = sem_matrices(theta);
  Eigen::LLT<MatrixXd> llt(s.Sigma);
  if (llt.info() != Eigen::Success) return VectorXd::Constant(39, NA_REAL);
  VectorXd g = VectorXd::Zero(39);

  // Log-likelihood, through G = d loglik / d Sigma
  int n = Y.rows();
  MatrixXd Omega = llt.solve(MatrixXd::Identity(11, 11));
  MatrixXd OD = Omega * (Y.rowwise() - s.nu.transpose()).transpose();
  MatrixXd G = 0.5 * (OD * OD.transpose() - n * Omega);

  g.head(11) = OD.rowwise().sum();
  MatrixXd dLambda = 2 * G * s.Lambda * s.M;
  for (int k = 0; k < 8; k++) {
    g[lambda_par[k]] += dLambda(lambda_row[k], lambda_col[k]);
  }
  MatrixXd LGL = s.Lambda.transpose() * G * s.Lambda;
  MatrixXd dB = 2 * s.A.transpose() * LGL * s.M;
  for (int k = 0; k < 3; k++) g[16 + k] = dB(beta_row[k], beta_col[k]);
  MatrixXd dPsi = s.A.transpose() * LGL * s.A;
  for (int k = 0; k < 3; k++) g[19 + k] = dPsi(k, k) * s.psi[k];
  for (int i = 0; i < 11; i++) g[22 + i] = G(i, i) * s.v[i];
  for (int k = 0; k < 6; k++) {
    int i = cov_i[k], j = cov_j[k];
    g[22 + i] += G(i, j) * s.Theta(i, j);
    g[22 + j] += G(i, j) * s.Theta(i, j);
    g[33 + k] = 2 * G(i, j) * (1 - s.rho[k] * s.rho[k]) *
      std::sqrt(s.v[i] * s.v[j]);
  }

  // Log-prior, with the Jacobians
  for (int i = 0; i < 11; i++) g[i] -= theta[i] / (32.0 * 32.0);
  for (int i = 11; i < 19; i++) g[i] -= theta[i] / (10.0 * 10.0);
  for (int i = 19; i < 33; i++) g[i] += 0.5 - std::exp(theta[i] / 2) / 4;
  for (int k = 0; k < 6; k++) g[33 + k] -= 2 * s.rho[k];
  return g;
}
)"
)
```

## The log-posterior in RTMB

`logpost_rtmb()` is `logpost_r()` with three changes:

- The densities, [`diag()`](https://rdrr.io/r/base/diag.html) and
  [`solve()`](https://rdrr.io/r/base/solve.html) are RTMB’s, which also
  take its numbers. We call them with `RTMB::` rather than attaching the
  package: [`library(RTMB)`](https://github.com/kaskr/RTMB) masks the
  base and stats functions of the same names, and slows ordinary R code
  that calls them, `logpost_r()` among it.
- `"[<-" <- RTMB::ADoverload("[<-")` lets RTMB’s numbers be assigned
  into the ordinary matrices \\\Lambda\\, \\B\\ and \\\Theta\\.
- [`RTMB::dmvnorm()`](https://rdrr.io/pkg/RTMB/man/MVgauss.html)
  replaces the Cholesky factor and its
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html). The tape keeps
  the path taken when it was recorded, so a branch on the parameters,
  such as returning `-Inf` when \\\Sigma\\ is not positive definite,
  cannot be recorded. Where \\\Sigma\\ is not positive definite the tape
  gives `NaN`, which
  [`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
  treats as `-Inf`.

``` r

logpost_rtmb <- function(theta, Y) {
  "[<-" <- RTMB::ADoverload("[<-")
  Lambda <- lambda_fixed
  Lambda[lambda_at] <- theta[lambda_par]
  B <- matrix(0, 3, 3)
  B[beta_at] <- theta[17:19]
  psi <- exp(theta[20:22])
  v <- exp(theta[23:33])
  Theta <- RTMB::diag(v)
  Theta[cov_at] <- tanh(theta[34:39]) * sqrt(v[cov_at[, 1]] * v[cov_at[, 2]])
  Theta[cov_at[, 2:1]] <- Theta[cov_at]
  LA <- Lambda %*% RTMB::solve(diag(3) - B)
  Sigma <- LA %*% (psi * t(LA)) + Theta
  loglik <- RTMB::dmvnorm(t(t(Y) - theta[1:11]), Sigma = Sigma, log = TRUE)

  sd <- exp(theta[20:33] / 2)
  rho <- tanh(theta[34:39])
  logprior <- sum(RTMB::dnorm(theta[1:11], 0, 32, log = TRUE)) +
    sum(RTMB::dnorm(theta[12:19], 0, 10, log = TRUE)) +
    sum(RTMB::dgamma(sd, shape = 1, rate = 0.5, log = TRUE) + log(sd / 2)) +
    sum(RTMB::dbeta((rho + 1) / 2, 1, 1, log = TRUE) - log(2) + log1p(-rho^2))
  sum(loglik) + logprior
}
```

`RTMB::MakeADFun()` takes a function of a list of parameters. TMB’s
objectives are usually negative log-likelihoods, to be minimised, but
the sign matters only when random effects are integrated out. There are
none here, so the function returns the log-posterior itself, and
`obj$fn`, `obj$gr` and `obj$he` go to
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
as they are.

``` r

obj <- RTMB::MakeADFun(
  function(p) logpost_rtmb(p$theta, Y),
  list(theta = unname(init)),
  silent = TRUE
)
```

`obj$he` works out the Hessian from the recorded function at each call.
A *Hessian tape* is faster:
[`RTMB::MakeTape()`](https://rdrr.io/pkg/RTMB/man/Tape.html) records the
function, which here takes the parameter vector itself, `$jacfun()`
records the computation of its gradient, and a second `$jacfun()` that
of its Hessian, so that each call replays the Hessian’s own computation.

``` r

hess_tape <- RTMB::MakeTape(
  function(theta) logpost_rtmb(theta, Y),
  unname(init)
)$jacfun()$jacfun()
```

At 200 points around the starting values, the value matches the R
log-posterior, and the gradient the hand-written one, to rounding error:

``` r

set.seed(1)
points <- replicate(200, init + rnorm(39, sd = 0.05), simplify = FALSE)
value_error <- function(th) abs(obj$fn(th) - logpost_r(th, Y))
grad_error <- function(th) max(abs(obj$gr(th) - grad_cpp(th, Y)))
c(
  value = max(vapply(points, value_error, 0)),
  gradient = max(vapply(points, grad_error, 0))
)
#>        value     gradient 
#> 1.364242e-12 1.136868e-13
```

## Speed

The time per call, in microseconds, over the same 200 points. (Repeated
calls at one point come back faster than calls at new points, as in a
fit.)

``` r

time_per_call <- function(f, reps = 50) {
  elapsed <- system.time(
    for (r in seq_len(reps)) {
      for (th in points) f(th)
    }
  )[["elapsed"]]
  1e6 * elapsed / (reps * length(points))
}
us <- rbind(
  logpost = c(
    R = time_per_call(function(th) logpost_r(th, Y)),
    RTMB = time_per_call(obj$fn),
    `C++` = time_per_call(function(th) logpost_cpp(th, Y))
  ),
  gradient = c(
    R = NA,
    RTMB = time_per_call(obj$gr),
    `C++` = time_per_call(function(th) grad_cpp(th, Y))
  ),
  hessian = c(R = NA, RTMB = time_per_call(hess_tape, reps = 1), `C++` = NA)
)
round(us, 1)
#>              R    RTMB C++
#> logpost  108.1   151.3 8.7
#> gradient    NA   251.4 9.3
#> hessian     NA 16780.0  NA
```

Replaying the tape takes 1.4 times as long as a call of `logpost_r()`,
and RTMB’s gradient costs as much as 1.7 calls of its log-posterior. C++
is 12 times faster than R, and the hand-written gradient costs as much
as 1.1 calls of its log-posterior. RTMB’s Hessian tape costs as much as
67 gradients, about one for each of the 39 parameters, and `obj$he` 39.
The five fits, with the Hessian tape in the third:

``` r

fits <- list(
  R = inla_vanilla(logpost_r, init, Y = Y, silent = TRUE),
  RTMB = inla_vanilla(obj$fn, init, grad = obj$gr, silent = TRUE),
  `RTMB and Hessian` = inla_vanilla(
    obj$fn,
    init,
    grad = obj$gr,
    hess = hess_tape,
    silent = TRUE
  ),
  `C++` = inla_vanilla(logpost_cpp, init, Y = Y, silent = TRUE),
  `C++ and gradient` = inla_vanilla(
    logpost_cpp,
    init,
    grad = grad_cpp,
    Y = Y,
    silent = TRUE
  )
)
seconds <- vapply(fits, function(f) f$timing[["total"]], 0)
cbind(t(sapply(fits, function(f) f$counts)), seconds = round(seconds, 2))
#>                  logpost grad hess seconds
#> R                  13582    0    0    1.72
#> RTMB                1291  483    0    0.48
#> RTMB and Hessian     953  330   12    0.57
#> C++                13433    0    0    0.33
#> C++ and gradient    1290  484    0    0.18
```

Without a gradient the fit calls `logpost` 13,582 times. With one it
calls `logpost` 1,291 times and the gradient 483 times. Without a
gradient the count grows with the square of the number of parameters,
and with one in proportion to it. Here a gradient saves time when it
costs less than 25 calls of `logpost`. RTMB’s costs far less than that,
and its fit is 3.6 times as fast as R’s. A hand-written gradient in C++
is the fastest of all.

## Accuracy

Each fit against the one with the hand-written gradient: the skew-normal
marginals (means, SDs and quantiles) and the VB mean in posterior SDs,
and the log marginal likelihood in nats.

``` r

ref <- fits[["C++ and gradient"]]
sds <- ref$summary[, "SD"]
signif(
  t(sapply(fits, function(f) {
    c(
      marginals = max(abs(f$summary - ref$summary) / sds),
      vb_mean = max(abs(f$mean_vb - ref$mean_vb) / sds),
      logml = f$logml - ref$logml
    )
  })),
  2
)
#>                  marginals vb_mean   logml
#> R                  2.9e-03 0.03300 2.2e-02
#> RTMB               9.9e-05 0.00018 7.6e-05
#> RTMB and Hessian   2.0e-04 0.00035 1.6e-04
#> C++                3.0e-03 0.03300 2.2e-02
#> C++ and gradient   0.0e+00 0.00000 0.0e+00
```

RTMB’s fit is the hand-written gradient’s fit, since the two gradients
agree to rounding error, and so is the fit with RTMB’s Hessian. The
skew-normal marginals are nearly the same on every route: they are
scanned around the mode and do not depend on the VB correction. Their
differences, a few ten-thousandths of a posterior SD, come from the
volume correction, a central difference with a gradient and a forward
one otherwise. The gradient changes the VB correction, which then takes
its step from the gradient at 2m + 1 = 79 points instead of from
`logpost` at (m + 1)(m + 2) = 1,640 points (see the section on the VB
correction in
[`?inla_vanilla`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)).
Here the fits without a gradient move the VB mean by up to 0.03
posterior SDs and the log marginal likelihood by 0.02. That is small,
but it carries over to everything built on the VB mean: the log marginal
likelihood, [`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md) and
[`waic()`](http://haziqj.ml/inlavanilla/reference/loo.md), and the
Gaussian marginals of `marginals = "gaussian"`.

## The Hessian

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
takes the Hessian as `hess`, after `grad`, as
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html) does, and uses it in
two places. [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) takes
Newton steps to the mode. And the Hessian at the mode, which gives the
covariance of the Laplace approximation, is RTMB’s rather than central
differences of the gradient. The default marginals do not use it: their
volume correction needs only two curvatures per parameter, and 4m + 1 =
157 gradients give them for every marginal at once, where a Hessian each
would cost far more. The time each stage took, in seconds:

``` r

rtmb_fits <- fits[c("RTMB", "RTMB and Hessian")]
t(sapply(rtmb_fits, function(f) round(f$timing, 2)))
#>                  init laplace   vb marginals diagnostics total
#> RTMB             0.00    0.09 0.06      0.33           0  0.48
#> RTMB and Hessian 0.02    0.18 0.06      0.31           0  0.57
```

The marginals and the VB correction take the same time with and without
it. The mode stage takes a little longer, because a Hessian costs about
as much as m gradients, and the whole fit takes 1.19 times as long. The
gap is smaller where the Hessian is sparse, as when each of many
parameters enters only a few terms of the log-posterior: a Hessian tape
then costs far fewer than m gradients.

The Hessian also finds the mode more precisely.
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html) stops when its steps
become small. The whitened gradient where it stopped in each fit is
about how far it was from the mode, in posterior SDs:

``` r

signif(sapply(fits, function(f) f$diagnostics$optim$grad_mode), 1)
#>                R             RTMB RTMB and Hessian              C++ 
#>            1e-04            8e-05            4e-07            1e-04 
#> C++ and gradient 
#>            1e-04
```

Without the Hessian, [`nlminb()`](https://rdrr.io/r/stats/nlminb.html)
stops about a ten-thousandth of a posterior SD from the mode. With the
Hessian it stops within a millionth. Either is far inside the check of
the mode, which asks for a tenth. A supplied `hess` is itself checked at
the mode against central differences of `grad`, with a warning if they
do not match. With random effects integrated out, TMB gives no Hessian.

## Using several cores

The skew-normal scans can run in parallel with `cores`, which needs the
optional future package. With no random effects, as here, RTMB needs
nothing special:

``` r

cores <- future::availableCores(omit = 1)
fit_par <- inla_vanilla(
  obj$fn, init, grad = obj$gr, width = "tk", cores = cores
)
```

Random effects are the exception. Forked workers, which `cores` uses in
a terminal on macOS or Linux, share the tape. Separate R sessions, which
it uses in RStudio, on Windows and in Positron, whose console refuses to
fork, do not: once `MakeADFun()` integrates random effects out, a copy
of `obj$fn` in another session returns `NaN`, with no error. The fit
checks a worker’s value at the mode, finds the difference, and runs
serially with a message.

To run in parallel there, we build the tape on each worker: the call of
`MakeADFun()` goes in a function, which is called again when the process
is not the one that started the fit.

``` r

build <- function() {
  RTMB::MakeADFun(model, parameters, random = "u", silent = TRUE)
}
obj <- build()
pid <- Sys.getpid()
rebuild <- function() {
  if (Sys.getpid() != pid) {
    obj <<- build()
    pid <<- Sys.getpid()
  }
}
logpost <- function(theta) {
  rebuild()
  -obj$fn(theta)
}
grad <- function(theta) {
  rebuild()
  -obj$gr(theta)
}
fit_par <- inla_vanilla(
  logpost, init, grad = grad, width = "tk", cores = cores
)
```

Whatever `build()` uses, here `model` and `parameters`, is found and
sent to the workers. On a Poisson model with 21 parameters and 60 random
intercepts integrated out, a fit with `width = "tk"` took 32 s serially
and 4.4 s on 12 cores from a terminal. In Positron’s console it took 9.7
s: each worker is a new R session that has to start up, and the fit
chose seven of the 12.

## Choosing a route

In order of the work they take, with the first two about equal:

| Route | What you write | Speed | Accuracy |
|:---|:---|:--:|:--:|
| R | `logpost` in R | ✗ | ✗ |
| RTMB | `logpost` in R, with RTMB’s functions, which give `grad` and `hess` | ✓ | ✓ |
| Rcpp | `logpost` in C++ | ✓ | ✗ |
| R and a gradient | `logpost` and its gradient in R | ✗ | ✓ |
| Rcpp and a gradient | `logpost` and its gradient in C++ | ✓ | ✓ |

Start with R. If RTMB has the functions your log-posterior needs, it
gives the gradient’s accuracy for the same work, and the fit is 3.6
times as fast as R’s. A Hessian finds the mode more precisely for a
little more time, so it is worth passing when it comes free, as from
RTMB, but seldom worth deriving by hand. C++ buys speed, which matters
when each call is expensive, as with a large data set, or when there are
many, as with many parameters or many fits in a simulation study. A
gradient by hand is the most work, since it has to be derived, checked
and kept in step with the log-posterior, and only in C++ does it buy
both.

## References

Kristensen, K. (2026). *RTMB: R bindings for TMB* (Version 2.0)
\[Computer software\]. <https://doi.org/10.32614/CRAN.package.RTMB>

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M.
(2016). TMB: Automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, *70*(5), 1–21.
<https://doi.org/10.18637/jss.v070.i05>

Neudecker, H., & Satorra, A. (1991). Linear structural relations:
Gradient and Hessian of the fitting function. *Statistics & Probability
Letters*, *11*(1), 57–61. <https://doi.org/10.1016/0167-7152(91)90178-T>
