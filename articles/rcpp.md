# Speeding up with Rcpp

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
calls `logpost` many times, and the count grows quickly with the number
of parameters $`m`$. Without a gradient, the VB correction alone needs
$`(m + 1)(2m^2 + 1)`$ calls per sweep. That is 76 calls for 3
parameters, and 121,720 for 39. Two things help: a faster `logpost`, and
a gradient.

This article does both for a structural equation model (SEM) with 39
parameters. We write the log-posterior in plain R first, then in C++
with [Rcpp](https://www.rcpp.org/) and
[RcppEigen](https://cran.r-project.org/package=RcppEigen), and then add
an analytic gradient in C++. The article [Using
TMB](http://haziqj.ml/inlavanilla/articles/tmb.md) fits the same model
with a TMB template, which computes the gradient by automatic
differentiation.

``` r

library(inlavanilla)
library(lavaan)
```

## The model

The `PoliticalDemocracy` data in lavaan come from Bollen (1989). For 75
developing countries they hold three indicators of industrialisation in
1960 (`x1`–`x3`), and the same four indicators of political democracy in
1960 (`y1`–`y4`) and 1965 (`y5`–`y8`). In lavaan syntax, the model is:

``` r

model <- "
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + a*y2 + b*y3 + c*y4
  dem65 =~ y5 + a*y6 + b*y7 + c*y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8
"
```

Three latent variables, industrialisation in 1960 (`ind60`) and
democracy in 1960 and 1965 (`dem60`, `dem65`), are each measured by
their indicators. The first loading of each is fixed at 1, and the
labels `a`, `b` and `c` hold the 1960 and 1965 loadings equal. Democracy
in 1960 depends on industrialisation, and democracy in 1965 on both. The
residuals of the same indicator in the two years are correlated, as are
those of `y2` with `y4` and `y6` with `y8`.

With a mean structure, the 11 indicators of each country are normal with
mean $`\nu`$ and covariance

``` math
\Sigma = \Lambda (I - B)^{-1} \Psi (I - B)^{-\top} \Lambda^\top + \Theta,
```

where $`\nu`$ holds the 11 intercepts, $`\Lambda`$ ($`11 \times 3`$) the
loadings, $`B`$ ($`3 \times 3`$) the regressions among the latent
variables, $`\Psi`$ the diagonal matrix of the variance of `ind60` and
the disturbance variances of `dem60` and `dem65`, and $`\Theta`$
($`11 \times 11`$) the residual covariances. The log-likelihood is the
sum of 75 multivariate normal log-densities.

## Parameters and priors

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
needs every parameter on an unconstrained scale. We use the
parameterisation and default priors of
[INLAvaan](https://inlavaan.haziqj.ml/), so the results can be compared
with `INLAvaan::asem(model, PoliticalDemocracy, meanstructure = TRUE)`:

| Parameters | Number | Unconstrained scale | Prior |
|:---|---:|:---|:---|
| Intercepts $`\nu`$ | 11 | as they are | $`N(0, 32^2)`$ |
| Loadings (`x2`, `x3`, `a`, `b`, `c`) | 5 | as they are | $`N(0, 10^2)`$ |
| Regressions | 3 | as they are | $`N(0, 10^2)`$ |
| Variances (3 latent, 11 residual) | 14 | $`\ell = \log \sigma^2`$ | $`\text{Gamma}(1, \text{rate } 0.5)`$ on $`\sigma`$ |
| Residual covariances | 6 | $`z = \operatorname{atanh} \rho`$ | $`\text{Beta}(1, 1)`$ on $`(\rho + 1)/2`$ |

A residual covariance enters through its correlation $`\rho`$, as
$`\theta_{ij} = \rho \sqrt{\theta_{ii} \theta_{jj}}`$, and both its
parameter $`z`$ and its prior refer to $`\rho`$.

Two of the priors are on a different scale from the parameter, so
`logpost` needs their log-Jacobians:

- A variance’s prior is on its standard deviation
  $`\sigma = e^{\ell / 2}`$, and $`d\sigma / d\ell = \sigma / 2`$. Its
  term in the log-posterior is $`\log p(\sigma) + \log(\sigma / 2)`$.
- A correlation’s prior is on $`\rho = \tanh z`$, and
  $`d\rho / dz = 1 - \rho^2`$. The Beta(1, 1) prior on $`(\rho + 1)/2`$
  is the uniform density $`1/2`$ on $`(-1, 1)`$, so the term is
  $`\log(1/2) + \log(1 - \rho^2)`$.

The normal priors are on the parameters themselves and need none.

## The log-posterior in R

The data, ordered as the rows of $`\Lambda`$, and the names of the 39
parameters. `logv_` marks a log variance and `z_` a residual correlation
on the atanh scale.

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

parnames <- c(
  paste0("nu_", ov), # 1-11
  c("lam_x2", "lam_x3", "a", "b", "c"), # 12-16
  c("beta_dem60_ind60", "beta_dem65_ind60", "beta_dem65_dem60"), # 17-19
  paste0("logv_", c("ind60", "dem60", "dem65")), # 20-22
  paste0("logv_", ov), # 23-33
  paste0("z_", cov_pairs[, 1], "_", cov_pairs[, 2]) # 34-39
)
```

The free entries of $`\Lambda`$, $`B`$ and $`\Theta`$, as matrix
indices, and the element of the parameter vector that fills each. The
loadings `a`, `b` and `c` fill two entries each.

``` r

lambda_fixed <- matrix(0, 11, 3)
lambda_fixed[cbind(c(1, 4, 8), 1:3)] <- 1
lambda_at <- cbind(c(2, 3, 5, 6, 7, 9, 10, 11), rep(1:3, c(2, 3, 3)))
lambda_par <- c(12, 13, 14, 15, 16, 14, 15, 16)
beta_at <- cbind(c(2, 3, 3), c(1, 1, 2))
cov_at <- matrix(match(cov_pairs, ov), ncol = 2)
```

The model-implied covariance matrix, the 75 pointwise log-likelihoods,
and the log-prior with its Jacobians. The log-posterior is their total.
Keeping the pointwise log-likelihood as a function of its own makes it
available for leave-one-out cross-validation.

``` r

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
```

Starting values: the sample means for the intercepts, loadings of 1, no
regressions, half of each sample variance for the residual variances,
and no residual correlations.

``` r

init <- c(
  colMeans(Y),
  rep(1, 5),
  rep(0, 3),
  rep(0, 3),
  log(apply(Y, 2, var) / 2),
  rep(0, 6)
)
names(init) <- parnames
```

### A check against lavaan

At lavaan’s maximum likelihood estimates, our log-likelihood should
equal lavaan’s. We map the estimates to the unconstrained scale:

``` r

fit_ml <- sem(model, PoliticalDemocracy, meanstructure = TRUE)
est <- coef(fit_ml)
est <- est[!duplicated(names(est))] # a, b and c appear twice

lv <- c("ind60", "dem60", "dem65")
var_ml <- est[paste0(c(lv, ov), "~~", c(lv, ov))]
sd_ml <- sqrt(var_ml[-(1:3)])
names(sd_ml) <- ov
cov_ml <- est[paste0(cov_pairs[, 1], "~~", cov_pairs[, 2])]

theta_ml <- c(
  est[paste0(ov, "~1")],
  est[c("ind60=~x2", "ind60=~x3", "a", "b", "c")],
  est[c("dem60~ind60", "dem65~ind60", "dem65~dem60")],
  log(var_ml),
  atanh(cov_ml / (sd_ml[cov_pairs[, 1]] * sd_ml[cov_pairs[, 2]]))
)
names(theta_ml) <- parnames

c(ours = sum(loglik_r(theta_ml, Y)), lavaan = logLik(fit_ml))
#>      ours    lavaan 
#> -1548.818 -1548.818
```

### How long R takes

A helper to time one call, in microseconds:

``` r

time_per_call <- function(f, n) {
  elapsed <- system.time(
    for (i in seq_len(n)) {
      f()
    }
  )[["elapsed"]]
  1e6 * elapsed / n
}
us_r <- time_per_call(function() logpost_r(init, Y), n = 1000)
us_r
#> [1] 139
```

With `vb = FALSE`,
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
skips the VB correction and runs only the mode search and the
skew-normal marginals:

``` r

fit_r <- inla_vanilla(logpost_r, init, Y = Y, vb = FALSE)
#> ✔ Found the posterior mode in 83 iterations [867ms]
#> ✔ Fitted 39 skew-normal marginals [1.0s]
fit_r$counts
#> logpost    grad 
#>   14295       0
```

These two stages take 14,295 calls. The VB correction would add far
more: it evaluates the gradient at $`2m^2 + 1 = 3{,}043`$ points per
sweep, and without `grad` each gradient costs $`m + 1 = 40`$ calls, over
a few sweeps. We come back to the cost of a full fit below.

## The log-posterior in C++

The C++ version mirrors the R one. RcppEigen provides the matrix algebra
and converts R’s numeric vectors and matrices to Eigen’s `VectorXd` and
`MatrixXd` (an `Eigen::Map` reuses R’s memory instead of copying). The
C++ file starts with the layout of the parameter vector, now with
indices from 0, and a function that builds the model matrices.

``` r

cpp_file <- file.path(tempdir(), "sem.cpp")
```

``` cpp
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
```

### The gradient

The gradient follows from the chain rule, and costs about as much as the
log-likelihood itself. Write $`D`$ for the $`n \times 11`$ matrix of
data minus $`\nu`$ and $`\Omega = \Sigma^{-1}`$. The derivative of the
log-likelihood is $`\Omega D^\top \mathbf{1}`$ with respect to $`\nu`$,
and

``` math
G = \tfrac{1}{2} \left( \Omega D^\top D \, \Omega - n \Omega \right)
```

with respect to $`\Sigma`$. Through
$`\Sigma = \Lambda M \Lambda^\top + \Theta`$, with $`A = (I - B)^{-1}`$
and $`M = A \Psi A^\top`$, the derivatives with respect to the model
matrices are

``` math
\frac{\partial \ell}{\partial \Lambda} = 2 G \Lambda M, \qquad
\frac{\partial \ell}{\partial B} = 2 A^\top \Lambda^\top G \Lambda M, \qquad
\frac{\partial \ell}{\partial \Psi} = A^\top \Lambda^\top G \Lambda A,
```

and $`G_{ii}`$ for a residual variance $`\theta_{ii}`$ and $`2 G_{ij}`$
for a residual covariance. A loading held equal across years sums its
two entries. On the unconstrained scale, a log variance multiplies its
derivative by the variance. A residual covariance
$`\theta_{ij} = \rho \sqrt{\theta_{ii} \theta_{jj}}`$ moves with $`z`$,
through
$`d\theta_{ij} / dz = (1 - \rho^2) \sqrt{\theta_{ii} \theta_{jj}}`$, and
with both log variances, through
$`d\theta_{ij} / d\ell_i = \theta_{ij} / 2`$.

The log-prior adds $`-\vartheta / s^2`$ for each normal prior,
$`1/2 - \sigma / 4`$ for each log variance and $`-2\rho`$ for each
correlation. We append the gradient to the same file:

``` cpp
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
```

### Compiling and checking

[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html)
compiles the file and defines the four exported functions in R.
(`Rcpp::sourceCpp(code = ...)` takes the same code as a string.)

``` r

Rcpp::sourceCpp(cpp_file)
```

The C++ log-posterior should match the R one to rounding error, at the
starting values, at lavaan’s estimates and at some random points:

``` r

set.seed(1)
points <- c(
  list(init, theta_ml),
  replicate(4, init + rnorm(39, sd = 0.3), simplify = FALSE)
)
vapply(points, function(th) logpost_cpp(th, Y) - logpost_r(th, Y), 0)
#> [1] -4.547474e-13  2.273737e-13  0.000000e+00 -2.273737e-13  0.000000e+00
#> [6] -4.547474e-13
max(abs(loglik_cpp(theta_ml, Y) - loglik_r(theta_ml, Y)))
#> [1] 7.105427e-15
```

The gradient should match central differences of the R log-posterior, up
to their own error, about $`10^{-8}`$ relative:

``` r

fd_grad <- function(f, x, h = 1e-5) {
  vapply(
    seq_along(x),
    function(j) {
      e <- replace(numeric(length(x)), j, h)
      (f(x + e) - f(x - e)) / (2 * h)
    },
    0
  )
}
vapply(
  points,
  function(th) {
    g_fd <- fd_grad(function(x) logpost_r(x, Y), th)
    max(abs(grad_cpp(th, Y) - g_fd) / (1 + abs(g_fd)))
  },
  0
)
#> [1] 1.461419e-08 1.228863e-08 9.214642e-09 7.379461e-09 3.295428e-09
#> [6] 6.734287e-09
```

## Speed

Per call, in microseconds:

``` r

us_cpp <- time_per_call(function() logpost_cpp(init, Y), n = 10000)
us_grad <- time_per_call(function() grad_cpp(init, Y), n = 10000)
c(R = us_r, `C++` = us_cpp, `C++ gradient` = us_grad)
#>            R          C++ C++ gradient 
#>        139.0          9.0         11.7
```

The C++ log-posterior is about 15 times faster than the R one, and the
analytic gradient costs about the same as one log-posterior, where
finite differences need 40.

Now the full fit, with the VB correction, first without a gradient:

``` r

fit_cpp <- inla_vanilla(logpost_cpp, init, Y = Y)
#> ✔ Found the posterior mode in 79 iterations [104ms]
#> ✔ VB mean correction, max shift 0.43σ [7.0s]
#> ✔ Fitted 39 skew-normal marginals [236ms]
```

It took 500,872 calls, nearly all in the VB correction. In R, at 139
microseconds a call, these would take about 70 seconds. With the
gradient:

``` r

fit <- inla_vanilla(logpost_cpp, init, grad = grad_cpp, Y = Y)
#> ✔ Found the posterior mode in 91 iterations [10ms]
#> ✔ VB mean correction, max shift 0.43σ [266ms]
#> ✔ Fitted 39 skew-normal marginals [187ms]
```

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
checks the gradient against finite differences at `init`, and its
Jacobian at the mode for symmetry, and warns if either check fails. Here
neither does. The calls and times (in seconds) of the three fits:

``` r

fits <- list(`R, no VB` = fit_r, `C++` = fit_cpp, `C++ and gradient` = fit)
t(sapply(fits, function(f) c(f$counts, round(f$timing, 2))))
#>                  logpost  grad laplace   vb marginals total
#> R, no VB           14295     0    0.87 0.00      1.01  1.88
#> C++               500872     0    0.10 7.00      0.24  7.34
#> C++ and gradient    1537 15541    0.01 0.27      0.19  0.47
```

The three fits agree. The skew-normal marginals do not depend on the VB
correction, so even `fit_r` gives the same table. The largest difference
in a posterior mean or quantile, in posterior SDs:

``` r

sds <- fit$summary[, "SD"]
c(
  R = max(abs(fit_r$summary - fit$summary) / sds),
  `C++` = max(abs(fit_cpp$summary - fit$summary) / sds)
)
#>            R          C++ 
#> 5.689665e-05 4.961973e-05
```

## Results

``` r

summary(fit)
#> ── INLAvanilla summary ─────────────────────────────────────────────────────────
#> 
#> Posterior marginals
#>                     Mean      SD    2.5%     25%     50%     75%   97.5%    Mode
#> nu_x1             5.0538  0.0846  4.8880  4.9967  5.0538  5.1108  5.2195  5.0538
#> nu_x2              4.791   0.174   4.450   4.673   4.791   4.908   5.132   4.791
#> nu_x3              3.557   0.162   3.239   3.447   3.557   3.666   3.874   3.557
#> nu_y1              5.461   0.299   4.875   5.259   5.461   5.663   6.048   5.461
#> nu_y2              4.252   0.438   3.394   3.957   4.252   4.547   5.110   4.252
#> nu_y3              6.559   0.394   5.787   6.293   6.559   6.824   7.331   6.559
#> nu_y4              4.448   0.380   3.704   4.192   4.448   4.704   5.192   4.448
#> nu_y5              5.133   0.306   4.534   4.927   5.133   5.339   5.732   5.133
#> nu_y6              2.974   0.392   2.206   2.710   2.974   3.238   3.742   2.974
#> nu_y7              6.192   0.365   5.477   5.946   6.192   6.438   6.907   6.192
#> nu_y8              4.039   0.375   3.303   3.786   4.039   4.292   4.775   4.039
#> lam_x2             2.219   0.147   1.952   2.117   2.212   2.314   2.530   2.196
#> lam_x3             1.840   0.155   1.550   1.734   1.835   1.941   2.157   1.825
#> a                  1.211   0.145   0.944   1.111   1.204   1.304   1.512   1.191
#> b                  1.189   0.122   0.965   1.105   1.183   1.267   1.443   1.172
#> c                  1.276   0.128   1.045   1.187   1.269   1.358   1.548   1.253
#> beta_dem60_ind60   1.468   0.391   0.717   1.202   1.463   1.729   2.252   1.451
#> beta_dem65_ind60   0.587   0.242   0.120   0.423   0.585   0.749   1.068   0.580
#> beta_dem65_dem60  0.8679  0.0770  0.7216  0.8153  0.8662  0.9186  1.0238  0.8627
#> logv_ind60        -0.794   0.192  -1.166  -0.925  -0.796  -0.666  -0.412  -0.800
#> logv_dem60         1.359   0.227   0.909   1.208   1.361   1.513   1.797   1.366
#> logv_dem65        -1.598   0.861  -3.598  -2.109  -1.460  -0.959  -0.302  -1.097
#> logv_x1           -2.438   0.238  -2.926  -2.594  -2.431  -2.274  -1.989  -2.416
#> logv_x2           -2.165   0.566  -3.461  -2.504  -2.086  -1.752  -1.271  -1.896
#> logv_x3           -0.695   0.194  -1.072  -0.827  -0.697  -0.565  -0.310  -0.700
#> logv_y1            0.685   0.246   0.181   0.524   0.692   0.853   1.146   0.708
#> logv_y2            2.065   0.178   1.719   1.944   2.063   2.184   2.418   2.060
#> logv_y3            1.655   0.196   1.277   1.522   1.653   1.785   2.044   1.649
#> logv_y4            1.207   0.231   0.732   1.056   1.214   1.365   1.641   1.229
#> logv_y5            0.902   0.211   0.487   0.760   0.903   1.045   1.314   0.904
#> logv_y6            1.643   0.181   1.290   1.520   1.642   1.765   2.001   1.640
#> logv_y7            1.320   0.209   0.909   1.179   1.320   1.460   1.727   1.321
#> logv_y8            1.226   0.215   0.792   1.083   1.230   1.372   1.636   1.238
#> z_y1_y5            0.290   0.154  -0.021   0.188   0.293   0.395   0.582   0.300
#> z_y2_y4            0.283   0.123   0.041   0.200   0.283   0.366   0.524   0.283
#> z_y2_y6            0.361   0.108   0.153   0.288   0.360   0.433   0.577   0.358
#> z_y3_y7            0.184   0.136  -0.085   0.093   0.185   0.277   0.449   0.187
#> z_y4_y8            0.108   0.133  -0.158   0.019   0.109   0.198   0.365   0.112
#> z_y6_y8            0.328   0.123   0.086   0.245   0.328   0.411   0.569   0.329
#> 
#> Diagnostics
#>                   VB shift   SN - VB  Shape     NMAD   Bend
#> nu_x1             -0.00014  -0.00036  -0.14  4.1e-06  0.054
#> nu_x2             -0.00025   -0.0004  -0.14  4.4e-06  0.054
#> nu_x3             -0.00022  -0.00036  -0.13  3.9e-06  0.054
#> nu_y1              -0.0004  -0.00032  -0.15  6.9e-06   0.05
#> nu_y2             -0.00049  -0.00032  -0.14  4.7e-06  0.042
#> nu_y3              -0.0005  -0.00033  -0.15  5.9e-06  0.047
#> nu_y4             -0.00057  -0.00037  -0.15  6.5e-06  0.047
#> nu_y5             -0.00038  -0.00029  -0.15  6.4e-06  0.049
#> nu_y6             -0.00051  -0.00033  -0.15  6.5e-06  0.051
#> nu_y7             -0.00051  -0.00031  -0.15  6.5e-06  0.046
#> nu_y8             -0.00057  -0.00034  -0.15  6.6e-06  0.048
#> lam_x2                0.33     0.047    1.5   0.0067  0.043
#> lam_x3                0.18   -0.0032    1.1   0.0043  0.079
#> a                     0.16     0.012    1.4   0.0043  0.079
#> b                     0.12    0.0087    1.4   0.0048   0.11
#> c                     0.22     0.017    1.6   0.0069  0.081
#> beta_dem60_ind60    0.0019   -0.0049   0.82   0.0019  0.082
#> beta_dem65_ind60    -0.074    -0.011   0.69  0.00038  0.075
#> beta_dem65_dem60      0.12     0.014   0.98   0.0038  0.042
#> logv_ind60          -0.037   -0.0022   0.69    0.003  0.038
#> logv_dem60          0.0018   -0.0072  -0.69   0.0022  0.059
#> logv_dem65           -0.43       0.2   -4.2    0.045   0.44
#> logv_x1               0.26     0.044   -1.1    0.006   0.14
#> logv_x2               -0.3     0.031   -3.3    0.033   0.11
#> logv_x3               0.29     0.036   0.65    0.003  0.015
#> logv_y1               0.29     0.011   -1.1   0.0085  0.033
#> logv_y2               0.28   0.00062   0.63  0.00018  0.017
#> logv_y3               0.32     0.012   0.71  0.00072  0.018
#> logv_y4               0.18   0.00025   -1.1   0.0074  0.041
#> logv_y5               0.31    0.0093  -0.45    0.005  0.028
#> logv_y6                0.3    0.0032   0.52   0.0014  0.024
#> logv_y7                0.3    0.0084  -0.42   0.0057  0.018
#> logv_y8               0.21    0.0037  -0.92   0.0051  0.034
#> z_y1_y5             0.0079    -0.012  -0.97   0.0026   0.05
#> z_y2_y4             -0.074    0.0031  -0.25    0.007  0.094
#> z_y2_y6             -0.062    0.0097   0.74    0.012   0.12
#> z_y3_y7              0.051   -0.0042  -0.58   0.0055  0.061
#> z_y4_y8             -0.094   -0.0018  -0.74   0.0041   0.13
#> z_y6_y8             -0.085    0.0053  -0.35   0.0053   0.08
#> 
#> ✔ Mode found in 91 iterations (gradient 7.8e-05 SDs).
#> ✔ VB correction in 4 iterations: max shift 0.43 SDs.
#> ! The straight scan misses up to 0.44 nats at +-2 SDs for b, logv_dem65, logv_x1, logv_x2, z_y2_y6, and z_y4_y8: the marginal may be too narrow.
#> ✔ Skew-normal fits: max NMAD 0.045.
#> ! Skew-normal and VB means differ by up to 0.2 SDs, for logv_dem65.
#> 
#> Log marginal likelihood: -1686.847 (Laplace -1687.471, VB term 0.62)
#> 1537 logpost and 15541 grad calls in 470ms (mode 10ms, VB 267ms, marginals 188ms)
```

The intercepts are close to symmetric, and the loadings and regressions
are mildly right-skewed. The most skewed marginal is that of the
disturbance variance of `dem65`. Democracy in 1960 and industrialisation
explain most of the variation in democracy in 1965 (see $`R^2`$ below),
so this variance is small and poorly determined, and its log has a long
left tail. [`summary()`](https://rdrr.io/r/base/summary.html) flags it
and a few others: there, the straight scan through the mode misses some
log-density, and the marginal may come out too narrow. The skew-normal
and VB means of `logv_dem65` also differ by 0.2 posterior SDs, another
sign that its marginal is hard to approximate.

Here the warning is borne out. Checked against MCMC, the marginals of
`logv_dem65` and `logv_x2` are less than half as wide as the posterior,
with the 2.5% quantile of `logv_dem65` at −3.6 against −7.3. The data
cannot rule out a variance near zero, so on the log scale the posterior
has a long exponential left tail, held back only by the prior, and a
skew-normal cannot follow it. On the variance scale, where the tail near
zero is compressed, the error is much smaller: the posterior means of
the two variances are 0.2 posterior SDs too high. The other 37 marginals
agree with MCMC to within 7% in SD.

``` r

plot(fit, pars = c("beta_dem65_ind60", "beta_dem65_dem60", "logv_dem65"))
```

![Posterior densities of three parameters: the two regressions of 1965
democracy, both slightly right-skewed, and the log disturbance variance
of 1965 democracy, strongly left-skewed. Dashed Gaussian curves are
symmetric.](rcpp_files/figure-html/plot-1.png)

### Derived quantities

[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) gives
draws of any function of the parameters. Industrialisation affects
democracy in 1965 directly, and indirectly through democracy in 1960.
The indirect effect is a product of two regressions:

``` r

set.seed(2026)
effects <- sampling(fit, nsamp = 10000, fun = function(theta) {
  direct <- theta[["beta_dem65_ind60"]]
  indirect <- theta[["beta_dem60_ind60"]] * theta[["beta_dem65_dem60"]]
  c(direct = direct, indirect = indirect, total = direct + indirect)
})
round(t(apply(effects, 2, quantile, c(0.025, 0.5, 0.975))), 3)
#>           2.5%   50% 97.5%
#> direct   0.117 0.587 1.068
#> indirect 0.613 1.257 2.011
#> total    1.132 1.860 2.613
```

The share of the variance of `dem65` explained by industrialisation and
democracy in 1960, $`R^2`$, needs the covariance matrix of the latent
variables, $`A \Psi A^\top`$:

``` r

r2_dem65 <- sampling(fit, nsamp = 10000, fun = function(theta) {
  B <- matrix(0, 3, 3)
  B[beta_at] <- theta[17:19]
  A <- solve(diag(3) - B)
  M <- A %*% diag(exp(theta[20:22])) %*% t(A)
  1 - exp(theta[["logv_dem65"]]) / M[3, 3]
})
round(quantile(r2_dem65, c(0.025, 0.5, 0.975)), 3)
#>  2.5%   50% 97.5% 
#> 0.847 0.952 0.994
```

### Compared with maximum likelihood

The posterior modes are close to lavaan’s maximum likelihood estimates,
but not equal to them, because of the priors:

``` r

pars <- c(parnames[12:19], "logv_dem65")
round(
  cbind(
    ML = theta_ml[pars],
    `posterior mode` = fit$mode[pars],
    `posterior mean` = fit$summary[pars, "Mean"],
    `posterior SD` = fit$summary[pars, "SD"]
  ),
  3
)
#>                      ML posterior mode posterior mean posterior SD
#> lam_x2            2.180          2.168          2.219        0.147
#> lam_x3            1.818          1.814          1.840        0.155
#> a                 1.191          1.187          1.211        0.145
#> b                 1.175          1.174          1.189        0.122
#> c                 1.251          1.248          1.276        0.128
#> beta_dem60_ind60  1.471          1.470          1.468        0.391
#> beta_dem65_ind60  0.600          0.608          0.587        0.242
#> beta_dem65_dem60  0.865          0.858          0.868        0.077
#> logv_dem65       -1.805         -1.399         -1.598        0.861
```

The largest gap is again for the log disturbance variance of `dem65`. On
the log scale, its prior (a Gamma(1, 0.5) on the SD, times the Jacobian)
falls to zero for small variances, which lifts the posterior mode above
the maximum likelihood estimate. The long left tail then puts the
posterior mean below the mode.

## References

Bollen, K. A. (1989). *Structural equations with latent variables*.
Wiley. <https://doi.org/10.1002/9781118619179>
