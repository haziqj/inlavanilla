# Speeding up with Rcpp

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
calls `logpost` many times, and without a gradient the count grows with
the square of the number of parameters \\m\\: the Hessian by finite
differences takes about \\m^2\\ calls, the skew-normal scans about
\\2m^2\\, and the VB correction \\(m + 1)(m + 2)\\ per sweep, over a few
sweeps. For 39 parameters that is over 20,000 calls. Two things help: a
faster `logpost`, and a gradient.

This article does both for a structural equation model (SEM) with 39
parameters. We write the log-posterior in plain R first, then in C++
with [Rcpp](https://www.rcpp.org/) and
[RcppEigen](https://cran.r-project.org/package=RcppEigen), and then add
an analytic gradient in C++. The article [Gradients with
RTMB](http://haziqj.ml/inlavanilla/articles/rtmb.md) writes the same
model in R for RTMB, which computes the gradient by automatic
differentiation, and compares the routes for speed and accuracy.

``` r

library(inlavanilla)
data(PoliticalDemocracy, package = "lavaan")
```

## The model

The `PoliticalDemocracy` data in lavaan come from Bollen (1989). For 75
developing countries they hold three indicators of industrialisation in
1960 (`x1`–`x3`), and the same four indicators of political democracy in
1960 (`y1`–`y4`) and 1965 (`y5`–`y8`). In the lavaan syntax that
INLAvaan takes, the model is:

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
mean \\\nu\\ and covariance

\\ \Sigma = \Lambda (I - B)^{-1} \Psi (I - B)^{-\top} \Lambda^\top +
\Theta, \\

where \\\nu\\ holds the 11 intercepts, \\\Lambda\\ (\\11 \times 3\\) the
loadings, \\B\\ (\\3 \times 3\\) the regressions among the latent
variables, \\\Psi\\ the diagonal matrix of the variance of `ind60` and
the disturbance variances of `dem60` and `dem65`, and \\\Theta\\ (\\11
\times 11\\) the residual covariances. The log-likelihood is the sum of
75 multivariate normal log-densities.

## Parameters and priors

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
needs every parameter on an unconstrained scale. We use the
parameterisation and default priors of
[INLAvaan](https://inlavaan.haziqj.ml/), so the results can be compared
with `INLAvaan::asem(model, PoliticalDemocracy, meanstructure = TRUE)`:

| Parameters | Number | Unconstrained scale | Prior |
|:---|---:|:---|:---|
| Intercepts \\\nu\\ | 11 | as they are | \\N(0, 32^2)\\ |
| Loadings (`x2`, `x3`, `a`, `b`, `c`) | 5 | as they are | \\N(0, 10^2)\\ |
| Regressions | 3 | as they are | \\N(0, 10^2)\\ |
| Variances (3 latent, 11 residual) | 14 | \\\ell = \log \sigma^2\\ | \\\text{Gamma}(1, \text{rate } 0.5)\\ on \\\sigma\\ |
| Residual covariances | 6 | \\z = \operatorname{atanh} \rho\\ | \\\text{Beta}(1, 1)\\ on \\(\rho + 1)/2\\ |

A residual covariance enters through its correlation \\\rho\\, as
\\\theta\_{ij} = \rho \sqrt{\theta\_{ii} \theta\_{jj}}\\, and both its
parameter \\z\\ and its prior refer to \\\rho\\.

Two of the priors are on a different scale from the parameter, so
`logpost` needs their log-Jacobians:

- A variance’s prior is on its standard deviation \\\sigma = e^{\ell /
  2}\\, and \\d\sigma / d\ell = \sigma / 2\\. Its term in the
  log-posterior is \\\log p(\sigma) + \log(\sigma / 2)\\.
- A correlation’s prior is on \\\rho = \tanh z\\, and \\d\rho / dz = 1 -
  \rho^2\\. The Beta(1, 1) prior on \\(\rho + 1)/2\\ is the uniform
  density \\1/2\\ on \\(-1, 1)\\, so the term is \\\log(1/2) + \log(1 -
  \rho^2)\\.

The normal priors are on the parameters themselves and need none.

## The log-posterior in R

The data, ordered as the rows of \\\Lambda\\, and the names of the 39
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

The free entries of \\\Lambda\\, \\B\\ and \\\Theta\\, as matrix
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

### A check against INLAvaan

INLAvaan fits the same model, with the same parameterisation and priors:

``` r

fit_inlavaan <- INLAvaan::asem(
  model,
  PoliticalDemocracy,
  meanstructure = TRUE,
  verbose = FALSE
)
int <- INLAvaan::get_inlavaan_internal(fit_inlavaan)
```

Its parameters are ours, in another order and under lavaan’s names. Here
is its posterior mode in our order:

``` r

lv <- c("ind60", "dem60", "dem65")
inlavaan_names <- c(
  paste0(ov, "~1"),
  "ind60=~x2",
  "ind60=~x3",
  "a",
  "b",
  "c",
  "dem60~ind60",
  "dem65~ind60",
  "dem65~dem60",
  paste0(c(lv, ov), "~~", c(lv, ov)),
  paste0(cov_pairs[, 1], "~~", cov_pairs[, 2])
)
at <- match(inlavaan_names, names(int$coefficients))
theta_inlavaan <- setNames(int$theta_star_novbc[at], parnames)
```

There, our log-posterior should equal INLAvaan’s, which its optimiser
minimised as the negative:

``` r

c(ours = logpost_r(theta_inlavaan, Y), INLAvaan = -int$opt$objective)
#>      ours  INLAvaan 
#> -1655.078 -1655.078
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
#> [1] 68
```

With `vb = FALSE`,
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
skips the VB correction and runs only the mode search and the
skew-normal marginals:

``` r

fit_r <- inla_vanilla(logpost_r, init, Y = Y, vb = FALSE)
#> ✔ Found the posterior mode in 80 iterations [537ms]
#> ✔ Fitted 39 skew-normal marginals [674ms]
fit_r$counts
#> logpost    grad 
#>   14055       0
```

These two stages take 14,055 calls. The VB correction would add \\(m +
1)(m + 2) = 1{,}640\\ calls per sweep, over a few sweeps: without a
gradient it takes the expected gradient from values of `logpost` alone.
We come back to the cost of a full fit below.

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

The gradient of the log-likelihood has a closed form (Neudecker &
Satorra, 1991), and costs about as much as the log-likelihood itself.
Write \\D\\ for the \\n \times 11\\ matrix of data minus \\\nu\\ and
\\\Omega = \Sigma^{-1}\\. The derivative of the log-likelihood is
\\\Omega D^\top \mathbf{1}\\ with respect to \\\nu\\, and

\\ G = \tfrac{1}{2} \left( \Omega D^\top D \\ \Omega - n \Omega \right)
\\

with respect to \\\Sigma\\. Through \\\Sigma = \Lambda M \Lambda^\top +
\Theta\\, with \\A = (I - B)^{-1}\\ and \\M = A \Psi A^\top\\, the
derivatives with respect to the model matrices are

\\ \frac{\partial \ell}{\partial \Lambda} = 2 G \Lambda M, \qquad
\frac{\partial \ell}{\partial B} = 2 A^\top \Lambda^\top G \Lambda M,
\qquad \frac{\partial \ell}{\partial \Psi} = A^\top \Lambda^\top G
\Lambda A, \\

and \\G\_{ii}\\ for a residual variance \\\theta\_{ii}\\ and \\2
G\_{ij}\\ for a residual covariance. A loading held equal across years
sums its two entries. On the unconstrained scale, a log variance
multiplies its derivative by the variance. A residual covariance
\\\theta\_{ij} = \rho \sqrt{\theta\_{ii} \theta\_{jj}}\\ moves with
\\z\\, through \\d\theta\_{ij} / dz = (1 - \rho^2) \sqrt{\theta\_{ii}
\theta\_{jj}}\\, and with both log variances, through \\d\theta\_{ij} /
d\ell_i = \theta\_{ij} / 2\\.

The log-prior adds \\-\vartheta / s^2\\ for each normal prior, \\1/2 -
\sigma / 4\\ for each log variance and \\-2\rho\\ for each correlation.
We append the gradient to the same file:

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
starting values, at INLAvaan’s posterior mode and at some random points:

``` r

set.seed(1)
points <- c(
  list(init, theta_inlavaan),
  replicate(4, init + rnorm(39, sd = 0.3), simplify = FALSE)
)
vapply(points, function(th) logpost_cpp(th, Y) - logpost_r(th, Y), 0)
#> [1] -4.547474e-13  0.000000e+00  0.000000e+00 -2.273737e-13  0.000000e+00
#> [6] -4.547474e-13
max(abs(loglik_cpp(theta_inlavaan, Y) - loglik_r(theta_inlavaan, Y)))
#> [1] 7.105427e-15
```

The gradient should match central differences of the R log-posterior, up
to their own error, about \\10^{-8}\\ relative:

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
#> [1] 1.461419e-08 2.310907e-08 9.214642e-09 7.379461e-09 5.716762e-09
#> [6] 6.734287e-09
```

## Speed

Per call, in microseconds:

``` r

us_cpp <- time_per_call(function() logpost_cpp(init, Y), n = 10000)
us_grad <- time_per_call(function() grad_cpp(init, Y), n = 10000)
c(R = us_r, `C++` = us_cpp, `C++ gradient` = us_grad)
#>            R          C++ C++ gradient 
#>         68.0          7.3          9.0
```

The C++ log-posterior is about 9 times faster than the R one, and the
analytic gradient costs about the same as one log-posterior, where
finite differences need 40.

Now the full fit, with the VB correction, first without a gradient:

``` r

fit_cpp <- inla_vanilla(logpost_cpp, init, Y = Y)
#> ✔ Found the posterior mode in 82 iterations [76ms]
#> ✔ VB mean correction, max shift 0.39σ [137ms]
#> ✔ Fitted 39 skew-normal marginals [209ms]
```

It took 22,437 calls. In R, at 68 microseconds a call, these would take
about 1.5 seconds. With the gradient:

``` r

fit <- inla_vanilla(logpost_cpp, init, grad = grad_cpp, Y = Y)
#> ✔ Found the posterior mode in 89 iterations [8ms]
#> ✔ VB mean correction, max shift 0.43σ [122ms]
#> ✔ Fitted 39 skew-normal marginals [133ms]
```

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
checks the gradient against finite differences at `init`, and its
Jacobian at the mode for symmetry, and warns if either check fails. Here
neither does. The calls and times (in seconds) of the three fits:

``` r

fits <- list(`R, no VB` = fit_r, `C++` = fit_cpp, `C++ and gradient` = fit)
t(sapply(fits, function(f) c(f$counts, round(f$timing, 2))))
#>                  logpost grad laplace   vb marginals total
#> R, no VB           14055    0    0.54 0.00      0.67  1.22
#> C++                22437    0    0.08 0.14      0.21  0.43
#> C++ and gradient    1535 9931    0.01 0.12      0.13  0.27
```

The VB correction takes the same nodes either way: the gradient at
\\(m + 1)(m + 2) + 1 = 1{,}641\\ points per sweep, or `logpost` at all
but the centre. The gradient buys accuracy there: its VB mean is exact
when the log-posterior is a polynomial of degree six, against degree
four from `logpost` alone.

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
#> 4.488435e-05 2.977629e-05
```

## Results

``` r

summary(fit)
#> ── INLAvanilla summary ───────────────────────────────── Skew-normal + VB ──
#> 
#> Log marginal likelihood: -1686.850 · Converged in 89 iterations (269ms)
#> Laplace -1687.471 + VB term 0.62
#> 1535 logpost and 9931 grad calls: mode 8ms, VB 122ms, marginals 133ms
#> 
#> ── Posterior marginals ──
#> 
#>                         Mean      SD     2.5%      50%    97.5%     Mode
#>  1 nu_x1              5.0538  0.0846   4.8880   5.0538   5.2195   5.0538
#>  2 nu_x2              4.791   0.174    4.450    4.791    5.132    4.791 
#>  3 nu_x3              3.557   0.162    3.239    3.557    3.874    3.557 
#>  4 nu_y1              5.461   0.299    4.875    5.461    6.048    5.461 
#>  5 nu_y2              4.252   0.438    3.394    4.252    5.110    4.252 
#>  6 nu_y3              6.559   0.394    5.787    6.559    7.331    6.559 
#>  7 nu_y4              4.448   0.380    3.704    4.448    5.192    4.448 
#>  8 nu_y5              5.133   0.306    4.534    5.133    5.732    5.133 
#>  9 nu_y6              2.974   0.392    2.206    2.974    3.742    2.974 
#> 10 nu_y7              6.192   0.365    5.477    6.192    6.907    6.192 
#> 11 nu_y8              4.039   0.375    3.303    4.039    4.775    4.039 
#> 12 lam_x2             2.219   0.147    1.952    2.212    2.530    2.196 
#> 13 lam_x3             1.840   0.155    1.550    1.835    2.157    1.825 
#> 14 a                  1.211   0.145    0.944    1.204    1.512    1.191 
#> 15 b                  1.189   0.122    0.965    1.183    1.443    1.172 
#> 16 c                  1.276   0.128    1.045    1.269    1.548    1.253 
#> 17 beta_dem60_ind60   1.468   0.391    0.717    1.463    2.252    1.451 
#> 18 beta_dem65_ind60   0.587   0.242    0.120    0.585    1.068    0.580 
#> 19 beta_dem65_dem60   0.8679  0.0770   0.7216   0.8662   1.0238   0.8627
#> 20 logv_ind60        -0.794   0.192   -1.166   -0.796   -0.412   -0.800 
#> 21 logv_dem60         1.359   0.227    0.909    1.361    1.797    1.366 
#> 22 logv_dem65        -1.598   0.861   -3.598   -1.460   -0.302   -1.097 
#> 23 logv_x1           -2.438   0.238   -2.926   -2.431   -1.989   -2.416 
#> 24 logv_x2           -2.165   0.566   -3.461   -2.086   -1.271   -1.896 
#> 25 logv_x3           -0.695   0.194   -1.072   -0.697   -0.310   -0.700 
#> 26 logv_y1            0.685   0.246    0.181    0.692    1.146    0.708 
#> 27 logv_y2            2.065   0.178    1.719    2.063    2.418    2.060 
#> 28 logv_y3            1.655   0.196    1.277    1.653    2.044    1.649 
#> 29 logv_y4            1.207   0.231    0.733    1.214    1.641    1.229 
#> 30 logv_y5            0.902   0.211    0.487    0.903    1.314    0.904 
#> 31 logv_y6            1.643   0.181    1.290    1.642    2.001    1.640 
#> 32 logv_y7            1.320   0.209    0.909    1.320    1.727    1.321 
#> 33 logv_y8            1.226   0.215    0.792    1.230    1.636    1.238 
#> 34 z_y1_y5            0.290   0.154   -0.021    0.293    0.582    0.300 
#> 35 z_y2_y4            0.283   0.123    0.041    0.283    0.524    0.283 
#> 36 z_y2_y6            0.361   0.108    0.153    0.360    0.577    0.358 
#> 37 z_y3_y7            0.184   0.136   -0.085    0.185    0.449    0.187 
#> 38 z_y4_y8            0.108   0.133   -0.158    0.109    0.365    0.112 
#> 39 z_y6_y8            0.328   0.123    0.086    0.328    0.569    0.329 
#> 
#> ── Diagnostics ──
#> 
#>                       VB shift       NMAD  Unscanned  Ridge bend  Mean gap
#>  1 nu_x1             -0.000080  0.0000042   0.000063       0.054  -0.00036
#>  2 nu_x2             -0.00024   0.0000044   0.000063       0.054  -0.00037
#>  3 nu_x3             -0.00018   0.0000039   0.000063       0.054  -0.00034
#>  4 nu_y1             -0.00041   0.0000068   0.000063       0.050  -0.00032
#>  5 nu_y2             -0.00050   0.0000047   0.000063       0.042  -0.00031
#>  6 nu_y3             -0.00050   0.0000059   0.000063       0.047  -0.00033
#>  7 nu_y4             -0.00056   0.0000065   0.000063       0.047  -0.00035
#>  8 nu_y5             -0.00052   0.0000063   0.000063       0.049  -0.00019
#>  9 nu_y6             -0.00061   0.0000065   0.000063       0.051  -0.00023
#> 10 nu_y7             -0.00062   0.0000064   0.000063       0.046  -0.00020
#> 11 nu_y8             -0.00069   0.0000066   0.000063       0.048  -0.00022
#> 12 lam_x2             0.33      0.0067      0.0014         0.043   0.049  
#> 13 lam_x3             0.18      0.0043      0.00032        0.079  -0.0027 
#> 14 a                  0.16      0.0043      0.00048        0.079   0.012  
#> 15 b                  0.12      0.0048      0.00042        0.11    0.0083 
#> 16 c                  0.22      0.0069      0.00090        0.081   0.016  
#> 17 beta_dem60_ind60   0.0018    0.0019      0.000093       0.082  -0.0047 
#> 18 beta_dem65_ind60  -0.074     0.00038     0.000071       0.075  -0.011  
#> 19 beta_dem65_dem60   0.12      0.0038      0.00024        0.042   0.016  
#> 20 logv_ind60        -0.037     0.0030      0.000069       0.038  -0.0028 
#> 21 logv_dem60         0.0018    0.0022      0.000078       0.059  -0.0073 
#> 22 logv_dem65        -0.43      0.045       0.0019         0.44    0.19   
#> 23 logv_x1            0.26      0.0060      0.000072       0.14    0.045  
#> 24 logv_x2           -0.30      0.033       0.0031         0.11    0.031  
#> 25 logv_x3            0.29      0.0030      0.00023        0.015   0.036  
#> 26 logv_y1            0.29      0.0085      0.000081       0.033   0.011  
#> 27 logv_y2            0.28      0.00018     0.00018        0.017   0.0010 
#> 28 logv_y3            0.32      0.00072     0.00025        0.018   0.012  
#> 29 logv_y4            0.17      0.0074      0.00011        0.041   0.0012 
#> 30 logv_y5            0.31      0.0050      0.00013        0.028   0.0094 
#> 31 logv_y6            0.30      0.0014      0.00016        0.024   0.0042 
#> 32 logv_y7            0.30      0.0057      0.00013        0.018   0.0093 
#> 33 logv_y8            0.21      0.0051      0.000084       0.034   0.0054 
#> 34 z_y1_y5            0.0089    0.0026      0.00013        0.050  -0.013  
#> 35 z_y2_y4           -0.075     0.0070      0.000052       0.094   0.0038 
#> 36 z_y2_y6           -0.062     0.012       0.000051       0.12    0.0095 
#> 37 z_y3_y7            0.052     0.0055      0.000057       0.061  -0.0056 
#> 38 z_y4_y8           -0.094     0.0041      0.00010        0.13   -0.0022 
#> 39 z_y6_y8           -0.087     0.0053      0.000058       0.080   0.0068 
#> 
#> ✔ Mode found in 89 iterations (Newton step 9.1e-05σ).
#> ✔ Hessian condition number 65, in SD units.
#> ✔ VB correction settled in 4 iterations: max shift 0.43σ.
#> ✔ Skew-normal fits: max NMAD 0.045, unscanned mass 0.0031.
#> ! Ridge bend > 0.25 for `logv_dem65` (0.44): the marginal SD may be >10%
#>   too small.
#> ! Mean gap > 0.1σ for `logv_dem65` (0.19): the posterior is far from
#>   Gaussian there, and either mean may be off.
```

The intercepts are close to symmetric, and the loadings and regressions
are mildly right-skewed. The most skewed marginal is that of the
disturbance variance of `dem65`. Democracy in 1960 and industrialisation
explain most of the variation in democracy in 1965 (see \\R^2\\ below),
so this variance is small and poorly determined, and its log has a long
left tail. [`summary()`](https://rdrr.io/r/base/summary.html) flags it
twice: the straight scan through the mode misses enough log-density that
its marginal may come out too narrow, and its skew-normal and VB means
differ by 0.2 posterior SDs, another sign that the marginal is hard to
approximate.

Here the flags are borne out. Checked against MCMC, the marginals of
`logv_dem65` and `logv_x2` are less than half as wide as the posterior,
with the 2.5% quantile of `logv_dem65` at −3.6 against −7.3. The data
cannot rule out a variance near zero, so on the log scale the posterior
has a long exponential left tail, held back only by the prior, and a
skew-normal cannot follow it. On the variance scale, where the tail near
zero is compressed, the error is much smaller: the posterior means of
the two variances are 0.2 posterior SDs too high. No check flags
`logv_x2`: its ridge bends less, and none of the checks is aimed at a
tail like this. The other 37 marginals agree with MCMC to within 7% in
SD.

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
democracy in 1960, \\R^2\\, needs the covariance matrix of the latent
variables, \\A \Psi A^\top\\:

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

### Compared with INLAvaan

INLAvaan runs the same approximation on the same log-posterior, so the
two fits should agree. The modes and the Laplace covariances, in
posterior SDs and relative to the largest covariance:

``` r

sds <- sqrt(diag(fit$cov))
c(
  mode = max(abs(theta_inlavaan - fit$mode) / sds),
  cov = max(abs(int$Sigma_theta[at, at] - fit$cov)) / max(abs(fit$cov))
)
#>         mode          cov 
#> 9.982413e-05 1.559542e-04
```

Both agree to a few parts in \\10^{4}\\, the precision of the two
optimisers and finite-difference Hessians.

INLAvaan fits a skew-normal to each parameter as well. Its means and
SDs, on the same scale, next to ours for the regressions and the
variance of `dem65`:

``` r

sn <- int$approx_data[inlavaan_names, ]
delta <- sn[, "alpha"] / sqrt(1 + sn[, "alpha"]^2)
sn_mean <- sn[, "xi"] + sn[, "omega"] * delta * sqrt(2 / pi)
sn_sd <- sn[, "omega"] * sqrt(1 - 2 * delta^2 / pi)
names(sn_mean) <- names(sn_sd) <- parnames

pars <- c(parnames[17:19], "logv_dem65")
round(
  cbind(
    `mean, INLAvaan` = sn_mean[pars],
    `mean, ours` = fit$summary[pars, "Mean"],
    `SD, INLAvaan` = sn_sd[pars],
    `SD, ours` = fit$summary[pars, "SD"]
  ),
  4
)
#>                  mean, INLAvaan mean, ours SD, INLAvaan SD, ours
#> beta_dem60_ind60         1.4685     1.4685       0.3915   0.3915
#> beta_dem65_ind60         0.5872     0.5872       0.2417   0.2417
#> beta_dem65_dem60         0.8679     0.8679       0.0770   0.0770
#> logv_dem65              -1.5981    -1.5980       0.8610   0.8610
max(abs(sn_mean - fit$summary[, "Mean"]) / fit$summary[, "SD"])
#> [1] 6.996214e-05
```

Over all 39 parameters, the means differ by at most about \\10^{-4}\\
posterior SDs. The two packages share the skew-normal marginals, and
with them the limitation for `logv_dem65` and `logv_x2` described above.

The log marginal likelihoods share their Laplace part. The VB terms
differ a little, because INLAvaan takes the expectation in the VB
correction with 64 Sobol points and inlavanilla with a deterministic
cubature rule:

``` r

round(
  rbind(
    INLAvaan = c(
      Laplace = int$mloglik - int$vb$kld_global,
      VB = int$vb$kld_global,
      total = int$mloglik
    ),
    ours = c(
      Laplace = fit$logml_laplace,
      VB = fit$logml - fit$logml_laplace,
      total = fit$logml
    )
  ),
  3
)
#>            Laplace    VB     total
#> INLAvaan -1687.471 0.687 -1686.784
#> ours     -1687.471 0.621 -1686.850
```

The Laplace parts agree to the last digit. The VB terms differ by 0.07,
from the two integration rules: 64 Sobol points estimate the expectation
with Monte Carlo error, which the deterministic cubature rule avoids.

## References

Bollen, K. A. (1989). *Structural equations with latent variables*.
Wiley. <https://doi.org/10.1002/9781118619179>

Neudecker, H., & Satorra, A. (1991). Linear structural relations:
Gradient and Hessian of the fitting function. *Statistics & Probability
Letters*, *11*(1), 57–61. <https://doi.org/10.1016/0167-7152(91)90178-T>
