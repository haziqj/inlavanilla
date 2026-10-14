# Using TMB

[TMB](https://github.com/kaskr/adcomp) (Kristensen et al., 2016)
compiles a log-density written as a C++ template, and computes its
gradient by automatic differentiation. The gradient is exact, and
nothing has to be derived by hand.
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
takes TMB’s objective and gradient as they are, apart from the sign: TMB
minimises.

This article fits the industrialisation and political democracy model of
Bollen (1989), a structural equation model with 39 parameters. The
article [Speeding up with
Rcpp](http://haziqj.ml/inlavanilla/articles/rcpp.md) describes the
model, its unconstrained parameterisation, the priors and their
Jacobians, and writes the same log-posterior in R and in C++. Here we
write it as a TMB template.

``` r

library(inlavanilla)
data(PoliticalDemocracy, package = "lavaan")
```

## Data and parameters

The 11 indicators, in the order `x1`, `x2`, `x3`, `y1`, …, `y8`, and the
six pairs with correlated residuals:

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
```

The 39 parameters, all on an unconstrained scale, in six blocks:

| Block | Number | Parameters | Prior |
|:---|---:|:---|:---|
| `nu` | 11 | intercepts | \\N(0, 32^2)\\ |
| `lam` | 5 | loadings of `x2`, `x3`, and `a`, `b`, `c` | \\N(0, 10^2)\\ |
| `beta` | 3 | regressions | \\N(0, 10^2)\\ |
| `logv_lv` | 3 | log variance of `ind60`, log disturbance variances of `dem60` and `dem65` | \\\text{Gamma}(1, \text{rate } 0.5)\\ on the SD |
| `logv_ov` | 11 | log residual variances | \\\text{Gamma}(1, \text{rate } 0.5)\\ on the SD |
| `z` | 6 | residual correlations \\\rho\\ on the atanh scale | \\\text{Beta}(1, 1)\\ on \\(\rho + 1)/2\\ |

We give them readable names, and the same starting values as in the Rcpp
article:

``` r

parnames <- c(
  paste0("nu_", ov),
  c("lam_x2", "lam_x3", "a", "b", "c"),
  c("beta_dem60_ind60", "beta_dem65_ind60", "beta_dem65_dem60"),
  paste0("logv_", c("ind60", "dem60", "dem65")),
  paste0("logv_", ov),
  paste0("z_", cov_pairs[, 1], "_", cov_pairs[, 2])
)
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

## The template

A TMB template declares its data and parameters, and returns the value
to minimise, here the negative log-posterior. `density::MVNORM_t` is
TMB’s multivariate normal; it returns the negative log-density,
normalising constant included. The 75 pointwise log-likelihoods are kept
in a vector, and `REPORT()` makes them available from R. The log-prior
adds the log-Jacobians of the standard deviations and correlations, as
in the Rcpp article.

``` r

tmb_file <- file.path(tempdir(), "sem_tmb.cpp")
```

``` cpp
#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()() {
  DATA_MATRIX(Y);            // 75 x 11: x1, x2, x3, y1, ..., y8

  PARAMETER_VECTOR(nu);      // 11 intercepts
  PARAMETER_VECTOR(lam);     // loadings: x2, x3, then a, b, c
  PARAMETER_VECTOR(beta);    // dem60 ~ ind60, dem65 ~ ind60, dem65 ~ dem60
  PARAMETER_VECTOR(logv_lv); // log variances of ind60, dem60, dem65
  PARAMETER_VECTOR(logv_ov); // log residual variances, x1 to y8
  PARAMETER_VECTOR(z);       // residual correlations, atanh scale

  int n = Y.rows(), p = Y.cols();

  // Loadings, with the first of each latent variable fixed at 1, and a, b, c
  // shared between 1960 and 1965
  matrix<Type> Lambda(p, 3);
  Lambda.setZero();
  Lambda(0, 0) = 1;
  Lambda(1, 0) = lam(0);
  Lambda(2, 0) = lam(1);
  Lambda(3, 1) = 1;
  Lambda(4, 1) = lam(2);
  Lambda(5, 1) = lam(3);
  Lambda(6, 1) = lam(4);
  Lambda(7, 2) = 1;
  Lambda(8, 2) = lam(2);
  Lambda(9, 2) = lam(3);
  Lambda(10, 2) = lam(4);

  matrix<Type> B(3, 3);
  B.setZero();
  B(1, 0) = beta(0);
  B(2, 0) = beta(1);
  B(2, 1) = beta(2);
  matrix<Type> I(3, 3);
  I.setIdentity();
  matrix<Type> A = atomic::matinv(matrix<Type>(I - B));

  vector<Type> psi = exp(logv_lv);
  vector<Type> v = exp(logv_ov);
  vector<Type> rho = tanh(z);
  matrix<Type> Psi = psi.matrix().asDiagonal();
  matrix<Type> Theta = v.matrix().asDiagonal();
  // Residual covariances: y1-y5, y2-y4, y2-y6, y3-y7, y4-y8, y6-y8
  int ci[6] = {3, 4, 4, 5, 6, 8};
  int cj[6] = {7, 6, 8, 9, 10, 10};
  for (int k = 0; k < 6; k++) {
    Theta(ci[k], cj[k]) = rho(k) * sqrt(v(ci[k]) * v(cj[k]));
    Theta(cj[k], ci[k]) = Theta(ci[k], cj[k]);
  }

  matrix<Type> LA = Lambda * A;
  matrix<Type> Sigma = LA * Psi * LA.transpose() + Theta;

  // Pointwise log-likelihood
  density::MVNORM_t<Type> mvn(Sigma);
  vector<Type> loglik(n);
  for (int i = 0; i < n; i++) {
    vector<Type> r = vector<Type>(Y.row(i)) - nu;
    loglik(i) = -mvn(r);
  }
  REPORT(loglik);

  // Log-prior on the unconstrained scale, with the Jacobians
  Type logprior = sum(dnorm(nu, Type(0), Type(32), true)) +
    sum(dnorm(lam, Type(0), Type(10), true)) +
    sum(dnorm(beta, Type(0), Type(10), true));
  vector<Type> logv(14);
  logv << logv_lv, logv_ov;
  for (int k = 0; k < 14; k++) {
    // Gamma(1, rate 0.5) on the SD (TMB's dgamma() takes the scale), and
    // d SD / d logv = SD / 2
    Type sd = exp(logv(k) / Type(2));
    logprior += dgamma(sd, Type(1), Type(1 / 0.5), true) + log(sd / Type(2));
  }
  for (int k = 0; k < 6; k++) {
    // Beta(1, 1) on (rho + 1) / 2, and d rho / d z = 1 - rho^2
    Type u = (rho(k) + Type(1)) / Type(2);
    logprior += dbeta(u, Type(1), Type(1), true) - log(Type(2)) +
      log(Type(1) - rho(k) * rho(k));
  }

  return -(sum(loglik) + logprior);
}
```

## Compiling

[`TMB::compile()`](https://rdrr.io/pkg/TMB/man/compile.html) builds a
shared library from the template, which takes a minute or so.
[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html) then
records the computation and returns an object whose `fn` and `gr`
evaluate the objective and its gradient. The parameters go in as a list
of the six vectors, and `obj$par` holds them in one vector, in the same
order as `init`.

``` r

invisible(TMB::compile(tmb_file))
dyn.load(TMB::dynlib(sub("\\.cpp$", "", tmb_file)))

blocks <- rep(
  c("nu", "lam", "beta", "logv_lv", "logv_ov", "z"),
  c(11, 5, 3, 3, 11, 6)
)
obj <- TMB::MakeADFun(
  data = list(Y = Y),
  parameters = split(unname(init), factor(blocks, unique(blocks))),
  DLL = "sem_tmb",
  silent = TRUE
)
all.equal(unname(obj$par), unname(init))
#> [1] TRUE
```

### A check against INLAvaan

INLAvaan fits the same model, with the same parameterisation and priors.
As in the Rcpp article, we put its posterior mode in our order:

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
fit_inlavaan <- INLAvaan::asem(
  model,
  PoliticalDemocracy,
  meanstructure = TRUE,
  verbose = FALSE
)
int <- INLAvaan::get_inlavaan_internal(fit_inlavaan)

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

There, the template’s objective should equal the one INLAvaan minimised:

``` r

c(TMB = obj$fn(theta_inlavaan), INLAvaan = int$opt$objective)
#>      TMB INLAvaan 
#> 1655.078 1655.078
```

## Fitting

TMB’s `fn` returns the negative log-posterior and `gr` its gradient as a
one-row matrix, so we negate both and flatten the gradient. With `init`
as the starting values, the parameters keep our names; `obj$par` would
also do, and its repeated names (`nu`, `nu`, …) would be numbered.

``` r

fit <- inla_vanilla(
  function(theta) -obj$fn(theta),
  init,
  grad = function(theta) -as.vector(obj$gr(theta))
)
#> ✔ Found the posterior mode in 91 iterations [54ms]
#> ✔ VB mean correction, max shift 0.43σ [2.1s]
#> ✔ Fitted 39 skew-normal marginals [683ms]
```

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
checks the gradient against finite differences at `init`, and would warn
if it disagreed. The calls, and the time per call in microseconds:

``` r

time_per_call <- function(f, n) {
  elapsed <- system.time(
    for (i in seq_len(n)) {
      f()
    }
  )[["elapsed"]]
  1e6 * elapsed / n
}
fit$counts
#> logpost    grad 
#>    1537   15541
c(
  fn = time_per_call(function() obj$fn(init), n = 2000),
  gr = time_per_call(function() obj$gr(init), n = 1000)
)
#>    fn    gr 
#>  59.5 146.0
```

A gradient costs a few evaluations of the objective, against 40 for
finite differences. Without it, the VB correction alone would take over
120,000 evaluations per sweep. A hand-written gradient, as in the Rcpp
article, can be faster still, but it has to be derived and checked;
TMB’s needs neither.

## Results

``` r

fit
#> ── INLAvanilla ───────────────────────────────────────────── Skew-normal + VB ──
#> 
#>                     Mean      SD    2.5%     50%   97.5%    Mode
#> nu_x1             5.0538  0.0846  4.8880  5.0538  5.2195  5.0538
#> nu_x2              4.791   0.174   4.450   4.791   5.132   4.791
#> nu_x3              3.557   0.162   3.239   3.557   3.874   3.557
#> nu_y1              5.461   0.299   4.875   5.461   6.048   5.461
#> nu_y2              4.252   0.438   3.394   4.252   5.110   4.252
#> nu_y3              6.559   0.394   5.787   6.559   7.331   6.559
#> nu_y4              4.448   0.380   3.704   4.448   5.192   4.448
#> nu_y5              5.133   0.306   4.534   5.133   5.732   5.133
#> nu_y6              2.974   0.392   2.206   2.974   3.742   2.974
#> nu_y7              6.192   0.365   5.477   6.192   6.907   6.192
#> nu_y8              4.039   0.375   3.303   4.039   4.775   4.039
#> lam_x2             2.219   0.147   1.952   2.212   2.530   2.196
#> lam_x3             1.840   0.155   1.550   1.835   2.157   1.825
#> a                  1.211   0.145   0.944   1.204   1.512   1.191
#> b                  1.189   0.122   0.965   1.183   1.443   1.172
#> c                  1.276   0.128   1.045   1.269   1.548   1.253
#> beta_dem60_ind60   1.468   0.391   0.717   1.463   2.252   1.451
#> beta_dem65_ind60   0.587   0.242   0.120   0.585   1.068   0.580
#> beta_dem65_dem60  0.8679  0.0770  0.7216  0.8662  1.0238  0.8627
#> logv_ind60        -0.794   0.192  -1.166  -0.796  -0.412  -0.800
#> logv_dem60         1.359   0.227   0.909   1.361   1.797   1.366
#> logv_dem65        -1.598   0.861  -3.598  -1.460  -0.302  -1.097
#> logv_x1           -2.438   0.238  -2.926  -2.431  -1.989  -2.416
#> logv_x2           -2.165   0.566  -3.461  -2.086  -1.271  -1.896
#> logv_x3           -0.695   0.194  -1.072  -0.697  -0.310  -0.700
#> logv_y1            0.685   0.246   0.181   0.692   1.146   0.708
#> logv_y2            2.065   0.178   1.719   2.063   2.418   2.060
#> logv_y3            1.655   0.196   1.277   1.653   2.044   1.649
#> logv_y4            1.207   0.231   0.732   1.214   1.641   1.229
#> logv_y5            0.902   0.211   0.487   0.903   1.314   0.904
#> logv_y6            1.643   0.181   1.290   1.642   2.001   1.640
#> logv_y7            1.320   0.209   0.909   1.320   1.727   1.321
#> logv_y8            1.226   0.215   0.792   1.230   1.636   1.238
#> z_y1_y5            0.290   0.154  -0.021   0.293   0.582   0.300
#> z_y2_y4            0.283   0.123   0.041   0.283   0.524   0.283
#> z_y2_y6            0.361   0.108   0.153   0.360   0.577   0.358
#> z_y3_y7            0.184   0.136  -0.085   0.185   0.449   0.187
#> z_y4_y8            0.108   0.133  -0.158   0.109   0.365   0.112
#> z_y6_y8            0.328   0.123   0.086   0.328   0.569   0.329
#> 
#> Log marginal likelihood: -1686.85
#> ! 2 checks flagged by `summary()`.
```

This agrees with the Rcpp article’s fit, which uses the same
log-posterior, to within about \\10^{-4}\\ posterior SDs. The two checks
flagged by [`summary()`](https://rdrr.io/r/base/summary.html) are the
scan-path bend, for several parameters, and the gap between the
skew-normal and VB means of `logv_dem65`; the [Rcpp
article](http://haziqj.ml/inlavanilla/articles/rcpp.md) discusses both.

``` r

plot(fit, pars = paste0("z_", cov_pairs[, 1], "_", cov_pairs[, 2]))
```

![Posterior densities of the six residual correlations on the atanh
scale, each close to a Gaussian, with the Laplace Gaussians
dashed.](tmb_files/figure-html/plot-1.png)

### On the natural scale

The parameters are on the unconstrained scale, and
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) with
`fun` maps its draws back. Here are the residual correlations \\\rho =
\tanh z\\ and the disturbance variance of `dem65`, next to the posterior
medians from INLAvaan’s own copula draws, which are on the natural
scale.

``` r

set.seed(2026)
natural <- sampling(fit, nsamp = 10000, fun = function(theta) {
  rho <- tanh(theta[paste0("z_", cov_pairs[, 1], "_", cov_pairs[, 2])])
  names(rho) <- paste0("rho_", cov_pairs[, 1], "_", cov_pairs[, 2])
  c(rho, psi_dem65 = exp(theta[["logv_dem65"]]))
})
draws <- INLAvaan::sampling(fit_inlavaan, nsamp = 10000)
v <- function(x) draws[, paste0(x, "~~", x)]
rho <- draws[, paste0(cov_pairs[, 1], "~~", cov_pairs[, 2])] /
  sqrt(v(cov_pairs[, 1]) * v(cov_pairs[, 2]))
round(
  cbind(
    t(apply(natural, 2, quantile, c(0.025, 0.5, 0.975))),
    INLAvaan = c(apply(rho, 2, median), median(v("dem65")))
  ),
  3
)
#>             2.5%   50% 97.5% INLAvaan
#> rho_y1_y5 -0.024 0.288 0.528    0.286
#> rho_y2_y4  0.039 0.276 0.475    0.274
#> rho_y2_y6  0.153 0.345 0.520    0.345
#> rho_y3_y7 -0.084 0.184 0.421    0.184
#> rho_y4_y8 -0.158 0.109 0.347    0.111
#> rho_y6_y8  0.080 0.315 0.509    0.317
#> psi_dem65  0.028 0.234 0.739    0.230
```

The medians agree with INLAvaan’s, up to the Monte Carlo error of the
two sets of draws. The disturbance variance of `dem65` is small and
poorly determined, with a long right tail on this scale.

### Leave-one-out cross-validation

`obj$report()` evaluates the template at any parameter value without
recording derivatives, and returns the 75 pointwise log-likelihoods: the
`loglik` that [`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md)
needs. It leaves out one country at a time, from this one fit:

``` r

loo_sem <- loo(fit, function(theta) obj$report(theta)$loglik)
loo_sem
#> ── Leave-one-out cross-validation ────────────────────────────────── 75 units ──
#> 
#>           Estimate    SE
#> elpd_loo   -1588.8  20.5
#> p_loo         38.6   3.0
#> looic       3177.7  40.9
#> 
#> ✔ Second-order expansion for every unit. Largest leverage k = 0.25, for unit
#>   35.
```

INLAvaan’s [`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md) uses
the same deterministic expansion (Alhyari et al., 2026), with its own
per-case log-likelihoods:

``` r

loo_inlavaan <- INLAvaan::loo(fit_inlavaan)
round(
  rbind(
    ours = loo_sem$estimates[, "Estimate"],
    INLAvaan = loo_inlavaan$estimates[, "Estimate"]
  ),
  2
)
#>          elpd_loo p_loo   looic
#> ours     -1588.83 38.60 3177.66
#> INLAvaan -1588.96 38.77 3177.92
```

The two agree to well within the standard error. They differ a little
because each expands about its own VB-corrected mean, and the two VB
corrections take their expectations with different rules.

With 39 parameters this takes \\39^2 + 39 + 1\\, or 1,561, evaluations
of the template, each giving all 75 terms. See [Model comparison with
LOO and WAIC](http://haziqj.ml/inlavanilla/articles/model-comparison.md)
for what the numbers mean and how to compare models with them.

## References

Alhyari, M., Jamil, H., Montcho, H., & Rue, H. (2026). *Deterministic
leave-one-cluster-out cross-validation for multilevel Bayesian models*
(arXiv:2609.00670 \[stat.ME\]). arXiv.
<https://doi.org/10.48550/arXiv.2609.00670>

Bollen, K. A. (1989). *Structural equations with latent variables*.
Wiley. <https://doi.org/10.1002/9781118619179>

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M.
(2016). TMB: Automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, *70*(5), 1–21.
<https://doi.org/10.18637/jss.v070.i05>
