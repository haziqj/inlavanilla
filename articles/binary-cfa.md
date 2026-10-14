# Binary factor analysis

A factor model for binary items has a latent factor score for each
person, and the likelihood integrates it out. Here a Gauss-Hermite
quadrature inside the log-posterior does the integral, so
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
sees only the item parameters.

``` r

library(inlavanilla)
```

## The model

Person \\i\\ answers item \\j\\ with \\y\_{ij} = 1\\ when a latent
response \\y^\*\_{ij} = \lambda_j \eta_i + \epsilon\_{ij}\\ exceeds the
threshold \\\tau_j\\, with factor \\\eta_i \sim N(0, 1)\\. We use the
*delta* parameterisation, as lavaan (Rosseel, 2012) calls it, which
fixes \\\operatorname{Var}(y^\*\_{ij}) = 1\\. Then \\\lambda_j\\ is the
correlation of \\y^\*\_{ij}\\ with the factor, \\\tau_j =
-\Phi^{-1}\\P(y\_{ij} = 1)\\\\, and

\\ P(y\_{ij} = 1 \mid \eta_i) = \Phi\left( \frac{\lambda_j \eta_i -
\tau_j}{\sqrt{1 - \lambda_j^2}} \right). \\

Each person’s likelihood is an integral over their factor, which an
\\n\\-point Gauss-Hermite rule with nodes \\z_q\\ and weights \\w_q\\
turns into a sum:

\\ p(y_i) = \int \prod_j P(y\_{ij} \mid \eta) \\ \phi(\eta) \\ d\eta
\approx \sum\_{q=1}^{n} w_q \prod_j P(y\_{ij} \mid z_q). \\

The priors are \\\lambda_j \sim U(-1, 1)\\ and \\\tau_j \sim N(0, 1)\\,
the latter a uniform prior on each item’s proportion of ones. The
parameters are \\\operatorname{atanh} \lambda_j\\ and \\\tau_j\\, and
the uniform prior on \\\lambda_j\\ becomes \\(1 - \lambda_j^2)/2\\ on
the atanh scale.

## The data

Simulated answers of 500 people to six items:

``` r

n <- 500
lambda <- c(0.5, 0.6, 0.7, 0.8, 0.85, 0.55)
tau <- c(-0.9, -0.4, 0, 0.2, 0.45, 1)
set.seed(1)
eta <- rnorm(n)
y_star <- outer(eta, lambda) +
  rnorm(n * 6) * rep(sqrt(1 - lambda^2), each = n)
y <- 1 * (y_star > rep(tau, each = n))
```

The likelihood depends on the data only through the counts of the
distinct response patterns, of which there are at most \\2^6 = 64\\:

``` r

tab <- table(apply(y, 1, paste, collapse = ""))
patterns <- 1 * (do.call(rbind, strsplit(names(tab), "")) == "1")
count <- as.vector(tab)
nrow(patterns)
#> [1] 48
```

## The quadrature

The Golub and Welsch (1969) algorithm gives the nodes and weights for
\\N(0, 1)\\ as the eigenvalues and squared first components of the
eigenvectors of a tridiagonal matrix. With 21 nodes, the rule integrates
polynomials up to degree 41 exactly, among them the variance:

``` r

gh_normal <- function(n) {
  jacobi <- matrix(0, n, n)
  off <- cbind(1:(n - 1), 2:n)
  jacobi[off] <- jacobi[off[, 2:1]] <- sqrt(1:(n - 1))
  e <- eigen(jacobi, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}
nodes <- gh_normal(21)
c(sum(nodes$w), sum(nodes$w * nodes$z^2))
#> [1] 1 1
```

## The log-posterior

With patterns in rows and nodes in columns, the log-probabilities of all
the patterns at all the nodes take two matrix products. A log-sum-exp
over the nodes then gives each pattern’s log-likelihood:

``` r

logpost <- function(theta, patterns, count, nodes) {
  p <- ncol(patterns)
  lambda <- tanh(theta[1:p])
  tau <- theta[p + 1:p]
  a <- (outer(lambda, nodes$z) - tau) / sqrt(1 - lambda^2) # item x node
  ell <- patterns %*% pnorm(a, log.p = TRUE) +
    (1 - patterns) %*% pnorm(-a, log.p = TRUE)
  ell <- t(t(ell) + log(nodes$w))
  ell_max <- apply(ell, 1, max)
  loglik <- sum(count * (ell_max + log(rowSums(exp(ell - ell_max)))))
  loglik + sum(log1p(-lambda^2) - log(2)) + sum(dnorm(tau, log = TRUE))
}
```

## The fit

Flipping the signs of all the loadings and of the factor leaves the
likelihood unchanged, so the posterior has a mirror image. Starting at
positive loadings, the fit stays on that side. No gradient is needed:

``` r

init <- c(
  setNames(rep(atanh(0.5), 6), paste0("atanh_lambda", 1:6)),
  setNames(rep(0, 6), paste0("tau", 1:6))
)
fit <- inla_vanilla(
  logpost,
  init,
  patterns = patterns,
  count = count,
  nodes = nodes
)
#> ✔ Found the posterior mode in 35 iterations [137ms]
#> ✔ VB mean correction, max shift 0.17σ [111ms]
#> ✔ Fitted 12 skew-normal marginals [134ms]
fit
#> ── INLAvanilla ───────────────────────────────────────── Skew-normal + VB ──
#> 
#> Log marginal likelihood: -1683.00 · Converged in 35 iterations (397ms)
#> 
#>                      Mean      SD     2.5%      50%    97.5%     Mode
#>  1 atanh_lambda1   0.4016  0.0832   0.2402   0.4009   0.5665   0.3996
#>  2 atanh_lambda2   0.5854  0.0819   0.4279   0.5843   0.7491   0.5821
#>  3 atanh_lambda3   0.940   0.104    0.747    0.936    1.155    0.927 
#>  4 atanh_lambda4   0.8750  0.0952   0.6955   0.8724   1.0695   0.8670
#>  5 atanh_lambda5   1.377   0.188    1.058    1.358    1.794    1.317 
#>  6 atanh_lambda6   0.6151  0.0916   0.4396   0.6135   0.7992   0.6104
#>  7 tau1           -0.8395  0.0637  -0.9652  -0.8392  -0.7153  -0.8386
#>  8 tau2           -0.4339  0.0579  -0.5477  -0.4338  -0.3207  -0.4335
#>  9 tau3            0.0050  0.0559  -0.1046   0.0050   0.1146   0.0050
#> 10 tau4            0.1462  0.0562   0.0362   0.1461   0.2563   0.1461
#> 11 tau5            0.4523  0.0579   0.3390   0.4522   0.5662   0.4519
#> 12 tau6            0.8627  0.0641   0.7379   0.8624   0.9893   0.8618
```

Quantiles carry over through \\\tanh\\, so the loadings’ 95% intervals
and medians are those of \\\operatorname{atanh} \lambda_j\\,
transformed. Next to the true values:

``` r

q <- fit$summary[, c("2.5%", "50%", "97.5%")]
q[1:6, ] <- tanh(q[1:6, ])
rownames(q) <- c(paste0("lambda", 1:6), paste0("tau", 1:6))
round(cbind(true = c(lambda, tau), q), 3)
#>          true   2.5%    50%  97.5%
#> lambda1  0.50  0.236  0.381  0.513
#> lambda2  0.60  0.404  0.526  0.635
#> lambda3  0.70  0.633  0.733  0.819
#> lambda4  0.80  0.602  0.703  0.789
#> lambda5  0.85  0.785  0.876  0.946
#> lambda6  0.55  0.413  0.547  0.664
#> tau1    -0.90 -0.965 -0.839 -0.715
#> tau2    -0.40 -0.548 -0.434 -0.321
#> tau3     0.00 -0.105  0.005  0.115
#> tau4     0.20  0.036  0.146  0.256
#> tau5     0.45  0.339  0.452  0.566
#> tau6     1.00  0.738  0.862  0.989
```

Ten of the twelve intervals cover the true values. The two misses come
from the simulated data, not from the approximation. Item 6, for one,
got 97 ones where 79 were expected, and the thresholds follow the data:
each posterior median is close to \\-\Phi^{-1}\\ of the item’s
proportion of ones.

``` r

round(rbind(data = -qnorm(colMeans(y)), posterior = q[7:12, "50%"]), 3)
#>             tau1   tau2  tau3  tau4  tau5  tau6
#> data      -0.842 -0.434 0.005 0.146 0.451 0.863
#> posterior -0.839 -0.434 0.005 0.146 0.452 0.862
```

## Why the delta parameterisation

The *theta* parameterisation fixes the residual variance at one instead,
with loadings \\\lambda_j / \sqrt{1 - \lambda_j^2}\\ and thresholds
\\\tau_j / \sqrt{1 - \lambda_j^2}\\. The model is the same, but its
thresholds scale with the loadings, which bends the posterior along
curved ridges. `lab/binary-cfa-delta.R`, in the package source, fits
both to 500 simulated data sets, with normal priors on the theta scale,
so the same posterior in both. Both covered the true values in about 95%
of the data sets. But [`summary()`](https://rdrr.io/r/base/summary.html)
flagged 37% to 93% of the theta fits, mostly for the ridge bend, against
under 7% of the delta fits, and the delta fits were closer to Stan.

For two correlated factors, the nodes become a product grid, mapped
through the Cholesky factor of the factor correlation matrix;
`lab/binary-cfa-model.R` has the code.

## References

Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss quadrature
rules. *Mathematics of Computation*, *23*(106), 221–230.
<https://doi.org/10.1090/S0025-5718-69-99647-1>

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software*, *48*(2), 1–36.
<https://doi.org/10.18637/jss.v048.i02>
