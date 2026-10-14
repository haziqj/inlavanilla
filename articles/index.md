# Articles

### Get started

- [Get started with
  inlavanilla](http://haziqj.ml/inlavanilla/articles/inlavanilla.md):

### Checking a fit

- [Diagnostics](http://haziqj.ml/inlavanilla/articles/diagnostics.md):
- [Model comparison with LOO and
  WAIC](http://haziqj.ml/inlavanilla/articles/model-comparison.md):

### Examples

- [Multilevel
  models](http://haziqj.ml/inlavanilla/articles/multilevel.md):

  A linear mixed model with the random effects integrated out and drawn
  back afterwards, checked against brms.

- [Diagnostic classification
  models](http://haziqj.ml/inlavanilla/articles/dcm.md):

  A DINA model for binary items, with the latent classes summed out,
  checked against Stan.

- [Correlation
  matrices](http://haziqj.ml/inlavanilla/articles/correlation.md):

  A multivariate normal with an LKJ prior on its correlation matrix,
  written in terms of canonical partial correlations, checked against
  Stan.

- [Binary factor
  analysis](http://haziqj.ml/inlavanilla/articles/binary-cfa.md):

  A one-factor model for binary items, with the factor integrated out by
  Gauss-Hermite quadrature inside the log-posterior.

### Speed and gradients

Faster log-posteriors and exact gradients for larger models. Both
articles fit the same structural equation model, with 39 parameters.

- [Speeding up with
  Rcpp](http://haziqj.ml/inlavanilla/articles/rcpp.md):

  A structural equation model with 39 parameters, its log-posterior
  written first in R and then in C++ with Rcpp and RcppEigen, with an
  analytic gradient.

- [Gradients with RTMB](http://haziqj.ml/inlavanilla/articles/rtmb.md):

  The structural equation model of the Rcpp article written for RTMB,
  which differentiates R code automatically, and compared with R and C++
  for speed and accuracy.
