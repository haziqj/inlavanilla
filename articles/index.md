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

  A linear mixed model with the random effects integrated out, checked
  against brms.

- [Diagnostic classification
  models](http://haziqj.ml/inlavanilla/articles/dcm.md):

  A DINA model for binary items, with the latent classes summed out,
  checked against Stan.

- [Correlation
  matrices](http://haziqj.ml/inlavanilla/articles/correlation.md):

  A multivariate normal with an LKJ prior on its correlation matrix,
  written in terms of canonical partial correlations, checked against
  Stan.

### Compiled code

Faster log-posteriors for larger models. Both articles fit the same
structural equation model, with 39 parameters.

- [Speeding up with
  Rcpp](http://haziqj.ml/inlavanilla/articles/rcpp.md):

  A structural equation model with 39 parameters, its log-posterior
  written first in R and then in C++ with Rcpp and RcppEigen, with an
  analytic gradient.

- [Using TMB](http://haziqj.ml/inlavanilla/articles/tmb.md):

  The structural equation model of the Rcpp article as a TMB template,
  with the gradient from automatic differentiation.
