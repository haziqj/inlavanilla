# Articles

### Get started

- [Get started with
  inlavanilla](https://haziqj.ml/inlavanilla/articles/inlavanilla.md):

### Examples

- [Multilevel
  models](https://haziqj.ml/inlavanilla/articles/multilevel.md):

  A linear mixed model with the random effects integrated out and drawn
  back afterwards, checked against brms.

- [Cognitive diagnostic
  models](https://haziqj.ml/inlavanilla/articles/cdm.md):

  A DINA model for binary items, with the latent classes summed out,
  checked against Stan.

- [Correlation
  matrices](https://haziqj.ml/inlavanilla/articles/correlation.md):

  A multivariate normal with an LKJ prior on its correlation matrix,
  written in terms of canonical partial correlations, checked against
  Stan.

- [Binary factor
  analysis](https://haziqj.ml/inlavanilla/articles/binary-cfa.md):

  A one-factor model for binary items, with the factor integrated out by
  Gauss-Hermite quadrature inside the log-posterior.

- [Diagnostic tests without a gold
  standard](https://haziqj.ml/inlavanilla/articles/diagnostic-tests.md):

  Three latent class models for five schistosomiasis tests in Ghana,
  with the infection status summed out and a random effect integrated by
  quadrature.

### Articles

- [Natural-scale
  summaries](https://haziqj.ml/inlavanilla/articles/natural-scale.md):

- [Diagnostics](https://haziqj.ml/inlavanilla/articles/diagnostics.md):

- [Model comparison with LOO and
  WAIC](https://haziqj.ml/inlavanilla/articles/model-comparison.md):

- [Latent variables](https://haziqj.ml/inlavanilla/articles/latents.md):

  Random effects, class memberships and other latent variables are
  integrated out of the log-posterior. Four ways to draw them back after
  the fit.

- [Speeding up with
  Rcpp](https://haziqj.ml/inlavanilla/articles/rcpp.md):

  A structural equation model with 39 parameters, its log-posterior
  written first in R and then in C++ with Rcpp and RcppEigen, with an
  analytic gradient.

- [Gradients and Hessians with
  RTMB](https://haziqj.ml/inlavanilla/articles/rtmb.md):

  The structural equation model of the Rcpp article written for RTMB,
  which differentiates R code automatically to give the gradient and the
  Hessian, and compared with R and C++ for speed and accuracy.
