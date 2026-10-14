# Summarise an INLAvanilla fit

Posterior summaries from the skew-normal marginals, on the unconstrained
scale of `logpost`, together with diagnostics of each stage of the
approximation.

## Usage

``` r
# S3 method for class 'inla_vanilla'
summary(object, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- object:

  An object of class `inla_vanilla`.

- probs:

  Probabilities for the reported quantiles.

- ...:

  Unused.

## Value

An object of class `summary.inla_vanilla`: a list with the summary table
`table`, the log marginal likelihood `logml` (and its Laplace part
`logml_laplace`), the per-parameter diagnostics `diagnostics`, and
details of the optimiser, the VB correction, the number of calls and the
timing.

## Details

Most of the checks are those of INLAvaan's
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md),
and
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
returns them under the same names. Jamil and Rue (2026) use the VB shift
and NMAD as fit diagnostics (Sections 3.1 and 4.4), and the gap between
the skew-normal and VB means as an internal check (Section 5.1). The
mode shift, the Hessian condition number and the scan-end mass come from
INLAvaan's
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md),
not from the paper. `tails` and `bend` are inlavanilla's own.

A posterior SD here is the Laplace SD at the mode,
\\\sqrt{\Sigma\_{jj}}\\, as in INLAvaan's
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md);
it is not the SD column of the table, which is the skew-normal's. The
diagnostics, with the level at which each is flagged where it is, are:

- **Mode.** `mode_shift_sigma`, the Newton step at the mode in posterior
  SDs: how far one step of Newton's method would still move each
  parameter. It is flagged above 0.1, whatever the return code of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), whose
  message is then shown: without a gradient
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) can report failure
  at an accurate mode, and success short of one.

- **Hessian.** `hess_cond`, the condition number of the Hessian at the
  mode with each parameter in posterior-SD units (the eigenvalue ratio
  of the posterior correlation matrix). It is flagged above \\10^8\\:
  some combination of the parameters is then weakly identified. Its line
  is shown when there is more than one parameter.

- **VB correction.** `vb_shift_sigma`, the shift of each mean from the
  mode in posterior SDs (`vb_shift` is the same shift on the scale of
  `logpost`). The VB correction is the shift of the Laplace Gaussian
  that minimises its Kullback-Leibler divergence from the posterior, and
  a large shift is a sign of a posterior far from that Gaussian (Jamil &
  Rue, 2026, Section 3.1). It is flagged above 1 SD. `kld`, the
  Kullback-Leibler divergence between each Laplace Gaussian marginal and
  its VB-shifted copy, is `vb_shift_sigma^2 / 2`, as in INLAvaan; it is
  not flagged. When the fixed-point iteration hands over to the fallback
  optimiser, the VB line says so, but, as in INLAvaan, that is not
  flagged: it is a route of the solver, not a sign of a poor
  approximation.

- **Skew-normal fits.** `nmad`, the largest gap between the fitted and
  the scanned density, relative to its peak (Jamil & Rue, 2026, Section
  4.4), flagged above 0.1. `alpha` (the Shape column) is the shape of
  the fitted skew-normal marginal, on the unconstrained scale: zero is a
  Gaussian marginal, and the sign gives the direction of the skew. It is
  not flagged.

- **Volume correction.** `gamma1`, the slope added to the scanned
  log-heights, per posterior SD along the scan, for the changing volume
  of the other parameters (Jamil & Rue, 2026, Corollary 3.3). It is in
  `summary(fit)$diagnostics` and in the verbose output of the fit, not
  in the printed table, and it is not flagged.

- **Scan-end mass.** `scan_end_mass`, the mass each fitted skew-normal
  puts outside the scanned window (4 posterior SDs either side of the
  mode), where the marginal is extrapolated, as INLAvaan's
  [`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
  reports it. A Gaussian marginal gives \\2\Phi(-4) = 6.3 \times
  10^{-5}\\; the larger it is, the more of the marginal is extrapolated.
  It is flagged above 0.05.

- **Tails.** `tails` is `FALSE` (and flagged) when the density at an end
  of the scan is still over 5% of its peak: the Laplace SD badly
  understates a heavy tail there, and the skew-normal extrapolates it.

- **Scan path.** Each parameter is scanned along the conditional mode
  path of the Laplace Gaussian (Jamil & Rue, 2026, Lemma 3.1), a
  straight line that is a first-order approximation to a possibly curved
  ridge (Section 3.2.1). `bend` is the log-density, in nats, that the
  straight scan misses at \\\pm 2\\ posterior SDs: the gain from one
  Newton step towards the conditional mode of the other parameters. It
  is flagged above 0.1. A large value means the other parameters'
  conditional mode bends away from the line, as a poorly determined
  noise scale does for a regression coefficient at small n, and the
  marginal comes out too narrow, by roughly 0.35 to 0.5 times `bend` in
  relative SD (see the Accuracy section of
  [`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)).

- **Skew-normal and VB means.** `sn_vs_vb`, the gap between the
  skew-normal mean and the VB-corrected mean in posterior SDs. Both
  means target the posterior mean, so they should agree closely (to
  \\O(n^{-2})\\; Jamil & Rue, 2026, Proposition 3.4), and a gap signals
  higher-order effects or numerical error (Section 5.1). A gap above 0.1
  SDs is flagged: the posterior is far from Gaussian along that
  parameter, and either approximation may be off.

Every check flagged here but the scan path and the SN-VB gap also raises
one warning at the end of the fit, of class
`inlavanilla_diagnostics_warning`, at the same level. The levels of the
mode shift, NMAD, scan-end mass, VB shift and Hessian condition number
are those of INLAvaan's fit-time warning; the heavy-tails check is
inlavanilla's own, as are the scan-path and SN-VB flags. Every level is
a package default, chosen so that a healthy fit stays silent, and not a
calibrated cut-off: a tripped check is a prompt to look again, and
silence is not a certificate of accuracy. In the simulations of Jamil
and Rue (2026, Section 4.4) the VB shift and NMAD track the loss of
accuracy, and the paper reports a calibration of the expected error
against them in its Supplementary Material.

The status lines list every flagged parameter, worst first, with its
value (the heavy-tails line, a yes/no check, lists names only, in
parameter order); the warning lists at most three, as INLAvaan's does.
With `silent = TRUE` the warning is only recorded in `fit$warnings`.
Printing the fit ends with a line for its warnings and one counting the
checks flagged here, each shown only when there are any.
[`vignette("diagnostics", package = "inlavanilla")`](http://haziqj.ml/inlavanilla/articles/diagnostics.md)
explains every check with examples.

With `vb = FALSE` there is no VB mean: the VB columns of
`summary(fit)$diagnostics` are `NA`, as INLAvaan's
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
leaves them, and are left out of the table. With
`marginals = "gaussian"` there are no scans, so the scan columns are
`NA` and only the VB shift is shown.

## References

Jamil, H., & Rue, H. (2026). *Approximate Bayesian inference for
structural equation models using integrated nested Laplace
approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
<https://doi.org/10.48550/arXiv.2603.25690>

## See also

[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
for the diagnostics as numbers, under INLAvaan's names, and
[`plot.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/plot.inla_vanilla.md)
with `type = "sn_fit"` for the skew-normal fits.
