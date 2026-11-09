#' INLA-style posterior approximation for any log-posterior
#'
#' @description
#' Approximates the posterior of a parameter vector from its log density
#' alone, using the approximation pipeline of INLAvaan (Jamil & Rue, 2026a,
#' 2026b) with no latent field and nothing model-specific. The pipeline has
#' four stages:
#'
#' 1. **Joint Gaussian (Laplace).** The posterior mode \eqn{\vartheta^*} and
#'    the inverse Hessian \eqn{\Omega} there give a Gaussian approximation
#'    \eqn{N(\vartheta^*, \Omega)} and the Laplace estimate of the log marginal
#'    likelihood.
#' 2. **VB mean correction.** The Gaussian is recentred at
#'    \eqn{\tilde\vartheta = \vartheta^* + \hat\delta}, the mean that maximises
#'    its expected log-posterior, with the covariance unchanged.
#' 3. **Skew-normal marginals.** For each parameter, the log-posterior is
#'    scanned along the conditional mode path of the Gaussian, corrected for the
#'    changing volume of the other parameters, and fitted by a skew-normal
#'    \eqn{SN(\hat\xi_j, \hat\omega_j, \hat\alpha_j)}.
#' 4. **Gaussian copula draws.** [sampling()] couples the skew-normal marginals
#'    through the correlations of \eqn{\Omega}.
#'
#' Stages 2 and 3 are optional and independent of each other. `vb = "off"`
#' skips the VB correction, which leaves the marginals unchanged (they are
#' scanned around the mode) but gives the Laplace log marginal likelihood.
#' `marginals = "gaussian"` skips the scans and uses the Gaussian's own
#' marginals \eqn{N(\tilde\vartheta_j, \Omega_{jj})}, centred at the VB mean
#' (or at the mode when `vb = "off"`); [sampling()] then draws from the joint
#' Gaussian. So `vb = "off", marginals = "gaussian"` is the Laplace
#' approximation alone, and `marginals = "gaussian"` alone is the Laplace
#' approximation with the VB correction.
#'
#' @section Writing `logpost`:
#' The approximation works on an unconstrained scale: every parameter must be
#' free to take any real value. Reparameterise anything constrained, for
#' example a standard deviation \eqn{\sigma > 0} as \eqn{\log\sigma}, a
#' probability as its logit, or a correlation as its Fisher \eqn{z}. Means and
#' regression coefficients need no change.
#'
#' Return `-Inf` wherever the density is zero or cannot be evaluated, such as
#' outside the support or where a matrix is not positive definite (`NA` and
#' `NaN` count the same). The scans and the VB nodes reach several posterior
#' SDs from the mode, and an error raised there stops the fit.
#'
#' `logpost` is then the log density **of the parameters you pass**. When a
#' prior is stated on the natural scale, include the log-Jacobian of the
#' transformation: with \eqn{\vartheta = \log\sigma} and a prior
#' \eqn{p(\sigma)}, the prior term is \eqn{\log p(e^\vartheta) + \vartheta}.
#' Leaving it out targets a different posterior, not a worse approximation of
#' the same one, and nothing downstream can detect it.
#'
#' The choice of parameterisation also sets the quality of the approximation,
#' since each parameter gets its own skew-normal on the scale you choose.
#' Scales on which the posterior is close to Gaussian (log standard
#' deviations, Fisher-\eqn{z} correlations, logits) work best. Natural-scale or
#' derived quantities are available through `sampling(fit, fun = ...)`.
#'
#' The log marginal likelihood is an estimate of \eqn{\log p(y)} only when
#' `logpost` returns the full \eqn{\log p(y, \vartheta)}, including every
#' normalising constant of the likelihood and the prior. Otherwise it is off by
#' the missing constants.
#'
#' @section Accuracy:
#' The skew-normal marginals capture skewness that the Laplace approximation
#' misses, and they are accurate to the extent that a skew-normal fits each
#' marginal (see the diagnostics in [summary.inla_vanilla()]). Like INLAvaan,
#' each marginal is scanned along a straight line through the mode, and the
#' changing volume of the other parameters enters only to first order. When a
#' parameter's spread depends strongly on another, poorly determined
#' parameter, its marginal comes out too narrow, by a factor that vanishes as
#' O(1/n). Examples are a regression coefficient or a mean when the noise SD
#' is uncertain. For a normal mean with unknown variance the SD is about 11%
#' too small at n = 12, 3% at n = 50, and under 1% at n = 200. The
#' `ridge_bend` diagnostic of [summary.inla_vanilla()] flags this narrowing
#' where it predicts about 10% or more. By default (`width = "trace"`) this
#' is the scan of INLAvaan; `width = "tk"`, or `tk_fallback = TRUE` for the
#' flagged parameters only, removes the narrowing at a higher cost (see the
#' section on the width term).
#'
#' A skew-normal also cannot follow a tail heavier than a Gaussian's. In a
#' structural equation model with 39 parameters and n = 75, the data could not
#' keep two residual variances away from zero, so on the log scale their
#' posteriors had a long, exponential left tail, set by the prior. Checked
#' against MCMC, those two skew-normal marginals were less than half as wide
#' as the posterior, with the 2.5% quantile about 2 posterior SDs too high;
#' the other 37 were within 7% in SD. On the variance scale, where the tail
#' near zero is compressed, the two were much closer: means 0.2 posterior SDs
#' too high and SDs 10 to 15% too small. INLAvaan's marginals of this model
#' are the same as inlavanilla's: Jamil and Rue (2026a, Section 4.4) find
#' the same miss in the lower tail of a variance near zero, the boundary
#' pile-up that a unimodal skew-normal cannot capture. No check is aimed at
#' this. The `ridge_bend` diagnostic of [summary.inla_vanilla()] flagged one
#' of the two (0.44), but the other's ridge bend (0.11) was no larger than
#' that of marginals within 7% in SD, and the `tails` check, which looks only
#' 4 posterior SDs either side of the mode (the Laplace SDs), missed both.
#'
#' @section The width term:
#' Each skew-normal scan follows the conditional mode path (CMP) of the
#' Laplace Gaussian, \eqn{\vartheta(t) = \vartheta^* + t v_j} with
#' \eqn{v_j = \Omega_{\cdot j} / \sqrt{\Omega_{jj}}}{v_j = Omega[, j] /
#' sqrt(Omega[j, j])}, so that \eqn{t} counts posterior SDs. The log marginal
#' is then about the height of the log-posterior along the path plus the
#' width term \eqn{\gamma_j(t) = -\frac12 \log|H_{-j}(t)|}{gamma_j(t) =
#' -1/2 log|H_-j(t)|}, where \eqn{H_{-j}}{H_-j} is the negative Hessian of
#' the log-posterior with row and column j deleted: the log-volume of the
#' other parameters. `width` sets how it is taken, in three levels that each
#' remove one error:
#' * `"trace"`, the default and INLAvaan's, linearises it,
#'   \eqn{\gamma_j(t) \approx \gamma_j' t}{gamma_j(t) ~ gamma_j' t}, with
#'   the slope from whitened differences of the gradient (Jamil & Rue, 2026a,
#'   Corollary 3.3).
#' * `"cmp"` takes \eqn{\gamma_j(t)}{gamma_j(t)} exactly along the same
#'   path. Its difference from `"trace"` is the error of the linearisation.
#' * `"tk"` moves the other parameters to their conditional mode at each
#'   point of the scan and takes the height and the width term there: the
#'   Laplace approximation of the marginal (Tierney & Kadane, 1986), INLA's
#'   "Laplace strategy". Its difference from `"cmp"` is the error of the
#'   straight path, which the `ridge_bend` diagnostic of
#'   [summary.inla_vanilla()] predicts.
#'
#' The skew-normal is fitted to the corrected heights in the same way for
#' all three, on the same grid. `tk_fallback = TRUE` keeps the cheap scan
#' for every parameter but those whose ridge bend is flagged, and rescans
#' those with `"tk"`. In `fit$diagnostics$marginals`, `width` gives the
#' method each parameter used and `gamma1` the slope of its width term at
#' the mode (for `"cmp"` and `"tk"`, a central difference across the mode).
#' The ridge bend of a parameter scanned with `"tk"` is still reported, but
#' not flagged.
#'
#' @section Gradients and Hessians:
#' `grad` is optional. Without it, derivatives come from finite differences
#' of `logpost`, and the number of calls grows with the square of the number
#' of parameters m. The mode and Hessian take about m^2 + 6m calls, and the
#' marginals about m(`n_grid` + 2m + 6), with 2(m - 1) more for each parameter
#' whose ridge bend is above 0.1. The VB correction takes (m + 1)(m + 2)
#' calls, or that many per sweep with `vb = "full"`: it needs no derivatives,
#' because Stein's identity gives the expected gradient from values of
#' `logpost` alone. With `vb = "full"`, the VB mean is then exact when the
#' log-posterior is a polynomial of degree four in whitened coordinates, and up
#' to degree six with `grad`. On test models it was about 0.01 posterior SDs from exact VB
#' (0.003 with `grad`), and up to 0.1 on a strongly skewed posterior (0.03
#' with `grad`), against a VB error of 0.01 to 0.2 SDs from the posterior mean
#' itself. On the most skewed posteriors the VB correction may then not
#' settle, and the fit says it has not converged (see
#' [summary.inla_vanilla()]); `grad` is the remedy.
#'
#' With `grad`, every stage is more accurate. The fit is also faster when a
#' gradient call costs about as much as a `logpost` call, as a hand-written
#' gradient usually does, since it then makes far fewer calls in all. A
#' gradient by automatic differentiation (RTMB, TMB) costs about two, and
#' beyond a handful of parameters the fit is still several times faster than
#' without one (four times at m = 39). The default VB correction,
#' `vb = "onestep"`, takes 2m + 1 gradient calls, about as many as one numerical
#' Hessian, and `vb = "full"` takes at most (m + 1)(m + 2) + 1 per sweep, usually
#' about four sweeps. Beyond a handful of parameters the marginals are the
#' slowest stage, along with the VB correction under `vb = "full"`. `vb = "off"`
#' skips the VB correction (and its term in the log marginal likelihood), and
#' the skew-normal marginals do not depend on it;
#' `marginals = "gaussian"` saves the calls of the scans.
#'
#' `hess` is optional too, and needs `grad`, as in [stats::nlminb()]. With
#' it, [stats::nlminb()] takes Newton steps to the mode, and the Hessian at
#' the mode is the supplied one. The VB correction does not use it, nor do the
#' default marginals, whose volume correction takes 4m + 1 `grad` calls for all
#' the parameters together, where a Hessian each would cost far more.
#' `width = "cmp"` and `width = "tk"` do use it. The fit is the one `grad`
#' alone gives, to the error of its finite differences (about 1e-8, relative),
#' except that the mode is found more precisely. A dense Hessian costs about m
#' gradients, so this takes a little longer: a structural equation model with
#' 39 parameters and an automatic-differentiation Hessian (RTMB) took 0.18 s
#' with it and 0.17 s without.
#'
#' A supplied `grad` is checked against finite differences at `init`, and
#' its Jacobian at the mode is checked for symmetry, with a warning if either
#' check fails. A supplied `hess` is checked at the mode against central
#' differences of `grad`, along two directions, with a warning if it does
#' not match. Any R function works for `logpost`, `grad` and `hess`,
#' including compiled code (Rcpp, TMB) and RTMB, which differentiates R code.
#' An RTMB function can return the log-posterior itself, and its `obj$fn`,
#' `obj$gr` and `obj$he` are then passed as they are. For a TMB template that
#' returns the negative log-posterior, pass `function(x) -obj$fn(x)`,
#' `function(x) -as.vector(obj$gr(x))` and `function(x) -obj$he(x)`, with
#' `init = obj$par`; its repeated names (one per element of a parameter
#' vector) are numbered, as `beta[1]`, `beta[2]`, and so on. TMB gives no
#' Hessian for a model with random effects.
#'
#' @section The VB correction:
#' The VB mean correction moves the centre of the Laplace Gaussian towards the
#' posterior, keeping its covariance. It maximises the expected log-posterior
#' under that Gaussian, so that at the solution the gradient, averaged over the
#' Gaussian, is zero.
#'
#' `vb = "onestep"` takes the first Newton step towards that solution, with the
#' Laplace curvature, and stops. With a gradient, the average is taken over the
#' mode and \eqn{\pm\sqrt 3}{+/- sqrt(3)} posterior SDs along each principal
#' axis of the covariance. Without one, it is one pass of Stein's identity on
#' the rule of `vb = "full"`. Most of the movement is in that step, and exact VB
#' is itself an approximation to the posterior mean with an error of the same
#' order. If the first step is over 2.5 posterior SDs, or a node fails, the
#' iteration carries on with capped steps, within the number of calls that the
#' mode and Hessian took. `diagnostics()` then has `vb_step_sigma`, the step
#' still to go, which is `NA` after a single step. After a single step the fit
#' checks that it delivered what it promised: `vb_gain_ratio` in `diagnostics()`
#' is the rise in the expected log-posterior over the rise the step assumed, and
#' a ratio below 0.5 raises a warning. A large first step can also
#' overstate the VB term of the log marginal likelihood, so after a shift above
#' 1 posterior SD, which is flagged, compare models with care: `vb = "full"` may
#' be closer.
#'
#' `vb = "full"` iterates to the solution, on a degree-5 rule with at most
#' (m + 1)(m + 2) + 1 nodes, and hands over to [stats::nlminb()] if it does not
#' converge.
#'
#' @section Control:
#' `control` is a list with any of:
#' * `n_grid`: number of grid points per skew-normal scan, spread over
#'   \eqn{\pm 4} posterior SDs (default `21`).
#' * `log_threshold`: grid points whose log-posterior falls more than this
#'   below the maximum are ignored in the skew-normal fit (default `-6`).
#' * `vb_maxit`, `vb_tol`: iteration limit (default `25`) and convergence
#'   tolerance in posterior SDs (default `1e-3`) for the VB correction's
#'   fixed-point iteration, with `vb = "full"` and when a one-step correction
#'   has to iterate. A step over 1.5 posterior SDs is cut to 1.5. With
#'   `vb = "full"`, if two steps in a row are cut, or it has not converged
#'   after `vb_maxit` steps, [stats::nlminb()] takes over.
#' * `nlminb`: a list of [stats::nlminb()] controls for the mode search.
#'   `iter.max` and `eval.max` default to `1000` and `2000`, and `trace = 1`
#'   prints its iterations.
#'
#' @section Output:
#' By default each stage prints one line as it finishes, with the time it
#' took, and shows a progress bar while it runs for longer than half a
#' second, or than the option `cli.progress_show_after` when it is set. In
#' knitr or Quarto documents, where a bar cannot redraw in place, a larger
#' `cli.progress_show_after` in the setup chunk keeps bars out of the output.
#' Warnings, such as a gradient that does not match `logpost`, are signalled
#' as usual. A fit whose checks look off, such as a mode that was not
#' reached, ends with one warning that lists them, as INLAvaan's does (see
#' [summary.inla_vanilla()] and [diagnostics()]).
#'
#' `verbose = TRUE` adds a full diagnostic of the fit, in grey lines under
#' each stage. It reports what the fit computes anyway, with no further calls
#' of `logpost` or `grad`:
#' * Before the mode search: the parameters, where the derivatives come from,
#'   `logpost` at `init`, and for a supplied `grad` its largest relative
#'   discrepancy from finite differences there.
#' * The mode: the message, iterations and evaluations of [stats::nlminb()],
#'   and whether it was restarted with the parameters rescaled; the
#'   log-posterior at the mode, the whitened gradient there (which decides
#'   the restart) and the largest Newton step in posterior SDs (which the
#'   fit-time warning checks); how the Hessian was taken (and, for a
#'   supplied `grad`, the asymmetry of its Jacobian, or, for a supplied
#'   `hess`, how far it is from central differences of `grad`); the largest
#'   posterior
#'   correlation and the smallest eigenvalue of the correlation matrix, which
#'   is near zero when some combination of the parameters is poorly
#'   identified, and the condition number of the Hessian in posterior-SD
#'   units; and the Laplace log marginal likelihood.
#' * The VB correction: the largest step, in posterior SDs, at each
#'   fixed-point iteration, whether Anderson acceleration or the fallback
#'   optimiser was used (and, after the fallback, the step still to go where
#'   it landed, and whether it was rerun from there), the largest shift, and
#'   the term it adds to the log marginal likelihood.
#' * The marginals: the shape, NMAD, ridge bend, tails, unscanned mass and
#'   volume correction `gamma1` of each parameter, with the values that
#'   [summary.inla_vanilla()] flags marked.
#' * The calls of `logpost`, `grad` and `hess` in each stage, and in all,
#'   with the total time.
#'
#' For the iteration trace of [stats::nlminb()] itself, pass
#' `control = list(nlminb = list(trace = 1))`.
#'
#' `silent = TRUE` prints nothing, for simulations and other repeated fits:
#' no progress bars and no stage lines, and every warning and message raised
#' during the fit, including those of `logpost`, `grad` and `hess`, is
#' muffled. Errors are still raised. Output printed with [cat()] or
#' [print()], by `logpost`, `grad`, `hess` or the trace of [stats::nlminb()],
#' is not captured.
#'
#' Whatever the level, the warnings raised during the fit are recorded in
#' `fit$warnings`, so that a simulation can check them afterwards.
#'
#' @param logpost A function `logpost(theta, ...)` returning the log-posterior
#'   density (up to a constant) at the named numeric vector `theta`, on the
#'   unconstrained scale. See the section on writing `logpost`.
#' @param init A named numeric vector of starting values. Its length sets the
#'   number of parameters and its names label them; unnamed values are called
#'   `theta[1]`, `theta[2]`, ..., and a repeated name is numbered in the same
#'   way.
#' @param grad Optional function `grad(theta, ...)` returning the gradient of
#'   `logpost` as a numeric vector. When `NULL` (the default), finite
#'   differences are used.
#' @param hess Optional function `hess(theta, ...)` returning the Hessian of
#'   `logpost`, the m x m matrix of its second derivatives (a sparse matrix
#'   from the Matrix package, and a number when m = 1, work too). It needs
#'   `grad`. When `NULL` (the default), finite differences of `grad` are used.
#' @param ... Further arguments passed to `logpost`, `grad` and `hess`, such
#'   as the data. As in [stats::nlminb()], they come after the gradient and
#'   the Hessian, and the arguments after them (`vb`, `marginals`, `width`,
#'   `tk_fallback`, `verbose`, `silent`, `control`) must be given by their
#'   full names, so data arguments cannot partially match them. A data
#'   argument cannot take one of those names. The arguments before them can be
#'   matched, as in `nlminb()`: a data argument named `g` or `gr` would be
#'   taken as `grad`, and one named `h` or `he` as `hess`, so give it a
#'   longer name.
#' @param vb The VB mean correction: `"onestep"` (the default) takes one Newton
#'   step from the mode, `"full"` iterates to the VB solution, and `"off"` (or
#'   `FALSE`) leaves it out. `TRUE` means the default, and abbreviations such as
#'   `"one"` work too. See the section on the VB correction. The correction
#'   enters the log marginal likelihood, and the skew-normal means are checked
#'   against it (see [summary.inla_vanilla()]).
#' @param marginals The marginal of each parameter: `"skewnormal"` (the
#'   default) fits a skew-normal to a scan of the log-posterior; `"gaussian"`
#'   uses the Gaussian marginal of the Laplace approximation, centred at the
#'   VB mean unless `vb = "off"`, with no scans. Abbreviations such as `"skew"`
#'   and `"gauss"` work too.
#' @param width How each skew-normal scan accounts for the changing volume of
#'   the other parameters, the width term
#'   \eqn{\gamma_j(t) = -\frac12 \log|H_{-j}|}{gamma_j(t) = -1/2 log|H_-j|}
#'   (see the section on the width term). With m parameters and
#'   `n_grid` scan points, per parameter:
#'   * `"trace"` (the default): its linear term \eqn{\gamma_j' t}{gamma_j' t},
#'     as in INLAvaan, from 4 `grad` calls, and one more for all the parameters
#'     together (2m + 5 `logpost` calls without `grad`), on top of the
#'     `n_grid - 1` `logpost` calls of the scan.
#'   * `"cmp"`: the exact width term along the same straight scan, from the
#'     Hessian of the other parameters at each scan point: about
#'     `n_grid - 1` times 2(m - 1) `grad` calls (`n_grid - 1` `hess` calls
#'     with `hess`, or `n_grid - 1` times (m - 1)m `logpost` calls without
#'     `grad`).
#'   * `"tk"`: the Laplace marginal of Tierney and Kadane (1986). At each scan
#'     point the other parameters are moved to their conditional mode by
#'     [stats::nlminb()], started from the straight scan, and the height and
#'     the width term are taken there: the calls of `"cmp"` plus those of
#'     `n_grid - 1` conditional optimisations, typically 5 to 20 `logpost`
#'     and `grad` calls each with `grad`, and many more without.
#'
#'   With one parameter the three are the same. Abbreviations work. Ignored
#'   for Gaussian marginals.
#' @param tk_fallback Logical. Rescan with `width = "tk"` every parameter
#'   whose `ridge_bend` is flagged (above 0.25, see [summary.inla_vanilla()])?
#'   Default `FALSE`, so that a fit costs what `width` says and gives the same
#'   marginals whatever its diagnostics. The parameters rescanned are listed
#'   in `fit$fallback` and by [summary.inla_vanilla()]. Ignored for Gaussian
#'   marginals and with `width = "tk"`.
#' @param verbose Logical. Print a full diagnostic of each stage under its
#'   progress line? Default `FALSE`. See the section on output.
#' @param silent Logical. Print nothing at all, and muffle the warnings and
#'   messages raised during the fit? Default `FALSE`. The warnings are still
#'   recorded in `fit$warnings`. Cannot be combined with `verbose = TRUE`.
#' @param control A list of tuning parameters; see the section on control.
#'
#' @returns An object of class `inla_vanilla`, with methods for [print()],
#'   [summary()], [coef()], [plot()], [sampling()], [diagnostics()] and
#'   [timing()]. It is a list
#'   holding, among other things:
#'   * `summary`: the posterior summary table (mean, SD, quantiles, mode).
#'   * `mode`, `cov`: the Laplace mode and covariance.
#'   * `mean_vb`: the VB-corrected mean (the mode when `vb = "off"`).
#'   * `sn`: the skew-normal parameters `xi`, `omega` and `alpha` per
#'     parameter (with `alpha = 0` for Gaussian marginals).
#'   * `vb`, `marginals`, `width`, `tk_fallback`, `verbose`, `silent`: the
#'     settings the fit used (`width` is `NA` for Gaussian marginals).
#'   * `fallback`: the parameters that `tk_fallback` rescanned with
#'     `width = "tk"` (empty when there were none).
#'   * `logml`: the log marginal likelihood, and `logml_laplace` its Laplace
#'     part.
#'   * `diagnostics`: the checks of each stage, which [diagnostics()]
#'     returns under INLAvaan's names and [summary.inla_vanilla()] explains.
#'     `optim` holds the return code, message and iterations of
#'     [stats::nlminb()], the gradient of minus the log-posterior at the mode
#'     (`grad`), the whitened gradient there (`grad_mode`), and the condition
#'     number of the Hessian in posterior-SD units (`hess_cond`). `vb` holds
#'     the iterations of the VB correction, whether it stopped at the first
#'     step (`one_pass`), whether it fell back to
#'     [stats::nlminb()] (`fallback`), and the Kullback-Leibler divergence of
#'     the shift (`kld_global`); it is `NULL` when `vb = "off"`. `marginals`
#'     is a data frame with a row per parameter: `mode_shift_sigma`, `kld`,
#'     `vb_shift`, `vb_shift_sigma`, `nmad`, `alpha`, `unscanned_mass`,
#'     `vb_step_sigma`, `mean_gap`, `gamma1`, `tails` and `ridge_bend`, with
#'     `NA` for the checks the fit did not run, and `width`, the width method
#'     of each scan (see the section on the width term).
#'   * `scan`: the scans the skew-normals were fitted to, for
#'     `plot(fit, type = "sn_fit")`: the grid `z` in posterior SDs from the
#'     mode, the log-posterior at each point (`logpost`, a column per
#'     parameter: along the CMP, or at the conditional modes for
#'     `width = "tk"`), the width term added to it (`width`, zero at the
#'     mode), and the log scale factor `logC` of each fit on that grid.
#'     It is `NULL` for Gaussian marginals.
#'   * `timing`: the elapsed time of each stage of the fit, in seconds, which
#'     [timing()] returns.
#'   * `warnings`: the messages of the warnings raised during the fit, as a
#'     character vector of plain text, each message once (empty when there
#'     were none). They are recorded whether or not `silent` is `TRUE`.
#'
#' @references
#' Jamil, H., & Rue, H. (2026a). *Approximate Bayesian inference for
#' structural equation models using integrated nested Laplace
#' approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
#' <https://doi.org/10.48550/arXiv.2603.25690>
#'
#' Jamil, H., & Rue, H. (2026b). *Implementation and workflows for INLA-based
#' approximate Bayesian structural equation modelling* (arXiv:2604.00671
#' \[stat.CO\]). arXiv. <https://doi.org/10.48550/arXiv.2604.00671>
#'
#' Tierney, L., & Kadane, J. B. (1986). Accurate approximations for posterior
#' moments and marginal densities. *Journal of the American Statistical
#' Association*, *81*(393), 82–86.
#' <https://doi.org/10.1080/01621459.1986.10478240>
#'
#' @seealso [sampling()] for posterior draws and derived quantities,
#'   [summary.inla_vanilla()] and [diagnostics()] for diagnostics.
#'
#' @examples
#' # Normal data with unknown mean and standard deviation, on the scale
#' # (mu, log sigma)
#' set.seed(1)
#' y <- rnorm(30, mean = 5, sd = 2)
#'
#' logpost <- function(theta, y) {
#'   mu <- theta[["mu"]]
#'   sigma <- exp(theta[["log_sigma"]])
#'   sum(dnorm(y, mu, sigma, log = TRUE)) +
#'     dnorm(mu, 0, 10, log = TRUE) + # prior on mu
#'     dexp(sigma, rate = 1, log = TRUE) + # prior on sigma ...
#'     theta[["log_sigma"]] # ... and its Jacobian
#' }
#'
#' fit <- inla_vanilla(logpost, init = c(mu = 0, log_sigma = 0), y = y)
#' fit
#'
#' # The posterior of sigma itself, through the draws
#' sigma <- sampling(fit, fun = function(theta) exp(theta[["log_sigma"]]))
#' quantile(sigma, c(0.025, 0.5, 0.975))
#'
#' # A full diagnostic of each stage
#' fit <- inla_vanilla(
#'   logpost,
#'   init = c(mu = 0, log_sigma = 0),
#'   y = y,
#'   verbose = TRUE
#' )
#'
#' # Nothing printed, as in a simulation; any warnings are kept in the fit
#' fit <- inla_vanilla(
#'   logpost,
#'   init = c(mu = 0, log_sigma = 0),
#'   y = y,
#'   silent = TRUE
#' )
#' fit$warnings
#' @export
inla_vanilla <- function(
  logpost,
  init,
  grad = NULL,
  hess = NULL,
  ...,
  vb = c("onestep", "full", "off"),
  marginals = c("skewnormal", "gaussian"),
  width = c("trace", "cmp", "tk"),
  tk_fallback = FALSE,
  verbose = FALSE,
  silent = FALSE,
  control = list()
) {
  t_start <- elapsed()
  call <- match.call()

  check_function(logpost, "logpost")
  if (!is.null(grad)) {
    check_function(grad, "grad")
  }
  if (!is.null(hess)) {
    check_function(hess, "hess")
    if (is.null(grad)) {
      abort(c(
        "{.arg hess} needs {.arg grad}.",
        "i" = "Supply the gradient too, or leave out {.arg hess}. As in
               {.fn stats::nlminb}, a Hessian is used only with a gradient."
      ))
    }
  }
  vb <- check_vb(vb)
  marginals <- check_choice(marginals, c("skewnormal", "gaussian"), "marginals")
  width <- check_choice(width, c("trace", "cmp", "tk"), "width")
  check_flag(tk_fallback, "tk_fallback")
  check_flag(verbose, "verbose")
  check_flag(silent, "silent")
  if (verbose && silent) {
    abort(c(
      "{.arg verbose} and {.arg silent} cannot both be {.code TRUE}.",
      "i" = "{.code verbose = TRUE} prints a full diagnostic of the fit, and
             {.code silent = TRUE} prints nothing."
    ))
  }
  init <- check_init(init, show = !silent)
  ctrl <- vanilla_control(control)

  # The user's functions with `...` bound once here: forwarding `...` to
  # internal helpers would let a data argument such as `g` partially match
  # one of their formals.
  logpost_x <- function(x) logpost(x, ...)
  grad_x <- if (is.null(grad)) NULL else function(x) grad(x, ...)
  hess_x <- if (is.null(hess)) NULL else function(x) hess(x, ...)

  # Every warning raised during the fit, by the package or by logpost, grad
  # and hess, is recorded once in fit$warnings. With silent = TRUE it is then
  # muffled, and so is every message; otherwise it is signalled as usual.
  # Errors pass through either way.
  prog <- progress_new(show = !silent, detail = verbose)
  withCallingHandlers(
    vanilla_fit(
      logpost_x,
      grad_x,
      hess_x,
      init,
      vb = vb,
      marginals = marginals,
      width = width,
      tk_fallback = tk_fallback,
      ctrl = ctrl,
      prog = prog,
      call = call,
      t_start = t_start
    ),
    warning = function(w) {
      msg <- warning_text(w)
      if (!msg %in% prog$warnings) {
        prog$warnings <- c(prog$warnings, msg)
      }
      if (silent) {
        # tryInvokeRestart(), as suppressWarnings() does: a condition
        # signalled with signalCondition() has no muffle restart
        tryInvokeRestart("muffleWarning")
      }
    },
    message = function(m) {
      if (silent) {
        tryInvokeRestart("muffleMessage")
      }
    }
  )
}

# The wrappers the stages call in place of the user's functions, which have
# `...` already bound. They name the parameters, and count the calls and record
# the last parameter value in the environment `state`. `on_call` runs at each
# call, for the progress bar. `gr` and `he` are NULL when the user gave none.
wrap_fns <- function(
  logpost_x,
  grad_x,
  hess_x,
  init,
  on_call = function() NULL
) {
  m <- length(init)
  parnames <- names(init)
  state <- new.env(parent = emptyenv())
  state$n_lp <- 0L
  state$n_gr <- 0L
  state$n_he <- 0L
  state$last_x <- init
  lp <- function(x) {
    names(x) <- parnames
    state$last_x <- x
    state$n_lp <- state$n_lp + 1L
    on_call()
    # [[1L]] also accepts a 1 x 1 matrix, such as a quadratic form
    v <- logpost_x(x)[[1L]]
    if (is.na(v)) -Inf else v
  }
  gr <- if (is.null(grad_x)) {
    NULL
  } else {
    function(x) {
      names(x) <- parnames
      state$last_x <- x
      state$n_gr <- state$n_gr + 1L
      on_call()
      as.numeric(grad_x(x))
    }
  }
  # A plain m x m matrix, whatever the user's function returns (a sparse
  # matrix from the Matrix package, or a number when m = 1)
  he <- if (is.null(hess_x)) {
    NULL
  } else {
    function(x) {
      names(x) <- parnames
      state$last_x <- x
      state$n_he <- state$n_he + 1L
      on_call()
      matrix(as.numeric(as.matrix(hess_x(x))), m, m)
    }
  }
  list(lp = lp, gr = gr, he = he, state = state)
}

# The stages of the fit, given the user's functions with `...` already bound,
# the checked arguments, and the progress reporter `prog`, which also holds
# the record of warnings.
vanilla_fit <- function(
  logpost_x,
  grad_x,
  hess_x,
  init,
  vb,
  marginals,
  width,
  tk_fallback,
  ctrl,
  prog,
  call,
  t_start
) {
  m <- length(init)
  parnames <- names(init)
  has_grad <- !is.null(grad_x)
  has_hess <- !is.null(hess_x)

  # The engine calls the wrappers of wrap_fns(), and the progress bar ticks
  # from them. The spinner counts the stage's logpost, grad and hess calls
  # together, from `stage_calls()`, so stages driven by grad tick too.
  n_base <- 0L
  stage_calls <- function() {
    n_base <<- st$n_lp + st$n_gr + st$n_he
  }
  tick <- function() {
    n <- st$n_lp + st$n_gr + st$n_he
    if (prog$calls && n %% 16L == 0L) {
      progress_tick(prog, set = n - n_base)
    }
  }
  fns <- wrap_fns(logpost_x, grad_x, hess_x, init, on_call = tick)
  lp <- fns$lp
  gr <- fns$gr
  he <- fns$he
  st <- fns$state
  # The calls since the last stage, for verbose output
  calls_seen <- c(0L, 0L, 0L)
  calls_since <- function() {
    now <- c(st$n_lp, st$n_gr, st$n_he)
    out <- now - calls_seen
    calls_seen <<- now
    out
  }
  # An error inside logpost, grad or hess stops the fit; say where it happened
  # and at which parameter value.
  with_context <- function(stage, expr) {
    tryCatch(expr, error = function(e) {
      if (inherits(e, "inlavanilla_error")) {
        stop(e)
      }
      at <- paste(parnames, "=", signif(st$last_x, 4), collapse = ", ")
      abort(
        c(
          "The fit stopped with an error while {stage}.",
          "i" = "The last parameter value tried was {at}.",
          "i" = "If {.arg logpost} cannot be evaluated somewhere, return
                 {.code -Inf} there instead of failing."
        ),
        parent = e,
        call = quote(inla_vanilla())
      )
    })
  }

  lp_init <- with_context(
    "evaluating logpost at init",
    check_logpost_at_init(logpost_x, init)
  )
  grad_check <- if (has_grad) {
    with_context(
      "evaluating grad at init",
      check_grad_at_init(grad_x, lp, init)
    )
  }
  if (has_hess) {
    with_context(
      "evaluating hess at init",
      check_hess_at_init(hess_x, init)
    )
  }
  if (prog$detail) {
    report_start(prog, init, lp_init, grad_check, calls_since(), has_hess)
  }

  # 1a. Joint Gaussian -------------------------------------------------------
  # Each segment of fit$timing runs from its start to the start of the next,
  # so the segments add up to the total. The progress lines time the work of
  # each stage alone, from `t0`.
  t0 <- elapsed()
  t_starts <- c(init = t_start, laplace = t0)
  stage_calls()
  progress_start(prog, "Finding the posterior mode, {cli::pb_current} calls")
  lap <- with_context(
    "finding the posterior mode",
    laplace_fit(lp, gr, init, ctrl, he)
  )
  if (has_grad && !has_hess && lap$asymmetry > 1e-3) {
    warn(c(
      "{.arg grad} does not look like the gradient of {.arg logpost}.",
      "x" = "Its Jacobian at the mode is not symmetric (relative asymmetry
             {signif(lap$asymmetry, 2)}), as the Hessian of a function must
             be.",
      "i" = "A wrong gradient gives a wrong mode and wrong marginals."
    ))
  }
  if (has_hess && !isTRUE(lap$hess_error <= 1e-3)) {
    warn(c(
      "{.arg hess} does not look like the Hessian of {.arg logpost}.",
      "x" = "At the mode it differs from central differences of {.arg grad}
             by {signif(lap$hess_error, 2)}, in units where the Hessian is
             the identity.",
      "i" = "A wrong Hessian gives a wrong covariance and wrong marginals."
    ))
  }
  # INLAvaan's Newton step at the mode, in posterior SDs, decides whether the
  # mode is marked here and reported in the fit-time warning, not nlminb()'s
  # return code: without a gradient nlminb() can report "false convergence"
  # at an accurate mode, and it can report success short of one.
  mode_shift_max <- max(lap$mode_shift)
  mode_ok <- is.finite(mode_shift_max) &&
    mode_shift_max <= warn_levels$mode_shift
  t_laplace <- elapsed() - t0
  progress_done(
    prog,
    if (mode_ok) {
      paste0(
        "Found the posterior mode in ",
        hl(lap$opt$iterations),
        " iteration",
        if (lap$opt$iterations == 1L) "" else "s"
      )
    } else {
      paste0(
        "Stopped near the posterior mode (Newton step ",
        hl(with_sd(mode_shift_max)),
        ")"
      )
    },
    t_laplace,
    ok = mode_ok
  )
  if (prog$detail) {
    report_mode(prog, lap, calls_since())
  }

  sds <- sqrt(diag(lap$Sigma))

  # 1b. VB mean correction ---------------------------------------------------
  # With the user's gradient when there is one; without, from log-posterior
  # values alone (see vb_correct())
  t0 <- elapsed()
  t_starts[["vb"]] <- t0
  vbc <- NULL
  use_vb <- vb != "off"
  if (use_vb) {
    stage_calls()
    progress_start(
      prog,
      "Applying the VB mean correction, {cli::pb_current} calls"
    )
    vbc <- with_context(
      "applying the VB mean correction",
      vb_correct(
        lp,
        gr,
        lap$mode,
        lap$Sigma,
        lap$R,
        lap$L,
        ctrl,
        method = vb,
        budget = st$n_lp + st$n_gr + st$n_he
      )
    )
    # The fallback is a route of the solver, reported but not marked, as in
    # INLAvaan. A correction that did not settle is marked.
    vb_max <- max(abs(vbc$shift) / sds)
    gain_ok <- !isTRUE(below_level(vbc$gain_ratio, warn_levels$vb_gain))
    settled <- !any(
      check_ran(vbc$step_sigma) &
        off_level(vbc$step_sigma, warn_levels$vb_step)
    )
    progress_done(
      prog,
      paste0(
        "VB mean correction, max shift ",
        hl(paste0(format(signif(vb_max, 2)), sd_unit())),
        if (vbc$fallback) " (fallback optimiser)" else "",
        if (!settled) ", not settled" else "",
        if (!gain_ok) ", gain below its promise" else ""
      ),
      elapsed() - t0,
      ok = settled &&
        gain_ok &&
        is.finite(vb_max) &&
        vb_max <= warn_levels$vb_shift
    )
    if (prog$detail) {
      report_vb(prog, vbc, sds, lap$logml, calls_since())
    }
  } else {
    progress_info(prog, "No VB correction (vb = \"off\")")
  }
  shift <- if (use_vb) vbc$shift else stats::setNames(numeric(m), parnames)
  mean_vb <- lap$mode + shift

  # 2. Skew-normal marginals -------------------------------------------------
  # With Gaussian marginals there are no scans: each marginal is the
  # Gaussian's own, N(mean_vb_j, Sigma_jj), a skew-normal with shape 0.
  t0 <- elapsed()
  t_starts[["marginals"]] <- t0
  if (marginals == "skewnormal") {
    start_marginals <- function(step) {
      progress_close(prog)
      if (step == "shared") {
        stage_calls()
        progress_start(
          prog,
          "Computing the volume corrections, {cli::pb_current} calls"
        )
      } else {
        progress_start(
          prog,
          paste0(
            "Fitting ",
            pb_item,
            " skew-normal marginal",
            if (m == 1L) "" else "s"
          ),
          total = m
        )
      }
    }
    mar <- with_context(
      "fitting the skew-normal marginals",
      fit_marginals(
        lp,
        gr,
        lap,
        ctrl,
        tick = function() progress_tick(prog),
        start = start_marginals,
        he = he,
        width = width,
        tk_fallback = tk_fallback,
        par = list(
          logpost = logpost_x,
          grad = grad_x,
          hess = hess_x,
          init = init,
          # Seconds per call so far, to judge how long the scans will take.
          # The clock has a resolution of a millisecond
          call_s = max(elapsed() - t_start, 1e-3) /
            (st$n_lp + st$n_gr + st$n_he),
          add_calls = function(n) {
            st$n_lp <- st$n_lp + n[[1L]]
            st$n_gr <- st$n_gr + n[[2L]]
            st$n_he <- st$n_he + n[[3L]]
          },
          # A message cannot share the console with a running progress line
          note = function(text) {
            progress_close(prog)
            progress_note(prog, text)
          }
        )
      )
    )
    workers <- mar$workers
    sn <- cbind(
      mar$sn,
      unscanned_mass = unscanned_mass(mar$sn, lap$mode, sds)
    )
    # The scans, kept for plot(type = "sn_fit"): no further calls
    scan <- list(
      z = mar$z,
      logpost = mar$logpost,
      width = mar$width,
      logC = sn[, "logC"]
    )
    method <- mar$method
    progress_done(
      prog,
      paste0(
        "Fitted ",
        hl(m),
        " skew-normal marginal",
        if (m == 1L) "" else "s"
      ),
      elapsed() - t0
    )
    if (prog$detail) {
      report_marginals(prog, sn, calls_since(), method, width)
      if (workers > 1L) {
        progress_info(prog, paste0("The scans ran on ", workers, " workers"))
      }
    }
  } else {
    workers <- 1L
    na <- rep(NA_real_, m)
    sn <- cbind(
      xi = mean_vb,
      omega = sds,
      alpha = 0,
      gamma1 = na,
      nmad = na,
      tails = na,
      ridge_bend = na,
      logC = na,
      unscanned_mass = na
    )
    rownames(sn) <- parnames
    scan <- NULL
    method <- stats::setNames(rep(NA_character_, m), parnames)
    progress_info(
      prog,
      "Gaussian marginals (marginals = \"gaussian\"): no scans"
    )
  }
  t_starts[["diagnostics"]] <- elapsed()

  logml <- lap$logml + if (use_vb) vbc$kld_global else 0
  mom <- sn_moments(sn[, "xi"], sn[, "omega"], sn[, "alpha"])

  # The diagnostics, in the names and order of INLAvaan's diagnostics() where
  # it has the quantity. Checks that did not run are NA: the VB columns with
  # vb = "off", the scan columns with Gaussian marginals.
  skew <- marginals == "skewnormal"
  na <- rep(NA_real_, m)
  diagnostics <- list(
    optim = list(
      convergence = lap$opt$convergence,
      message = lap$opt$message,
      iterations = lap$opt$iterations,
      # INLAvaan's sign: the gradient of nlminb()'s objective, -logpost
      grad = -lap$grad,
      grad_mode = lap$grad_mode,
      hess_cond = lap$hess_cond
    ),
    vb = if (use_vb) {
      list(
        iterations = vbc$iterations,
        fallback = vbc$fallback,
        one_pass = vbc$one_pass,
        gain_ratio = vbc$gain_ratio,
        kld_global = vbc$kld_global
      )
    },
    marginals = data.frame(
      mode_shift_sigma = lap$mode_shift,
      kld = if (use_vb) 0.5 * (shift / sds)^2 else na,
      vb_shift = if (use_vb) shift else na,
      vb_shift_sigma = if (use_vb) shift / sds else na,
      nmad = sn[, "nmad"],
      alpha = if (skew) sn[, "alpha"] else na,
      unscanned_mass = sn[, "unscanned_mass"],
      vb_step_sigma = if (use_vb) vbc$step_sigma else na,
      mean_gap = if (use_vb && skew) (mom$mean - mean_vb) / sds else na,
      gamma1 = sn[, "gamma1"],
      tails = as.logical(sn[, "tails"]),
      ridge_bend = sn[, "ridge_bend"],
      width = unname(method),
      row.names = parnames,
      stringsAsFactors = FALSE
    )
  )

  # One warning for every check that looks off, before the fit is built so
  # that fit$warnings records it
  warn_diagnostics(
    diagnostics_internal(
      diagnostics$optim,
      diagnostics$vb,
      diagnostics$marginals,
      pars = mean_vb
    ),
    gaussian = !skew,
    has_grad = has_grad
  )

  t_end <- elapsed()
  t_total <- t_end - t_start
  if (prog$detail) {
    report_end(prog, st$n_lp, st$n_gr, st$n_he, t_total)
  }
  new_inla_vanilla(
    sn = sn[, c("xi", "omega", "alpha"), drop = FALSE],
    mode = lap$mode,
    cov = lap$Sigma,
    mean_vb = mean_vb,
    logml = logml,
    logml_laplace = lap$logml,
    diagnostics = diagnostics,
    scan = scan,
    counts = c(logpost = st$n_lp, grad = st$n_gr, hess = st$n_he),
    workers = workers,
    timing = c(
      stats::setNames(diff(c(t_starts, t_end)), names(t_starts)),
      total = t_total
    ),
    vb = vb,
    marginals = marginals,
    width = if (skew) width else NA_character_,
    tk_fallback = tk_fallback,
    fallback = if (skew) parnames[method != width] else character(),
    has_grad = has_grad,
    verbose = prog$detail,
    silent = !prog$show,
    warnings = prog$warnings,
    call = call
  )
}

new_inla_vanilla <- function(
  sn,
  mode,
  cov,
  mean_vb,
  logml,
  logml_laplace,
  diagnostics,
  scan,
  counts,
  workers,
  timing,
  vb,
  marginals,
  width,
  tk_fallback,
  fallback,
  has_grad,
  verbose,
  silent,
  warnings,
  call
) {
  structure(
    list(
      summary = sn_summary(sn),
      sn = sn,
      mode = mode,
      cov = cov,
      mean_vb = mean_vb,
      logml = logml,
      logml_laplace = logml_laplace,
      diagnostics = diagnostics,
      scan = scan,
      counts = counts,
      workers = workers,
      timing = timing,
      vb = vb,
      marginals = marginals,
      width = width,
      tk_fallback = tk_fallback,
      fallback = fallback,
      has_grad = has_grad,
      verbose = verbose,
      silent = silent,
      warnings = warnings,
      call = call
    ),
    class = "inla_vanilla"
  )
}

# Posterior summary table from the skew-normal marginals: moments in closed
# form, quantiles by inversion, and the mode numerically.
sn_summary <- function(sn, probs = c(0.025, 0.5, 0.975)) {
  m <- nrow(sn)
  mom <- sn_moments(sn[, "xi"], sn[, "omega"], sn[, "alpha"])
  qq <- matrix(NA_real_, m, length(probs))
  md <- numeric(m)
  for (j in seq_len(m)) {
    qq[j, ] <- qsnorm(probs, sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"])
    md[j] <- sn_mode(sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"])
  }
  tab <- cbind(mom$mean, mom$sd, qq, md)
  dimnames(tab) <- list(
    rownames(sn),
    c("Mean", "SD", format_probs(probs), "Mode")
  )
  tab
}

format_probs <- function(probs) {
  # as stats::quantile() names them, with no padding
  sprintf("%s%%", formatC(100 * probs, format = "fg", width = 1, digits = 7))
}

vanilla_control <- function(control) {
  defaults <- list(
    n_grid = 21L,
    log_threshold = -6,
    vb_maxit = 25L,
    vb_tol = 1e-3,
    nlminb = list()
  )
  if (!is.list(control)) {
    abort("{.arg control} must be a list.")
  }
  nms <- names(control)
  if (length(control) > 0L && (is.null(nms) || any(is.na(nms) | nms == ""))) {
    abort(c(
      "Every {.arg control} entry must be named.",
      "i" = "Valid entries are {.val {names(defaults)}}."
    ))
  }
  unknown <- setdiff(nms, names(defaults))
  if (length(unknown) > 0L) {
    abort(c(
      "Unknown {.arg control} entr{?y/ies}: {.val {unknown}}.",
      "i" = "Valid entries are {.val {names(defaults)}}."
    ))
  }
  ctrl <- utils::modifyList(defaults, control)
  is_num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is_num(ctrl$n_grid) || ctrl$n_grid < 7) {
    abort("{.code control$n_grid} must be a single number of at least 7.")
  }
  ctrl$n_grid <- as.integer(round(ctrl$n_grid))
  if (!is_num(ctrl$log_threshold) || ctrl$log_threshold >= 0) {
    abort("{.code control$log_threshold} must be a single negative number.")
  }
  if (!is_num(ctrl$vb_maxit) || ctrl$vb_maxit < 1) {
    abort("{.code control$vb_maxit} must be a single number of at least 1.")
  }
  ctrl$vb_maxit <- as.integer(round(ctrl$vb_maxit))
  if (!is_num(ctrl$vb_tol) || ctrl$vb_tol <= 0) {
    abort("{.code control$vb_tol} must be a single positive number.")
  }
  if (!is.list(ctrl$nlminb)) {
    abort(
      "{.code control$nlminb} must be a list of {.fn stats::nlminb} controls."
    )
  }
  ctrl$nlminb <- utils::modifyList(
    list(iter.max = 1000L, eval.max = 2000L),
    ctrl$nlminb
  )
  ctrl
}

check_function <- function(f, arg) {
  if (!is.function(f)) {
    abort(
      "{.arg {arg}} must be a function, not {.obj_type_friendly {f}}."
    )
  }
}

check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    abort("{.arg {arg}} must be {.code TRUE} or {.code FALSE}.")
  }
}

# The VB method: one of "onestep", "full" and "off", where TRUE is the default
# and FALSE is "off"
check_vb <- function(vb) {
  if (isTRUE(vb)) {
    return("onestep")
  }
  if (isFALSE(vb)) {
    return("off")
  }
  check_choice(vb, c("onestep", "full", "off"), "vb")
}

# One of `choices` (lower case), or an unambiguous abbreviation of one in any
# case, as with match.arg(): "skew" and "Gauss" both work. The whole vector
# (the default) means the first. The package's own error otherwise.
check_choice <- function(x, choices, arg) {
  if (identical(x, choices)) {
    return(choices[[1L]])
  }
  i <- NA_integer_
  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
    i <- pmatch(tolower(x), choices)
  }
  if (is.na(i)) {
    abort("{.arg {arg}} must be one of {.or {.val {choices}}}.")
  }
  choices[[i]]
}

check_init <- function(init, show = TRUE) {
  if (!is.numeric(init) || length(init) == 0L) {
    abort(
      "{.arg init} must be a non-empty numeric vector, not {.obj_type_friendly {init}}."
    )
  }
  if (!all(is.finite(init))) {
    abort("{.arg init} must contain only finite values.")
  }
  nms <- names(init)
  if (is.null(nms)) {
    nms <- paste0("theta[", seq_along(init), "]")
  } else if (any(is.na(nms) | nms == "")) {
    abort("{.arg init} must be fully named, or not named at all.")
  } else if (anyDuplicated(nms)) {
    # As TMB names obj$par: one repeated name per element of a vector
    dup <- unique(nms[duplicated(nms)])
    renamed <- nms %in% dup
    for (d in dup) {
      i <- which(nms == d)
      nms[i] <- paste0(d, "[", seq_along(i), "]")
    }
    if (anyDuplicated(nms)) {
      abort("{.arg init} has names that cannot be made unique.")
    }
    if (show) {
      message(cli::format_message(c(
        "i" = "Numbered the repeated names in {.arg init}: {.val {nms[renamed]}}."
      )))
    }
  }
  stats::setNames(as.numeric(init), nms)
}

check_logpost_at_init <- function(logpost, init) {
  v <- logpost(init)
  if (!is.numeric(v) || length(v) != 1L) {
    abort(
      c(
        "{.arg logpost} must return a single number.",
        "x" = "{.code logpost(init)} returned {.obj_type_friendly {v}} of length {length(v)}."
      )
    )
  }
  if (!is.finite(v)) {
    abort(
      c(
        "{.arg logpost} is not finite at {.arg init}.",
        "x" = "{.code logpost(init)} returned {.val {v}}.",
        "i" = "Start inside the support of the posterior, and check that every
             parameter is on an unconstrained scale."
      )
    )
  }
  invisible(v)
}

# A Hessian of the right shape, finite at init. Whether it is the Hessian of
# `logpost` is checked at the mode, against the gradient (hess_discrepancy()).
check_hess_at_init <- function(hess, init) {
  m <- length(init)
  h <- hess(init)
  H <- tryCatch(as.matrix(h), error = function(e) NULL)
  if (!is.numeric(H) || !identical(dim(H), c(m, m))) {
    shape <- if (is.null(dim(h))) {
      paste("of length", length(h))
    } else {
      paste("with dimensions", paste(dim(h), collapse = " x "))
    }
    abort(c(
      "{.arg hess} must return a numeric {m} x {m} matrix.",
      "x" = "{.code hess(init)} returned {.obj_type_friendly {h}} {shape}."
    ))
  }
  if (!all(is.finite(H))) {
    abort("{.arg hess} is not finite at {.arg init}.")
  }
  invisible(H)
}

check_grad_at_init <- function(grad, lp, init) {
  g <- grad(init)
  if (!is.numeric(g) || length(g) != length(init)) {
    abort(
      c(
        "{.arg grad} must return a numeric vector of length {length(init)}.",
        "x" = "{.code grad(init)} returned {.obj_type_friendly {g}} of length {length(g)}."
      )
    )
  }
  if (!all(is.finite(g))) {
    abort(
      "{.arg grad} is not finite at {.arg init}.",
    )
  }
  # Steps from the local curvature scale, so the check holds for posterior
  # SDs far from 1; a fixed step where the curvature search fails
  scales <- fd_scales(lp, init, lp(init))
  h <- ifelse(attr(scales, "ok"), 1e-3 * scales, 1e-5 * pmax(1, abs(init)))
  g_fd <- fd_grad(lp, init, h = h)
  # Compare per curvature scale, the change in log-posterior over about one
  # posterior SD, so a gradient that is wrong by a factor is caught whatever
  # the scale of the parameter
  u <- ifelse(attr(scales, "ok"), scales, 1)
  score <- function(g, g_fd) abs(g - g_fd) * u / (1 + abs(g_fd) * u)
  err <- score(g, g_fd)
  where <- "at {.arg init}"
  # Near the mode both gradients are close to zero, so a gradient that is
  # wrong by a factor can still pass there. Check again one curvature scale
  # away, where the gradient is about one per scale; skip the point if
  # logpost or grad cannot be evaluated there.
  x1 <- init + u
  if (is.finite(lp(x1))) {
    g1 <- tryCatch(as.numeric(grad(x1)), error = function(e) NULL)
    if (length(g1) == length(init) && all(is.finite(g1))) {
      g1_fd <- fd_grad(lp, x1, h = h)
      err1 <- score(g1, g1_fd)
      if (all(is.finite(err1)) && max(err1) > max(err)) {
        err <- err1
        g <- g1
        g_fd <- g1_fd
        where <- "one curvature scale from {.arg init}"
      }
    }
  }
  if (any(err > 1e-3)) {
    worst <- which.max(err)
    warn(c(
      paste0(
        "{.arg grad} does not match finite differences of {.arg logpost} ",
        where,
        "."
      ),
      "x" = "For {.field {names(init)[worst]}}: {.arg grad} gives
             {.val {signif(g[worst], 4)}}, finite differences give
             {.val {signif(g_fd[worst], 4)}}.",
      "i" = "A wrong gradient gives a wrong mode and wrong marginals."
    ))
  }
  # The largest discrepancy, for verbose output
  invisible(max(err))
}
