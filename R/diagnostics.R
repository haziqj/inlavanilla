# The diagnostics() generic and its methods, and the fit-time warning. The
# checks themselves are in R/vanilla-diagnostics.R.

# The wording of the diagnostics for inla_vanilla(), in its own arguments (see
# advice_entries in R/vanilla-diagnostics.R), and the vignette that explains
# them
inlavanilla_advice <- list(
  grad_gap = "{.arg grad} does not match {.arg logpost}",
  no_grad = "supply {.arg grad} if you can",
  vb_full = "{.code vb_correction = \"full\"}",
  vb_none = "{.code vb_correction = \"none\"}",
  mode_limit = paste0(
    "Raise {.code mode_maxit} or {.code mode_maxeval} in ",
    "{.fn vanilla_control}."
  ),
  mode_retry = "Try other starting values, or supply {.arg grad}.",
  funnel = paste0(
    "Random effects in {.arg logpost} whose scale can run to zero do this ",
    "(a funnel): integrate them out of {.arg logpost}."
  ),
  vignette = "{.code vignette(\"diagnostics\", package = \"inlavanilla\")}"
)

# One warning, of class inlavanilla_diagnostics_warning, when any check is
# off; FALSE (invisibly) when none is.
warn_diagnostics <- function(dg, gaussian, has_grad = TRUE) {
  problems <- diagnostic_problems(dg, gaussian, has_grad, inlavanilla_advice)
  if (length(problems) == 0L) {
    return(invisible(FALSE))
  }
  env <- attr(problems, "envir")
  env$n_issues <- sum(names(problems) == "x")
  warn(
    c(
      "Fit diagnostics flagged {n_issues} potential issue{?s}:",
      problems,
      "i" = paste0(
        "Inspect with {.code diagnostics(fit)} and ",
        "{.code diagnostics(fit, type = \"param\")}; ",
        inlavanilla_advice$vignette,
        " explains each check."
      )
    ),
    class = "inlavanilla_diagnostics_warning",
    .envir = env
  )
  invisible(TRUE)
}

#' Diagnostics of an INLAvanilla fit
#'
#' @description
#' The checks of each stage of the approximation, as numbers: a summary of the
#' whole fit, or a table with one row per parameter. They are the checks of
#' INLAvaan's `diagnostics()`, in the same order and under the same names, two
#' of them under the names INLAvaan will use from its version 0.4.0
#' (`hess_cond_cor`, and `unscanned_mass` for its `scan_end_mass`), followed by
#' inlavanilla's own. The fit-time warning and the flags of
#' [summary.inla_vanilla()] are computed from the same numbers, and
#' [summary.inla_vanilla()] explains each check and the level at which it is
#' flagged.
#'
#' Jamil and Rue (2026) use the VB shift and NMAD as fit diagnostics
#' (Sections 3.1 and 4.4), and the gap between the skew-normal and VB means
#' (`mean_gap`) as an internal check (Section 5.1). The mode shift, the
#' Hessian condition number and the unscanned mass come from INLAvaan's
#' `diagnostics()`, not from the paper. `grad_gap`, `vb_step_sigma`, `tails`
#' and `ridge_bend` are inlavanilla's own.
#'
#' @details
#' A posterior SD here is the Laplace SD at the mode,
#' \eqn{s_j = \sqrt{\Sigma_{jj}}}{s_j = sqrt(Sigma_jj)}, as in INLAvaan; it
#' is not the SD column of the summary table, which is the skew-normal's.
#' Values that a fit did not compute are `NA`: the VB quantities with
#' `vb_correction = "none"`, the scan quantities with `marginals = "gaussian"`,
#' and `grad_gap` without `grad`.
#'
#' **Global diagnostics** (`type = "global"`), a named numeric vector:
#' * `npar`: the number of parameters.
#' * `converged`: 1 if [stats::nlminb()] reported convergence (return code
#'   0), 0 otherwise, as INLAvaan defines it. No check uses it: without a
#'   gradient `nlminb()` can report false convergence at an accurate mode.
#'   `mode_shift_max` checks the mode instead.
#' * `iterations`: the iterations of `nlminb()`, the total over both runs if
#'   it was restarted.
#' * `grad_inf`, `grad_inf_rel`, `grad_l2`: the largest absolute value, the
#'   largest relative value and the Euclidean norm of `grad` (see below).
#'   They depend on the scale of each parameter.
#' * `mode_shift_max`: the largest `mode_shift_sigma`. Unlike the gradient
#'   norms it is scale-free. It should be close to 0.
#' * `hess_cond`: the condition number of the Hessian at the mode, the ratio
#'   of the largest to the smallest eigenvalue of the Laplace covariance, as
#'   INLAvaan defines it. It depends on the units of the parameters: posterior
#'   SDs of 1e-3 and 1e3 alone make it 1e12.
#' * `hess_cond_cor`: the same with each parameter in its own posterior-SD
#'   units, the ratio of the largest to the smallest eigenvalue of the
#'   posterior correlation matrix. It is free of units. Large values mean that
#'   some combination of the parameters is weakly identified, and the fit
#'   warns above \eqn{10^8}.
#' * `vb_applied`: 1 if the VB correction ran, 0 otherwise.
#' * `vb_shift_max`: the largest absolute `vb_shift_sigma`.
#' * `vb_kld_global`: the Kullback-Leibler divergence between the Laplace
#'   Gaussian and its VB-shifted copy,
#'   \eqn{\frac12 \delta^\top H \delta}{delta' H delta / 2}, which is also
#'   the VB term of the log marginal likelihood.
#' * `kld_max`, `kld_mean`: the largest and the mean `kld`.
#' * `nmad_max`, `nmad_mean`: the largest and the mean `nmad`.
#' * `unscanned_mass_max`: the largest `unscanned_mass`.
#'
#' and then inlavanilla's own:
#' * `grad_gap_max`: the largest `grad_gap`. It should be close to 0.
#' * `vb_fallback`: 1 if the VB fixed-point iteration handed over to the
#'   fallback optimiser, 0 otherwise. Like INLAvaan, inlavanilla records it
#'   but flags nothing: the fallback is a route of the solver, not a measure
#'   of the approximation.
#' * `vb_step_max`: the largest `vb_step_sigma`. Like `mode_shift_max` it
#'   should be close to 0.
#' * `vb_gain_ratio`: how much of the promised rise in the expected
#'   log-posterior a one-step correction delivered: the rise under the
#'   Gaussian, over the \eqn{\frac12 \delta^\top H \delta}{1/2 d' H d} that
#'   the step assumes. Near 1 is healthy, and it is flagged below 0.5. It
#'   checks a one-step correction without a gradient, and is `NA` with one,
#'   when the correction iterated, or when the step was too small to check.
#' * `vb_lhat`, `vb_j0`, `vb_j1`: the check of a one-step correction with a
#'   gradient, from a second sweep of gradients at the step. One step is the
#'   first iterate of a fixed-point map, and holds only where that map is a
#'   contraction. `vb_lhat` is the size of the residual of the VB equation
#'   at the step over its size at the mode, which is the contraction ratio of
#'   the map over the step. `vb_j0` and `vb_j1` are the slope of the map
#'   along the step, at the mode and at the step. The check passes when
#'   `vb_lhat` is below 1 and both slopes lie between -1 and 1. A value that
#'   is not finite fails it. `NA` without a gradient, with
#'   `vb_correction = "full"`, and when the step was below `vb_tol` in every
#'   parameter, which is too small to check.
#' * `vb_solve`: 1 if that check failed, 0 if it passed, `NA` if it did not
#'   run. After a failed check the fit solves the VB equation on the same
#'   nodes, within the calls that the mode and Hessian took, and raises a
#'   warning of its own: the posterior is far from Gaussian along the mean
#'   correction, whether or not the solve converged.
#' * `mean_gap_max`: the largest absolute `mean_gap`.
#' * `ridge_bend_max`: the largest `ridge_bend`.
#' * `n_heavy_tails`: the number of parameters whose `tails` is `FALSE`.
#'
#' **Per-parameter diagnostics** (`type = "param"`), a data frame:
#' * `names`: the parameter.
#' * `grad`: the gradient of minus the log-posterior at the mode, as INLAvaan
#'   reports it: from `grad` when the fit has one, and otherwise from central
#'   differences of `logpost` a thousandth of a posterior SD apart. It
#'   should be close to 0.
#' * `grad_abs`: its absolute value.
#' * `grad_rel`: `grad_abs / (abs(par) + 1e-6)`, with `par` the VB mean (the
#'   mode when `vb_correction = "none"`), as in INLAvaan.
#' * `mode_shift_sigma`: the Newton step at the mode, in posterior SDs,
#'   \eqn{|(\Sigma g)_j| / s_j} with \eqn{g} the gradient. It should be
#'   close to 0.
#' * `kld`: the Kullback-Leibler divergence between the Laplace Gaussian
#'   marginal and its VB-shifted copy,
#'   \eqn{\delta_j^2 / (2\Sigma_{jj})}{delta_j^2 / (2 Sigma_jj)} =
#'   `vb_shift_sigma^2 / 2`. It is a function of the VB shift alone, and
#'   says nothing about the skew-normal fit.
#' * `vb_shift`: the VB shift of the mean, \eqn{\delta_j}, the VB mean minus
#'   the mode, on the scale of `logpost`.
#' * `vb_shift_sigma`: the same in posterior SDs.
#' * `nmad`: the NMAD of the skew-normal fit.
#' * `alpha`: the shape of the fitted skew-normal marginal.
#' * `unscanned_mass`: the mass the fitted skew-normal puts beyond the
#'   scanned window.
#'
#' and then inlavanilla's own:
#' * `grad_gap`: how far a supplied `grad` is from `logpost` at the mode: the
#'   larger of the difference in slope per conditional posterior SD (about how
#'   far it moves the mode, in SDs) and half the relative difference in
#'   curvature (about the relative error it makes in an SD), the curvature being
#'   that of the Hessian at the mode, the supplied `hess` when there is one. It
#'   comes from central differences of `logpost` at the mode, up to a posterior
#'   SD or so apart: without `hess` the fit takes them anyway, so it costs no
#'   calls, and with `hess` it takes them for this check alone, about 6m to 8m
#'   `logpost` calls (see "Gradients and Hessians" in [inla_vanilla()]). It
#'   should be close to 0, and is flagged above 0.1.
#' * `vb_step_sigma`: the step still to go at the VB mean, in posterior SDs:
#'   the fixed-point step there, the Newton step of the VB stationarity
#'   condition, as `mode_shift_sigma` is for the mode. After the fallback
#'   optimiser it is taken where the optimiser landed. On the fixed-point
#'   path it is the last step taken, which bounds the next, unless that step
#'   was 0.1 or more (a `vb_tol` above 0.1), when it is taken afresh. It
#'   should be close to 0 (`Inf` where the step is not finite). It is `NA`
#'   after a one-step correction, which has no step still to go. When the
#'   check of a one-step correction failed (`vb_solve`), it is the step still
#'   to go where the solve stopped, at its convergence or at its budget of
#'   calls. That value is recorded and not flagged, since the budget is not a
#'   setting.
#' * `mean_gap`, `gamma1`, `tails` and `ridge_bend`, which
#'   [summary.inla_vanilla()] describes.
#'
#' Some of INLAvaan's diagnostics are left out. `nsamp` counts posterior
#' draws, which [sampling()] makes on demand rather than during the fit.
#' `grad_num` and `grad_diff` compare an analytic gradient with finite
#' differences at the mode, which would take further calls; `grad_gap` makes
#' that comparison from differences the fit has already taken. `vb_mcse_max`,
#' `vb_mcse_mean` and `vb_mcse_sigma` estimate the quadrature error of
#' INLAvaan's Sobol points, and the cubature rule of inlavanilla is
#' deterministic (INLAvaan gives `NA` for its own deterministic rule).
#'
#' INLAvaan also exports a `diagnostics()` generic, an S4 one. When both
#' packages are attached, the one attached last masks the other. With
#' inlavanilla attached last, its `diagnostics()` passes an INLAvaan fit on
#' to INLAvaan's. With INLAvaan attached last, call
#' `inlavanilla::diagnostics(fit)` for an inlavanilla fit.
#'
#' @param object An object of class `inla_vanilla`.
#' @param type `"global"` (the default) for the diagnostics of the whole fit,
#'   or `"param"` for those of each parameter.
#' @param ... Unused.
#'
#' @returns For `type = "global"`, a named numeric vector of class
#'   `inla_vanilla_diagnostics`. For `type = "param"`, a data frame of class
#'   `inla_vanilla_diagnostics_param`. Their `print()` methods use
#'   INLAvaan's formats.
#'
#' @references
#' Jamil, H., & Rue, H. (2026). *Approximate Bayesian inference for
#' structural equation models using integrated nested Laplace
#' approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
#' <https://doi.org/10.48550/arXiv.2603.25690>
#'
#' @seealso [summary.inla_vanilla()] for what each check means and when it is
#'   flagged, and [plot.inla_vanilla()] with `type = "sn_fit"` for the
#'   skew-normal fits themselves.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(30, mean = 5, sd = 2)
#' logpost <- function(theta) {
#'   sum(dnorm(y, theta[["mu"]], exp(theta[["log_sigma"]]), log = TRUE)) +
#'     dnorm(theta[["mu"]], 0, 10, log = TRUE) +
#'     dnorm(theta[["log_sigma"]], 0, 2, log = TRUE)
#' }
#' fit <- inla_vanilla(logpost, c(mu = 0, log_sigma = 0), silent = TRUE)
#'
#' # The whole fit
#' diagnostics(fit)
#'
#' # One row per parameter
#' diagnostics(fit, type = "param")
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.default <- function(object, ...) {
  # INLAvaan's diagnostics() is an S4 generic: pass its fits on to it when it
  # is loaded, so that attaching inlavanilla after INLAvaan does not break it
  inlavaan <- if (isS4(object)) inlavaan_diagnostics()
  if (!is.null(inlavaan)) {
    return(inlavaan(object, ...))
  }
  abort(
    "There is no {.fn diagnostics} method for {.obj_type_friendly {object}}."
  )
}

# INLAvaan's diagnostics() generic when INLAvaan is loaded, or NULL.
inlavaan_diagnostics <- function() {
  if (isNamespaceLoaded("INLAvaan")) {
    getExportedValue("INLAvaan", "diagnostics")
  }
}

#' @rdname diagnostics
#' @export
diagnostics.inla_vanilla <- function(object, type = c("global", "param"), ...) {
  type <- check_choice(type, c("global", "param"), "type")
  dg <- diagnostics_internal(
    fit_optim(object),
    object$diagnostics$vb,
    object$diagnostics$marginals,
    pars = object$mean_vb
  )
  dg[[type]]
}

# The diagnostics of the mode search of a fit. A fit saved before 26.11.24
# has one condition number, the scale-free one, under the name `hess_cond`:
# it becomes `hess_cond_cor`, and `hess_cond` is the raw one of the
# covariance, as in newer fits.
fit_optim <- function(object) {
  optim <- object$diagnostics$optim
  if (!is.null(optim) && is.null(optim$hess_cond_cor)) {
    optim$hess_cond_cor <- optim$hess_cond
    ev <- eigen(object$cov, symmetric = TRUE, only.values = TRUE)$values
    optim$hess_cond <- if (min(ev) > 0) max(ev) / min(ev) else Inf
  }
  optim
}

# Printed as INLAvaan prints its diagnostics: counts as integers, the
# gradient norms and the quantities that span orders of magnitude in
# scientific notation, and the rest to four decimals.
#' @export
print.inla_vanilla_diagnostics <- function(x, ...) {
  nm <- names(x)
  int_names <- c(
    "npar",
    "converged",
    "iterations",
    "vb_applied",
    "vb_fallback",
    "vb_solve",
    "n_heavy_tails"
  )
  sci_names <- c(
    "hess_cond",
    "hess_cond_cor",
    "mode_shift_max",
    "unscanned_mass_max",
    "vb_step_max"
  )
  formatted <- vapply(
    seq_along(x),
    function(i) {
      val <- unclass(x)[[i]]
      name <- nm[i]
      if (is.na(val)) {
        "NA"
      } else if (name %in% int_names) {
        as.character(as.integer(round(val)))
      } else if (startsWith(name, "grad_") || name %in% sci_names) {
        formatC(val, digits = 2, format = "e")
      } else {
        formatC(val, digits = 4, format = "f", drop0trailing = FALSE)
      }
    },
    character(1)
  )
  names(formatted) <- nm
  print(formatted, quote = FALSE, right = TRUE)
  invisible(x)
}

# A healthy unscanned mass is below the last of four decimals, so it keeps
# three significant figures instead, as in INLAvaan.
#' @export
print.inla_vanilla_diagnostics_param <- function(x, digits = 4, ...) {
  out <- x
  for (nm in names(out)[vapply(out, is.numeric, logical(1))]) {
    out[[nm]] <- if (nm == "unscanned_mass") {
      signif(out[[nm]], 3)
    } else {
      round(out[[nm]], digits)
    }
  }
  print.data.frame(out, ...)
  invisible(x)
}
