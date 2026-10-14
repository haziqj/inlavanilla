#' @export
print.inla_vanilla <- function(x, digits = 3, ...) {
  cli::cat_rule(
    width = out_width(),
    left = cli::style_bold("INLAvanilla"),
    right = cli::col_grey(method_label(is_gaussian(x), isTRUE(x$vb))),
    line_col = vanilla_hex
  )
  cli::cat_line()
  cat_table(x$summary, digits = digits, scale = x$summary[, "SD"])
  cli::cat_line()
  cli::cat_line(
    cli::style_bold("Log marginal likelihood: "),
    col_vanilla(format(round(x$logml, 2), nsmall = 2))
  )
  # Converged means the mode check of summary() passes, whatever the return
  # code of nlminb(): the Newton step at the mode is within its level
  dg <- diagnostics_internal(
    x$diagnostics$optim,
    x$diagnostics$vb,
    x$diagnostics$marginals
  )
  flags <- check_flags(dg, is_gaussian(x))
  took <- paste0(
    plural(x$diagnostics$optim$iterations, "iteration"),
    " (",
    format_time(x$timing[["total"]]),
    ")"
  )
  if (flags[["mode"]]) {
    cat_status("warn", "Not converged after {took}.")
  } else {
    cli::cat_line("Converged in ", took, ".")
  }
  # Whether to trust the table: shown only when something is off
  n_warnings <- length(x$warnings)
  if (n_warnings > 0L) {
    cat_status(
      "warn",
      "{n_warnings} warning{?s} during the fit, kept in {.code $warnings}."
    )
  }
  n_flagged <- sum(flags)
  if (n_flagged > 0L) {
    cat_status("warn", "{n_flagged} check{?s} flagged by {.fn summary}.")
  }
  invisible(x)
}

#' Summarise an INLAvanilla fit
#'
#' Posterior summaries from the skew-normal marginals, on the unconstrained
#' scale of `logpost`, together with diagnostics of each stage of the
#' approximation.
#'
#' Printed, the summary gives the fit first (the log marginal likelihood,
#' and the calls and time each stage took), then the posterior marginals,
#' then the diagnostics.
#'
#' Most of the checks are those of INLAvaan's `diagnostics()`, and
#' [diagnostics()] returns them under the same names. Jamil and Rue (2026) use
#' the VB shift and NMAD as fit diagnostics (Sections 3.1 and 4.4), and the
#' gap between the skew-normal and VB means as an internal check (Section
#' 5.1). The mode shift, the Hessian condition number and the unscanned mass
#' come from INLAvaan's `diagnostics()`, not from the paper. `vb_step_sigma`,
#' `tails` and `ridge_bend` are inlavanilla's own.
#'
#' A posterior SD here is the Laplace SD at the mode,
#' \eqn{\sqrt{\Sigma_{jj}}}{sqrt(Sigma_jj)}, as in INLAvaan's
#' `diagnostics()`; it is not the SD column of the table, which is the
#' skew-normal's. The diagnostics, with the level at which each is flagged
#' where it is, are:
#' * **Mode.** `mode_shift_sigma`, the Newton step at the mode in posterior
#'   SDs: how far one step of Newton's method would still move each
#'   parameter. It is flagged above 0.1, whatever the return code of
#'   [stats::nlminb()], whose message is then shown: without a gradient
#'   `nlminb()` can report failure at an accurate mode, and success short of
#'   one. With the VB check below, it is what makes the fit converged.
#' * **Hessian.** `hess_cond`, the condition number of the Hessian at the
#'   mode with each parameter in posterior-SD units (the eigenvalue ratio of
#'   the posterior correlation matrix). It is flagged above \eqn{10^8}: some
#'   combination of the parameters is then weakly identified. Its line is
#'   shown when there is more than one parameter.
#' * **VB correction.** `vb_shift_sigma`, the shift of each mean from the
#'   mode in posterior SDs (`vb_shift` is the same shift on the scale of
#'   `logpost`). The VB correction is the shift of the Laplace Gaussian that
#'   minimises its Kullback-Leibler divergence from the posterior, and a
#'   large shift is a sign of a posterior far from that Gaussian (Jamil &
#'   Rue, 2026, Section 3.1). It is flagged above 1 SD. `kld`, the
#'   Kullback-Leibler divergence between each Laplace Gaussian marginal and
#'   its VB-shifted copy, is `vb_shift_sigma^2 / 2`, as in INLAvaan; it is
#'   not flagged. When the fixed-point iteration hands over to the fallback
#'   optimiser, the VB line says so, but, as in INLAvaan, that is not
#'   flagged: it is a route of the solver, not a sign of a poor
#'   approximation. Where the solver lands is checked, as the mode is:
#'   `vb_step_sigma` is the step still to go at the VB mean, in posterior
#'   SDs (the fixed-point step there, which is the Newton step of the VB
#'   stationarity condition). It is flagged above 0.1, the level of the mode
#'   check, as a correction that did not settle. That check is inlavanilla's
#'   own.
#' * **Skew-normal fits.** `nmad`, the largest gap between the fitted and
#'   the scanned density, relative to its peak (Jamil & Rue, 2026, Section
#'   4.4), flagged above 0.1. `alpha` is the shape of the fitted skew-normal
#'   marginal, on the unconstrained scale: zero is a Gaussian marginal, and
#'   the sign gives the direction of the skew. It describes the marginal
#'   rather than checking it, so it is not in the printed table: it is in
#'   `summary(fit)$diagnostics` and [diagnostics()], and `fit$sn` holds it
#'   with the other skew-normal parameters.
#' * **Volume correction.** `gamma1`, the slope added to the scanned
#'   log-heights, per posterior SD along the scan, for the changing volume of
#'   the other parameters (Jamil & Rue, 2026, Corollary 3.3). It is in
#'   `summary(fit)$diagnostics` and in the verbose output of the fit, not in
#'   the printed table, and it is not flagged.
#' * **Unscanned mass.** `unscanned_mass`, the mass each fitted skew-normal
#'   puts outside the scanned window (4 posterior SDs either side of the
#'   mode), where the marginal is extrapolated, as INLAvaan's `diagnostics()`
#'   reports it. A Gaussian marginal gives \eqn{2\Phi(-4) = 6.3 \times
#'   10^{-5}}; the larger it is, the more of the marginal is extrapolated.
#'   It is flagged above 0.05.
#' * **Tails.** `tails` is `FALSE` (and flagged) when the density at an end
#'   of the scan is still over 5% of its peak: the Laplace SD badly
#'   understates a heavy tail there, and the skew-normal extrapolates it.
#' * **Ridge bend.** Each parameter is scanned along the conditional mode
#'   path (CMP) of the Laplace Gaussian (Jamil & Rue, 2026, Lemma 3.1), a
#'   straight line that is a first-order approximation to a possibly curved
#'   ridge of high posterior density (Section 3.2.1). `ridge_bend` is the
#'   log-density, in nats, that the straight scan misses at \eqn{\pm 2}
#'   posterior SDs: the gain from one Newton step towards the conditional
#'   mode of the other parameters. A large value means the ridge bends away
#'   from the CMP, as it does for a regression coefficient when the noise
#'   scale is poorly determined at small n, and the marginal comes out too
#'   narrow. For a normal mean with an unknown SD, a ridge bend \eqn{b} makes
#'   the marginal SD too small by
#'   \eqn{1 - \sqrt{1 - 3b/4}}{1 - sqrt(1 - 3b/4)}, about \eqn{3b/8}. In
#'   other models this is a rough guide: on test models it was mostly within
#'   a factor of 1.5 of the true shortfall, above it in logistic and Poisson
#'   regressions and below it in a nonlinear regression. `ridge_bend` is
#'   flagged above 0.25, a shortfall of 10% by that formula (see the Accuracy
#'   section of [inla_vanilla()]).
#' * **Mean gap.** `mean_gap`, the gap between the skew-normal mean and the
#'   VB-corrected mean, in posterior SDs. Both recover the leading,
#'   \eqn{O(n^{-1})} shift from the mode to the posterior mean, so they
#'   agree to \eqn{O(n^{-2})} (Jamil & Rue, 2026, Proposition 3.4). A gap
#'   signals that the higher-order terms of that expansion about the mode
#'   matter, or numerical error (Section 5.1). A gap above 0.1 SDs is
#'   flagged: the posterior is far from Gaussian along that parameter, and
#'   either approximation may be off.
#'
#' Every check flagged here but the ridge bend and the mean gap also raises
#' one warning at the end of the fit, of class
#' `inlavanilla_diagnostics_warning`, at the same level. The levels of the
#' mode shift, NMAD, unscanned mass, VB shift and Hessian condition number
#' are those of INLAvaan's fit-time warning; the VB and heavy-tails checks
#' are inlavanilla's own, as are the ridge bend and mean gap flags. Every
#' level is a package default, chosen so that a healthy fit stays silent, and
#' not a calibrated cut-off: a tripped check is a prompt to look again, and
#' silence is not a certificate of accuracy. In the simulations of Jamil and
#' Rue (2026, Section 4.4) the VB shift and NMAD track the loss of accuracy,
#' and the paper reports a calibration of the expected error against them in
#' its Supplementary Material.
#'
#' The table gives each value to two significant figures, as a tibble
#' prints numbers, with the flagged values in red where the console shows
#' colour. The status lines list every flagged parameter, worst first, with
#' its value (the heavy-tails line, a yes/no check, lists names only, in
#' parameter order); the warning lists at most three, as INLAvaan's does. With
#' `silent = TRUE` the warning is only recorded in `fit$warnings`. Printing
#' the fit ends by saying whether it converged, which is the mode check
#' here, with the iterations of [stats::nlminb()] and the time the whole fit
#' took. A line for its warnings and one counting the checks flagged here
#' follow, each shown only when there are any.
#' `vignette("diagnostics", package = "inlavanilla")` explains every check
#' with examples.
#'
#' With `vb = FALSE` there is no VB mean: the VB columns of
#' `summary(fit)$diagnostics` are `NA`, as INLAvaan's `diagnostics()` leaves
#' them, and are left out of the table. With `marginals = "gaussian"` there
#' are no scans, so the scan columns are `NA` and only the VB shift is shown.
#'
#' @param object An object of class `inla_vanilla`.
#' @param probs Probabilities for the reported quantiles.
#' @param ... Unused.
#'
#' @returns An object of class `summary.inla_vanilla`: a list with the summary
#'   table `table`, the log marginal likelihood `logml` (and its Laplace part
#'   `logml_laplace`), the per-parameter diagnostics `diagnostics`, and
#'   details of the optimiser, the VB correction, the number of calls and the
#'   timing.
#'
#' @references
#' Jamil, H., & Rue, H. (2026). *Approximate Bayesian inference for
#' structural equation models using integrated nested Laplace
#' approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
#' <https://doi.org/10.48550/arXiv.2603.25690>
#'
#' @seealso [diagnostics()] for the diagnostics as numbers, under INLAvaan's
#'   names, and [plot.inla_vanilla()] with `type = "sn_fit"` for the
#'   skew-normal fits.
#' @export
summary.inla_vanilla <- function(
  object,
  probs = c(0.025, 0.5, 0.975),
  ...
) {
  if (!is.numeric(probs) || any(!is.finite(probs) | probs <= 0 | probs >= 1)) {
    abort(
      "{.arg probs} must be probabilities strictly between 0 and 1."
    )
  }
  structure(
    list(
      table = sn_summary(object$sn, probs = probs),
      logml = object$logml,
      logml_laplace = object$logml_laplace,
      diagnostics = object$diagnostics$marginals,
      optim = object$diagnostics$optim,
      vb = object$diagnostics$vb,
      counts = object$counts,
      timing = object$timing,
      marginals = object$marginals,
      has_grad = object$has_grad
    ),
    class = "summary.inla_vanilla"
  )
}

#' @export
print.summary.inla_vanilla <- function(x, digits = 3, ...) {
  cli::cat_rule(
    width = out_width(),
    left = cli::style_bold("INLAvanilla summary"),
    right = cli::col_grey(method_label(is_gaussian(x), !is.null(x$vb))),
    line_col = vanilla_hex
  )
  cli::cat_line()
  # The fit, then its parameters, then how far to trust them
  vb_term <- x$logml - x$logml_laplace
  cat_wrap(paste0(
    cli::style_bold("Log marginal likelihood: "),
    col_vanilla(format(round(x$logml, 3), nsmall = 3)),
    if (!is.null(x$vb)) {
      cli::col_grey(
        " (Laplace ",
        format(round(x$logml_laplace, 3), nsmall = 3),
        ", VB term ",
        format(signif(vb_term, 2)),
        ")"
      )
    }
  ))
  calls <- paste0(x$counts[["logpost"]], " logpost")
  if (x$has_grad) {
    calls <- paste0(calls, " and ", x$counts[["grad"]], " grad")
  }
  tm <- x$timing
  cat_wrap(
    cli::col_grey(
      calls,
      " calls in ",
      format_time(tm[["total"]]),
      " (mode ",
      format_time(tm[["laplace"]]),
      ", VB ",
      format_time(tm[["vb"]]),
      ", marginals ",
      format_time(tm[["marginals"]]),
      ")"
    )
  )
  cli::cat_line()

  cat_heading("Posterior marginals")
  cat_table(x$table, digits = digits, scale = x$table[, "SD"])
  cli::cat_line()

  cat_heading("Diagnostics")
  diag <- x$diagnostics
  gaussian <- is_gaussian(x)
  has_vb <- !is.null(x$vb)
  opt <- x$optim
  dg <- diagnostics_internal(opt, x$vb, diag)
  glob <- dg$global
  # VB first, then the skew-normal checks in the order of their status
  # lines, ending with the comparison of the two means
  cols <- c(
    "vb_shift_sigma",
    "nmad",
    "unscanned_mass",
    "ridge_bend",
    "mean_gap"
  )
  heads <- c("VB shift", "NMAD", "Unscanned", "Ridge bend", "Mean gap")
  keep <- (has_vb | c(FALSE, TRUE, TRUE, TRUE, FALSE)) &
    (!gaussian | c(TRUE, FALSE, FALSE, FALSE, FALSE))
  if (any(keep)) {
    dtab <- as.matrix(diag[, cols[keep], drop = FALSE])
    colnames(dtab) <- heads[keep]
    flag <- param_flags(dg, gaussian)[, cols[keep], drop = FALSE]
    cat_table(dtab, digits = 2, flag = flag)
    cli::cat_line()
  }

  # A status line per stage, and one "!" line for each flagged check, so that
  # the count in print() matches. Flagged lines list every flagged parameter,
  # worst first, through variables in `env` (see list_template()).
  flags <- check_flags(dg, gaussian)
  env <- new.env(parent = environment())
  flagged_line <- function(check) {
    cat_status("warn", problem_text(check, dg, env, n = Inf), .envir = env)
  }

  if (!flags[["mode"]]) {
    msm <- with_sd(glob[["mode_shift_max"]])
    cat_status(
      "ok",
      "Mode found in {opt$iterations} iteration{?s} (Newton step {msm})."
    )
  } else {
    listed <- flagged_list(
      env,
      "mode",
      dg,
      "mode_shift_sigma",
      warn_levels$mode_shift,
      Inf
    )
    cat_status(
      "warn",
      paste0(
        "Mode search stopped short ({.val {opt$message}}): Newton step > ",
        with_sd(warn_levels$mode_shift),
        " for ",
        listed,
        "."
      ),
      .envir = env
    )
  }

  if (flags[["hess_cond"]]) {
    flagged_line("hess_cond")
  } else if (glob[["npar"]] > 1) {
    hc <- glob[["hess_cond"]]
    cat_status("ok", "Hessian condition number {fmt_num(hc)}, in SD units.")
  }

  if (has_vb) {
    how <- if (isTRUE(x$vb$fallback)) " (via the fallback optimiser)" else ""
    # Whether it settled, the VB half of converged, then the size of the
    # shift: a line for each that is flagged, or one line when neither is
    if (flags[["vb_step"]]) {
      listed <- flagged_list(
        env,
        "vbstep",
        dg,
        "vb_step_sigma",
        warn_levels$vb_step,
        Inf
      )
      cat_status(
        "warn",
        paste0(
          "VB correction did not settle in {x$vb$iterations} ",
          "iteration{?s}{how}: step > ",
          with_sd(warn_levels$vb_step),
          " for ",
          listed,
          "."
        ),
        .envir = env
      )
    }
    if (!flags[["vb"]] && !flags[["vb_step"]]) {
      max_shift <- with_sd(glob[["vb_shift_max"]])
      cat_status(
        "ok",
        "VB correction settled in {x$vb$iterations} iteration{?s}{how}: max shift {max_shift}."
      )
    } else if (flags[["vb"]]) {
      listed <- flagged_list(
        env,
        "vbs",
        dg,
        "vb_shift_sigma",
        warn_levels$vb_shift,
        Inf
      )
      cat_status(
        "warn",
        paste0(
          "VB correction in {x$vb$iterations} iteration{?s}{how}: shift > ",
          with_sd(warn_levels$vb_shift),
          " for ",
          listed,
          "."
        ),
        .envir = env
      )
    }
  } else {
    cat_status("info", "No VB correction (vb = FALSE).")
  }

  if (gaussian) {
    cat_status(
      "info",
      "Gaussian marginals: no scans, so no skew-normal checks."
    )
  } else {
    # One line for the checks that pass, with the largest value of each as
    # the table shows it, then a line for each that does not
    fmt_cell <- function(v) trimws(format_sigfig(v, 2)$plain)
    healthy <- c(
      if (!flags[["nmad"]]) {
        paste("NMAD", fmt_cell(glob[["nmad_max"]]))
      },
      if (!flags[["unscanned_mass"]]) {
        paste("unscanned mass", fmt_cell(glob[["unscanned_mass_max"]]))
      },
      if (!flags[["ridge_bend"]]) {
        paste("ridge bend", fmt_cell(glob[["ridge_bend_max"]]))
      },
      if (has_vb && !flags[["mean_gap"]]) {
        paste0("mean gap ", fmt_cell(glob[["mean_gap_max"]]), sd_unit())
      }
    )
    if (length(healthy) > 0L) {
      sn_ok <- paste(healthy, collapse = ", ")
      cat_status("ok", "Skew-normal fits: max {sn_ok}.")
    }
    for (check in c(
      "nmad",
      "unscanned_mass",
      "tails",
      "ridge_bend",
      "mean_gap"
    )) {
      if (flags[[check]]) {
        flagged_line(check)
      }
    }
  }
  invisible(x)
}

#' Plot an INLAvanilla fit
#'
#' By default (`type = "marg_pdf"`), draws the marginal density of each
#' parameter, on the unconstrained scale, with the Gaussian marginal of the
#' Laplace approximation for comparison. The vertical line marks the
#' posterior mean. For a fit with `marginals = "gaussian"` the marginals are
#' Gaussian: centred at the VB mean, or, with `vb = FALSE`, the Laplace
#' marginals themselves, drawn without a second copy.
#'
#' `type = "sn_fit"` shows each skew-normal fit against the scan it was
#' fitted to, as INLAvaan's `plot(type = "sn_fit")` does. The x axis is the
#' scan coordinate, in posterior SDs from the mode (the Laplace SDs, from -4
#' to 4). "Raw" is the log-posterior along the scan and "Corrected" the same
#' with the volume correction `gamma1` added, each as a density scaled to a
#' peak of 1, and "Skew-normal fit" is the fitted skew-normal at the height
#' it was fitted at. The largest gap between the last two, as a fraction of
#' the peak, is the NMAD of [summary.inla_vanilla()]. `type = "sn_fit_log"`
#' draws the same on a base-10 log scale, which shows the tails. Both use the
#' scans kept in the fit, with no further calls of `logpost`, and need
#' skew-normal marginals.
#'
#' @param x An object of class `inla_vanilla`.
#' @param pars Parameters to plot, by name or position. Defaults to all.
#' @param type The plot: `"marg_pdf"` (the default) for the marginal
#'   densities, `"sn_fit"` for each skew-normal fit against its scan, or
#'   `"sn_fit_log"` for the same on a log scale, as INLAvaan names them.
#' @param laplace Logical. Overlay the Laplace (Gaussian) marginals? Default
#'   `TRUE`; ignored when the marginals are the Laplace ones, and by the
#'   `"sn_fit"` plots.
#' @param ... Further graphical parameters passed to [graphics::plot()] for
#'   each panel, such as `xlab` or `main`; they replace the defaults (the
#'   parameter name as the title, no axis labels, and for `"marg_pdf"` no y
#'   axis).
#'
#' @returns `x`, invisibly.
#' @export
plot.inla_vanilla <- function(
  x,
  pars = NULL,
  type = c("marg_pdf", "sn_fit", "sn_fit_log"),
  laplace = TRUE,
  ...
) {
  type <- check_choice(type, c("marg_pdf", "sn_fit", "sn_fit_log"), "type")
  check_flag(laplace, "laplace")
  sn <- x$sn
  idx <- if (is.null(pars)) seq_len(nrow(sn)) else pars
  if (is.character(idx)) {
    bad <- setdiff(idx, rownames(sn))
    if (length(bad) > 0L) {
      abort("Unknown parameter{?s}: {.val {bad}}.")
    }
    idx <- match(idx, rownames(sn))
  }
  if (
    !is.numeric(idx) ||
      length(idx) == 0L ||
      anyNA(idx) ||
      any(idx < 1 | idx > nrow(sn))
  ) {
    abort("{.arg pars} must name or index parameters of the fit.")
  }
  if (type == "marg_pdf") {
    plot_marginals(x, idx, laplace, ...)
  } else {
    if (is.null(x$scan)) {
      abort(c(
        "There are no scans to plot: the marginals are Gaussian
         ({.code marginals = \"gaussian\"}).",
        "i" = "Use {.code type = \"marg_pdf\"}, or fit with
               {.code marginals = \"skewnormal\"}."
      ))
    }
    plot_scans(x, idx, log = type == "sn_fit_log", ...)
  }
  invisible(x)
}

# Panels side by side, then in rows, one per parameter in `idx`, each drawn
# by `draw(j)`, with one legend for all of them in the outer margin at the
# top when `legend` (a list of arguments to legend()) is given.
plot_panels <- function(idx, draw, legend = NULL, mar = c(2.5, 1, 2, 1)) {
  # Setting `mfrow` resets `cex`, and so does restoring it: keep the caller's
  old_cex <- graphics::par("cex")
  # Wider than tall: panels side by side, then rows
  op <- graphics::par(
    mfrow = rev(grDevices::n2mfrow(length(idx))),
    mar = mar,
    oma = c(0, 0, if (is.null(legend)) 0 else 1.5, 0),
    mgp = c(2, 0.6, 0),
    las = 1
  )
  on.exit({
    graphics::par(op)
    graphics::par(cex = old_cex)
  })
  for (j in idx) {
    draw(j)
  }
  if (!is.null(legend)) {
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    graphics::par(new = TRUE)
    graphics::plot.new()
    do.call(
      graphics::legend,
      c(list("top", horiz = TRUE, bty = "n", cex = 0.85, xpd = NA), legend)
    )
  }
}

# The marginal densities, with the Laplace marginals overlaid.
plot_marginals <- function(x, idx, laplace, ...) {
  sn <- x$sn
  col_sn <- "#C2410C"
  col_lap <- "#0F766E"
  label <- if (!is_gaussian(x)) {
    "Skew-normal"
  } else if (isTRUE(x$vb)) {
    "VB Gaussian"
  } else {
    # The marginals are the Laplace ones: nothing to overlay
    laplace <- FALSE
    "Laplace"
  }
  mom <- sn_moments(sn[, "xi"], sn[, "omega"], sn[, "alpha"])
  draw <- function(j) {
    lap_sd <- sqrt(x$cov[j, j])
    q_sn <- qsnorm(
      c(0.0005, 0.9995),
      sn[j, "xi"],
      sn[j, "omega"],
      sn[j, "alpha"]
    )
    q_lap <- stats::qnorm(c(0.0005, 0.9995), x$mode[j], lap_sd)
    lims <- if (laplace) range(q_sn, q_lap) else q_sn
    tt <- seq(lims[1], lims[2], length.out = 400)
    f_sn <- dsnorm(tt, sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"])
    f_lap <- stats::dnorm(tt, x$mode[j], lap_sd)
    args <- utils::modifyList(
      list(
        x = tt,
        y = f_sn,
        type = "n",
        ylim = c(0, max(f_sn, if (laplace) f_lap) * 1.05),
        main = rownames(sn)[j],
        font.main = 1,
        xlab = "",
        ylab = "",
        yaxt = "n",
        bty = "n"
      ),
      list(...)
    )
    do.call(graphics::plot, args)
    graphics::polygon(
      c(tt, rev(tt)),
      c(f_sn, numeric(length(tt))),
      col = grDevices::adjustcolor(col_sn, 0.15),
      border = NA
    )
    graphics::lines(tt, f_sn, lwd = 2, col = col_sn)
    if (laplace) {
      graphics::lines(tt, f_lap, lty = 2, lwd = 1.5, col = col_lap)
    }
    graphics::segments(
      mom$mean[j],
      0,
      mom$mean[j],
      dsnorm(mom$mean[j], sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"]),
      col = col_sn,
      lty = 3
    )
  }
  legend <- if (laplace) {
    list(
      legend = c(label, "Laplace"),
      col = c(col_sn, col_lap),
      lty = c(1, 2),
      lwd = c(2, 1.5)
    )
  }
  plot_panels(idx, draw, legend)
}

# Each skew-normal fit against its scan, as INLAvaan's visual_debug() draws
# it: the raw and the corrected scan, each scaled to a peak of 1, and the
# fitted skew-normal on the scan coordinate, at the height it was fitted at
# (its logC). Scan points where the log-posterior is not finite were left
# out of the fit, and are left out here. On the log scale, values that are
# not positive are dropped.
plot_scans <- function(x, idx, log = FALSE, ...) {
  z <- x$scan$z
  sds <- sqrt(diag(x$cov))
  gamma1 <- x$diagnostics$marginals$gamma1
  labels <- c("Raw", "Corrected", "Skew-normal fit")
  cols <- c("gray60", "black", "red2")
  ltys <- c(1, 1, 2)
  scale_y <- function(v) {
    if (log) {
      v[v <= 0] <- NA
      log10(v)
    } else {
      v
    }
  }
  draw <- function(j) {
    raw <- x$scan$logpost[, j]
    ok <- is.finite(raw)
    zs <- z[ok]
    raw <- raw[ok]
    y <- raw + gamma1[j] * zs
    zz <- seq(min(z), max(z), length.out = 201)
    fitted <- dsnorm(
      zz,
      (x$sn[j, "xi"] - x$mode[[j]]) / sds[[j]],
      x$sn[j, "omega"] / sds[[j]],
      x$sn[j, "alpha"],
      logC = x$scan$logC[[j]]
    )
    curves <- list(
      scale_y(exp(raw - max(raw))),
      scale_y(exp(y - max(y))),
      scale_y(fitted)
    )
    ylim <- range(unlist(curves), na.rm = TRUE)
    if (!log) {
      ylim <- c(0, ylim[2] * 1.05)
    }
    args <- utils::modifyList(
      list(
        x = range(z),
        y = ylim,
        type = "n",
        main = rownames(x$sn)[j],
        font.main = 1,
        xlab = "",
        ylab = "",
        bty = "n"
      ),
      list(...)
    )
    do.call(graphics::plot, args)
    xs <- list(zs, zs, zz)
    for (k in seq_along(curves)) {
      graphics::lines(
        xs[[k]],
        curves[[k]],
        col = cols[k],
        lty = ltys[k],
        lwd = 1.3
      )
    }
  }
  plot_panels(
    idx,
    draw,
    legend = list(legend = labels, col = cols, lty = ltys, lwd = 1.3),
    mar = c(2.5, 2.5, 2, 1)
  )
}

# Print a numeric matrix as an aligned table: bold column headers, plain row
# names. With `scale` (one value per row, such as the posterior SD), each row
# is printed to `digits` significant figures of its scale, so parameters on
# different scales each get a sensible number of decimals. Without it, each
# column is formatted as a tibble prints it (see format_sigfig()), to
# `digits` significant figures. The cells where `flag` (a logical matrix
# like `tab`) is TRUE are red. Widths are worked out on the plain strings
# before any styling.
cat_table <- function(tab, digits = 3, scale = NULL, flag = NULL) {
  cells <- styled <- matrix("", nrow(tab), ncol(tab))
  if (is.null(scale)) {
    for (j in seq_len(ncol(tab))) {
      col <- format_sigfig(tab[, j], digits)
      cells[, j] <- col$plain
      styled[, j] <- col$styled
    }
  } else {
    for (i in seq_len(nrow(tab))) {
      s <- abs(scale[i])
      dec <- if (is.finite(s) && s > 0) digits - 1 - floor(log10(s)) else digits
      cells[i, ] <- formatC(
        tab[i, ],
        digits = min(max(dec, 0), 8),
        format = "f"
      )
    }
    styled <- cells
  }
  if (!is.null(flag)) {
    styled[flag] <- cli::col_red(cells[flag])
  }
  heads <- colnames(tab)
  widths <- pmax(nchar(heads), apply(nchar(cells), 2, max))
  rn <- rownames(tab)
  rw <- max(nchar(rn), 0L)
  pad <- function(s, w) {
    vapply(seq_along(s), function(i) formatC(s[i], width = w[i]), "")
  }
  # Right-align a styled cell by the width of its plain text
  pad_cells <- function(i, cols) {
    w <- widths[cols]
    paste0(strrep(" ", w - nchar(cells[i, cols])), styled[i, cols])
  }
  # Split the columns into blocks that fit the console width, as print() does
  # for a wide matrix
  avail <- out_width() - rw
  block <- integer(length(widths))
  b <- 1L
  used <- 0
  for (j in seq_along(widths)) {
    if (used > 0 && used + widths[j] + 2 > avail) {
      b <- b + 1L
      used <- 0
    }
    block[j] <- b
    used <- used + widths[j] + 2
  }
  for (b in unique(block)) {
    cols <- block == b
    if (b > 1) {
      cli::cat_line()
    }
    cli::cat_line(
      strrep(" ", rw),
      paste0(
        "  ",
        cli::style_bold(pad(heads[cols], widths[cols])),
        collapse = ""
      )
    )
    for (i in seq_len(nrow(cells))) {
      cli::cat_line(
        formatC(rn[i], width = -rw),
        paste0("  ", pad_cells(i, cols), collapse = "")
      )
    }
  }
}

# Format a numeric column as the pillar package prints a tibble's: `sigfig`
# significant figures for each number, every digit before the decimal point
# kept, and a decimal digit dropped only where the number has no more (so 1.5
# and 0.33 but 0.30 for 0.3001). A number rounded to a whole one keeps its
# point ("123."). The decimal points line up. A column wider than 13
# characters is printed in scientific notation instead. Returns the plain
# cells, padded to one width, and the same cells styled: the leading zeros
# of a number below 1 and the digits past the significant ones are grey.
format_sigfig <- function(x, sigfig = 3) {
  n <- length(x)
  neg <- ifelse(!is.na(x) & x < 0, "-", "")
  lhs <- rhs <- lhs_grey <- character(n)
  for (k in seq_len(n)) {
    ax <- abs(x[k])
    if (!is.finite(ax)) {
      neg[k] <- ""
      lhs[k] <- format(x[k])
      next
    }
    if (ax == 0) {
      lhs[k] <- "0"
      next
    }
    expo <- floor(log10(ax))
    # Decimals for `sigfig` figures, fewer where the number has no more
    rhs_digits <- max(sigfig - 1 - expo, 0)
    while (rhs_digits > 0) {
      v <- ax * 10^(rhs_digits - 1)
      if (abs(v - round(v)) > sqrt(.Machine$double.eps) * max(1, v)) {
        break
      }
      rhs_digits <- rhs_digits - 1
    }
    r <- signif(ax, max(sigfig, expo + 1))
    s <- formatC(r, format = "f", digits = rhs_digits)
    if (rhs_digits == 0 && r != ax) {
      s <- paste0(s, ".")
    }
    parts <- regmatches(s, regexpr("\\.", s), invert = TRUE)[[1]]
    lhs[k] <- parts[1]
    rhs[k] <- if (length(parts) > 1) paste0(".", parts[2]) else ""
  }
  plain_num <- paste0(neg, lhs, rhs)
  if (max(nchar(plain_num)) > 13) {
    expo <- ifelse(is.finite(x) & x != 0, floor(log10(abs(x))), 0)
    mant <- trimws(formatC(x / 10^expo, format = "fg", digits = sigfig))
    plain_num <- ifelse(is.finite(x), paste0(mant, "e", expo), format(x))
    plain <- formatC(plain_num, width = max(nchar(plain_num)))
    return(list(plain = plain, styled = plain))
  }
  # Line up the decimal points: integer parts right-aligned, decimals left
  int <- paste0(neg, lhs)
  int_w <- max(nchar(int))
  dec_w <- max(nchar(rhs))
  pad_l <- strrep(" ", int_w - nchar(int))
  pad_r <- strrep(" ", dec_w - nchar(rhs))
  plain <- paste0(pad_l, int, rhs, pad_r)
  # Grey: "0.00" before the first significant digit of a number below 1, and
  # the integer digits past `sigfig`
  grey <- cli::col_grey
  styled <- character(n)
  for (k in seq_len(n)) {
    if (lhs[k] == "0" && nzchar(rhs[k])) {
      lead <- regmatches(rhs[k], regexpr("^\\.0*", rhs[k]))
      body <- paste0(
        neg[k],
        grey(paste0("0", lead)),
        substring(rhs[k], nchar(lead) + 1)
      )
    } else if (grepl("^[0-9]+$", lhs[k]) && nchar(lhs[k]) > sigfig) {
      body <- paste0(
        neg[k],
        substring(lhs[k], 1, sigfig),
        grey(substring(lhs[k], sigfig + 1)),
        rhs[k]
      )
    } else {
      body <- paste0(neg[k], lhs[k], rhs[k])
    }
    styled[k] <- paste0(pad_l[k], body, pad_r[k])
  }
  list(plain = plain, styled = styled)
}

# A section heading, as cli's level-two headings look: "── Title ──", and
# a blank line under it.
cat_heading <- function(title) {
  line <- col_vanilla(strrep(cli::symbol$line, 2))
  cli::cat_line(line, " ", cli::style_bold(title), " ", line)
  cli::cat_line()
}

cat_status <- function(type, msg, .envir = parent.frame()) {
  sym <- switch(
    type,
    ok = cli::col_green(cli::symbol$tick),
    warn = cli::col_yellow("!"),
    fail = cli::col_red(cli::symbol$cross),
    info = col_vanilla(cli::symbol$info)
  )
  cat_wrap(paste(sym, cli::format_inline(msg, .envir = .envir)), exdent = 2)
}

# Write a line wrapped to out_width(), with continuation lines indented by
# `exdent` spaces.
cat_wrap <- function(text, exdent = 0) {
  cli::cat_line(cli::ansi_strwrap(text, width = out_width(), exdent = exdent))
}

# Whether a fit has Gaussian marginals (marginals = "gaussian").
is_gaussian <- function(x) {
  identical(x$marginals, "gaussian")
}

# The approximation, as the headers of print() and summary() name it: the
# marginals, and "+ VB" when the VB correction ran. Gaussian marginals are
# those of the Laplace approximation (centred at the VB mean with VB).
method_label <- function(gaussian, vb) {
  paste0(if (gaussian) "Laplace" else "Skew-normal", if (vb) " + VB")
}

# The checks that summary() flags, each TRUE when flagged: the mode, the
# Hessian condition number, the VB shift, and, for skew-normal marginals,
# NMAD, the unscanned mass, the tails, the ridge bend and the mean gap.
# `dg` is from diagnostics_internal(). The levels are those of the fit-time
# warning (warn_levels), which raises every one of these but the ridge bend
# and the mean gap (summary_levels). A check that did not run is not
# flagged; a value that is not finite is. print() counts the flags, and
# summary() prints one "!" line for each and colours the flagged cells of
# its table, all from param_flags().
check_flags <- function(dg, gaussian) {
  pf <- param_flags(dg, gaussian)
  c(
    mode = any(pf[, "mode_shift_sigma"]),
    vb_step = any(pf[, "vb_step_sigma"]),
    hess_cond = off_level(dg$global[["hess_cond"]], warn_levels$hess_cond),
    vb = any(pf[, "vb_shift_sigma"]),
    nmad = any(pf[, "nmad"]),
    unscanned_mass = any(pf[, "unscanned_mass"]),
    tails = any(pf[, "tails"]),
    ridge_bend = any(pf[, "ridge_bend"]),
    mean_gap = any(pf[, "mean_gap"])
  )
}

# The per-parameter flags behind check_flags(): a logical matrix with a row
# per parameter and a column per check, named as in dg$param.
param_flags <- function(dg, gaussian) {
  p <- dg$param
  sn <- !gaussian
  vb <- isTRUE(dg$global[["vb_applied"]] == 1)
  cbind(
    mode_shift_sigma = off_level(p$mode_shift_sigma, warn_levels$mode_shift),
    vb_step_sigma = vb & off_level(p$vb_step_sigma, warn_levels$vb_step),
    vb_shift_sigma = vb &
      off_level(abs(p$vb_shift_sigma), warn_levels$vb_shift),
    nmad = sn & off_level(p$nmad, warn_levels$nmad),
    unscanned_mass = sn &
      off_level(p$unscanned_mass, warn_levels$unscanned_mass),
    tails = sn & !p$tails %in% TRUE,
    ridge_bend = sn & off_level(p$ridge_bend, summary_levels$ridge_bend),
    mean_gap = sn & vb & off_level(abs(p$mean_gap), summary_levels$mean_gap)
  )
}

# The posterior-SD unit, in plain ASCII where the output cannot show UTF-8.
sd_unit <- function() {
  if (cli::is_utf8_output()) "\u03c3" else " SD"
}

fmt_num <- function(x) {
  if (!is.finite(x)) {
    return(format(x))
  }
  format(signif(x, 2), scientific = x != 0 && (abs(x) < 1e-3 || abs(x) >= 1e5))
}
