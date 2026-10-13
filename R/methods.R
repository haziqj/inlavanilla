#' @export
print.inla_vanilla <- function(x, digits = 3, ...) {
  m <- nrow(x$sn)
  cli::cat_rule(
    left = cli::style_bold("INLAvanilla"),
    right = cli::col_grey(format_time(x$timing[["total"]])),
    line_col = "cyan"
  )
  cli::cat_line(
    cli::col_grey(
      m,
      if (m == 1L) " parameter" else " parameters",
      sep_dot(),
      if (is_gaussian(x)) "Gaussian marginals" else "skew-normal marginals",
      if (isTRUE(x$vb)) paste0(sep_dot(), "VB-corrected") else ""
    )
  )
  cli::cat_line()
  cat_table(x$summary, digits = digits, scale = x$summary[, "SD"])
  cli::cat_line()
  cli::cat_line(
    cli::style_bold("Log marginal likelihood: "),
    cli::col_cyan(format(round(x$logml, 2), nsmall = 2))
  )
  invisible(x)
}

#' Summarise an INLAvanilla fit
#'
#' Posterior summaries from the skew-normal marginals, on the unconstrained
#' scale of `logpost`, together with diagnostics of each stage of the
#' approximation.
#'
#' The diagnostics, and the levels at which they are flagged, are:
#' * **Mode.** The largest gradient at the mode in whitened units: how far, in
#'   posterior SDs, a Newton step would still move it. The mode is accepted
#'   when this is below 0.01, whatever the return code of [stats::nlminb()],
#'   whose message is shown otherwise.
#' * **VB correction.** The shift of each mean from the mode, in posterior
#'   SDs (`vb_shift`). It is flagged when it exceeds 1 SD, or when the
#'   fixed-point iteration failed and the fallback optimiser ran: both mean
#'   the posterior is far from the Laplace Gaussian.
#' * **Tails.** Each parameter is scanned over \eqn{\pm 4} Laplace SDs.
#'   `tails` is `FALSE` (and flagged) when the density at an end of the scan
#'   is still over 5% of its peak: the Laplace SD badly understates a heavy
#'   tail there, and the skew-normal extrapolates it.
#' * **Scan path.** `bend` is the log-density, in nats, that the straight
#'   scan misses at \eqn{\pm 2} posterior SDs: the gain from one Newton step
#'   towards the conditional mode of the other parameters. It is flagged above
#'   0.1. A large value means the other parameters' conditional mode bends away
#'   from the line, as a poorly determined noise scale does for a regression
#'   coefficient at small n, and the marginal comes out too narrow, by roughly
#'   0.35 to 0.5 times `bend` in relative SD (see the Accuracy section of
#'   [inla_vanilla()]).
#' * **Skew-normal fits.** `nmad`, the largest gap between the fitted and the
#'   scanned density relative to its peak (flagged above 0.05), and
#'   `sn_vs_vb`, the gap between the
#'   skew-normal mean and the VB-corrected mean in posterior SDs. Both means
#'   target the posterior mean, so they should agree closely (to
#'   \eqn{O(n^{-2})}; Proposition 3.4 of Jamil and Rue, 2026a). A gap above
#'   0.1 SDs is flagged: the posterior is far from Gaussian along that
#'   parameter, and either approximation may be off.
#'
#' * **Scan-end mass.** `scan_end_mass`, the mass each fitted skew-normal
#'   puts outside the scanned window (the mode \eqn{\pm 4} Laplace SDs),
#'   where the marginal is extrapolated, as INLAvaan's `diagnostics()`
#'   reports it. A Gaussian marginal gives \eqn{2\Phi(-4) = 6.3 \times
#'   10^{-5}}, and well-fitted marginals usually sit between \eqn{10^{-3}}
#'   and \eqn{10^{-2}}. It is flagged above 0.05.
#'
#' Some problems also raise one warning at the end of the fit, of class
#' `inlavanilla_diagnostics_warning`: heavy tails, an NMAD above 0.1, a
#' scan-end mass above 0.05, a VB shift above 1 SD or the VB fallback. These
#' levels match INLAvaan's fit-time warning; the flags here are more
#' sensitive. With `silent = TRUE` the warning is only recorded in
#' `fit$warnings`. `vignette("diagnostics", package = "inlavanilla")`
#' explains every check with examples.
#'
#' With `vb = FALSE` there is no VB mean, and the VB columns are left out.
#' With `marginals = "gaussian"` there are no scans, so only the VB shift is
#' reported.
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
#' Jamil, H. and Rue, H. (2026a). Approximate Bayesian inference for
#' structural equation models using integrated nested Laplace approximations.
#' \doi{10.48550/arXiv.2603.25690}
#' @export
summary.inla_vanilla <- function(
  object,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
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
    left = cli::style_bold("INLAvanilla summary"),
    line_col = "cyan"
  )
  cli::cat_line()
  cli::cat_line(cli::style_bold("Posterior marginals"))
  cat_table(x$table, digits = digits, scale = x$table[, "SD"])
  cli::cat_line()

  cli::cat_line(cli::style_bold("Diagnostics"))
  diag <- x$diagnostics
  gaussian <- is_gaussian(x)
  cols <- c("vb_shift", "sn_vs_vb", "alpha", "nmad", "bend")
  heads <- c("VB shift", "SN - VB", "Shape", "NMAD", "Bend")
  keep <- (!is.null(x$vb) | c(FALSE, FALSE, TRUE, TRUE, TRUE)) &
    (!gaussian | c(TRUE, FALSE, FALSE, FALSE, FALSE))
  if (any(keep)) {
    dtab <- as.matrix(diag[, cols[keep], drop = FALSE])
    colnames(dtab) <- heads[keep]
    cat_table(dtab, digits = 2, by_cell = TRUE)
    cli::cat_line()
  }

  opt <- x$optim
  if (opt$grad_mode < 1e-2) {
    cat_status(
      "ok",
      "Mode found in {opt$iterations} iteration{?s} (gradient {fmt_num(opt$grad_mode)} SDs)."
    )
  } else {
    cat_status(
      "warn",
      "Mode search stopped with a gradient of {fmt_num(opt$grad_mode)} SDs ({opt$message})."
    )
  }

  if (!is.null(x$vb)) {
    max_shift <- max(abs(diag$vb_shift))
    how <- if (x$vb$fallback) " (via the fallback optimiser)" else ""
    cat_status(
      if (x$vb$fallback || max_shift > 1) "warn" else "ok",
      "VB correction in {x$vb$iterations} iteration{?s}{how}: max shift {fmt_num(max_shift)} SDs."
    )
  } else {
    cat_status("info", "No VB correction (vb = FALSE).")
  }

  if (gaussian) {
    cat_status(
      "info",
      "Gaussian marginals: no scans, so no skew-normal checks."
    )
  }
  sem <- diag$scan_end_mass
  if (!gaussian && !is.null(sem) && any(sem > 0.05, na.rm = TRUE)) {
    beyond <- rownames(diag)[!is.na(sem) & sem > 0.05]
    cat_status(
      "warn",
      "The skew-normal for {.field {beyond}} puts up to {fmt_num(max(sem, na.rm = TRUE))} of its mass beyond the scanned window: its outer quantiles are extrapolated."
    )
  }
  if (!gaussian && !all(diag$tails)) {
    cat_status(
      "warn",
      "Heavy tails for {.field {rownames(diag)[!diag$tails]}}: the density is still over 5% of its peak 4 Laplace SDs out, and the skew-normal extrapolates beyond."
    )
  }
  bent <- !gaussian & diag$bend > 0.1
  if (any(bent)) {
    cat_status(
      "warn",
      "The straight scan misses up to {fmt_num(max(diag$bend))} nats at +-2 SDs for {.field {rownames(diag)[bent]}}: the marginal may be too narrow."
    )
  }
  max_nmad <- max(diag$nmad)
  max_gap <- max(abs(diag$sn_vs_vb))
  if (gaussian) {
    # nothing was scanned or fitted
  } else if (max_nmad > 0.05) {
    cat_status(
      "warn",
      "Skew-normal fit is poor for {.field {rownames(diag)[which.max(diag$nmad)]}} (NMAD {fmt_num(max_nmad)})."
    )
  } else {
    cat_status("ok", "Skew-normal fits: max NMAD {fmt_num(max_nmad)}.")
  }
  if (!is.null(x$vb) && !gaussian) {
    if (max_gap > 0.1) {
      cat_status(
        "warn",
        "Skew-normal and VB means differ by up to {fmt_num(max_gap)} SDs, for {.field {rownames(diag)[which.max(abs(diag$sn_vs_vb))]}}."
      )
    } else {
      cat_status(
        "ok",
        "Skew-normal and VB means agree to within {fmt_num(max_gap)} SDs."
      )
    }
  }
  cli::cat_line()

  vb_term <- x$logml - x$logml_laplace
  cli::cat_line(
    cli::style_bold("Log marginal likelihood: "),
    cli::col_cyan(format(round(x$logml, 3), nsmall = 3)),
    if (!is.null(x$vb)) {
      cli::col_grey(
        " (Laplace ",
        format(round(x$logml_laplace, 3), nsmall = 3),
        ", VB term ",
        format(signif(vb_term, 2)),
        ")"
      )
    }
  )
  calls <- paste0(x$counts[["logpost"]], " logpost")
  if (x$has_grad) {
    calls <- paste0(calls, " and ", x$counts[["grad"]], " grad")
  }
  tm <- x$timing
  cli::cat_line(
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
  invisible(x)
}

#' Plot the posterior marginals of an INLAvanilla fit
#'
#' Draws the marginal density of each parameter, on the unconstrained scale,
#' with the Gaussian marginal of the Laplace approximation for comparison. The
#' vertical line marks the posterior mean. For a fit with
#' `marginals = "gaussian"` the marginals are Gaussian: centred at the VB mean,
#' or, with `vb = FALSE`, the Laplace marginals themselves, drawn without a
#' second copy.
#'
#' @param x An object of class `inla_vanilla`.
#' @param pars Parameters to plot, by name or position. Defaults to all.
#' @param laplace Logical. Overlay the Laplace (Gaussian) marginals? Default
#'   `TRUE`; ignored when the marginals are the Laplace ones.
#' @param ... Further graphical parameters passed to [graphics::plot()] for
#'   each panel, such as `xlab` or `main`; they replace the defaults (the
#'   parameter name as the title, no axis labels, no y axis).
#'
#' @returns `x`, invisibly.
#' @export
plot.inla_vanilla <- function(x, pars = NULL, laplace = TRUE, ...) {
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
  k <- length(idx)
  # Setting `mfrow` resets `cex`, and so does restoring it: keep the caller's
  old_cex <- graphics::par("cex")
  # Wider than tall: panels side by side, then rows
  op <- graphics::par(
    mfrow = rev(grDevices::n2mfrow(k)),
    mar = c(2.5, 1, 2, 1),
    oma = c(0, 0, if (laplace) 1.5 else 0, 0),
    mgp = c(2, 0.6, 0),
    las = 1
  )
  on.exit({
    graphics::par(op)
    graphics::par(cex = old_cex)
  })
  mom <- sn_moments(sn[, "xi"], sn[, "omega"], sn[, "alpha"])
  for (j in idx) {
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
  if (laplace) {
    # One legend for all panels, in the outer margin at the top
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    graphics::par(new = TRUE)
    graphics::plot.new()
    graphics::legend(
      "top",
      legend = c(label, "Laplace"),
      col = c(col_sn, col_lap),
      lty = c(1, 2),
      lwd = c(2, 1.5),
      horiz = TRUE,
      bty = "n",
      cex = 0.85,
      xpd = NA
    )
  }
  invisible(x)
}

# Print a numeric matrix as an aligned table: bold column headers, cyan row
# names. With `scale` (one value per row, such as the posterior SD), each row
# is printed to `digits` significant figures of its scale, so parameters on
# different scales each get a sensible number of decimals. With `by_cell`,
# each number is formatted on its own. Widths are worked out on the plain
# strings before any styling.
cat_table <- function(tab, digits = 3, scale = NULL, by_cell = FALSE) {
  cells <- matrix("", nrow(tab), ncol(tab))
  for (i in seq_len(nrow(tab))) {
    cells[i, ] <- if (by_cell) {
      formatC(tab[i, ], digits = digits, format = "g")
    } else {
      s <- abs(scale[i])
      dec <- if (is.finite(s) && s > 0) digits - 1 - floor(log10(s)) else digits
      formatC(tab[i, ], digits = min(max(dec, 0), 8), format = "f")
    }
  }
  heads <- colnames(tab)
  widths <- pmax(nchar(heads), apply(nchar(cells), 2, max))
  rn <- rownames(tab)
  rw <- max(nchar(rn), 0L)
  pad <- function(s, w) {
    vapply(seq_along(s), function(i) formatC(s[i], width = w[i]), "")
  }
  header <- paste0(
    strrep(" ", rw),
    paste0("  ", cli::style_bold(pad(heads, widths)), collapse = "")
  )
  cli::cat_line(header)
  for (i in seq_len(nrow(cells))) {
    cli::cat_line(
      cli::col_cyan(formatC(rn[i], width = -rw)),
      paste0("  ", pad(cells[i, ], widths), collapse = "")
    )
  }
}

cat_status <- function(type, msg, .envir = parent.frame()) {
  sym <- switch(
    type,
    ok = cli::col_green(cli::symbol$tick),
    warn = cli::col_yellow("!"),
    fail = cli::col_red(cli::symbol$cross),
    info = cli::col_cyan(cli::symbol$info)
  )
  cli::cat_line(sym, " ", cli::format_inline(msg, .envir = .envir))
}

# Whether a fit has Gaussian marginals (marginals = "gaussian").
is_gaussian <- function(x) {
  identical(x$marginals, "gaussian")
}

# Separator and the posterior-SD unit, in plain ASCII where the output cannot
# show UTF-8.
sep_dot <- function() {
  if (cli::is_utf8_output()) " \u00b7 " else ", "
}

sd_unit <- function() {
  if (cli::is_utf8_output()) "\u03c3" else " SD"
}

fmt_num <- function(x) {
  format(signif(x, 2), scientific = abs(x) < 1e-3 && x != 0)
}
