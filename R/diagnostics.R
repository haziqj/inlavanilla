# Fit-time diagnostics warning.
#
# One consolidated warning at the end of a fit when a marginal is clearly
# unreliable. The levels match INLAvaan's warn_fit_diagnostics() where the
# check is the same (NMAD, scan-end mass, VB shift), and are looser than the
# flags in summary(): a healthy fit stays silent, and a warning means that a
# marginal is clearly off. `tails` is inlavanilla's own check (see
# fit_marginals()). The mode and the user's gradient have their own warnings.
warn_levels <- list(nmad = 0.1, scan_end_mass = 0.05, vb_shift = 1)

# `diag` is fit$diagnostics$marginals; `vbd` is fit$diagnostics$vb (NULL when
# vb = FALSE). Returns the bullets of the warning, or none, as cli templates
# with the environment that fills them in attribute "envir". The templates
# are filled in only when the warning is raised, so every value they show
# needs a variable of its own.
diagnostic_problems <- function(diag, vbd) {
  rn <- rownames(diag)
  out <- character()
  heavy <- rn[!is.na(diag$tails) & !diag$tails]
  if (length(heavy) > 0L) {
    out <- c(
      out,
      x = paste0(
        "Heavy tails for {.field {heavy}}: the density is still over 5% of ",
        "its peak at the ends of the scan, 4 Laplace SDs out, and the ",
        "skew-normal extrapolates beyond them."
      )
    )
  }
  nmad <- diag$nmad
  bad <- !is.na(nmad) & nmad > warn_levels$nmad
  if (any(bad)) {
    worst_nmad <- signif(max(nmad[bad]), 2)
    poor <- rn[bad]
    out <- c(
      out,
      x = paste0(
        "Poor skew-normal fit for {.field {poor}}: largest NMAD {worst_nmad}, ",
        "above ",
        warn_levels$nmad,
        "."
      )
    )
  }
  sem <- diag$scan_end_mass
  bad <- !is.na(sem) & sem > warn_levels$scan_end_mass
  if (any(bad)) {
    worst_sem <- signif(max(sem[bad]), 2)
    beyond <- rn[bad]
    out <- c(
      out,
      x = paste0(
        "The skew-normal for {.field {beyond}} puts up to {worst_sem} of its ",
        "mass beyond the scanned window (above ",
        warn_levels$scan_end_mass,
        "), so its outer quantiles are extrapolated."
      )
    )
  }
  if (!is.null(vbd)) {
    shift <- abs(diag$vb_shift)
    if (isTRUE(vbd$fallback)) {
      out <- c(
        out,
        x = paste0(
          "The VB fixed-point iteration failed and the fallback optimiser ",
          "ran: the posterior is far from the Laplace Gaussian."
        )
      )
    } else if (any(shift > warn_levels$vb_shift)) {
      worst_shift <- signif(max(shift), 2)
      moved <- rn[shift > warn_levels$vb_shift]
      out <- c(
        out,
        x = paste0(
          "The VB correction moved the mean of {.field {moved}} by up to ",
          "{worst_shift} posterior SDs from the mode: the posterior is far from ",
          "the Laplace Gaussian."
        )
      )
    }
  }
  attr(out, "envir") <- environment()
  out
}

warn_diagnostics <- function(diag, vbd) {
  problems <- diagnostic_problems(diag, vbd)
  if (length(problems) == 0L) {
    return(invisible(FALSE))
  }
  warn(
    c(
      "The approximation may be poor for some parameters.",
      problems,
      "i" = paste0(
        "See {.code summary()} of the fit for every diagnostic, and ",
        "{.code vignette(\"diagnostics\", package = \"inlavanilla\")} for ",
        "what they mean."
      )
    ),
    class = "inlavanilla_diagnostics_warning",
    .envir = attr(problems, "envir")
  )
  invisible(TRUE)
}
