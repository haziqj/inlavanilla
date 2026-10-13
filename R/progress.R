# Progress reporting, at three levels: `silent = TRUE` shows nothing, the
# default shows one line per stage, and `verbose = TRUE` adds indented grey
# detail lines under each one. Each stage gets a cli spinner or progress bar,
# which appears only when the stage runs for longer than the user's
# `cli.progress_show_after` seconds, or half a second when that option is
# unset (cli's own default is 2), and ends with one coloured line. The lines
# are built directly from cli's styles rather than through its glue-based
# formatters, which cost a few milliseconds a call: as long as a whole fit of
# a small model. The detail lines are only built when they are shown.

# `show`: print anything at all; `detail`: add the detail lines. The warnings
# raised during the fit are recorded in `warnings`.
progress_new <- function(show, detail = FALSE, frame = parent.frame()) {
  p <- new.env(parent = emptyenv())
  p$show <- show
  p$detail <- show && detail
  p$frame <- frame
  p$id <- NULL
  p$calls <- FALSE # whether the log-posterior wrapper should tick the bar
  p$warnings <- character()
  p
}

# Start a stage. With `total`, a bar over that many items; without, a spinner
# that counts log-posterior calls. The delay is set while the bar is created
# (cli reads it then) and the option is restored afterwards. A user's own
# setting is kept, so a knitr or Quarto setup chunk can keep bars, which
# cannot redraw in place there, out of the document.
progress_start <- function(p, msg, total = NA) {
  if (!p$show) {
    return(invisible(p))
  }
  old <- options(
    cli.progress_show_after = getOption("cli.progress_show_after", 0.5)
  )
  on.exit(options(old))
  fmt <- if (is.na(total)) {
    paste0(
      "{cli::pb_spin} ",
      msg,
      " {cli::col_grey(paste(cli::pb_current, 'calls'))}"
    )
  } else {
    paste0(
      "{cli::pb_spin} ",
      msg,
      " {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
      " {cli::col_grey(cli::pb_eta_str)}"
    )
  }
  p$id <- cli::cli_progress_bar(
    format = fmt,
    total = total,
    clear = TRUE,
    .envir = p$frame
  )
  p$calls <- is.na(total)
  invisible(p)
}

progress_tick <- function(p, set = NULL) {
  if (!is.null(p$id)) {
    if (is.null(set)) {
      cli::cli_progress_update(id = p$id, .envir = p$frame)
    } else {
      cli::cli_progress_update(id = p$id, set = set, .envir = p$frame)
    }
  }
}

# End a stage with a tick (or a yellow ! when `ok` is FALSE), `text`, and the
# time it took.
progress_done <- function(p, text, secs, ok = TRUE) {
  if (!p$show) {
    return(invisible(p))
  }
  if (!is.null(p$id)) {
    cli::cli_progress_done(id = p$id, .envir = p$frame)
  }
  p$id <- NULL
  p$calls <- FALSE
  message(
    if (ok) cli::col_green(cli::symbol$tick) else cli::col_yellow("!"),
    " ",
    text,
    " ",
    cli::col_grey(paste0("[", format_time(secs), "]"))
  )
  invisible(p)
}

# Highlight a value in a progress line.
hl <- function(x) {
  cli::col_cyan(x)
}

# Verbose detail ---------------------------------------------------------------
#
# The report_*() functions are called only when p$detail is TRUE, and use
# only what the fit has already computed, with no further calls of logpost or
# grad.

# A line of its own, marked with a cyan i.
progress_info <- function(p, text) {
  if (p$detail) {
    message(cli::col_cyan(cli::symbol$info), " ", text)
  }
}

# Detail lines, already styled (grey, with any flags in yellow), indented
# under the stage line, one message each.
progress_detail <- function(p, lines) {
  if (p$detail) {
    for (line in lines) {
      message("  ", line)
    }
  }
}

grey <- function(...) {
  cli::col_grey(paste0(...))
}

flag <- function() {
  cli::col_yellow("!")
}

with_sd <- function(x) {
  paste0(fmt_num(x), sd_unit())
}

plural <- function(n, word) {
  paste0(n, " ", word, if (n == 1) "" else "s")
}

# "12 logpost calls", "3 logpost and 40 grad calls"
format_calls <- function(n_lp, n_gr) {
  parts <- c(
    if (n_lp > 0 || n_gr == 0) paste(n_lp, "logpost"),
    if (n_gr > 0) paste(n_gr, "grad")
  )
  paste0(
    paste(parts, collapse = " and "),
    if (n_lp + n_gr == 1) " call" else " calls"
  )
}

# "a", "a and b", "a, b and c", "a, b, c, d, e and 25 more"
format_names <- function(x, max = 6L) {
  n <- length(x)
  if (n > max) {
    x <- c(x[seq_len(max - 1L)], paste(n - max + 1L, "more"))
    n <- max
  }
  if (n == 1L) x else paste(paste(x[-n], collapse = ", "), "and", x[n])
}

format_value <- function(x) {
  format(round(x, 2), nsmall = 2)
}

# Before the mode search: the parameters, where the derivatives come from,
# logpost at init, and for a user gradient its check at init.
report_start <- function(p, init, lp_init, grad_check, calls) {
  m <- length(init)
  progress_info(
    p,
    paste0(plural(m, "parameter"), ": ", format_names(names(init)))
  )
  lines <- if (is.null(grad_check)) {
    grey("Derivatives by finite differences of logpost")
  } else {
    paste0(
      grey(
        "Gradient from grad, checked at init: relative discrepancy ",
        fmt_num(grad_check)
      ),
      if (grad_check > 1e-3) flag() else "",
      grey(
        " (warning above 0.001), ",
        format_calls(calls[[1]], calls[[2]])
      )
    )
  }
  lines <- c(lines, grey("logpost at init ", format_value(lp_init[[1L]])))
  progress_detail(p, lines)
}

# The mode and the Laplace approximation.
report_mode <- function(p, lap, calls) {
  opt <- lap$opt
  ev <- opt$evaluations
  lines <- grey(
    "nlminb: ",
    opt$message,
    ", after ",
    plural(opt$iterations, "iteration"),
    " (",
    ev[["function"]],
    " function and ",
    ev[["gradient"]],
    " gradient evaluations)"
  )
  rs <- lap$restart
  if (!is.null(rs)) {
    first <- if (is.finite(rs$grad_mode)) {
      paste0("whitened gradient ", with_sd(rs$grad_mode))
    } else {
      "Hessian not positive definite"
    }
    lines <- c(
      lines,
      grey(
        "Restarted nlminb with each parameter scaled by its curvature (first ",
        "run: ",
        plural(rs$iterations, "iteration"),
        ", ",
        first,
        "); ",
        if (rs$kept) "kept the restart" else "kept the first run"
      )
    )
  }
  # The whitened gradient decides the restart, and the Newton step (INLAvaan's
  # mode shift, in posterior SDs) the warning
  msm <- max(lap$mode_shift)
  lines <- c(
    lines,
    paste0(
      grey(
        "logpost at the mode ",
        format_value(lap$lp_max),
        ", whitened gradient ",
        with_sd(lap$grad_mode),
        ", Newton step ",
        with_sd(msm)
      ),
      if (any(off_level(lap$mode_shift, warn_levels$mode_shift))) {
        flag()
      } else {
        ""
      },
      grey(" (warning above ", with_sd(warn_levels$mode_shift), ")")
    )
  )
  hs <- lap$hessian
  lines <- c(
    lines,
    if (hs$by == "grad") {
      c(
        grey(
          "Hessian by central differences of grad: ",
          plural(hs$calls, "call"),
          ", and ",
          hs$step_calls,
          " logpost calls to set the steps"
        ),
        paste0(
          grey("Jacobian asymmetry ", fmt_num(lap$asymmetry)),
          if (lap$asymmetry > 1e-3) flag() else "",
          grey(" (warning above 0.001)")
        )
      )
    } else {
      grey(
        "Hessian by second differences of logpost: ",
        plural(hs$calls, "call"),
        ", and ",
        hs$step_calls,
        " more to set the steps"
      )
    }
  )
  # The smallest eigenvalue of the correlation matrix is near 0 when some
  # combination of the parameters is poorly identified; the ratio of the
  # largest to it is the condition number of the Hessian in posterior-SD
  # units, which the fit-time warning checks
  if (length(lap$mode) > 1L) {
    rho <- stats::cov2cor(lap$Sigma)
    upper <- which(upper.tri(rho), arr.ind = TRUE)
    k <- which.max(abs(rho[upper]))
    pair <- rownames(rho)[upper[k, ]]
    ev_min <- min(eigen(rho, symmetric = TRUE, only.values = TRUE)$values)
    lines <- c(
      lines,
      grey(
        "Posterior correlation matrix: largest correlation ",
        fmt_num(rho[upper][k]),
        " (",
        pair[1],
        ", ",
        pair[2],
        "), smallest eigenvalue ",
        fmt_num(ev_min),
        if (ev_min < 0.01) {
          ": some combination of the parameters is poorly identified"
        } else {
          ""
        }
      ),
      paste0(
        grey("Hessian condition number ", fmt_num(lap$hess_cond)),
        if (off_level(lap$hess_cond, warn_levels$hess_cond)) flag() else "",
        grey(
          ", in posterior-SD units (warning above ",
          fmt_num(warn_levels$hess_cond),
          ")"
        )
      )
    )
  }
  lines <- c(
    lines,
    grey("Laplace log marginal likelihood ", format_value(lap$logml)),
    grey(format_calls(calls[[1]], calls[[2]]))
  )
  progress_detail(p, lines)
}

# The VB correction: the fixed-point steps, the fallback, the term it adds to
# the log marginal likelihood and the largest shift.
report_vb <- function(p, vbc, sds, logml_laplace, calls) {
  h <- vbc$history
  steps <- vapply(h$steps, fmt_num, "")
  n <- length(steps)
  if (n == 0L) {
    steps <- "none"
  } else if (n > 8L) {
    steps <- c(steps[1:3], "...", steps[(n - 2L):n])
  }
  lines <- grey(
    "Fixed-point steps, largest over the parameters (",
    trimws(sd_unit()),
    "): ",
    paste(steps, collapse = ", ")
  )
  n_acc <- sum(h$accelerated)
  lines <- c(
    lines,
    if (vbc$fallback) {
      paste0(
        flag(),
        grey(
          " Fell back to nlminb, as ",
          h$reason,
          ": ",
          h$nlminb,
          ", after ",
          plural(vbc$iterations, "iteration")
        )
      )
    } else {
      grey(
        "Converged in ",
        plural(vbc$iterations, "iteration"),
        if (n_acc > 0) {
          paste0(", Anderson acceleration at ", n_acc, " of them")
        } else {
          ", without Anderson acceleration"
        },
        "; no fallback"
      )
    }
  )
  shift_sd <- vbc$shift / sds
  j <- which.max(abs(shift_sd))
  lines <- c(
    lines,
    grey(
      "Largest shift ",
      with_sd(shift_sd[[j]]),
      ", for ",
      names(vbc$shift)[j]
    ),
    grey(
      "KL term ",
      fmt_num(vbc$kld_global),
      ", giving a log marginal likelihood of ",
      format_value(logml_laplace + vbc$kld_global)
    ),
    grey(format_calls(calls[[1]], calls[[2]]))
  )
  progress_detail(p, lines)
}

# The skew-normal marginals: one aligned line per parameter, with the values
# that summary() flags marked, at its levels.
report_marginals <- function(p, sn, calls) {
  rn <- rownames(sn)
  cols <- list(
    Shape = formatC(sn[, "alpha"], digits = 2, format = "f"),
    NMAD = formatC(sn[, "nmad"], digits = 2, format = "g"),
    Bend = formatC(sn[, "bend"], digits = 2, format = "g"),
    Tails = ifelse(sn[, "tails"] == 1, "ok", "heavy"),
    `End mass` = formatC(sn[, "scan_end_mass"], digits = 2, format = "g"),
    gamma1 = formatC(sn[, "gamma1"], digits = 2, format = "g")
  )
  none <- logical(length(rn))
  flags <- list(
    Shape = none,
    NMAD = off_level(sn[, "nmad"], warn_levels$nmad),
    Bend = off_level(sn[, "bend"], summary_levels$bend),
    Tails = !sn[, "tails"] %in% 1,
    `End mass` = off_level(sn[, "scan_end_mass"], warn_levels$scan_end_mass),
    gamma1 = none
  )
  rw <- max(nchar(rn))
  widths <- vapply(
    names(cols),
    function(k) max(nchar(k), nchar(cols[[k]])),
    numeric(1)
  )
  # Right-aligned, each value followed by its flag or a space; the last
  # column is never flagged
  header <- grey(
    formatC("", width = rw),
    paste0("  ", mapply(formatC, names(cols), width = widths), collapse = "")
  )
  rows <- vapply(
    seq_along(rn),
    function(i) {
      cells <- vapply(
        names(cols),
        function(k) {
          paste0(
            grey(formatC(cols[[k]][i], width = widths[[k]])),
            if (isTRUE(flags[[k]][i])) flag() else " "
          )
        },
        ""
      )
      sub(
        " $",
        "",
        paste0(
          grey(formatC(rn[i], width = -rw)),
          " ",
          paste0(" ", cells, collapse = "")
        )
      )
    },
    ""
  )
  progress_detail(
    p,
    c(
      header,
      rows,
      paste0(
        flag(),
        grey(
          " marks what summary() flags: NMAD above ",
          warn_levels$nmad,
          ", bend above ",
          summary_levels$bend,
          ", heavy tails, scan-end mass above ",
          warn_levels$scan_end_mass
        )
      ),
      grey(format_calls(calls[[1]], calls[[2]]))
    )
  )
}

# The end of the fit: the calls and time in all, and the warnings.
report_end <- function(p, n_lp, n_gr, secs) {
  n_warn <- length(p$warnings)
  progress_info(
    p,
    paste0(
      "Done: ",
      format_calls(n_lp, n_gr),
      " in ",
      format_time(secs),
      if (n_warn > 0) {
        paste0(", with ", cli::col_yellow(plural(n_warn, "warning")))
      } else {
        ""
      }
    )
  )
}
