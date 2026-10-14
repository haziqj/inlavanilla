# Finite differences. Every step `h` is per coordinate and recycled to the
# length of `x`, so callers can scale the steps by the posterior SDs once those
# are known.

# Central-difference gradient of a scalar function: 2m calls of `fn`.
fd_grad <- function(fn, x, h) {
  m <- length(x)
  h <- rep_len(h, m)
  g <- numeric(m)
  for (j in seq_len(m)) {
    xp <- xm <- x
    xp[j] <- x[j] + h[j]
    xm[j] <- x[j] - h[j]
    g[j] <- (fn(xp) - fn(xm)) / (2 * h[j])
  }
  g
}

# Central-difference Jacobian of a gradient function: 2m calls of `gr`.
fd_jacobian <- function(gr, x, h) {
  m <- length(x)
  h <- rep_len(h, m)
  H <- matrix(0, m, m)
  for (j in seq_len(m)) {
    xp <- xm <- x
    xp[j] <- x[j] + h[j]
    xm[j] <- x[j] - h[j]
    H[, j] <- (gr(xp) - gr(xm)) / (2 * h[j])
  }
  H
}

# Hessian of a scalar function from function values alone. The diagonal uses
# the usual three-point second difference. Each off-diagonal entry reuses the
# axial points and needs only two new calls,
#   f(x + hi ei + hj ej) + f(x - hi ei - hj ej) - f(x +- hi ei) - f(x +- hj ej)
#     + 2 f(x) = 2 hi hj H_ij + O(h^4),
# so the error is O(h^2), as for the four-point formula, at m^2 + m + 1 calls
# of `fn` rather than 2m^2 + 1.
fd_hessian <- function(fn, x, h, f0 = fn(x)) {
  m <- length(x)
  h <- rep_len(h, m)
  fp <- fm <- numeric(m)
  H <- matrix(0, m, m)
  for (i in seq_len(m)) {
    xp <- xm <- x
    xp[i] <- x[i] + h[i]
    xm[i] <- x[i] - h[i]
    fp[i] <- fn(xp)
    fm[i] <- fn(xm)
    H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h[i]^2
  }
  if (m > 1) {
    for (i in seq_len(m - 1)) {
      for (j in (i + 1):m) {
        xpp <- xmm <- x
        xpp[c(i, j)] <- x[c(i, j)] + h[c(i, j)]
        xmm[c(i, j)] <- x[c(i, j)] - h[c(i, j)]
        H[i, j] <- H[j, i] <- (fn(xpp) +
          fn(xmm) -
          fp[i] -
          fm[i] -
          fp[j] -
          fm[j] +
          2 * f0) /
          (2 * h[i] * h[j])
      }
    }
  }
  H
}

# Second derivative of a scalar function along direction `v` at `x`, from
# three calls (two when f(x) is supplied).
fd_curvature <- function(fn, x, v, h, f0 = fn(x)) {
  (fn(x + h * v) - 2 * f0 + fn(x - h * v)) / h^2
}

# Format a number of seconds for messages.
format_time <- function(s) {
  if (s < 1) {
    paste0(round(1000 * s), "ms")
  } else if (s < 60) {
    paste0(format(round(s, 1), nsmall = 1), "s")
  } else {
    paste0(floor(s / 60), "m ", round(s %% 60), "s")
  }
}

elapsed <- function() {
  proc.time()[["elapsed"]]
}

# The package's colour in the console: vanilla, the amber of a cured vanilla
# pod, between the bean and the orchid's throat of the logo
# (lab/logo/logo.R). It marks the rules and the values to read first. cli
# shows the nearest colour a console has, or none.
vanilla_hex <- "#B5772E"

col_vanilla <- function(x) {
  cli::make_ansi_style(vanilla_hex)(x)
}

# The width of everything the package prints: the console's, but never more
# than 80 columns. cli::console_width() follows options(width) in knitr, so a
# document can set a narrower width for its output.
out_width <- function() {
  min(cli::console_width(), 80L)
}

# Errors and warnings, formatted by cli but raised with base R. cli_abort()
# and cli_warn() go through rlang, which cli only suggests, so they could
# fail where rlang is not installed. Conditions carry the class
# "inlavanilla_error" or "inlavanilla_warning". An error raised while another
# is being handled keeps it as `parent`, and its message is appended.
abort <- function(
  message,
  class = NULL,
  parent = NULL,
  call = NULL,
  .envir = parent.frame()
) {
  old <- options(cli.condition_width = out_width())
  on.exit(options(old))
  msg <- cli::format_error(message, .envir = .envir)
  if (!is.null(parent)) {
    msg <- paste0(msg, "\nCaused by error: ", conditionMessage(parent))
  }
  stop(structure(
    class = c(class, "inlavanilla_error", "error", "condition"),
    list(message = msg, call = call, parent = parent)
  ))
}

# A warning is wrapped to out_width(), like the rest of the output. It also
# carries its message as plain text, with no colours and no line wrapping,
# for the record of warnings in a fit (see warning_text()).
warn <- function(message, class = NULL, .envir = parent.frame()) {
  force(.envir)
  old <- options(cli.condition_width = out_width())
  on.exit(options(old))
  warning(structure(
    class = c(class, "inlavanilla_warning", "warning", "condition"),
    list(
      message = cli::format_warning(message, .envir = .envir),
      call = NULL,
      text = format_plain(message, .envir)
    )
  ))
}

format_plain <- function(message, .envir) {
  old <- options(
    cli.condition_width = Inf,
    cli.width = Inf,
    cli.num_colors = 1
  )
  on.exit(options(old))
  cli::format_warning(message, .envir = .envir)
}

# The message of a warning as plain text: the package's own from warn(), any
# other with its colours stripped.
warning_text <- function(w) {
  if (inherits(w, "inlavanilla_warning") && is.character(w[["text"]])) {
    w[["text"]]
  } else {
    cli::ansi_strip(conditionMessage(w))
  }
}
