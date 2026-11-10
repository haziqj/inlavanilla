# Natural-scale reporting. The engine works on the unconstrained scale of
# `logpost` and knows nothing of this file. A parameter whose name starts with
# a prefix of `transform_dict` ("log_sigma") is reported on its natural scale
# ("sigma") by print(), summary(), coef(), plot() and sampling(), each a map of
# the fitted skew-normal marginal through a one-to-one, increasing function.

# The transforms recognised in parameter names, keyed by the prefix without
# its underscore. `forward` takes the natural scale to the unconstrained one
# and `inverse` takes it back, `log_jac` is the log of the derivative of
# `inverse`, and `call` is how the note under a table writes the inverse.
# These three cover the internal scales of nearly every hyperparameter in
# R-INLA: the log (precisions, shapes, ranges), the logit (probabilities) and
# the correlation map 2 plogis(x) - 1 = tanh(x / 2), Fisher's z doubled.
# Every inverse must be increasing, as natural_summary() maps the quantiles
# straight through it.
transform_dict <- list(
  log = list(
    forward = log,
    inverse = exp,
    log_jac = function(theta) theta,
    call = "exp"
  ),
  logit = list(
    forward = stats::qlogis,
    inverse = stats::plogis,
    log_jac = function(theta) {
      stats::plogis(theta, log.p = TRUE) + stats::plogis(-theta, log.p = TRUE)
    },
    call = "plogis"
  ),
  atanh = list(
    forward = atanh,
    inverse = tanh,
    # log(1 - tanh(theta)^2), written so that it does not cancel in the tails
    log_jac = function(theta) {
      a <- abs(theta)
      log(4) - 2 * a - 2 * log1p(exp(-2 * a))
    },
    call = "tanh"
  )
)

# The transform of each parameter, read from its name: a data frame with a row
# per parameter, holding its name on the unconstrained scale (`parameter`),
# the name it is reported under (`name`) and the transform (`transform`,
# "identity" for none). A prefix counts only when it is followed by an
# underscore and something more, so "log_sigma" is a log and "logistic",
# "log_" and "Log_sigma" are not. A name that would clash once its prefix is
# gone, as "log_sigma" does beside a parameter called "sigma", is left alone,
# with a warning. With `detect = FALSE` nothing is transformed.
detect_transforms <- function(parnames, detect = TRUE) {
  type <- rep("identity", length(parnames))
  natural <- parnames
  if (detect) {
    for (key in names(transform_dict)) {
      prefix <- paste0(key, "_")
      hit <- type == "identity" &
        startsWith(parnames, prefix) &
        nchar(parnames) > nchar(prefix)
      type[hit] <- key
      natural[hit] <- substring(parnames[hit], nchar(prefix) + 1L)
    }
    # Undoing one clash can make another ("log_log_x" beside "log_x" and "x"),
    # so repeat until the names are unique. Each pass undoes at least one.
    clashed <- logical(length(parnames))
    repeat {
      dup <- natural %in% natural[duplicated(natural)] & type != "identity"
      if (!any(dup)) {
        break
      }
      clashed <- clashed | dup
      type[dup] <- "identity"
      natural[dup] <- parnames[dup]
    }
    if (any(clashed)) {
      left <- parnames[clashed]
      warn(c(
        "Left {length(left)} parameter{?s} on the unconstrained scale:
         {.val {left}}.",
        "i" = "Without its prefix the name is that of another parameter.
               Rename one of them to have it reported on the natural scale."
      ))
    }
  }
  data.frame(
    parameter = parnames,
    name = natural,
    transform = type,
    row.names = parnames,
    stringsAsFactors = FALSE
  )
}

# The transforms a method applies: those of the fit, or NULL when there are
# none, when the method was asked for the unconstrained scale, or when the
# fit predates them.
active_transforms <- function(x, transform = TRUE) {
  tr <- x$transform
  if (!transform || is.null(tr) || all(tr$transform == "identity")) {
    return(NULL)
  }
  tr
}

# The mean, SD, quantiles and mode of psi = inverse(theta) for theta ~
# SN(xi, omega, alpha). The quantiles map straight through the increasing
# inverse. The mean and SD of exp(theta) come from the moment generating
# function of the skew-normal, and the others from natural_moments(). The
# mode is that of the density of psi, not the image of the mode of theta.
natural_summary <- function(xi, omega, alpha, type, probs) {
  tf <- transform_dict[[type]]
  mom <- if (type == "log") {
    lognormal_moments(xi, omega, alpha)
  } else {
    natural_moments(xi, omega, alpha, tf$inverse)
  }
  c(
    mom$mean,
    mom$sd,
    tf$inverse(qsnorm(probs, xi, omega, alpha)),
    natural_mode(xi, omega, alpha, tf)
  )
}

# Mean and SD of exp(theta) for theta ~ SN(xi, omega, alpha), from
# E exp(t theta) = 2 exp(xi t + omega^2 t^2 / 2) Phi(delta omega t). The
# variance is the squared mean times expm1() of the log ratio of the second
# moment to the squared mean, which keeps its accuracy when omega is small.
lognormal_moments <- function(xi, omega, alpha) {
  d <- alpha / sqrt(1 + alpha^2) * omega
  log_mean <- log(2) + xi + omega^2 / 2 + stats::pnorm(d, log.p = TRUE)
  log_ratio <- omega^2 -
    log(2) +
    stats::pnorm(2 * d, log.p = TRUE) -
    2 * stats::pnorm(d, log.p = TRUE)
  mean <- exp(log_mean)
  list(mean = mean, sd = mean * sqrt(expm1(max(log_ratio, 0))))
}

# Mean and SD of inverse(theta) for theta ~ SN(xi, omega, alpha), for a
# bounded inverse, by the trapezoid rule on 10 SDs either side of xi. The
# integrand is smooth and dies off at both ends, where the rule converges
# faster than any power of the step. The step resolves the inverse, which
# varies over about one unit of theta, and the edge of a strongly skewed
# density, which is 1 / |alpha| wide.
natural_moments <- function(xi, omega, alpha, inverse) {
  step <- min(0.005, 0.1 / omega, 0.1 / (1 + abs(alpha)))
  z <- seq(-10, 10, length.out = min(ceiling(20 / step), 2e5) + 1)
  w <- dsnorm(z, 0, 1, alpha)
  w <- w / sum(w)
  psi <- inverse(xi + omega * z)
  mean <- sum(w * psi)
  list(mean = mean, sd = sqrt(sum(w * (psi - mean)^2)))
}

# Mode of the density of psi = inverse(theta) for theta ~ SN(xi, omega,
# alpha). On the scale of theta that density is the skew-normal's over the
# derivative of the inverse, and the log of that derivative has a slope of at
# most 2 in the dictionary, so its peak is within 2 omega SDs of the
# skew-normal's. A grid finds the highest peak (a wide posterior on the logit
# scale has one at each end of the unit interval), and optimize() refines it.
natural_mode <- function(xi, omega, alpha, tf) {
  zmax <- 10 + 2 * omega
  z <- seq(-zmax, zmax, length.out = min(ceiling(200 * zmax), 2e5) + 1)
  log_dens <- function(z) {
    dsnorm(z, 0, 1, alpha, log = TRUE) - tf$log_jac(xi + omega * z)
  }
  k <- which.max(log_dens(z))
  peak <- stats::optimize(
    log_dens,
    interval = z[c(max(k - 1L, 1L), min(k + 1L, length(z)))],
    maximum = TRUE,
    tol = 1e-10
  )$maximum
  tf$inverse(xi + omega * peak)
}

# Density of psi = inverse(theta) at `psi`, given the density of theta as a
# function `dens`: the change of variables f(theta) / (d psi / d theta).
natural_density <- function(psi, dens, tf) {
  theta <- tf$forward(psi)
  dens(theta) * exp(-tf$log_jac(theta))
}

# Draws on the natural scale: each transformed column of `theta` mapped
# through its inverse, and the columns renamed.
natural_draws <- function(theta, tr) {
  for (j in which(tr$transform != "identity")) {
    theta[, j] <- transform_dict[[tr$transform[j]]]$inverse(theta[, j])
  }
  colnames(theta) <- tr$name
  theta
}

# The line under a table that has natural-scale rows, saying which they are
# and how each was obtained: each one in full when there are at most `n`, and
# otherwise a count for each transform. The names are pasted as plain text and
# never pass through a cli template.
cat_transform_note <- function(tr, n = 4L) {
  rows <- which(tr$transform != "identity")
  calls <- vapply(
    names(transform_dict),
    function(k) transform_dict[[k]]$call,
    ""
  )
  items <- if (length(rows) <= n) {
    paste0(
      tr$name[rows],
      " = ",
      calls[tr$transform[rows]],
      "(",
      tr$parameter[rows],
      ")"
    )
  } else {
    count <- table(factor(tr$transform[rows], names(transform_dict)))
    count <- count[count > 0]
    paste0(
      calls[names(count)],
      "() of ",
      count,
      " ",
      names(count),
      "_ parameter",
      ifelse(count == 1, "", "s")
    )
  }
  cat_wrap(
    paste0(
      col_vanilla(cli::symbol$info),
      cli::col_grey(" Natural scale: ", paste(items, collapse = ", "), ".")
    ),
    exdent = 2
  )
}
