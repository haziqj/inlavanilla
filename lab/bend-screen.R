################################################################################
#
# The ridge bend without a gradient (November 2026). The full check takes a
# one-sided gradient in the other parameters at +-2 SDs, 2(m - 1) calls per
# parameter. scan_bend() now tries first the direction that gamma1_lp()'s calls
# already give, 2 calls per parameter, and takes the full gradient only when
# the gain along it is above `bend_screen`. This script fits each model three
# ways, by setting that level: 0 (the full check on every parameter), 0.1 (the
# default) and Inf (the direction alone), and compares the ridge bends and the
# calls of logpost.
#
# Result: the flag agrees on 350 of 351 parameters of 28 models, with no false
# flags, and the marginals do not change. The one miss is the intercept of a
# random-intercept model whose joint mode has collapsed (log tau = -29.5),
# where the other nine parameters are flagged. Where the full value is above
# 0.02 the median ratio to it is 1.000 and the 5th percentile 0.965. The
# direction alone misses three flags, among them the only one on the SEM.
#
# Takes under a minute, and needs lavaan, Rcpp, RcppEigen and knitr for the
# SEM. Run from the package root:
#   Rscript lab/bend-screen.R
#
################################################################################

## ----- Configuration ----------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE))
options(width = 160)
flag_level <- summary_levels$ridge_bend

set_screen <- function(value) {
  ns <- asNamespace("inlavanilla")
  unlockBinding("bend_screen", ns)
  assign("bend_screen", value, envir = ns)
}

## =============================================================================
## Models
## =============================================================================

## ----- The small models of lab/bend-calibration.R ------------------------------

script <- readLines("lab/bend-calibration.R")
from <- grep("^## Small models", script) + 1L
to <- grep("^## Political Democracy SEM", script) - 2L
compare_grid <- function(case, logpost, init) {
  list(name = case, logpost = logpost, init = init)
}
eval(parse(text = script[from:to]))
cases <- unname(res)
add <- function(name, logpost, init) {
  cases[[length(cases) + 1L]] <<- list(
    name = name,
    logpost = logpost,
    init = init
  )
}

## ----- The Political Democracy SEM ---------------------------------------------

article <- "vignettes/articles/rcpp.qmd"
purled <- tempfile(fileext = ".R")
invisible(knitr::purl(article, purled, documentation = 0, quiet = TRUE))
code <- readLines(purled)
code <- code[seq_len(grep("^fit_inlavaan <- ", code)[1] - 1)]
code <- code[!grepl("^library\\(inlavanilla\\)", code)]
eval(parse(text = code))
qmd <- readLines(article)
fences <- grep("^```", qmd)
cpp_starts <- grep("^```\\{cat, engine.opts", qmd)
cpp <- unlist(lapply(cpp_starts, function(s) {
  qmd[(s + 1):(fences[fences > s][1] - 1)]
}))
cpp_file <- tempfile(fileext = ".cpp")
writeLines(cpp, cpp_file)
Rcpp::sourceCpp(cpp_file)
add("PolDem SEM", function(theta) logpost_cpp(theta, Y), init)

## ----- Larger models, with one or many scale parameters ------------------------

linreg <- function(p, n) {
  set.seed(1)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  y <- drop(X %*% rnorm(p, 0, 1)) + rnorm(n, 0, 1.5)
  add(
    sprintf("linear regression, p = %d, n = %d", p, n),
    function(theta) {
      b <- theta[seq_len(p)]
      sum(dnorm(y, drop(X %*% b), exp(theta[[p + 1L]]), log = TRUE)) +
        sum(dnorm(b, 0, 10, log = TRUE)) +
        dnorm(theta[[p + 1L]], 0, 5, log = TRUE)
    },
    setNames(c(rep(0, p), 0), c(paste0("b", seq_len(p)), "log_sigma"))
  )
}
linreg(10, 20)
linreg(20, 40)
linreg(40, 90)

hetero <- function(p, g, n) {
  set.seed(2)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  grp <- rep_len(seq_len(g), n)
  y <- drop(X %*% rnorm(p, 0, 1)) +
    rnorm(n, 0, exp(seq(-0.5, 0.5, length.out = g))[grp])
  add(
    sprintf("regression, %d noise scales, p = %d, n = %d", g, p, n),
    function(theta) {
      b <- theta[seq_len(p)]
      ls <- theta[p + seq_len(g)]
      sum(dnorm(y, drop(X %*% b), exp(ls[grp]), log = TRUE)) +
        sum(dnorm(b, 0, 10, log = TRUE)) +
        sum(dnorm(ls, 0, 3, log = TRUE))
    },
    setNames(
      numeric(p + g),
      c(paste0("b", seq_len(p)), paste0("ls", seq_len(g)))
    )
  )
}
hetero(6, 6, 48)

pois_ri <- function(J, nj) {
  set.seed(3)
  u <- rnorm(J, 0, 0.7)
  grp <- rep(seq_len(J), each = nj)
  y <- rpois(length(grp), exp(0.8 + u[grp]))
  add(
    sprintf("Poisson random intercepts, J = %d, nj = %d", J, nj),
    function(theta) {
      uu <- theta[-(1:2)]
      sum(dpois(y, exp(theta[[1L]] + uu[grp]), log = TRUE)) +
        sum(dnorm(uu, 0, exp(theta[[2L]]), log = TRUE)) +
        dnorm(theta[[1L]], 0, 5, log = TRUE) +
        dnorm(theta[[2L]], 0, 2, log = TRUE)
    },
    setNames(numeric(J + 2L), c("b0", "log_tau", paste0("u", seq_len(J))))
  )
}
pois_ri(8, 4)
pois_ri(20, 5)
pois_ri(40, 3)

logit <- function(m, n) {
  set.seed(1)
  X <- cbind(1, matrix(rnorm(n * (m - 1)), n, m - 1))
  y <- rbinom(n, 1, plogis(drop(X %*% c(-0.5, rnorm(m - 1, 0, 0.5)))))
  add(
    sprintf("logistic, m = %d, n = %d", m, n),
    function(theta) {
      eta <- drop(X %*% theta)
      sum(y * eta - log1p(exp(eta))) + sum(dnorm(theta, 0, 5, log = TRUE))
    },
    setNames(numeric(m), paste0("b", seq_len(m)))
  )
}
logit(5, 100)
logit(20, 200)
logit(40, 500)

negbin <- function(p, n) {
  set.seed(4)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  y <- rnbinom(n, size = 2, mu = exp(drop(X %*% c(1, rnorm(p - 1, 0, 0.3)))))
  add(
    sprintf("negative binomial, p = %d, n = %d", p, n),
    function(theta) {
      b <- theta[seq_len(p)]
      sum(dnbinom(
        y,
        size = exp(theta[[p + 1L]]),
        mu = exp(drop(X %*% b)),
        log = TRUE
      )) +
        sum(dnorm(b, 0, 5, log = TRUE)) +
        dnorm(theta[[p + 1L]], 0, 2, log = TRUE)
    },
    setNames(numeric(p + 1L), c(paste0("b", seq_len(p)), "log_size"))
  )
}
negbin(8, 40)
negbin(20, 120)

student <- function(p, n) {
  set.seed(5)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  y <- drop(X %*% rnorm(p, 0, 0.7)) + 1.2 * rt(n, 4)
  add(
    sprintf("Student-t regression, p = %d, n = %d", p, n),
    function(theta) {
      b <- theta[seq_len(p)]
      s <- exp(theta[[p + 1L]])
      nu <- 1 + exp(theta[[p + 2L]])
      sum(dt((y - drop(X %*% b)) / s, nu, log = TRUE)) -
        n * theta[[p + 1L]] +
        sum(dnorm(b, 0, 10, log = TRUE)) +
        dnorm(theta[[p + 1L]], 0, 3, log = TRUE) +
        dnorm(theta[[p + 2L]], 1, 1, log = TRUE)
    },
    setNames(
      numeric(p + 2L),
      c(paste0("b", seq_len(p)), "log_sigma", "log_nu1")
    )
  )
}
student(6, 30)
student(15, 60)

anova_vc <- function(J, nj) {
  set.seed(6)
  a <- rnorm(J, 0, 1)
  grp <- rep(seq_len(J), each = nj)
  y <- 5 + a[grp] + rnorm(length(grp), 0, 2)
  ybar <- tapply(y, grp, mean)
  S <- sum((y - ybar[grp])^2)
  N <- length(y)
  add(
    sprintf("variance components, J = %d, nj = %d", J, nj),
    function(theta) {
      s2 <- exp(2 * theta[[2L]])
      t2 <- exp(2 * theta[[3L]])
      v <- t2 + s2 / nj
      -0.5 *
        (N - J) *
        log(s2) -
        0.5 * S / s2 -
        0.5 * J * log(nj) -
        0.5 * J * log(v) -
        0.5 * sum((ybar - theta[[1L]])^2) / v -
        0.5 * N * log(2 * pi) +
        dnorm(theta[[1L]], 0, 20, log = TRUE) +
        dnorm(theta[[2L]], 0, 2, log = TRUE) +
        dnorm(theta[[3L]], 0, 2, log = TRUE)
    },
    c(mu = mean(y), log_sigma = 0.5, log_tau = 0)
  )
}
anova_vc(6, 4)
anova_vc(12, 3)

## =============================================================================
## Fits
## =============================================================================

fit_at <- function(case, level) {
  set_screen(level)
  fit <- inla_vanilla(case$logpost, case$init, silent = TRUE)
  fit
}

rows <- lapply(cases, function(case) {
  full <- fit_at(case, 0)
  scr <- fit_at(case, 0.1)
  dir <- fit_at(case, Inf)
  data.frame(
    case = case$name,
    m = length(case$init),
    par = names(case$init),
    full = full$diagnostics$marginals$ridge_bend,
    screen = scr$diagnostics$marginals$ridge_bend,
    dir = dir$diagnostics$marginals$ridge_bend,
    sum_diff = max(abs(full$summary - scr$summary)),
    calls_full = full$counts[["logpost"]],
    calls_screen = scr$counts[["logpost"]],
    calls_dir = dir$counts[["logpost"]],
    log_tau_mode = if ("log_tau" %in% names(case$init)) {
      full$mode[["log_tau"]]
    } else {
      NA_real_
    }
  )
})
d <- do.call(rbind, rows)
set_screen(0.1)

## =============================================================================
## Results
## =============================================================================

cat("Parameters:", nrow(d), " models:", length(unique(d$case)), "\n")
cat(
  "Marginal summaries unchanged, largest difference:",
  max(d$sum_diff),
  "\n\n"
)

for (v in c("screen", "dir")) {
  cat(v, ", flag at ", flag_level, ":\n", sep = "")
  print(table(full = d$full > flag_level, new = d[[v]] > flag_level))
  big <- d$full > 0.02
  cat("ratio to the full value where it is above 0.02:\n")
  print(round(
    quantile(
      d[[v]][big] / d$full[big],
      c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
    ),
    3
  ))
  cat("\n")
}

models <- do.call(
  rbind,
  lapply(split(d, d$case), function(x) {
    data.frame(
      case = x$case[1L],
      m = x$m[1L],
      flag_full = sum(x$full > flag_level),
      flag_screen = sum(x$screen > flag_level),
      refined = sum(x$screen != x$dir),
      calls_full = x$calls_full[1L],
      calls_screen = x$calls_screen[1L],
      saved_pct = round(100 * (1 - x$calls_screen[1L] / x$calls_full[1L])),
      row.names = NULL
    )
  })
)
print(models[order(models$m, models$case), ], row.names = FALSE)

cat("\nDisagreements at the flag level:\n")
print(
  d[
    (d$full > flag_level) != (d$screen > flag_level),
    c("case", "par", "full", "screen", "log_tau_mode")
  ],
  row.names = FALSE
)
cat("\nThe random-intercept fits, joint mode of log tau:\n")
print(
  unique(d[!is.na(d$log_tau_mode), c("case", "log_tau_mode")]),
  row.names = FALSE
)
