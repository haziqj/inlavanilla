# The marginals on workers, through `cores` or a future::plan(). Either way a
# fit must give the numbers of a serial fit: the same summaries, diagnostics,
# counts and warnings. A worker is only ever a saving of time, so anything that
# goes wrong on one, from a NaN to a killed process, must end in the serial
# answer.
#
# The fits here are tiny, so `force_workers()` says the scans will be long, to
# get past the checks that keep short scans serial and size the jobs. Forks
# (`multicore`) share the package under test, so they run under
# devtools::test() too. Separate sessions (`multisession`) load the installed
# package, so those tests only run on an installed copy.

skip_if_not_installed("future")

## ----- A small model ----------------------------------------------------------

set.seed(11)
reg_n <- 40
reg_x <- cbind(1, rnorm(reg_n), rnorm(reg_n))
reg_y <- drop(reg_x %*% c(1, 0.5, -0.5)) + rnorm(reg_n, 0, 0.8)
lp_reg <- function(theta) {
  eta <- drop(reg_x %*% theta[1:3])
  sum(dnorm(reg_y, eta, exp(theta[[4]]), log = TRUE)) +
    sum(dnorm(theta[1:3], 0, 10, log = TRUE)) +
    dnorm(theta[[4]], 0, 2, log = TRUE)
}
gr_reg <- function(theta) {
  eta <- drop(reg_x %*% theta[1:3])
  s2 <- exp(2 * theta[[4]])
  res <- reg_y - eta
  c(
    drop(crossprod(reg_x, res)) / s2 - theta[1:3] / 100,
    -reg_n + sum(res^2) / s2 - theta[[4]] / 4
  )
}
init_reg <- c(b0 = 0, b1 = 0, b2 = 0, ls = 0)
parent_pid <- Sys.getpid()

fit_reg <- function(lp = lp_reg, gr = NULL, ..., silent = TRUE) {
  inla_vanilla(lp, init_reg, gr, ..., silent = silent)
}

# Run `code` on a plan, and leave the plan as it was. Skips where the plan is
# not available and future would run sequentially instead (Windows, RStudio).
on_plan <- function(strategy, workers, code) {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(strategy, workers = workers)
  skip_if(future::nbrOfWorkers() < workers, "this plan is not available here")
  force(code)
}
on_forks <- function(code) on_plan(future::multicore, 2L, code)

# As if the scans of these tiny fits would take a long time
force_workers <- function(env = parent.frame()) {
  local_mocked_bindings(scan_calls = function(...) 1e12, .env = env)
}

expect_same_fit <- function(fit, ref) {
  expect_identical(fit$summary, ref$summary)
  expect_identical(fit$diagnostics$marginals, ref$diagnostics$marginals)
  expect_identical(fit$scan$logpost, ref$scan$logpost)
  expect_identical(fit$counts, ref$counts)
}

## ----- With no workers --------------------------------------------------------

test_that("with no plan, or a sequential one, the fit is a plain loop", {
  expect_identical(fit_reg()$workers, 1L)
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  expect_identical(par_workers(), 1L)
  force_workers()
  expect_identical(fit_reg()$workers, 1L)
})

test_that("a short scan stays serial, whatever the plan", {
  on_forks({
    fit <- fit_reg()
    expect_identical(fit$workers, 1L)
  })
})

test_that("Gaussian marginals have no scans to send anywhere", {
  force_workers()
  on_forks({
    fit <- fit_reg(marginals = "gaussian")
    expect_identical(fit$workers, 1L)
  })
})

## ----- On workers -------------------------------------------------------------

test_that("a fit on workers is the same as a serial fit", {
  force_workers()
  grads <- list(no_grad = NULL, grad = gr_reg)
  for (g in names(grads)) {
    for (w in c("trace", "cmp", "tk")) {
      ref <- fit_reg(gr = grads[[g]], width = w)
      fit <- on_forks(fit_reg(gr = grads[[g]], width = w))
      expect_identical(fit$workers, 2L, label = paste(g, w))
      expect_same_fit(fit, ref)
    }
  }
})

test_that("tk_fallback gives the same rescans on workers", {
  force_workers()
  set.seed(4)
  y <- rnorm(12, 10, 2)
  lp <- function(theta) {
    sum(dnorm(y, theta[["mu"]], exp(theta[["ls"]]), log = TRUE))
  }
  init <- c(mu = mean(y), ls = log(sd(y)))
  ref <- inla_vanilla(lp, init, tk_fallback = TRUE, silent = TRUE)
  fit <- on_forks(inla_vanilla(lp, init, tk_fallback = TRUE, silent = TRUE))
  expect_identical(fit$workers, 2L)
  expect_identical(fit$fallback, ref$fallback)
  expect_same_fit(fit, ref)
})

test_that("the progress line ticks once per parameter", {
  force_workers()
  ctrl <- vanilla_control(list())
  lap <- laplace_fit(lp_reg, NULL, init_reg, ctrl)
  ticks <- 0L
  par <- list(
    logpost = lp_reg,
    grad = NULL,
    hess = NULL,
    init = init_reg,
    call_s = 1,
    add_calls = function(n) NULL,
    note = function(text) FALSE
  )
  out <- on_forks(fit_marginals(
    lp_reg,
    NULL,
    lap,
    ctrl,
    tick = function() ticks <<- ticks + 1L,
    par = par
  ))
  expect_identical(out$workers, 2L)
  expect_identical(ticks, length(init_reg))
})

test_that("a warning on a worker is recorded as in a serial fit", {
  force_workers()
  edge <- fit_reg()$mode[[1]]
  lp_warn <- function(theta) {
    if (theta[[1]] > edge + 0.3) {
      warning("a warning from logpost far out")
    }
    lp_reg(theta)
  }
  ref <- fit_reg(lp_warn)
  expect_true(any(grepl("far out", ref$warnings)))
  fit <- on_forks(fit_reg(lp_warn))
  expect_identical(fit$workers, 2L)
  expect_identical(fit$warnings, ref$warnings)
})

test_that("verbose = TRUE says how many workers the scans used", {
  force_workers()
  msgs <- on_forks(capture_messages(fit_reg(verbose = TRUE, silent = FALSE)))
  expect_true(any(grepl("ran on 2 workers", msgs)))
})

## ----- When a worker cannot be trusted ----------------------------------------

test_that("a function that differs on a worker sends the fit back to serial", {
  force_workers()
  # What RTMB with random effects does in a separate session: NaN, no error
  lp_nan <- function(theta) {
    if (Sys.getpid() != parent_pid) NaN else lp_reg(theta)
  }
  ref <- fit_reg()
  msgs <- character()
  fit <- on_forks(withCallingHandlers(
    fit_reg(lp_nan, silent = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  expect_identical(fit$workers, 1L)
  expect_true(any(grepl("will run serially", msgs)))
  expect_same_fit(fit, ref)
})

test_that("a function that fails on a worker sends the fit back to serial", {
  force_workers()
  lp_err <- function(theta) {
    if (Sys.getpid() != parent_pid) {
      stop("no such function here")
    }
    lp_reg(theta)
  }
  fit <- on_forks(fit_reg(lp_err))
  expect_identical(fit$workers, 1L)
  expect_same_fit(fit, fit_reg())
})

test_that("a job that errors is done again in the parent", {
  force_workers()
  edge <- fit_reg()$mode[[1]]
  # Fine at the mode, where the check looks, and on a worker only
  lp_far <- function(theta) {
    if (Sys.getpid() != parent_pid && theta[[1]] > edge + 0.3) {
      stop("only on a worker, and only far out")
    }
    lp_reg(theta)
  }
  ref <- fit_reg()
  fit <- on_forks(fit_reg(lp_far))
  expect_identical(fit$workers, 2L)
  expect_identical(fit$summary, ref$summary)
  expect_identical(fit$diagnostics$marginals, ref$diagnostics$marginals)
})

test_that("a worker that dies is replaced by the parent", {
  force_workers()
  edge <- fit_reg()$mode[[1]]
  lp_kill <- function(theta) {
    if (Sys.getpid() != parent_pid && theta[[1]] > edge + 0.3) {
      tools::pskill(Sys.getpid(), tools::SIGKILL)
      Sys.sleep(5)
    }
    lp_reg(theta)
  }
  ref <- fit_reg()
  fit <- on_forks(fit_reg(lp_kill))
  expect_identical(fit$summary, ref$summary)
})

test_that("an error everywhere reads as in a serial fit", {
  force_workers()
  edge <- fit_reg()$mode[[1]]
  lp_bad <- function(theta) {
    if (theta[[1]] > edge + 0.3) {
      stop("logpost failed far out")
    }
    lp_reg(theta)
  }
  ref <- tryCatch(fit_reg(lp_bad), error = function(e) e)
  err <- tryCatch(on_forks(fit_reg(lp_bad)), error = function(e) e)
  expect_s3_class(ref, "inlavanilla_error")
  expect_s3_class(err, "inlavanilla_error")
  expect_identical(conditionMessage(err), conditionMessage(ref))
})

## ----- loo() and waic() ------------------------------------------------------

ll_reg <- function(theta) {
  dnorm(reg_y, drop(reg_x %*% theta[1:3]), exp(theta[[4]]), log = TRUE)
}

# As if the calls of loglik were slow
force_loo_workers <- function(env = parent.frame()) {
  local_mocked_bindings(taylor_est_s = function(...) 1e6, .env = env)
}

test_that("loo() and waic() on workers are the same as serially", {
  fit <- fit_reg()
  ref_loo <- loo(fit, ll_reg)
  ref_waic <- waic(fit, ll_reg)
  force_loo_workers()
  res_loo <- on_forks(loo(fit, ll_reg))
  res_waic <- on_forks(waic(fit, ll_reg))
  expect_identical(res_loo$estimates, ref_loo$estimates)
  expect_identical(res_loo$pointwise, ref_loo$pointwise)
  expect_identical(res_waic$estimates, ref_waic$estimates)
  expect_identical(res_waic$pointwise, ref_waic$pointwise)
})

test_that("the Taylor terms use the workers, and count every call", {
  fit <- fit_reg()
  ref <- taylor_pointwise(fit, ll_reg)
  expect_identical(ref$workers, 1L)
  force_loo_workers()
  tp <- on_forks(taylor_pointwise(fit, ll_reg))
  expect_identical(tp$workers, 2L)
  expect_identical(tp$calls, ref$calls)
  m <- length(init_reg)
  expect_identical(tp$calls, m * m + m + 1L)
  expect_identical(tp$l0, ref$l0)
  expect_identical(tp$lpd_2, ref$lpd_2)
})

test_that("the points of the Taylor terms are the ones the differences need", {
  U <- taylor_points(3L, 0.5)
  expect_identical(dim(U), c(3L, 13L))
  expect_identical(U[, 1], c(0, 0, 0))
  expect_identical(U[, 2], c(0.5, 0, 0))
  expect_identical(U[, 3], c(-0.5, 0, 0))
  expect_identical(U[, 6], c(0, 0, 0.5))
  expect_identical(U[, 7], c(0, 0, -0.5))
  expect_identical(U[, 8], c(0.5, 0.5, 0))
  expect_identical(U[, 13], c(0, -0.5, -0.5))
})

test_that("a loglik that differs on a worker runs serially", {
  fit <- fit_reg()
  force_loo_workers()
  ll_nan <- function(theta) {
    if (Sys.getpid() != parent_pid) NaN * reg_y else ll_reg(theta)
  }
  tp <- on_forks(taylor_pointwise(fit, ll_nan))
  expect_identical(tp$workers, 1L)
  expect_identical(tp$l0, taylor_pointwise(fit, ll_reg)$l0)
})

test_that("an unusable loglik on a worker gets the usual error", {
  fit <- fit_reg()
  ref_err <- function(f) tryCatch(f, error = function(e) conditionMessage(e))
  # Fine at the centre, where the check looks
  centre <- fit$mean_vb
  ll_bad <- function(theta) {
    v <- ll_reg(theta)
    if (
      !isTRUE(all.equal(unname(theta), unname(centre))) && v[[1]] < v[[1]] + 1
    ) {
      v[[2]] <- NaN
    }
    v
  }
  serial <- ref_err(taylor_pointwise(fit, ll_bad))
  force_loo_workers()
  on_workers <- ref_err(on_forks(taylor_pointwise(fit, ll_bad)))
  expect_match(serial, "not finite")
  expect_identical(on_workers, serial)
})

test_that("a warning from loglik on a worker is recorded", {
  fit <- fit_reg()
  ll_warn <- function(theta) {
    warning("a warning from loglik")
    ll_reg(theta)
  }
  force_loo_workers()
  tp <- on_forks(taylor_pointwise(fit, ll_warn))
  expect_identical(tp$workers, 2L)
  expect_true(any(grepl("a warning from loglik", tp$loglik_warnings)))
})

## ----- The cores argument -----------------------------------------------------

test_that("cores = n runs on workers, the same as serially, and leaves no plan", {
  force_workers()
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  skip_if(!par_can_fork(), "needs forks to run under devtools::test()")
  ref <- fit_reg()
  fit <- fit_reg(cores = 2)
  expect_identical(fit$workers, 2L)
  expect_identical(fit$cores, 2L)
  expect_same_fit(fit, ref)
  expect_identical(future::nbrOfWorkers(), 1L)
})

test_that("a short scan with cores stays serial and starts no workers", {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  local_mocked_bindings(
    par_local_plan = function(cores) stop("no workers should be started")
  )
  expect_identical(fit_reg(cores = 2)$workers, 1L)
})

test_that("cores = 1 is serial whatever the plan", {
  force_workers()
  fit <- on_forks(fit_reg(cores = 1))
  expect_identical(fit$workers, 1L)
  expect_identical(fit$cores, 1L)
})

test_that("cores puts the user's own plan back", {
  force_workers()
  skip_if(!par_can_fork(), "needs forks to run under devtools::test()")
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  # A different kind of plan from the forks that `cores` chooses here
  future::plan(future::multisession, workers = 2L)
  fit <- fit_reg(cores = 2)
  expect_identical(fit$workers, 2L)
  expect_s3_class(future::plan(), "multisession")
})

test_that("cores never asks for more workers than the machine offers", {
  force_workers()
  skip_if(!par_can_fork(), "needs forks to run under devtools::test()")
  local_mocked_bindings(par_available = function() 2L)
  expect_identical(fit_reg(cores = 100)$workers, 2L)
})

test_that("cores must be a positive whole number or NULL", {
  for (bad in list(0, -1, 1.5, "a", c(1, 2), NA, Inf)) {
    expect_error(fit_reg(cores = bad), class = "inlavanilla_error")
  }
  expect_null(check_cores(NULL))
  expect_identical(check_cores(4), 4L)
})

test_that("cores takes what future::availableCores() returns", {
  n <- future::availableCores(omit = 1)
  expect_identical(check_cores(n), as.integer(n))
  expect_identical(fit_reg(cores = n)$workers, 1L) # a short scan, so serial
})

test_that("cores without the future package runs serially and says so", {
  force_workers()
  local_mocked_bindings(par_has_future = function() FALSE)
  msgs <- capture_messages(fit <- fit_reg(cores = 2, silent = FALSE))
  expect_identical(fit$workers, 1L)
  expect_true(any(grepl("needs the future package", gsub("\\s+", " ", msgs))))
})

test_that("cores uses separate sessions where forks are not safe", {
  # Positron's console refuses to fork
  local_mocked_bindings(par_program = function() "ark")
  expect_false(par_can_fork())
  # RStudio's console is unstable when forked, by future's own check
  Sys.setenv(RSTUDIO = "1", RSTUDIO_TERM = "")
  on.exit(Sys.unsetenv(c("RSTUDIO", "RSTUDIO_TERM")), add = TRUE)
  local_mocked_bindings(par_program = function() "r")
  expect_false(par_can_fork())
})

test_that("the POSITRON variable alone does not stop forks", {
  skip_if(!future::supportsMulticore(), "no forks here")
  local_mocked_bindings(par_program = function() "r")
  Sys.setenv(POSITRON = "1")
  on.exit(Sys.unsetenv("POSITRON"), add = TRUE)
  expect_true(par_can_fork())
})

test_that("loo() and waic() use the cores of the fit", {
  force_workers()
  force_loo_workers()
  skip_if(!par_can_fork(), "needs forks to run under devtools::test()")
  ref <- taylor_pointwise(fit_reg(), ll_reg)
  fit <- fit_reg(cores = 2)
  tp <- taylor_pointwise(fit, ll_reg)
  expect_identical(tp$workers, 2L)
  expect_identical(tp$lpd_2, ref$lpd_2)
  expect_identical(future::nbrOfWorkers(), 1L)
})

## ----- Choosing the number of jobs ---------------------------------------------

test_that("the jobs grow with the scan and stop at the workers and the items", {
  # The time with w jobs is est_s / w + cost * w, least at sqrt(est_s / cost)
  jobs <- function(est_s, cost, workers = 12L, n = 40L) {
    now <- 100
    clock <- function() {
      now <<- now + cost
      now
    }
    local_mocked_bindings(
      par_workers = function() workers,
      par_globals = function(functions) {
        list(found = list(), packages = character())
      },
      par_future = function(job, ids, data, setup) "a future",
      elapsed = clock
    )
    local_mocked_bindings(
      value = function(x, ...) 1,
      .package = "future"
    )
    par_setup(
      list(),
      function(ids, data) 1,
      list(),
      est_s = est_s,
      n = n,
      note = function(text) FALSE,
      label = "",
      what = ""
    )$jobs
  }
  expect_identical(jobs(est_s = 9, cost = 1), 3L)
  expect_identical(jobs(est_s = 2.3, cost = 0.14), 5L)
  expect_identical(jobs(est_s = 1e4, cost = 0.01), 12L)
  expect_identical(jobs(est_s = 1e4, cost = 0.01, n = 5L), 5L)
  expect_identical(jobs(est_s = 1, cost = 1), 1L)
})

test_that("the time a worker spends on its job does not reduce the jobs", {
  # Handing a job over takes 0.01 s and the job 5 s, as when a worker builds a
  # compiled model again. Workers do that side by side, so only the handing
  # over counts.
  now <- 100
  local_mocked_bindings(
    par_workers = function() 12L,
    par_globals = function(functions) {
      list(found = list(), packages = character())
    },
    par_future = function(job, ids, data, setup) {
      now <<- now + 0.01
      "a future"
    },
    elapsed = function() now
  )
  local_mocked_bindings(
    value = function(x, ...) {
      now <<- now + 5
      1
    },
    .package = "future"
  )
  setup <- par_setup(
    list(),
    function(ids, data) 1,
    list(),
    est_s = 100,
    n = 40L,
    note = function(text) FALSE,
    label = "",
    what = ""
  )
  expect_identical(setup$jobs, 12L)
})

test_that("the items are dealt out in turn, with none left out", {
  g <- par_groups(10L, 3L)
  expect_length(g, 3L)
  expect_identical(g[[1]], c(1L, 4L, 7L, 10L))
  expect_identical(sort(unlist(g)), 1:10)
})

test_that("the estimate of the calls follows the width and the route", {
  expect_identical(scan_calls(5, 21L, "trace", TRUE, FALSE), 22)
  expect_identical(scan_calls(5, 21L, "trace", FALSE, FALSE), 37)
  expect_gt(
    scan_calls(5, 21L, "cmp", TRUE, FALSE),
    scan_calls(5, 21L, "trace", TRUE, FALSE)
  )
  expect_gt(
    scan_calls(5, 21L, "tk", FALSE, FALSE),
    scan_calls(5, 21L, "tk", TRUE, FALSE)
  )
  expect_lt(
    scan_calls(5, 21L, "cmp", TRUE, TRUE),
    scan_calls(5, 21L, "cmp", TRUE, FALSE)
  )
})

## ----- The hint ---------------------------------------------------------------

test_that("a long interactive fit with no plan says once that a plan would help", {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  the$hinted <- FALSE
  on.exit(the$hinted <- FALSE, add = TRUE)
  local_mocked_bindings(is_interactive = function() TRUE, par_hint_s = 0)
  first <- capture_messages(fit_reg(silent = FALSE))
  expect_true(any(grepl("`cores`", first)))
  second <- capture_messages(fit_reg(silent = FALSE))
  expect_false(any(grepl("`cores`", second)))
})

test_that("the hint waits for a fit that can show it", {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  the$hinted <- FALSE
  on.exit(the$hinted <- FALSE, add = TRUE)
  local_mocked_bindings(is_interactive = function() TRUE, par_hint_s = 0)
  fit_reg()
  expect_false(the$hinted)
  first <- capture_messages(fit_reg(silent = FALSE))
  expect_true(any(grepl("`cores`", first)))
})

test_that("rendered documents and short fits get no hint", {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  the$hinted <- FALSE
  on.exit(the$hinted <- FALSE, add = TRUE)
  local_mocked_bindings(is_interactive = function() FALSE, par_hint_s = 0)
  expect_false(any(grepl(
    "future::plan",
    capture_messages(fit_reg(silent = FALSE))
  )))
  local_mocked_bindings(is_interactive = function() TRUE, par_hint_s = 1e6)
  expect_false(any(grepl(
    "future::plan",
    capture_messages(fit_reg(silent = FALSE))
  )))
})

## ----- Separate sessions ----------------------------------------------------------

test_that("cores runs separate sessions where forks are not safe", {
  skip_if(
    exists(
      ".__DEVTOOLS__",
      envir = asNamespace("inlavanilla"),
      inherits = FALSE
    ),
    "needs an installed copy"
  )
  force_workers()
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  local_mocked_bindings(par_program = function() "ark")
  ref <- fit_reg()
  fit <- fit_reg(cores = 2)
  expect_identical(fit$workers, 2L)
  expect_same_fit(fit, ref)
  expect_identical(future::nbrOfWorkers(), 1L)
})

test_that("separate sessions find global helpers and data passed in the dots", {
  # pkgload marks the namespace of a package it has loaded from source
  skip_if(
    exists(
      ".__DEVTOOLS__",
      envir = asNamespace("inlavanilla"),
      inherits = FALSE
    ),
    "needs an installed copy"
  )
  force_workers()
  # As a user's script has them: in the global environment
  assign(".par_test_sd", 10, envir = globalenv())
  assign(
    ".par_test_prior",
    function(b) sum(dnorm(b, 0, .par_test_sd, log = TRUE)),
    envir = globalenv()
  )
  on.exit(rm(.par_test_sd, .par_test_prior, envir = globalenv()), add = TRUE)
  lp <- function(theta, x, y) {
    sum(dnorm(y, drop(x %*% theta[1:3]), exp(theta[[4]]), log = TRUE)) +
      .par_test_prior(theta[1:3]) +
      dnorm(theta[[4]], 0, 2, log = TRUE)
  }
  environment(lp) <- globalenv()
  ref <- inla_vanilla(lp, init_reg, x = reg_x, y = reg_y, silent = TRUE)
  fit <- on_plan(
    future::multisession,
    2L,
    inla_vanilla(lp, init_reg, x = reg_x, y = reg_y, silent = TRUE)
  )
  expect_identical(fit$workers, 2L)
  expect_same_fit(fit, ref)
})
