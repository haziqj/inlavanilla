# The times of the stages of a fit (R/timing.R): timing() and its print()
# method. The times themselves vary from run to run, so the tests check the
# segments, their sum and the format, not the values.

lp_timing <- function(theta) -0.5 * sum((theta - 1)^2)
fit_timing <- inla_vanilla(lp_timing, c(a = 0, b = 0), silent = TRUE)
segments <- c("init", "laplace", "vb", "marginals", "diagnostics", "total")

test_that("timing() gives the total by default", {
  tt <- timing(fit_timing)
  expect_s3_class(tt, "inla_vanilla_timing")
  expect_identical(names(tt), "total")
  expect_identical(unclass(tt), fit_timing$timing["total"])
})

test_that("timing(what = \"all\") gives every segment, which add up to the total", {
  tt <- timing(fit_timing, what = "all")
  expect_identical(names(tt), segments)
  expect_true(all(is.finite(tt) & tt >= 0))
  expect_equal(sum(tt[segments != "total"]), tt[["total"]])

  # The same segments without VB and with Gaussian marginals
  fit <- inla_vanilla(
    lp_timing,
    c(a = 0, b = 0),
    vb = FALSE,
    marginals = "gaussian",
    silent = TRUE
  )
  tt <- timing(fit, what = "all")
  expect_identical(names(tt), segments)
  expect_equal(sum(tt[segments != "total"]), tt[["total"]])
})

test_that("timing() gives the requested segments, in the order asked", {
  tt <- timing(fit_timing, what = c("marginals", "vb"))
  expect_identical(names(tt), c("marginals", "vb"))
  expect_identical(unclass(tt), fit_timing$timing[c("marginals", "vb")])
})

test_that("timing() rejects unknown segments and a bad what", {
  expect_error(
    timing(fit_timing, what = c("vb", "optim", "loo")),
    "Unknown timing segments: \"optim\" and \"loo\"",
    class = "inlavanilla_error"
  )
  expect_error(
    timing(fit_timing, what = "optim"),
    "any of \"init\", \"laplace\"",
    class = "inlavanilla_error"
  )
  # "all" only on its own
  expect_error(
    timing(fit_timing, what = c("all", "vb")),
    "Unknown timing segment: \"all\"",
    class = "inlavanilla_error"
  )
  for (bad in list(1, character(), NA_character_, NULL)) {
    expect_error(
      timing(fit_timing, what = bad),
      "`what` must name timing segments",
      class = "inlavanilla_error"
    )
  }
})

test_that("timing() prints each time as summary() does", {
  tt <- timing(fit_timing, what = "all")
  out <- capture_output(res <- withVisible(print(tt)))
  expect_false(res$visible)
  expect_identical(res$value, tt)
  lines <- strsplit(out, "\n")[[1]]
  expect_identical(strsplit(trimws(lines[1]), " +")[[1]], segments)
  expect_identical(
    strsplit(trimws(lines[2]), " +")[[1]],
    unname(vapply(unclass(tt), format_time, character(1)))
  )
  # Milliseconds, seconds, then minutes and seconds
  tt[] <- c(0.0123, 1.23, 75, 0, 0, 76.2533)
  expect_match(capture_output(print(tt)), "12ms +1\\.2s +1m 15s")
})

test_that("timing() of other objects is an error, or INLAvaan's", {
  expect_error(
    timing(1:3),
    "There is no `timing\\(\\)` method for an integer vector",
    class = "inlavanilla_error"
  )
  # An S4 object is passed on to INLAvaan's timing() when INLAvaan is loaded,
  # with its arguments, and anything else is not
  s4 <- asS4(structure(list(), class = "not_inlavaan"))
  local({
    local_mocked_bindings(
      inlavaan_timing = function() function(object, ...) list(object, ...)
    )
    expect_identical(timing(s4, what = "all"), list(s4, what = "all"))
    expect_error(timing(1:3), class = "inlavanilla_error")
  })
  # Without INLAvaan, an S4 object is an error too
  skip_if(isNamespaceLoaded("INLAvaan"))
  expect_error(
    timing(s4),
    "There is no `timing\\(\\)` method",
    class = "inlavanilla_error"
  )
})
