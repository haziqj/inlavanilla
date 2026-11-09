# Parallel evaluation through the future package, which is optional. The user's
# future::plan() is the switch, as in blavaan: with no plan, or a sequential
# one, nothing here runs and the fit is a plain loop. Whatever the plan, a fit
# gives the numbers of a serial fit. The work divides into items that do not
# depend on each other, and anything that goes wrong on a worker is done again
# in the parent, so parallelism only ever saves time.
#
# The items are spread over the jobs at the start, at most one job per worker,
# and each job runs to its end. A job signals the parent after each item, so
# the parent's progress line moves one item at a time, as in a serial fit. See
# lab/DESIGN.md.

# Below this estimate of the serial time, in seconds, the work is not worth
# sending to workers: a job costs tens of milliseconds at least, and on the
# fits tried the gain was a few hundred milliseconds at best.
par_min_s <- 1

# Above this estimate, with no plan set, an interactive fit says once per
# session that a plan would help.
par_hint_s <- 10

# What is kept for the session
the <- new.env(parent = emptyenv())
the$hinted <- FALSE

# A function of our own, so that the tests can replace it
is_interactive <- function() interactive()

# The workers the user's plan gives. Without the future namespace loaded nobody
# has set a plan, so it is not loaded here, which would cost a whole small fit.
par_workers <- function() {
  if (!isNamespaceLoaded("future")) {
    return(1L)
  }
  n <- tryCatch(future::nbrOfWorkers(), error = function(e) 1L)
  if (is.finite(n)) as.integer(n) else .Machine$integer.max
}

# What the items need besides their arguments: anything the user's functions
# use from the global environment, and the packages they use. `functions`
# holds the user's functions, NULL for those not given. future finds these by
# reading the functions, once, because reading them again for each job is slow.
par_globals <- function(functions) {
  functions <- functions[!vapply(functions, is.null, NA)]
  names(functions) <- paste0(".par_fn", seq_along(functions))
  expr <- as.call(c(list(as.name("list")), lapply(names(functions), as.name)))
  gp <- future::getGlobalsAndPackages(
    expr,
    envir = list2env(functions, parent = globalenv())
  )
  found <- setdiff(
    names(gp$globals),
    c(names(functions), "...", "future.call.arguments", "list")
  )
  list(
    found = stats::setNames(lapply(found, function(k) gp$globals[[k]]), found),
    packages = unique(c(gp$packages, "inlavanilla"))
  )
}

# One job as a future: `job(ids, data)` on a worker. The function and its data
# travel as one object each, so an environment they share is sent once. The
# names are ours so that they cannot mask a global the user's function needs.
par_future <- function(job, ids, data, setup) {
  future::future(
    quote(.par_job(.par_ids, .par_data)),
    substitute = FALSE,
    globals = c(
      setup$found,
      list(.par_job = job, .par_ids = ids, .par_data = data)
    ),
    packages = setup$packages,
    seed = NULL
  )
}

# Decide how to run `n` items of about `est_s` serial seconds in all. `probe`
# is run with `data` on a worker and here, and the two values must agree: the
# user's function may not survive the trip to a separate R session, and may
# then fail or, worse, return NaN with no error (compiled code does, depending
# on how it was built). `label` and `what` word the messages. Returns `jobs`,
# 1 to stay serial, and what par_future() needs.
#
# The elapsed time of the probe is the real cost of reaching one worker, with
# the real data. The time with w jobs is about est_s / w + cost * w, which is
# least at w = sqrt(est_s / cost), so a long scan uses every worker of the plan
# and a shorter one fewer.
par_setup <- function(functions, probe, data, est_s, n, note, label, what) {
  serial <- list(jobs = 1L)
  workers <- par_workers()
  if (workers < 2L) {
    par_hint(est_s, note, label)
    return(serial)
  }
  if (est_s < par_min_s || n < 2L) {
    return(serial)
  }
  # Only the values matter here. Anything the function says while being
  # checked it has said already, or will say again in the real run
  checked <- suppressMessages(suppressWarnings(tryCatch(
    {
      setup <- par_globals(functions)
      t0 <- elapsed()
      got <- future::value(par_future(probe, 0L, data, setup))
      cost <- elapsed() - t0
      isTRUE(all.equal(
        got,
        probe(0L, data),
        tolerance = 1e-6,
        check.attributes = FALSE
      ))
    },
    error = function(e) FALSE
  )))
  if (!checked) {
    note(paste0(
      label,
      " will run serially, because `",
      what,
      "` did not give the same value on a parallel worker. Compiled code ",
      "(Rcpp, TMB, RTMB) cannot be sent to a separate R session, but ",
      "future::plan(future::multicore) can use it where forking is available."
    ))
    return(serial)
  }
  jobs <- min(workers, n, ceiling(sqrt(est_s / max(cost, 1e-3))))
  c(list(jobs = if (jobs < 2L) 1L else as.integer(jobs)), setup)
}

# Say once per session, at the console, that a plan would help. `note` returns
# whether the message was shown, so a silent fit does not use the hint up.
par_hint <- function(est_s, note, label) {
  if (!the$hinted && est_s >= par_hint_s && is_interactive()) {
    the$hinted <- isTRUE(note(paste0(
      label,
      " should take about ",
      format_clock(est_s),
      ". To use more cores, call future::plan(future::multisession) first."
    )))
  }
}

# Signal the parent from a worker that an item is done. future relays an
# immediate condition while the parent polls, and not only at the end. Some
# backends relay it again with the job's result, so the parent counts each item
# once.
par_tick <- function(item) {
  signalCondition(structure(
    class = c("inlavanilla_tick", "immediateCondition", "condition"),
    list(message = "", call = NULL, item = item)
  ))
  invisible()
}

# The items split evenly over the jobs, in turn, so that items that are alike
# and next to each other, as parameters often are, do not land on one worker.
par_groups <- function(n, jobs) {
  unname(split(seq_len(n), rep_len(seq_len(jobs), n)))
}

# Run one job per group and wait for all of them, calling `on_item` as items
# finish. Returns a list with the value of each job, or NULL for one that
# failed, for any reason, which the caller then does itself. Warnings raised on
# the workers are signalled here by value().
par_run <- function(groups, job, data, setup, on_item) {
  withCallingHandlers(
    {
      running <- lapply(groups, function(ids) {
        tryCatch(par_future(job, ids, data, setup), error = function(e) NULL)
      })
      started <- running[!vapply(running, is.null, NA)]
      tryCatch(
        while (!all(future::resolved(started))) {
          Sys.sleep(0.02)
        },
        error = function(e) NULL
      )
      lapply(running, function(f) {
        if (is.null(f)) {
          return(NULL)
        }
        tryCatch(future::value(f), error = function(e) NULL)
      })
    },
    inlavanilla_tick = function(cond) on_item(cond$item)
  )
}
