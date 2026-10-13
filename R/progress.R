# Progress reporting. Each stage gets a cli spinner or progress bar that
# appears only when the stage runs for longer than half a second, and ends
# with one coloured line. The lines are built directly from cli's styles
# rather than through its glue-based formatters, which cost a few
# milliseconds a call: as long as a whole fit of a small model.

progress_new <- function(verbose, frame = parent.frame()) {
  p <- new.env(parent = emptyenv())
  p$verbose <- verbose
  p$frame <- frame
  p$id <- NULL
  p$calls <- FALSE # whether the log-posterior wrapper should tick the bar
  p
}

# Start a stage. With `total`, a bar over that many items; without, a spinner
# that counts log-posterior calls.
progress_start <- function(p, msg, total = NA) {
  if (!p$verbose) {
    return(invisible(p))
  }
  old <- options(cli.progress_show_after = 0.5)
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
  if (!p$verbose) {
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
