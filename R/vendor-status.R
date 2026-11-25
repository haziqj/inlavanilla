#' Check a copy of the engine of inlavanilla
#'
#' @description
#' Checks the copy of the engine that [use_vanilla()] put into the package
#' whose root is `path`, and prints one line per check: a tick when it passes,
#' a cross when it fails, and an exclamation mark for something to look at,
#' each with what to do. The checks, in this order:
#' * **Lock**: `inst/inlavanilla.lock` exists and can be read.
#' * **Copies**: each engine file has the MD5 hash in the lock, and every file
#'   `R/vanilla-*.R` is in the lock. The hash is taken after carriage returns
#'   are removed, so that a checkout with Windows line endings passes.
#' * **Names**: no object of an engine file is also defined in another file of
#'   `R/`, or imported by name in `NAMESPACE`. R keeps one of the two with no
#'   warning, so a clash breaks the engine or your code silently.
#' * **Dependencies**: `DESCRIPTION` has `cli (>= 3.6.0)`, `stats` and `utils`
#'   in `Imports` (or `Depends`), `future` and `parallelly` in `Suggests` (or
#'   above), and asks for R 4.1.0 or later.
#' * **Authors**: `Authors@R` names the copyright holders of the engine, with
#'   the role `"cph"`.
#' * **License**: the package's `License` allows code under the GNU General
#'   Public License, version 3 or later.
#' * **Release**: the release in the lock is a tag on GitHub at the commit in
#'   the lock, and each engine file at that tag has the hash in the lock.
#' * **Update**: no newer release of inlavanilla changes a copied file. A
#'   newer release that does is reported, with the command that updates the
#'   copy, and is not a failure.
#'
#' The last two checks need GitHub, or a local clone of inlavanilla given as
#' `clone`. Without either they are skipped with a note, and the others still
#' run.
#'
#' @param path The root of the package, or a folder inside it. By default the
#'   working directory.
#' @param clone The path of a local clone of inlavanilla, for the last two
#'   checks without the internet. It needs git, and its tags must be up to date
#'   (`git fetch --tags`). By default they use GitHub.
#' @param error If `TRUE`, stop with an error when a check fails, as continuous
#'   integration needs. A skipped check or an exclamation mark is never an
#'   error.
#'
#' @returns Invisibly, a data frame with a row per check, in the order above:
#'   `check` (`"lock"`, `"copies"`, `"names"`, `"dependencies"`, `"authors"`,
#'   `"license"`, `"release"` or `"update"`), `status` (`"pass"`, `"warn"`,
#'   `"fail"` or `"skip"`), `message`, the line printed for it, and `details`,
#'   the lines printed under it, as plain text.
#'
#' @seealso [use_vanilla()] to copy and update the engine, and
#'   [use_vanilla_ci()] to run these checks on GitHub.
#'
#' @examples
#' \dontrun{
#' # From the root of a package that holds a copy of the engine
#' vanilla_status()
#'
#' # In continuous integration: an error when a check fails
#' vanilla_status(error = TRUE)
#' }
#' @export
vanilla_status <- function(path = ".", clone = NULL, error = FALSE) {
  check_flag(error, "error")
  clone <- vendor_check_clone(clone)
  root <- vendor_root(path)
  desc <- vendor_desc(root)
  pkg <- vendor_package(root, desc)
  vendor_start()
  rows <- vendor_status_rows()

  if (identical(pkg, "inlavanilla")) {
    vendor_say(list(vendor_msg(
      "i",
      "This is inlavanilla, which holds the originals of the engine files:
       nothing to check."
    )))
    return(invisible(vendor_status_frame(rows)))
  }

  cli::cli_h2("Engine copy in {pkg}")

  # The lock
  lock_path <- file.path(root, vendor_lock_file)
  lock <- NULL
  if (!file.exists(lock_path)) {
    rows$lock <- vendor_row(
      "fail",
      vendor_msg(
        "x",
        "There is no lock {.file inst/inlavanilla.lock}, so {pkg} holds no copy
         of the engine that {.fn use_vanilla} made."
      ),
      vendor_msg(">", "Copy the engine with {.run inlavanilla::use_vanilla()}.")
    )
  } else {
    lock <- tryCatch(
      vendor_read_lock(lock_path),
      inlavanilla_lock_error = function(e) e
    )
    if (inherits(lock, "error")) {
      reason <- lock$reason
      rows$lock <- vendor_row(
        "fail",
        vendor_msg(
          "x",
          "The lock {.file inst/inlavanilla.lock} is malformed: {reason}."
        ),
        vendor_msg(
          ">",
          "Restore it from git: {.fn use_vanilla} does not change a lock that
           it cannot read."
        )
      )
      lock <- NULL
    }
  }
  if (is.null(lock)) {
    vendor_status_print(rows, NULL)
    return(vendor_status_end(rows, error))
  }

  h <- lock$header
  release <- h[["Release"]]
  sha7 <- vendor_short_sha(h[["Commit"]])
  n <- nrow(lock$files)
  synced <- h[["Synced"]]
  header <- list(vendor_msg(
    "i",
    "inlavanilla {release}, commit {sha7}: {n} file{?s}, copied on {synced}."
  ))
  installed <- vendor_installed_version()
  installed_v <- tryCatch(numeric_version(installed), error = function(e) NULL)
  lock_v <- tryCatch(numeric_version(h[["Version"]]), error = function(e) NULL)
  if (!is.null(installed_v) && !is.null(lock_v) && installed_v < lock_v) {
    header <- c(
      header,
      list(vendor_msg(
        "!",
        "inlavanilla {installed}, installed here, is older than the copy,
       {release}: install a newer inlavanilla for checks that know that
       release."
      ))
    )
  }
  upstream <- h[["Upstream"]]
  known <- identical(upstream, vendor_upstream)
  rows$lock <- if (known) {
    vendor_row(
      "pass",
      vendor_msg("v", "The lock {.file inst/inlavanilla.lock} is well formed.")
    )
  } else {
    vendor_row(
      "warn",
      vendor_msg(
        "!",
        "The lock copies from {.url {upstream}}, which these checks do not
         know."
      )
    )
  }

  # Fetch first: the tags, a snapshot of the lock's commit, and one of the
  # newest release
  src <- NULL
  snap <- NULL
  snap_failed <- FALSE
  if (known) {
    src <- vendor_tags(clone)
  }
  if (!is.null(src)) {
    snap <- tryCatch(
      vendor_snapshot(h[["Commit"]], release, clone),
      inlavanilla_source_error = function(e) NULL
    )
    snap_failed <- is.null(snap)
    if (!is.null(snap)) {
      on.exit(unlink(snap$stage, recursive = TRUE), add = TRUE)
    }
  }

  # The copies
  problems <- vendor_check_copies(root, lock)
  rows$copies <- if (nrow(problems) == 0L) {
    vendor_row(
      "pass",
      vendor_msg(
        "v",
        "The {n} engine file{?s} match{?es/} the lock."
      )
    )
  } else {
    k <- nrow(problems)
    details <- lapply(seq_len(k), function(i) {
      f <- problems$file[i]
      switch(
        problems$problem[i],
        differs = vendor_msg(" ", "{.file {f}} differs from the lock."),
        missing = vendor_msg(" ", "{.file {f}} is in the lock but missing."),
        unlisted = vendor_msg(" ", "{.file {f}} is not in the lock.")
      )
    })
    do.call(
      vendor_row,
      c(
        list(
          "fail",
          vendor_msg(
            "x",
            "{k} engine file{?s} do{?es/} not match the lock:"
          )
        ),
        details,
        list(vendor_msg(
          ">",
          "Never edit an engine file. Restore it from git, or delete it and run
         {.code inlavanilla::use_vanilla(release = \"{release}\")}. Give your
         own files names that do not start with {.file vanilla-}."
        ))
      )
    )
  }

  # The names
  rows$names <- vendor_names_row(vendor_check_names(root, snap))

  # Dependencies and authors
  rows$dependencies <- vendor_deps_row(desc, vendor_copies(root), release)
  rows$authors <- vendor_authors_row(desc, release)

  # The license
  lic <- vendor_license_msgs(desc[["License"]])
  grade <- vendor_license_grade(desc[["License"]])
  rows$license <- do.call(
    vendor_row,
    c(
      list(switch(grade, pass = "pass", fail = "fail", "warn")),
      lic
    )
  )

  # The release and the update
  if (!known) {
    note <- vendor_msg(
      " ",
      "Skipped: the lock does not copy from {.url {vendor_upstream}}."
    )
    rows$release <- vendor_row("skip", note)
    rows$update <- vendor_row("skip", note)
  } else if (is.null(src)) {
    rows$release <- vendor_row(
      "skip",
      vendor_msg(" ", "Skipped: GitHub could not be reached.")
    )
    update_note <- list(vendor_msg(
      " ",
      "Skipped: GitHub could not be reached."
    ))
    if (!is.null(installed_v) && !is.null(lock_v) && installed_v > lock_v) {
      update_note <- c(
        update_note,
        list(vendor_msg(
          " ",
          "inlavanilla {installed} is installed here, newer than the copy,
         {release}."
        ))
      )
    }
    rows$update <- do.call(vendor_row, c(list("skip"), update_note))
  } else {
    rows$release <- vendor_release_row(lock, src, snap, snap_failed, clone)
    update <- vendor_update_row(lock, src, clone)
    if (!is.null(update$stage)) {
      on.exit(unlink(update$stage, recursive = TRUE), add = TRUE)
    }
    rows$update <- update$row
  }

  vendor_status_print(rows, header, offline = known && is.null(src))
  vendor_status_end(rows, error)
}

## ----- Rows ------------------------------------------------------------------

vendor_status_checks <- c(
  "lock",
  "copies",
  "names",
  "dependencies",
  "authors",
  "license",
  "release",
  "update"
)

# The rows, all skipped until a check runs
vendor_status_rows <- function() {
  skip <- vendor_row("skip", list(type = " ", cli = "", plain = ""))
  stats::setNames(
    rep(list(skip), length(vendor_status_checks)),
    vendor_status_checks
  )
}

vendor_status_frame <- function(rows) {
  data.frame(
    check = names(rows),
    status = vapply(rows, `[[`, "", "status", USE.NAMES = FALSE),
    message = vapply(
      rows,
      function(r) r$msgs[[1L]]$plain,
      "",
      USE.NAMES = FALSE
    ),
    details = vapply(
      rows,
      function(r) {
        paste(vapply(r$msgs[-1L], `[[`, "", "plain"), collapse = "\n")
      },
      "",
      USE.NAMES = FALSE
    )
  )
}

# Prints the header and the rows. A skipped row prints its note on an i line,
# and the two network rows skipped offline print one line together.
vendor_status_print <- function(rows, header, offline = FALSE) {
  vendor_say(header)
  for (check in vendor_status_checks) {
    row <- rows[[check]]
    if (row$status != "skip") {
      vendor_say(row$msgs)
    } else if (
      !(offline && check %in% c("release", "update")) &&
        nzchar(row$msgs[[1L]]$cli)
    ) {
      row$msgs[[1L]]$type <- "i"
      vendor_say(row$msgs)
    }
  }
  if (offline) {
    vendor_say(c(
      list(vendor_msg(
        "i",
        "Skipped 2 checks: GitHub could not be reached. Give {.arg clone}, a
         local clone of inlavanilla, to run them without the internet."
      )),
      rows$update$msgs[-1L]
    ))
  }
}

# The annotations on GitHub Actions, the error of `error = TRUE`, and the data
# frame
vendor_status_end <- function(rows, error) {
  frame <- vendor_status_frame(rows)
  if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
    for (i in which(frame$status %in% c("fail", "warn"))) {
      text <- paste(
        c(frame$message[i], frame$details[i][nzchar(frame$details[i])]),
        collapse = "\n"
      )
      text <- gsub("%", "%25", text, fixed = TRUE)
      text <- gsub("\r", "%0D", text, fixed = TRUE)
      text <- gsub("\n", "%0A", text, fixed = TRUE)
      level <- if (frame$status[i] == "fail") "error" else "warning"
      cat(paste0("::", level, "::", text, "\n"))
    }
  }
  failed <- frame$check[frame$status == "fail"]
  if (error && length(failed) > 0L) {
    n <- length(failed)
    failed <- vendor_and(failed)
    abort(
      c(
        "{n} check{?s} of the engine copy failed: {.val {failed}}.",
        i = "The lines marked with a cross say what to do."
      ),
      class = "inlavanilla_status_error"
    )
  }
  invisible(frame)
}

# The release row: the lock against its release, on GitHub or in the clone
vendor_release_row <- function(lock, src, snap, snap_failed, clone) {
  h <- lock$header
  release <- h[["Release"]]
  commit <- h[["Commit"]]
  sha7 <- vendor_short_sha(commit)
  where <- if (is.null(clone)) "on GitHub" else "in the clone"
  details <- list()
  bad <- function(template, .envir = parent.frame()) {
    details[[length(details) + 1L]] <<- vendor_msg(" ", template, .envir)
  }
  extra <- list()
  if (!release %in% names(src$tags)) {
    bad("There is no tag {release} {where}.")
  } else if (!identical(unname(src$tags[[release]]), commit)) {
    tag7 <- vendor_short_sha(src$tags[[release]])
    bad("The tag {release} is at commit {tag7}, but the lock says {sha7}.")
  }
  if (!is.null(clone) && release %in% names(src$tags)) {
    clone_commit <- unname(src$tags[[release]])
    clone7 <- vendor_short_sha(clone_commit)
    state <- vendor_github_state(release, clone_commit, src$github)
    if (state == "moved") {
      gh7 <- vendor_short_sha(src$github[[release]])
      bad("{release} is at {clone7} in the clone but at {gh7} on GitHub.")
    } else if (state == "missing") {
      bad("{release} is not on GitHub.")
    } else if (state == "offline") {
      extra <- list(vendor_msg(
        "!",
        "Could not reach GitHub to check that {release} is at {clone7} there."
      ))
    }
  }
  if (snap_failed) {
    note <- if (is.null(clone)) {
      vendor_msg(" ", "Skipped: could not download commit {sha7}.")
    } else {
      vendor_msg(" ", "Skipped: could not read commit {sha7} in the clone.")
    }
    if (length(details) == 0L) {
      return(do.call(vendor_row, c(list("skip", note), extra)))
    }
  } else {
    version <- h[["Version"]]
    upstream_v <- unname(snap$desc["Version"])
    if (!identical(version, upstream_v)) {
      bad(
        "Version is {version}, but the DESCRIPTION of inlavanilla at {sha7}
         says {upstream_v}."
      )
    }
    for (i in seq_len(nrow(lock$files))) {
      f <- lock$files$file[i]
      if (!f %in% snap$files) {
        bad("{.file {f}} is not in inlavanilla at {sha7}.")
      } else if (vendor_md5(file.path(snap$dir, f)) != lock$files$md5[i]) {
        bad("{.file {f}} has another MD5 at {release} than in the lock.")
      }
    }
  }
  if (length(details) == 0L) {
    line <- if (is.null(clone)) {
      vendor_msg("v", "The lock matches release {release} on GitHub.")
    } else {
      vendor_msg(
        "v",
        "The lock matches release {release} in the clone {.path {clone}}."
      )
    }
    return(do.call(vendor_row, c(list("pass", line), extra)))
  }
  line <- if (is.null(clone)) {
    vendor_msg("x", "The lock does not match release {release} on GitHub:")
  } else {
    vendor_msg(
      "x",
      "The lock does not match release {release} in the clone {.path {clone}}:"
    )
  }
  fix <- vendor_msg(
    ">",
    "Copy the release again with
     {.code inlavanilla::use_vanilla(release = \"{release}\")}."
  )
  do.call(vendor_row, c(list("fail", line), details, list(fix), extra))
}

# The update row: whether the newest release changes a copied file. Only the
# newest release is judged, so a change undone between releases does not
# count. Never a failure. Returns list(row, stage).
vendor_update_row <- function(lock, src, clone) {
  h <- lock$header
  release <- h[["Release"]]
  releases <- vendor_release_tags(src$tags)
  newer <- vendor_newer_tags(names(releases), release)
  if (length(newer) == 0L) {
    return(list(
      row = vendor_row(
        "pass",
        vendor_msg("v", "No release of inlavanilla is newer than {release}.")
      )
    ))
  }
  newest <- newer[length(newer)]
  if (identical(unname(releases[[newest]]), h[["Commit"]])) {
    return(list(
      row = vendor_row(
        "pass",
        vendor_msg(
          "v",
          "{newest}, the newest release, is the same commit as {release}."
        )
      )
    ))
  }
  snap <- tryCatch(
    vendor_snapshot(releases[[newest]], newest, clone),
    inlavanilla_source_error = function(e) NULL
  )
  if (is.null(snap)) {
    return(list(
      row = vendor_row(
        "skip",
        vendor_msg(" ", "Skipped: could not get {newest}, the newest release.")
      )
    ))
  }
  files <- lock$files$file
  gone <- setdiff(files, snap$files)
  kept <- intersect(files, snap$files)
  md5 <- vendor_md5(file.path(snap$dir, kept))
  changed <- kept[md5 != lock$files$md5[match(kept, files)]]
  not_copied <- setdiff(snap$files, files)
  details <- list()
  if (length(gone) > 0L) {
    gone <- vendor_and(gone)
    details <- c(
      details,
      list(vendor_msg(
        " ",
        "{newest} no longer has {.file {gone}}."
      ))
    )
  }
  if (length(not_copied) > 0L) {
    not_copied <- vendor_and(not_copied)
    details <- c(
      details,
      list(vendor_msg(
        " ",
        "{newest} also has {.file {not_copied}}, which {?is/are} not copied."
      ))
    )
  }
  if (length(changed) == 0L && length(gone) == 0L) {
    row <- do.call(
      vendor_row,
      c(
        list("pass", vendor_msg("v", "{newest} changes no copied file.")),
        details
      )
    )
    return(list(row = row, stage = snap$stage))
  }
  compare <- paste0(vendor_upstream, "/compare/", release, "...", newest)
  line <- if (length(changed) > 0L) {
    changed <- vendor_and(changed)
    vendor_msg(
      "!",
      "A newer release changes copied files: {newest} changes
       {.file {changed}}."
    )
  } else {
    vendor_msg(
      "!",
      "A newer release changes the set of engine files: {newest}."
    )
  }
  row <- do.call(
    vendor_row,
    c(
      list("warn", line),
      details,
      list(
        vendor_msg(" ", "Changes: {.url {compare}}"),
        vendor_msg(
          ">",
          "Update with {.code inlavanilla::use_vanilla(release = \"{newest}\")},
         then run your tests."
        )
      )
    )
  )
  list(row = row, stage = snap$stage)
}
