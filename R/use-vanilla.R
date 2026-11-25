#' Copy the engine of inlavanilla into a package
#'
#' @description
#' Copies the engine of inlavanilla into the package whose root is `path`, so
#' that the package can fit models with `vanilla_core()` and does not need
#' inlavanilla to run. The engine is the code behind [inla_vanilla()]: the
#' files `R/vanilla-*.R` of an inlavanilla release. They are copied byte for
#' byte, under the same names, into the package's `R/` folder, and listed with
#' their MD5 hashes in the lock `inst/inlavanilla.lock`.
#'
#' `use_vanilla()` also adds to the package's `DESCRIPTION` what the engine
#' needs, where it is missing:
#' * the copyright holders of the engine in `Authors@R`, with the role
#'   `"cph"`, as the `DESCRIPTION` of inlavanilla at that release lists them;
#' * `cli (>= 3.6.0)`, `stats` and `utils` in `Imports`;
#' * `future` and `parallelly` in `Suggests`, for the optional parallel scans;
#' * `R (>= 4.1.0)` in `Depends`, where the package asks for an older R or for
#'   none.
#'
#' It never removes anything, and changes no other file of the package. It
#' prints each change, and ends with a list of what is still to do: any name
#' clashes, the function that calls `vanilla_core()`, the documentation, the
#' tests, R CMD check and a NEWS bullet.
#'
#' Run it again with a newer `release` to update the copy.
#'
#' @details
#' **Where the files come from.** By default from the release on GitHub,
#' <https://github.com/haziqj/inlavanilla>, which needs the internet but not
#' git. With `clone`, from the release tag in a local clone of inlavanilla,
#' which needs git but not the internet. The lock records the release, its
#' commit and the hash of each file, so that [vanilla_status()] can check the
#' copy later.
#'
#' **Which files.** A first copy takes every engine file of the release, or the
#' files in `files`. An update copies the files in the lock, and adds the files
#' in `files`, so that every copy comes from one release. It reports the
#' engine files of the release that are not copied, and the ones that the
#' copies need.
#'
#' **A copy edited by hand.** When an engine file of the package differs from
#' the lock, `use_vanilla()` stops and changes nothing, so that the edit is not
#' lost. It also stops for a file `R/vanilla-*.R` that is not in the lock and
#' differs from the engine file of that name. Never edit an engine file: the
#' next update would undo the edit. Report the problem at
#' <https://github.com/haziqj/inlavanilla/issues> instead. The names
#' `R/vanilla-*.R` are kept for the engine, so give your own files other names,
#' and leave the engine files out of any code formatter or linter.
#'
#' **A rights holder already listed.** A copyright holder of the engine who is
#' already in `Authors@R` is not added again. When that entry lacks the role
#' `"cph"`, `use_vanilla()` leaves it as it is, and lists the role in what is
#' still to do.
#'
#' **License.** The engine is under the GNU General Public License, version 3
#' or later, so a package that holds it must be too. `use_vanilla()` warns when
#' the package's `License` does not allow that.
#'
#' `vignette("engine", package = "inlavanilla")`, also at
#' <https://haziqj.ml/inlavanilla/articles/engine.html>, describes the
#' interface of `vanilla_core()`, with an example.
#'
#' @param release The release of inlavanilla to copy: its tag on GitHub, such
#'   as `"v26.11.24"`, or its version, `"26.11.24"`. By default the version of
#'   inlavanilla installed here, or, when that version is not a release, the
#'   newest release before it.
#' @param files The engine files to copy, as `"vanilla-vb.R"` or
#'   `"R/vanilla-vb.R"`. By default every engine file of the release on a first
#'   copy, and the files in the lock on an update (see Details). The engine
#'   files call each other, so copy only some of them when you know which you
#'   need.
#' @param clone The path of a local clone of inlavanilla, to copy from without
#'   the internet. It needs git, and the release tag must be in the clone
#'   (`git fetch --tags`). By default the files come from GitHub.
#' @param path The root of the package, the folder that holds its
#'   `DESCRIPTION`, or a folder inside it. By default the working directory.
#'
#' @returns Invisibly, a data frame with a row per engine file in the lock:
#'   `file`, its `md5` hash, and its `state`, `"new"`, `"restored"`,
#'   `"changed"` or `"unchanged"`. Its attributes `release`, `commit` and
#'   `previous` give the release copied, its commit and the release before it
#'   (`NA` for a first copy), and `description` and `open` the changes made to
#'   `DESCRIPTION` and what is still to do, as plain text.
#'
#' @seealso [vanilla_status()] to check the copy, and [use_vanilla_test()] and
#'   [use_vanilla_ci()] to check it in the tests and on GitHub.
#'
#' @examples
#' \dontrun{
#' # From the root of your package: copy the engine from the release of the
#' # inlavanilla installed here
#' use_vanilla()
#'
#' # Update the copy to a later release, for example
#' use_vanilla("v26.12.01")
#'
#' # Without the internet, from a clone of inlavanilla
#' use_vanilla(clone = "~/src/inlavanilla")
#' }
#' @export
use_vanilla <- function(
  release = NULL,
  files = NULL,
  clone = NULL,
  path = "."
) {
  error_class <- "inlavanilla_vendor_error"

  ## ----- Arguments and the package -------------------------------------------
  ok <- is.null(release) ||
    (is.character(release) &&
      length(release) == 1L &&
      !is.na(release) &&
      nzchar(release))
  if (!ok) {
    abort(
      "{.arg release} must be a single string, such as {.val v26.11.24}, or
       {.code NULL}.",
      class = error_class
    )
  }
  ok <- is.null(files) ||
    (is.character(files) &&
      length(files) > 0L &&
      !anyNA(files) &&
      all(nzchar(files)))
  if (!ok) {
    abort(
      "{.arg files} must be the names of engine files, such as
       {.val vanilla-vb.R}, or {.code NULL}.",
      class = error_class
    )
  }
  clone <- vendor_check_clone(clone)
  root <- vendor_root(path)
  desc <- vendor_desc(root)
  pkg <- vendor_package(root, desc)
  if (identical(pkg, "inlavanilla")) {
    abort(
      "{.path {root}} is inlavanilla itself. {.fn use_vanilla} copies its
       engine into another package.",
      class = error_class
    )
  }
  vendor_start()

  ## ----- The old lock and the copies -----------------------------------------
  lock_path <- file.path(root, vendor_lock_file)
  lock <- NULL
  if (file.exists(lock_path)) {
    lock <- tryCatch(
      vendor_read_lock(lock_path),
      inlavanilla_lock_error = function(e) {
        reason <- e$reason
        abort(
          "{.file inst/inlavanilla.lock} is not a lock that {.fn use_vanilla}
           can read: {reason}.",
          class = error_class
        )
      }
    )
    upstream <- lock$header[["Upstream"]]
    if (!identical(upstream, vendor_upstream)) {
      abort(
        "The lock copies from {.url {upstream}}, not from
         {.url {vendor_upstream}}.",
        class = error_class
      )
    }
    problems <- vendor_check_copies(root, lock)
    differ <- problems$file[problems$problem == "differs"]
    if (length(differ) > 0L) {
      n <- length(differ)
      differ <- vendor_and(differ)
      issues <- vendor_issues_url
      them <- if (n == 1L) "it" else "them"
      abort(
        c(
          "{n} engine file{?s} differ{?s/} from the lock: {.file {differ}}.",
          i = "{.fn use_vanilla} would overwrite your changes, so it changed
               nothing. Never edit an engine file: report the problem at
               {.url {issues}}. To discard the changes, restore {them} from
               git, or delete {them}, then run {.fn use_vanilla} again."
        ),
        class = error_class
      )
    }
  }
  previous <- if (is.null(lock)) NA_character_ else lock$header[["Release"]]
  unlisted <- setdiff(vendor_copies(root), lock$files$file)

  ## ----- The release ---------------------------------------------------------
  src <- vendor_tags(clone)
  if (is.null(src)) {
    abort(
      c(
        "Could not reach GitHub to read the releases of inlavanilla.",
        i = "Check the connection, or copy from a local clone of inlavanilla
             with {.arg clone}."
      ),
      class = error_class
    )
  }
  res <- vendor_resolve_release(
    release,
    src$tags,
    vendor_installed_version(),
    lock,
    clone
  )
  release <- res$release
  commit <- res$commit
  sha7 <- vendor_short_sha(commit)
  github_msgs <- list()
  if (!is.null(clone)) {
    state <- vendor_github_state(release, commit, src$github)
    if (state == "moved") {
      gh7 <- vendor_short_sha(src$github[[release]])
      abort(
        "{release} is at {sha7} in the clone but at {gh7} on GitHub.",
        class = error_class
      )
    }
    if (state == "missing") {
      abort(
        c(
          "{release} is not on GitHub.",
          i = "The lock names a release on GitHub, so copy a release that is
               there."
        ),
        class = error_class
      )
    }
    if (state == "offline") {
      github_msgs <- list(vendor_msg(
        "!",
        "Could not reach GitHub to check that {release} is at {sha7} there."
      ))
    }
  }

  ## ----- The files of the release --------------------------------------------
  snap <- vendor_snapshot(commit, release, clone)
  on.exit(unlink(snap$stage, recursive = TRUE), add = TRUE)
  if (length(snap$files) == 0L) {
    abort(
      "{release} has no engine files {.file R/vanilla-*.R}.",
      class = error_class
    )
  }
  want <- NULL
  if (!is.null(files)) {
    want <- ifelse(startsWith(files, "R/"), files, file.path("R", files))
    unknown <- setdiff(want, snap$files)
    if (length(unknown) > 0L) {
      unknown <- vendor_and(unknown)
      all_files <- vendor_and(snap$files)
      abort(
        c(
          "Not an engine file of {release}: {.file {unknown}}.",
          i = "Its engine files are {.file {all_files}}."
        ),
        class = error_class
      )
    }
  }
  set <- if (is.null(lock)) {
    if (is.null(want)) snap$files else unique(want)
  } else {
    union(lock$files$file, want)
  }
  gone <- setdiff(set, snap$files)
  still <- gone[file.exists(file.path(root, gone))]
  if (length(still) > 0L) {
    still <- vendor_and(still)
    abort(
      c(
        "{release} no longer has {.file {still}}.",
        i = "Delete {.file {still}}, which match{?es/} the lock and so
             hold{?s/} no edits of yours, then run {.fn use_vanilla} again."
      ),
      class = error_class
    )
  }
  set <- sort(setdiff(set, gone), method = "radix")
  for (f in unlisted) {
    same <- f %in%
      set &&
      identical(
        vendor_md5(file.path(root, f)),
        vendor_md5(file.path(snap$dir, f))
      )
    if (!same) {
      abort(
        c(
          "{.file {f}} has the name of an engine file but is not in the lock,
           so it may hold your own code.",
          i = "Rename it, or delete it if it is an older copy of the engine,
               then run {.fn use_vanilla} again."
        ),
        class = error_class
      )
    }
  }
  for (f in set) {
    lines <- readLines(file.path(snap$dir, f), warn = FALSE)
    tags <- trimws(grep(vendor_roxygen_ns, lines, value = TRUE))
    if (length(tags) > 0L) {
      tags <- vendor_and(tags)
      abort(
        "{.file {f}} at {release} has roxygen tags that would change
         {.file NAMESPACE}: {.code {tags}}.",
        class = error_class
      )
    }
  }

  ## ----- The plan ------------------------------------------------------------
  md5 <- vendor_md5(file.path(snap$dir, set))
  old_md5 <- lock$files$md5[match(set, lock$files$file)]
  if (is.null(old_md5)) {
    old_md5 <- rep(NA_character_, length(set))
  }
  state <- ifelse(
    is.na(old_md5),
    "new",
    ifelse(
      !file.exists(file.path(root, set)),
      "restored",
      ifelse(md5 != old_md5, "changed", "unchanged")
    )
  )
  holders <- vendor_holders(unname(snap$desc["Authors@R"]))
  plan <- vendor_plan_description(root, holders)

  ## ----- Write ---------------------------------------------------------------
  copy <- set[state != "unchanged"]
  if (length(copy) > 0L) {
    r_dir <- file.path(root, "R")
    dir.create(r_dir, showWarnings = FALSE)
    copied <- file.copy(
      file.path(snap$dir, copy),
      file.path(root, copy),
      overwrite = TRUE
    )
    if (!all(copied)) {
      abort(
        "Could not copy the engine files into {.path {r_dir}}.",
        class = error_class
      )
    }
  }
  same_lock <- !is.null(lock) &&
    identical(lock$header[["Release"]], release) &&
    identical(lock$header[["Commit"]], commit) &&
    setequal(lock$files$file, set) &&
    all(state == "unchanged")
  if (!same_lock) {
    header <- c(
      Upstream = vendor_upstream,
      Release = release,
      Commit = commit,
      Version = unname(snap$desc["Version"]),
      Synced = format(vendor_today()),
      Hash = vendor_hash
    )
    vendor_write_lock(lock_path, header, data.frame(file = set, md5 = md5))
  }
  desc_written <- length(plan$changes) > 0L
  if (desc_written) {
    vendor_desc_write(plan$d, file.path(root, "DESCRIPTION"))
  }
  written <- length(copy) > 0L || !same_lock || desc_written

  ## ----- After ---------------------------------------------------------------
  new_lock <- vendor_read_lock(lock_path)
  post <- vendor_check_copies(root, new_lock)
  names <- vendor_check_names(root, snap)
  missing <- names[names$kind == "missing", , drop = FALSE]
  names <- names[names$kind != "missing", , drop = FALSE]
  desc <- vendor_desc(root)
  grade <- vendor_license_grade(desc[["License"]])
  collate <- length(vendor_collate_lacks(desc, set)) > 0L
  releases <- vendor_release_tags(src$tags)
  newest <- names(releases)[length(releases)]

  ## ----- The report ----------------------------------------------------------
  first <- is.null(lock)
  same_release <- !first && identical(previous, release)
  if (first) {
    cli::cli_h2("Copying the engine of inlavanilla into {pkg}")
  } else if (same_release) {
    cli::cli_h2("The engine of inlavanilla in {pkg}")
  } else {
    cli::cli_h2("Updating the engine of inlavanilla in {pkg}")
  }
  done <- list()
  say <- function(type, template, .envir = parent.frame()) {
    done[[length(done) + 1L]] <<- vendor_msg(type, template, .envir)
  }
  from <- if (is.null(clone)) {
    "from GitHub."
  } else {
    "from the clone {.path {clone}}."
  }
  if (first || same_release) {
    say("i", paste("Release {release}, commit {sha7},", from))
  } else {
    say(
      "i",
      paste(
        "Release {release}, commit {sha7},",
        from,
        "The lock was at {previous}."
      )
    )
  }
  done <- c(done, res$notes)
  n <- length(set)
  n_changed <- sum(state %in% c("changed", "restored"))
  n_unchanged <- sum(state == "unchanged")
  n_new <- sum(state == "new")
  if (first) {
    say("v", "Copied {n} engine file{?s} to {.file R/}.")
  } else if (n_unchanged == n) {
    say("v", "The {n} engine file{?s} {?is/are} those of {release} already.")
  } else {
    parts <- c(
      if (n_changed > 0L) paste(n_changed, "changed"),
      if (n_unchanged > 0L) paste(n_unchanged, "unchanged"),
      if (n_new > 0L) paste(n_new, "new")
    )
    parts <- paste(parts, collapse = ", ")
    say("v", "Copied {n} engine file{?s} to {.file R/}: {parts}.")
    changed <- set[state %in% c("changed", "restored")]
    if (length(changed) > 0L) {
      changed <- vendor_and(changed)
      say(" ", "Changed: {.file {changed}}.")
    }
  }
  if (!same_lock) {
    if (first) {
      say("v", "Wrote the lock {.file inst/inlavanilla.lock}.")
    } else {
      say("v", "Updated the lock {.file inst/inlavanilla.lock}.")
    }
  }
  for (f in gone) {
    say(
      "v",
      "Dropped {.file {f}} from the lock: {release} no longer has it, and the
       copy was already deleted."
    )
  }
  done <- c(done, vendor_change_msgs(plan$changes))
  if (!isTRUE(plan$safe) && length(plan$changes) == 0L) {
    say(
      "!",
      "{.file DESCRIPTION} could not be edited safely, so it was left as it
       was: make the changes below by hand."
    )
  } else if (identical(plan$left, "Authors@R")) {
    say(
      "!",
      "{.field Authors@R} could not be edited safely, so it was left as it
       was: add the copyright holders by hand, as below."
    )
  } else if (identical(plan$left, "dependencies")) {
    say(
      "!",
      "The dependencies in {.file DESCRIPTION} could not be edited safely, so
       they were left as they were: make the changes below by hand."
    )
  } else if (length(plan$changes) == 0L && length(plan$open) == 0L) {
    say("v", "{.file DESCRIPTION} already has what the engine needs.")
  }
  if (nrow(post) > 0L) {
    issues <- vendor_issues_url
    say(
      "x",
      "The copies do not match the new lock, which is a bug: report it at
       {.url {issues}}."
    )
  }
  names_row <- vendor_names_row(names)
  clashes <- names_row$status == "fail"
  names_msgs <- Filter(function(m) m$type != ">", names_row$msgs)
  done <- c(done, names_msgs)
  uncopied <- setdiff(snap$files, set)
  if (length(uncopied) > 0L) {
    k <- length(uncopied)
    uncopied <- vendor_and(uncopied)
    say(
      if (nrow(missing) > 0L) "!" else "i",
      "{release} has {k} engine file{?s} that {?is/are} not copied:
       {.file {uncopied}}."
    )
    for (i in seq_len(nrow(missing))) {
      done[[length(done) + 1L]] <- list(
        type = " ",
        cli = vendor_names_item(missing[i, ]),
        plain = missing$text[i]
      )
    }
  }
  if (grade != "pass") {
    done <- c(done, vendor_license_msgs(desc[["License"]], pkg)[1L])
  }
  done <- c(done, github_msgs)
  if (
    !is.na(newest) && vendor_tag_version(newest) > vendor_tag_version(release)
  ) {
    say(
      "i",
      "{newest}, a newer release, exists:
       {.run inlavanilla::use_vanilla(\"{newest}\")} updates to it."
    )
  }
  vendor_say(done)

  todo <- list()
  open <- list()
  add_todo <- function(template, .envir = parent.frame(), is_open = FALSE) {
    msg <- vendor_msg("_", template, .envir)
    todo[[length(todo) + 1L]] <<- msg
    if (is_open) {
      open[[length(open) + 1L]] <<- msg
    }
  }
  if (clashes) {
    add_todo(
      "Rename or remove your own objects named above. Never rename an engine
       object. Then run {.run inlavanilla::vanilla_status()}.",
      is_open = TRUE
    )
  }
  if (nrow(missing) > 0L) {
    need <- basename(missing$name)
    files_code <- if (length(need) == 1L) {
      paste0("\"", need, "\"")
    } else {
      paste0("c(", paste0("\"", need, "\"", collapse = ", "), ")")
    }
    add_todo(
      "Copy the engine files that the copies need:
       {.run inlavanilla::use_vanilla(release = \"{release}\", files =
       {files_code})}.",
      is_open = TRUE
    )
  }
  if (collate) {
    add_todo(
      "{.file DESCRIPTION} has a {.field Collate} field, which must list every
       file in {.file R/}: add the engine files to it.",
      is_open = TRUE
    )
  }
  todo <- c(todo, plan$open)
  open <- c(open, plan$open)
  if (grade != "pass") {
    vignette <- "vignette(\"engine\", package = \"inlavanilla\")"
    add_todo(
      "License {pkg} under GPL version 3 or later before you release it. See
       \"Licensing\" in {.code {vignette}}.",
      is_open = TRUE
    )
  }
  formatter <- c("air.toml", ".air.toml", ".lintr", ".styler")
  if (written && any(file.exists(file.path(root, formatter)))) {
    add_todo(
      "Leave {.file R/vanilla-*.R} out of your formatter and linter: the
       copies must stay byte for byte."
    )
  }
  if (first) {
    vignette <- "vignette(\"engine\", package = \"inlavanilla\")"
    url <- vendor_vignette_url
    add_todo(
      "Write the function that calls {.fn vanilla_core} with your model's
       log-posterior and turns the list it returns into your results.
       {.code {vignette}} shows how, and so does {.url {url}}."
    )
  }
  if (!first && !same_release) {
    news <- paste0(vendor_upstream, "/blob/", release, "/NEWS.md")
    compare <- paste0(vendor_upstream, "/compare/", previous, "...", release)
    add_todo(
      "Read what changed since {previous}: the NEWS of inlavanilla at
       {.url {news}}, and the engine's code at {.url {compare}}."
    )
  }
  authors_changed <- any(vapply(
    plan$changes,
    function(ch) identical(ch$kind, "authors"),
    NA
  ))
  if (first || clashes || authors_changed) {
    add_todo("Run {.run devtools::document()}.")
  }
  if (written) {
    add_todo("Run your tests and R CMD check: {.run devtools::check()}.")
    add_todo(
      "Add a NEWS bullet that names the release, such as \"Uses the engine of
       inlavanilla {release}.\""
    )
    if (file.exists(file.path(root, ".git"))) {
      add_todo(
        "Commit {.file R/vanilla-*.R}, {.file inst/inlavanilla.lock} and
         {.file DESCRIPTION} together."
      )
    }
  }
  if (length(todo) > 0L) {
    cli::cli_h2("Still to do")
    vendor_say(todo)
  }

  optional <- list()
  if (!file.exists(file.path(root, vendor_test_file))) {
    optional <- c(
      optional,
      list(vendor_msg(
        "*",
        "{.run inlavanilla::use_vanilla_test()} writes a test that fails when an
       engine file differs from the lock."
      ))
    )
  }
  if (!file.exists(file.path(root, vendor_ci_file))) {
    optional <- c(
      optional,
      list(vendor_msg(
        "*",
        "{.run inlavanilla::use_vanilla_ci()} adds a GitHub Actions workflow that
       runs {.fn vanilla_status} on every push."
      ))
    )
  }
  if (length(optional) > 0L) {
    cli::cli_h2("Optional")
    vendor_say(optional)
  }

  out <- data.frame(file = set, md5 = md5, state = state)
  attr(out, "release") <- release
  attr(out, "commit") <- commit
  attr(out, "previous") <- previous
  attr(out, "description") <- vapply(
    vendor_change_msgs(plan$changes),
    `[[`,
    "",
    "plain"
  )
  attr(out, "open") <- vapply(open, `[[`, "", "plain")
  invisible(out)
}

#' Add a test that checks the engine copy
#'
#' @description
#' Writes the test file `tests/testthat/test-vanilla-engine.R` in the package
#' whose root is `path`. The test fails when an engine file differs from the
#' lock `inst/inlavanilla.lock`, when a file in the lock is missing, or when a
#' file `R/vanilla-*.R` is not in the lock, so that `devtools::test()` catches
#' an edit to a copy. It uses base R and testthat only: the package does not
#' need inlavanilla to run it.
#'
#' The test reads the source files in `R/`, which an installed package does not
#' hold. So it runs with `devtools::test()` from the source tree, and skips
#' under R CMD check.
#'
#' It is optional: [use_vanilla()] does not add it, and offers it at the end of
#' its list. The package must already use testthat, as
#' `usethis::use_testthat()` sets up.
#'
#' @param path The root of the package, or a folder inside it. By default the
#'   working directory.
#' @param overwrite If `TRUE`, replace a test file of that name that differs
#'   from the one this version of inlavanilla writes. By default such a file is
#'   left alone, with a message.
#'
#' @returns Invisibly, the path of the test file.
#'
#' @seealso [vanilla_status()], which also checks the lock against its release,
#'   the names and `DESCRIPTION`.
#'
#' @examples
#' \dontrun{
#' # From the root of a package that holds a copy of the engine
#' use_vanilla_test()
#' }
#' @export
use_vanilla_test <- function(path = ".", overwrite = FALSE) {
  check_flag(overwrite, "overwrite")
  root <- vendor_root(path)
  desc <- vendor_desc(root)
  pkg <- vendor_package(root, desc)
  vendor_start()
  if (!dir.exists(file.path(root, "tests", "testthat"))) {
    abort(
      c(
        "{pkg} has no testthat tests ({.file tests/testthat/}).",
        i = "Set them up with {.code usethis::use_testthat(3)}, then run
             {.fn use_vanilla_test} again."
      ),
      class = "inlavanilla_vendor_error"
    )
  }
  if (!file.exists(file.path(root, vendor_lock_file))) {
    abort(
      "There is no lock {.file inst/inlavanilla.lock}: copy the engine with
       {.fn use_vanilla} first.",
      class = "inlavanilla_vendor_error"
    )
  }
  dest <- vendor_template(
    root,
    "test-vanilla-engine.R",
    vendor_test_file,
    overwrite
  )
  invisible(dest)
}

#' Check the engine copy on GitHub Actions
#'
#' @description
#' Writes the GitHub Actions workflow `.github/workflows/vanilla-status.yaml`
#' in the package whose root is `path`. On every push to `main` or `master`,
#' and on every pull request, it installs the latest release of inlavanilla
#' from GitHub and runs `vanilla_status(error = TRUE)`, which fails the run
#' when a check fails. A newer release of inlavanilla is reported in the log,
#' and does not fail the run.
#'
#' It also adds `^\.github$` to `.Rbuildignore` where that line is missing,
#' so that R CMD build leaves the folder out.
#'
#' It is optional: [use_vanilla()] does not add it, and offers it at the end of
#' its list.
#'
#' @param path The root of the package, or a folder inside it. By default the
#'   working directory.
#' @param overwrite If `TRUE`, replace a workflow file of that name that
#'   differs from the one this version of inlavanilla writes. By default such
#'   a file is left alone, with a message.
#'
#' @returns Invisibly, the path of the workflow file.
#'
#' @seealso [vanilla_status()] for the checks.
#'
#' @examples
#' \dontrun{
#' use_vanilla_ci()
#' }
#' @export
use_vanilla_ci <- function(path = ".", overwrite = FALSE) {
  check_flag(overwrite, "overwrite")
  root <- vendor_root(path)
  desc <- vendor_desc(root)
  vendor_package(root, desc)
  vendor_start()
  dest <- vendor_template(
    root,
    "vanilla-status.yaml",
    vendor_ci_file,
    overwrite
  )
  ignore <- file.path(root, ".Rbuildignore")
  line <- "^\\.github$"
  lines <- if (file.exists(ignore)) {
    readLines(ignore, warn = FALSE)
  } else {
    character()
  }
  if (!line %in% lines) {
    if (file.exists(ignore)) {
      d <- vendor_desc_read(ignore)
      d$lines <- c(d$lines, line)
      d$final_nl <- TRUE
    } else {
      d <- list(lines = line, eol = "\n", final_nl = TRUE)
    }
    vendor_desc_write(d, ignore)
    vendor_say(list(vendor_msg(
      "v",
      "Added {.code {line}} to {.file .Rbuildignore}."
    )))
  }
  if (!file.exists(file.path(root, vendor_lock_file))) {
    vendor_say(list(vendor_msg(
      "!",
      "There is no lock yet: the workflow fails until {.fn use_vanilla} has
       copied the engine."
    )))
  }
  invisible(dest)
}

# Writes a template of inst/templates to `target` in the package, unless a file
# that differs is there and `overwrite` is FALSE. Carriage returns do not
# count, so that a checkout with CRLF line endings is up to date. Returns the
# path.
vendor_template <- function(root, template, target, overwrite) {
  src <- system.file(
    "templates",
    template,
    package = "inlavanilla",
    mustWork = TRUE
  )
  bytes <- readBin(src, "raw", n = file.size(src))
  dest <- file.path(root, target)
  shown <- gsub("\\", "/", target, fixed = TRUE)
  if (!file.exists(dest)) {
    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    writeBin(bytes, dest)
    vendor_say(list(vendor_msg("v", "Wrote {.file {shown}}.")))
  } else if (identical(vendor_lf(dest), bytes[bytes != as.raw(13L)])) {
    vendor_say(list(vendor_msg("v", "{.file {shown}} is up to date.")))
  } else if (overwrite) {
    writeBin(bytes, dest)
    vendor_say(list(vendor_msg("v", "Replaced {.file {shown}}.")))
  } else {
    vendor_say(list(vendor_msg(
      "!",
      "{.file {shown}} differs from the file this version of inlavanilla
       writes, so it was left alone. Give {.code overwrite = TRUE} to replace
       it."
    )))
  }
  dest
}

# The bytes of a file with its carriage returns removed
vendor_lf <- function(path) {
  bytes <- readBin(path, "raw", n = file.size(path))
  bytes[bytes != as.raw(13L)]
}
