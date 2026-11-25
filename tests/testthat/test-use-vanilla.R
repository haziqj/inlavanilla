# use_vanilla(), use_vanilla_test() and use_vanilla_ci() on host packages, with
# a fake GitHub

quiet <- function(expr) suppressMessages(expr)

test_that("a first copy takes every engine file of the release", {
  gh <- local_vendor()
  root <- local_host()
  namespace <- read_text(file.path(root, "NAMESPACE"))
  expect_snapshot(out <- use_vanilla(path = root))
  release <- fake_releases()[["v26.11.24-1"]]
  for (f in c("R/vanilla-a.R", "R/vanilla-b.R", "R/vanilla-c.R")) {
    expect_identical(read_text(file.path(root, f)), release[[f]])
  }
  expect_false(file.exists(file.path(root, "R", "methods.R")))
  lock <- vendor_read_lock(file.path(root, "inst", "inlavanilla.lock"))
  expect_identical(
    lock$header,
    c(
      Upstream = "https://github.com/haziqj/inlavanilla",
      Release = "v26.11.24-1",
      Commit = gh$commits[["v26.11.24-1"]],
      Version = "26.11.24-1",
      Synced = "2026-11-25",
      Hash = "md5 of the text with LF line endings"
    )
  )
  expect_identical(nrow(vendor_check_copies(root, lock)), 0L)
  expect_identical(read_text(file.path(root, "NAMESPACE")), namespace)
  expect_identical(read_text(file.path(root, "DESCRIPTION")), minimal_expected)
  expect_identical(out$file, lock$files$file)
  expect_identical(out$md5, lock$files$md5)
  expect_identical(out$state, rep("new", 3))
  expect_identical(attr(out, "release"), "v26.11.24-1")
  expect_identical(attr(out, "commit"), gh$commits[["v26.11.24-1"]])
  expect_identical(attr(out, "previous"), NA_character_)
  expect_length(attr(out, "description"), 4L)
  expect_identical(attr(out, "open"), character())
  # Each URL is fetched once
  expect_identical(anyDuplicated(gh$urls()), 0L)
})

test_that("the default release falls back to the newest before the installed", {
  local_vendor(installed = "26.11.25")
  root <- local_host()
  expect_snapshot(out <- use_vanilla(path = root))
  expect_identical(attr(out, "release"), "v26.11.24-1")
})

test_that("a subset of files is copied, and the set grows with files", {
  local_vendor()
  root <- local_host()
  expect_snapshot(out <- use_vanilla(files = "vanilla-a.R", path = root))
  expect_identical(out$file, "R/vanilla-a.R")
  expect_match(
    attr(out, "open"),
    "use_vanilla(release = \"v26.11.24-1\", files = \"vanilla-b.R\")",
    fixed = TRUE
  )
  out <- quiet(use_vanilla(files = "R/vanilla-b.R", path = root))
  expect_identical(out$file, c("R/vanilla-a.R", "R/vanilla-b.R"))
  expect_identical(out$state, c("unchanged", "new"))
  expect_identical(attr(out, "open"), character())
  expect_false(file.exists(file.path(root, "R", "vanilla-c.R")))
})

test_that("the same release again changes nothing", {
  local_vendor()
  root <- local_host()
  quiet(use_vanilla(path = root))
  before <- folder_state(root)
  local_vendor(today = "2026-12-31")
  expect_snapshot(out <- use_vanilla(path = root))
  expect_identical(folder_state(root), before)
  expect_identical(out$state, rep("unchanged", 3))
})

test_that("an update copies the changed files and rewrites the lock", {
  gh <- local_vendor(installed = "26.11.24")
  root <- local_host()
  quiet(use_vanilla(path = root))
  local_vendor(gh, installed = "26.11.24-1", today = "2026-11-26")
  expect_snapshot(out <- use_vanilla(path = root))
  expect_identical(out$state, c("unchanged", "changed", "unchanged"))
  expect_identical(attr(out, "previous"), "v26.11.24")
  lock <- vendor_read_lock(file.path(root, "inst", "inlavanilla.lock"))
  expect_identical(lock$header[["Release"]], "v26.11.24-1")
  expect_identical(lock$header[["Synced"]], "2026-11-26")
  expect_identical(
    read_text(file.path(root, "R", "vanilla-b.R")),
    engine_b("26.11.24-1")
  )
})

test_that("refusals change nothing in the host", {
  refused <- function(root, pattern, ...) {
    before <- folder_state(root)
    expect_error(
      quiet(use_vanilla(..., path = root)),
      pattern,
      class = "inlavanilla_vendor_error"
    )
    expect_identical(folder_state(root), before)
  }
  gh <- local_vendor()
  root <- local_host()
  quiet(use_vanilla(path = root))
  # An edited copy
  write_text(file.path(root, "R", "vanilla-a.R"), paste0(engine_a, "# mine\n"))
  refused(root, "1 engine file differs from the lock")
  write_text(file.path(root, "R", "vanilla-a.R"), engine_a)
  # A malformed lock, and a lock of another upstream
  lock_path <- file.path(root, "inst", "inlavanilla.lock")
  lock_text <- read_text(lock_path)
  write_text(lock_path, sub("Commit: c", "Commit: x", lock_text))
  refused(root, "is not a lock that")
  write_text(lock_path, sub("haziqj/inlavanilla", "someone/fork", lock_text))
  refused(root, "The lock copies from")
  write_text(lock_path, lock_text)
  # Releases that cannot be copied
  refused(root, "no release \"v26.12.01\" on GitHub", release = "v26.12.01")
  refused(root, "is not a release of inlavanilla", release = "main")
  refused(root, "older than v26.11.24", release = "v26.11.23")
  refused(root, "Not an engine file of v26.11.24-1", files = "vanilla-z.R")
  # Offline
  local_vendor(fake_github(offline = TRUE))
  refused(root, "Could not reach GitHub")
  # An unlisted file R/vanilla-*.R that is not the release's
  local_vendor()
  root <- local_host(files = list("R/vanilla-a.R" = "mine <- 1\n"))
  refused(root, "has the name of an engine file but is not in the lock")
  # Inside inlavanilla itself
  root <- local_host(desc = upstream_desc("26.11.24"))
  refused(root, "is inlavanilla itself")
})

test_that("a release with roxygen tags or without a locked file is refused", {
  releases <- fake_releases()
  releases[["v26.11.25"]] <- fake_release(
    "26.11.25",
    extra = list("R/vanilla-a.R" = paste0("#' @export\n", engine_a))
  )
  releases[["v26.11.26"]] <- fake_release("26.11.26", c = FALSE)
  local_vendor(fake_github(releases))
  root <- local_host()
  quiet(use_vanilla(release = "v26.11.24", path = root))
  before <- folder_state(root)
  expect_error(
    quiet(use_vanilla(release = "v26.11.25", path = root)),
    "has roxygen tags",
    class = "inlavanilla_vendor_error"
  )
  expect_error(
    quiet(use_vanilla(release = "v26.11.26", path = root)),
    "v26.11.26 no longer has 'R/vanilla-c.R'",
    class = "inlavanilla_vendor_error"
  )
  expect_identical(folder_state(root), before)
  # Once the copy is deleted, it is dropped from the lock
  unlink(file.path(root, "R", "vanilla-c.R"))
  expect_snapshot(out <- use_vanilla(release = "v26.11.26", path = root))
  expect_identical(out$file, c("R/vanilla-a.R", "R/vanilla-b.R"))
  lock <- vendor_read_lock(file.path(root, "inst", "inlavanilla.lock"))
  expect_identical(lock$files$file, c("R/vanilla-a.R", "R/vanilla-b.R"))
})

test_that("an unlisted copy identical to the release's file is accepted", {
  local_vendor()
  root <- local_host(
    files = list("R/vanilla-a.R" = gsub("\n", "\r\n", engine_a))
  )
  out <- quiet(use_vanilla(path = root))
  expect_identical(out$state, rep("new", 3))
  root <- local_host(files = list("R/vanilla-a.R" = engine_a))
  expect_error(
    quiet(use_vanilla(files = "vanilla-b.R", path = root)),
    "has the name of an engine file",
    class = "inlavanilla_vendor_error"
  )
})

test_that("a failed safety check leaves DESCRIPTION and copies the engine", {
  local_vendor()
  local_mocked_bindings(vendor_desc_safe = function(...) "forced")
  root <- local_host()
  expect_snapshot(out <- use_vanilla(path = root))
  expect_identical(read_text(file.path(root, "DESCRIPTION")), host_desc_minimal)
  expect_length(attr(out, "open"), 7L)
  expect_true(file.exists(file.path(root, "inst", "inlavanilla.lock")))
})

test_that("an Authors@R that fails the check leaves the other edits", {
  local_vendor()
  real_authors <- vendor_plan_authors
  local_mocked_bindings(
    vendor_plan_authors = function(d, holders) {
      out <- real_authors(d, holders)
      out$d$lines <- sub("\"Ann\"", "\"Anne\"", out$d$lines, fixed = TRUE)
      out
    }
  )
  root <- local_host()
  text <- messages_of(out <- use_vanilla(path = root))
  expect_match(
    text,
    "Authors@R could not be edited safely, so it was left as it was",
    fixed = TRUE
  )
  desc <- read_text(file.path(root, "DESCRIPTION"))
  expect_match(desc, "Imports: cli (>= 3.6.0), stats, utils\n", fixed = TRUE)
  expect_match(
    desc,
    "Authors@R: person(\"Ann\", \"Example\", email",
    fixed = TRUE
  )
  expect_identical(
    attr(out, "open"),
    "Add Haziq Jamil and Thomas Luu to Authors@R, with the role \"cph\"."
  )
})

test_that("clashes, a license and a Collate field are left to do", {
  local_vendor()
  desc <- sub(
    "GPL (>= 3)",
    "MIT + file LICENSE",
    host_desc_minimal,
    fixed = TRUE
  )
  desc <- paste0(desc, "Collate: 'hello.R' 'zzz.R'\n")
  root <- local_host(
    desc = desc,
    files = list("R/zzz.R" = "b_fun <- function() 0\n", "air.toml" = ""),
    namespace = "importFrom(rlang, abort)\nexport(hello)\n"
  )
  dir.create(file.path(root, ".git"))
  expect_snapshot(out <- use_vanilla(path = root))
  expect_length(attr(out, "open"), 3L)
})

test_that("the clone route copies the release from the clone", {
  skip_if_not(has_git(), "needs git")
  clone <- local_clone()
  root <- local_host()
  local_vendor(fake_github_for(clone))
  out <- quiet(use_vanilla(
    release = "v26.11.24",
    clone = clone$dir,
    path = root
  ))
  expect_identical(attr(out, "commit"), clone$commits[["v26.11.24"]])
  expect_identical(
    read_text(file.path(root, "R", "vanilla-b.R")),
    engine_b("26.11.24")
  )
  # Offline, with a line that GitHub was not checked
  local_vendor(fake_github(offline = TRUE), installed = "26.11.24-1")
  text <- messages_of(out <- use_vanilla(clone = clone$dir, path = root))
  expect_match(text, "Could not reach GitHub to check that v26.11.24-1 is at")
  expect_identical(attr(out, "release"), "v26.11.24-1")
  # A tag at another commit on GitHub is refused
  moved <- clone
  moved$commits[["v26.11.24-1"]] <- strrep("0", 40)
  local_vendor(fake_github_for(moved))
  expect_error(
    quiet(use_vanilla(release = "v26.11.24-1", clone = clone$dir, path = root)),
    "in the clone but at 0000000 on GitHub",
    class = "inlavanilla_vendor_error"
  )
  expect_error(
    quiet(use_vanilla(release = "v26.11.30", clone = clone$dir, path = root)),
    "git fetch --tags",
    class = "inlavanilla_vendor_error"
  )
})

test_that("use_vanilla_test() writes the test file once", {
  local_vendor()
  root <- local_host()
  expect_error(
    quiet(use_vanilla_test(path = root)),
    "has no testthat tests",
    class = "inlavanilla_vendor_error"
  )
  dir.create(file.path(root, "tests", "testthat"), recursive = TRUE)
  expect_error(
    quiet(use_vanilla_test(path = root)),
    "There is no lock",
    class = "inlavanilla_vendor_error"
  )
  quiet(use_vanilla(path = root))
  expect_snapshot({
    path <- use_vanilla_test(path = root)
    use_vanilla_test(path = root)
  })
  template <- system.file(
    "templates",
    "test-vanilla-engine.R",
    package = "inlavanilla"
  )
  expect_identical(read_text(path), read_text(template))
  write_text(path, "# mine\n")
  expect_snapshot(use_vanilla_test(path = root))
  expect_identical(read_text(path), "# mine\n")
  quiet(use_vanilla_test(path = root, overwrite = TRUE))
  expect_identical(read_text(path), read_text(template))
})

test_that("the written test passes on good copies and fails on bad ones", {
  local_vendor()
  root <- local_host()
  dir.create(file.path(root, "tests", "testthat"), recursive = TRUE)
  quiet(use_vanilla(path = root))
  path <- quiet(use_vanilla_test(path = root))
  run <- function(file = path) {
    res <- as.data.frame(testthat::test_file(file, reporter = "silent"))
    c(
      failed = sum(res$failed),
      skipped = sum(res$skipped),
      passed = sum(res$nb) - sum(res$failed)
    )
  }
  good <- run()
  expect_identical(good[["failed"]], 0L)
  expect_identical(good[["skipped"]], 0L)
  write_text(file.path(root, "R", "vanilla-a.R"), paste0(engine_a, "# mine\n"))
  expect_gt(run()[["failed"]], 0L)
  write_text(file.path(root, "R", "vanilla-a.R"), engine_a)
  write_text(file.path(root, "R", "vanilla-x.R"), "x <- 1\n")
  expect_gt(run()[["failed"]], 0L)
  unlink(file.path(root, "R", "vanilla-x.R"))
  unlink(file.path(root, "inst", "inlavanilla.lock"))
  expect_gt(run()[["failed"]], 0L)
  # As R CMD check runs it, with no source tree two levels up
  alone <- tempfile()
  dir.create(file.path(alone, "tests", "testthat"), recursive = TRUE)
  on.exit(unlink(alone, recursive = TRUE), add = TRUE)
  file.copy(path, file.path(alone, "tests", "testthat"))
  skipped <- run(file.path(alone, "tests", "testthat", basename(path)))
  expect_identical(skipped[["failed"]], 0L)
  expect_identical(skipped[["skipped"]], 1L)
})

test_that("use_vanilla_ci() writes the workflow and .Rbuildignore", {
  root <- local_host()
  expect_snapshot({
    path <- use_vanilla_ci(path = root)
    use_vanilla_ci(path = root)
  })
  template <- system.file(
    "templates",
    "vanilla-status.yaml",
    package = "inlavanilla"
  )
  expect_identical(read_text(path), read_text(template))
  expect_identical(read_text(file.path(root, ".Rbuildignore")), "^\\.github$\n")
  write_text(file.path(root, ".Rbuildignore"), "^README\\.Rmd$")
  write_text(path, "on: push\n")
  expect_snapshot(use_vanilla_ci(path = root))
  expect_identical(read_text(path), "on: push\n")
  expect_identical(
    read_text(file.path(root, ".Rbuildignore")),
    "^README\\.Rmd$\n^\\.github$\n"
  )
  quiet(use_vanilla_ci(path = root, overwrite = TRUE))
  expect_identical(read_text(path), read_text(template))
  expect_identical(
    read_text(file.path(root, ".Rbuildignore")),
    "^README\\.Rmd$\n^\\.github$\n"
  )
})
