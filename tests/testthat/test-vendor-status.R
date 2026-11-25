# vanilla_status() on host packages, with a fake GitHub

quiet <- function(expr) suppressMessages(expr)

# A host with a copy of the fake engine at `release`
copied_host <- function(
  release = "v26.11.24",
  gh = NULL,
  files = NULL,
  envir = parent.frame()
) {
  if (is.null(gh)) {
    gh <- fake_github(envir = envir)
  }
  local_vendor(gh, installed = sub("^v", "", release), envir = envir)
  root <- local_host(envir = envir)
  quiet(use_vanilla(release = release, files = files, path = root))
  root
}

test_that("a healthy copy passes, with a newer release reported", {
  gh <- local_vendor()
  root <- copied_host("v26.11.24", gh)
  expect_snapshot(out <- vanilla_status(path = root))
  expect_identical(
    out$check,
    c(
      "lock",
      "copies",
      "names",
      "dependencies",
      "authors",
      "license",
      "release",
      "update"
    )
  )
  expect_identical(out$status, c(rep("pass", 7), "warn"))
  expect_identical(
    out$message[[8]],
    "A newer release changes copied files: v26.11.24-1 changes 'R/vanilla-b.R'."
  )
  expect_match(
    out$details[[8]],
    "compare/v26.11.24...v26.11.24-1",
    fixed = TRUE
  )
  expect_match(
    out$details[[8]],
    "inlavanilla::use_vanilla(release = \"v26.11.24-1\")",
    fixed = TRUE
  )
  # At the newest release
  root <- copied_host("v26.11.24-1", gh)
  out <- quiet(vanilla_status(path = root))
  expect_identical(out$status, rep("pass", 8))
  expect_identical(
    out$message[[8]],
    "No release of inlavanilla is newer than v26.11.24-1."
  )
})

test_that("a newer release that changes only an uncopied file passes", {
  releases <- list(
    "v26.11.24" = fake_release("26.11.24", b = "1"),
    "v26.11.24-1" = fake_release(
      "26.11.24-1",
      b = "1",
      extra = list("R/vanilla-c.R" = "c_fun <- function() \"c2\"\n")
    )
  )
  gh <- fake_github(releases)
  root <- copied_host("v26.11.24", gh, files = c("vanilla-a.R", "vanilla-b.R"))
  out <- quiet(vanilla_status(path = root))
  expect_identical(out$status[[8]], "pass")
  expect_identical(out$message[[8]], "v26.11.24-1 changes no copied file.")
  expect_identical(
    out$details[[8]],
    "v26.11.24-1 also has 'R/vanilla-c.R', which is not copied."
  )
})

test_that("offline, the release and update checks are skipped", {
  root <- copied_host("v26.11.24")
  local_vendor(fake_github(offline = TRUE), installed = "26.11.24-1")
  expect_snapshot(out <- vanilla_status(path = root, error = TRUE))
  expect_identical(out$status[7:8], c("skip", "skip"))
  expect_identical(out$status[1:6], rep("pass", 6))
  expect_match(out$details[[8]], "inlavanilla 26.11.24-1 is installed here")
})

test_that("error = TRUE stops after printing when a check fails", {
  root <- copied_host("v26.11.24-1")
  write_text(file.path(root, "R", "vanilla-a.R"), paste0(engine_a, "# mine\n"))
  text <- NULL
  err <- expect_error(
    text <- messages_of(vanilla_status(path = root, error = TRUE)),
    class = "inlavanilla_status_error"
  )
  expect_match(conditionMessage(err), "1 check of the engine copy failed")
  out <- quiet(vanilla_status(path = root))
  expect_identical(out$status[[2]], "fail")
  expect_identical(out$message[[2]], "1 engine file does not match the lock:")
  expect_match(
    out$details[[2]],
    "'R/vanilla-a.R' differs from the lock.",
    fixed = TRUE
  )
  expect_snapshot(vanilla_status(path = root))
})

test_that("without a lock, or with a malformed one, only the lock is checked", {
  local_vendor()
  root <- local_host()
  expect_snapshot(out <- vanilla_status(path = root))
  expect_identical(out$status, c("fail", rep("skip", 7)))
  write_text(file.path(root, "inst", "inlavanilla.lock"), "Upstream: x\n")
  out <- quiet(vanilla_status(path = root))
  expect_identical(out$status, c("fail", rep("skip", 7)))
  expect_match(out$message[[1]], "is malformed: its header lacks the fields")
})

test_that("the release check finds a lock that its release does not match", {
  gh <- local_vendor()
  root <- copied_host("v26.11.24", gh)
  lock_path <- file.path(root, "inst", "inlavanilla.lock")
  lock_text <- read_text(lock_path)
  release_row <- function() {
    out <- quiet(vanilla_status(path = root))
    out[out$check == "release", ]
  }
  # The tag moved
  write_text(
    lock_path,
    sub(gh$commits[["v26.11.24"]], gh$commits[["v26.11.24-1"]], lock_text)
  )
  row <- release_row()
  expect_identical(row$status, "fail")
  expect_match(
    row$details,
    "The tag v26.11.24 is at commit c010000, but the lock says c020000.",
    fixed = TRUE
  )
  # No such tag
  write_text(
    lock_path,
    sub("Release: v26.11.24", "Release: v26.11.30", lock_text)
  )
  row <- release_row()
  expect_match(
    row$details,
    "There is no tag v26.11.30 on GitHub.",
    fixed = TRUE
  )
  # Another version
  write_text(lock_path, sub("Version: 26.11.24", "Version: 9.9.9", lock_text))
  row <- release_row()
  expect_match(
    row$details,
    "Version is 9.9.9, but the DESCRIPTION of inlavanilla",
    fixed = TRUE
  )
  # A file not in the release, and an MD5 that is not the release's
  edited <- paste0(engine_a, "# mine\n")
  write_text(file.path(root, "R", "vanilla-a.R"), edited)
  write_text(file.path(root, "R", "vanilla-z.R"), "z <- 1\n")
  md5 <- vendor_md5(file.path(root, "R", c("vanilla-a.R", "vanilla-z.R")))
  lock <- vendor_read_lock(lock_path)
  files <- lock$files
  files$md5[files$file == "R/vanilla-a.R"] <- md5[[1]]
  files <- rbind(files, data.frame(file = "R/vanilla-z.R", md5 = md5[[2]]))
  vendor_write_lock(lock_path, lock$header, files)
  out <- quiet(vanilla_status(path = root))
  expect_identical(out$status[[2]], "pass")
  row <- out[out$check == "release", ]
  expect_identical(row$status, "fail")
  expect_match(
    row$details,
    "'R/vanilla-a.R' has another MD5 at v26.11.24 than in the lock.",
    fixed = TRUE
  )
  expect_match(
    row$details,
    "'R/vanilla-z.R' is not in inlavanilla at c010000.",
    fixed = TRUE
  )
})

test_that("a commit that cannot be downloaded skips the release check", {
  gh <- local_vendor()
  root <- copied_host("v26.11.24", gh)
  refs_only <- gh
  refs_only$fetch <- function(url, dest, timeout = 60) {
    if (identical(url, vendor_refs_url())) {
      gh$fetch(url, dest, timeout)
    } else {
      FALSE
    }
  }
  local_vendor(refs_only, installed = "26.11.24")
  text <- messages_of(out <- vanilla_status(path = root))
  expect_identical(out$status[7:8], c("skip", "skip"))
  expect_match(
    text,
    "Skipped: could not download commit c010000.",
    fixed = TRUE
  )
  expect_match(text, "Skipped: could not get v26.11.24-1", fixed = TRUE)
})

test_that("inside inlavanilla, nothing is checked", {
  root <- local_host(desc = upstream_desc("26.11.24"))
  expect_snapshot(out <- vanilla_status(path = root))
  expect_identical(out$status, rep("skip", 8))
})

test_that("on GitHub Actions, failures and warnings are annotations", {
  gh <- local_vendor()
  root <- copied_host("v26.11.24", gh)
  write_text(file.path(root, "R", "vanilla-a.R"), paste0(engine_a, "# mine\n"))
  old <- Sys.getenv("GITHUB_ACTIONS", unset = NA)
  Sys.setenv(GITHUB_ACTIONS = "true")
  on.exit(
    if (is.na(old)) {
      Sys.unsetenv("GITHUB_ACTIONS")
    } else {
      Sys.setenv(GITHUB_ACTIONS = old)
    },
    add = TRUE
  )
  out <- utils::capture.output(quiet(vanilla_status(path = root)))
  expect_match(
    out,
    "^::error::1 engine file does not match the lock:%0A",
    all = FALSE
  )
  expect_match(
    out,
    "^::warning::A newer release changes copied files",
    all = FALSE
  )
  expect_length(out, 2L)
})

test_that("an installed inlavanilla older than the copy is pointed out", {
  root <- copied_host("v26.11.24-1")
  local_vendor(installed = "26.11.24")
  text <- messages_of(vanilla_status(path = root))
  expect_match(
    text,
    "inlavanilla 26.11.24, installed here, is older than the copy, v26.11.24-1",
    fixed = TRUE
  )
})

test_that("the clone route checks the release in the clone", {
  skip_if_not(has_git(), "needs git")
  clone <- local_clone()
  local_vendor(fake_github_for(clone), installed = "26.11.24")
  root <- local_host()
  quiet(use_vanilla(clone = clone$dir, path = root))
  out <- quiet(vanilla_status(path = root, clone = clone$dir))
  expect_identical(out$status, c(rep("pass", 7), "warn"))
  local_vendor(fake_github(offline = TRUE), installed = "26.11.24")
  out <- quiet(vanilla_status(path = root, clone = clone$dir))
  expect_identical(out$status[[7]], "pass")
  expect_match(
    out$details[[7]],
    "Could not reach GitHub to check that v26.11.24"
  )
})

test_that("the DESCRIPTION rows name use_vanilla() only for what it fixes", {
  gh <- local_vendor()
  root <- copied_host("v26.11.24", gh)
  desc_path <- file.path(root, "DESCRIPTION")
  desc <- read_text(desc_path)
  desc <- sub("Depends: R (>= 4.1.0)", "Depends: R (> 4.0)", desc, fixed = TRUE)
  desc <- sub("Suggests: future, ", "Suggests: ", desc, fixed = TRUE)
  write_text(desc_path, paste0(desc, "Collate: 'hello.R'\n"))
  row <- function() {
    out <- quiet(vanilla_status(path = root))
    out[out$check == "dependencies", ]
  }
  dep <- row()
  expect_identical(dep$status, "fail")
  expect_identical(
    strsplit(dep$details, "\n")[[1]],
    c(
      "`R (>= 4.1.0)` in Depends: it has `R (> 4.0)`",
      "future in Suggests",
      paste(
        "the engine files in Collate: it lacks 'R/vanilla-a.R',",
        "'R/vanilla-b.R' and 'R/vanilla-c.R'"
      ),
      paste(
        "`inlavanilla::use_vanilla(release = \"v26.11.24\")` adds future in",
        "Suggests, and copies nothing new."
      ),
      paste(
        "Change `R (> 4.0)` in Depends to `R (>= 4.1.0)` by hand:",
        "`use_vanilla()` does not change it."
      ),
      paste(
        "Add the engine files to Collate by hand: `use_vanilla()` does not",
        "change it."
      )
    )
  )
  # After that fix, only the fixes by hand are left
  quiet(use_vanilla(release = "v26.11.24", path = root))
  dep <- row()
  expect_identical(dep$status, "fail")
  details <- strsplit(dep$details, "\n")[[1]]
  expect_length(details, 4L)
  expect_false(any(grepl("use_vanilla(release", details, fixed = TRUE)))
  # When use_vanilla() fixes every item, it says so in one line
  write_text(desc_path, host_desc_minimal)
  dep <- row()
  expect_identical(
    strsplit(dep$details, "\n")[[1]][7],
    paste(
      "`inlavanilla::use_vanilla(release = \"v26.11.24\")` adds them, and",
      "copies nothing new."
    )
  )
})

test_that("the authors row names use_vanilla() only when it can add them", {
  gh <- local_vendor()
  root <- copied_host("v26.11.24", gh)
  desc <- sub(
    "Encoding: UTF-8\n",
    paste0(
      "Depends: R (>= 4.1.0)\n",
      "Imports: cli (>= 3.6.0), stats, utils\n",
      "Suggests: future, parallelly\n",
      "Encoding: UTF-8\n"
    ),
    host_desc_minimal,
    fixed = TRUE
  )
  write_text(file.path(root, "DESCRIPTION"), desc)
  authors_row <- function() {
    out <- quiet(vanilla_status(path = root))
    out[out$check == "authors", ]
  }
  row <- authors_row()
  expect_identical(row$status, "fail")
  expect_identical(
    row$details,
    paste(
      "`inlavanilla::use_vanilla(release = \"v26.11.24\")` adds the copyright",
      "holders, and copies nothing new."
    )
  )
  # An edit of Authors@R that would fail the safety check
  real_authors <- vendor_plan_authors
  local_mocked_bindings(
    vendor_plan_authors = function(d, holders) {
      out <- real_authors(d, holders)
      out$d$lines <- sub("\"Ann\"", "\"Anne\"", out$d$lines, fixed = TRUE)
      out
    }
  )
  row <- authors_row()
  expect_identical(row$status, "fail")
  expect_identical(
    row$details,
    paste(
      "Add Haziq Jamil and Thomas Luu to Authors@R by hand, with the role",
      "\"cph\": `use_vanilla()` cannot edit this Authors@R safely."
    )
  )
})

test_that("a lock of another upstream skips the network checks in one line", {
  gh <- local_vendor()
  root <- copied_host("v26.11.24", gh)
  lock_path <- file.path(root, "inst", "inlavanilla.lock")
  write_text(
    lock_path,
    sub(
      "https://github.com/haziqj/inlavanilla",
      "https://github.com/someone/fork",
      read_text(lock_path),
      fixed = TRUE
    )
  )
  text <- messages_of(out <- vanilla_status(path = root))
  expect_identical(out$status[[1]], "warn")
  expect_identical(out$status[7:8], c("skip", "skip"))
  expect_identical(
    lengths(regmatches(text, gregexpr("Skipped", text, fixed = TRUE))),
    1L
  )
  expect_match(
    text,
    paste(
      "Skipped 2 checks: the lock does not copy from",
      "<https://github.com/haziqj/inlavanilla>."
    ),
    fixed = TRUE
  )
})
