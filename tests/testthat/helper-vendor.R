# Fixtures for the tests of use_vanilla() and vanilla_status(): a fake engine
# of a few small files, a fake GitHub that serves releases of it without the
# network, host packages in temporary folders, and clones under git.

## ----- Files -----------------------------------------------------------------

write_text <- function(path, text) {
  force(text)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeBin(charToRaw(enc2utf8(text)), con)
  invisible(path)
}

read_text <- function(path) {
  text <- rawToChar(readBin(path, "raw", n = file.size(path)))
  Encoding(text) <- "UTF-8"
  text
}

# Runs `expr` when the test, or the function, whose frame is `envir` ends
defer_in <- function(expr, envir) {
  thunk <- as.call(list(function() expr))
  do.call(base::on.exit, list(thunk, TRUE, FALSE), envir = envir)
}

# Every file of a folder with its MD5, to show that a refusal changed nothing
folder_state <- function(root) {
  files <- sort(list.files(root, recursive = TRUE, all.files = TRUE))
  stats::setNames(unname(tools::md5sum(file.path(root, files))), files)
}

## ----- The fake engine and its releases --------------------------------------

engine_a <- "a_fun <- function() b_fun()\nabort <- function(...) stop(...)\n"
engine_b <- function(v) paste0("b_fun <- function() \"", v, "\"\n")
engine_c <- "c_fun <- function() \"c\"\n"

# inlavanilla's own Authors@R: three persons, two of them rights holders
upstream_authors <- paste0(
  "Authors@R: c(\n",
  "    person(\"Haziq\", \"Jamil\", , \"haziq.jamil@gmail.com\", ",
  "role = c(\"aut\", \"cre\", \"cph\"),\n",
  "           comment = c(ORCID = \"0000-0003-3298-1010\")),\n",
  "    person(\"H",
  intToUtf8(229L),
  "vard\", \"Rue\", role = \"ctb\",\n",
  "           comment = c(ORCID = \"0000-0002-0222-1881\", ",
  "\"Statistical and computational methodology\")),\n",
  "    person(\"Thomas\", \"Luu\", role = \"cph\",\n",
  "           comment = \"Skew-normal quantile approximation ",
  "(Apache License 2.0)\")\n",
  "  )\n"
)

upstream_desc <- function(version) {
  paste0(
    "Package: inlavanilla\n",
    "Title: INLA-Style Skew-Normal Marginals for Any Log-Posterior\n",
    "Version: ",
    version,
    "\n",
    upstream_authors,
    "License: GPL (>= 3)\n",
    "Encoding: UTF-8\n"
  )
}

# A release of the fake engine: DESCRIPTION and the engine files
fake_release <- function(version, b = version, c = TRUE, extra = list()) {
  files <- list(
    "DESCRIPTION" = upstream_desc(version),
    "R/vanilla-a.R" = engine_a,
    "R/vanilla-b.R" = engine_b(b),
    "R/methods.R" = "print_fit <- function(x) x\n"
  )
  if (c) {
    files[["R/vanilla-c.R"]] <- engine_c
  }
  utils::modifyList(files, extra)
}

# The usual releases: v26.11.24 and v26.11.24-1, which changes vanilla-b.R
fake_releases <- function() {
  list(
    "v26.11.24" = fake_release("26.11.24"),
    "v26.11.24-1" = fake_release("26.11.24-1")
  )
}

# A SHA whose first seven characters differ between releases
fake_sha <- function(i, prefix = "") {
  head <- paste0(prefix, sprintf("%02d", i))
  paste0(head, strrep("0", 40L - nchar(head)))
}

# One pkt-line of git's smart-HTTP protocol: four hexadecimal digits that give
# the length, those four included, then the text. A NUL separates the
# capabilities from the first ref.
pkt_line <- function(text, capabilities = NULL) {
  bytes <- charToRaw(text)
  if (!is.null(capabilities)) {
    bytes <- c(bytes, as.raw(0L), charToRaw(capabilities))
  }
  bytes <- c(bytes, charToRaw("\n"))
  c(charToRaw(sprintf("%04x", length(bytes) + 4L)), bytes)
}

# A fake GitHub: `releases` is a named list, tag to a named list of files
# (paths to text). Each tag gets a commit and an archive in the layout of
# GitHub's archive of a commit, and the ref list in git's smart-HTTP format,
# with the first tag annotated. fetch() replaces vendor_fetch(): it serves
# the ref list and the archives, and fails for any other URL. With offline =
# TRUE it fails for every URL.
fake_github <- function(
  releases = fake_releases(),
  offline = FALSE,
  envir = parent.frame()
) {
  if (offline) {
    return(list(fetch = function(url, dest, timeout = 60) FALSE))
  }
  dir <- tempfile("fake-github-")
  dir.create(file.path(dir, "src"), recursive = TRUE)
  defer_in(unlink(dir, recursive = TRUE), envir)
  tags <- names(releases)
  commits <- stats::setNames(
    vapply(seq_along(tags), function(i) fake_sha(i, "c"), ""),
    tags
  )
  for (tag in tags) {
    top <- paste0("inlavanilla-", commits[[tag]])
    for (f in names(releases[[tag]])) {
      write_text(file.path(dir, "src", top, f), releases[[tag]][[f]])
    }
    old <- setwd(file.path(dir, "src"))
    utils::tar(
      file.path(dir, paste0(commits[[tag]], ".tar.gz")),
      top,
      compression = "gzip",
      tar = "internal"
    )
    setwd(old)
  }
  head <- if (length(tags) > 0L) commits[[length(tags)]] else fake_sha(0L)
  refs <- c(
    pkt_line("# service=git-upload-pack"),
    charToRaw("0000"),
    pkt_line(paste(head, "HEAD"), "multi_ack thin-pack side-band"),
    pkt_line(paste(head, "refs/heads/main")),
    pkt_line(paste(fake_sha(99L, "b"), "refs/heads/gh-pages"))
  )
  for (i in seq_along(tags)) {
    if (i == 1L) {
      refs <- c(
        refs,
        pkt_line(paste(fake_sha(i, "a"), paste0("refs/tags/", tags[i]))),
        pkt_line(paste(commits[[i]], paste0("refs/tags/", tags[i], "^{}")))
      )
    } else {
      refs <- c(
        refs,
        pkt_line(paste(commits[[i]], paste0("refs/tags/", tags[i])))
      )
    }
  }
  refs <- c(refs, charToRaw("0000"))
  urls <- character()
  fetch <- function(url, dest, timeout = 60) {
    urls <<- c(urls, url)
    if (identical(url, vendor_refs_url())) {
      writeBin(refs, dest)
      return(TRUE)
    }
    for (sha in commits) {
      if (identical(url, vendor_archive_url(sha))) {
        return(file.copy(
          file.path(dir, paste0(sha, ".tar.gz")),
          dest,
          overwrite = TRUE
        ))
      }
    }
    FALSE
  }
  list(
    fetch = fetch,
    commits = commits,
    refs = refs,
    dir = dir,
    urls = function() urls
  )
}

# Replaces the network, the date and the installed version in the test that
# calls it, and hides GitHub Actions from it (see local_github_actions()).
local_vendor <- function(
  gh = fake_github(envir = envir),
  installed = "26.11.24-1",
  today = "2026-11-25",
  envir = parent.frame()
) {
  testthat::local_mocked_bindings(
    vendor_fetch = gh$fetch,
    vendor_today = function() as.Date(today),
    vendor_installed_version = function() installed,
    .package = "inlavanilla",
    .env = envir
  )
  local_github_actions(NA, envir = envir)
  invisible(gh)
}

# Sets the environment variable GITHUB_ACTIONS to `value`, or unsets it when
# `value` is NA, until the test whose frame is `envir` ends. When it is
# "true", vanilla_status() also prints its failures and warnings as
# annotations, which would change the snapshots of its output on CI.
local_github_actions <- function(value = NA, envir = parent.frame()) {
  old <- Sys.getenv("GITHUB_ACTIONS", unset = NA)
  set <- function(v) {
    if (is.na(v)) {
      Sys.unsetenv("GITHUB_ACTIONS")
    } else {
      Sys.setenv(GITHUB_ACTIONS = v)
    }
  }
  set(value)
  defer_in(set(old), envir)
  invisible(old)
}

## ----- Host packages ---------------------------------------------------------

host_desc_minimal <- paste0(
  "Package: myhost\n",
  "Title: A Host of the Engine\n",
  "Version: 0.1.0\n",
  "Authors@R: person(\"Ann\", \"Example\", email = \"ann@example.org\", ",
  "role = c(\"aut\", \"cre\"))\n",
  "Description: Fits models with a copy of the engine of inlavanilla.\n",
  "License: GPL (>= 3)\n",
  "Encoding: UTF-8\n"
)

# The minimal DESCRIPTION after use_vanilla() edited it
minimal_expected <- paste0(
  "Package: myhost\n",
  "Title: A Host of the Engine\n",
  "Version: 0.1.0\n",
  "Authors@R: c(\n",
  "    person(\"Ann\", \"Example\", email = \"ann@example.org\", ",
  "role = c(\"aut\", \"cre\")),\n",
  "    person(\"Haziq\", \"Jamil\", role = \"cph\",\n",
  "           comment = c(ORCID = \"0000-0003-3298-1010\", ",
  "\"Code in R/vanilla-*.R, from inlavanilla\")),\n",
  "    person(\"Thomas\", \"Luu\", role = \"cph\",\n",
  "           comment = \"Skew-normal quantile approximation ",
  "(Apache License 2.0)\")\n",
  "  )\n",
  "Description: Fits models with a copy of the engine of inlavanilla.\n",
  "License: GPL (>= 3)\n",
  "Depends: R (>= 4.1.0)\n",
  "Imports: cli (>= 3.6.0), stats, utils\n",
  "Suggests: future, parallelly\n",
  "Encoding: UTF-8\n"
)

# A host package in a temporary folder whose path has a space, removed when
# the test ends.
local_host <- function(
  desc = host_desc_minimal,
  files = list(),
  namespace = "export(hello)\n",
  envir = parent.frame()
) {
  base <- tempfile("host-")
  root <- file.path(base, "my host")
  dir.create(file.path(root, "R"), recursive = TRUE)
  defer_in(unlink(base, recursive = TRUE), envir)
  write_text(file.path(root, "DESCRIPTION"), desc)
  if (!is.null(namespace)) {
    write_text(file.path(root, "NAMESPACE"), namespace)
  }
  write_text(file.path(root, "R", "hello.R"), "hello <- function() \"hi\"\n")
  for (f in names(files)) {
    write_text(file.path(root, f), files[[f]])
  }
  normalizePath(root, winslash = "/")
}

# The messages of `expr`, as one string with the lines joined by spaces, so
# that a pattern can match across the wrapping.
messages_of <- function(expr) {
  text <- character()
  withCallingHandlers(
    expr,
    message = function(m) {
      text <<- c(text, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  gsub("[[:space:]]+", " ", paste(text, collapse = " "))
}

## ----- Clones under git ------------------------------------------------------

has_git <- function() {
  nzchar(Sys.which("git"))
}

git_run <- function(dir, ...) {
  args <- c(
    "-C",
    dir,
    "-c",
    "user.name=Fixture",
    "-c",
    "user.email=fixture@example.org",
    "-c",
    "commit.gpgsign=false",
    "-c",
    "tag.gpgsign=false",
    "-c",
    "core.autocrlf=false",
    ...
  )
  out <- suppressWarnings(system2(
    "git",
    shQuote(args),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("git failed: ", paste(out, collapse = "\n"))
  }
  out
}

# A clone of inlavanilla with one commit per release, each tagged, the first
# with an annotated tag. Returns the folder and the commit of each tag.
local_clone <- function(releases = fake_releases(), envir = parent.frame()) {
  dir <- tempfile("clone ")
  dir.create(dir)
  defer_in(unlink(dir, recursive = TRUE), envir)
  git_run(dir, "init", "--quiet")
  commits <- character()
  for (i in seq_along(releases)) {
    tag <- names(releases)[i]
    unlink(file.path(dir, "R"), recursive = TRUE)
    for (f in names(releases[[i]])) {
      write_text(file.path(dir, f), releases[[i]][[f]])
    }
    git_run(dir, "add", "--all")
    git_run(dir, "commit", "--quiet", "-m", paste("Release", tag))
    if (i == 1L) {
      git_run(dir, "tag", "-a", "-m", tag, tag)
    } else {
      git_run(dir, "tag", tag)
    }
    commits[[tag]] <- git_run(dir, "rev-parse", "HEAD")
  }
  list(dir = normalizePath(dir, winslash = "/"), commits = commits)
}

# A fake GitHub whose releases are at the commits of a clone
fake_github_for <- function(
  clone,
  releases = fake_releases(),
  envir = parent.frame()
) {
  gh <- fake_github(releases, envir = envir)
  refs <- c(pkt_line("# service=git-upload-pack"), charToRaw("0000"))
  first <- TRUE
  for (tag in names(clone$commits)) {
    sha <- clone$commits[[tag]]
    if (first) {
      refs <- c(refs, pkt_line(paste(sha, "HEAD"), "multi_ack"))
      first <- FALSE
    }
    refs <- c(refs, pkt_line(paste(sha, paste0("refs/tags/", tag))))
  }
  refs <- c(refs, charToRaw("0000"))
  gh$fetch <- function(url, dest, timeout = 60) {
    if (identical(url, vendor_refs_url())) {
      writeBin(refs, dest)
      return(TRUE)
    }
    FALSE
  }
  gh
}
