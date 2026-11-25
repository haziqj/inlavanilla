# Releases on GitHub and in a clone: the ref list, the release to copy, and
# snapshots of a commit.

test_that("the ref list gives the tags and their commits", {
  gh <- fake_github()
  tags <- vendor_parse_refs(gh$refs)
  # The annotated first tag is peeled to its commit. Branches and the
  # capabilities after the NUL are not tags.
  expect_identical(tags, gh$commits)
  # The ref list as GitHub sends it, with its real line lengths
  live <- c(
    pkt_line("# service=git-upload-pack"),
    charToRaw("0000"),
    pkt_line(
      paste(strrep("0", 40), "HEAD"),
      "multi_ack symref=HEAD:refs/heads/main"
    ),
    pkt_line(paste(strrep("1", 40), "refs/pull/1/head")),
    pkt_line(paste(strrep("2", 40), "refs/tags/v26.11.24")),
    charToRaw("0000")
  )
  expect_identical(
    vendor_parse_refs(live),
    c("v26.11.24" = strrep("2", 40))
  )
  expect_null(vendor_parse_refs(charToRaw("<!DOCTYPE html><html></html>")))
  expect_null(vendor_parse_refs(raw()))
  expect_null(vendor_parse_refs(c(charToRaw("00ff"), charToRaw("short"))))
  expect_identical(
    vendor_parse_refs(charToRaw("0000")),
    stats::setNames(character(), character())
  )
})

test_that("the tags come from GitHub, or NULL offline", {
  gh <- local_vendor()
  expect_identical(vendor_github_tags(), gh$commits)
  local_vendor(fake_github(offline = TRUE))
  expect_null(vendor_github_tags())
  expect_null(vendor_tags(NULL))
})

test_that("the release to copy is resolved from the tags", {
  tags <- c(
    "v26.11.23" = strrep("1", 40),
    "v26.11.24" = strrep("2", 40),
    "v26.11.24-1" = strrep("3", 40),
    "v0.1.1" = strrep("4", 40)
  )
  res <- vendor_resolve_release("26.11.24", tags, "26.11.24-1", NULL)
  expect_identical(res$release, "v26.11.24")
  expect_identical(res$commit, strrep("2", 40))
  expect_true(res$explicit)
  expect_length(res$notes, 0L)
  expect_error(
    vendor_resolve_release("main", tags, "26.11.24-1", NULL),
    "is not a release of inlavanilla",
    class = "inlavanilla_vendor_error"
  )
  expect_error(
    vendor_resolve_release("v26.11.23", tags, "26.11.24-1", NULL),
    "older than v26.11.24, the first release",
    class = "inlavanilla_vendor_error"
  )
  err <- expect_error(
    vendor_resolve_release("v26.12.01", tags, "26.11.24-1", NULL),
    "no release",
    class = "inlavanilla_vendor_error"
  )
  expect_match(conditionMessage(err), "The newest release is \"v26.11.24-1\"")
  # By default the installed version, when it is a release
  res <- vendor_resolve_release(NULL, tags, "26.11.24", NULL)
  expect_identical(res$release, "v26.11.24")
  expect_false(res$explicit)
  # Else the newest release before it, with a note
  res <- vendor_resolve_release(NULL, tags, "26.11.25", NULL)
  expect_identical(res$release, "v26.11.24-1")
  expect_length(res$notes, 1L)
  expect_identical(
    res$notes[[1]]$plain,
    paste(
      "inlavanilla 26.11.25, installed here, is not a release on GitHub, so",
      "this copies v26.11.24-1, the newest release before it."
    )
  )
  expect_error(
    vendor_resolve_release(NULL, tags, "0.1.1.9003", NULL),
    "installed here, is older than v26.11.24",
    class = "inlavanilla_vendor_error"
  )
  # Going back: refused by default, allowed when asked
  lock <- list(header = c(Release = "v26.11.24-1"))
  expect_error(
    vendor_resolve_release(NULL, tags, "26.11.24", lock),
    "The copy is at v26.11.24-1, newer than v26.11.24",
    class = "inlavanilla_vendor_error"
  )
  res <- vendor_resolve_release("v26.11.24", tags, "26.11.24", lock)
  expect_identical(
    res$notes[[1]]$plain,
    "Going back from v26.11.24-1 to v26.11.24."
  )
  # In a clone
  expect_error(
    vendor_resolve_release("v26.12.01", tags, "26.11.24", NULL, clone = "x"),
    "git fetch --tags",
    class = "inlavanilla_vendor_error"
  )
})

test_that("a snapshot holds the DESCRIPTION and the engine files of a commit", {
  gh <- local_vendor()
  sha <- gh$commits[["v26.11.24"]]
  snap <- vendor_snapshot(sha, "v26.11.24")
  on.exit(unlink(snap$stage, recursive = TRUE), add = TRUE)
  expect_identical(
    snap$files,
    c("R/vanilla-a.R", "R/vanilla-b.R", "R/vanilla-c.R")
  )
  expect_identical(unname(snap$desc["Version"]), "26.11.24")
  expect_identical(read_text(file.path(snap$dir, "R", "vanilla-a.R")), engine_a)
  # Only the DESCRIPTION and R/vanilla-*.R are extracted
  expect_false(file.exists(file.path(snap$dir, "R", "methods.R")))
  expect_error(
    vendor_snapshot(fake_sha(7L), "v26.11.25"),
    "Could not download release v26.11.25",
    class = "inlavanilla_source_error"
  )
})

test_that("an archive of another commit, or with other members, is refused", {
  dir <- tempfile()
  dir.create(file.path(dir, "src"), recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  sha <- fake_sha(1L, "c")
  make <- function(top, files) {
    for (f in names(files)) {
      write_text(file.path(dir, "src", top, f), files[[f]])
    }
    archive <- tempfile(tmpdir = dir, fileext = ".tar.gz")
    old <- setwd(file.path(dir, "src"))
    on.exit(setwd(old))
    utils::tar(archive, top, compression = "gzip", tar = "internal")
    archive
  }
  other <- make(paste0("inlavanilla-", fake_sha(2L, "c")), fake_release("1"))
  stage <- tempfile(tmpdir = dir)
  dir.create(stage)
  expect_error(
    vendor_extract(other, sha, stage),
    "is not that commit",
    class = "inlavanilla_source_error"
  )
  nested <- make(
    paste0("inlavanilla-", sha),
    c(fake_release("1"), list("R/sub/vanilla-x.R" = "x <- 1\n"))
  )
  out <- vendor_extract(nested, sha, stage)
  expect_identical(
    sort(list.files(out, recursive = TRUE)),
    c("DESCRIPTION", "R/vanilla-a.R", "R/vanilla-b.R", "R/vanilla-c.R")
  )
  no_desc <- make(
    paste0("inlavanilla-", fake_sha(3L, "c")),
    list("R/vanilla-a.R" = engine_a)
  )
  expect_error(
    vendor_extract(no_desc, fake_sha(3L, "c"), stage),
    "has no 'DESCRIPTION'",
    class = "inlavanilla_source_error"
  )
})

test_that("a clone gives the tags and the same snapshot as GitHub", {
  skip_if_not(has_git(), "needs git")
  clone <- local_clone()
  tags <- vendor_clone_tags(clone$dir)
  # The annotated tag is peeled
  expect_identical(tags[names(clone$commits)], clone$commits)
  snap <- vendor_snapshot(clone$commits[["v26.11.24"]], "v26.11.24", clone$dir)
  on.exit(unlink(snap$stage, recursive = TRUE), add = TRUE)
  expect_identical(
    snap$files,
    c("R/vanilla-a.R", "R/vanilla-b.R", "R/vanilla-c.R")
  )
  expect_identical(
    read_text(file.path(snap$dir, "R", "vanilla-b.R")),
    engine_b("26.11.24")
  )
  not_git <- tempfile()
  dir.create(not_git)
  on.exit(unlink(not_git, recursive = TRUE), add = TRUE)
  expect_error(
    vendor_clone_tags(not_git),
    "is not a git repository",
    class = "inlavanilla_source_error"
  )
  # GitHub's tags beside a clone's
  local_vendor(fake_github_for(clone))
  src <- vendor_tags(clone$dir)
  expect_identical(src$github, clone$commits)
  expect_identical(
    vendor_github_state("v26.11.24", clone$commits[["v26.11.24"]], src$github),
    "same"
  )
  expect_identical(
    vendor_github_state("v26.11.24", strrep("0", 40), src$github),
    "moved"
  )
  expect_identical(vendor_github_state("v9.9.9", "x", src$github), "missing")
  expect_identical(vendor_github_state("v26.11.24", "x", NULL), "offline")
})

test_that("GitHub serves the releases of inlavanilla", {
  skip_on_cran()
  skip_if_offline("github.com")
  # skip_if_offline() checks only that the name resolves, so a proxy or a
  # firewall that blocks GitHub skips here.
  tags <- vendor_github_tags()
  skip_if(is.null(tags), "GitHub cannot be reached")
  expect_identical(
    tags[["v26.11.24"]],
    "88056cfa545cee5a601af4485bbc2a3a5443e0fc"
  )
  snap <- tryCatch(
    vendor_snapshot(tags[["v26.11.24"]], "v26.11.24"),
    inlavanilla_source_error = function(e) NULL
  )
  skip_if(is.null(snap), "GitHub's archive of v26.11.24 cannot be reached")
  on.exit(unlink(snap$stage, recursive = TRUE), add = TRUE)
  expect_length(snap$files, 13L)
  expect_identical(
    vendor_md5(file.path(
      snap$dir,
      c("R/vanilla-skew-normal.R", "R/vanilla-qsnorm-fast.R")
    )),
    c("faa33742c788625510ca1b332aae47d5", "000ab78ae2b4b2c2a8d7222704be8617")
  )
})
