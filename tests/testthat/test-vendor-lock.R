# The lock of an engine copy, the MD5 of LF text, and the copies check

lock_header <- function(release = "v26.11.24", commit = NULL, version = NULL) {
  c(
    Upstream = "https://github.com/haziqj/inlavanilla",
    Release = release,
    Commit = if (is.null(commit)) {
      "88056cfa545cee5a601af4485bbc2a3a5443e0fc"
    } else {
      commit
    },
    Version = if (is.null(version)) sub("^v", "", release) else version,
    Synced = "2026-11-24",
    Hash = "md5 of the text with LF line endings"
  )
}

test_that("a lock has the format of INLAvaan's, byte for byte", {
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  files <- data.frame(
    file = c("R/vanilla-skew-normal.R", "R/vanilla-qsnorm-fast.R"),
    md5 = c(
      "faa33742c788625510ca1b332aae47d5",
      "000ab78ae2b4b2c2a8d7222704be8617"
    )
  )
  vendor_write_lock(path, lock_header(), files)
  expected <- paste0(
    "Upstream: https://github.com/haziqj/inlavanilla\n",
    "Release: v26.11.24\n",
    "Commit: 88056cfa545cee5a601af4485bbc2a3a5443e0fc\n",
    "Version: 26.11.24\n",
    "Synced: 2026-11-24\n",
    "Hash: md5 of the text with LF line endings\n",
    "\n",
    "File: R/vanilla-qsnorm-fast.R\n",
    "MD5: 000ab78ae2b4b2c2a8d7222704be8617\n",
    "\n",
    "File: R/vanilla-skew-normal.R\n",
    "MD5: faa33742c788625510ca1b332aae47d5\n"
  )
  expect_identical(read_text(path), expected)
})

test_that("a lock reads back as written, and with CRLF line endings", {
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  files <- data.frame(
    file = c("R/vanilla-a.R", "R/vanilla-b.R"),
    md5 = c(strrep("a", 32), strrep("b", 32))
  )
  vendor_write_lock(path, lock_header(), files)
  lock <- vendor_read_lock(path)
  expect_identical(lock$header, lock_header())
  expect_identical(lock$files, files)
  write_text(path, gsub("\n", "\r\n", read_text(path), fixed = TRUE))
  expect_identical(vendor_read_lock(path), lock)
  # A header field it does not know is dropped
  write_text(
    path,
    sub("Hash:", "Extra: x\nHash:", read_text(path), fixed = TRUE)
  )
  expect_identical(vendor_read_lock(path), lock)
  # A lock with no files
  vendor_write_lock(path, lock_header(), files[0, ])
  expect_identical(nrow(vendor_read_lock(path)$files), 0L)
})

test_that("a malformed lock is refused with its reason", {
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  good <- paste0(
    paste0(names(lock_header()), ": ", lock_header(), collapse = "\n"),
    "\n\nFile: R/vanilla-a.R\nMD5: ",
    strrep("a", 32),
    "\n"
  )
  refused <- function(text, reason) {
    write_text(path, text)
    err <- expect_error(
      vendor_read_lock(path),
      class = "inlavanilla_lock_error"
    )
    expect_match(err$reason, reason, fixed = TRUE)
  }
  refused(sub("Version: 26.11.24\n", "", good), "lacks the field Version")
  refused(sub("Hash:", "File: R/vanilla-x.R\nHash:", good), "File or MD5")
  refused(paste0(good, "Release: v1\n"), "a file record has a header field")
  refused(sub("88056cf", "88056CF", good), "Commit is not a SHA")
  refused(sub(strrep("a", 32), strrep("a", 31), good), "not 32 hexadecimal")
  refused(
    paste0(good, "\nFile: R/vanilla-a.R\nMD5: ", strrep("b", 32), "\n"),
    "lists a file twice"
  )
  refused(sub("R/vanilla-a.R", "R/a.R", good), "not R/vanilla-*.R")
  refused(sub("md5 of the text", "sha256 of the text", good), "its Hash is not")
  refused(sub("v26.11.24", "main", good), "Release is not a release tag")
  refused(sub("MD5: a+\n$", "", good), "lacks File or MD5")
  refused(
    "this is\nnot: DCF\n  at all\n\n  \nreally not\n",
    "not in DCF format"
  )
  refused("", "it is empty")
  err <- expect_error(
    vendor_read_lock(file.path(tempdir(), "no-such-lock")),
    class = "inlavanilla_lock_error"
  )
  expect_match(err$reason, "no such file")
})

test_that("the MD5 is that of the text with carriage returns removed", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- function(name, text) write_text(file.path(dir, name), text)
  md5 <- function(name) unname(tools::md5sum(file.path(dir, name)))
  lf <- path("lf", "a <- 1\nb <- 2\n")
  crlf <- path("crlf", "a <- 1\r\nb <- 2\r\n")
  no_final <- path("no-final", "a <- 1\nb <- 2")
  no_final_lf <- path("no-final-lf", "a <- 1\nb <- 2")
  lone_cr <- path("lone-cr", "a\rb")
  lone_cr_lf <- path("lone-cr-lf", "ab")
  empty <- path("empty", "")
  utf8 <- path("utf8", paste0("x <- \"", intToUtf8(c(229L, 963L)), "\"\r\n"))
  utf8_lf <- path(
    "utf8-lf",
    paste0("x <- \"", intToUtf8(c(229L, 963L)), "\"\n")
  )
  expect_identical(vendor_md5(lf), md5("lf"))
  expect_identical(vendor_md5(crlf), md5("lf"))
  expect_identical(vendor_md5(no_final), md5("no-final-lf"))
  expect_identical(vendor_md5(lone_cr), md5("lone-cr-lf"))
  expect_identical(vendor_md5(empty), "d41d8cd98f00b204e9800998ecf8427e")
  expect_identical(vendor_md5(utf8), md5("utf8-lf"))
  expect_identical(vendor_md5(c(lf, crlf)), rep(md5("lf"), 2))
  expect_identical(vendor_md5(character()), character())
})

test_that("the copies check finds missing, differing and unlisted files", {
  root <- local_host()
  write_text(file.path(root, "R", "vanilla-a.R"), engine_a)
  write_text(file.path(root, "R", "vanilla-b.R"), engine_b("1"))
  files <- data.frame(
    file = c("R/vanilla-a.R", "R/vanilla-b.R"),
    md5 = vendor_md5(file.path(root, c("R/vanilla-a.R", "R/vanilla-b.R")))
  )
  lock <- list(header = lock_header(), files = files)
  expect_identical(nrow(vendor_check_copies(root, lock)), 0L)
  # A checkout with CRLF line endings passes
  write_text(file.path(root, "R", "vanilla-a.R"), gsub("\n", "\r\n", engine_a))
  expect_identical(nrow(vendor_check_copies(root, lock)), 0L)
  write_text(file.path(root, "R", "vanilla-a.R"), paste0(engine_a, "# edit\n"))
  unlink(file.path(root, "R", "vanilla-b.R"))
  write_text(file.path(root, "R", "vanilla-mine.R"), "mine <- 1\n")
  problems <- vendor_check_copies(root, lock)
  expect_identical(
    problems$file,
    c("R/vanilla-a.R", "R/vanilla-b.R", "R/vanilla-mine.R")
  )
  expect_identical(problems$problem, c("differs", "missing", "unlisted"))
  expect_identical(problems$lock_md5[1:2], files$md5)
  expect_identical(is.na(problems$md5), c(FALSE, TRUE, FALSE))
})

test_that("newer tags are versions above the release, oldest first", {
  tags <- c("v26.12.1", "v26.11.24-1", "v26.11.24", "main", "v0.1.1", "x1")
  expect_identical(
    vendor_newer_tags(tags, "v26.11.24"),
    c("v26.11.24-1", "v26.12.1")
  )
  expect_identical(vendor_newer_tags(tags, "v26.12.1"), character())
  expect_true(
    vendor_tag_version("v26.11.24-1") > vendor_tag_version("v26.11.24")
  )
  commits <- c(
    "v0.1.1" = "a",
    "v26.11.24-1" = "c",
    "v26.11.24" = "b",
    main = "d"
  )
  expect_identical(
    names(vendor_release_tags(commits)),
    c("v26.11.24", "v26.11.24-1")
  )
})
