# The engine of inlavanilla: the files R/vanilla-*.R are byte copies of the
# engine files of an inlavanilla release, listed with their MD5 hashes in
# inst/inlavanilla.lock. This test fails when a copy differs from the lock,
# when a listed file is missing, or when a file R/vanilla-*.R is not in the
# lock. Never edit a copy: update it with inlavanilla::use_vanilla(). The test
# reads the sources in R/, so it skips under R CMD check, which installs none.
# Written by inlavanilla::use_vanilla_test().

test_that("the engine files match inst/inlavanilla.lock", {
  root <- test_path("..", "..")
  skip_if_not(
    file.exists(file.path(root, "DESCRIPTION")) &&
      dir.exists(file.path(root, "R")),
    "needs the source tree of the package"
  )
  lock_path <- file.path(root, "inst", "inlavanilla.lock")
  if (!file.exists(lock_path)) {
    fail("inst/inlavanilla.lock is missing")
  } else {
    records <- read.dcf(lock_path)[-1, , drop = FALSE]
    has <- function(field) field %in% colnames(records)
    listed <- if (has("File")) unname(records[, "File"]) else character()
    md5 <- if (has("MD5")) unname(records[, "MD5"]) else character()
    # The MD5 of the text with carriage returns removed, so that a checkout
    # with Windows line endings passes.
    md5_lf <- function(path) {
      bytes <- readBin(path, "raw", n = file.size(path))
      tmp <- tempfile()
      on.exit(unlink(tmp))
      writeBin(bytes[bytes != as.raw(13L)], tmp)
      unname(tools::md5sum(tmp))
    }
    for (i in seq_along(listed)) {
      path <- file.path(root, listed[[i]])
      expect_true(file.exists(path), label = paste(listed[[i]], "exists"))
      if (file.exists(path)) {
        expect_identical(
          md5_lf(path),
          md5[[i]],
          label = paste("the MD5 of", listed[[i]])
        )
      }
    }
    copies <- list.files(file.path(root, "R"), "^vanilla-.+[.]R$")
    expect_identical(
      setdiff(file.path("R", copies), listed),
      character(),
      label = "the files R/vanilla-*.R that are not in the lock"
    )
  }
})
