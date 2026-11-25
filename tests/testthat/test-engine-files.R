# The engine files R/vanilla-*.R are copied byte for byte into INLAvaan, so
# they hold only internal code and call nothing from the other files of R/.
# The CI workflow source-checks.yaml runs the first checks on every push; these
# tests run them from the source tree, and skip without it, as in R CMD check.

engine_dir <- function() {
  dir <- testthat::test_path("..", "..", "R")
  skip_if_not(
    length(list.files(dir, "^vanilla-.*[.]R$")) > 0L,
    "needs the source tree"
  )
  dir
}

test_that("the engine files hold no roxygen, generic or package name", {
  dir <- engine_dir()
  for (f in list.files(dir, "^vanilla-.*[.]R$", full.names = TRUE)) {
    lines <- readLines(f)
    expect_false(any(grepl("^\\s*#'", lines)), label = basename(f))
    expect_false(
      any(grepl(
        "(UseMethod|NextMethod|registerS3method|setGeneric|setMethod)\\s*\\(|S3method",
        lines
      )),
      label = basename(f)
    )
    # The name in code, quoted or not, but not in comments or as the prefix
    # of a condition class (inlavanilla_error)
    expect_false(
      any(grepl("^[^#]*inlavanilla([^_]|$)", lines)),
      label = basename(f)
    )
  }
})

test_that("the engine files use nothing defined in the other files", {
  dir <- engine_dir()
  files <- list.files(dir, "[.]R$", full.names = TRUE)
  engine <- grepl("^vanilla-", basename(files))
  top_names <- function(f) {
    ex <- parse(f, keep.source = FALSE)
    unlist(lapply(ex, function(e) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) {
        as.character(e[[2]])
      }
    }))
  }
  furniture <- setdiff(
    unlist(lapply(files[!engine], top_names)),
    unlist(lapply(files[engine], top_names))
  )
  # inla_vanilla() names the call of the engine's errors, quoted, and
  # `diagnostics` is a list in vanilla_core()
  allowed <- c("inla_vanilla", "diagnostics")
  for (f in files[engine]) {
    used <- unique(all.names(parse(f, keep.source = FALSE)))
    expect_identical(
      setdiff(intersect(used, furniture), allowed),
      character(),
      label = basename(f)
    )
  }
})

test_that("the vendoring code reuses no engine name", {
  dir <- engine_dir()
  items <- vendor_check_names(dirname(dir))
  expect_identical(items$name[items$kind == "duplicate"], character())
})

test_that("vendor_needs covers the engine's packages", {
  dir <- engine_dir()
  used <- unlist(lapply(
    list.files(dir, "^vanilla-.*[.]R$", full.names = TRUE),
    function(f) {
      pd <- utils::getParseData(parse(f, keep.source = TRUE))
      pd$text[pd$token == "SYMBOL_PACKAGE"]
    }
  ))
  allowed <- c("base", names(vendor_needs$imports), vendor_needs$suggests)
  expect_identical(setdiff(unique(used), allowed), character())
})
