# The names check of an engine copy

engine_host <- function(
  files = list(),
  namespace = "export(hello)\n",
  envir = parent.frame()
) {
  local_host(
    files = utils::modifyList(
      list("R/vanilla-a.R" = engine_a, "R/vanilla-b.R" = engine_b("1")),
      files
    ),
    namespace = namespace,
    envir = envir
  )
}

test_that("an engine name defined in another file of R/ is a failure", {
  root <- engine_host(list(
    "R/zzz.R" = "b_fun <- function() 2\n",
    "R/unix/x.R" = "abort <- function(...) NULL\n"
  ))
  items <- vendor_check_names(root)
  expect_identical(items$kind, c("duplicate", "duplicate"))
  expect_identical(items$level, c("fail", "fail"))
  expect_identical(items$name, c("abort", "b_fun"))
  expect_identical(items$files[[1]], c("R/vanilla-a.R", "R/unix/x.R"))
  expect_identical(
    items$text[[2]],
    "`b_fun`: 'R/vanilla-b.R' and 'R/zzz.R'"
  )
  expect_identical(nrow(vendor_check_names(engine_host())), 0L)
})

test_that("bindings in blocks, if bodies, assign() and setGeneric() count", {
  code <- paste(
    "{ x1 <- 1; x2 = 2 }",
    "if (TRUE) x3 <- 3 else x4 <- 4",
    "assign(\"x5\", 5)",
    "delayedAssign(x = \"x6\", 6)",
    "methods::setGeneric(\"x7\", function(x) standardGeneric(\"x7\"))",
    "x8 <- x9 <- function() NULL",
    "f(x10 <- 1)",
    sep = "\n"
  )
  path <- tempfile(fileext = ".R")
  on.exit(unlink(path), add = TRUE)
  write_text(path, code)
  expect_setequal(vendor_top_names(path)$names, paste0("x", 1:9))
})

test_that("imports that mask engine names are found", {
  median_engine <- "median <- function(x) 0\n"
  root <- engine_host(
    list("R/vanilla-c.R" = median_engine),
    namespace = "importFrom(rlang, abort)\nimport(stats)\n"
  )
  items <- vendor_check_names(root)
  expect_identical(items$kind, c("import_from", "import_all"))
  expect_identical(items$level, c("fail", "warn"))
  expect_identical(items$pkg, c("rlang", "stats"))
  expect_identical(items$names[[2]], "median")
  expect_match(
    items$text[[1]],
    "`abort`: 'R/vanilla-a.R', and imported from rlang in 'NAMESPACE'",
    fixed = TRUE
  )
  root <- engine_host(
    list("R/vanilla-c.R" = median_engine),
    namespace = paste0(
      "import(stats, except = c(\"median\"))\n",
      "import(notinstalledpackage.xyz)\n"
    )
  )
  items <- vendor_check_names(root)
  expect_identical(items$kind, "import_unchecked")
  expect_identical(items$level, "note")
  expect_match(items$text, "notinstalledpackage.xyz is not installed")
})

test_that("an export pattern of engine names warns, an export() is a note", {
  root <- engine_host(namespace = "exportPattern(\"^[[:alpha:]]+\")\n")
  items <- vendor_check_names(root)
  expect_identical(items$kind, "export_pattern")
  expect_identical(items$level, "warn")
  expect_identical(items$count, 3L)
  root <- engine_host(namespace = "export(a_fun)\nexport(hello)\n")
  items <- vendor_check_names(root)
  expect_identical(items$kind, "export")
  expect_identical(items$level, "note")
  expect_identical(items$names[[1]], "a_fun")
})

test_that("a file that does not parse is a warning, not an error", {
  root <- engine_host(list("R/broken.R" = "f <- function( {\n"))
  items <- vendor_check_names(root)
  expect_identical(items$kind, "parse")
  expect_identical(items$level, "warn")
  expect_identical(items$files[[1]], "R/broken.R")
})

test_that("with a snapshot, a copy that needs an uncopied engine file warns", {
  root <- local_host(files = list("R/vanilla-a.R" = engine_a))
  snap_dir <- tempfile()
  on.exit(unlink(snap_dir, recursive = TRUE), add = TRUE)
  write_text(file.path(snap_dir, "R", "vanilla-a.R"), engine_a)
  write_text(file.path(snap_dir, "R", "vanilla-b.R"), engine_b("1"))
  write_text(file.path(snap_dir, "R", "vanilla-c.R"), engine_c)
  snap <- list(
    dir = snap_dir,
    files = c("R/vanilla-a.R", "R/vanilla-b.R", "R/vanilla-c.R")
  )
  items <- vendor_check_names(root, snap)
  expect_identical(items$kind, "missing")
  expect_identical(items$level, "warn")
  expect_identical(items$name, "R/vanilla-b.R")
  expect_identical(items$names[[1]], "b_fun")
  expect_identical(
    items$text,
    "'R/vanilla-b.R' binds `b_fun`, which the copies use."
  )
  # Not when the package defines the name itself
  write_text(file.path(root, "R", "mine.R"), "b_fun <- function() 1\n")
  expect_identical(nrow(vendor_check_names(root, snap)), 0L)
})
