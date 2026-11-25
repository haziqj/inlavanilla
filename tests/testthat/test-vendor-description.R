# The edits of a host's DESCRIPTION, and the checks of vanilla_status() on it

fixture_holders <- function() {
  vendor_holders(trimws(sub("^Authors@R:", "", upstream_authors)))
}

plan_for <- function(
  desc,
  holders = fixture_holders(),
  envir = parent.frame()
) {
  root <- local_host(desc = desc, envir = envir)
  plan <- vendor_plan_description(root, holders)
  plan$text <- vendor_desc_text(plan$d)
  plan$open_text <- vapply(plan$open, `[[`, "", "plain")
  plan$change_text <- vapply(
    vendor_change_msgs(plan$changes),
    `[[`,
    "",
    "plain"
  )
  plan
}

# The persons of the Authors@R field of a DESCRIPTION text
persons_of <- function(text) {
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  write_text(path, text)
  value <- read.dcf(path, keep.white = "Authors@R")[1, "Authors@R"]
  Encoding(value) <- "UTF-8"
  vendor_read_persons(value)
}

hv <- intToUtf8(229L)

# The Authors@R line of host_desc_minimal
ann_line <- paste0(
  "Authors@R: person(\"Ann\", \"Example\", email = \"ann@example.org\", ",
  "role = c(\"aut\", \"cre\"))"
)

inlavaan_like <- paste0(
  "Type: Package\n",
  "Package: INLAvaan\n",
  "Title: Approximate Bayesian Latent Variable Analysis\n",
  "Version: 0.3.99.9000\n",
  "Authors@R: c(\n",
  "    person(\"Haziq\", \"Jamil\", , \"haziq.jamil@gmail.com\", ",
  "role = c(\"aut\", \"cre\", \"cph\"),\n",
  "           comment = c(ORCID = \"0000-0003-3298-1010\")),\n",
  "    person(\"H",
  hv,
  "vard\", \"Rue\", role = \"ctb\"),\n",
  "    person(\"Thomas\", \"Luu\", role = \"cph\",\n",
  "           comment = \"Skew-normal quantile approximation ",
  "(Apache License 2.0)\")\n",
  "  )\n",
  "Description: Implements approximate Bayesian inference for Structural\n",
  "    Equation Models.\n",
  "License: GPL (>= 3)\n",
  "Depends: \n",
  "    R (>= 4.1.0)\n",
  "Imports:\n",
  "    cli,\n",
  "    graphics,\n",
  "    lavaan (>= 0.7-2),\n",
  "    methods,\n",
  "    stats,\n",
  "    utils\n",
  "Suggests: \n",
  "    blavaan,\n",
  "    ggplot2,\n",
  "    numDeriv,\n",
  "    quarto,\n",
  "    testthat (>= 3.2.0)\n",
  "Encoding: UTF-8\n",
  "LazyData: true\n"
)

old_style <- paste0(
  "Package: oldhost\n",
  "Type: Package\n",
  "Title: An Old Host\n",
  "Version: 1.0\n",
  "Author: Ann Example\n",
  "Maintainer: Ann Example <ann@example.org>\n",
  "Description: An old package.\n",
  "License: GPL-3\n",
  "Depends: R (>= 3.5.0)\n",
  "Imports: cli (>= 3.0.0), stats\n"
)

test_that("an INLAvaan-like DESCRIPTION gains only what it lacks", {
  plan <- plan_for(inlavaan_like)
  expected <- sub(
    "    cli,\n",
    "    cli (>= 3.6.0),\n",
    inlavaan_like,
    fixed = TRUE
  )
  expected <- sub(
    "    blavaan,\n",
    "    blavaan,\n    future,\n",
    expected,
    fixed = TRUE
  )
  expected <- sub(
    "    numDeriv,\n",
    "    numDeriv,\n    parallelly,\n",
    expected,
    fixed = TRUE
  )
  expect_identical(plan$text, expected)
  expect_true(plan$safe)
  expect_identical(plan$open_text, character())
  expect_identical(
    plan$change_text,
    c(
      "Raised cli in Imports to `(>= 3.6.0)`, from no version.",
      "Added future and parallelly to Suggests in 'DESCRIPTION'."
    )
  )
  # The persons are kept as they were, Haziq Jamil and Thomas Luu included
  expect_identical(
    format(persons_of(plan$text)),
    format(persons_of(inlavaan_like))
  )
})

test_that("a minimal DESCRIPTION gains every field, with CRLF kept", {
  plan <- plan_for(host_desc_minimal)
  expect_identical(plan$text, minimal_expected)
  expect_identical(
    plan$change_text,
    c(
      paste(
        "Added the engine's copyright holders to Authors@R in 'DESCRIPTION',",
        "with the role \"cph\": Haziq Jamil and Thomas Luu."
      ),
      "Added cli (>= 3.6.0), stats and utils to Imports in 'DESCRIPTION'.",
      "Added future and parallelly to Suggests in 'DESCRIPTION'.",
      "Added `R (>= 4.1.0)` to Depends."
    )
  )
  persons <- persons_of(plan$text)
  expect_length(persons, 3L)
  expect_identical(format(persons[1]), format(persons_of(host_desc_minimal)))
  expect_identical(
    format(persons[2:3], include = c("given", "family", "role")),
    c("Haziq Jamil [cph]", "Thomas Luu [cph]")
  )
  plan <- plan_for(gsub("\n", "\r\n", host_desc_minimal, fixed = TRUE))
  expect_identical(
    plan$text,
    gsub("\n", "\r\n", minimal_expected, fixed = TRUE)
  )
  # Applied again, nothing changes
  again <- plan_for(minimal_expected)
  expect_identical(again$text, minimal_expected)
  expect_length(again$changes, 0L)
  expect_length(again$open, 0L)
})

test_that("a single person over several lines is wrapped in c()", {
  desc <- sub(
    paste0(ann_line, "\n"),
    paste0(
      "Authors@R: person(\"Ann\", \"Example\", email = \"ann@example.org\",\n",
      "    role = c(\"aut\", \"cre\"))\n"
    ),
    host_desc_minimal,
    fixed = TRUE
  )
  plan <- plan_for(desc)
  expect_match(
    plan$text,
    paste0(
      "Authors@R: c(\n",
      "    person(\"Ann\", \"Example\", email = \"ann@example.org\",\n",
      "        role = c(\"aut\", \"cre\")),\n",
      "    person(\"Haziq\", \"Jamil\", role = \"cph\",\n"
    ),
    fixed = TRUE
  )
  expect_identical(
    format(persons_of(plan$text)[1]),
    format(persons_of(desc))
  )
  # Arguments aligned under the opening parenthesis move left with it
  desc <- sub(
    paste0(ann_line, "\n"),
    paste0(
      "Authors@R: person(\"Ann\", \"Example\", email = \"ann@example.org\",\n",
      "                  role = c(\"aut\", \"cre\"))\n"
    ),
    host_desc_minimal,
    fixed = TRUE
  )
  plan <- plan_for(desc)
  expect_match(
    plan$text,
    paste0(
      "Authors@R: c(\n",
      "    person(\"Ann\", \"Example\", email = \"ann@example.org\",\n",
      "           role = c(\"aut\", \"cre\")),\n",
      "    person(\"Haziq\", \"Jamil\", role = \"cph\",\n"
    ),
    fixed = TRUE
  )
  expect_identical(
    format(persons_of(plan$text)[1]),
    format(persons_of(desc))
  )
  expect_identical(
    vendor_shift_lines(c("                  b = 2,", ")"), 11L),
    c("           b = 2,", ")")
  )
  expect_identical(
    vendor_shift_lines(c("  b = 2,", "                  c = 3)"), 11L),
    c("      b = 2,", "                      c = 3)")
  )
})

test_that("a comment in Authors@R does not take the comma that is added", {
  ann <- "person(\"Ann\", \"B\", role = c(\"aut\", \"cre\"))"
  with_authors <- function(field) {
    sub(
      paste0(ann_line, "\n"),
      field,
      host_desc_minimal,
      fixed = TRUE
    )
  }
  check <- function(field, expected) {
    desc <- with_authors(field)
    plan <- plan_for(desc)
    expect_true(plan$safe)
    expect_match(plan$text, expected, fixed = TRUE)
    persons <- persons_of(plan$text)
    expect_length(persons, 3L)
    expect_identical(format(persons[1]), format(persons_of(desc)))
  }
  # A comment after the last person, and the closing parenthesis on a line
  # of its own.
  check(
    paste0("Authors@R: c(\n    ", ann, " # the maintainer\n  )\n"),
    paste0(
      "    ",
      ann,
      ", # the maintainer\n    person(\"Haziq\", \"Jamil\", role = \"cph\","
    )
  )
  # A line with a comment only, before the closing parenthesis
  check(
    paste0("Authors@R: c(\n    ", ann, "\n    # more to come\n  )\n"),
    paste0(
      "    ",
      ann,
      ",\n    # more to come\n    person(\"Haziq\", \"Jamil\""
    )
  )
  # The closing parenthesis on the line of the last person, and a comment
  check(
    paste0("Authors@R: c(", ann, ") # me\n"),
    paste0("Authors@R: c(", ann, ", # me\n    person(\"Haziq\"")
  )
  # A single person with a comment
  check(
    paste0("Authors@R: ", ann, " # me\n"),
    paste0("Authors@R: c(\n    ", ann, ", # me\n    person(\"Haziq\"")
  )
  # A line with a comment only, after the closing parenthesis
  check(
    paste0("Authors@R: c(\n    ", ann, "\n  )\n  # the end\n"),
    paste0(
      "    ",
      ann,
      ",\n    person(\"Haziq\", \"Jamil\""
    )
  )
  check(
    paste0("Authors@R: c(", ann, ")\n  # the end\n"),
    paste0("Authors@R: c(", ann, ",\n    person(\"Haziq\"")
  )
  # A # inside a string is not a comment
  sharp <- "person(\"Ann\", \"B\", role = \"aut\", comment = \"C# and R\")"
  check(
    paste0("Authors@R: c(\n    ", sharp, "\n  )\n"),
    paste0("    ", sharp, ",\n    person(\"Haziq\"")
  )
})

test_that("an old-style DESCRIPTION is raised in place, authors left to do", {
  plan <- plan_for(old_style)
  expected <- sub("R (>= 3.5.0)", "R (>= 4.1.0)", old_style, fixed = TRUE)
  expected <- sub(
    "Imports: cli (>= 3.0.0), stats\n",
    "Imports: cli (>= 3.6.0), stats, utils\nSuggests: future, parallelly\n",
    expected,
    fixed = TRUE
  )
  expect_identical(plan$text, expected)
  expect_identical(
    plan$open_text,
    paste(
      "'DESCRIPTION' has no Authors@R field: add Haziq Jamil and Thomas Luu to",
      "its Author field, as copyright holders."
    )
  )
  expect_true(
    "Raised R in Depends to `(>= 4.1.0)`, from `(>= 3.5.0)`." %in%
      plan$change_text
  )
  expect_true(
    "Raised cli in Imports to `(>= 3.6.0)`, from `(>= 3.0.0)`." %in%
      plan$change_text
  )
})

test_that("an entry goes last in an unsorted field and on a wrapped line", {
  desc <- sub(
    "License: GPL (>= 3)\n",
    "License: GPL (>= 3)\nImports: utils, methods\n",
    host_desc_minimal,
    fixed = TRUE
  )
  plan <- plan_for(desc)
  expect_match(
    plan$text,
    "Imports: utils, methods, cli (>= 3.6.0), stats\n",
    fixed = TRUE
  )
  desc <- sub(
    "License: GPL (>= 3)\n",
    "License: GPL (>= 3)\nImports: methods, utils,\n    graphics\n",
    host_desc_minimal,
    fixed = TRUE
  )
  plan <- plan_for(desc)
  expect_match(
    plan$text,
    "Imports: methods, utils,\n    graphics, cli (>= 3.6.0), stats\n",
    fixed = TRUE
  )
  # A new field one entry per line when another one is
  desc <- sub(
    "License: GPL (>= 3)\n",
    "License: GPL (>= 3)\nImports:\n    stats,\n    utils\n",
    host_desc_minimal,
    fixed = TRUE
  )
  plan <- plan_for(desc)
  expect_match(
    plan$text,
    paste0(
      "Imports:\n    cli (>= 3.6.0),\n    stats,\n    utils\n",
      "Suggests:\n    future,\n    parallelly\n"
    ),
    fixed = TRUE
  )
  expect_match(plan$text, "Depends:\n    R (>= 4.1.0)\nImports:", fixed = TRUE)
})

test_that("R in Depends is raised, or left alone with an item to do", {
  with_depends <- function(depends) {
    sub(
      "License: GPL (>= 3)\n",
      paste0("License: GPL (>= 3)\nDepends: ", depends, "\n"),
      host_desc_minimal,
      fixed = TRUE
    )
  }
  plan <- plan_for(with_depends("R (>= 3.5)"))
  expect_match(plan$text, "Depends: R (>= 4.1.0)\n", fixed = TRUE)
  plan <- plan_for(with_depends("R"))
  expect_match(plan$text, "Depends: R (>= 4.1.0)\n", fixed = TRUE)
  plan <- plan_for(with_depends("R (> 4.0)"))
  expect_match(plan$text, "Depends: R (> 4.0)\n", fixed = TRUE)
  expect_match(
    plan$open_text[[1]],
    "asks for `R (>= 4.1.0)`: it has `R (> 4.0)`",
    fixed = TRUE
  )
  plan <- plan_for(with_depends("R (>= 4.2.0)"))
  expect_match(plan$text, "Depends: R (>= 4.2.0)\n", fixed = TRUE)
  expect_false(any(grepl("Depends", plan$change_text)))
  plan <- plan_for(with_depends("methods"))
  expect_match(plan$text, "Depends: R (>= 4.1.0), methods\n", fixed = TRUE)
})

test_that("package names are escaped in the patterns that find entries", {
  expect_identical(
    vendor_regex_escape(c("R.utils", "a+b", "x[1]", "a\\b", "c{2}|d")),
    c("R\\.utils", "a\\+b", "x\\[1\\]", "a\\\\b", "c\\{2\\}\\|d")
  )
  for (name in c("R.utils", "a+b", "x[1]", "a\\b")) {
    pattern <- paste0("^", vendor_regex_escape(name), "$")
    expect_true(grepl(pattern, name), label = name)
    expect_true(grepl(pattern, name, perl = TRUE), label = name)
  }
  expect_false(grepl(
    paste0("^", vendor_regex_escape("R.utils"), "$"),
    "RXutils"
  ))
  # An entry whose name only matches the unescaped pattern is not taken
  desc <- sub(
    "License: GPL (>= 3)\n",
    "License: GPL (>= 3)\nImports:\n    cliX,\n    stats\n",
    host_desc_minimal,
    fixed = TRUE
  )
  d <- vendor_desc_read(write_text(tempfile(), desc))
  rows <- which(startsWith(d$lines, "    "))
  expect_identical(vendor_entry_row(d$lines, rows, "stats"), rows[2])
  expect_identical(vendor_entry_row(d$lines, rows, "cl.X"), rows[2])
})

test_that("cli in Suggests only is added to Imports, and left in Suggests", {
  desc <- sub(
    "License: GPL (>= 3)\n",
    "License: GPL (>= 3)\nSuggests: cli, testthat\n",
    host_desc_minimal,
    fixed = TRUE
  )
  plan <- plan_for(desc)
  expect_match(
    plan$text,
    "Imports: cli (>= 3.6.0), stats, utils\n",
    fixed = TRUE
  )
  expect_match(
    plan$text,
    "Suggests: cli, future, parallelly, testthat\n",
    fixed = TRUE
  )
  expect_true(
    "Remove cli from Suggests: it is now in Imports." %in% plan$open_text
  )
})

test_that("a rights holder already listed is not added again", {
  desc <- sub(
    paste0(ann_line, "\n"),
    paste0(
      "Authors@R: c(\n",
      "    person(\"Haziq\", \"Jamil\", role = c(\"aut\", \"cre\")),\n",
      "    person(\"T.\", \"Luu\", role = \"cph\",\n",
      "           comment = c(ORCID = \"0000-0000-0000-0001\")),\n",
      "    person(\"Tom\", \"Luu\", role = \"ctb\")\n",
      "  )\n"
    ),
    host_desc_minimal,
    fixed = TRUE
  )
  holders <- fixture_holders()
  # The same person by ORCID, under another spelling
  holders[2] <- utils::person(
    "Thomas",
    "Luu",
    role = "cph",
    comment = c(ORCID = "0000-0000-0000-0001")
  )
  plan <- plan_for(desc, holders)
  expect_identical(format(persons_of(plan$text)), format(persons_of(desc)))
  expect_false(any(grepl("Authors@R", plan$change_text)))
  expect_identical(
    plan$open_text,
    "Add the role \"cph\" to Haziq Jamil in Authors@R."
  )
})

test_that("an Authors@R that is not plain person() calls is not evaluated", {
  for (value in c(
    "c(person(\"A\", \"B\"), system(\"x\"))",
    "as.person(\"A B\")",
    "authors"
  )) {
    desc <- sub(
      ann_line,
      paste("Authors@R:", value),
      host_desc_minimal,
      fixed = TRUE
    )
    plan <- plan_for(desc)
    expect_match(plan$text, paste("Authors@R:", value), fixed = TRUE)
    expect_match(plan$open_text, "could not be read safely", all = FALSE)
  }
  expect_null(vendor_read_persons("c(person(\"A\", \"B\"), system(\"x\"))"))
  expect_null(vendor_read_persons("utils::as.person(\"A B\")"))
  expect_s3_class(
    vendor_read_persons(
      "c(utils::person(\"A\", \"B\", , \"a@b.org\", role = \"aut\"))"
    ),
    "person"
  )
})

test_that("a DESCRIPTION in latin1 is edited and written back in latin1", {
  rene <- paste0("\"Ren", intToUtf8(233L), "\"")
  as_latin1 <- function(text) {
    text <- sub("Encoding: UTF-8", "Encoding: latin1", text, fixed = TRUE)
    text <- sub("\"Ann\"", rene, text, fixed = TRUE)
    charToRaw(iconv(text, from = "UTF-8", to = "latin1"))
  }
  root <- local_host()
  path <- file.path(root, "DESCRIPTION")
  writeBin(as_latin1(host_desc_minimal), path)
  plan <- vendor_plan_description(root, fixture_holders())
  expect_true(plan$safe)
  expect_identical(plan$d$encoding, "latin1")
  vendor_desc_write(plan$d, path)
  expect_identical(
    readBin(path, "raw", n = file.size(path)),
    as_latin1(minimal_expected)
  )
  persons <- vendor_read_persons(vendor_desc(root)[["Authors@R"]])
  expect_identical(
    vendor_person_name(persons[1]),
    paste0("Ren", intToUtf8(233L), " Example")
  )
  # An encoding that cannot be read stops with a clear error
  bad <- sub(
    "latin1",
    "no-such-encoding",
    rawToChar(as_latin1(host_desc_minimal)),
    useBytes = TRUE
  )
  writeBin(charToRaw(bad), path)
  expect_error(
    vendor_plan_description(root, fixture_holders()),
    "cannot be read in its encoding",
    class = "inlavanilla_vendor_error"
  )
})

test_that("an ORCID that is not valid is read without a warning", {
  text <- paste0(
    "person(\"First\", \"Last\", , \"first.last@example.com\", ",
    "role = c(\"aut\", \"cre\"), comment = c(ORCID = \"YOUR-ORCID-ID\"))"
  )
  expect_no_warning(persons <- vendor_read_persons(text))
  expect_identical(vendor_person_name(persons), "First Last")
  desc <- sub(
    ann_line,
    paste("Authors@R:", text),
    host_desc_minimal,
    fixed = TRUE
  )
  expect_no_warning(plan <- plan_for(desc))
  expect_true(plan$safe)
  expect_length(persons_of(plan$text), 3L)
})

test_that("a holder with a non-ASCII name is written in ASCII", {
  holders <- vendor_holders(paste0(
    "c(person(\"H",
    hv,
    "vard\", \"Rue\", role = c(\"ctb\", \"cph\"), ",
    "comment = c(ORCID = \"0000-0002-0222-1881\")))"
  ))
  plan <- plan_for(host_desc_minimal, holders)
  expect_match(
    plan$text,
    "person(\"H\\u00e5vard\", \"Rue\", role = \"cph\",",
    fixed = TRUE
  )
  expect_true(all(utf8ToInt(plan$text) < 128L))
  persons <- persons_of(plan$text)
  expect_identical(vendor_person_name(persons[2]), paste0("H", hv, "vard Rue"))
})

test_that("a failed safety check leaves DESCRIPTION as it was", {
  local_mocked_bindings(vendor_desc_safe = function(...) "forced")
  plan <- plan_for(host_desc_minimal)
  expect_identical(plan$text, host_desc_minimal)
  expect_identical(plan$safe, "forced")
  expect_length(plan$changes, 0L)
  expect_identical(
    plan$open_text,
    c(
      "Add `R (>= 4.1.0)` to Depends.",
      "Add cli (>= 3.6.0) to Imports.",
      "Add stats to Imports.",
      "Add utils to Imports.",
      "Add future to Suggests.",
      "Add parallelly to Suggests.",
      "Add Haziq Jamil and Thomas Luu to Authors@R, with the role \"cph\"."
    )
  )
  expect_identical(plan$left, c("dependencies", "Authors@R"))
})

test_that("the safety check catches an edit that went wrong", {
  # Authors@R: the role of the old person is changed. The dependencies are
  # still edited.
  real_authors <- vendor_plan_authors
  local_mocked_bindings(
    vendor_plan_authors = function(d, holders) {
      out <- real_authors(d, holders)
      out$d$lines <- sub(
        "role = c(\"aut\", \"cre\")",
        "role = \"aut\"",
        out$d$lines,
        fixed = TRUE
      )
      out
    }
  )
  plan <- plan_for(host_desc_minimal)
  expect_identical(
    plan$safe,
    "Authors@R would not hold the old persons and then the new ones"
  )
  expect_identical(plan$left, "Authors@R")
  expect_identical(
    format(persons_of(plan$text)),
    format(persons_of(host_desc_minimal))
  )
  expect_match(
    plan$text,
    "Imports: cli (>= 3.6.0), stats, utils\n",
    fixed = TRUE
  )
  expect_identical(
    plan$open_text,
    "Add Haziq Jamil and Thomas Luu to Authors@R, with the role \"cph\"."
  )
  # The dependencies: another field changes. Authors@R is still edited.
  local_mocked_bindings(vendor_plan_authors = real_authors)
  real_add <- vendor_add_dep
  local_mocked_bindings(
    vendor_add_dep = function(d, field, pkg, entry, after) {
      d <- real_add(d, field, pkg, entry, after)
      d$lines <- sub("A Host", "A Ghost", d$lines, fixed = TRUE)
      d
    }
  )
  plan <- plan_for(host_desc_minimal)
  expect_identical(plan$safe, "a field that should not change would change")
  expect_identical(plan$left, "dependencies")
  expect_match(plan$text, "Title: A Host of the Engine\n", fixed = TRUE)
  expect_false(grepl("Imports", plan$text, fixed = TRUE))
  expect_length(persons_of(plan$text), 3L)
  expect_identical(
    plan$open_text,
    c(
      "Add `R (>= 4.1.0)` to Depends.",
      "Add cli (>= 3.6.0) to Imports.",
      "Add stats to Imports.",
      "Add utils to Imports.",
      "Add future to Suggests.",
      "Add parallelly to Suggests."
    )
  )
})

test_that("the safety check refuses each kind of bad edit", {
  as_d <- function(text) {
    list(lines = strsplit(text, "\n")[[1]], eol = "\n", final_nl = TRUE)
  }
  old <- as_d(minimal_expected)
  edit <- function(from, to) as_d(sub(from, to, minimal_expected, fixed = TRUE))
  add_x <- list(list(kind = "add", field = "Imports", entry = "x"))
  raise_cli <- list(
    list(kind = "raise", field = "Imports", pkg = "cli", new = "(>= 9)")
  )
  expect_true(vendor_desc_safe(old, old, list(), NULL))
  # A field that the changes do not edit
  expect_identical(
    vendor_desc_safe(
      old,
      edit("Version: 0.1.0", "Version: 0.2.0"),
      add_x,
      NULL
    ),
    "a field that should not change would change"
  )
  expect_identical(
    vendor_desc_safe(
      old,
      edit("stats, utils", "stats, utils, x"),
      list(),
      NULL
    ),
    "a field that should not change would change"
  )
  expect_identical(
    vendor_desc_safe(old, edit("cph", "ctb"), add_x, NULL),
    "a field that should not change would change"
  )
  expect_identical(
    vendor_desc_safe(
      old,
      edit("Encoding", "LazyData: true\nEncoding"),
      add_x,
      NULL
    ),
    "a field would be added that should not be"
  )
  # The entries of an edited field
  expect_true(
    vendor_desc_safe(old, edit("stats, utils", "stats, utils, x"), add_x, NULL)
  )
  expect_identical(
    vendor_desc_safe(old, edit("stats, utils", "x, utils"), add_x, NULL),
    "the entries of Imports would not be as planned"
  )
  expect_true(
    vendor_desc_safe(old, edit("cli (>= 3.6.0)", "cli (>= 9)"), raise_cli, NULL)
  )
  expect_identical(
    vendor_desc_safe(
      old,
      edit("cli (>= 3.6.0)", "cli (>= 3.0)"),
      raise_cli,
      NULL
    ),
    "an entry of Imports would change other than by a raise"
  )
  expect_identical(
    vendor_desc_safe(
      old,
      edit("cli (>= 3.6.0)", "cli (== 9)"),
      raise_cli,
      NULL
    ),
    "an entry of Imports would change other than by a raise"
  )
  # Text that does not read as DCF
  expect_identical(
    vendor_desc_safe(
      old,
      edit("Version: 0.1.0\n", "Version: 0.1.0\n  \n "),
      list(),
      NULL
    ),
    "the new text does not read as a DESCRIPTION file"
  )
})

test_that("the license is graded against GPL (>= 3)", {
  grades <- c(
    "GPL-3" = "pass",
    "GPL (>= 3)" = "pass",
    "GPL (>=3.0)" = "pass",
    "GPL-2 | GPL-3" = "pass",
    "GPL-3 + file LICENSE" = "warn",
    "GPL" = "warn",
    "GPL (>= 2)" = "warn",
    "GPL (>= 2.0)" = "warn",
    "AGPL-3" = "warn",
    "AGPL (>= 3)" = "warn",
    "file LICENSE" = "file",
    "MIT + file LICENSE" = "fail",
    "BSD_3_clause + file LICENSE" = "fail",
    "Apache License 2.0" = "fail",
    "LGPL-3" = "fail",
    "GPL-2" = "fail",
    "CC BY 4.0" = "fail",
    # Other spellings of GPL version 3 that R accepts
    "GPL version 3" = "pass",
    "GNU General Public License version 3" = "pass",
    "GNU General Public License (>= 3)" = "pass",
    "GPL (>= 3.0.0)" = "pass",
    "GPL (> 2)" = "pass",
    "GPL-3 | file LICENSE" = "pass",
    "GNU General Public License" = "warn",
    "GNU Affero General Public License version 3" = "warn",
    "GPL (<= 3)" = "warn",
    "GPL (< 3)" = "fail",
    "GPL version 2" = "fail",
    "GNU Lesser General Public License (>= 3)" = "fail",
    "Unknown License version 3" = "fail"
  )
  for (license in names(grades)) {
    expect_identical(
      vendor_license_grade(license),
      grades[[license]],
      label = license
    )
  }
  expect_identical(vendor_license_grade(NA_character_), "fail")
})

test_that("the dependency and authors checks of vanilla_status() work", {
  root <- local_host(desc = host_desc_minimal)
  check <- vendor_check_deps(vendor_desc(root), character())
  expect_identical(
    vapply(check$fail, `[[`, "", "plain"),
    c(
      "`R (>= 4.1.0)` in Depends",
      "cli (>= 3.6.0) in Imports",
      "stats in Imports",
      "utils in Imports",
      "future in Suggests",
      "parallelly in Suggests"
    )
  )
  root <- local_host(desc = minimal_expected)
  check <- vendor_check_deps(vendor_desc(root), "R/vanilla-a.R")
  expect_length(check$fail, 0L)
  expect_length(check$warn, 0L)
  expect_identical(vendor_authors_row(vendor_desc(root), "v1")$status, "pass")
  root <- local_host(desc = old_style)
  check <- vendor_check_deps(vendor_desc(root), "R/vanilla-a.R")
  expect_identical(
    vapply(check$fail, `[[`, "", "plain"),
    c(
      "`R (>= 4.1.0)` in Depends: it has `R (>= 3.5.0)`",
      "cli (>= 3.6.0) in Imports: it has `cli (>= 3.0.0)`",
      "utils in Imports",
      "future in Suggests",
      "parallelly in Suggests"
    )
  )
  expect_identical(vendor_authors_row(vendor_desc(root), "v1")$status, "warn")
  root <- local_host(desc = host_desc_minimal)
  row <- vendor_authors_row(vendor_desc(root), "v1")
  expect_identical(row$status, "fail")
  expect_identical(
    row$msgs[[1]]$plain,
    paste(
      "Authors@R does not name Haziq Jamil and Thomas Luu, who hold rights in",
      "the engine."
    )
  )
  desc <- paste0(minimal_expected, "Collate: 'hello.R'\n")
  root <- local_host(desc = desc)
  check <- vendor_check_deps(vendor_desc(root), "R/vanilla-a.R")
  expect_identical(
    vapply(check$fail, `[[`, "", "plain"),
    "the engine files in Collate: it lacks 'R/vanilla-a.R'"
  )
})
