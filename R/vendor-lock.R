# Copies of the engine in other packages. use_vanilla() copies the engine files
# R/vanilla-*.R of an inlavanilla release into another package (the host), and
# lists them with their MD5 hashes in the lock inst/inlavanilla.lock, in the
# format that INLAvaan's inst/dev/engine-tools.R reads and writes.
# vanilla_status() checks the copy. The helpers of both are in the files
# R/vendor-*.R, and every one of their names starts with vendor_, so that none
# can clash with a name of the engine.

## ----- Constants -------------------------------------------------------------

vendor_upstream <- "https://github.com/haziqj/inlavanilla"
vendor_issues_url <- "https://github.com/haziqj/inlavanilla/issues"
vendor_lock_file <- file.path("inst", "inlavanilla.lock")
vendor_file_pattern <- "^R/vanilla-[^/]+[.]R$"
vendor_hash <- "md5 of the text with LF line endings"
vendor_header_fields <- c(
  "Upstream",
  "Release",
  "Commit",
  "Version",
  "Synced",
  "Hash"
)
vendor_tag_pattern <- "^v[0-9]+([.-][0-9]+)*$"

# The first release with the 13 engine files. Change it only with the engine.
vendor_first_release <- "v26.11.24"

# Roxygen tags that would change NAMESPACE, or the collation order through
# @include, if a copied file held them.
vendor_roxygen_ns <- paste0(
  "^\\s*#'.*@",
  "(export|import|useDynLib|rawNamespace|S3method|include)"
)

# What the engine needs in the DESCRIPTION of a package that holds it. Change
# it only with the engine: a test in tests/testthat/test-engine-files.R checks
# it against every pkg:: in R/vanilla-*.R.
vendor_needs <- list(
  r = "4.1.0",
  imports = c(cli = "3.6.0", stats = NA, utils = NA),
  suggests = c("future", "parallelly")
)

vendor_test_file <- file.path("tests", "testthat", "test-vanilla-engine.R")
vendor_ci_file <- file.path(".github", "workflows", "vanilla-status.yaml")
vendor_vignette_url <- "https://haziqj.ml/inlavanilla/articles/engine.html"

## ----- The lock --------------------------------------------------------------

# The lock is DCF: a header record, then one record per copied file with the
# fields File and MD5. Returns list(header = <the six fields, named, in
# order>, files = data.frame(file, md5)). A lock that cannot be read stops
# with the class inlavanilla_lock_error, whose `reason` says why in plain
# text. Extra header fields are dropped, so that a later field does not break
# an older reader.
vendor_read_lock <- function(path) {
  if (!file.exists(path)) {
    vendor_lock_error("there is no such file")
  }
  dcf <- tryCatch(
    read.dcf(path),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(dcf)) {
    vendor_lock_error("it is not in DCF format")
  }
  if (nrow(dcf) == 0L) {
    vendor_lock_error("it is empty")
  }
  header <- dcf[1, ]
  header <- header[!is.na(header)]
  missing <- setdiff(vendor_header_fields, names(header))
  if (length(missing) > 0L) {
    vendor_lock_error(paste0(
      "its header lacks the field",
      if (length(missing) > 1L) "s",
      " ",
      paste(missing, collapse = ", ")
    ))
  }
  if (any(c("File", "MD5") %in% names(header))) {
    vendor_lock_error("its header has a File or MD5 field")
  }
  if (!grepl("^[0-9a-f]{40}$", header[["Commit"]])) {
    vendor_lock_error("its Commit is not a SHA of 40 hexadecimal characters")
  }
  if (!identical(unname(header[["Hash"]]), vendor_hash)) {
    vendor_lock_error(paste0("its Hash is not \"", vendor_hash, "\""))
  }
  if (!grepl(vendor_tag_pattern, header[["Release"]])) {
    vendor_lock_error("its Release is not a release tag such as v26.11.24")
  }
  files <- data.frame(file = character(), md5 = character())
  records <- dcf[-1L, , drop = FALSE]
  if (nrow(records) > 0L) {
    if (!all(c("File", "MD5") %in% colnames(records))) {
      vendor_lock_error("a file record lacks File or MD5")
    }
    in_header <- intersect(vendor_header_fields, colnames(records))
    if (any(!is.na(records[, in_header]))) {
      vendor_lock_error("a file record has a header field")
    }
    files <- data.frame(
      file = unname(records[, "File"]),
      md5 = unname(records[, "MD5"])
    )
    if (anyNA(files$file) || anyNA(files$md5)) {
      vendor_lock_error("a file record lacks File or MD5")
    }
    if (!all(grepl(vendor_file_pattern, files$file))) {
      vendor_lock_error("it lists a file that is not R/vanilla-*.R")
    }
    if (!all(grepl("^[0-9a-f]{32}$", files$md5))) {
      vendor_lock_error("it has an MD5 that is not 32 hexadecimal characters")
    }
    if (anyDuplicated(files$file) > 0L) {
      vendor_lock_error("it lists a file twice")
    }
  }
  list(header = vendor_chr(header[vendor_header_fields]), files = files)
}

vendor_lock_error <- function(reason) {
  stop(structure(
    class = c(
      "inlavanilla_lock_error",
      "inlavanilla_error",
      "error",
      "condition"
    ),
    list(message = reason, call = NULL, reason = reason)
  ))
}

# A named character vector with plain values: unname() would drop the names
vendor_chr <- function(x) {
  stats::setNames(as.character(x), names(x))
}

# Written with LF line endings on every platform, the files in radix order
vendor_write_lock <- function(path, header, files) {
  files <- files[order(files$file, method = "radix"), , drop = FALSE]
  lines <- paste0(vendor_header_fields, ": ", header[vendor_header_fields])
  for (i in seq_len(nrow(files))) {
    lines <- c(
      lines,
      "",
      paste0("File: ", files$file[i]),
      paste0("MD5: ", files$md5[i])
    )
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(lines, con, sep = "\n")
  invisible(path)
}

vendor_short_sha <- function(sha) {
  substr(sha, 1L, 7L)
}

## ----- Hashes of LF text -----------------------------------------------------

# The MD5 of each file with every carriage return removed, as the shell step
# `tr -d '\r' < file | md5sum` computes it, so that a checkout with CRLF line
# endings passes. tools::md5sum() reads files only (its `bytes` argument needs
# R 4.5.0), hence the temporary copy.
vendor_md5 <- function(paths) {
  vapply(paths, vendor_md5_one, character(1), USE.NAMES = FALSE)
}

vendor_md5_one <- function(path) {
  bytes <- readBin(path, "raw", n = file.size(path))
  tmp <- tempfile("md5-lf-")
  on.exit(unlink(tmp))
  writeBin(bytes[bytes != as.raw(13L)], tmp)
  unname(tools::md5sum(tmp))
}

## ----- The package -----------------------------------------------------------

# The root of the package that holds `path`: the first folder at or above it
# with a DESCRIPTION.
vendor_root <- function(path) {
  error_class <- "inlavanilla_vendor_error"
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    abort(
      "{.arg path} must be a single string, the root of a package.",
      class = error_class
    )
  }
  if (!file.exists(path)) {
    abort(
      "{.arg path} is {.path {path}}, which does not exist.",
      class = error_class
    )
  }
  dir <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "DESCRIPTION"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      abort(
        c(
          "There is no {.file DESCRIPTION} in {.path {path}} or a folder
           above it.",
          i = "Give the root of a package as {.arg path}."
        ),
        class = error_class
      )
    }
    dir <- parent
  }
}

vendor_desc_fields_read <- c(
  "Package",
  "Version",
  "License",
  "Depends",
  "Imports",
  "Suggests",
  "Collate",
  "Encoding",
  "Authors@R",
  "Author"
)

# The fields of a DESCRIPTION that the checks read, NA where absent
vendor_desc <- function(root) {
  path <- file.path(root, "DESCRIPTION")
  dcf <- tryCatch(
    read.dcf(path, fields = vendor_desc_fields_read, keep.white = "Authors@R"),
    error = function(e) {
      abort(
        "{.file {path}} cannot be read as a {.file DESCRIPTION} file.",
        parent = e
      )
    }
  )
  desc <- vendor_chr(dcf[1L, ])
  # The fields of a DESCRIPTION in another encoding, such as latin1, as UTF-8
  encoding <- desc[["Encoding"]]
  if (!is.na(encoding) && !toupper(encoding) %in% c("UTF-8", "UTF8")) {
    utf8 <- tryCatch(
      iconv(desc, from = encoding, to = "UTF-8"),
      error = function(e) desc
    )
    desc[!is.na(utf8)] <- utf8[!is.na(utf8)]
  }
  utf8 <- !is.na(desc) & validUTF8(desc)
  Encoding(desc[utf8]) <- "UTF-8"
  desc
}

# The name of the package at `root`, from its DESCRIPTION
vendor_package <- function(root, desc) {
  pkg <- unname(desc[["Package"]])
  if (is.na(pkg)) {
    abort(
      "{.path {root}} is not the root of a package: its {.file DESCRIPTION}
       has no {.field Package} field.",
      class = "inlavanilla_vendor_error"
    )
  }
  pkg
}

# The files R/vanilla-*.R of a package, as "R/vanilla-x.R"
vendor_copies <- function(root) {
  files <- list.files(file.path(root, "R"), pattern = "^vanilla-.+[.]R$")
  file.path("R", sort(files, method = "radix"))
}

## ----- The copies against the lock -------------------------------------------

# The problems of the copies: a data frame with a row per problem, `file`,
# `problem` ("missing", "differs" or "unlisted"), the `md5` of the copy and
# the `lock_md5`. No rows when the copies match the lock.
vendor_check_copies <- function(root, lock) {
  listed <- lock$files$file
  present <- file.exists(file.path(root, listed))
  md5 <- rep(NA_character_, length(listed))
  md5[present] <- vendor_md5(file.path(root, listed[present]))
  differs <- present & md5 != lock$files$md5
  unlisted <- setdiff(vendor_copies(root), listed)
  out <- data.frame(
    file = c(listed[!present], listed[differs], unlisted),
    problem = rep(
      c("missing", "differs", "unlisted"),
      c(sum(!present), sum(differs), length(unlisted))
    ),
    md5 = c(
      rep(NA_character_, sum(!present)),
      md5[differs],
      vendor_md5(file.path(root, unlisted))
    ),
    lock_md5 = c(
      lock$files$md5[!present],
      lock$files$md5[differs],
      rep(NA_character_, length(unlisted))
    )
  )
  out[order(out$file, method = "radix"), , drop = FALSE]
}

## ----- Releases --------------------------------------------------------------

# A tag as a version: v26.11.24-1 is 26.11.24.1, above 26.11.24
vendor_tag_version <- function(tag) {
  numeric_version(sub("^v", "", tag))
}

# The tags among `tags` (a named vector of commits) that are releases with the
# engine files: version tags at or above vendor_first_release, oldest first.
vendor_release_tags <- function(tags) {
  tags <- tags[grepl(vendor_tag_pattern, names(tags))]
  if (length(tags) == 0L) {
    return(tags)
  }
  v <- vendor_tag_version(names(tags))
  keep <- v >= vendor_tag_version(vendor_first_release)
  tags <- tags[keep]
  tags[order(v[keep])]
}

# The version tags among `tags` (names) newer than `release`, oldest first
vendor_newer_tags <- function(tags, release) {
  tags <- tags[grepl(vendor_tag_pattern, tags)]
  if (length(tags) == 0L) {
    return(character())
  }
  v <- vendor_tag_version(tags)
  newer <- v > vendor_tag_version(release)
  tags <- tags[newer]
  tags[order(v[newer])]
}

# Wrappers that tests replace
vendor_today <- function() {
  Sys.Date()
}

vendor_installed_version <- function() {
  utils::packageDescription("inlavanilla", fields = "Version")
}
