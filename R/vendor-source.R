# Where the engine files come from: the releases of inlavanilla on GitHub, or a
# local clone of inlavanilla. On GitHub, the tags come from git's smart-HTTP
# ref list and the files from GitHub's archive of a commit, so that neither
# git, a token nor the GitHub API is needed. A clone needs git: the tags come
# from `git for-each-ref`, and the files from `git archive`, which gives the
# same layout as GitHub's archive.

## ----- The network -----------------------------------------------------------

# Downloads `url` to `dest`. TRUE on success, FALSE on any failure (a 404, no
# connection, a proxy, a timeout), and never an error. `mode = "wb"` keeps the
# bytes as they are on Windows. The method is R's default, which follows the
# proxy settings.
vendor_fetch <- function(url, dest, timeout = 60) {
  old <- options(timeout = max(timeout, getOption("timeout")))
  on.exit(options(old))
  ok <- tryCatch(
    suppressWarnings(utils::download.file(
      url,
      dest,
      mode = "wb",
      quiet = TRUE
    )),
    error = function(e) 1L
  )
  identical(as.integer(ok), 0L) && file.exists(dest) && file.size(dest) > 0
}

vendor_refs_url <- function() {
  paste0(vendor_upstream, ".git/info/refs?service=git-upload-pack")
}

vendor_archive_url <- function(commit) {
  paste0(vendor_upstream, "/archive/", commit, ".tar.gz")
}

# The tags of inlavanilla on GitHub, named by tag, with the commit each names,
# or NULL when GitHub cannot be reached.
vendor_github_tags <- function() {
  tmp <- tempfile("inlavanilla-refs-")
  on.exit(unlink(tmp))
  if (!vendor_fetch(vendor_refs_url(), tmp, timeout = 15)) {
    return(NULL)
  }
  vendor_parse_refs(readBin(tmp, "raw", n = file.size(tmp)))
}

# The tags in a ref list of git's smart-HTTP protocol. Each pkt-line starts
# with four hexadecimal digits that give its length, those four included, and
# 0000 is a flush. The first ref line has the capabilities after a NUL byte.
# The peeled entry <tag>^{} of an annotated tag gives the commit of the tag.
# NULL when the bytes are not such a list, such as the HTML page of a portal.
vendor_parse_refs <- function(bytes) {
  lines <- character()
  i <- 1L
  n <- length(bytes)
  if (n < 4L) {
    return(NULL)
  }
  while (i + 3L <= n) {
    len <- strtoi(rawToChar(bytes[i:(i + 3L)]), 16L)
    if (is.na(len) || (len > 0L && len < 4L) || i + len - 1L > n) {
      return(NULL)
    }
    if (len == 0L) {
      i <- i + 4L
      next
    }
    line <- if (len > 4L) bytes[(i + 4L):(i + len - 1L)] else raw()
    nul <- which(line == as.raw(0L))
    if (length(nul) > 0L) {
      line <- line[seq_len(nul[1L] - 1L)]
    }
    lines <- c(lines, sub("\n$", "", rawToChar(line)))
    i <- i + len
  }
  refs <- grep("^[0-9a-f]{40} refs/tags/", lines, value = TRUE)
  sha <- substr(refs, 1L, 40L)
  tag <- substring(refs, 52L)
  peeled <- endsWith(tag, "^{}")
  tags <- stats::setNames(sha[!peeled], tag[!peeled])
  tags[sub("\\^\\{\\}$", "", tag[peeled])] <- sha[peeled]
  tags
}

## ----- Snapshots of a commit -------------------------------------------------

# The DESCRIPTION and the engine files of a commit of inlavanilla, in a
# temporary folder: list(dir, files, desc, stage). `files` are the engine files
# as "R/vanilla-x.R", in radix order. `desc` holds the fields of the
# DESCRIPTION, and `stage` is the folder to delete afterwards. Stops with the
# class inlavanilla_source_error.
vendor_snapshot <- function(commit, release, clone = NULL) {
  stage <- tempfile("inlavanilla-snapshot-")
  dir.create(stage)
  ext <- if (is.null(clone)) ".tar.gz" else ".tar"
  archive <- tempfile("inlavanilla-", fileext = ext)
  on.exit(unlink(archive))
  if (is.null(clone)) {
    if (!vendor_fetch(vendor_archive_url(commit), archive, timeout = 300)) {
      unlink(stage, recursive = TRUE)
      abort(
        c(
          "Could not download release {release} of inlavanilla from GitHub.",
          i = "Check the connection, or copy from a local clone of
               inlavanilla with {.arg clone}."
        ),
        class = c("inlavanilla_source_error", "inlavanilla_vendor_error")
      )
    }
  } else {
    vendor_clone_archive(clone, commit, archive)
  }
  dir <- tryCatch(
    vendor_extract(archive, commit, stage),
    error = function(e) {
      unlink(stage, recursive = TRUE)
      stop(e)
    }
  )
  files <- list.files(file.path(dir, "R"), pattern = "^vanilla-.+[.]R$")
  files <- file.path("R", sort(files, method = "radix"))
  desc <- read.dcf(file.path(dir, "DESCRIPTION"), keep.white = "Authors@R")
  desc <- vendor_chr(desc[1L, ])
  utf8 <- !is.na(desc) & validUTF8(desc)
  Encoding(desc[utf8]) <- "UTF-8"
  if (!identical(unname(desc["Package"]), "inlavanilla")) {
    unlink(stage, recursive = TRUE)
    pkg <- unname(desc["Package"])
    abort(
      "{.path {clone}} is not a clone of inlavanilla: its {.file DESCRIPTION}
       at {release} is that of {.pkg {pkg}}.",
      class = c("inlavanilla_source_error", "inlavanilla_vendor_error")
    )
  }
  list(dir = dir, files = files, desc = desc, stage = stage, commit = commit)
}

# Extracts the DESCRIPTION and the engine files R/vanilla-*.R of an archive
# whose members all sit in the folder inlavanilla-<commit>/, as in GitHub's
# archive of a commit. Only those members are extracted, so that no other
# path, such as ../x, can be written. Returns the folder.
vendor_extract <- function(archive, commit, stage) {
  sha7 <- vendor_short_sha(commit)
  error_class <- c("inlavanilla_source_error", "inlavanilla_vendor_error")
  top <- paste0("inlavanilla-", commit)
  members <- tryCatch(
    suppressWarnings(utils::untar(archive, list = TRUE, tar = "internal")),
    error = function(e) character()
  )
  if (length(members) == 0L || !all(startsWith(members, paste0(top, "/")))) {
    abort("The archive of {sha7} is not that commit.", class = error_class)
  }
  pattern <- paste0("^", top, "/(DESCRIPTION|R/vanilla-[^/]+[.]R)$")
  wanted <- members[grepl(pattern, members)]
  if (!paste0(top, "/DESCRIPTION") %in% wanted) {
    abort(
      "The archive of {sha7} has no {.file DESCRIPTION}.",
      class = error_class
    )
  }
  suppressWarnings(utils::untar(
    archive,
    files = wanted,
    exdir = stage,
    tar = "internal"
  ))
  dir <- file.path(stage, top)
  paths <- file.path(stage, wanted)
  if (!all(file.exists(paths)) || any(nzchar(Sys.readlink(paths)))) {
    abort(
      "The archive of {sha7} could not be extracted as plain files.",
      class = error_class
    )
  }
  dir
}

## ----- A local clone ---------------------------------------------------------

# Runs git in `dir` and returns its output. A failure stops, or returns NULL
# with error = FALSE.
vendor_git <- function(dir, args, error = TRUE) {
  error_class <- c("inlavanilla_source_error", "inlavanilla_vendor_error")
  if (!nzchar(Sys.which("git"))) {
    abort(
      c(
        "{.arg clone} needs git on the PATH.",
        i = "Leave out {.arg clone} to copy from GitHub, which needs no git."
      ),
      class = error_class
    )
  }
  err <- tempfile("git-stderr-")
  on.exit(unlink(err))
  out <- suppressWarnings(system2(
    "git",
    shQuote(c("-C", dir, args)),
    stdout = TRUE,
    stderr = err
  ))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    if (!error) {
      return(NULL)
    }
    command <- paste(c("git", args), collapse = " ")
    message <- paste(readLines(err, warn = FALSE), collapse = " ")
    abort(
      "{.code {command}} failed in {.path {dir}}: {message}",
      class = error_class
    )
  }
  out
}

# The tags of a clone, named by tag, with the commit each names
vendor_clone_tags <- function(clone) {
  if (is.null(vendor_git(clone, c("rev-parse", "--git-dir"), error = FALSE))) {
    abort(
      "{.path {clone}} is not a clone of inlavanilla: it is not a git
       repository.",
      class = c("inlavanilla_source_error", "inlavanilla_vendor_error")
    )
  }
  out <- vendor_git(
    clone,
    c(
      "for-each-ref",
      "--format=%(refname) %(objectname) %(*objectname)",
      "refs/tags"
    )
  )
  out <- out[nzchar(out)]
  parts <- strsplit(out, " ", fixed = TRUE)
  tag <- sub("^refs/tags/", "", vapply(parts, `[`, "", 1L))
  sha <- vapply(
    parts,
    function(p) if (length(p) >= 3L && nzchar(p[3L])) p[3L] else p[2L],
    ""
  )
  stats::setNames(sha, tag)
}

# The DESCRIPTION and R/ of a commit of a clone, as a tar archive with the
# layout of GitHub's archive.
vendor_clone_archive <- function(clone, commit, dest) {
  vendor_git(
    clone,
    c(
      "archive",
      "--format=tar",
      paste0("--prefix=inlavanilla-", commit, "/"),
      "-o",
      dest,
      commit,
      "--",
      "DESCRIPTION",
      "R"
    )
  )
  invisible(dest)
}

# The path of a local clone, checked
vendor_check_clone <- function(clone) {
  if (is.null(clone)) {
    return(NULL)
  }
  if (!is.character(clone) || length(clone) != 1L || is.na(clone)) {
    abort(
      "{.arg clone} must be the path of a local clone of inlavanilla, or
       {.code NULL}.",
      class = "inlavanilla_vendor_error"
    )
  }
  if (!dir.exists(clone)) {
    abort(
      "{.arg clone} is {.path {clone}}, which is not a folder.",
      class = "inlavanilla_vendor_error"
    )
  }
  normalizePath(clone, winslash = "/")
}

## ----- Which release ---------------------------------------------------------

# The tags of the source, and how the release on GitHub was checked: list(tags,
# github, where), where `github` is GitHub's tags when the source is a clone
# (NULL when GitHub cannot be reached) and `where` names the source in
# messages. Without a clone, NULL when GitHub cannot be reached.
vendor_tags <- function(clone) {
  if (is.null(clone)) {
    tags <- vendor_github_tags()
    if (is.null(tags)) {
      return(NULL)
    }
    return(list(tags = tags, github = tags, where = "on GitHub"))
  }
  list(
    tags = vendor_clone_tags(clone),
    github = vendor_github_tags(),
    where = "in the clone"
  )
}

# The release to copy and its commit: list(release, commit, explicit, notes),
# where `notes` are vendor_msg() lines to print. `tags` are the tags of the
# source, `installed` the version of the inlavanilla installed here and `lock`
# the old lock, or NULL.
vendor_resolve_release <- function(
  release,
  tags,
  installed,
  lock,
  clone = NULL
) {
  error_class <- "inlavanilla_vendor_error"
  first <- vendor_first_release
  releases <- vendor_release_tags(tags)
  newest <- if (length(releases) > 0L) names(releases)[length(releases)] else NA
  where <- if (is.null(clone)) "on GitHub" else "in the clone"
  notes <- list()
  explicit <- !is.null(release)
  if (explicit) {
    given <- release
    if (!startsWith(release, "v")) {
      release <- paste0("v", release)
    }
    if (!grepl(vendor_tag_pattern, release)) {
      abort(
        "{.val {given}} is not a release of inlavanilla: give a tag such as
         {.val {first}}.",
        class = error_class
      )
    }
    if (vendor_tag_version(release) < vendor_tag_version(first)) {
      abort(
        "{release} is older than {first}, the first release with the engine
         files.",
        class = error_class
      )
    }
  } else {
    release <- paste0("v", installed)
    ok <- grepl(vendor_tag_pattern, release) &&
      vendor_tag_version(release) >= vendor_tag_version(first)
    if (!ok) {
      abort(
        c(
          "inlavanilla {installed}, installed here, is older than {first}, the
           first release with the engine files.",
          i = "Install a newer inlavanilla, or give {.arg release}."
        ),
        class = error_class
      )
    }
    if (!release %in% names(releases)) {
      older <- releases[
        vendor_tag_version(names(releases)) <= vendor_tag_version(release)
      ]
      if (length(older) == 0L) {
        abort(
          c(
            "inlavanilla {installed}, installed here, is not a release
             {where}, and no release before it is.",
            i = "Install a newer inlavanilla, or give {.arg release}."
          ),
          class = error_class
        )
      }
      release <- names(older)[length(older)]
      notes <- c(
        notes,
        list(vendor_msg(
          "i",
          "inlavanilla {installed}, installed here, is not a release {where}, so
         this copies {release}, the newest release before it."
        ))
      )
    }
  }
  if (!release %in% names(tags)) {
    if (is.null(clone)) {
      abort(
        c(
          "inlavanilla has no release {.val {release}} on GitHub.",
          i = if (!is.na(newest)) "The newest release is {.val {newest}}."
        ),
        class = error_class
      )
    }
    abort(
      c(
        "There is no release {.val {release}} in the clone {.path {clone}}.",
        i = "Run {.code git fetch --tags} in the clone.",
        i = if (!is.na(newest)) "Its newest release is {.val {newest}}."
      ),
      class = error_class
    )
  }
  if (!is.null(lock)) {
    old <- lock$header[["Release"]]
    if (vendor_tag_version(release) < vendor_tag_version(old)) {
      if (!explicit) {
        abort(
          c(
            "The copy is at {old}, newer than {release}, the release of the
             inlavanilla installed here.",
            i = "Give {.code release = \"{release}\"} to go back, or install a
                 newer inlavanilla."
          ),
          class = error_class
        )
      }
      notes <- c(
        notes,
        list(vendor_msg(
          "!",
          "Going back from {old} to {release}."
        ))
      )
    }
  }
  list(
    release = release,
    commit = unname(tags[[release]]),
    explicit = explicit,
    notes = notes
  )
}

# How a release of a clone stands on GitHub: "same", "moved" (another commit),
# "missing" (no such tag) or "offline".
vendor_github_state <- function(release, commit, github) {
  if (is.null(github)) {
    return("offline")
  }
  if (!release %in% names(github)) {
    return("missing")
  }
  if (!identical(unname(github[[release]]), commit)) "moved" else "same"
}
