# The DESCRIPTION of a package that holds the engine. R's write.dcf() and
# format() of persons reflow fields, drop trailing spaces, turn CRLF into LF,
# collapse dependencies written one per line, rewrite every person and can
# corrupt non-ASCII text, so use_vanilla() changes only the lines of the fields
# it edits, and keeps the file's line endings.

## ----- Reading and writing ---------------------------------------------------

# The lines of a DESCRIPTION, with its line ending, its final newline and its
# encoding. A file that is not valid UTF-8 is read in the encoding of its
# Encoding field, such as latin1, edited as UTF-8 and written back in its own
# encoding (`encoding`, NULL for UTF-8).
vendor_desc_read <- function(path) {
  bytes <- readBin(path, "raw", n = file.size(path))
  eol <- if (any(bytes == as.raw(13L))) "\r\n" else "\n"
  text <- rawToChar(bytes[bytes != as.raw(13L)])
  encoding <- NULL
  if (validUTF8(text)) {
    Encoding(text) <- "UTF-8"
  } else {
    encoding <- vendor_declared_encoding(text)
    utf8 <- tryCatch(
      iconv(text, from = encoding, to = "UTF-8"),
      error = function(e) NA_character_
    )
    if (is.na(utf8)) {
      abort(
        "{.file {path}} is not valid UTF-8, and cannot be read in its
         encoding, {.val {encoding}}.",
        class = "inlavanilla_vendor_error"
      )
    }
    text <- utf8
  }
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  list(
    lines = lines,
    eol = eol,
    final_nl = endsWith(text, "\n"),
    encoding = encoding
  )
}

# The Encoding field of a text that is not valid UTF-8, read byte by byte, or
# latin1 when there is none, in which any bytes can be read.
vendor_declared_encoding <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE, useBytes = TRUE)[[1L]]
  field <- lines[grepl("^Encoding:", lines, useBytes = TRUE)]
  if (length(field) == 0L) {
    return("latin1")
  }
  value <- sub("^Encoding:[[:space:]]*", "", field[1L], useBytes = TRUE)
  value <- sub("[[:space:]]*$", "", value, useBytes = TRUE)
  if (nzchar(value)) value else "latin1"
}

vendor_desc_text <- function(d) {
  text <- paste(d$lines, collapse = d$eol)
  if (d$final_nl) paste0(text, d$eol) else text
}

vendor_desc_write <- function(d, path) {
  text <- vendor_desc_text(d)
  if (!is.null(d$encoding)) {
    text <- iconv(text, from = "UTF-8", to = d$encoding)
    if (is.na(text)) {
      abort(
        "The new text of {.file {path}} cannot be written in its encoding,
         {.val {d$encoding}}.",
        class = "inlavanilla_vendor_error"
      )
    }
  }
  writeBin(charToRaw(text), path)
}

# The first and last line of each field. A field starts at a line with a name
# and a colon. Trailing blank lines belong to no field.
vendor_desc_fields <- function(lines) {
  start <- grep("^[^[:space:]:]+:", lines, useBytes = TRUE)
  end <- c(start[-1L] - 1L, length(lines))
  for (i in seq_along(end)) {
    while (end[i] > start[i] && !nzchar(trimws(lines[end[i]]))) {
      end[i] <- end[i] - 1L
    }
  }
  data.frame(
    field = sub(":.*$", "", lines[start], useBytes = TRUE),
    start = start,
    end = end
  )
}

vendor_field_value <- function(lines, f) {
  v <- lines[f$start:f$end]
  v[1L] <- sub("^[^:]*:", "", v[1L], useBytes = TRUE)
  trimws(paste(v, collapse = "\n"))
}

## ----- Dependencies ----------------------------------------------------------

# The entries of a dependency field: name, operator, version and the text of
# each entry
vendor_parse_deps <- function(value) {
  if (is.null(value) || is.na(value)) {
    value <- ""
  }
  entries <- trimws(strsplit(gsub("\n", " ", value), ",", fixed = TRUE)[[1L]])
  entries <- entries[nzchar(entries)]
  has <- grepl("(", entries, fixed = TRUE)
  data.frame(
    name = sub("^([^[:space:](]+).*$", "\\1", entries),
    op = ifelse(
      has,
      sub("^.*\\([[:space:]]*([<>=]+).*$", "\\1", entries),
      NA_character_
    ),
    version = ifelse(
      has,
      sub(
        "^.*\\([[:space:]]*[<>=]+[[:space:]]*([^)[:space:]]+).*$",
        "\\1",
        entries
      ),
      NA_character_
    ),
    entry = entries
  )
}

# Whether an entry asks for at least version `min`. NA when its operator says
# nothing either way, such as == 3.0 or < 4.
vendor_dep_ok <- function(op, version, min) {
  if (is.na(min)) {
    return(TRUE)
  }
  if (is.na(op)) {
    return(FALSE)
  }
  v <- tryCatch(numeric_version(version), error = function(e) NULL)
  if (is.null(v)) {
    return(NA)
  }
  if (op %in% c(">=", ">", "==") && v >= numeric_version(min)) {
    return(TRUE)
  }
  if (op == ">=") FALSE else NA
}

vendor_dep_text <- function(pkg, version) {
  if (is.na(version)) pkg else paste0(pkg, " (>= ", version, ")")
}

vendor_field_entries <- function(d, field) {
  fields <- vendor_desc_fields(d$lines)
  i <- match(field, fields$field)
  if (is.na(i)) {
    return(vendor_parse_deps(""))
  }
  vendor_parse_deps(vendor_field_value(d$lines, fields[i, ]))
}

# Adds `entry` (for package `pkg`) to the dependency field `field`, where it
# is not. Where the field is missing, it goes after the last present of the
# fields `after`, else at the end of the file.
vendor_add_dep <- function(d, field, pkg, entry, after) {
  lines <- d$lines
  fields <- vendor_desc_fields(lines)
  i <- match(field, fields$field)
  dep_fields <- c("Depends", "Imports", "Suggests", "LinkingTo", "Enhances")
  if (is.na(i)) {
    anchor <- fields$end[fields$field %in% after]
    at <- if (length(anchor) > 0L) max(anchor) else max(fields$end)
    per_line <- any(vapply(
      which(fields$field %in% dep_fields),
      function(k) vendor_one_per_line(lines, fields[k, ]),
      NA
    ))
    add <- if (per_line) {
      c(paste0(field, ":"), paste0("    ", entry))
    } else {
      paste0(field, ": ", entry)
    }
    d$lines <- append(lines, add, after = at)
    return(d)
  }
  f <- fields[i, ]
  deps <- vendor_parse_deps(vendor_field_value(lines, f))
  sorted <- !is.unsorted(tolower(deps$name))
  pos <- if (sorted) sum(tolower(deps$name) < tolower(pkg)) else nrow(deps)
  if (vendor_one_per_line(lines, f)) {
    rows <- (f$start + 1L):f$end
    indent <- sub("^([[:space:]]*).*$", "\\1", lines[rows[1L]])
    if (pos == nrow(deps)) {
      last <- rows[length(rows)]
      if (!grepl(",[[:space:]]*$", lines[last])) {
        lines[last] <- sub("[[:space:]]*$", ",", lines[last])
      }
      lines <- append(lines, paste0(indent, entry), after = last)
    } else {
      # Before the line of the entry at `pos + 1`
      at <- vendor_entry_row(lines, rows, deps$name[pos + 1L])
      lines <- append(lines, paste0(indent, entry, ","), after = at - 1L)
    }
  } else if (f$end == f$start) {
    entries <- append(deps$entry, entry, after = pos)
    lines[f$start] <- paste0(field, ": ", paste(entries, collapse = ", "))
  } else {
    lines[f$end] <- paste0(
      sub("[[:space:],]*$", "", lines[f$end]),
      ", ",
      entry
    )
  }
  d$lines <- lines
  d
}

# Whether a field has one entry per line: nothing after the colon, and no
# line with more than one comma
vendor_one_per_line <- function(lines, f) {
  if (f$end <= f$start) {
    return(FALSE)
  }
  body <- lines[(f$start + 1L):f$end]
  first_empty <- grepl("^[^:]*:[[:space:]]*$", lines[f$start], useBytes = TRUE)
  first_empty &&
    all(lengths(regmatches(body, gregexpr(",", body, fixed = TRUE))) <= 1L)
}

vendor_entry_row <- function(lines, rows, name) {
  pat <- paste0("^[[:space:]]*", vendor_regex_escape(name), "([[:space:](,]|$)")
  hit <- rows[grepl(pat, lines[rows])]
  if (length(hit) > 0L) hit[1L] else rows[length(rows)]
}

vendor_regex_escape <- function(x) {
  gsub("([.|()^${}+*?\\[\\]\\\\])", "\\\\\\1", x)
}

# Replaces the entry of `pkg` in `field` by `entry`, on its own line. NULL
# when the entry cannot be found on one line.
vendor_raise_dep <- function(d, field, pkg, entry) {
  fields <- vendor_desc_fields(d$lines)
  f <- fields[match(field, fields$field), ]
  pat <- paste0(
    "(^",
    field,
    ":[[:space:]]*|^[[:space:]]+|,[[:space:]]*)",
    vendor_regex_escape(pkg),
    "([[:space:]]*\\([^)]*\\))?(?=[[:space:]]*(,|$))"
  )
  for (j in f$start:f$end) {
    if (grepl(pat, d$lines[j], perl = TRUE)) {
      d$lines[j] <- sub(pat, paste0("\\1", entry), d$lines[j], perl = TRUE)
      return(d)
    }
  }
  NULL
}

## ----- Authors@R -------------------------------------------------------------

# The persons of an Authors@R field, evaluated only when the field is made of
# c(), list(), person() and utils::person() calls of constants. NULL when it
# is anything else, or does not evaluate to persons.
vendor_read_persons <- function(text) {
  if (is.null(text) || is.na(text)) {
    return(NULL)
  }
  encoding <- if (validUTF8(text)) "UTF-8" else "unknown"
  expr <- tryCatch(
    parse(text = text, keep.source = FALSE, encoding = encoding),
    error = function(e) NULL
  )
  if (length(expr) != 1L || !vendor_safe_call(expr[[1L]])) {
    return(NULL)
  }
  env <- new.env(parent = baseenv())
  env$person <- utils::person
  out <- tryCatch(eval(expr[[1L]], env), error = function(e) NULL)
  if (inherits(out, "person")) out else NULL
}

vendor_safe_call <- function(e) {
  if (is.null(e) || (is.atomic(e) && length(e) == 1L)) {
    return(is.null(e) || is.character(e) || is.numeric(e) || is.logical(e))
  }
  if (is.name(e)) {
    # An empty argument, as in person("A", "B", , "email")
    return(identical(e, quote(expr = )))
  }
  if (!is.call(e)) {
    return(FALSE)
  }
  f <- e[[1L]]
  ok <- (is.name(f) && as.character(f) %in% c("c", "list", "person")) ||
    identical(f, quote(utils::person))
  ok && all(vapply(as.list(e)[-1L], vendor_safe_call, NA))
}

vendor_person_part <- function(p, part) {
  unclass(p)[[1L]][[part]]
}

vendor_person_name <- function(p) {
  trimws(paste(
    c(vendor_person_part(p, "given"), vendor_person_part(p, "family")),
    collapse = " "
  ))
}

vendor_names_of <- function(persons) {
  vapply(seq_along(persons), function(i) vendor_person_name(persons[i]), "")
}

vendor_orcid <- function(p) {
  cm <- vendor_person_part(p, "comment")
  if (!is.null(names(cm)) && "ORCID" %in% names(cm)) {
    trimws(cm[["ORCID"]])
  } else {
    NA_character_
  }
}

# The same person: the ORCIDs match, or the given and family names match,
# ignoring case and surrounding spaces
vendor_same_person <- function(a, b) {
  oa <- vendor_orcid(a)
  ob <- vendor_orcid(b)
  if (!is.na(oa) && !is.na(ob) && identical(oa, ob)) {
    return(TRUE)
  }
  norm <- function(x) tolower(trimws(paste(x, collapse = " ")))
  identical(
    norm(vendor_person_part(a, "given")),
    norm(vendor_person_part(b, "given"))
  ) &&
    identical(
      norm(vendor_person_part(a, "family")),
      norm(vendor_person_part(b, "family"))
    )
}

# The persons whose roles include "cph", from the text of an Authors@R field,
# or NULL when it cannot be read
vendor_holders <- function(text) {
  persons <- vendor_read_persons(text)
  if (is.null(persons)) {
    return(NULL)
  }
  cph <- vapply(
    seq_along(persons),
    function(i) "cph" %in% vendor_person_part(persons[i], "role"),
    NA
  )
  persons[cph]
}

# Where each holder stands in `persons`: "listed" with the role "cph",
# "no_cph" when listed without it, or "absent"
vendor_holder_state <- function(holders, persons) {
  vapply(
    seq_along(holders),
    function(i) {
      for (k in seq_along(persons)) {
        if (vendor_same_person(holders[i], persons[k])) {
          role <- vendor_person_part(persons[k], "role")
          return(if ("cph" %in% role) "listed" else "no_cph")
        }
      }
      "absent"
    },
    ""
  )
}

# A string as R code, ASCII only: other characters as \uXXXX
vendor_r_string <- function(x) {
  chars <- strsplit(enc2utf8(x), "")[[1L]]
  code <- vapply(chars, utf8ToInt, 1L, USE.NAMES = FALSE)
  chars[chars == "\\"] <- "\\\\"
  chars[chars == "\""] <- "\\\""
  chars[code > 127L] <- sprintf("\\u%04x", code[code > 127L])
  paste0("\"", paste(chars, collapse = ""), "\"")
}

# A rights holder as lines of Authors@R, in usethis's layout: name and role,
# then the comment. The ORCID is copied, never the email. The comment is the
# holder's own text, else a note on the files.
vendor_person_lines <- function(p, indent = "    ") {
  head <- paste0(
    "person(",
    vendor_r_string(paste(vendor_person_part(p, "given"), collapse = " ")),
    ", ",
    vendor_r_string(paste(vendor_person_part(p, "family"), collapse = " ")),
    ", role = \"cph\""
  )
  cm <- vendor_person_part(p, "comment")
  nm <- names(cm)
  if (is.null(nm)) {
    nm <- rep("", length(cm))
  }
  text <- unname(cm[nm == ""])
  if (length(text) == 0L) {
    text <- "Code in R/vanilla-*.R, from inlavanilla"
  }
  orcid <- vendor_orcid(p)
  parts <- c(
    if (!is.na(orcid)) paste0("ORCID = ", vendor_r_string(orcid)),
    vapply(text, vendor_r_string, "", USE.NAMES = FALSE)
  )
  comment <- if (length(parts) == 1L && is.na(orcid)) {
    parts
  } else {
    paste0("c(", paste(parts, collapse = ", "), ")")
  }
  c(
    paste0(indent, head, ","),
    paste0(indent, strrep(" ", 7L), "comment = ", comment, ")")
  )
}

## ----- The plan of the edits -------------------------------------------------

# The edits of a DESCRIPTION that the engine needs, in this order: R in
# Depends, the Imports, the Suggests, and the rights holders in Authors@R.
# Nothing is written here. Returns list(d, changes, open, safe, left): the new
# lines, the changes made (a list of records), the items left for the user to
# do (a list of vendor_msg()), TRUE when every edit passed the safety checks,
# else the first reason, and the parts left as they were, "dependencies" or
# "Authors@R". The edits of the dependencies and the edit of Authors@R are
# checked apart. The edits of a part that fails are not made, and become
# items left to do.
vendor_plan_description <- function(root, holders) {
  path <- file.path(root, "DESCRIPTION")
  old <- vendor_desc_read(path)
  d <- old
  changes <- list()
  open <- list()
  note <- function(change) {
    changes[[length(changes) + 1L]] <<- change
  }
  todo <- function(template, .envir = parent.frame()) {
    open[[length(open) + 1L]] <<- vendor_msg("_", template, .envir)
  }

  # R in Depends
  r_min <- vendor_needs$r
  depends <- vendor_field_entries(d, "Depends")
  k <- match("R", depends$name)
  if (is.na(k)) {
    d <- vendor_add_r(d, r_min)
    note(list(kind = "r_add"))
  } else {
    ok <- vendor_dep_ok(depends$op[k], depends$version[k], r_min)
    old_r <- if (is.na(depends$op[k])) {
      NA_character_
    } else {
      paste0("(", depends$op[k], " ", depends$version[k], ")")
    }
    if (isFALSE(ok)) {
      new_d <- vendor_raise_dep(d, "Depends", "R", vendor_dep_text("R", r_min))
      if (is.null(new_d)) {
        todo("Raise R in {.field Depends} to {.code (>= {r_min})}.")
      } else {
        d <- new_d
        note(list(kind = "r_raise", old = old_r))
      }
    } else if (is.na(ok)) {
      todo(
        "Make sure that {.field Depends} asks for {.code R (>= {r_min})}: it
         has {.code R {old_r}}, which {.fn use_vanilla} does not change."
      )
    }
  }

  # Imports
  imports_after <- c("Depends", "License", "URL", "BugReports")
  for (pkg in names(vendor_needs$imports)) {
    min <- vendor_needs$imports[[pkg]]
    entry <- vendor_dep_text(pkg, min)
    found <- FALSE
    for (field in c("Depends", "Imports")) {
      deps <- vendor_field_entries(d, field)
      k <- match(pkg, deps$name)
      if (is.na(k)) {
        next
      }
      found <- TRUE
      ok <- vendor_dep_ok(deps$op[k], deps$version[k], min)
      old_v <- if (is.na(deps$op[k])) {
        NA_character_
      } else {
        paste0("(", deps$op[k], " ", deps$version[k], ")")
      }
      if (isFALSE(ok)) {
        new_d <- vendor_raise_dep(d, field, pkg, entry)
        if (is.null(new_d)) {
          todo("Raise {.pkg {pkg}} in {.field {field}} to {.code {entry}}.")
        } else {
          d <- new_d
          note(list(
            kind = "raise",
            field = field,
            pkg = pkg,
            old = old_v,
            new = paste0("(>= ", min, ")")
          ))
        }
      } else if (is.na(ok)) {
        todo(
          "Make sure that {.field {field}} asks for {.code {entry}}: it has
           {.code {deps$entry[k]}}, which {.fn use_vanilla} does not change."
        )
      }
      break
    }
    if (found) {
      next
    }
    d <- vendor_add_dep(d, "Imports", pkg, entry, imports_after)
    note(list(kind = "add", field = "Imports", entry = entry))
    if (pkg %in% vendor_field_entries(d, "Suggests")$name) {
      todo(
        "Remove {.pkg {pkg}} from {.field Suggests}: it is now in
         {.field Imports}."
      )
    }
  }

  # Suggests
  suggests_after <- c("Imports", "Depends", "License", "URL", "BugReports")
  for (pkg in vendor_needs$suggests) {
    present <- vapply(
      c("Depends", "Imports", "Suggests"),
      function(field) pkg %in% vendor_field_entries(d, field)$name,
      NA
    )
    if (!any(present)) {
      d <- vendor_add_dep(d, "Suggests", pkg, pkg, suggests_after)
      note(list(kind = "add", field = "Suggests", entry = pkg))
    }
  }

  # The edits of the dependencies are checked apart from the edit of
  # Authors@R, so that a problem with one does not stop the other.
  safe <- TRUE
  left <- character()
  if (length(changes) > 0L) {
    ok <- vendor_desc_safe(old, d, changes, NULL)
    if (!isTRUE(ok)) {
      safe <- ok
      left <- "dependencies"
      d <- old
      open <- c(lapply(changes, vendor_change_todo), open)
      changes <- list()
    }
  }

  # Authors@R
  authors <- vendor_plan_authors(d, holders)
  open <- c(open, authors$open)
  if (!is.null(authors$change)) {
    ok <- vendor_desc_safe(d, authors$d, list(), authors$added)
    if (isTRUE(ok)) {
      d <- authors$d
      note(authors$change)
    } else {
      if (isTRUE(safe)) {
        safe <- ok
      }
      left <- c(left, "Authors@R")
      open <- c(open, list(vendor_change_todo(authors$change)))
    }
  }
  list(d = d, changes = changes, open = open, safe = safe, left = left)
}

# R (>= r_min) as the first entry of Depends, or a new Depends field
vendor_add_r <- function(d, r_min) {
  entry <- vendor_dep_text("R", r_min)
  fields <- vendor_desc_fields(d$lines)
  i <- match("Depends", fields$field)
  if (is.na(i)) {
    return(vendor_add_dep(
      d,
      "Depends",
      "R",
      entry,
      c("License", "URL", "BugReports")
    ))
  }
  f <- fields[i, ]
  first <- d$lines[f$start]
  if (grepl("^Depends:[[:space:]]*$", first) && f$end > f$start) {
    indent <- sub("^([[:space:]]*).*$", "\\1", d$lines[f$start + 1L])
    d$lines <- append(d$lines, paste0(indent, entry, ","), after = f$start)
  } else {
    d$lines[f$start] <- sub(
      "^Depends:[[:space:]]*",
      paste0("Depends: ", entry, ", "),
      first
    )
  }
  d
}

# The rights holders that Authors@R lacks, added after its last person.
# Returns list(d, change, open, added).
vendor_plan_authors <- function(d, holders) {
  out <- list(d = d, change = NULL, open = list(), added = NULL)
  todo <- function(template, .envir = parent.frame()) {
    out$open[[length(out$open) + 1L]] <<- vendor_msg("_", template, .envir)
  }
  if (is.null(holders)) {
    todo(
      "The {.field Authors@R} of inlavanilla could not be read, so its
       copyright holders were not added: add them to {.field Authors@R}
       with the role {.val cph}."
    )
    return(out)
  }
  if (length(holders) == 0L) {
    return(out)
  }
  fields <- vendor_desc_fields(d$lines)
  i <- match("Authors@R", fields$field)
  all_names <- vendor_and(vendor_names_of(holders))
  if (is.na(i)) {
    todo(
      "{.file DESCRIPTION} has no {.field Authors@R} field: add {all_names}
       to its {.field Author} field, as copyright holders."
    )
    return(out)
  }
  f <- fields[i, ]
  value <- vendor_field_value(d$lines, f)
  persons <- vendor_read_persons(value)
  if (is.null(persons)) {
    todo(
      "{.field Authors@R} could not be read safely: add {all_names} to it by
       hand, with the role {.val cph}."
    )
    return(out)
  }
  state <- vendor_holder_state(holders, persons)
  for (k in which(state == "no_cph")) {
    name <- vendor_person_name(holders[k])
    todo("Add the role {.val cph} to {name} in {.field Authors@R}.")
  }
  add <- holders[state == "absent"]
  if (length(add) == 0L) {
    return(out)
  }
  add_names <- vendor_and(vendor_names_of(add))
  expr <- parse(text = value, keep.source = FALSE)[[1L]]
  form <- if (is.call(expr)) vendor_call_symbol(expr) else ""
  lines <- d$lines
  body <- lines[f$start:f$end]
  starts <- grep("^[[:space:]]*(utils::)?person\\(", body, value = TRUE)
  indent <- if (form == "c" && length(starts) > 0L) {
    sub("^([[:space:]]*).*$", "\\1", starts[length(starts)])
  } else {
    "    "
  }
  if (!nzchar(indent)) {
    indent <- "    "
  }
  new <- lapply(seq_along(add), function(k) vendor_person_lines(add[k], indent))
  for (k in seq_along(new)[-length(new)]) {
    n <- length(new[[k]])
    new[[k]][n] <- paste0(new[[k]][n], ",")
  }
  new <- unlist(new)
  if (form == "c") {
    # The closing parenthesis of c() is on the last line that has code. Lines
    # after it hold comments only.
    rows <- f$start:f$end
    ends <- vendor_code_ends(lines[rows])
    last <- max(which(ends > 0L))
    close <- rows[last]
    code <- substr(lines[close], 1L, ends[last])
    if (grepl("^[[:space:]]*\\)$", code)) {
      # The closing parenthesis on a line of its own: a comma after the code
      # of the line before it that has code, before any comment there, and
      # the new persons just before the closing line.
      k <- max(which(ends[seq_len(last - 1L)] > 0L))
      lines[rows[k]] <- vendor_add_comma(lines[rows[k]], ends[k])
      lines <- append(lines, new, after = close - 1L)
    } else if (endsWith(code, ")")) {
      # The closing parenthesis on the line of the last person: it moves
      # after the new persons, and a comment on that line stays there.
      n <- nchar(code)
      head <- sub("[[:space:]]*$", "", substr(code, 1L, n - 1L))
      lines[close] <- paste0(head, ",", substring(lines[close], n + 1L))
      new[length(new)] <- paste0(new[length(new)], ")")
      lines <- append(lines, new, after = close)
    } else {
      todo(
        "{.field Authors@R} does not end with the closing parenthesis of
         {.code c()}: add {add_names} to it by hand, with the role {.val cph}."
      )
      return(out)
    }
  } else if (form %in% c("person", "utils::person")) {
    first <- lines[f$start]
    if (grepl("^Authors@R:[[:space:]]*$", first)) {
      lines[f$start] <- "Authors@R: c("
    } else {
      # The person moves to a line of its own, at column 4
      prefix <- sub("^(Authors@R:[[:space:]]*).*$", "\\1", first)
      lines[f$start] <- paste0("    ", substring(first, nchar(prefix) + 1L))
      if (f$end > f$start) {
        rest <- (f$start + 1L):f$end
        lines[rest] <- vendor_shift_lines(lines[rest], nchar(prefix))
      }
      lines <- append(lines, "Authors@R: c(", after = f$start - 1L)
      f$end <- f$end + 1L
    }
    rows <- f$start:f$end
    ends <- vendor_code_ends(lines[rows])
    k <- max(which(ends > 0L))
    lines[rows[k]] <- vendor_add_comma(lines[rows[k]], ends[k])
    lines <- append(lines, c(new, "  )"), after = f$end)
  } else {
    todo(
      "{.field Authors@R} is not a list of {.code person()} calls: add
       {add_names} to it by hand, with the role {.val cph}."
    )
    return(out)
  }
  out$d$lines <- lines
  out$change <- list(kind = "authors", names = vendor_names_of(add))
  out$added <- add
  out
}

# The other lines of a single person whose first line moves from column `from`
# to column 4, when it is wrapped in c(). Lines aligned at or beyond `from`,
# such as arguments under the opening parenthesis, move left with the person.
# Lines with a hanging indent keep their place relative to the start of the
# line, which moves 4 columns right. Blank lines and lines that start with a
# closing parenthesis do not decide which, and a line is never moved left by
# more than its leading spaces.
vendor_shift_lines <- function(rest, from) {
  spaces <- nchar(sub("^( *).*$", "\\1", rest))
  judge <- nzchar(trimws(rest)) & !grepl("^[[:space:]]*\\)", rest)
  if (!all(spaces[judge] >= from)) {
    return(paste0("    ", rest))
  }
  drop <- from - 4L
  ifelse(spaces >= drop, substring(rest, drop + 1L), rest)
}

# Where the R code of each line ends: the number of characters up to the last
# one that is code, with any comment and the spaces before it left out. A #
# inside a string does not start a comment, and a string can run on over
# several lines, so the lines are read in order from the start of the field.
vendor_code_ends <- function(lines) {
  ends <- integer(length(lines))
  quote <- ""
  escaped <- FALSE
  for (i in seq_along(lines)) {
    chars <- strsplit(lines[i], "", fixed = TRUE)[[1L]]
    end <- 0L
    for (k in seq_along(chars)) {
      ch <- chars[k]
      if (nzchar(quote)) {
        if (escaped) {
          escaped <- FALSE
        } else if (ch == "\\") {
          escaped <- TRUE
        } else if (ch == quote) {
          quote <- ""
        }
        end <- k
      } else if (ch == "#") {
        break
      } else {
        if (ch %in% c("\"", "'", "`")) {
          quote <- ch
        }
        if (!ch %in% c(" ", "\t")) {
          end <- k
        }
      }
    }
    ends[i] <- end
  }
  ends
}

# The line with a comma after its code, which ends at character `end`, unless
# the code already ends with one.
vendor_add_comma <- function(line, end) {
  code <- substr(line, 1L, end)
  if (endsWith(code, ",")) {
    return(line)
  }
  paste0(code, ",", substring(line, end + 1L))
}

# The checks of the new lines before they are written. They read as DCF. The
# fields that the planned changes (`changes`, and the persons `added` to
# Authors@R) do not edit are unchanged. Each edited dependency field holds the
# old entries plus the additions, with versions only raised. Authors@R holds
# the old persons unchanged, then the new ones. TRUE, or the reason in plain
# text.
vendor_desc_safe <- function(old, new, changes, added) {
  read <- function(d) {
    tmp <- tempfile("DESCRIPTION-")
    on.exit(unlink(tmp))
    vendor_desc_write(d, tmp)
    m <- tryCatch(
      read.dcf(tmp, keep.white = c("Authors@R", "Description")),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(m) && !is.null(d$encoding)) {
      m[] <- iconv(m, from = d$encoding, to = "UTF-8")
    }
    m
  }
  a <- read(old)
  b <- read(new)
  if (is.null(a) || is.null(b) || nrow(b) != 1L) {
    return("the new text does not read as a DESCRIPTION file")
  }
  edited <- unique(unlist(lapply(changes, function(ch) {
    if (ch$kind %in% c("r_add", "r_raise")) "Depends" else ch$field
  })))
  if (!is.null(added)) {
    edited <- c(edited, "Authors@R")
  }
  kept <- setdiff(colnames(a), edited)
  if (!all(kept %in% colnames(b)) || !identical(a[1L, kept], b[1L, kept])) {
    return("a field that should not change would change")
  }
  if (!all(setdiff(colnames(b), colnames(a)) %in% edited)) {
    return("a field would be added that should not be")
  }
  value <- function(m, field) {
    if (field %in% colnames(m)) m[1L, field] else ""
  }
  for (field in intersect(c("Depends", "Imports", "Suggests"), edited)) {
    was <- vendor_parse_deps(value(a, field))
    now <- vendor_parse_deps(value(b, field))
    adds <- unlist(lapply(changes, function(ch) {
      if (identical(ch$kind, "add") && identical(ch$field, field)) {
        vendor_parse_deps(ch$entry)$name
      } else if (identical(ch$kind, "r_add") && field == "Depends") {
        "R"
      }
    }))
    if (!setequal(now$name, c(was$name, adds)) || anyDuplicated(now$name)) {
      return(paste0("the entries of ", field, " would not be as planned"))
    }
    for (k in seq_len(nrow(was))) {
      j <- match(was$name[k], now$name)
      if (identical(was$entry[k], now$entry[j])) {
        next
      }
      raised <- tryCatch(
        identical(now$op[j], ">=") &&
          (is.na(was$version[k]) ||
            numeric_version(now$version[j]) > numeric_version(was$version[k])),
        error = function(e) FALSE
      )
      if (!isTRUE(raised)) {
        return(paste0(
          "an entry of ",
          field,
          " would change other than by a raise"
        ))
      }
    }
  }
  if (!is.null(added)) {
    was <- vendor_read_persons(value(a, "Authors@R"))
    now <- vendor_read_persons(value(b, "Authors@R"))
    n <- length(was)
    if (
      is.null(now) ||
        length(now) != n + length(added) ||
        !identical(format(now[seq_len(n)]), format(was)) ||
        !identical(
          format(now[n + seq_along(added)], include = c("given", "family")),
          format(added, include = c("given", "family"))
        )
    ) {
      return("Authors@R would not hold the old persons and then the new ones")
    }
  }
  TRUE
}

# A planned change as an item left for the user to do
vendor_change_todo <- function(change) {
  field <- change$field
  entry <- change$entry
  pkg <- change$pkg
  new <- change$new
  r_min <- vendor_needs$r
  names <- if (is.null(change$names)) "" else vendor_and(change$names)
  template <- switch(
    change$kind,
    add = "Add {.pkg {entry}} to {.field {field}}.",
    raise = "Raise {.pkg {pkg}} in {.field {field}} to {.code {new}}.",
    r_add = "Add {.code R (>= {r_min})} to {.field Depends}.",
    r_raise = "Raise R in {.field Depends} to {.code (>= {r_min})}.",
    authors = "Add {names} to {.field Authors@R}, with the role {.val cph}."
  )
  vendor_msg("_", template)
}

# The changes as lines of the report of use_vanilla()
vendor_change_msgs <- function(changes) {
  kinds <- vapply(changes, `[[`, "", "kind")
  fields <- vapply(
    changes,
    function(ch) if (is.null(ch$field)) "" else ch$field,
    ""
  )
  msgs <- list()
  r_min <- vendor_needs$r
  for (ch in changes[kinds == "authors"]) {
    names <- vendor_and(ch$names)
    msgs <- c(
      msgs,
      list(vendor_msg(
        "v",
        "Added the engine's copyright holder{?s} to {.field Authors@R} in
       {.file DESCRIPTION}, with the role {.val cph}: {names}."
      ))
    )
  }
  for (field in c("Imports", "Suggests")) {
    deps <- vendor_and(vapply(
      changes[kinds == "add" & fields == field],
      `[[`,
      "",
      "entry"
    ))
    if (length(deps) > 0L) {
      msgs <- c(
        msgs,
        list(vendor_msg(
          "v",
          "Added {.pkg {deps}} to {.field {field}} in {.file DESCRIPTION}."
        ))
      )
    }
    for (ch in changes[kinds == "raise" & fields == field]) {
      pkg <- ch$pkg
      new <- ch$new
      old <- ch$old
      msgs <- c(
        msgs,
        list(
          if (is.na(old)) {
            vendor_msg(
              "v",
              "Raised {.pkg {pkg}} in {.field {field}} to {.code {new}}, from no
           version."
            )
          } else {
            vendor_msg(
              "v",
              "Raised {.pkg {pkg}} in {.field {field}} to {.code {new}}, from
           {.code {old}}."
            )
          }
        )
      )
    }
  }
  for (ch in changes[kinds == "raise" & fields == "Depends"]) {
    pkg <- ch$pkg
    new <- ch$new
    old <- if (is.na(ch$old)) "no version" else ch$old
    msgs <- c(
      msgs,
      list(vendor_msg(
        "v",
        "Raised {.pkg {pkg}} in {.field Depends} to {.code {new}}, from
       {.code {old}}."
      ))
    )
  }
  for (ch in changes[kinds == "r_raise"]) {
    old <- ch$old
    msgs <- c(
      msgs,
      list(
        if (is.na(old)) {
          vendor_msg(
            "v",
            "Raised R in {.field Depends} to {.code (>= {r_min})}, from no
         version."
          )
        } else {
          vendor_msg(
            "v",
            "Raised R in {.field Depends} to {.code (>= {r_min})}, from
         {.code {old}}."
          )
        }
      )
    )
  }
  if (any(kinds == "r_add")) {
    msgs <- c(
      msgs,
      list(vendor_msg(
        "v",
        "Added {.code R (>= {r_min})} to {.field Depends}."
      ))
    )
  }
  msgs
}

## ----- The checks of vanilla_status() ----------------------------------------

# What the DESCRIPTION lacks of what the engine needs, from vendor_desc():
# list(fail, warn), each a list of vendor_msg() details. `plan` is what
# use_vanilla() would change, from vendor_plan_description(), or NULL. Each
# fail item says whether that plan fixes it (`fixed`), and holds the short
# text of the item for the line of use_vanilla() (`short`, a list of a cli and
# a plain text) and the line of a fix by hand (`manual`).
vendor_check_deps <- function(desc, copies, plan = NULL) {
  fail <- list()
  warn <- list()
  changes <- if (is.null(plan)) list() else plan$changes
  add_fail <- function(
    template,
    short,
    manual,
    fixed,
    .envir = parent.frame()
  ) {
    item <- vendor_msg(" ", template, .envir)
    item$fixed <- fixed
    item$short <- list(
      cli = cli::format_inline(vendor_squash(short), .envir = .envir),
      plain = vendor_plain(short, .envir)
    )
    item$manual <- vendor_msg(">", manual, .envir)
    fail[[length(fail) + 1L]] <<- item
  }
  depends <- vendor_parse_deps(desc[["Depends"]])
  imports <- vendor_parse_deps(desc[["Imports"]])
  suggests <- vendor_parse_deps(desc[["Suggests"]])
  r_min <- vendor_needs$r
  k <- match("R", depends$name)
  r_short <- "{.code R (>= {r_min})} in {.field Depends}"
  if (is.na(k)) {
    add_fail(
      r_short,
      r_short,
      "Add {.code R (>= {r_min})} to {.field Depends} by hand.",
      vendor_plan_fixes(changes, "R")
    )
  } else if (!isTRUE(vendor_dep_ok(depends$op[k], depends$version[k], r_min))) {
    has <- depends$entry[k]
    add_fail(
      "{.code R (>= {r_min})} in {.field Depends}: it has {.code {has}}",
      r_short,
      "Change {.code {has}} in {.field Depends} to {.code R (>= {r_min})} by
       hand: {.fn use_vanilla} does not change it.",
      vendor_plan_fixes(changes, "R")
    )
  }
  both <- rbind(depends, imports)
  for (pkg in names(vendor_needs$imports)) {
    min <- vendor_needs$imports[[pkg]]
    entry <- vendor_dep_text(pkg, min)
    short <- "{.pkg {entry}} in {.field Imports}"
    k <- match(pkg, both$name)
    if (is.na(k)) {
      add_fail(
        short,
        short,
        "Add {.pkg {entry}} to {.field Imports} by hand.",
        vendor_plan_fixes(changes, pkg)
      )
    } else if (!isTRUE(vendor_dep_ok(both$op[k], both$version[k], min))) {
      has <- both$entry[k]
      field <- if (k <= nrow(depends)) "Depends" else "Imports"
      add_fail(
        "{.pkg {entry}} in {.field Imports}: it has {.code {has}}",
        short,
        "Change {.code {has}} in {.field {field}} to {.code {entry}} by hand:
         {.fn use_vanilla} does not change it.",
        vendor_plan_fixes(changes, pkg)
      )
    }
    if (pkg %in% suggests$name && pkg %in% both$name) {
      warn[[length(warn) + 1L]] <- vendor_msg(
        " ",
        "{.pkg {pkg}} is also in {.field Suggests}: remove it from there."
      )
    }
  }
  for (pkg in vendor_needs$suggests) {
    if (!pkg %in% c(both$name, suggests$name)) {
      short <- "{.pkg {pkg}} in {.field Suggests}"
      add_fail(
        short,
        short,
        "Add {.pkg {pkg}} to {.field Suggests} by hand.",
        vendor_plan_fixes(changes, pkg)
      )
    }
  }
  lacks <- vendor_collate_lacks(desc, copies)
  if (length(lacks) > 0L) {
    lacks <- vendor_and(lacks)
    add_fail(
      "the engine files in {.field Collate}: it lacks {.file {lacks}}",
      "the engine files in {.field Collate}",
      "Add the engine files to {.field Collate} by hand: {.fn use_vanilla}
       does not change it.",
      FALSE
    )
  }
  list(fail = fail, warn = warn)
}

# Whether the changes that use_vanilla() would make add or raise `name`, R or
# a package.
vendor_plan_fixes <- function(changes, name) {
  any(vapply(
    changes,
    function(ch) {
      switch(
        ch$kind,
        r_add = ,
        r_raise = identical(name, "R"),
        raise = identical(ch$pkg, name),
        add = identical(vendor_parse_deps(ch$entry)$name, name),
        FALSE
      )
    },
    NA
  ))
}

# The engine files that a Collate field lacks: none when there is no Collate
vendor_collate_lacks <- function(desc, files) {
  if (is.na(desc[["Collate"]])) {
    return(character())
  }
  collate <- scan(
    text = desc[["Collate"]],
    what = "",
    quiet = TRUE,
    quote = "'\""
  )
  setdiff(files, file.path("R", collate))
}

# The license against the engine's, GPL (>= 3): "pass" when an alternative
# is GPL version 3 or later only, "warn" when it allows GPL version 3 or later
# only for the package as a whole (GPL with no version, GPL (>= 2), AGPL
# version 3, or GPL version 3 with + file LICENSE), "file" for file LICENSE
# alone, which cannot be checked, and "fail" for anything else. GPL is
# matched by its abbreviation or its name in R's license database, in each
# form that R accepts: GPL-3, GPL version 3, GPL (>= 3) and so on.
vendor_license_grade <- function(license) {
  if (is.null(license) || is.na(license) || !nzchar(trimws(license))) {
    return("fail")
  }
  alts <- trimws(strsplit(license, "|", fixed = TRUE)[[1L]])
  file_re <- "[[:space:]]*\\+[[:space:]]*file[[:space:]]+LICEN[CS]E$"
  grades <- vapply(
    alts,
    function(alt) {
      if (grepl("^file[[:space:]]+LICEN[CS]E$", alt)) {
        return("file")
      }
      file <- grepl(file_re, alt)
      spec <- vendor_license_spec(sub(file_re, "", alt))
      if (is.null(spec)) {
        return("fail")
      }
      v3 <- vendor_license_allows(spec, "3")
      v2 <- vendor_license_allows(spec, "2")
      if (spec$abbrev == "GPL" && v3 && !v2 && !file) {
        "pass"
      } else if (spec$abbrev %in% c("GPL", "AGPL") && v3) {
        "warn"
      } else {
        "fail"
      }
    },
    "",
    USE.NAMES = FALSE
  )
  for (g in c("pass", "warn", "file")) {
    if (g %in% grades) {
      return(g)
    }
  }
  "fail"
}

# A license as list(abbrev, op, version), from the forms that R accepts: the
# abbreviation or name of a license in R's license database, alone, with
# "-<version>" or " version <version>", or with a version range in
# parentheses. NULL for any other text.
vendor_license_spec <- function(text) {
  text <- gsub("[[:space:]]+", " ", trimws(text))
  db <- vendor_license_db()
  abbrev_of <- function(name) {
    k <- match(name, c(db[, "Abbrev"], db[, "Name"]))
    if (is.na(k)) NA_character_ else c(db[, "Abbrev"], db[, "Abbrev"])[k]
  }
  forms <- list(
    c("^(.+) ?\\( ?(<=|>=|<|>) ?([0-9.]+) ?\\)$", "range"),
    c("^(.+) version ([0-9.]+)$", "version"),
    c("^([^ ]+)-([0-9.]+)$", "version")
  )
  for (form in forms) {
    if (!grepl(form[1L], text)) {
      next
    }
    m <- regmatches(text, regexec(form[1L], text))[[1L]]
    abbrev <- abbrev_of(trimws(m[2L]))
    if (is.na(abbrev)) {
      return(NULL)
    }
    if (form[2L] == "range") {
      return(list(abbrev = abbrev, op = m[3L], version = m[4L]))
    }
    return(list(abbrev = abbrev, op = "==", version = m[3L]))
  }
  abbrev <- abbrev_of(text)
  if (is.na(abbrev)) NULL else list(abbrev = abbrev, op = NA, version = NA)
}

# Whether a license spec allows version `v` of its license, such as "3". The
# version is compared with as many parts as the spec gives, so that 3 is
# 3.0.0 in GPL (>= 3.0.0).
vendor_license_allows <- function(spec, v) {
  if (is.na(spec$op)) {
    return(TRUE)
  }
  have <- tryCatch(numeric_version(spec$version), error = function(e) NULL)
  if (is.null(have)) {
    return(FALSE)
  }
  n <- length(unclass(have)[[1L]])
  v <- numeric_version(paste(c(v, rep("0", n - 1L)), collapse = "."))
  switch(
    spec$op,
    "==" = v == have,
    ">=" = v >= have,
    ">" = v > have,
    "<=" = v <= have,
    "<" = v < have,
    FALSE
  )
}

# The Name and Abbrev of R's license database, with the GNU licenses as a
# fallback when the database cannot be read
vendor_license_db <- function() {
  path <- file.path(R.home("share"), "licenses", "license.db")
  db <- tryCatch(
    read.dcf(path, fields = c("Name", "Abbrev")),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(db)) {
    db <- cbind(
      Name = c(
        "GNU General Public License",
        "GNU Affero General Public License",
        "GNU Lesser General Public License"
      ),
      Abbrev = c("GPL", "AGPL", "LGPL")
    )
  }
  db[!is.na(db[, "Abbrev"]) & !is.na(db[, "Name"]), , drop = FALSE]
}

# The license line of vanilla_status() and use_vanilla(), with its fix
vendor_license_msgs <- function(license, pkg = NULL) {
  grade <- vendor_license_grade(license)
  gpl3 <- vendor_gpl3()
  if (is.na(license)) {
    license <- ""
  }
  vignette <- "vignette(\"engine\", package = \"inlavanilla\")"
  of <- if (is.null(pkg)) "" else paste0(" of ", pkg)
  switch(
    grade,
    pass = list(vendor_msg(
      "v",
      "The license{of}, {.val {license}}, allows the engine's license,
       {.val {gpl3}}."
    )),
    warn = list(vendor_msg(
      "!",
      "The license{of}, {.val {license}}, allows the engine only if the
       package as a whole is distributed under GPL version 3 or later: say
       so in {.field License}."
    )),
    file = list(vendor_msg(
      "!",
      "The license{of}, {.val {license}}, cannot be checked: make sure it
       allows code under {.val {gpl3}}."
    )),
    fail = list(
      vendor_msg(
        if (is.null(pkg)) "x" else "!",
        "The license{of}, {.val {license}}, does not allow the engine's
         license, {.val {gpl3}}."
      ),
      vendor_msg(
        ">",
        "A package that holds the engine must be licensed under GPL version
         3 or later. See \"Licensing\" in {.code {vignette}}."
      )
    )
  )
}

# The dependency row of vanilla_status(). Its fixes say which items
# use_vanilla() adds, as `plan` shows, and which must be made by hand.
vendor_deps_row <- function(desc, copies, release, plan = NULL) {
  check <- vendor_check_deps(desc, copies, plan)
  r_min <- vendor_needs$r
  warns <- lapply(check$warn, function(m) {
    m$type <- "!"
    m
  })
  if (length(check$fail) > 0L) {
    fixed <- vapply(check$fail, `[[`, NA, "fixed")
    items <- lapply(check$fail, `[`, c("type", "cli", "plain"))
    command <- "{.code inlavanilla::use_vanilla(release = \"{release}\")}"
    fixes <- list()
    if (all(fixed)) {
      fixes <- list(vendor_msg(
        ">",
        paste(command, "adds them, and copies nothing new.")
      ))
    } else if (any(fixed)) {
      template <- paste(command, "adds {what}, and copies nothing new.")
      shorts <- lapply(check$fail[fixed], `[[`, "short")
      env <- new.env(parent = environment())
      env$what <- vendor_and(vapply(shorts, `[[`, "", "cli"))
      line_cli <- vendor_inline(template, env)
      env$what <- vendor_and(vapply(shorts, `[[`, "", "plain"))
      fixes <- list(list(
        type = ">",
        cli = line_cli,
        plain = vendor_plain(template, env)
      ))
    }
    fixes <- c(fixes, lapply(check$fail[!fixed], `[[`, "manual"))
    return(list(
      status = "fail",
      msgs = c(
        list(vendor_msg(
          "x",
          "{.file DESCRIPTION} lacks some of what the engine needs:"
        )),
        items,
        fixes,
        warns
      )
    ))
  }
  if (length(warns) > 0L) {
    return(list(status = "warn", msgs = warns))
  }
  vendor_row(
    "pass",
    vendor_msg(
      "v",
      "{.file DESCRIPTION} has the engine's dependencies and
       {.code R (>= {r_min})}."
    )
  )
}

# The rights holders of the engine, from the Authors@R of the inlavanilla
# installed here, or NULL when it cannot be read.
vendor_installed_holders <- function() {
  vendor_holders(utils::packageDescription("inlavanilla")[["Authors@R"]])
}

# The authors row of vanilla_status(): the rights holders against Authors@R.
# Its fix names use_vanilla() only when `plan`, what use_vanilla() would
# change, adds them.
vendor_authors_row <- function(
  desc,
  release,
  plan = NULL,
  holders = vendor_installed_holders()
) {
  if (is.null(holders) || length(holders) == 0L) {
    return(vendor_row(
      "skip",
      vendor_msg(
        " ",
        "Skipped: the copyright holders of inlavanilla could not be read."
      )
    ))
  }
  names <- vendor_and(vendor_names_of(holders))
  if (is.na(desc[["Authors@R"]])) {
    return(vendor_row(
      "warn",
      vendor_msg(
        "!",
        "{.file DESCRIPTION} has no {.field Authors@R}: make sure that its
         {.field Author} field names {names} as copyright holders."
      )
    ))
  }
  persons <- vendor_read_persons(desc[["Authors@R"]])
  if (is.null(persons)) {
    return(vendor_row(
      "warn",
      vendor_msg(
        "!",
        "{.field Authors@R} could not be read safely, so it was not checked:
         make sure that it names {names} with the role {.val cph}."
      )
    ))
  }
  state <- vendor_holder_state(holders, persons)
  msgs <- list()
  status <- "pass"
  if (any(state == "absent")) {
    status <- "fail"
    absent <- vendor_and(vendor_names_of(holders[state == "absent"]))
    adds <- !is.null(plan) &&
      any(vapply(plan$changes, `[[`, "", "kind") == "authors")
    fix <- if (adds) {
      vendor_msg(
        ">",
        "{.code inlavanilla::use_vanilla(release = \"{release}\")} adds the
         copyright holders, and copies nothing new."
      )
    } else {
      vendor_msg(
        ">",
        "Add {absent} to {.field Authors@R} by hand, with the role {.val cph}:
         {.fn use_vanilla} cannot edit this {.field Authors@R} safely."
      )
    }
    msgs <- list(
      vendor_msg(
        "x",
        "{.field Authors@R} does not name {absent}, who hold{?s/} rights in the
         engine."
      ),
      fix
    )
  }
  if (any(state == "no_cph")) {
    if (status == "pass") {
      status <- "warn"
    }
    no_cph <- vendor_and(vendor_names_of(holders[state == "no_cph"]))
    msgs <- c(
      msgs,
      list(vendor_msg(
        "!",
        "{.field Authors@R} names {no_cph}, but without the role {.val cph}: add
       it by hand."
      ))
    )
  }
  if (status == "pass") {
    msgs <- list(vendor_msg(
      "v",
      "{.field Authors@R} names the engine's copyright holders: {names}."
    ))
  }
  list(status = status, msgs = msgs)
}
