# The names check of an engine copy. Without a Collate field, R loads R/ in
# alphabetical order, and an object defined in two files silently takes the
# later definition. An object of the package also silently masks an import
# of the same name. So an engine object that has the name of another object
# of the package breaks the engine or the package's own code, with no warning.

## ----- Names bound by R code -------------------------------------------------

# The names bound at the top level of an R file, and the parse error of a file
# that does not parse: list(names, error). The file is read in `encoding`, the
# encoding of its package, and parsed as UTF-8, since parse() ignores an
# encoding of latin1 in a UTF-8 locale.
vendor_top_names <- function(path, encoding = "UTF-8") {
  exprs <- tryCatch(
    parse(
      text = vendor_read_source(path, encoding),
      keep.source = FALSE,
      srcfile = path,
      encoding = "UTF-8"
    ),
    error = function(e) e
  )
  if (inherits(exprs, "error")) {
    return(list(names = character(), error = conditionMessage(exprs)))
  }
  names <- unique(as.character(unlist(lapply(exprs, vendor_bound_names))))
  list(names = names, error = NA_character_)
}

# The lines of an R file as UTF-8, from the encoding of its package
vendor_read_source <- function(path, encoding) {
  lines <- readLines(path, warn = FALSE)
  if (!toupper(encoding) %in% c("UTF-8", "UTF8")) {
    lines <- iconv(lines, from = encoding, to = "UTF-8")
    if (anyNA(lines)) {
      stop("it cannot be read in its encoding, ", encoding, call. = FALSE)
    }
  }
  lines
}

# The names that a top-level expression binds: assignments with <-, = or <<-
# (chains such as a <- b <- f included), setGeneric("name"), and assign() and
# delayedAssign() with a literal name, also inside a { block or an if or else
# body. A call such as methods::setGeneric counts by its last symbol.
vendor_bound_names <- function(e) {
  if (!is.call(e)) {
    return(character())
  }
  fun <- vendor_call_symbol(e)
  args <- as.list(e)[-1L]
  if (fun %in% c("<-", "=", "<<-")) {
    lhs <- args[[1L]]
    here <- if (is.name(lhs) || (is.character(lhs) && length(lhs) == 1L)) {
      as.character(lhs)
    }
    return(c(here, vendor_bound_names(args[[2L]])))
  }
  if (fun %in% c("{", "(")) {
    return(unlist(lapply(args, vendor_bound_names)))
  }
  if (fun == "if") {
    return(unlist(lapply(args[-1L], vendor_bound_names)))
  }
  if (fun %in% c("setGeneric", "assign", "delayedAssign")) {
    name <- vendor_call_arg(args, if (fun == "setGeneric") "name" else "x")
    if (is.character(name) && length(name) == 1L) {
      return(name)
    }
  }
  character()
}

# The symbol of the function that a call calls: f for f(), and f for pkg::f()
vendor_call_symbol <- function(e) {
  f <- e[[1L]]
  if (
    is.call(f) && is.name(f[[1L]]) && as.character(f[[1L]]) %in% c("::", ":::")
  ) {
    f <- f[[3L]]
  }
  if (is.name(f)) as.character(f) else ""
}

# The argument `name` of a call, given by name or as the first unnamed one
vendor_call_arg <- function(args, name) {
  tags <- names(args)
  if (is.null(tags)) {
    tags <- rep("", length(args))
  }
  if (name %in% tags) {
    return(args[[name]])
  }
  unnamed <- which(tags == "")
  if (length(unnamed) > 0L) args[[unnamed[1L]]] else NULL
}

## ----- NAMESPACE -------------------------------------------------------------

# The imports and exports of a package's NAMESPACE: the names imported one by
# one (importFrom) with the package each comes from, the whole-package imports
# with their exceptions, the exported names and the export patterns. Empty
# parts when there is no NAMESPACE, and `error` when it cannot be read.
vendor_namespace <- function(root) {
  out <- list(
    import_from = data.frame(name = character(), from = character()),
    import_all = list(),
    exports = character(),
    export_patterns = character(),
    error = NA_character_
  )
  if (!file.exists(file.path(root, "NAMESPACE"))) {
    return(out)
  }
  ns <- tryCatch(
    parseNamespaceFile(basename(root), dirname(root)),
    error = function(e) e
  )
  if (inherits(ns, "error")) {
    out$error <- conditionMessage(ns)
    return(out)
  }
  from <- list()
  for (x in ns$imports) {
    if (is.character(x)) {
      out$import_all <- c(
        out$import_all,
        list(list(pkg = x[[1L]], except = NULL))
      )
    } else if (is.list(x) && !is.null(x$except)) {
      out$import_all <- c(
        out$import_all,
        list(list(pkg = x[[1L]], except = as.character(x$except)))
      )
    } else if (is.list(x) && length(x) >= 2L) {
      nm <- x[[2L]]
      local <- if (is.null(names(nm))) {
        nm
      } else {
        ifelse(names(nm) == "", nm, names(nm))
      }
      from <- c(
        from,
        list(data.frame(name = as.character(local), from = x[[1L]]))
      )
    }
  }
  if (length(from) > 0L) {
    out$import_from <- do.call(rbind, from)
  }
  out$exports <- as.character(ns$exports)
  out$export_patterns <- as.character(ns$exportPatterns)
  out
}

## ----- The check -------------------------------------------------------------

# The names check of the engine copies in the package at `root`. Returns a data
# frame with a row per item: `kind`, `level` ("fail", "warn" or "note"),
# `name` (the object, or the package or file the item is about), `files` and
# `names` (list columns), `pkg`, `count`, and `text`, the item as plain text.
# The kinds:
# * duplicate (fail): an engine object is bound in another file of R/;
# * import_from (fail): NAMESPACE imports an engine name with importFrom();
# * import_all (warn): a whole-package import() brings engine names, which
#   the package's own code then does not see;
# * import_unchecked (note): such an import of a package that is not
#   installed here, which cannot be checked;
# * export_pattern (warn): an exportPattern() exports engine objects;
# * export (note): export() of an engine object;
# * parse (warn): a file of R/ that does not parse;
# * namespace (warn): a NAMESPACE that cannot be read;
# * missing (warn, with a snapshot of the release only): names that the
#   copies use and that only an engine file of the release that is not copied
#   binds. all.names() over-counts, so this is not a failure.
vendor_check_names <- function(root, snapshot = NULL) {
  desc <- vendor_desc(root)
  encoding <- if (is.na(desc[["Encoding"]])) "UTF-8" else desc[["Encoding"]]
  dirs <- c("R", file.path("R", c("unix", "windows")))
  files <- unlist(lapply(dirs, function(d) {
    f <- list.files(file.path(root, d), pattern = "[.][RrSsQq]$")
    if (length(f) > 0L) file.path(d, f) else character()
  }))
  parsed <- lapply(
    file.path(root, files),
    vendor_top_names,
    encoding = encoding
  )
  where <- data.frame(
    name = as.character(unlist(lapply(parsed, `[[`, "names"))),
    file = rep(files, vapply(parsed, function(p) length(p$names), 1L))
  )
  engine <- grepl(vendor_file_pattern, where$file)
  engine_names <- sort(unique(where$name[engine]), method = "radix")
  items <- list()
  add <- function(
    kind,
    level,
    name,
    files = character(),
    names = character(),
    pkg = NA_character_,
    count = NA_integer_
  ) {
    items[[length(items) + 1L]] <<- list(
      kind = kind,
      level = level,
      name = name,
      files = files,
      names = names,
      pkg = pkg,
      count = count
    )
  }
  for (name in engine_names) {
    in_files <- where$file[where$name == name]
    if (length(in_files) > 1L) {
      add("duplicate", "fail", name, files = in_files)
    }
  }
  ns <- vendor_namespace(root)
  if (!is.na(ns$error)) {
    add("namespace", "warn", "NAMESPACE", names = ns$error)
  }
  masked <- ns$import_from[
    ns$import_from$name %in% engine_names,
    ,
    drop = FALSE
  ]
  masked <- masked[order(masked$name, method = "radix"), , drop = FALSE]
  for (i in seq_len(nrow(masked))) {
    in_files <- where$file[where$name == masked$name[i] & engine]
    add(
      "import_from",
      "fail",
      masked$name[i],
      files = in_files,
      pkg = masked$from[i]
    )
  }
  for (imp in ns$import_all) {
    if (!requireNamespace(imp$pkg, quietly = TRUE)) {
      add("import_unchecked", "note", imp$pkg, pkg = imp$pkg)
      next
    }
    clash <- setdiff(
      intersect(engine_names, getNamespaceExports(imp$pkg)),
      imp$except
    )
    if (length(clash) > 0L) {
      add(
        "import_all",
        "warn",
        imp$pkg,
        names = sort(clash, method = "radix"),
        pkg = imp$pkg,
        count = length(clash)
      )
    }
  }
  for (pattern in ns$export_patterns) {
    hit <- engine_names[grepl(pattern, engine_names)]
    if (length(hit) > 0L) {
      add(
        "export_pattern",
        "warn",
        pattern,
        names = hit,
        count = length(hit)
      )
    }
  }
  exported <- intersect(engine_names, ns$exports)
  if (length(exported) > 0L) {
    add("export", "note", "export", names = exported, count = length(exported))
  }
  for (i in which(!is.na(vapply(parsed, `[[`, "", "error")))) {
    add(
      "parse",
      "warn",
      files[[i]],
      files = files[[i]],
      names = parsed[[i]]$error
    )
  }
  if (!is.null(snapshot)) {
    for (m in vendor_missing_engine(
      root,
      files[!engine],
      where,
      snapshot,
      encoding
    )) {
      add(
        "missing",
        "warn",
        m$file,
        files = m$file,
        names = m$names,
        count = length(m$names)
      )
    }
  }
  vendor_names_frame(items)
}

# The engine files of the release that the copies need: for each engine file
# of the snapshot that is not copied, the names that a copy uses and that only
# that file binds. Names bound by the package's own files are left out.
vendor_missing_engine <- function(root, own, where, snapshot, encoding) {
  copied <- unique(where$file[grepl(vendor_file_pattern, where$file)])
  uncopied <- setdiff(snapshot$files, copied)
  if (length(uncopied) == 0L || length(copied) == 0L) {
    return(list())
  }
  used <- unique(unlist(lapply(file.path(root, copied), function(f) {
    exprs <- tryCatch(
      parse(f, keep.source = FALSE, encoding = encoding),
      error = function(e) NULL
    )
    all.names(as.expression(exprs))
  })))
  used <- setdiff(used, where$name)
  out <- list()
  for (f in uncopied) {
    bound <- vendor_top_names(file.path(snapshot$dir, f))$names
    need <- sort(intersect(used, bound), method = "radix")
    if (length(need) > 0L) {
      out[[length(out) + 1L]] <- list(file = f, names = need)
    }
  }
  out
}

vendor_names_frame <- function(items) {
  if (length(items) == 0L) {
    out <- data.frame(
      kind = character(),
      level = character(),
      name = character(),
      pkg = character(),
      count = integer(),
      text = character()
    )
    out$files <- list()
    out$names <- list()
    return(out)
  }
  out <- data.frame(
    kind = vapply(items, `[[`, "", "kind"),
    level = vapply(items, `[[`, "", "level"),
    name = vapply(items, `[[`, "", "name"),
    pkg = vapply(items, function(x) as.character(x$pkg), ""),
    count = vapply(items, function(x) as.integer(x$count), 1L)
  )
  out$files <- lapply(items, `[[`, "files")
  out$names <- lapply(items, `[[`, "names")
  out$text <- vapply(
    seq_len(nrow(out)),
    function(i) vendor_names_item(out[i, ], plain = TRUE),
    ""
  )
  out
}

# One item of the names check as the text of a cli bullet (see vendor_inline()),
# or as plain text
vendor_names_item <- function(item, plain = FALSE) {
  name <- item$name
  files <- vendor_and(item$files[[1L]])
  names <- vendor_and(item$names[[1L]])
  pkg <- item$pkg
  count <- item$count
  template <- switch(
    item$kind,
    duplicate = "{.code {name}}: {.file {files}}",
    import_from = "{.code {name}}: {.file {files}}, and imported from
      {.pkg {pkg}} in {.file NAMESPACE}",
    import_all = "{.code import({pkg})} in {.file NAMESPACE} brings
      {count} name{?s} that the engine also defines, so your own code gets
      the engine's object{?s}, not {.pkg {pkg}}'s: {.code {names}}",
    import_unchecked = "Could not check {.code import({pkg})} in
      {.file NAMESPACE}: {.pkg {pkg}} is not installed here.",
    export_pattern = "{.code exportPattern(\"{name}\")} in {.file NAMESPACE}
      exports {count} engine object{?s}: export your own functions by name
      instead.",
    export = "{.file NAMESPACE} exports {count} engine object{?s}:
      {.code {names}}. Document {?it/them} in a file of your own.",
    parse = "{.file {files}} does not parse, so its names were not checked:
      {names}",
    namespace = "{.file NAMESPACE} cannot be read, so its imports were not
      checked: {names}",
    missing = "{.file {name}} binds {.code {names}}, which the copies use."
  )
  if (plain) vendor_plain(template) else vendor_inline(template)
}

# The names check as lines of the report: list(status, msgs). Failures go
# under a cross line with the fix, warnings on lines of their own, and notes
# under the line of the check.
vendor_names_row <- function(items) {
  fails <- items[items$level == "fail", , drop = FALSE]
  warns <- items[items$level == "warn", , drop = FALSE]
  notes <- items[items$level == "note", , drop = FALSE]
  item_msgs <- function(x, type) {
    lapply(seq_len(nrow(x)), function(i) {
      text <- vendor_names_item(x[i, ])
      plain <- x$text[i]
      list(type = type, cli = text, plain = plain)
    })
  }
  msgs <- list()
  status <- "pass"
  if (nrow(fails) > 0L) {
    status <- "fail"
    k <- length(unique(fails$name))
    msgs <- c(
      list(vendor_msg(
        "x",
        "{k} engine object{?s} {?has/have} the name of another object, and R
         keeps only one of each, with no warning:"
      )),
      item_msgs(fails, " "),
      list(vendor_msg(
        ">",
        "Rename or remove your own objects named above. Never rename an
         engine object."
      ))
    )
  }
  if (nrow(warns) > 0L) {
    if (status == "pass") {
      status <- "warn"
    }
    msgs <- c(msgs, item_msgs(warns, "!"))
  }
  if (status == "pass") {
    msgs <- list(vendor_msg(
      "v",
      "No engine object has the name of another object in {.file R/} or of an
       import in {.file NAMESPACE}."
    ))
  }
  msgs <- c(msgs, item_msgs(notes, " "))
  list(status = status, msgs = msgs)
}
