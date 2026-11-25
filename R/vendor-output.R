# The output of use_vanilla() and vanilla_status(): cli messages, as usethis
# prints them, so that suppressMessages() silences them. They are wrapped to
# out_width(), and the i symbol has the package's colour, vanilla_hex.

# The cli theme: the i symbol in vanilla, a checkbox for the "_" bullets (the
# list of what is still to do), a plain bullet for "*", and level-two
# headings drawn as cat_heading() draws them.
vendor_theme <- function() {
  list(
    ".bullets .bullet-i" = list(
      before = function(x) paste0(col_vanilla(cli::symbol$info), " ")
    ),
    ".bullets .bullet-_" = list(
      "text-exdent" = 2,
      before = function(x) paste0(cli::symbol$checkbox_off, " ")
    ),
    ".bullets .bullet-*" = list(
      before = function(x) paste0(cli::symbol$bullet, " ")
    ),
    h2 = list(
      fmt = function(x) {
        line <- col_vanilla(strrep(cli::symbol$line, 2))
        paste0(line, " ", cli::style_bold(x), " ", line)
      },
      "margin-top" = 1,
      "margin-bottom" = 1
    )
  )
}

# Sets the width of the output and the theme until `.envir` exits
vendor_start <- function(.envir = parent.frame()) {
  old <- options(cli.width = out_width())
  vendor_defer(options(old), .envir)
  cli::cli_div(theme = vendor_theme(), .envir = .envir)
}

# Runs `expr` when the function whose frame is `envir` exits
vendor_defer <- function(expr, envir) {
  thunk <- as.call(list(function() expr))
  do.call(base::on.exit, list(thunk, TRUE, FALSE), envir = envir)
}

# A message: its bullet type ("v", "x", "!", "i", " ", ">", "_" or "*"), and
# its text formatted for cli_bullets() and as plain text.
vendor_msg <- function(type, template, .envir = parent.frame()) {
  list(
    type = type,
    cli = vendor_inline(template, .envir),
    plain = vendor_plain(template, .envir)
  )
}

# Prints a list of messages from vendor_msg()
vendor_say <- function(msgs) {
  if (length(msgs) == 0L) {
    return(invisible())
  }
  text <- vapply(msgs, `[[`, "", "cli")
  cli::cli_bullets(stats::setNames(text, vapply(msgs, `[[`, "", "type")))
}

# A cli template formatted, with its braces escaped, so that cli_bullets()
# prints it as it is. Values reach a template only through variables, so a
# file name with a brace in it cannot break it.
vendor_inline <- function(template, .envir = parent.frame()) {
  text <- cli::format_inline(vendor_squash(template), .envir = .envir)
  gsub("}", "}}", gsub("{", "{{", text, fixed = TRUE), fixed = TRUE)
}

# The plain text of a cli template, for the data frames the functions return
vendor_plain <- function(template, .envir = parent.frame()) {
  old <- options(cli.num_colors = 1L)
  on.exit(options(old))
  text <- cli::format_inline(vendor_squash(template), .envir = .envir)
  text <- cli::ansi_strip(text)
  gsub("\u00a0", " ", text, fixed = TRUE)
}

# A template written over several lines of code, as one line
vendor_squash <- function(template) {
  gsub("[[:space:]]*\n[[:space:]]*", " ", template)
}

# A vector that cli lists as "a, b and c"
vendor_and <- function(x) {
  cli::cli_vec(x, list("vec-last" = " and "))
}

# The engine's license, with non-breaking spaces where the output can show
# them, so that it does not break across lines.
vendor_gpl3 <- function() {
  if (cli::is_utf8_output()) "GPL\u00a0(>=\u00a03)" else "GPL (>= 3)"
}

# A row: its status and its messages, the first the line of the check
vendor_row <- function(status, ...) {
  list(status = status, msgs = list(...))
}
