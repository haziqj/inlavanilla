# Make cli optional

For the agent picking this up. Line numbers refer to branch
`claude/dazzling-ritchie-cj1doo` at commit `ac551d7`. Grep again before
editing, because they drift.

## Goal

Move cli from `Imports` to `Suggests`, so that INLAvanilla has no hard
dependencies beyond base R. The README badge goes from 1/1 to 0/0. When cli is
installed, the package uses it for colour, symbols, condition layout and
progress bars, much as now. When it is not, everything still runs, with plain
uncoloured output and no live progress bars.

## The rule that makes this honest

The user considers it dishonest to list a package in `Suggests` while the code
relies on it. `Suggests` is acceptable here only because the package must work
completely without cli. Every exported function, every error, warning and
message, and every `print()`, `summary()` and `plot()` method must run and give
correct, readable text when cli is absent. cli adds polish, nothing more. If
some feature cannot work without cli, stop and ask the user rather than gating
it on cli.

This reverses `lab/DESIGN.md` line 12 ("cli is a required Imports entry") and
the statements in `.claude/CLAUDE.md` that cli is required. Update both (see
Housekeeping).

## Where cli is used now

About 75 `cli::` calls in five files (`grep -rn "cli::" R/`):

- `R/utils.R:109-163`: `abort()`, `warn()`, `format_plain()` and
  `warning_text()` format with `cli::format_error()` and
  `cli::format_warning()`, and strip with `cli::ansi_strip()`.
- `R/progress.R`: `progress_start()`, `progress_tick()` and `progress_done()`
  drive `cli::cli_progress_bar()`. `hl()`, `grey()`, `flag()` and
  `progress_info()` use `col_*()` and `symbol`.
- `R/loo.R`: its own progress bar (about lines 266-364), `format_inline()` and
  `ansi_strip()` to record warnings (about 436 and 489), and `cat_rule()` in
  the print method (about 557).
- `R/methods.R`: `cat_rule()`, `cat_line()`, `style_bold()` and `col_*()`
  throughout `print()` and `summary()`. `cat_status()` (about 445) renders cli
  markup with `format_inline()`. `sd_unit()` (about 491) uses
  `is_utf8_output()`.
- `R/inla-vanilla.R:827`: a `cli::format_message()` info message about
  renamed `init` names.

Message markup: 47 `abort()` and `warn()` call sites, plus the `cat_status()`
messages, are written in cli markup. Find them with
`grep -rnE "\{\.[a-z_]+|\{\?" R/`. They use glue interpolation (`{x}`), the
inline classes `.arg`, `.code`, `.val`, `.field`, `.fn`, `.obj_type_friendly`
and `.or`, and plurals (`{?s}`, `{?has/have}`, `{?y/ies}`, `{?it/them}`,
`cli::qty()`).

Tests:

- 20 `expect_snapshot()` calls (`tests/testthat/_snaps/interface.md`,
  `loo.md`, `marginals-option.md`) capture cli-formatted output.
- `cli::ansi_strip()` in `test-diagnostics.R:41`, `:129` and
  `test-interface.R:1452`.
- `test-interface.R:1177-1207` mocks `cli::cli_progress_bar()` to check that
  the user's `cli.progress_show_after` is honoured and restored.
- `test-interface.R:1394-1399` checks the plain `text` of warnings under
  `cli.condition_width = 30`.

Docs: the roxygen of `inla_vanilla()` (`R/inla-vanilla.R:125-127`) and
`NEWS.md:31` mention `cli.progress_show_after`. Three vignettes set
`options(cli.num_colors = 1, cli.dynamic = FALSE, cli.progress_show_after =
Inf)` in their setup chunk. Setting an option needs no cli, so all of these
stay valid. Only reword the docs to say "when cli is installed".

## Design

### 1. One switch: `use_cli()`

```r
use_cli <- function() {
  requireNamespace("cli", quietly = TRUE)
}
```

Every `cli::` call goes behind this. Do not cache the answer at load time:
tests force the fallback with
`testthat::local_mocked_bindings(use_cli = function() FALSE)`, which only works
if the function is called each time. It costs about half a microsecond once
cli is loaded, but keep it out of the per-call path of the `lp` and `gr`
closures anyway: `progress_new()` calls it once and stores the answer in the
`prog` environment.

### 2. Styling shim (new file `R/style.R`)

Small wrappers that return their input unchanged without cli:

- `col_grey()`, `col_cyan()`, `col_yellow()`, `col_green()`, `col_red()` and
  `style_bold()`: `if (use_cli()) cli::col_grey(...) else paste0(...)`.
- `utf8_output()`: `cli::is_utf8_output()` with cli, otherwise
  `isTRUE(l10n_info()[["UTF-8"]])`. `sd_unit()` uses it.
- Symbols (`sym_tick()`, `sym_info()`, `sym_cross()`): `cli::symbol` with cli,
  otherwise `"✔"`, `"ℹ"` and `"✖"` when `utf8_output()`, and
  cli's own ASCII fallbacks `"v"`, `"i"` and `"x"` when not.
- `strip_ansi()`: always base R, `gsub("\033\\[[0-9;]*m", "", x)`. Use it in
  `warning_text()` and in the tests, in place of `cli::ansi_strip()`.
- `cat_line()`: always base R, `cat(paste0(..., "\n"), sep = "")`, which is
  what cli does.
- `cat_rule(left, right, line_col)`: `cli::cat_rule()` with cli. Without cli,
  draw `-- left ---- right --` to `getOption("width")`, measuring the labels
  with `nchar(strip_ansi(x), type = "width")`.

No colour at all without cli. Deciding when ANSI codes are safe (`NO_COLOR`,
RStudio, Positron, Windows Rgui, Emacs, knitr, `sink()`, CI logs) is the part
of cli that is not worth copying. Users who want colour install cli.

### 3. Conditions and messages

Rewrite every message as plain text built in base R. The text is the same with
or without cli. Do not write an interpreter for cli markup: it would be about
200 lines, with nested braces (`{.or {.val {choices}}}`) and cli's plural
rules, and two renderers of the same markup would drift apart.

Helpers, next to `abort()` in `R/utils.R`. Match cli's plain rendering, which
was checked with `cli.num_colors = 1`:

| cli markup | cli's plain output | Helper |
|---|---|---|
| `{.arg init}`, `{.code -Inf}` | `` `init` ``, `` `-Inf` `` | `fmt_arg()`, `fmt_code()` |
| `{.fn foo}` | `` `foo()` `` | `fmt_fn()` |
| `{.field mean}` | `mean` | `fmt_field()` |
| `{.val {x}}` | `"a", "b", and "c"`, `1, 2, and 3`, `"a" and "b"` | `fmt_val()` |
| `{.or {.val {x}}}` | `"a", "b", or "c"` | `fmt_or()` |
| `{.val {1:30}}` | the first 18 values, then `…, 29, and 30` | `fmt_val()` truncates too |
| `{.obj_type_friendly {x}}` | `a character vector`, `a number`, `NULL` | `type_friendly()` |
| `{?s}`, `{?has/have}` | | `ngettext()` or the existing `plural()` (`R/progress.R:132`) |

`.obj_type_friendly` has three uses (`R/sampling.R:98`, `R/loo.R:179`,
`R/loo.R:300`), so `type_friendly()` only needs to cover the types that can
reach them. When translating a plural, check which count it refers to. cli
takes the nearest interpolated value before it, or the one set by `qty()`.

The helpers may colour their output through the shim, as cli does now. The
plain text of a warning then comes from `strip_ansi()`.

`abort()`, `warn()` and a new `inform()` (for `R/inla-vanilla.R:827`) keep
their message form: a character vector whose names are `""`, `"i"`, `"x"`,
`"!"`, `"*"` or `" "`. Now the elements are finished text, not markup, and the
`.envir` argument goes. Rendering:

1. Squash whitespace in each element with `gsub("\\s+", " ", x)`. Many
   messages are string literals that span several indented source lines (for
   example `R/inla-vanilla.R:408-409`), and cli collapses that whitespace
   today. Without this step, the indentation leaks into the message.
2. With cli, let it add the bullets and the wrapping, without reading braces in
   the text as markup. A user's error message or parameter name may contain
   `{`. Interpolated values are not parsed again, so pass each element as a
   value:

   ```r
   tpl <- stats::setNames(paste0("{msg[[", seq_along(msg), "]]}"), names(msg))
   cli::format_error(tpl, .envir = environment())
   ```

   Use `format_warning()` and `format_message()` in the same way.
3. Without cli, prefix each element by its name (`"i"` gives `sym_info()`,
   `"x"` gives `sym_cross()`, `"!"` gives `"!"`, `"*"` gives `"•"` or
   `"*"`, `" "` gives two spaces, and an unnamed element is left as it is), then
   join with `"\n"`. Do not wrap, as base R does not wrap condition messages.

The `text` field of a warning (recorded in `fit$warnings` and in the
`warnings` of a `loo` result) stays unwrapped and uncoloured. Build it with the
base renderer and `strip_ansi()` whether or not cli is installed, so that it
is identical in both cases. The `parent` chaining in `abort()` stays.

### 4. Progress

- With cli: no change. Keep the `cli.progress_show_after` handling in
  `progress_start()` and the bar in `loo()`.
- Without cli: no live indicator. `progress_start()` returns without creating a
  bar, so `p$id` stays `NULL`, `progress_tick()` is then already a no-op, and
  `progress_done()` prints its stage line as now. Handle `loo()`'s `pb` in the
  same way.
- The default path must not get slower. Time `progress_tick()` and a small fit
  (the README's beta-binomial) before and after.
- Out of scope: `utils::txtProgressBar()` as a fallback. It has no mode for an
  unknown total (the mode and VB stages count calls), no delayed start, and it
  does not clear itself. Ask the user before adding it.

### 5. Print and summary methods

Route `cat_line()`, `cat_rule()`, `style_bold()`, `col_*()` and the symbols
through the shim. `cat_status()` takes finished text, so rewrite its callers'
messages (plurals at about `R/methods.R:155` and `:169`, `{.field}` at about
197-229). Without cli, the output has the same text, no colour, and a rule
drawn with dashes.

## Testing

testthat imports cli, so cli is installed wherever the tests run. The fallback
can only be tested by mocking `use_cli()`:

```r
local_mocked_bindings(use_cli = function() FALSE)
```

- Add `tests/testthat/test-no-cli.R`, which runs the whole surface under the
  mock: `inla_vanilla()` at all three output levels, `print()`, `summary()`,
  `plot()`, `sampling()`, `loo()`, `waic()` and their print methods, and a
  selection of errors and warnings. Check the condition classes
  (`inlavanilla_error`, `inlavanilla_warning` and the specific ones such as
  `inlavanilla_diagnostics_warning`), that no output contains `\033`, and that
  no message contains leftover markup (`{.` or `{?`).
- Add snapshots of the fallback output for `print()`, `summary()` and a few
  conditions, as separate tests.
- Existing snapshots will change wherever message text changes. Review every
  diff with `testthat::snapshot_review()`. The text should differ by no more
  than small details of quoting and punctuation. Do not accept blindly.
- Replace `cli::ansi_strip()` in the tests with the internal `strip_ansi()`.
- Keep the `cli.progress_show_after` test (`test-interface.R:1177`). Add a case
  showing that no bar is created under the mock.
- Healthy fits in every document must stay silent, with and without cli.

An end-to-end check in a library without cli would be ideal, but on this
machine cli sits in R's system library
(`/Library/Frameworks/R.framework/Versions/4.6/Resources/library`), which
`.libPaths()` cannot drop. Rely on the mock locally. A CI step that removes cli
after the check and then runs an example with `Rscript` would cover it, but CI
is kept minimal on purpose, so ask the user first.

## Housekeeping

- `DESCRIPTION`: move `cli (>= 3.6.0)` from `Imports` to `Suggests`.
  `devtools::check()` must give no NOTE about undeclared or unconditional use.
- `README.Rmd:23`: change the badge from `dependencies-1%2F1` to
  `dependencies-0%2F0`, then run `devtools::build_readme()`.
- `NEWS.md`: one bullet. cli is now optional. Without it, output is plain and
  there are no progress bars.
- Roxygen of `inla_vanilla()` (`R/inla-vanilla.R:125-127`): progress bars
  appear when cli is installed. Run `devtools::document()`.
- `.claude/CLAUDE.md`: the Dependencies section ("plus cli, which the user
  requires"), Errors and warnings, Progress, Methods, and a line stating the
  honesty rule above.
- `lab/DESIGN.md` lines 9 and 12: record the new decision and the reason.
- Vignettes: keep the `cli.*` options in their setup chunks.
- Run `air format .`. R source stays ASCII, so write symbols as `\u` escapes.

## Suggested commits

One logical change each, subject line only. Commit locally and never push: the
user pushes.

1. Add `use_cli()` and the styling shim, and route print and progress styling
   through them.
2. Build condition messages in base R. This means the new `abort()`, `warn()`
   and `inform()`, the helpers, and all call sites at once, because the new
   `abort()` no longer reads markup. Update the snapshots in the same commit.
3. Skip progress bars when cli is not installed.
4. Add the tests without cli.
5. Move cli to Suggests and update the badge.
6. Update NEWS, CLAUDE.md and DESIGN.md.

## Done when

- `grep -rn "cli::" R/` finds only the shim, the renderers and the progress
  code, each behind `use_cli()`.
- `grep -rnE "\{\.[a-z_]+|\{\?" R/` finds no cli markup.
- `devtools::test()` passes, including `test-no-cli.R`.
- `devtools::check()` gives 0 errors, 0 warnings and 0 notes.
- A small fit is no slower than before.

## Decisions to confirm if unsure

The user agreed the outline: cli optional, and plainer output without it.
These details are recommendations that have not been confirmed one by one:

- Message text is built in base R even when cli is installed. cli then only
  adds the bullets, the wrapping and the colour.
- No colour at all without cli.
- No progress indicator without cli (no `txtProgressBar()`).
