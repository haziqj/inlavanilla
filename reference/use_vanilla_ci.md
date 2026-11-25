# Check the engine copy on GitHub Actions

Writes the GitHub Actions workflow
`.github/workflows/vanilla-status.yaml` in the package whose root is
`path`. On every push to `main` or `master`, and on every pull request,
it installs the latest release of inlavanilla from GitHub and runs
`vanilla_status(error = TRUE)`, which fails the run when a check fails.
A newer release of inlavanilla is reported in the log, and does not fail
the run.

It also adds `^\.github$` to `.Rbuildignore` where that line is missing,
so that R CMD build leaves the folder out.

It is optional:
[`use_vanilla()`](https://haziqj.ml/inlavanilla/reference/use_vanilla.md)
does not add it, and offers it at the end of its list.

## Usage

``` r
use_vanilla_ci(path = ".", overwrite = FALSE)
```

## Arguments

- path:

  The root of the package, or a folder inside it. By default the working
  directory.

- overwrite:

  If `TRUE`, replace a workflow file of that name that differs from the
  one this version of inlavanilla writes. By default such a file is left
  alone, with a message.

## Value

Invisibly, the path of the workflow file.

## See also

[`vanilla_status()`](https://haziqj.ml/inlavanilla/reference/vanilla_status.md)
for the checks.

## Examples

``` r
if (FALSE) { # \dontrun{
use_vanilla_ci()
} # }
```
