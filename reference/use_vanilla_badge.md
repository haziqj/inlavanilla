# Add a README badge of the engine release

Writes the file `.github/inlavanilla-badge.json` in the package whose
root is `path`, and prints the line to add to its README. The badge
reads "vanilla", with the release of inlavanilla that the lock
`inst/inlavanilla.lock` names, such as "v26.11.25-2", in vanilla
colours.

The badge is not a fixed image. The README line asks
[shields.io](https://shields.io/badges/endpoint-badge) to draw it from
the file on GitHub, so it shows the release of the copy on the branch
`branch` as soon as that file is pushed.
[`use_vanilla()`](https://haziqj.ml/inlavanilla/reference/use_vanilla.md)
rewrites the file whenever it updates the copy, so the badge follows the
lock. Commit the file with the lock.

It also adds `^\.github$` to `.Rbuildignore` where that line is missing,
so that R CMD build leaves the folder out.

It is optional:
[`use_vanilla()`](https://haziqj.ml/inlavanilla/reference/use_vanilla.md)
does not add it, and offers it at the end of its list.

## Usage

``` r
use_vanilla_badge(path = ".", branch = "main")
```

## Arguments

- path:

  The root of the package, or a folder inside it. By default the working
  directory.

- branch:

  The branch of the package on GitHub whose file the badge reads. By
  default `"main"`.

## Value

Invisibly, the README line, as Markdown. It takes the GitHub repository
of the package from the `URL` or `BugReports` field of its
`DESCRIPTION`. Where neither names one, the line holds `OWNER/REPO`, to
be replaced by hand.

## See also

[`use_vanilla()`](https://haziqj.ml/inlavanilla/reference/use_vanilla.md),
which keeps the badge in step with the lock.

## Examples

``` r
if (FALSE) { # \dontrun{
# From the root of a package that holds a copy of the engine
use_vanilla_badge()
} # }
```
