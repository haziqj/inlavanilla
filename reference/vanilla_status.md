# Check a copy of the engine of inlavanilla

Checks the copy of the engine that
[`use_vanilla()`](https://haziqj.ml/inlavanilla/reference/use_vanilla.md)
put into the package whose root is `path`, and prints one line per
check: a tick when it passes, a cross when it fails, and an exclamation
mark for something to look at, each with what to do. The checks, in this
order:

- **Lock**: `inst/inlavanilla.lock` exists and can be read.

- **Copies**: each engine file has the MD5 hash in the lock, and every
  file `R/vanilla-*.R` is in the lock. The hash is taken after carriage
  returns are removed, so that a checkout with Windows line endings
  passes.

- **Names**: no object of an engine file is also defined in another file
  of `R/`, or imported by name in `NAMESPACE`. R keeps one of the two
  with no warning, so a clash breaks the engine or your code silently.

- **Dependencies**: `DESCRIPTION` has `cli (>= 3.6.0)`, `stats` and
  `utils` in `Imports` (or `Depends`), `future` and `parallelly` in
  `Suggests` (or above), and asks for R 4.1.0 or later.

- **Authors**: `Authors@R` names the copyright holders of the engine,
  with the role `"cph"`.

- **License**: the package's `License` allows code under the GNU General
  Public License, version 3 or later.

- **Release**: the release in the lock is a tag on GitHub at the commit
  in the lock, and each engine file at that tag has the hash in the
  lock.

- **Update**: no newer release of inlavanilla changes a copied file. A
  newer release that does is reported, with the command that updates the
  copy, and is not a failure.

The last two checks need GitHub, or a local clone of inlavanilla given
as `clone`. Without either they are skipped with a note, and the others
still run.

## Usage

``` r
vanilla_status(path = ".", clone = NULL, error = FALSE)
```

## Arguments

- path:

  The root of the package, or a folder inside it. By default the working
  directory.

- clone:

  The path of a local clone of inlavanilla, for the last two checks
  without the internet. It needs git, and its tags must be up to date
  (`git fetch --tags`). By default they use GitHub.

- error:

  If `TRUE`, stop with an error when a check fails, as continuous
  integration needs. A skipped check or an exclamation mark is never an
  error.

## Value

Invisibly, a data frame with a row per check, in the order above:
`check` (`"lock"`, `"copies"`, `"names"`, `"dependencies"`, `"authors"`,
`"license"`, `"release"` or `"update"`), `status` (`"pass"`, `"warn"`,
`"fail"` or `"skip"`), `message`, the line printed for it, and
`details`, the lines printed under it, as plain text.

## See also

[`use_vanilla()`](https://haziqj.ml/inlavanilla/reference/use_vanilla.md)
to copy and update the engine, and
[`use_vanilla_ci()`](https://haziqj.ml/inlavanilla/reference/use_vanilla_ci.md)
to run these checks on GitHub.

## Examples

``` r
if (FALSE) { # \dontrun{
# From the root of a package that holds a copy of the engine
vanilla_status()

# In continuous integration: an error when a check fails
vanilla_status(error = TRUE)
} # }
```
