# Add a test that checks the engine copy

Writes the test file `tests/testthat/test-vanilla-engine.R` in the
package whose root is `path`. The test fails when an engine file differs
from the lock `inst/inlavanilla.lock`, when a file in the lock is
missing, or when a file `R/vanilla-*.R` is not in the lock, so that
`devtools::test()` catches an edit to a copy. It uses base R and
testthat only: the package does not need inlavanilla to run it.

The test reads the source files in `R/`, which an installed package does
not hold. So it runs with `devtools::test()` from the source tree, and
skips under R CMD check.

It is optional:
[`use_vanilla()`](https://haziqj.ml/inlavanilla/reference/use_vanilla.md)
does not add it, and offers it at the end of its list. The package must
already use testthat, as `usethis::use_testthat()` sets up.

## Usage

``` r
use_vanilla_test(path = ".", overwrite = FALSE)
```

## Arguments

- path:

  The root of the package, or a folder inside it. By default the working
  directory.

- overwrite:

  If `TRUE`, replace a test file of that name that differs from the one
  this version of inlavanilla writes. By default such a file is left
  alone, with a message.

## Value

Invisibly, the path of the test file.

## See also

[`vanilla_status()`](https://haziqj.ml/inlavanilla/reference/vanilla_status.md),
which also checks the lock against its release, the names and
`DESCRIPTION`.

## Examples

``` r
if (FALSE) { # \dontrun{
# From the root of a package that holds a copy of the engine
use_vanilla_test()
} # }
```
