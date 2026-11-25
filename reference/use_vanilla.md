# Copy the engine of inlavanilla into a package

Copies the engine of inlavanilla into the package whose root is `path`,
so that the package can fit models with `vanilla_core()` and does not
need inlavanilla to run. The engine is the code behind
[`inla_vanilla()`](https://haziqj.ml/inlavanilla/reference/inla_vanilla.md):
the files `R/vanilla-*.R` of an inlavanilla release. They are copied
byte for byte, under the same names, into the package's `R/` folder, and
listed with their MD5 hashes in the lock `inst/inlavanilla.lock`.

`use_vanilla()` also adds to the package's `DESCRIPTION` what the engine
needs, where it is missing:

- the copyright holders of the engine in `Authors@R`, with the role
  `"cph"`, as the `DESCRIPTION` of inlavanilla at that release lists
  them;

- `cli (>= 3.6.0)`, `stats` and `utils` in `Imports`;

- `future` and `parallelly` in `Suggests`, for the optional parallel
  scans;

- `R (>= 4.1.0)` in `Depends`, where the package asks for an older R or
  for none.

It never removes anything, and changes no other file of the package. It
prints each change, and ends with a list of what is still to do: any
name clashes, the function that calls `vanilla_core()`, the
documentation, the tests, R CMD check and a NEWS bullet.

Run it again with a newer `release` to update the copy.

## Usage

``` r
use_vanilla(release = NULL, files = NULL, clone = NULL, path = ".")
```

## Arguments

- release:

  The release of inlavanilla to copy: its tag on GitHub, such as
  `"v26.11.24"`, or its version, `"26.11.24"`. By default the version of
  inlavanilla installed here, or, when that version is not a release,
  the newest release before it.

- files:

  The engine files to copy, as `"vanilla-vb.R"` or `"R/vanilla-vb.R"`.
  By default every engine file of the release on a first copy, and the
  files in the lock on an update (see Details). The engine files call
  each other, so copy only some of them when you know which you need.

- clone:

  The path of a local clone of inlavanilla, to copy from without the
  internet. It needs git, and the release tag must be in the clone
  (`git fetch --tags`). By default the files come from GitHub.

- path:

  The root of the package, the folder that holds its `DESCRIPTION`, or a
  folder inside it. By default the working directory.

## Value

Invisibly, a data frame with a row per engine file in the lock: `file`,
its `md5` hash, and its `state`, `"new"`, `"restored"`, `"changed"` or
`"unchanged"`. Its attributes `release`, `commit` and `previous` give
the release copied, its commit and the release before it (`NA` for a
first copy), and `description` and `open` the changes made to
`DESCRIPTION` and what is still to do, as plain text.

## Details

**Where the files come from.** By default from the release on GitHub,
<https://github.com/haziqj/inlavanilla>, which needs the internet but
not git. With `clone`, from the release tag in a local clone of
inlavanilla, which needs git but not the internet. The lock records the
release, its commit and the hash of each file, so that
[`vanilla_status()`](https://haziqj.ml/inlavanilla/reference/vanilla_status.md)
can check the copy later.

**Which files.** A first copy takes every engine file of the release, or
the files in `files`. An update copies the files in the lock, and adds
the files in `files`, so that every copy comes from one release. It
reports the engine files of the release that are not copied, and the
ones that the copies need.

**A copy edited by hand.** When an engine file of the package differs
from the lock, `use_vanilla()` stops and changes nothing, so that the
edit is not lost. It also stops for a file `R/vanilla-*.R` that is not
in the lock and differs from the engine file of that name. Never edit an
engine file: the next update would undo the edit. Report the problem at
<https://github.com/haziqj/inlavanilla/issues> instead. The names
`R/vanilla-*.R` are kept for the engine, so give your own files other
names, and leave the engine files out of any code formatter or linter.

**A rights holder already listed.** A copyright holder of the engine who
is already in `Authors@R` is not added again. When that entry lacks the
role `"cph"`, `use_vanilla()` leaves it as it is, and lists the role in
what is still to do.

**License.** The engine is under the GNU General Public License, version
3 or later, so a package that holds it must be too. `use_vanilla()`
warns when the package's `License` does not allow that.

[`vignette("engine", package = "inlavanilla")`](https://haziqj.ml/inlavanilla/articles/engine.md),
also at <https://haziqj.ml/inlavanilla/articles/engine.html>, describes
the interface of `vanilla_core()`, with an example.

## See also

[`vanilla_status()`](https://haziqj.ml/inlavanilla/reference/vanilla_status.md)
to check the copy, and
[`use_vanilla_test()`](https://haziqj.ml/inlavanilla/reference/use_vanilla_test.md)
and
[`use_vanilla_ci()`](https://haziqj.ml/inlavanilla/reference/use_vanilla_ci.md)
to check it in the tests and on GitHub.

## Examples

``` r
if (FALSE) { # \dontrun{
# From the root of your package: copy the engine from the release of the
# inlavanilla installed here
use_vanilla()

# Update the copy to a later release, for example
use_vanilla("v26.12.01")

# Without the internet, from a clone of inlavanilla
use_vanilla(clone = "~/src/inlavanilla")
} # }
```
