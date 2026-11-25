# a first copy takes every engine file of the release

    Code
      out <- use_vanilla(path = root)
    Message
      
      -- Copying the engine of inlavanilla into myhost --
      
      i Release v26.11.24-1, commit c020000, from GitHub.
      v Copied 3 engine files to 'R/'.
      v Wrote the lock 'inst/inlavanilla.lock'.
      v Added the engine's copyright holders to Authors@R in 'DESCRIPTION', with the
        role "cph": Haziq Jamil and Thomas Luu.
      v Added cli (>= 3.6.0), stats and utils to Imports in 'DESCRIPTION'.
      v Added future and parallelly to Suggests in 'DESCRIPTION'.
      v Added `R (>= 4.1.0)` to Depends.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      
      -- Still to do --
      
      [ ] Write the function that calls `vanilla_core()` with your model's
        log-posterior and turns the list it returns into your results.
        `vignette("engine", package = "inlavanilla")` shows how, and so does
        <https://haziqj.ml/inlavanilla/articles/engine.html>.
      [ ] Run `devtools::document()`.
      [ ] Run your tests and R CMD check: `devtools::check()`.
      [ ] Add a NEWS bullet that names the release, such as "Uses the engine of
        inlavanilla v26.11.24-1."
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

# the default release falls back to the newest before the installed

    Code
      out <- use_vanilla(path = root)
    Message
      
      -- Copying the engine of inlavanilla into myhost --
      
      i Release v26.11.24-1, commit c020000, from GitHub.
      i inlavanilla 26.11.25, installed here, is not a release on GitHub, so this
        copies v26.11.24-1, the newest release before it.
      v Copied 3 engine files to 'R/'.
      v Wrote the lock 'inst/inlavanilla.lock'.
      v Added the engine's copyright holders to Authors@R in 'DESCRIPTION', with the
        role "cph": Haziq Jamil and Thomas Luu.
      v Added cli (>= 3.6.0), stats and utils to Imports in 'DESCRIPTION'.
      v Added future and parallelly to Suggests in 'DESCRIPTION'.
      v Added `R (>= 4.1.0)` to Depends.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      
      -- Still to do --
      
      [ ] Write the function that calls `vanilla_core()` with your model's
        log-posterior and turns the list it returns into your results.
        `vignette("engine", package = "inlavanilla")` shows how, and so does
        <https://haziqj.ml/inlavanilla/articles/engine.html>.
      [ ] Run `devtools::document()`.
      [ ] Run your tests and R CMD check: `devtools::check()`.
      [ ] Add a NEWS bullet that names the release, such as "Uses the engine of
        inlavanilla v26.11.24-1."
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

# a subset of files is copied, and the set grows with files

    Code
      out <- use_vanilla(files = "vanilla-a.R", path = root)
    Message
      
      -- Copying the engine of inlavanilla into myhost --
      
      i Release v26.11.24-1, commit c020000, from GitHub.
      v Copied 1 engine file to 'R/'.
      v Wrote the lock 'inst/inlavanilla.lock'.
      v Added the engine's copyright holders to Authors@R in 'DESCRIPTION', with the
        role "cph": Haziq Jamil and Thomas Luu.
      v Added cli (>= 3.6.0), stats and utils to Imports in 'DESCRIPTION'.
      v Added future and parallelly to Suggests in 'DESCRIPTION'.
      v Added `R (>= 4.1.0)` to Depends.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      ! v26.11.24-1 has 2 engine files that are not copied: 'R/vanilla-b.R' and
        'R/vanilla-c.R'.
        'R/vanilla-b.R' binds `b_fun`, which the copies use.
      
      -- Still to do --
      
      [ ] Copy the engine files that the copies need:
        `inlavanilla::use_vanilla(release = "v26.11.24-1", files = "vanilla-b.R")`.
      [ ] Write the function that calls `vanilla_core()` with your model's
        log-posterior and turns the list it returns into your results.
        `vignette("engine", package = "inlavanilla")` shows how, and so does
        <https://haziqj.ml/inlavanilla/articles/engine.html>.
      [ ] Run `devtools::document()`.
      [ ] Run your tests and R CMD check: `devtools::check()`.
      [ ] Add a NEWS bullet that names the release, such as "Uses the engine of
        inlavanilla v26.11.24-1."
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

# the same release again changes nothing

    Code
      out <- use_vanilla(path = root)
    Message
      
      -- The engine of inlavanilla in myhost --
      
      i Release v26.11.24-1, commit c020000, from GitHub.
      v The 3 engine files are those of v26.11.24-1 already.
      v 'DESCRIPTION' already has what the engine needs.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

# an update copies the changed files and rewrites the lock

    Code
      out <- use_vanilla(path = root)
    Message
      
      -- Updating the engine of inlavanilla in myhost --
      
      i Release v26.11.24-1, commit c020000, from GitHub. The lock was at v26.11.24.
      v Copied 3 engine files to 'R/': 1 changed, 2 unchanged.
        Changed: 'R/vanilla-b.R'.
      v Updated the lock 'inst/inlavanilla.lock'.
      v 'DESCRIPTION' already has what the engine needs.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      
      -- Still to do --
      
      [ ] Read what changed since v26.11.24: the NEWS of inlavanilla at
        <https://github.com/haziqj/inlavanilla/blob/v26.11.24-1/NEWS.md>, and the
        engine's code at
        <https://github.com/haziqj/inlavanilla/compare/v26.11.24...v26.11.24-1>.
      [ ] Run your tests and R CMD check: `devtools::check()`.
      [ ] Add a NEWS bullet that names the release, such as "Uses the engine of
        inlavanilla v26.11.24-1."
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

# a release with roxygen tags or without a locked file is refused

    Code
      out <- use_vanilla(release = "v26.11.26", path = root)
    Message
      
      -- Updating the engine of inlavanilla in myhost --
      
      i Release v26.11.26, commit c040000, from GitHub. The lock was at v26.11.24.
      v Copied 2 engine files to 'R/': 1 changed, 1 unchanged.
        Changed: 'R/vanilla-b.R'.
      v Updated the lock 'inst/inlavanilla.lock'.
      v Dropped 'R/vanilla-c.R' from the lock: v26.11.26 no longer has it, and the
        copy was already deleted.
      v 'DESCRIPTION' already has what the engine needs.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      
      -- Still to do --
      
      [ ] Read what changed since v26.11.24: the NEWS of inlavanilla at
        <https://github.com/haziqj/inlavanilla/blob/v26.11.26/NEWS.md>, and the
        engine's code at
        <https://github.com/haziqj/inlavanilla/compare/v26.11.24...v26.11.26>.
      [ ] Run your tests and R CMD check: `devtools::check()`.
      [ ] Add a NEWS bullet that names the release, such as "Uses the engine of
        inlavanilla v26.11.26."
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

# a failed safety check leaves DESCRIPTION and copies the engine

    Code
      out <- use_vanilla(path = root)
    Message
      
      -- Copying the engine of inlavanilla into myhost --
      
      i Release v26.11.24-1, commit c020000, from GitHub.
      v Copied 3 engine files to 'R/'.
      v Wrote the lock 'inst/inlavanilla.lock'.
      ! 'DESCRIPTION' could not be edited safely, so it was left as it was: make the
        changes below by hand.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      
      -- Still to do --
      
      [ ] Add `R (>= 4.1.0)` to Depends.
      [ ] Add cli (>= 3.6.0) to Imports.
      [ ] Add stats to Imports.
      [ ] Add utils to Imports.
      [ ] Add future to Suggests.
      [ ] Add parallelly to Suggests.
      [ ] Add Haziq Jamil and Thomas Luu to Authors@R, with the role "cph".
      [ ] Write the function that calls `vanilla_core()` with your model's
        log-posterior and turns the list it returns into your results.
        `vignette("engine", package = "inlavanilla")` shows how, and so does
        <https://haziqj.ml/inlavanilla/articles/engine.html>.
      [ ] Run `devtools::document()`.
      [ ] Run your tests and R CMD check: `devtools::check()`.
      [ ] Add a NEWS bullet that names the release, such as "Uses the engine of
        inlavanilla v26.11.24-1."
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

# clashes, a license and a Collate field are left to do

    Code
      out <- use_vanilla(path = root)
    Message
      
      -- Copying the engine of inlavanilla into myhost --
      
      i Release v26.11.24-1, commit c020000, from GitHub.
      v Copied 3 engine files to 'R/'.
      v Wrote the lock 'inst/inlavanilla.lock'.
      v Added the engine's copyright holders to Authors@R in 'DESCRIPTION', with the
        role "cph": Haziq Jamil and Thomas Luu.
      v Added cli (>= 3.6.0), stats and utils to Imports in 'DESCRIPTION'.
      v Added future and parallelly to Suggests in 'DESCRIPTION'.
      v Added `R (>= 4.1.0)` to Depends.
      x 2 engine objects have the name of another object, and R keeps only one of
        each, with no warning:
        `b_fun`: 'R/vanilla-b.R' and 'R/zzz.R'
        `abort`: 'R/vanilla-a.R', and imported from rlang in 'NAMESPACE'
      ! The license of myhost, "MIT + file LICENSE", does not allow the engine's
        license, "GPL (>= 3)".
      
      -- Still to do --
      
      [ ] Rename or remove your own objects named above. Never rename an engine
        object. Then run `inlavanilla::vanilla_status()`.
      [ ] 'DESCRIPTION' has a Collate field, which must list every file in 'R/': add
        the engine files to it.
      [ ] License myhost under GPL version 3 or later before you release it. See
        "Licensing" in `vignette("engine", package = "inlavanilla")`.
      [ ] Leave 'R/vanilla-*.R' out of your formatter and linter: the copies must
        stay byte for byte.
      [ ] Write the function that calls `vanilla_core()` with your model's
        log-posterior and turns the list it returns into your results.
        `vignette("engine", package = "inlavanilla")` shows how, and so does
        <https://haziqj.ml/inlavanilla/articles/engine.html>.
      [ ] Run `devtools::document()`.
      [ ] Run your tests and R CMD check: `devtools::check()`.
      [ ] Add a NEWS bullet that names the release, such as "Uses the engine of
        inlavanilla v26.11.24-1."
      [ ] Commit 'R/vanilla-*.R', 'inst/inlavanilla.lock' and 'DESCRIPTION' together.
      
      -- Optional --
      
      * `inlavanilla::use_vanilla_test()` writes a test that fails when an engine
        file differs from the lock.
      * `inlavanilla::use_vanilla_ci()` adds a GitHub Actions workflow that runs
        `vanilla_status()` on every push.

