# a healthy copy passes, with a newer release reported

    Code
      out <- vanilla_status(path = root)
    Message
      
      -- Engine copy in myhost --
      
      i inlavanilla v26.11.24, commit c010000: 3 files, copied on 2026-11-25.
      v The lock 'inst/inlavanilla.lock' is well formed.
      v The 3 engine files match the lock.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      v 'DESCRIPTION' has the engine's dependencies and `R (>= 4.1.0)`.
      v Authors@R names the engine's copyright holders: Haziq Jamil and Thomas Luu.
      v The license, "GPL (>= 3)", allows the engine's license, "GPL (>= 3)".
      v The lock matches release v26.11.24 on GitHub.
      ! A newer release changes copied files: v26.11.24-1 changes 'R/vanilla-b.R'.
        Changes:
        <https://github.com/haziqj/inlavanilla/compare/v26.11.24...v26.11.24-1>
      > Update with `inlavanilla::use_vanilla(release = "v26.11.24-1")`, then run
        your tests.

# offline, the release and update checks are skipped

    Code
      out <- vanilla_status(path = root, error = TRUE)
    Message
      
      -- Engine copy in myhost --
      
      i inlavanilla v26.11.24, commit c010000: 3 files, copied on 2026-11-25.
      v The lock 'inst/inlavanilla.lock' is well formed.
      v The 3 engine files match the lock.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      v 'DESCRIPTION' has the engine's dependencies and `R (>= 4.1.0)`.
      v Authors@R names the engine's copyright holders: Haziq Jamil and Thomas Luu.
      v The license, "GPL (>= 3)", allows the engine's license, "GPL (>= 3)".
      i Skipped 2 checks: GitHub could not be reached. Give `clone`, a local clone of
        inlavanilla, to run them without the internet.
        inlavanilla 26.11.24-1 is installed here, newer than the copy, v26.11.24.

# error = TRUE stops after printing when a check fails

    Code
      vanilla_status(path = root)
    Message
      
      -- Engine copy in myhost --
      
      i inlavanilla v26.11.24-1, commit c020000: 3 files, copied on 2026-11-25.
      v The lock 'inst/inlavanilla.lock' is well formed.
      x 1 engine file does not match the lock:
        'R/vanilla-a.R' differs from the lock.
      > Never edit an engine file. Restore it from git, or delete it and run
        `inlavanilla::use_vanilla(release = "v26.11.24-1")`. Give your own files
        names that do not start with 'vanilla-'.
      v No engine object has the name of another object in 'R/' or of an import in
        'NAMESPACE'.
      v 'DESCRIPTION' has the engine's dependencies and `R (>= 4.1.0)`.
      v Authors@R names the engine's copyright holders: Haziq Jamil and Thomas Luu.
      v The license, "GPL (>= 3)", allows the engine's license, "GPL (>= 3)".
      v The lock matches release v26.11.24-1 on GitHub.
      v No release of inlavanilla is newer than v26.11.24-1.

# without a lock, or with a malformed one, only the lock is checked

    Code
      out <- vanilla_status(path = root)
    Message
      
      -- Engine copy in myhost --
      
      x There is no lock 'inst/inlavanilla.lock', so myhost holds no copy of the
        engine that `use_vanilla()` made.
      > Copy the engine with `inlavanilla::use_vanilla()`.

# inside inlavanilla, nothing is checked

    Code
      out <- vanilla_status(path = root)
    Message
      i This is inlavanilla, which holds the originals of the engine files: nothing
        to check.

