## Status

Not submitted. This file records the pre-submission checks (#84) for the
v1.0.0 submission, which follows the maintainer's sign-off (#85).

## Test environments

* Local: Windows 11, R 4.6.1, `R CMD check --as-cran`.
* GitHub Actions: Windows, macOS, and Ubuntu (R release); Ubuntu (R devel,
  also with `--as-cran`; R oldrel-1).
* win-builder: R-release (R 4.6.1) and R-devel.
* mac-builder: macOS 14.4 (Apple M1), R 4.6.1 Patched.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Possibly misspelled words in DESCRIPTION: "Engelenburg" is the surname of
  the author of van Engelenburg (1999), cited in the description.

The note "Version contains large components" in development versions
(0.8.0.9000) does not apply to a release version.
