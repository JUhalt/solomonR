# Contributing to solomonR

Contributions, bug reports, methodological questions, and feature
proposals are welcome.

## Project philosophy

`solomonR` is designed to make Solomon four-group analyses easier to
conduct, interpret, teach, compare, and reproduce while preserving the
structure and logic of the design.

Contributions should preserve the package’s central principles:

- distinguish historically important Solomon procedures from
  contemporary recommendations and package-specific extensions;
- make the statistical estimand and inferential target explicit;
- preserve the distinction between treatment effects, pretest effects,
  and treatment-by-pretest sensitization;
- treat structurally absent pretests differently from incidental missing
  data;
- disclose assumptions, diagnostics, uncertainty, and unsupported cases;
- avoid presenting different analytical frameworks as interchangeable
  when their estimands or assumptions differ;
- make consequential researcher decisions explicit and reproducible.

## Development workflow

1.  Open or identify an issue describing the proposed change.
2.  Create a feature branch from the current default branch.
3.  Add or update tests alongside substantive code.
4.  Update documentation and examples when user-facing behavior changes.
5.  Run `devtools::document()`, `devtools::test()`, and
    `devtools::check()` before opening a pull request.
6.  For substantive computational changes, include known-answer,
    simulation/recovery, or direct-method regression tests where
    appropriate.
7.  Open a pull request and allow all GitHub Actions checks to complete.

## Statistical-method contributions

A proposed statistical feature should identify:

- the Solomon four-group question it supports;
- the estimand or statistical quantity being reported;
- assumptions and known limitations;
- how uncertainty is represented;
- how structural and incidental missingness are handled;
- how the method relates to historical Solomon procedures and modern
  alternatives;
- why the feature belongs in `solomonR`;
- relevant methodological references or an explicit derivation if the
  method is novel.

## Testing expectations

Tests should prioritize correctness and consequential behavior rather
than coverage percentages for their own sake.

Useful tests include known-answer comparisons, simulation/recovery
studies, boundary cases, malformed inputs, missing-data behavior,
agreement with established implementations, and regression tests for
interpretation and API behavior.

## Pull requests

Pull requests should describe what changed, why it changed, how it was
tested, and any methodological or API decisions that deserve review.

The current development specification is maintained in `ROADMAP.md`.
