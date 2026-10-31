# Simulation harness (issue #117)

This folder runs the simulation studies of the v1.0 validation protocols
(#116, and the items of #117) either on GitHub-hosted runners or on one
machine, with the same scripts and the same checks. It splits a study's
tasks among jobs, checks and merges what the jobs return, and computes the
summary only after every check has passed. It never chooses a seed: the
study script does, from its protocol.

`tools/` is in `.Rbuildignore`, and `tools/check-validation-scripts.R` looks
only at `vignettes/articles/*-validation/`, so nothing here is part of the
package or of the validation-scripts check.

## Files

| File | Role |
|---|---|
| `contract.R` | Sourced by every study script, from the package checkout: stages, stream assignment, task files, summary order |
| `studies.csv` | Registry: study key, issue, script, pins file, whether it is a protocol, commits the pinned commit must include |
| `assign-chunks.R` | Assigns tasks to chunks, longest processing time first (no random numbers) |
| `check-software.R` | R, the pinned package versions, BLAS and LAPACK, RNG kinds, repositories; writes `software.csv` |
| `job-info.R` | The machine of a chunk: runner image, CPUs, processor, memory, attempt |
| `merge-chunks.R` | The checks before any summary; manifest, completeness, replay list |
| `check-replay.R` | Compares the tasks run again on the merge machine with the chunks' |
| `check-results.R` | Required files; an MCSE or a note in every performance row; no performance files in a protocol's gate or test |
| `finalize-results.R` | `run-information.csv` and the audit files of the results folder |
| `scan-addresses.R` | Fails on any text shaped like an e-mail address before an upload |
| `rng-fingerprint.R`, `rng-fingerprint-windows-4.6.1.txt` | Platform fingerprint, and the laptop's output to compare it with |
| `compare-runs.R` | Compares two runs: task files, or summaries and manifests |
| `run-local.R` | Runs every stage on one machine |
| `harness-lib.R` | Helpers of the harness scripts |
| `studies/smoke/` | The harness's smoke study, with `fit_solomon_glm()`; not a protocol |

The workflow is `.github/workflows/simulation-run.yaml`, with the composite
action `.github/actions/sim-setup/action.yml`.

## How a run works

1. **Plan** (one job). The study script writes `tasks.csv`: one row per
   task (a scenario and a block of replications), with the 7 integers of
   the task's L'Ecuyer-CMRG substream (L'Ecuyer et al., 2002). Scenario 1
   takes `nextRNGStream()` of the state `set.seed(base)` leaves, scenario
   s the s-th stream, in the order of the scenario table; task k takes the
   stream after k - 1 calls of `nextRNGSubStream()`. Then the protocol's
   stops run (preflight), and `assign-chunks.R` divides the tasks.
2. **Chunks** (up to 256 jobs, at most 16 at a time, 4 workers each). Each
   job runs its tasks; each task sets its substream at its start and draws
   only from it, and is written as `tasks/<task id>.rds` with its envelope
   (the plan row, the package commit, the script's SHA-256), results,
   agreement checks, failures, and timing.
3. **Merge** (one job). `merge-chunks.R` requires every planned task exactly
   once, in its assigned chunk, with the plan's seed, replications, commit,
   and script; every scenario's replications 1 to R exactly once; and the
   same software in every job. Four tasks are run again on the merge machine
   and must be `identical()` to the chunks'. Only then does the study script
   summarize, after sorting every replication by scenario and replication
   number.

So a result cannot depend on the number of chunks, the job, the worker, or
the order of the tasks: 1 chunk and 3 chunks give identical task files
(except their timing) and identical summaries. `compare-runs.R` checks this.

The package commit fixes everything a task computes, because the study
script, its pins, and `contract.R` are read from the package checkout. The
harness checkout (the workflow's own commit) only runs the jobs and checks.

## Running in the cloud

The workflow must be on master before `workflow_dispatch` can start it; the
harness's own branch and pull request run the smoke study through their
`push` and `pull_request` triggers. Then, from any machine with `gh`:

```
gh workflow run simulation-run.yaml --repo JUhalt/solomonR --ref master \
  -f study=smoke -f package_sha=<40-character master SHA> \
  -f stage=run -f reps=20 -f chunks=3 -f smoke=true
gh run list --repo JUhalt/solomonR --workflow simulation-run.yaml --limit 3
gh run watch <run id> --repo JUhalt/solomonR
```

- `smoke=true` (the default) is a code-testing run: the study's test seed,
  at most 20 replications. A study run is `smoke=false`, `stage=run`, `reps`
  on the protocol's list, and a `package_sha` on master; a timing gate is
  `smoke=false`, `stage=gate`.
- Only one study run (`smoke=false`) runs at a time; others queue.
- A failed or timed-out chunk: `gh run rerun <run id> --failed --repo
  JUhalt/solomonR` runs it again from its plan (same tasks, same seeds,
  same files), then the merge.
- Results: `gh run download <run id> --repo JUhalt/solomonR -n
  results-<study>-<purpose> -D <folder>`. The task files: `gh run download
  <run id> --repo JUhalt/solomonR -p "chunk-*" -D <archive folder>`.
  Artifacts and run logs are deleted after 90 days, so download them when
  the run ends.
- The results folder holds the files that go into the study's folder
  (`performance.csv`, `agreement.csv`, and the others its protocol names),
  plus `run-information.csv`, `manifest.csv`, `software.csv`, and
  `timing/`. They go into the repository by a pull request; the workflow has
  read permission only. The task files are archived outside the repository,
  where `manifest.csv` gives the SHA-256 of each.

## Running on the laptop or the PC

From the package root:

```
Rscript tools/sim/run-local.R --study smoke --chunks 3 --reps 20 --out <empty folder>
Rscript tools/sim/run-local.R --out <same folder> --rerun-chunk 2
Rscript tools/sim/compare-runs.R results <folder A>/results <folder B>/results
Rscript tools/sim/compare-runs.R tasks <folder A>/chunks <folder B>/chunks
```

`--workers 0` (the default) runs a chunk's tasks in the stage's own process;
`--workers 10` starts 10 PSOCK workers. `--purpose`, `--design`, and
`--reps` are those of the workflow. The software check runs without its
repository check, so the local versions must match the pins.

A study runs entirely on one platform. The merge refuses chunks whose
`software.csv` differ, and the protocols' amendments name the platform: a
run on the laptop or the PC instead of the runners needs its own amendment
(platform, workers, BLAS and LAPACK) before it starts.

## Adding a study

In the study's own pull request, before its run (the order rule of #117):

1. Its script under `vignettes/articles/<name>-validation/` sources
   `tools/sim/contract.R`, defines `study` (the fields are listed in
   `sim_check_study()`), and ends with `sim_main(study)`. Its seeds are the
   protocol's: study, gate, and the code-testing seed.
2. Its `run-pins.csv` gives R, the dated snapshot, `libraries,blas,reference`,
   the extra packages, and every package version the protocol records.
3. A row in `studies.csv`, and the key in the workflow's `study` options.
4. Code-testing runs (`smoke=true`) before the gate; the gate's amendment
   and the run's amendment are posted on the issue before each.

## Safety

- `permissions: contents: read`, no secrets, checkouts without stored
  credentials; inputs reach scripts only through the environment.
- No step prints DESCRIPTION, `packageDescription()`, `citation()`,
  `maintainer()`, or git author fields, and package versions are read from
  the Version field only. Before every upload, `scan-addresses.R` fails on
  any text shaped like an e-mail address.
- R 4.6.1, the dated snapshot, and Ubuntu's reference BLAS and LAPACK are
  checked in every job before any task.

## Reference

L'Ecuyer, P., Simard, R., Chen, E. J., & Kelton, W. D. (2002). An
object-oriented random-number package with many long streams and substreams.
*Operations Research, 50*(6), 1073-1075.
https://doi.org/10.1287/opre.50.6.1073.358
