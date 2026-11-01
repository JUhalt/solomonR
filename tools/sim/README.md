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
| `contract.R` | Sourced by every study script, from the package checkout (contract version 2): stages, stream assignment, task files with their timing, the timing gate's projection, summary order, fault handling |
| `studies.csv` | Registry: study key, issue, script, pins file, whether it is a protocol, commits the pinned commit must include |
| `assign-chunks.R` | Assigns tasks to chunks, longest processing time first (no random numbers); checks the expected hours per chunk when the gate's seconds per unit of cost are given |
| `check-software.R` | R and its build, the pinned package and tool versions, BLAS and LAPACK, the processor flags, RNG kinds, repositories; writes `software.csv` |
| `job-info.R` | The machine of a job (a chunk, the plan, or the merge): runner image, CPUs, processor, memory, attempt, local resumptions |
| `merge-chunks.R` | The checks before any summary; manifest, completeness, replay list |
| `check-replay.R` | Compares the tasks run again on the merge machine with the chunks' |
| `check-results.R` | Required files; an MCSE above 0 or a note in every performance row; a basis in every decision row; in a protocol's gate or test, only the files agreement, run information, timing, and the harness's own (by name; contents are not inspected) |
| `finalize-results.R` | `run-information.csv` (the study's items, then the harness's) and the audit files of the results folder |
| `scan-addresses.R` | Fails on any text shaped like an e-mail address, in a file or in its path, before an upload or before a log is shown; fails closed (a missing path, a format other than .csv, .txt, .json, and .rds, a NUL byte or UTF-16 text, an object it cannot walk) |
| `rng-fingerprint.R`, `rng-fingerprint-windows-4.6.1.txt` | Platform fingerprint, and the laptop's output to compare it with |
| `compare-runs.R` | Compares two runs: task files, or summaries and manifests |
| `run-local.R` | Runs every stage on one machine; resumes an interrupted run |
| `harness-lib.R` | Helpers of the harness scripts |
| `tests/run-tests.R` | The harness's own tests: `check-results.R`, `scan-addresses.R`, `assign-chunks.R`, the workflow's check of the inputs, and `merge-chunks.R` on tampered copies of a run |
| `studies/smoke/` | The harness's smoke study; not a protocol, and it computes no measure for the package's methods |

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
   agreement checks, failures, and timing (start, finish, the job's start,
   and worker seconds for replications and for agreement checks apart).
3. **Merge** (one job). `merge-chunks.R` requires every planned task exactly
   once, in its assigned chunk, with the plan's seed, replications, commit,
   and script; exactly the planned rows of every replication, and the same
   columns in every task; every scenario's replications 1 to R exactly
   once; and the same software in every job. Up to four distinct tasks are
   run again on the merge machine and must be `identical()` to the chunks':
   the most costly task, the most costly task 1 (which carries a scenario's
   agreement checks), and the cheapest tasks of the first and last chunks,
   each from a chunk not yet chosen while one remains. `replay.csv` says
   which of them first ran on another processor model. Only then does the
   study script summarize, from the replications in the order of scenario
   and replication number.

So a result cannot depend on the number of chunks, the job, the worker, or
the order of the tasks, given the same arithmetic on every machine: 1 chunk
and 3 chunks give identical task files (every field except `timing`) and
identical summaries. `compare-runs.R` checks this, and the replay checks the
arithmetic across machines.

"The same software" is every row of `software.csv`: R's version, the r-builds
package and the SHA-256 of R's shared library; `libc6` (glibc),
`libgfortran5`, `libblas3`, and `liblapack3` with their package versions; the
BLAS and LAPACK files R uses, with their SHA-256; `La_version()` and the long
double digits; the processor's AVX2 and FMA flags (both required, so that
glibc picks the same variants of `exp()`, `log()`, and `pow()` everywhere);
the RNG kinds; the repositories of R and of pak; pak; and every installed
package with its version. The runner image's version may differ between
jobs; it is recorded in `timing/jobs.csv` (the chunk jobs),
`timing/plan-job.csv` and `timing/merge-job.csv` (the plan job, which runs
the preflight, and the merge job, which runs the replay and the summary),
and `run-information.csv` (`runner_images`, `plan_runner_image`, and
`merge_runner_image`, with the processor models beside them).

The package commit fixes everything a task computes, because the study
script, its pins, and `contract.R` are read from the package checkout. The
harness checkout (the workflow's own commit, `harness_sha`) only runs the
jobs and checks. On a `pull_request` run, that commit is GitHub's temporary
merge commit of the pull request into master (`refs/pull/<number>/merge`,
on no branch), while `package_sha` is the pull request's head; on a gate
or study run, both are commits on master.

## Running in the cloud

The workflow must be on master before `workflow_dispatch` can start it; the
pull request that changes the harness runs the smoke study through its
`pull_request` trigger. Then, from any machine with `gh`:

```
gh workflow run simulation-run.yaml --repo JUhalt/solomonR --ref master \
  -f study=smoke -f package_sha=<40-character master SHA> \
  -f stage=run -f reps=20 -f chunks=3 -f smoke=true
gh run list --repo JUhalt/solomonR --workflow simulation-run.yaml --limit 3
gh run watch <run id> --repo JUhalt/solomonR
```

The task files record the package commit, so two runs have the same
manifest only if they ran at the same `package_sha`; the summaries do not
depend on it. The harness's own tests after its merge are therefore
defined against each other: T3 (3 chunks, `fail_once=true`, then
`gh run rerun --failed`) and T4 (1 chunk) are dispatched with the same
`package_sha`, a master commit; the content hashes of their manifests must
be identical (`compare-runs.R results <T3 results> <T4 results>` prints
"Manifest content hashes: identical"), and both runs' summaries must
equal those of T1, the first cloud run (run 41255391619, at e13e8bb). The
two `manifest.csv` files themselves differ: their `chunk` and `attempt`
columns follow the division and the reruns, and `file_sha256` hashes the
whole task file, timing included. T1's task files differ from theirs in
the envelope only: the package commit, and the contract's checksum, since
`contract.R` has changed after e13e8bb.

- `smoke=true` (the default) is a code-testing run: the study's test seed,
  at most 20 replications. A study run is `smoke=false`, `stage=run`, `reps`
  on the protocol's list, a `package_sha` on master, and `cost_seconds`,
  the gate's `seconds_per_cost` in `timing/gate-summary.csv` as written
  there (fixed notation; the input also takes an exponent, such as
  `1.2e-05`); a timing gate is `smoke=false`, `stage=gate`. Both must be
  dispatched with `--ref master`, so that the harness that checks them is
  the merged one.
- The plan job multiplies each chunk's expected hours on 4 workers by an
  allowance, then stops a study run whose largest chunk would take more
  than 4 hours and warns above 2 hours; use more chunks. The allowance is
  the pins' `sizing,allowance` if given, otherwise 1.3 for a protocol study
  (the 30% allowance of #116's timing gate) and 1 for any other. 4 hours
  with 1.3 is 3.1 hours expected, against the 5.5 hours a job's tasks get.
- One study run (`smoke=false`) runs at a time. GitHub keeps at most one
  pending run per concurrency group and cancels an older pending one, so a
  third dispatch would cancel the second: start a study only when the
  previous one has finished.
- A failed or timed-out chunk: `gh run rerun <run id> --failed --repo
  JUhalt/solomonR` runs it again from its plan (same tasks, same seeds, same
  results), then the merge. This is the recovery path that test T3 checks,
  including that the rerun reads the first attempt's plan artifact; "Re-run
  all jobs" makes a new plan and recomputes every chunk, and has not been
  tested. GitHub allows at most 50 reruns of a run, within 30 days.
- A fault in the code (an error in a replication that is not a failed fit)
  fails the chunk at once: the chunk's diagnostic artifact holds
  `failed-tasks/<task id>.csv` with the replication and the message, and
  the job's remaining tasks are skipped. It fails the same way on a rerun.
  The remedy is a new commit, a new plan, and an amendment, never a rerun.
- Before the results of a gate or a study run are used, check where the
  run came from: `gh run view <run id> --repo JUhalt/solomonR --json
  event,headBranch,headSha,workflowName` must show
  `"event":"workflow_dispatch"`, `"headBranch":"master"`,
  `"workflowName":"simulation-run"`, and as `headSha` the master commit it
  was dispatched at. A pull request from a
  fork runs this workflow from its own merge commit and can upload
  artifacts with the same names, but it cannot start a `workflow_dispatch`
  run, which needs write access to the repository. `run-information.csv`
  repeats the event, the ref, and the workflow path (`event`, `ref`,
  `workflow_ref`) for the record; the run writes that file itself, so it is
  not the check.
- Results: `gh run download <run id> --repo JUhalt/solomonR -n
  results-<study>-<purpose>-<stage> -D <folder>`. The task files: `gh run
  download <run id> --repo JUhalt/solomonR -p "chunk-*" -D <archive
  folder>`; each artifact `chunk-NNNN` holds the folder `chunk-NNNN/` itself
  (so that the merge job's download gives `chunks/chunk-NNNN/` for 1 chunk as
  for 256), and `gh` puts each artifact in a folder of its name, so the
  archive holds `chunk-NNNN/chunk-NNNN/tasks/`. Artifacts and run logs are
  deleted after 90 days, so download them when the run ends.
- The results folder holds the files that go into the study's folder
  (`performance.csv`, `agreement.csv`, and the others its protocol names),
  plus `run-information.csv`, `manifest.csv`, `software.csv`, `replay.csv`,
  `completeness.csv`, `preflight.csv`, `timing/`, `logs/`, and
  `summary.rds` (the summary as R objects, which `compare-runs.R` reads;
  it does not go into the study's folder). They go into the repository by a
  pull request; the workflow has read permission only. The task files are
  archived outside the repository, where `manifest.csv` gives the SHA-256
  of each.

## Running on the laptop or the PC

From the package root:

```
Rscript tools/sim/run-local.R --study smoke --chunks 3 --reps 20 --out <empty folder>
Rscript tools/sim/run-local.R --out <same folder> --resume
Rscript tools/sim/run-local.R --out <same folder> --rerun-chunk 2
Rscript tools/sim/compare-runs.R results <folder A>/results <folder B>/results
Rscript tools/sim/compare-runs.R tasks <folder A>/chunks <folder B>/chunks
```

`--workers 0` (the default) runs a chunk's tasks in the stage's own process;
`--workers 10` starts 10 PSOCK workers. `--purpose`, `--design`, and
`--reps` are those of the workflow. Each stage writes its log as it runs
(`<folder>/logs/` and each chunk's `log.txt`), so a long chunk can be
followed with `Get-Content -Wait` or `tail -f`.

After a reboot, sleep, or crash, `--resume` keeps the plan and the finished
chunks, runs only the tasks whose files are missing or unreadable, records
the resumption in the chunk's `resumes.csv` (and in `run-information.csv`),
and merges afresh; it refuses to continue if `software.csv` has changed. A
gate or study run needs a clean checkout of a commit.

The harness's own tests need R alone and take a few minutes:

```
Rscript tools/sim/tests/run-tests.R
Rscript tools/sim/tests/run-tests.R --only merge --run <folder of a finished run>
```

Each case builds its input in a temporary folder (`--work` names another)
and runs one check in a new R process: the MCSE, decision, and file-name
rules of `check-results.R`; 47 folders that `scan-addresses.R` must refuse,
built around a placeholder address at example.org that is assembled at run
time, and 5 it must pass; the job sizes of `assign-chunks.R`; the
workflow's "Check the inputs" step, run with bash (`--bash <path>` if bash
is not on the PATH); and, with `--run`, `merge-chunks.R` on tampered copies
of a run's `plan/` and `chunks/`. Any case that ends otherwise than
expected makes the exit status 1. The workflow runs them in the pull
requests that change the harness: the job "Harness tests", and the merge
cases at the end of the merge job, on that run's own chunks. Run them
after any change to a check.

The software check runs without its repository and tool checks, so the
local package versions must match the pins. On Windows, R reports no BLAS or
LAPACK path; the check then requires R's own `Rblas.dll` and `Rlapack.dll`,
refuses them if they carry the name of an optimized library, records their
SHA-256, and compares it with the pins when the pins give
`libraries,rblas_sha256` and `libraries,rlapack_sha256`.

The design of the harness (1.4) had a study script run a whole study in one
process when `SOLOMONR_SIM_STAGE` is unset; the contract does not do that,
and `run-local.R` is the single-machine route instead, so that a local run
passes the same merge checks as a cloud run.

A study runs entirely on one platform. The merge refuses chunks whose
`software.csv` differ, and the protocols' amendments name the platform: a
run on the laptop or the PC instead of the runners needs its own amendment
(platform, workers, BLAS and LAPACK) before it starts.

## Adding a study

In the study's own pull request, before its run (the order rule of #117):

1. Its script under `vignettes/articles/<name>-validation/` sources
   `tools/sim/contract.R`, defines `study` (the fields are listed in
   `sim_check_study()`), and ends with `sim_main(study)`. Its seeds are the
   protocol's: study, gate, and the code-testing seed. Every task's results
   have the same columns and column classes; `cost()` is greater than 0.
2. Its agreement checks run inside `sim_time_agreement()`, so that the gate
   can add their time once instead of multiplying it by R.
3. Its `gate` gives the protocol's constants (`reps`, `task_size`,
   `ladder`, `divisor`, `cap_hours`, `allowance`, `warmup_minutes`), the
   grouping of scenarios the protocol averages over (`group`), and its
   fallback (`fallback`, `fallback_label`); or a `project` function when
   its rule differs from `sim_gate_default()`.
4. Its `summarize()` returns `agreement`, and in the study run
   `performance` (with `mcse` and `mcse_note`), `decisions` (with `rule`,
   `verdict`, and at least one of `mcse`, `threshold`, and `p_value`), and
   `run-information` (item, value: what the protocol lists that the harness
   does not, such as the number of permutations). In a protocol's smoke
   mode, test runs, and gate, it may return only elements named
   `agreement`, `run-information`, and `timing/<name>`; the contract stops
   on any other name, and `check-results.R` on any other file. Both check
   names only, not contents: a performance measure, decision, or failure
   count returned as a column of `agreement`, an item of `run-information`,
   or a `timing/` element would pass them. That these runs compute none
   rests on the study's `summarize()`, whose code its pull request shows.
   A `gate$project` may not
   return a name the harness or the summary writes (other than
   `gate-times`, `gate-projection`, and `gate-summary`).
5. Its `run-pins.csv` gives R, the dated snapshot, `libraries,blas,reference`,
   the extra packages, every package version the protocol records, and
   `tool,pak`; optionally `sizing,allowance`, the factor on a job's
   expected hours (default 1.3 for a protocol). A protocol study (registry
   flag `TRUE`) must also give `r,libR_sha256`, the SHA-256 of `libR.so` of
   the r-builds package that setup-r installs; every job on the runners
   stops if it differs. For `r-4.6.1 1`, the build that test T2 checked
   (run 41255391619), it is
   `70c4de5e29eef6269aec5e1841477c07e84839812aaca3b2f240941bfd223953`.
   Elsewhere (Windows) the library's checksum is recorded, not compared.
6. A row in `studies.csv`, and the key in the workflow's `study` options.
7. Code-testing runs (`smoke=true`) before the gate; the gate's amendment
   and the run's amendment are posted on the issue before each.

## Safety

- `permissions: contents: read`, no secrets, checkouts without stored
  credentials; inputs reach scripts only through the environment.
- No step prints DESCRIPTION, `packageDescription()`, `citation()`,
  `maintainer()`, or git author fields, and package versions are read from
  the Version field only. The study script's output goes to files, which
  the job prints only after `scan-addresses.R` has found no text shaped
  like an e-mail address in them; the same scan runs before every upload,
  and before `check-results.R` names any file of the results. It also
  scans the paths of the files, reads `%40` as `@`, walks every string of
  an `.rds` object (attributes, functions and their environments, promises
  as code, never evaluated), and fails on a path that does not exist and
  on any file it cannot read completely. A chunk job first stops every R
  process left from its tasks step (after a timeout, a PSOCK worker could
  otherwise write a task file after the scan).
- Every job installs its R library afresh from the dated snapshot (pak
  included): no Actions cache, which any workflow on the branch can write
  and to which pak only adds; no Bioconductor repository
  (`PKG_USE_BIOCONDUCTOR=false`); extra packages only from CRAN (a name,
  `any::`, or `cran::`); no pandoc or Quarto. R 4.6.1 (for a protocol, also
  the checksum of its shared library), the snapshot as the only repository
  of R and of pak, pak, the pinned package versions, and Ubuntu's reference
  BLAS and LAPACK are checked in every job before any task. The R package
  itself (the r-builds `.deb`) is installed by setup-r, which checks no
  checksum, before any of these checks runs.

## Reference

L'Ecuyer, P., Simard, R., Chen, E. J., & Kelton, W. D. (2002). An
object-oriented random-number package with many long streams and substreams.
*Operations Research, 50*(6), 1073-1075.
https://doi.org/10.1287/opre.50.6.1073.358
