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
| `job-info.R` | The machine of a chunk: runner image, CPUs, processor, memory, attempt, local resumptions |
| `merge-chunks.R` | The checks before any summary; manifest, completeness, replay list |
| `check-replay.R` | Compares the tasks run again on the merge machine with the chunks' |
| `check-results.R` | Required files; an MCSE above 0 or a note in every performance row; a basis in every decision row; no performance files in a protocol's gate or test |
| `finalize-results.R` | `run-information.csv` (the study's items, then the harness's) and the audit files of the results folder |
| `scan-addresses.R` | Fails on any text shaped like an e-mail address before an upload or before a log is shown |
| `rng-fingerprint.R`, `rng-fingerprint-windows-4.6.1.txt` | Platform fingerprint, and the laptop's output to compare it with |
| `compare-runs.R` | Compares two runs: task files, or summaries and manifests |
| `run-local.R` | Runs every stage on one machine; resumes an interrupted run |
| `harness-lib.R` | Helpers of the harness scripts |
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
the RNG kinds; the repositories; pak; and every installed package with its
version. The runner image's version may differ between jobs; it is recorded
in `timing/jobs.csv` and `run-information.csv`.

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
  on the protocol's list, a `package_sha` on master, and `cost_seconds` from
  the gate's `timing/gate-summary.csv`; a timing gate is `smoke=false`,
  `stage=gate`. Both must be dispatched with `--ref master`, so that the
  harness that checks them is the merged one.
- The plan job stops a study run whose largest chunk is expected to take
  more than 4 hours on 4 workers, and warns above 2 hours; use more chunks.
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
- Results: `gh run download <run id> --repo JUhalt/solomonR -n
  results-<study>-<purpose>-<stage> -D <folder>`. The task files: `gh run
  download <run id> --repo JUhalt/solomonR -p "chunk-*" -D <archive
  folder>`. Artifacts and run logs are deleted after 90 days, so download
  them when the run ends.
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
   does not, such as the number of permutations). A protocol's smoke mode,
   test runs, and gate may not return `performance`, `decisions`, or
   `failures`; the contract stops if they do.
5. Its `run-pins.csv` gives R, the dated snapshot, `libraries,blas,reference`,
   the extra packages, every package version the protocol records, and
   `tool,pak`.
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
  like an e-mail address in them; the same scan runs before every upload.
- R 4.6.1, the dated snapshot (pak included), and Ubuntu's reference BLAS
  and LAPACK are checked in every job before any task.

## Reference

L'Ecuyer, P., Simard, R., Chen, E. J., & Kelton, W. D. (2002). An
object-oriented random-number package with many long streams and substreams.
*Operations Research, 50*(6), 1073-1075.
https://doi.org/10.1287/opre.50.6.1073.358
