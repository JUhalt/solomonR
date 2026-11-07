# Simulation harness (issue #117)

This folder runs the simulation studies of the v1.0 validation protocols
(#116, and the items of #117) either on GitHub-hosted runners or on one
machine, with the same scripts and the same checks. It splits a study's
tasks among jobs, checks and merges what the jobs return, and computes the
summary only after every check has passed. It never chooses a task's
seed: the study script does, from its protocol. In a gate or study run,
its only random numbers are those of the RNG fingerprint
(`rng-fingerprint.R`), from `set.seed(116)` in a process of its own that
runs no task. Its own tests draw their own, in its pull requests only.

`tools/` is in `.Rbuildignore`, and `tools/check-validation-scripts.R` looks
only at `vignettes/articles/*-validation/`, so nothing here is part of the
package or of the validation-scripts check.

## Files

| File | Role |
|---|---|
| `contract.R` | Sourced by every study script, from the package checkout (contract version 3): stages, stream assignment, scenario subsets and test task sizes, the SHA-256 of the files a script sources, task files with their timing, the timing gate's projection, summary order, fault handling |
| `studies.csv` | Registry: study key, issue, script, pins file, whether it is a protocol, commits the pinned commit must include |
| `assign-chunks.R` | Assigns tasks to chunks, longest processing time first (no random numbers); checks the expected hours per chunk, with the study's allowance from the plan, when the gate's seconds per unit of cost are given |
| `check-software.R` | R and its build, the pinned package and tool versions, BLAS and LAPACK, the processor flags, RNG kinds, repositories; writes `software.csv` |
| `job-info.R` | The machine of a job (a chunk, the plan, or the merge): runner image, CPUs, processor, memory, attempt, local resumptions |
| `merge-chunks.R` | The checks before any summary; manifest, completeness, replay list |
| `check-replay.R` | Compares the tasks run again on the merge machine with the chunks' |
| `check-results.R` | Required files; an MCSE above 0 or a note in every performance row; a basis in every decision row; in every run, no file under a name the harness writes after it, in any case; in a protocol's gate or test, only agreement, run information, `summary.rds`, and timing (by name; contents are not inspected) |
| `finalize-results.R` | `run-information.csv` (the study's items, then the harness's) and the audit files of the results folder; refuses a file under one of their names, and a `logs/` folder, that is already there, in any case |
| `scan-addresses.R` | Fails on any text shaped like an e-mail address, in a file or in its path, before an upload or before a log is shown; fails closed (a missing path; a format other than .csv, .txt, .json, and .rds; a text file with a NUL byte, a UTF-16 byte-order mark at its start, a gzip header or one of the four zlib headers of a 32 KB window without a preset dictionary (78 01, 78 5e, 78 9c, 78 da) at its start (after a UTF-8 byte-order mark, if it has one), or bytes that are not valid UTF-8; an object it cannot walk); decodes nothing but `%40` (base64 text, a zlib stream with another header, raw deflate, and compressed bytes in an .rds are scanned as stored) |
| `rng-fingerprint.R`, `rng-fingerprint-windows-4.6.1.txt` | Platform fingerprint, and the laptop's output to compare it with |
| `compare-runs.R` | Compares two runs: task files, or summaries and manifests |
| `run-local.R` | Runs every stage on one machine; resumes an interrupted run |
| `harness-lib.R` | Helpers of the harness scripts |
| `tests/run-tests.R`, `tests/synthetic-study.R` | The harness's own tests: `check-results.R`, `scan-addresses.R`, `assign-chunks.R`, the workflow's checks of the inputs and of the package commit, the contract's plan stage on a synthetic study, `finalize-results.R`, `run-local.R` on the synthetic study, and `merge-chunks.R` on tampered copies of a run |
| `studies/smoke/` | The harness's smoke study (its script, and `smoke-summarize.R`, which the script sources); not a protocol, and it computes no measure for the package's methods |

The workflow is `.github/workflows/simulation-run.yaml`, with the composite
action `.github/actions/sim-setup/action.yml`.

## How a run works

1. **Plan** (one job). The study script writes `tasks.csv`: one row per
   task (a scenario and a block of replications), with the 7 integers of
   the task's L'Ecuyer-CMRG substream (L'Ecuyer et al., 2002). Scenario 1
   takes `nextRNGStream()` of the state `set.seed(base)` leaves, scenario
   s the s-th stream, in the order of the scenario table; task k takes the
   stream after k - 1 calls of `nextRNGSubStream()`. A plan of some
   scenarios only (the study's fallback, or a test's list) keeps these
   streams. Then the protocol's stops run (preflight), and
   `assign-chunks.R` divides the tasks.
2. **Chunks** (up to 256 jobs, at most 16 at a time, 4 workers each). Each
   job runs its tasks; each task sets its substream at its start and draws
   only from it, and is written as `tasks/<task id>.rds` with its envelope
   (the plan row, the package commit, the SHA-256 of the script, of
   `contract.R`, and of the files the script sources), results, agreement
   checks, failures, and timing (start, finish, the job's start, and worker
   seconds for replications and for agreement checks apart).
3. **Merge** (one job). `merge-chunks.R` requires every planned task exactly
   once, in its assigned chunk, with the plan's seed, replications, commit,
   script, and sourced files; exactly the planned rows of every
   replication, and the same columns in every task; every scenario's
   replications 1 to R exactly once; and the same software in every job.
   Up to four distinct tasks are run again on the merge machine and must
   be `identical()` to the chunks': the most costly task, the most costly
   task 1 (which carries a scenario's agreement checks), and the cheapest
   tasks of the first and last chunks, each from a chunk not yet chosen
   while one remains (once every chunk has been chosen, the cheapest tasks
   of the lowest- and highest-numbered chunks); ties go by task id, and no
   task is chosen twice. `replay.csv` says which of them first ran on
   another processor model. Only then does the study script summarize,
   from the replications in the order of scenario and replication number.

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
script, the files it sources, its pins, and `contract.R` are read from the
package checkout. The
harness checkout (the workflow's own commit, `harness_sha`) only runs the
jobs and checks. On a `pull_request` run, that commit is GitHub's temporary
merge commit of the pull request into master (`refs/pull/<number>/merge`,
on no branch), while `package_sha` is the pull request's head; on a gate
or study run, both are commits on master.

The contract reads the package commit with git only where the package copy
is the top level of its git repository (`git rev-parse --show-toplevel` is
the copy's folder), with `GIT_DIR` and `GIT_WORK_TREE` unset, as the
validation scripts do (#152, #162), and where `HEAD` is a commit of 40 hex
digits (a repository with no commit, or one with SHA-256 names, has none).
A copy unpacked inside another repository would otherwise record that
repository's commit. Where git is used, a commit given in
`SOLOMONR_SIM_PACKAGE_SHA` must equal `HEAD`. Elsewhere a test run records
`unknown`, or the given commit without a check, and a gate or study run
stops: its plan needs git to find the checkout clean, and its run stage
needs git to check the plan's commit against `HEAD` again. Every job of
the workflow checks out the package commit at the top level of its own
repository.

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
the envelope only: the package commit; the contract's version and
checksum, since `contract.R` has changed after e13e8bb; the script's
checksum, since the smoke study now sources `smoke-summarize.R`; and
`files_sha256`, which contract 3 added (`compare-runs.R tasks` names
every such field).

- `smoke=true` (the default) is a code-testing run: the study's test seed,
  at most 20 replications. A study run is `smoke=false`, `stage=run`, `reps`
  on the protocol's list, a `package_sha` on master, and `cost_seconds`,
  the gate's `seconds_per_cost` in `timing/gate-summary.csv` as written
  there (fixed notation; the input also takes an exponent, such as
  `1.2e-05`); a timing gate is `smoke=false`, `stage=gate`. Both must be
  dispatched with `--ref master`, so that the harness that checks them is
  the merged one.
- The plan job expects each chunk to take its cost times `cost_seconds`
  times the study's allowance, divided by its 4 workers, then stops a study
  run whose largest chunk would take more than 4 hours and warns above 2
  hours; use more chunks. The allowance is that of the study's timing gate
  (`allowance`, which `plan-info.csv` records: 1.3 for #116 and for the
  smoke study; 1 for a study without one), so that a run's jobs are sized
  with the allowance with which its gate chose R. 4 hours with 1.3 is 3.1
  hours expected, against the 5.5 hours a job's tasks get. A study that
  gives `agreement_cost` has task 1 of each scenario priced with its
  agreement checks, and its gate's `seconds_per_cost` then comes from
  replication time alone (`agreement_seconds_per_cost` beside it shows
  how well the agreement checks were priced).
- `scenarios=fallback` plans only the scenarios of the study's fallback
  (its gate's `fallback`), with the full table's streams: a study run after
  a gate that chose the fallback, or a test of that plan; a gate always
  runs every scenario. `task_size` (smoke runs of stage run only) gives a
  test several tasks per scenario, so that a test of at most 20
  replications reaches substreams 2 and later and the merge's checks
  across a scenario's tasks. The plan records both (`subset`,
  `task_size`), and the contract checks them again.
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
- A fault (an error in a replication that is not a failed fit) fails the
  chunk at once: the chunk's diagnostic artifact holds
  `failed-tasks/<task id>.csv` with the replication and the message, and
  the job's remaining tasks are skipped. Download that artifact, then run
  the job again once, on a new machine (`gh run rerun <run id> --failed`;
  every job gets a new virtual machine). A fault that recurs is a code
  fault: it needs a new commit, a new plan, and an amendment. A fault that
  does not recur is reported with the results, and the rerun's task files
  are used.
- The merge job is not run again after one of its checks failed: the
  completeness, envelope, shape, or software checks of `merge-chunks.R`,
  or the replay check of `check-replay.R`. The failure is reported, and an
  amendment says how the study continues (a replayed task that differs is
  not run again to choose between two results). A summary that stops on a
  fault in the study's code stops the same way on a rerun, since it runs
  the same code on the same task files. A fix needs a new commit, and the
  merge refuses task files of another commit.
- A merge job that fails uploads `failed-merge-attempt-<attempt>`, kept
  for 14 days, if its files pass the address scan: those of
  `merged/replay.csv`, `merged/job-info.csv` (the merging machine, recorded
  before the replay), and the logs of the replay and the summary that
  exist. After a replay that differs, `replay.csv` gives each replayed
  task's largest difference in the results, the columns that differ, and
  the processor models of its chunk and of the merging machine. Download
  it when the run ends.
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
`--reps` are those of the workflow. `--subset` plans only some scenarios
with the full table's streams (`fallback`, or in a test a list such as
`"1 18 131"`), and `--task-size` sets a test's replications per task. Each
stage writes its log as it runs (`<folder>/logs/` and each chunk's
`log.txt`), so a long chunk can be followed with `Get-Content -Wait` or
`tail -f`.

After a reboot, sleep, or crash, `--resume` keeps the plan and the finished
chunks, runs only the tasks whose files are missing or unreadable, records
the resumption in the chunk's `resumes.csv` (and in `run-information.csv`),
and merges afresh; it refuses to continue if `software.csv` has changed. A
gate or study run needs a clean checkout of a commit, at the top level of
its git repository.

A new run is of the study `--study` names (smoke without it). `--resume`
and `--rerun-chunk` take the study from the run's `plan/plan-info.csv`
(item `study`), and its script and pins from `studies.csv`, so the
commands above need no `--study` for any study; a `--study` that names
another study is refused before anything runs or changes. Each run first
prints the study it runs. `--rerun-chunk` takes a chunk number from 1 to
256 in at most four digits, with or without leading zeros (`2` and `0002`
are chunk 0002).

The harness's own tests need R alone (and git, for the clean-checkout
cases, and pkgload, for the cases of `run-local.R`) and take a few
minutes:

```
Rscript tools/sim/tests/run-tests.R
Rscript tools/sim/tests/run-tests.R --only merge --run <folder of a finished run>
```

Each case builds its input in a temporary folder (`--work` names another)
and runs one check in a new R process: the MCSE, decision, and file-name
rules of `check-results.R`, including the names the harness writes after
it, in any case; 62 folders that `scan-addresses.R` must refuse, built around a
placeholder address at example.org that is assembled at run time (17 of
them by a named rule: a NUL byte, a gzip or zlib header, or bytes that are
not valid UTF-8, among them compressed text without a NUL byte, a zlib
stream of stored blocks under 78 01, a gzip header after a UTF-8
byte-order mark, and the 7 zlib headers outside the zlib rule that cannot
start UTF-8 text, each refused by the UTF-8 rule), 11 it must pass (among
them UTF-8 text with the characters the package prints, a UTF-8
byte-order mark, empty files, text that starts with `x` but not with a
zlib header, and text that starts with a zlib header the scan does not
refuse), one case of the limit stated under "Safety" (the same stored
blocks under 48 4b, `HK`, pass), and two checks of the text, in the
tests' own process: the "Safety" section and the scanner's header must
say that the other zlib headers are "not in the zlib rule", and they and
the scanner's row under "Files" must say that the gzip and zlib rule
applies after a UTF-8 byte-order mark too; the job sizes and the number
formats of `assign-chunks.R`; the workflow's "Check the inputs" step, run
with bash (`--bash <path>` if bash is not on the PATH), the registry's
required commits, and, with git, the step "Check the package commit" on a
repository made there; the contract's plan
stage, its refusals of numbers it cannot read as whole numbers, and
`sim_verify_plan()` on the synthetic study `tests/synthetic-study.R`, which
loads no package there (with git, also the clean-checkout rule; a copy
inside another repository, whose commit is not read, whose gate and study
plans are refused, and where the run stage refuses a gate plan made at a
top level; and repositories with no commit or with SHA-256 names); the
replay list of `merge-chunks.R` and the items and refusals of
`finalize-results.R` on files made up there; `run-local.R` on the
synthetic study, as a second study beside smoke, in a package folder made
there (the study it takes for `--resume` and `--rerun-chunk`, and its
refusal of another `--study`); and, with `--run`, `merge-chunks.R` on
tampered copies of a run's `plan/` and `chunks/`, and the address scan of
the run's own folders. Any case that ends otherwise than expected makes
the exit status 1. The workflow runs them in the pull requests that change
the harness: the job "Harness tests", and the merge and `run-local.R` cases
at the end of the merge job, on that run's own chunks and with the study's
software. Run them after any change to a check.

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
   have the same columns and column classes; `cost()` is greater than 0,
   and `agreement_cost()`, if given, prices the agreement checks that task
   1 of a scenario carries, in the units of `cost()` (it changes no
   result). `files` lists the files the script sources besides
   `contract.R` (paths from the package root), so that their SHA-256 is in
   the plan and in every envelope, the merge compares it, and a gate or
   study run needs them committed.
2. Its agreement checks run inside `sim_time_agreement()`, so that the gate
   can add their time once instead of multiplying it by R.
3. Its `gate` gives the protocol's constants (`reps`, `task_size`,
   `ladder`, `divisor`, `cap_hours`, `allowance`, `warmup_minutes`), the
   grouping of scenarios the protocol averages over (`group`), and its
   fallback (`fallback`, `fallback_label`); or a `project` function when
   its rule differs from `sim_gate_default()`. Its `allowance` also sizes
   the run's jobs (through `plan-info.csv`).
4. Its `summarize()` returns `agreement`, and in the study run
   `performance` (with `mcse` and `mcse_note`), `decisions` (with `rule`,
   `verdict`, and at least one of `mcse`, `threshold`, and `p_value`), and
   `run-information` (item, value: what the protocol lists that the harness
   does not, such as the number of permutations). In a protocol's smoke
   mode, test runs, and gate, it may return only elements named
   `agreement`, `run-information`, and `timing/<name>`; the contract stops
   on any other name. In test and gate runs, `check-results.R` then allows
   in the results folder only `agreement.csv`, `run-information.csv`,
   `summary.rds`, and files directly under `timing/`. In every run it also
   refuses a file under a name that the harness writes after it
   (`manifest.csv`, `completeness.csv`, `software.csv`, `replay.csv`,
   `preflight.csv`, `timing/jobs.csv`, `timing/plan-job.csv`,
   `timing/merge-job.csv`, `timing/rng-fingerprint.txt`,
   `timing/sessionInfo-<chunk>.txt`, and anything under `logs/`), compared
   without regard to case, and anything under `timing/` spelled in another
   case (`Timing/`): a file system that ignores case, as Windows' and
   macOS' do, makes `Preflight.csv` the harness's `preflight.csv`, and on
   one that does not, the results would hold both. Only a `summarize()`
   that writes into `SOLOMONR_SIM_OUT` itself can have put such a file
   there, since both routes run the check after the summary and before
   `finalize-results.R`. That script refuses the same names, in the same
   way, before it writes anything, even a name it has nothing to put in
   (`preflight.csv` of a study without checks). These checks
   read names only, not contents: a performance measure, decision, or
   failure count returned as a column of `agreement`, an item of
   `run-information`, or a `timing/` element, or written as a file under
   `timing/` with a name the harness does not write, would pass them. They
   also look at the results folder only. `finalize-results.R` copies its
   audit files from the folders the summary read (`manifest.csv`,
   `completeness.csv`, `software.csv`, `replay.csv`, `jobs.csv`, and the
   merge job's `job-info.csv`, written before the replay, from `merged/`;
   `preflight.csv`, `rng-fingerprint.txt`, and the plan job's
   `job-info.csv` from `plan/`; `session-info.txt` from each chunk), and
   it reads most items of `run-information.csv` from them and from
   `plan/plan-info.csv` (`replay` and `reps` among them), as the summary
   left them. The summary runs on the same machine as
   those folders, the harness's scripts, and every check after it, and
   can write to all of them (`SOLOMONR_SIM_PLAN` names the plan folder,
   and `SOLOMONR_SIM_IN` a file in `merged/`), so no check after it can
   show that it wrote nowhere else. That these runs compute none, and that
   `summarize()` writes nowhere but `SOLOMONR_SIM_OUT`, rest on the
   study's code, which its pull request shows. A `gate$project` may not
   return a name the harness or the summary writes (other than
   `gate-times`, `gate-projection`, and `gate-summary`).
5. Its `run-pins.csv` gives R, the dated snapshot, `libraries,blas,reference`,
   the extra packages, every package version the protocol records, and
   `tool,pak`. Every pinned package must be installed, as a dependency of
   the package commit or an extra package; otherwise the software check
   stops the run in its Plan job (and `run-local.R` stops too). A package
   that leaves `Imports` leaves the pins too. A protocol study (registry flag `TRUE`) must also give
   `r,libR_sha256`, the SHA-256 of `libR.so` of the r-builds package that
   setup-r installs; every job on the runners stops if it differs. For
   `r-4.6.1 1`, the build that test T2 checked (run 41255391619), it is
   `70c4de5e29eef6269aec5e1841477c07e84839812aaca3b2f240941bfd223953`.
   Elsewhere (Windows) the library's checksum is recorded, not compared.
6. A row in `studies.csv`, and the key in the workflow's `study` options.
   The row's `required_ancestors` are the commits that the package commit
   of a gate or study run must include, as full SHAs separated by spaces;
   the Plan job checks them, and that the commit is on master, before
   anything is installed (a smoke run checks neither). The smoke study's
   row requires #132's merge (b4ade64), which master holds, so that a run
   of the smoke study with `smoke=false` takes that check too.
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
  as code, never evaluated), and fails on a path that does not exist, on
  any file it cannot read completely, and on a text file (.csv, .txt,
  .json) that holds a NUL byte, starts with a UTF-16 byte-order mark,
  starts (after a UTF-8 byte-order mark, if it has one) with a gzip header
  or with one of the four zlib headers of a 32 KB window without a preset
  dictionary (78 01, 78 5e, 78 9c, 78 da; zlib's default window), or is
  not valid UTF-8. Valid UTF-8 text that starts (after a UTF-8 byte-order
  mark, if it has one) with `x^`, with `x` and U+0001, or with `x` and a
  character from U+0680 to U+06BF is refused as a zlib header too. These
  rules catch, in a text file, a gzip stream or a zlib stream with one of
  those four headers at its start or after a UTF-8 byte-order mark, and
  bytes that are not UTF-8. Beyond `%40`, the scan decodes
  nothing: a zlib stream with another header (a smaller window, or a
  preset dictionary), raw deflate, and text compressed or encoded in
  another way (base64, for one) are not detected when their bytes stay
  valid UTF-8 and hold no NUL byte, as a stream of stored deflate blocks
  can, with an address split between two blocks. The other zlib headers
  are not in the zlib rule because many of them can start ordinary text:
  of the 28 others without a preset dictionary, 21 can start valid UTF-8
  text, 14 of them as two ASCII characters (such as `HK`, `(S`, or `X` and
  a tab). The other 7 cannot (their second byte is a UTF-8 continuation
  byte), so the UTF-8 rule refuses a file that starts with one of them.
  Raw deflate has no header at all.
  In an `.rds` the strings and raw vectors are scanned as they are
  stored, so compressed bytes there, in a raw vector or in a string, are
  not detected either. A chunk job first stops every R process left from
  its tasks step (after a timeout, a PSOCK worker could otherwise write a
  task file after the scan).
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
