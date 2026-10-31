# Records the machine a chunk ran on (issue #117). Usage:
#
#   Rscript tools/sim/job-info.R <chunk output dir> <chunk id>
#
# Writes job-info.csv (item, value): the runner image and its version, the
# number of CPUs, the processor model and whether it has AVX2 and FMA, the
# memory, the workers, the number of tasks, the run attempt, the number of
# local resumptions (tools/sim/run-local.R --resume), and the time.
# Only these named variables are read; nothing else of the environment.

local({
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(file), "harness-lib.R"))
})
args <- h_args(2L, "job-info.R <chunk output dir> <chunk id>")
out_dir <- args[[1]]
env <- function(name) { v <- Sys.getenv(name); if (nzchar(v)) v else NA_character_ }

cpu <- NA_character_
flags <- character()
memory_gb <- NA_real_
if (file.exists("/proc/cpuinfo")) {
  info <- readLines("/proc/cpuinfo", warn = FALSE)
  model <- grep("^model name", info, value = TRUE)
  if (length(model)) cpu <- trimws(sub("^[^:]*:", "", model[[1]]))
  fl <- grep("^flags", info, value = TRUE)
  if (length(fl)) flags <- strsplit(trimws(sub("^[^:]*:", "", fl[[1]])), " +")[[1]]
}
if (file.exists("/proc/meminfo")) {
  mem <- grep("^MemTotal:", readLines("/proc/meminfo", warn = FALSE), value = TRUE)
  if (length(mem)) memory_gb <- round(as.numeric(gsub("[^0-9]", "", mem)) / 1024^2, 1)
}
if (is.na(cpu)) cpu <- env("PROCESSOR_IDENTIFIER")
ids <- file.path(out_dir, "task-ids.txt")

values <- c(
  chunk = args[[2]],
  runner_os = if (is.na(env("RUNNER_OS"))) R.version$os else env("RUNNER_OS"),
  runner_arch = env("RUNNER_ARCH"),
  image_os = env("ImageOS"), image_version = env("ImageVersion"),
  cpus = parallel::detectCores(), cpu_model = cpu,
  avx2 = if (length(flags)) "avx2" %in% flags else NA,
  fma = if (length(flags)) "fma" %in% flags else NA,
  memory_gb = memory_gb,
  workers = if (is.na(env("SOLOMONR_SIM_WORKERS"))) env("SIM_WORKERS") else env("SOLOMONR_SIM_WORKERS"),
  tasks = if (file.exists(ids)) sum(nzchar(readLines(ids, warn = FALSE))) else NA,
  run_id = env("GITHUB_RUN_ID"), run_attempt = env("GITHUB_RUN_ATTEMPT"),
  resumes = if (file.exists(file.path(out_dir, "resumes.csv"))) nrow(h_read_csv(file.path(out_dir, "resumes.csv"))) else 0L,
  recorded_utc = h_now())
h_write_kv(values, file.path(out_dir, "job-info.csv"))
cat("Job: ", values[["cpus"]], " CPUs, ", values[["cpu_model"]], ", image ",
    values[["image_os"]], " ", values[["image_version"]], ".\n", sep = "")
