# Array task entry point for the fixed-K clustering blocks: pilot, core grid and
# variations.
#
#   Rscript R/run_cell.R <block> <task_id> [n_tasks] [outdir] [n_reps] [n_override]
#
#   block   pilot | core | variations
#
# The K-recovery, estimation and settings blocks have different output shapes
# and get their own runners; overloading this one would make every caller pay
# for the union of their arguments.
#
# WHY THE JOB LIST IS BUILT THE SAME WAY IN EVERY TASK
#
# Each task builds the FULL scenario x replicate list and then takes the slice
# whose index is congruent to its own task id. No task needs to know what any
# other task did, and the per-cell seed comes from the cell's own identity, so:
#
#   - re-running one failed task reproduces exactly the same data;
#   - the array can be resized without invalidating anything already computed;
#   - results are independent of the order tasks happen to execute in.
#
# That last property is what makes the array restartable. Do not replace the
# modulo assignment with a contiguous block split unless you also stop deriving
# seeds from cell identity.

source(file.path("R", "methods.R"))

args     <- commandArgs(trailingOnly = TRUE)
block    <- if (length(args) >= 1) args[1] else "pilot"
task_id  <- if (length(args) >= 2) as.integer(args[2]) else 1L
n_tasks  <- if (length(args) >= 3) as.integer(args[3]) else 1L
outdir   <- if (length(args) >= 4) args[4] else file.path("results", block)
n_reps   <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_
n_over   <- if (length(args) >= 6) as.integer(args[6]) else NA_integer_

cfg <- load_config()
K   <- gint(cfg, "K_true")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# The job list lives in R/jobs.R, shared with combine.R so that the
# completeness check can compute what SHOULD be present from the config rather
# than from the rows it happens to see.
source(file.path("R", "jobs.R"))

jobs <- build_jobs(cfg, block, n_reps)
if (!is.na(n_over)) {
  # SMOKE RUNS ONLY. This overwrites every row's n, including the per-scenario
  # values the sample-size variation rows carry, so a run using it cannot be
  # reported. Loud on purpose.
  cat(sprintf("*** n OVERRIDDEN to %d for every cell. SMOKE RUN -- not reportable. ***\n",
              n_over))
  jobs$n <- n_over
}
mine <- jobs[(seq_len(nrow(jobs)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]

inc     <- unique(jobs$include)
methods <- cfg$methods[as.logical(cfg$methods[[inc]]), , drop = FALSE]

cat(sprintf("block %s | task %d of %d | %d of %d cells | %d methods: %s\n",
            block, task_id, n_tasks, nrow(mine), nrow(jobs), nrow(methods),
            paste(methods$method_id, collapse = " ")))

# An empty slice is almost always a mistake, and a silent one: it writes a valid
# empty file and exits 0. The usual cause is re-running a failed task with
# n_tasks = 1, so the divisor no longer
# matches the array the cell was assigned under. Re-run with the ORIGINAL count:
#   Rscript R/run_cell.R core 57 100 <outdir>
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO cells (the block has %d).",
                     "If this is a re-run of one failed task, pass the ORIGINAL",
                     "number of tasks as n_tasks; with n_tasks = 1 every",
                     "cell is reassigned."),
               task_id, n_tasks, nrow(jobs)), call. = FALSE)

# --- run ----------------------------------------------------------------------

rows <- vector("list", nrow(mine) * nrow(methods))
k <- 0L
t_start <- Sys.time()

for (i in seq_len(nrow(mine))) {
  jb  <- mine[i, ]
  arm <- cfg$arms[cfg$arms$arm == jb$arm, ]
  seed <- cell_seed(gint(cfg, "seed_master"), jb$scenario_id, jb$rep)

  dat <- tryCatch(
    generate_replicate(
      cfg, arm, d = jb$d, rep = jb$rep, scenario_id = jb$scenario_id,
      contrast = jb$contrast, M_true = jb$M_true,
      n = if (is.na(jb$n)) NULL else jb$n,
      pi = if (is.na(jb$pi)) NULL else as.numeric(strsplit(jb$pi, ",")[[1]]),
      lambda = if (is.na(jb$lambda)) NULL else jb$lambda,
      error_sd = if (is.na(jb$error_sd)) NULL else jb$error_sd,
      window = jb$window, count_mode = jb$count_mode,
      obs_min_override = jb$obs_min_override),
    error = function(e) structure(conditionMessage(e), class = "genfail"))

  if (inherits(dat, "genfail")) {
    # A generation failure is fatal for the CELL, never for the task, and it is
    # recorded so the combined table shows a hole rather than silently omitting
    # the cell and shrinking the denominator.
    #
    # ONE ROW PER METHOD, with exactly the columns a successful row has. Two
    # reasons, both of which bit an earlier version. do.call(rbind, ...) on data
    # frames whose column sets differ is a hard error, so a short failure row
    # killed the whole task at the final rbind -- after every fit in the slice
    # had already been computed. And emitting one row instead of `nrow(methods)`
    # would make the completeness check report the task as short, conflating a
    # recorded generation failure with a task that died.
    for (m in seq_len(nrow(methods))) {
      k <- k + 1L
      rows[[k]] <- cbind(
        mine[i, ],
        empty_scores(methods$method_id[m], methods$label[m],
                     paste("GENERATION:", as.character(dat)))[
                       , setdiff(names(empty_scores()), "label")],
        n_obs_mean = NA_real_, n_obs_min = NA_integer_, n_obs_max = NA_integer_,
        row.names = NULL)
    }
    next
  }

  for (m in seq_len(nrow(methods))) {
    res <- score_method(methods[m, ], dat, cfg, K, seed = seed)
    k <- k + 1L
    rows[[k]] <- cbind(mine[i, ], res[, setdiff(names(res), "label")],
                       n_obs_mean = mean(dat$n_obs),
                       n_obs_min = min(dat$n_obs), n_obs_max = max(dat$n_obs),
                       row.names = NULL)
  }
  if (i %% 5L == 0L || i == nrow(mine))
    cat(sprintf("  %d/%d cells, %.1f min elapsed\n", i, nrow(mine),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

out <- do.call(rbind, rows[seq_len(k)])

# PROVENANCE. Two task files in one directory can differ in their slicing or in
# the code that produced them, at exactly 100 percent completeness, and nothing
# else here can see it. The combiner refuses to merge files that disagree.
out$n_tasks   <- n_tasks
out$code_hash <- source_hash(file.path("R", "run_cell.R"))
f <- file.path(outdir, sprintf("%s_task%04d.rds", block, task_id))
saveRDS(out, f)
cat(sprintf("wrote %s  (%d rows, %.1f min)\n", f, nrow(out),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
