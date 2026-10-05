# STAGE F0 (36 months) -- WHICH KENDALL-SURFACE BANDWIDTH FITS THE 36-MONTH WINDOW?
#
#   Rscript R/run_realdata36_f0.R <task_id> [n_tasks] [outdir] [n_sub]
#
#   n_sub  fit a random subsample of this many subjects -- SMOKE ONLY; the task
#          file is stamped and the reader says so.
#
# The 18-month decision rested on F0/F0B/F0C, all on the RAW scale. This sweep
# asks the question again on BOTH scales and with BOTH starts the analysis uses,
# one fit per cell:
#
#   scale {log, raw} x start {random, two-staged} x bwK {6, 8, 10, 12} x seed {1, 2, 3}
#   = 48 cells
#
# "random" is the start of F1 and of F4's 100 restarts; "two-staged" (k-means on
# the Kendall scores of a fit to all subjects, random if a class is too small) is
# the start of F4's 500 subsamples and its K-selection fits. They matter
# separately: on 2026-10-02, locally, log-scale random starts fitted 7 of 7 at
# bwK 8 and 10, while a two-staged K = 2 fit at bwK 8 stopped with "No enough
# points in local window".
#
# Everything else is fixed at the 36-month settings (realdata36_settings.R):
# window [0, 36], bw 3.6, two components, 51-point grid. A cell that stops is
# recorded, not lost, so the reader can report feasibility per bandwidth. Read the
# result with R/combine_realdata36_f0.R, choose bwK, and pass it to F1 and F4 as
# RKCFC36_BWK.
#
# Writes class labels in cohort order only (no CEA value, no identifier).

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[args != "--args"]
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 48L
outdir  <- if (length(args) >= 3) args[3] else file.path(REAL_RESULTS, "F0")
n_sub   <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!is.na(n_sub)) cat(sprintf("*** SUBSAMPLE OF %d SUBJECTS -- SMOKE RUN, NOT REPORTABLE ***\n", n_sub))

cfg <- load_config()
K   <- gint(cfg, "K_true")

grid <- expand.grid(scale = c("log", "raw"), init = c("random", "two-staged"),
                    bwK = c(6, 8, 10, 12), seed_idx = 1:3,
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
grid <- grid[order(grid$scale, grid$init, grid$bwK, grid$seed_idx), , drop = FALSE]
rownames(grid) <- NULL

mine <- grid[(seq_len(nrow(grid)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]
cat(sprintf("F0-36 | task %d of %d | %d of %d cells\n", task_id, n_tasks, nrow(mine), nrow(grid)))
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO cells (there are %d).",
                     "Pass the ORIGINAL array size as n_tasks when re-running one task."),
               task_id, n_tasks, nrow(grid)), call. = FALSE)

data_by_scale <- list()
res <- list(); t0 <- Sys.time()
for (i in seq_len(nrow(mine))) {
  sc <- mine$scale[i]; init <- mine$init[i]; bk <- mine$bwK[i]; si <- mine$seed_idx[i]
  if (is.null(data_by_scale[[sc]])) {
    dd <- truncate_cea36(load_cea36(scale = sc), REAL_TMAX)
    if (!is.na(n_sub)) {
      set.seed(20260909L)
      keep <- sort(sample.int(length(dd$Ly), min(n_sub, length(dd$Ly))))
      dd$Ly <- dd$Ly[keep]; dd$Lt <- dd$Lt[keep]; dd$ids <- dd$ids[keep]
      dd$base <- dd$base[keep, , drop = FALSE]
    }
    data_by_scale[[sc]] <- dd
  }
  d <- data_by_scale[[sc]]
  optns <- realdata_optns(cfg, interval = d$interval, bw = REAL_BW, bwK = bk)
  seed  <- cell_seed(gint(cfg, "seed_master"), sprintf("F0-36|%s|%s|%g", sc, init, bk), si)
  set.seed(seed)
  ts <- Sys.time()
  f <- tryCatch(
    RKCFC(Ly = d$Ly, Lt = d$Lt, K = K, robust = TRUE, optns = optns,
          init.method = init, maxIter = gint(cfg, "maxIter"),
          conv.prop = gnum(cfg, "conv_prop"), LOO = glog(cfg, "LOO"),
          parallel = FALSE, trace = 0),
    error = function(e) structure(conditionMessage(e), class = "f0fail"))
  secs <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
  r <- list(scale = sc, init = init, bwK = bk, seed_idx = si, seed = seed, seconds = secs,
            n_subjects = length(d$Ly), n_obs = length(unlist(d$Lt)))
  if (inherits(f, "f0fail")) {
    r$ok <- FALSE; r$error <- as.character(f)
  } else {
    r$ok <- TRUE; r$error <- NA_character_
    r$g <- as.integer(f$res_clust$g); r$iter <- f$iter; r$conv <- f$conv
    r$stop_reason <- f$stop_reason; r$objective <- f$objective_final
    r$sizes <- as.integer(table(factor(r$g, levels = seq_len(K))))
  }
  res[[i]] <- r
  cat(sprintf("  %d/%d  %s %-10s bwK=%-2g seed=%d  %s  (%.1f min)\n", i, nrow(mine), sc, init, bk, si,
              if (r$ok) sprintf("ok, %s, %d iter, sizes %s", r$stop_reason, r$iter,
                                paste(r$sizes, collapse = "/"))
              else sprintf("FAILED: %s", substr(r$error, 1, 50)), secs / 60))
}

saveRDS(list(cells = res, Tmax = REAL_TMAX, bw = REAL_BW, K = K, n_tasks = n_tasks,
             n_cells = nrow(grid), smoke = !is.na(n_sub),
             code_hash = source_hash(c(file.path("R", "run_realdata36_f0.R"), REAL36_SETTINGS_FILES))),
        file.path(outdir, sprintf("F0_task%04d.rds", task_id)))
cat(sprintf("wrote %s (%.1f min)\n", file.path(outdir, sprintf("F0_task%04d.rds", task_id)),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
