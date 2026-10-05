# THE 36-MONTH REAL-DATA ANALYSIS: its settings, stated once.
#
# Every *realdata36* script sources THIS file instead of realdata_load.R. It
# sources realdata_load.R for the loader and helpers (load_cea36, truncate_cea36,
# realdata_optns, peak_cea, cohort_descriptives) and then REPLACES the 18-month
# constants. realdata_load.R itself is not edited, so the 18-month pipeline and
# its results in results_real/ are untouched (frozen copies: backup_realdata_18m/).
#
# WHY 36 MONTHS IS BACK. The 18-month window was adopted after the F0/F0B/F0C
# sweeps found the 36-month window unfittable at every Kendall-surface bandwidth
# tried. All three sweeps were run on the RAW scale, where the K-centers loop
# isolates a handful of extreme patients and the small class has almost no
# observation pairs late in follow-up, so the smoother stops ("No enough points
# in local window"). The primary analysis is on the LOG scale, which was never
# tried on 36 months. The submitted analysis did run on 36 months, on rescaled
# log CEA, with a Kendall-surface bandwidth of 10 months (and four components,
# a 721-point grid; it stopped at 100 iterations without converging). Local
# checks on 2026-10-02, log scale, one random start each: bwK 10 fitted and
# converged 3 of 3 (classes about 600/800, starts agreeing at aRI 0.93-0.99);
# bwK 8 fitted 3 of 3 (1 converged, 2 cycled); bwK 6 converged 1 of 3, the other
# two still running after two hours. The two-staged start failed at once at bwK
# 8, 10 and 12. Stage F0 measures all of this properly on the cluster.
#
# THE SETTINGS, each overridable from the environment so a sweep or a re-run at
# another value needs no code change (for example RKCFC36_BWK=8; one
# value per variable, no commas):
#
#   REAL_TMAX = 36       the clustering window, as in the submitted analysis.
#   REAL_LANDMARK = 36   tied to the window: membership is not known before the
#                        trajectory is complete, so survival starts at its end.
#                        1259 of 1405 patients are alive at month 36,
#                        with 127 deaths after it (18 months: 1341 and 187).
#   REAL_BW = 3.6        subject-level Nadaraya-Watson bandwidth, 0.1 of the
#                        window -- the simulation's rule (bw = 0.1 on [0,1]) and
#                        the value first chosen for [0,36]. F0 showed it does not
#                        affect feasibility.            env: RKCFC36_BW
#   REAL_BWK = 24        Kendall-surface bandwidth, the parameter that binds.
#                        WAS 10 (the submitted fit's value; the log-scale run of
#                        2026-10-02 used it). 24 SINCE 2026-10-03: the longer
#                        original-scale runs (RAW2S, REALDATA36_RUNBOOK.md §8)
#                        showed that on the original scale the iteration does
#                        not settle at 10 or 12 months, cycles at 18, and
#                        converges at 24 (11 of 11 fits that ran, two partitions
#                        differing on 15 patients).        env: RKCFC36_BWK
#   REAL_PRIMARY         "raw_M1" (RKCFC on CEA on its original scale) SINCE
#                        2026-10-03; WAS "log_M1". The
#                        log-scale analysis at bwK 10 is now the sensitivity
#                        analysis (its results stay in results_real36 and
#                        results_real36_sel).           env: RKCFC36_PRIMARY
#   REAL_RESULTS         where every 36-month stage reads and writes, relative
#                        to study/ unless absolute. "results_real36_raw24" for
#                        the primary run; "results_real36" holds the log-scale
#                        run, the F0 sweep and RAW2S.   env: RKCFC36_RESULTS
#
#   To re-run or re-read the LOG-SCALE analysis, set all three:
#     RKCFC36_BWK=10 RKCFC36_PRIMARY=log_M1 RKCFC36_RESULTS=results_real36
#
# Settings that do not depend on the window (K = 2, two components, 51-point
# grid, B-spline basis, maxIter, convergence rule, LOO off, seeds) come from
# config/params_global.csv exactly as in the 18-month analysis.

source(file.path("R", "realdata_load.R"))

.num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (!nzchar(v)) return(default)
  x <- suppressWarnings(as.numeric(v))
  if (!is.finite(x) || x <= 0) stop(name, " must be a positive number, got '", v, "'", call. = FALSE)
  x
}

REAL_TMAX     <- 36
REAL_LANDMARK <- REAL_TMAX
REAL_BW       <- .num_env("RKCFC36_BW", 3.6)
REAL_BWK      <- .num_env("RKCFC36_BWK", 24)
REAL_PRIMARY  <- Sys.getenv("RKCFC36_PRIMARY", "raw_M1")
REAL_RESULTS  <- Sys.getenv("RKCFC36_RESULTS", "results_real36_raw24")
# THE TIME ORIGIN OF THE REPORTED SURVIVAL ESTIMATE, since 2026-10-03:
# "surgery" (time from surgery on all patients, the submitted design; the reported
# fully adjusted HR of the 24-month original-scale partition is 2.24, 1.67-3.01) with
# the month-36 landmark analysis as the sensitivity analysis (1.55, 1.02-2.37). Was
# "landmark" from 2026-09-09 to 2026-10-03, which removed the guarantee-time
# bias by design; both are now reported. F2 and F3 fit
# both designs whatever this says; it decides which one the tables and Figure 2
# lead with.                                        env: RKCFC36_COX_ORIGIN
REAL_COX_ORIGIN <- Sys.getenv("RKCFC36_COX_ORIGIN", "surgery")
if (!REAL_COX_ORIGIN %in% c("surgery", "landmark"))
  stop("RKCFC36_COX_ORIGIN must be surgery or landmark, got '", REAL_COX_ORIGIN, "'", call. = FALSE)

if (!REAL_PRIMARY %in% c("log_M1", "raw_M1"))
  stop("RKCFC36_PRIMARY must be log_M1 or raw_M1, got '", REAL_PRIMARY, "'", call. = FALSE)

cat(sprintf("[36-month settings] window [0, %g], landmark %g, bw %g, bwK %g, primary %s, survival origin %s, results %s\n",
            REAL_TMAX, REAL_LANDMARK, REAL_BW, REAL_BWK, REAL_PRIMARY, REAL_COX_ORIGIN, REAL_RESULTS))

# The files every 36-month result is stamped with, so a result records the
# settings that produced it.
REAL36_SETTINGS_FILES <- c(file.path("R", "realdata_load.R"),
                           file.path("R", "realdata36_settings.R"))
