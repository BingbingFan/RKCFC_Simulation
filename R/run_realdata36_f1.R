# STAGE F1 (36 months) -- ONE RANDOM START OF ONE FIT PER TASK.
#
#   Rscript R/run_realdata36_f1.R <task_id> [outdir] [n_sub]
#
#   task_id  1..80. Four fits -- raw_M1, raw_M2, log_M1, log_M2 (M1 = RKCFC,
#            M2 = KCFC) -- times NSTART = 20 random starts. Task t runs start
#            ((t-1) %% 20) + 1 of fit ((t-1) %/% 20) + 1.
#   outdir   default <REAL_RESULTS>/F1
#   n_sub    fit a random subsample of this many subjects -- SMOKE ONLY; every
#            output is stamped and the combiner and the tables refuse it.
#
# THE SAME PROCEDURE AS THE 18-MONTH F1, SPREAD OVER AN ARRAY. The 18-month F1
# ran RKCFC_restarts() over 20 starts for each fit in one ~19-hour job. Here each
# start is its own task, seeded exactly as RKCFC_restarts() seeds it (start b of
# a fit uses seed0 + b, seed0 from cell_seed("F1|<scale>|<method>")), and
# R/combine_realdata36_f1.R rebuilds what RKCFC_restarts() returns -- the best
# start by objective, prop_best, the per-start table -- and writes F1.rds in the
# 18-month format, so F2, F3 and F5 read it unchanged.
#
# Each task also scores its own start (in-sample and held-out relative
# reconstruction error, and for log fits the in-sample error back-transformed to
# the raw scale), because those use the fitted object; the combiner then reports
# the numbers of the best start. A start that throws is saved with its error, so
# the fit survives as long as any start succeeds.
#
# Writes class labels in cohort order only (no CEA value, no identifier).

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[args != "--args"]
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
outdir  <- if (length(args) >= 2) args[2] else file.path(REAL_RESULTS, "F1")
n_sub   <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

NSTART <- as.integer(Sys.getenv("RKCFC36_NSTART", "20"))   # 20 for the real run
FITS   <- c("raw_M1", "raw_M2", "log_M1", "log_M2")
if (is.na(task_id) || task_id < 1L || task_id > length(FITS) * NSTART)
  stop(sprintf("task_id must be 1..%d", length(FITS) * NSTART), call. = FALSE)
fit_key <- FITS[(task_id - 1L) %/% NSTART + 1L]
b       <- (task_id - 1L) %% NSTART + 1L
sc      <- sub("_.*$", "", fit_key)
m       <- sub("^.*_", "", fit_key)
robust  <- m == "M1"

cfg  <- load_config()
K    <- gint(cfg, "K_true")
NREG <- gint(cfg, "nRegGrid")

d_raw <- truncate_cea36(load_cea36(scale = "raw"), REAL_TMAX)
d_log <- truncate_cea36(load_cea36(scale = "log"), REAL_TMAX)
if (!is.na(n_sub)) {
  cat(sprintf("*** SUBSAMPLE OF %d SUBJECTS -- SMOKE RUN, NOT REPORTABLE ***\n", n_sub))
  set.seed(20260909L)
  keep <- sort(sample.int(length(d_raw$Ly), min(n_sub, length(d_raw$Ly))))
  sub <- function(d) { d$Ly <- d$Ly[keep]; d$Lt <- d$Lt[keep]
                       d$ids <- d$ids[keep]; d$base <- d$base[keep, , drop = FALSE]; d }
  d_raw <- sub(d_raw); d_log <- sub(d_log)
}
d  <- if (sc == "raw") d_raw else d_log
iv <- d_raw$interval

optns <- if (robust) realdata_optns(cfg, interval = iv, bw = REAL_BW, bwK = REAL_BWK) else
  list(dataType = "Sparse", methodMuCovEst = "smooth", nRegGrid = NREG,
       methodSelectK = gint(cfg, "nK_fit"), useBinnedData = "OFF")

seed0 <- cell_seed(gint(cfg, "seed_master"), paste0("F1|", sc, "|", m), 1L)
seed  <- as.integer(seed0) + b
cat(sprintf("F1-36 | task %d | %s start %d of %d | seed %d | n = %d | window [0, %g] | bw %g | bwK %g\n",
            task_id, fit_key, b, NSTART, seed, length(d$Ly), REAL_TMAX, REAL_BW, REAL_BWK))

set.seed(seed)
t0 <- Sys.time()
f <- tryCatch(
  RKCFC(Ly = d$Ly, Lt = d$Lt, K = K, robust = robust, optns = optns,
        maxIter = gint(cfg, "maxIter"), conv.prop = gnum(cfg, "conv_prop"),
        LOO = glog(cfg, "LOO"), parallel = FALSE, trace = 0, init.method = "random"),
  error = function(e) structure(conditionMessage(e), class = "ffail"))
fit_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

out <- list(fit = fit_key, scale = sc, method = m, start = b, seed = seed, nstart = NSTART,
            seconds = fit_secs, n_subjects = length(d$Ly), smoke = !is.na(n_sub),
            Tmax = REAL_TMAX, landmark = REAL_LANDMARK, bw_months = REAL_BW, bwK = REAL_BWK,
            K = K, primary = REAL_PRIMARY)
if (inherits(f, "ffail")) {
  out$ok <- FALSE; out$error <- as.character(f)
  cat(sprintf("  FAILED after %.1f min: %s\n", fit_secs / 60, out$error))
} else {
  g <- as.integer(f$res_clust$g)
  out$ok <- TRUE; out$error <- NA_character_; out$g <- g
  out$J <- f$objective_final; out$iter <- f$iter; out$conv <- f$conv
  out$stop_reason <- f$stop_reason
  out$sizes <- as.integer(table(factor(g, levels = seq_len(K))))
  out$rre_in <- tryCatch(reconstruction_errors(d, g, f$FPCA_C, robust, optns)$rre_in,
                         error = function(e) NA_real_)
  sp <- holdout_split(d, gnum(cfg, "holdout_frac"), gint(cfg, "holdout_min_ni"), 1L)
  out$rre_out <- tryCatch({
    fh <- univar_FPCA_k(sp$fit$Ly, sp$fit$Lt, K, g, robust, optns)
    holdout_error(sp, g, fh, robust, optns)
  }, error = function(e) NA_real_)
  # A log-scale fit scored on the RAW scale after back-transforming (yit_lg is
  # the natural log of yit), so all four fits can be compared on one scale.
  out$rre_in_rawback <- if (sc == "log")
    tryCatch(reconstruction_errors(d_raw, g, f$FPCA_C, robust, optns,
                                   back = exp, fwd = log)$rre_in,
             error = function(e) NA_real_) else NA_real_
  cat(sprintf("  ok after %.1f min: %s, %d iter, sizes %s, J %.4f, rre_in %.4f, rre_out %s\n",
              fit_secs / 60, out$stop_reason, out$iter, paste(out$sizes, collapse = "/"),
              out$J, out$rre_in, if (is.finite(out$rre_out)) sprintf("%.4f", out$rre_out) else "NA"))
}
out$code_hash <- source_hash(c(file.path("R", "run_realdata36_f1.R"), REAL36_SETTINGS_FILES))

fn <- file.path(outdir, sprintf("F1_%s_start%02d.rds", fit_key, b))
saveRDS(out, fn)
cat(sprintf("wrote %s (%.1f min in all)\n", fn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
