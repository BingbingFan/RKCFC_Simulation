# The CORRECTED silhouette for stored real-data partitions.
#
#   Rscript R/check_realdata_silhouette.R <scale: log|raw> [F4.rds]
#
# WHY THIS EXISTS. Neither F1 nor F4 ever produced a corrected silhouette on the
# real data: F1 left `sil_manuscript`/`sil_pairwise` at their NA template values,
# and the F4 runner called silhouette_score() with the wrong argument list
# (fixed 2026-09-11) inside a tryCatch that swallowed the message, so both F4
# panels printed "-" for every K. The partitions themselves were stored, and the
# coefficient is a deterministic function of a partition, the data and the
# smoothing options, so it is computed here post hoc under EXACTLY the options
# the fits used (realdata_optns: interval [0, 18], bw 3.6, bwK 6), which is the
# whole point of the correction: in release 0.1.x the robust branch omitted the
# covariance bandwidth. Class means are re-estimated by KFPCA per class.
#
# Inputs: F1.rds (the K = 2 fits on the requested scale) and, if given, F4.rds
# (its `kselect` items carry partitions at K = 2..5 where the fit succeeded).
# NO IDENTIFIERS: subjects are the integers 1..n of the analysis set.
#
# THE 36-MONTH ANALYSIS (2026-10-02). Set RKCFC_WINDOW=36 to use the 36-month
# settings (R/realdata36_settings.R: interval [0, 36], bw 3.6, bwK 10) and the
# selected partition in results_real36_sel/ (R/realdata36_select.R); RKCFC_REAL_ROOT
# overrides the results folder. Unset, the script is the 18-month one, unchanged.

source(file.path("R", "methods.R"))
WIN36 <- identical(Sys.getenv("RKCFC_WINDOW"), "36")
if (WIN36) source(file.path("R", "realdata36_settings.R")) else source(file.path("R", "realdata_load.R"))
ROOT  <- Sys.getenv("RKCFC_REAL_ROOT", if (WIN36) "results_real36_sel" else "results_real")

args  <- commandArgs(trailingOnly = TRUE); args <- args[args != "--args"]
scale <- if (length(args) >= 1) args[1] else "log"
f4p   <- if (length(args) >= 2) args[2] else NA_character_
stopifnot(scale %in% c("log", "raw"))

cfg   <- load_config()
d     <- truncate_cea36(load_cea36(scale = scale), REAL_TMAX)
optns <- realdata_optns(cfg, interval = d$interval, bw = REAL_BW, bwK = REAL_BWK)
# KCFC partitions are scored with the NON-robust branch and the fdapace options F1
# used for M2 (its optns_for(FALSE)): the help page says to use the same branch
# that produced the partition, and the corrected coefficient is about matching
# the smoothing of the model scored.
optns_kcfc <- list(dataType = "Sparse", methodMuCovEst = "smooth", nRegGrid = gint(cfg, "nRegGrid"),
                   methodSelectK = gint(cfg, "nK_fit"), useBinnedData = "OFF")
n     <- length(d$Ly)
cat(sprintf("scale %s | n = %d | interval [%g, %g] | bw %g | bwK %s\n", scale, n,
            d$interval[1], d$interval[2], optns$bw, format(optns$bwK)))

sil_of <- function(g, robust = TRUE) {
  g <- as.integer(factor(g))               # classes coded 1..K, as the function expects
  long <- data.frame(id = rep(seq_len(n), lengths(d$Ly)), g = rep(g, lengths(d$Ly)),
                     time = unlist(d$Lt), y = unlist(d$Ly))
  t0 <- Sys.time()
  r <- tryCatch(silhouette_score(long, class = "g", id = "id", y = "y", time = "time",
                                 robust = robust, option = if (robust) optns else optns_kcfc),
                error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error")) list(sil = NA_real_, err = conditionMessage(r), secs = secs)
  else list(sil = as.numeric(r$silhouette), err = NA_character_, secs = secs)
}

rows <- list()
f1 <- readRDS(file.path(ROOT, "F1", "F1.rds"))
for (m in c("M1", "M2")) {
  nm <- paste0(scale, "_", m); g <- f1$fits[[nm]]
  if (is.null(g)) next
  r <- sil_of(g, robust = (m == "M1"))
  rows[[length(rows) + 1]] <- data.frame(source = "F1", partition = nm, K = length(unique(g)),
    sizes = paste(as.integer(table(g)), collapse = "/"), silhouette = r$sil, error = r$err,
    secs = round(r$secs), stringsAsFactors = FALSE)
  cat(sprintf("  %-4s %-8s K=%d sizes %-22s silhouette %s  (%.0f s)%s\n", "F1", nm, length(unique(g)),
              paste(as.integer(table(g)), collapse = "/"),
              if (is.finite(r$sil)) sprintf("%.4f", r$sil) else "-", r$secs,
              if (is.na(r$err)) "" else paste0("  ERROR: ", r$err)))
}
if (!is.na(f4p)) {
  # The combined F4.rds written by combine_realdata_f4.R holds SUMMARIES only (no
  # partitions), so the second argument must be the directory of per-task files
  # F4_task*.rds, whose items carry `g`. A combined file is refused, not misread.
  if (dir.exists(f4p)) {
    tf <- sort(list.files(f4p, pattern = "^F4_task[0-9]+[.]rds$", full.names = TRUE))
    stopifnot(length(tf) > 0)
    parts <- lapply(tf, readRDS)
    sc <- unique(vapply(parts, function(z) as.character(z$scale), ""))
    stopifnot(length(sc) == 1, identical(sc, scale), all(vapply(parts, function(z) z$n, 1L) == n))
    f4 <- list(items = do.call(c, lapply(parts, `[[`, "items")), scale = sc, n = n)
    cat(sprintf("  read %d task files, %d items\n", length(tf), length(f4$items)))
  } else {
    f4 <- readRDS(f4p)
    if (is.null(f4$items)) stop("this is the combined summary file; pass the task-file directory", call. = FALSE)
    stopifnot(identical(as.character(f4$scale), scale), f4$n == n)
  }
  ks <- Filter(function(z) identical(z$kind, "kselect") && !is.null(z$g), f4$items)
  for (z in ks[order(vapply(ks, function(z) as.integer(z$idx)[1], 1L))]) {
    g <- z$g; r <- sil_of(g)
    rows[[length(rows) + 1]] <- data.frame(source = "F4", partition = sprintf("kselect_K%d", as.integer(z$idx)[1]),
      K = length(unique(g)), sizes = paste(as.integer(table(g)), collapse = "/"), silhouette = r$sil,
      error = r$err, secs = round(r$secs), stringsAsFactors = FALSE)
    cat(sprintf("  %-4s %-8s K=%d sizes %-22s silhouette %s  (%.0f s)%s\n", "F4", "kselect", length(unique(g)),
                paste(as.integer(table(g)), collapse = "/"),
                if (is.finite(r$sil)) sprintf("%.4f", r$sil) else "-", r$secs,
                if (is.na(r$err)) "" else paste0("  ERROR: ", r$err)))
  }
}
out <- do.call(rbind, rows)
dir.create(file.path(ROOT, "SIL"), showWarnings = FALSE)
saveRDS(list(scale = scale, table = out, bw = REAL_BW, bwK = REAL_BWK, Tmax = REAL_TMAX,
             f4 = f4p, date = Sys.Date()),
        file.path(ROOT, "SIL", sprintf("silhouette_%s.rds", scale)))
cat(sprintf("wrote %s/SIL/silhouette_%s.rds\n", ROOT, scale))
