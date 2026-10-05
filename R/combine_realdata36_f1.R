# Assemble stage F1 (36 months) from its per-start files into F1.rds.
#
#   Rscript R/combine_realdata36_f1.R [f1dir]
#
# Rebuilds, for each of the four fits, exactly what RKCFC_restarts() returns:
# the best start by objective J, the per-start table (start, seed, J, iter, conv,
# stop_reason, error, aRI_to_best), prop_best (share of successful starts whose
# partition matches the best, aRI > 0.999) and the number of starts that threw.
# Then writes F1.rds in the 18-month layout (table, descriptives, fits, starts,
# bw_months, bwK, Tmax, landmark, K, nstart, primary), so run_realdata36_f2/f3/f5
# and make_realdata_tables.R read it as they read the 18-month one.
#
# REFUSES to combine starts run under different settings (window, bandwidths,
# number of starts) or a mix of smoke and real starts, and reports any start
# that is missing rather than quietly using fewer.

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args  <- commandArgs(trailingOnly = TRUE)
args  <- args[args != "--args"]
f1dir <- if (length(args) >= 1) args[1] else file.path(REAL_RESULTS, "F1")

files <- sort(list.files(f1dir, pattern = "^F1_(raw|log)_M[12]_start[0-9]+[.]rds$", full.names = TRUE))
if (!length(files)) stop("no F1_*_start*.rds files in ", f1dir, call. = FALSE)
xs <- lapply(files, readRDS)

key_of <- function(x) paste(x$Tmax, x$landmark, x$bw_months, x$bwK, x$nstart, x$smoke, x$K, x$primary)
keys <- unique(vapply(xs, key_of, ""))
if (length(keys) > 1)
  stop("start files were run under different settings and cannot be combined:\n  ",
       paste(keys, collapse = "\n  "), "\n(fields: Tmax landmark bw bwK nstart smoke K primary)",
       call. = FALSE)
x1 <- xs[[1]]
NSTART <- x1$nstart; K <- x1$K
FITS <- c("raw_M1", "raw_M2", "log_M1", "log_M2")

missing <- character(0)
for (fk in FITS) {
  have <- vapply(Filter(function(x) x$fit == fk, xs), `[[`, 0L, "start")
  miss <- setdiff(seq_len(NSTART), have)
  if (length(miss)) missing <- c(missing, sprintf("%s starts %s", fk, paste(miss, collapse = ",")))
}
if (length(missing)) {
  cat("MISSING START FILES (re-run those tasks before using F1):\n  ", paste(missing, collapse = "\n  "), "\n", sep = "")
  cat("  task id = (fit index - 1) * ", NSTART, " + start, fit index: raw_M1 1, raw_M2 2, log_M1 3, log_M2 4\n", sep = "")
  if (!nzchar(Sys.getenv("RKCFC36_ALLOW_MISSING", "")))
    stop("refusing to combine an incomplete F1. Set RKCFC36_ALLOW_MISSING=1 to override.", call. = FALSE)
}

ari <- function(a, b) mclust::adjustedRandIndex(a, b)
row0 <- function() data.frame(
  scale = NA_character_, method = NA_character_, metric_scale = NA_character_,
  error = NA_character_, seconds = NA_real_, converged = NA, iter = NA_integer_,
  stop_reason = NA_character_, objective = NA_real_,
  n1 = NA_integer_, n2 = NA_integer_, smallest_frac = NA_real_,
  rre_in = NA_real_, rre_out = NA_real_,
  sil_manuscript = NA_real_, sil_pairwise = NA_real_,
  n_starts = NA_integer_, n_failed_starts = NA_integer_, prop_best = NA_real_,
  stringsAsFactors = FALSE)

rows <- list(); fits <- list(); starts_store <- list(); best_x <- list()
for (fk in FITS) {
  xx <- Filter(function(x) x$fit == fk, xs)
  if (!length(xx)) next
  xx <- xx[order(vapply(xx, `[[`, 0L, "start"))]
  st <- do.call(rbind, lapply(xx, function(x) data.frame(
    start = x$start, seed = x$seed, J = if (x$ok) x$J else NA_real_,
    iter = if (x$ok) x$iter else NA_integer_, conv = if (x$ok) x$conv else NA,
    stop_reason = if (x$ok) x$stop_reason else NA_character_,
    error = if (x$ok) NA_character_ else x$error, stringsAsFactors = FALSE)))
  ok <- which(vapply(xx, `[[`, TRUE, "ok"))
  r <- row0(); r$scale <- xx[[1]]$scale; r$method <- xx[[1]]$method; r$metric_scale <- xx[[1]]$scale
  r$seconds <- sum(vapply(xx, `[[`, 0, "seconds")); r$n_starts <- NSTART
  r$n_failed_starts <- length(xx) - length(ok)
  if (!length(ok)) {
    r$error <- "every start failed; see `starts` for causes"
    st$aRI_to_best <- NA_real_; starts_store[[fk]] <- st; rows[[length(rows) + 1L]] <- r
    cat(sprintf("  %s: FAILED -- all %d starts threw. First: %s\n", fk, length(xx), substr(xx[[1]]$error, 1, 60)))
    next
  }
  best <- ok[which.min(vapply(xx[ok], `[[`, 0, "J"))]
  xb <- xx[[best]]
  st$aRI_to_best <- NA_real_
  st$aRI_to_best[ok] <- vapply(ok, function(i) ari(xx[[i]]$g, xb$g), 0)
  prop_best <- mean(st$aRI_to_best[ok] > 0.999)
  r$converged <- xb$conv; r$iter <- xb$iter; r$stop_reason <- xb$stop_reason; r$objective <- xb$J
  r$n1 <- xb$sizes[1]; r$n2 <- xb$sizes[2]; r$smallest_frac <- min(xb$sizes) / sum(xb$sizes)
  r$rre_in <- xb$rre_in; r$rre_out <- xb$rre_out; r$prop_best <- prop_best
  rows[[length(rows) + 1L]] <- r
  fits[[fk]] <- xb$g; starts_store[[fk]] <- st; best_x[[fk]] <- xb
  cat(sprintf("  %s: best start %d, %s, iter %d, sizes %s, rre_in %.4f, rre_out %s\n", fk, xb$start,
              xb$stop_reason, xb$iter, paste(xb$sizes, collapse = "/"), xb$rre_in,
              if (is.finite(xb$rre_out)) sprintf("%.4f", xb$rre_out) else "NA"))
  cat(sprintf("        %d starts, %d threw, prop_best %.2f, J best %.1f / worst %.1f, %.1f core-h\n",
              NSTART, r$n_failed_starts, prop_best, min(st$J, na.rm = TRUE),
              suppressWarnings(max(st$J, na.rm = TRUE)), r$seconds / 3600))
}
# the cross-scale rows: each log fit scored on the raw scale after back-transforming
for (m in c("M1", "M2")) {
  fk <- paste0("log_", m); xb <- best_x[[fk]]
  if (is.null(xb)) next
  r <- row0(); r$scale <- "log"; r$method <- m; r$metric_scale <- "raw (back-transformed)"
  r$rre_in <- xb$rre_in_rawback
  rows[[length(rows) + 1L]] <- r
  cat(sprintf("  log / %s scored on the RAW scale: rre_in %s\n", m,
              if (is.finite(r$rre_in)) sprintf("%.4f", r$rre_in) else "NA"))
}

# descriptives of the analyzed cohort, from the data (the start files hold none)
d_raw <- truncate_cea36(load_cea36(scale = "raw"), x1$Tmax)
desc  <- cohort_descriptives(d_raw)
if (!x1$smoke && desc$n_subjects != x1$n_subjects)
  stop(sprintf("the start files fitted %d subjects but the cohort on [0, %g] has %d.",
               x1$n_subjects, x1$Tmax, desc$n_subjects), call. = FALSE)

out <- do.call(rbind, rows)
out$n_subjects <- x1$n_subjects
out$smoke      <- x1$smoke
out$code_hash  <- paste(unique(vapply(xs, `[[`, "", "code_hash")), collapse = "+")

saveRDS(list(table = out, descriptives = desc,
             fits = fits,              # partitions only, in cohort order, no PHI
             starts = starts_store,    # per-start J and agreement, no PHI
             bw_months = x1$bw_months, bwK = x1$bwK, Tmax = x1$Tmax,
             landmark = x1$landmark, K = K, nstart = NSTART, primary = x1$primary,
             missing_starts = missing),
        file.path(f1dir, "F1.rds"))
cat(sprintf("\nwrote %s\n", file.path(f1dir, "F1.rds")))
