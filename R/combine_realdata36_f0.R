# Read stage F0 (36 months): feasibility and partitions by scale, start and bandwidth.
#
#   Rscript R/combine_realdata36_f0.R [f0dir]
#
# Prints, per scale, start type and bwK: how many of the three fits ran and
# converged, the smallest class, the time per fit, and how well the fits agree
# (adjusted Rand index). Then, per scale, how the best-objective partitions at
# different bandwidths agree with each other.
#
# CHOOSE bwK as the smallest value at which, on the PRIMARY (log) scale, every
# fit with EITHER start runs, none returns a tiny class, and the partition does
# not move when the bandwidth moves. F1 needs the random start; F4 needs both.
# Writes F0.rds beside the task files.

source(file.path("R", "realdata36_settings.R"))

args  <- commandArgs(trailingOnly = TRUE)
args  <- args[args != "--args"]
f0dir <- if (length(args) >= 1) args[1] else file.path(REAL_RESULTS, "F0")

files <- sort(list.files(f0dir, pattern = "^F0_task[0-9]+[.]rds$", full.names = TRUE))
if (!length(files)) stop("no F0_task*.rds files in ", f0dir, call. = FALSE)
xs <- lapply(files, readRDS)
n_tasks <- unique(vapply(xs, `[[`, 0, "n_tasks"))
n_cells <- unique(vapply(xs, function(x) if (is.null(x$n_cells)) NA_real_ else x$n_cells, 0))
cells <- unlist(lapply(xs, `[[`, "cells"), recursive = FALSE)
cat(sprintf("F0-36: %d task files (array size %s), %d of %s cells\n",
            length(files), paste(n_tasks, collapse = ","), length(cells), paste(n_cells, collapse = ",")))
if (any(vapply(xs, function(x) isTRUE(x$smoke), TRUE)))
  cat("*** SMOKE RUN (subsample) -- NOT REPORTABLE ***\n")
hashes <- unique(vapply(xs, `[[`, "", "code_hash"))
if (length(hashes) > 1) cat("WARNING: task files come from", length(hashes), "different code versions\n")

ari <- function(a, b) mclust::adjustedRandIndex(a, b)
init_of <- function(r) if (is.null(r$init)) "random" else r$init
tab <- do.call(rbind, lapply(cells, function(r) data.frame(
  scale = r$scale, init = init_of(r), bwK = r$bwK, seed = r$seed_idx, ok = r$ok,
  minutes = round(r$seconds / 60, 1),
  stop = if (r$ok) r$stop_reason else "FAILED",
  iter = if (r$ok) r$iter else NA_integer_,
  sizes = if (r$ok) paste(r$sizes, collapse = "/") else "",
  smallest = if (r$ok) round(min(r$sizes) / sum(r$sizes), 3) else NA_real_,
  J = if (r$ok) r$objective else NA_real_,
  error = if (r$ok) "" else substr(r$error, 1, 40), stringsAsFactors = FALSE)))
tab <- tab[order(tab$scale, tab$init, tab$bwK, tab$seed), ]
cat("\n== every cell ==\n"); print(tab, row.names = FALSE)

cat("\n== by scale, start and bandwidth ==\n")
summ <- NULL
for (sc in unique(tab$scale)) for (it in unique(tab$init)) for (bk in sort(unique(tab$bwK))) {
  rs <- Filter(function(r) r$scale == sc && init_of(r) == it && r$bwK == bk, cells)
  if (!length(rs)) next
  ok <- Filter(function(r) r$ok, rs)
  agree <- if (length(ok) > 1) {
    pr <- utils::combn(length(ok), 2)
    range(apply(pr, 2, function(p) ari(ok[[p[1]]]$g, ok[[p[2]]]$g)))
  } else c(NA, NA)
  row <- data.frame(scale = sc, start = it, bwK = bk, ran = sprintf("%d/%d", length(ok), length(rs)),
                    converged = sum(vapply(ok, function(r) identical(r$stop_reason, "converged"), TRUE)),
                    smallest_class = if (length(ok)) round(min(vapply(ok, function(r) min(r$sizes) / sum(r$sizes), 0)), 3) else NA,
                    median_min = if (length(ok)) round(stats::median(vapply(ok, `[[`, 0, "seconds")) / 60, 1) else NA,
                    ari_between_fits = if (all(is.na(agree))) "" else sprintf("%.3f-%.3f", agree[1], agree[2]),
                    stringsAsFactors = FALSE)
  summ <- rbind(summ, row)
}
print(summ, row.names = FALSE)

cat("\n== does the partition move with the bandwidth? best-objective fit per bwK, either start ==\n")
for (sc in unique(tab$scale)) {
  best <- list()
  for (bk in sort(unique(tab$bwK))) {
    ok <- Filter(function(r) r$ok && r$scale == sc && r$bwK == bk, cells)
    if (length(ok)) best[[as.character(bk)]] <- ok[[which.min(vapply(ok, `[[`, 0, "objective"))]]$g
  }
  if (length(best) > 1) {
    M <- outer(seq_along(best), seq_along(best), Vectorize(function(i, j) ari(best[[i]], best[[j]])))
    dimnames(M) <- list(paste0("bwK", names(best)), paste0("bwK", names(best)))
    cat(sprintf("\n%s scale (adjusted Rand index):\n", sc)); print(round(M, 3))
  }
}

saveRDS(list(cells = tab, summary = summ, Tmax = REAL_TMAX, bw = REAL_BW),
        file.path(f0dir, "F0.rds"))
cat(sprintf("\nwrote %s\n", file.path(f0dir, "F0.rds")))
