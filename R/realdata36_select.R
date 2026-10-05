# SELECT THE REPORTED 36-MONTH PARTITION, by the rule of Section 2.3: the partition
# with the smallest objective J over every random start of the PRIMARY fit, reported
# with J and the share of starts that reached it.
#
#   Rscript R/realdata36_select.R [in_root] [out_root]
#     in_root   cluster output, default REAL_RESULTS (never modified)
#     out_root  default <in_root>_sel
#
# THE PRIMARY FIT comes from realdata36_settings.R (REAL_PRIMARY: "raw_M1" since
# 2026-10-03, RKCFC on CEA on its original scale at bwK 24; "log_M1" for the
# log-scale run of 2026-10-02, now the sensitivity analysis, which this script
# selected with the same rule). The primary fit is started 120 times with identical
# settings: the 20 starts of F1 and the 100 restarts of F4. On the log scale the
# smallest J among F1's 20 (1079.48, 822/583) was not the smallest overall: six F4
# restarts reached 1077.95 (579/826). Section 2.3 reports the smallest J, so that
# partition is the one reported.
#
# WHAT IT WRITES. out_root/F1/F1.rds: F1.rds of in_root with the primary partition
# replaced, its labels oriented to agree with F1's (class numbering unchanged), its
# table row rescored exactly as F1 scores a fit (rre_in, held-out rre_out, and for a
# log-scale primary rre_in back-transformed), and its starts table extended to all
# 120 starts. The other three fits are copied unchanged. Class FPCAs are refitted on
# the selected partition with the F1 settings; as a check, the same rescoring is
# applied to F1's own best start and compared with the values F1 stored. Class
# names are set here, keyed to the data so a relabelling cannot swap them.
#
# Prints summaries only; no identifier is printed or written.

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))
if (!requireNamespace("mclust", quietly = TRUE)) stop("package `mclust` is required.", call. = FALSE)

args <- commandArgs(trailingOnly = TRUE); args <- args[args != "--args"]
IN   <- if (length(args) >= 1) args[1] else REAL_RESULTS
OUT  <- if (length(args) >= 2) args[2] else paste0(IN, "_sel")
dir.create(file.path(OUT, "F1"), recursive = TRUE, showWarnings = FALSE)

PRIM   <- REAL_PRIMARY
SCALE  <- sub("_.*$", "", PRIM)          # "raw" or "log"
METHOD <- sub("^.*_", "", PRIM)          # "M1"
cfg   <- load_config()
K     <- gint(cfg, "K_true")
f1    <- readRDS(file.path(IN, "F1", "F1.rds"))
if (!identical(f1$primary, PRIM) || !isTRUE(f1$bwK == REAL_BWK))
  stop(sprintf("F1 under %s was run with primary %s at bwK %s; the settings say %s at bwK %s",
               IN, f1$primary, format(f1$bwK), PRIM, format(REAL_BWK)), call. = FALSE)
d_raw <- truncate_cea36(load_cea36(scale = "raw"), REAL_TMAX)
d_log <- truncate_cea36(load_cea36(scale = "log"), REAL_TMAX)
d_p   <- if (SCALE == "raw") d_raw else d_log
n     <- length(d_p$Ly)
optns <- realdata_optns(cfg, interval = d_raw$interval, bw = REAL_BW, bwK = REAL_BWK)
cat(sprintf("primary fit %s (scale %s), bwK %g, window [0, %g], %d patients\n", PRIM, SCALE, REAL_BWK, REAL_TMAX, n))

# --- every start of the primary fit ----------------------------------------------
s1 <- lapply(list.files(file.path(IN, "F1"), sprintf("^F1_%s_start[0-9]+\\.rds$", PRIM), full.names = TRUE), readRDS)
s1 <- Filter(function(z) isTRUE(z$ok) && length(z$g) == n, s1)
ff <- list.files(file.path(IN, "F4"), "^F4_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(ff)) stop("no F4 task files under ", file.path(IN, "F4"), call. = FALSE)
z4 <- readRDS(ff[1])
stopifnot(z4$scale == SCALE, z4$bwK == REAL_BWK, z4$bw == REAL_BW, z4$Tmax == REAL_TMAX, z4$K == K)
no_error <- function(z) is.null(z$error) || isTRUE(is.na(z$error))
rs <- Filter(function(z) identical(z$kind, "restart") && no_error(z) && length(z$g) == n,
             unlist(lapply(ff, function(f) readRDS(f)$items), recursive = FALSE))
pool <- c(lapply(s1, function(z) list(source = "F1", start = z$start, g = as.integer(z$g),
                                      J = as.numeric(z$J), stop = z$stop_reason)),
          lapply(rs, function(z) list(source = "F4", start = z$idx, g = as.integer(z$g),
                                      J = as.numeric(z$J), stop = z$stop_reason)))
if (!length(pool)) stop("no successful start of ", PRIM, " found", call. = FALSE)
J    <- vapply(pool, `[[`, 0, "J")
best <- which.min(J)
g    <- pool[[best]]$g
g_f1 <- as.integer(f1$fits[[PRIM]])
# LABEL ORIENTATION (K = 2). Original scale: class 2 is the HIGH-CEA class (the one
# with the higher median peak CEA), so that every "class 2 against class 1" hazard
# ratio downstream (F2, F3, Table S15) reads high-CEA against low-CEA; F1's own
# numbering put the 249-patient high class first. Log scale: as F1's, the published
# orientation (class 1 flat, class 2 falling).
if (SCALE == "raw") {
  pk_or <- vapply(d_raw$Ly, max, 0)
  if (which.max(tapply(pk_or, g, stats::median)) == 1L) g <- 3L - g
  if (which.max(tapply(pk_or, g_f1, stats::median)) == 1L) g_f1 <- 3L - g_f1
} else {
  if (mean(g == g_f1) < 0.5) g <- 3L - g
}
ari  <- vapply(pool, function(p) mclust::adjustedRandIndex(p$g, g), 0)
hit  <- ari > 0.999
stops <- vapply(pool, function(p) as.character(p$stop), "")
cat(sprintf("%s starts: %d (F1 %d, F4 %d); stop reasons: %s\n", PRIM, length(pool), length(s1), length(rs),
            paste(sprintf("%s %d", names(table(stops)), table(stops)), collapse = ", ")))
cat(sprintf("smallest J %.4f (%s start %s, %s), sizes %s; reached by %d of %d starts (%.0f%%)\n",
            J[best], pool[[best]]$source, pool[[best]]$start, pool[[best]]$stop,
            paste(table(g), collapse = "/"), sum(hit), length(pool), 100 * mean(hit)))
cat(sprintf("F1's best of its 20: J %.4f, sizes %s; adjusted Rand index with the selected %.3f\n",
            min(vapply(s1, function(z) as.numeric(z$J), 0)), paste(table(g_f1), collapse = "/"),
            mclust::adjustedRandIndex(g, g_f1)))
pa <- if (length(pool) > 1) utils::combn(length(pool), 2, function(ij) mclust::adjustedRandIndex(pool[[ij[1]]]$g, pool[[ij[2]]]$g)) else NA
cat(sprintf("pairwise adjusted Rand index over the %d starts: median %.3f (range %.3f to %.3f); %d distinct partitions\n",
            length(pool), stats::median(pa), min(pa), max(pa),
            length(unique(vapply(pool, function(p) paste(match(p$g, unique(p$g)), collapse = ""), "")))))

# --- score a partition exactly as F1 scores a fit -----------------------------
score <- function(gg) {
  fc  <- univar_FPCA_k(d_p$Ly, d_p$Lt, K, gg, TRUE, optns)
  rin <- reconstruction_errors(d_p, gg, fc, TRUE, optns)$rre_in
  sp  <- holdout_split(d_p, gnum(cfg, "holdout_frac"), gint(cfg, "holdout_min_ni"), 1L)
  fh  <- univar_FPCA_k(sp$fit$Ly, sp$fit$Lt, K, gg, TRUE, optns)
  rout <- holdout_error(sp, gg, fh, TRUE, optns)
  rbk <- if (SCALE == "log") reconstruction_errors(d_raw, gg, fc, TRUE, optns, back = exp, fwd = log)$rre_in else NA_real_
  c(rre_in = rin, rre_out = rout, rre_in_rawback = rbk)
}
tb <- f1$table
i  <- which(tb$scale == SCALE & tb$method == METHOD & tb$metric_scale == SCALE)
ib <- which(tb$scale == SCALE & tb$method == METHOD & tb$metric_scale != SCALE)   # the back-transformed row (log only)
stopifnot(length(i) == 1, length(ib) <= 1)
t0 <- Sys.time()
chk <- score(g_f1)
cat(sprintf("\ncheck on F1's own best start (refit vs stored): rre_in %.4f vs %.4f, rre_out %.4f vs %.4f%s\n",
            chk["rre_in"], tb$rre_in[i], chk["rre_out"], tb$rre_out[i],
            if (length(ib)) sprintf(", back-tr. %.4f vs %.4f", chk["rre_in_rawback"], tb$rre_in[ib]) else ""))
sc <- score(g)
cat(sprintf("selected partition: rre_in %.4f, rre_out %.4f%s (%.1f min)\n", sc["rre_in"], sc["rre_out"],
            if (length(ib)) sprintf(", back-transformed rre_in %.4f", sc["rre_in_rawback"]) else "",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# --- assignment margins of the selected partition -----------------------------
# For each patient, (reconstruction error under the other class - under its own) /
# under its own, from the class FPCAs refitted on the selected partition. On a
# partition the assignment rule leaves unchanged every margin is >= 0; a negative
# one would mean the refit does not reproduce the fit, and is counted.
fc_sel <- univar_FPCA_k(d_p$Ly, d_p$Lt, K, g, TRUE, optns)
err_k  <- sapply(seq_len(K), function(k) vapply(seq_len(n), function(ii) {
  p <- univar_predicted_yhat(fc_sel[[k]], d_p$Ly[[ii]], d_p$Lt[[ii]], TRUE, optns)
  sum((d_p$Ly[[ii]] - p$fittedX)^2) }, 0))
own    <- err_k[cbind(seq_len(n), g)]
oth    <- err_k[cbind(seq_len(n), 3L - g)]
margin <- (oth - own) / own
cat(sprintf("\nmargins: %d negative; median %.3f (IQR %.3f-%.3f); below 0.10: %.1f%% (class 1 %.1f%%, class 2 %.1f%%)\n",
            sum(margin < 0), stats::median(margin), stats::quantile(margin, .25), stats::quantile(margin, .75),
            100 * mean(margin < 0.10), 100 * mean(margin[g == 1] < 0.10), 100 * mean(margin[g == 2] < 0.10)))

# --- write the F1 record with the selected partition --------------------------
f1s <- f1
f1s$fits[[PRIM]] <- g
sz <- as.integer(table(factor(g, levels = seq_len(K))))
tb$seconds[i] <- NA; tb$converged[i] <- identical(pool[[best]]$stop, "converged"); tb$iter[i] <- NA
tb$stop_reason[i] <- pool[[best]]$stop; tb$objective[i] <- J[best]
tb$n1[i] <- sz[1]; tb$n2[i] <- sz[2]; tb$smallest_frac[i] <- min(sz) / n
tb$rre_in[i] <- sc["rre_in"]; tb$rre_out[i] <- sc["rre_out"]
if (length(ib)) tb$rre_in[ib] <- sc["rre_in_rawback"]
tb$n_starts[i] <- length(pool); tb$n_failed_starts[i] <- 0L; tb$prop_best[i] <- mean(hit)
f1s$table <- tb
s1t <- f1$starts[[PRIM]]
s1_starts <- vapply(s1, function(z) as.integer(z$start), 0L)
f1s$starts[[PRIM]] <- data.frame(
  source = vapply(pool, `[[`, "", "source"),
  start = vapply(pool, function(p) as.integer(p$start), 0L),
  seed = c(s1t$seed[match(s1_starts, s1t$start)], rep(NA_integer_, length(rs))),
  J = J, iter = c(s1t$iter[match(s1_starts, s1t$start)], rep(NA_integer_, length(rs))),
  conv = vapply(pool, function(p) identical(p$stop, "converged"), TRUE),
  stop_reason = stops, error = NA_character_, aRI_to_best = ari, stringsAsFactors = FALSE)
# CLASS NAMES, read by the table and figure scripts, keyed to the data so that a
# relabelling cannot swap them.
#   original scale: the partition is a level split (the 24-month fit of 2026-10-03: the
#   smaller class has median peak CEA 28 against 4 ng/mL and preoperative CEA 22 against
#   3), so the classes are "high-CEA" and "low-CEA", keyed to the median peak CEA.
#   log scale: one class starts high at surgery and falls below the other by month 2-3
#   (preoperative CEA reproduces the split at AUC 0.82), so "falling" and "flat", keyed
#   to the median preoperative CEA, as on 2026-10-02.
if (SCALE == "raw") {
  pk <- vapply(d_raw$Ly, max, 0)
  hi <- which.max(tapply(pk, g, stats::median))
  f1s$class_names <- ifelse(seq_len(K) == hi, "high-CEA", "low-CEA")
} else {
  pre <- vapply(d_raw$Ly, function(y) y[1], 0)
  hi  <- which.max(tapply(pre, g, stats::median))
  f1s$class_names <- ifelse(seq_len(K) == hi, "falling", "flat")
}
cat(sprintf("class names: %s\n", paste(sprintf("%d = %s (n %d)", seq_len(K), f1s$class_names, sz), collapse = ", ")))
f1s$selection <- list(rule = sprintf("smallest J over all %s starts (F1 %d + F4 %d)", PRIM, length(s1), length(rs)),
                      primary = PRIM, scale = SCALE, bwK = REAL_BWK,
                      n_starts = length(pool), J_best = J[best], n_reached = sum(hit),
                      pairwise_ari = c(median = stats::median(pa), min = min(pa), max = max(pa)),
                      source = pool[[best]]$source, start = pool[[best]]$start, stop = pool[[best]]$stop,
                      ari_with_f1_best = mclust::adjustedRandIndex(g, g_f1),
                      margin = margin, check_f1_refit = chk)
saveRDS(f1s, file.path(OUT, "F1", "F1.rds"))
cat(sprintf("\nwrote %s\n", file.path(OUT, "F1", "F1.rds")))
