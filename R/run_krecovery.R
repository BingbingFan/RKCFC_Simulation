# K recovery -- block K.
#
#   Rscript R/run_krecovery.R <task_id> [n_tasks] [outdir] [n_reps] [n_override]
#
# True K in {1, 2, 3}, each fitted at K = 1..K_fit_max, BOTH RKCFC and KCFC.
#
# WHY THIS BLOCK MATTERS MORE THAN ITS SIZE SUGGESTS. The submitted real-data
# silhouettes are
#
#          K = 1     2         3        4
#   RKCFC      0   0.0534   -0.0740  -0.1251
#   KCFC       0   0.6694    0.5491   0.4873
#
# Both select K = 2, but KCFC scores TWELVE AND A HALF TIMES higher on the
# criterion this paper proposes, so on the cohort the paper's own K-selection
# criterion strongly preferred the method the paper argues against.
#
# Simulation is the only place that can be answered, because it is the only place
# the true K is known. Either the silhouette tracks actual recovery accuracy -- in
# which case the real-data gap is a finding about the cohort -- or it rewards
# clusters separated in LEVEL, which is what KCFC produces, over clusters
# separated in SHAPE, which is what RKCFC produces. The second is the argument the
# paper needs, and it is far stronger made voluntarily than extracted under
# questioning.
#
# THREE REPORTING RULES, all of which the submitted paper broke:
#
#   1. Report the full DISTRIBUTION of the selected K, never the mode. A mode of
#      2 hides whether the runner-up was 1 or 5.
#   2. Report the silhouette VALUE, not only its argmax. The gap above is a
#      statement about values.
#   3. Report BOTH criteria. One weak criterion carrying the whole decision is
#      the softest point in the paper, so Calinski-Harabasz on the FPC scores
#      under a robust whitening runs beside it -- `select_K()`, already in the
#      package.
#
# AND THE SILHOUETTE IS REPORTED TWICE, ON PURPOSE. The package's curve-space
# `silhouette_score()` is the quantity the MANUSCRIPT reports, and it is a
# simplified, squared, self-inclusive variant rather than Rousseeuw's: it
# measures distance to the cluster MEAN FUNCTION, scores on squared L2
# deviations, and includes subject i in its own cluster's mean. The squaring is
# not innocuous -- it inflates the coefficient relative to a true distance. So
# this runner
# records the manuscript's quantity (to be comparable with the submitted table)
# AND Rousseeuw's pairwise form on the scores (to be defensible), and the paper
# must say which it is quoting.
#
# K = 1 IS A SPECIAL CASE. The silhouette is UNDEFINED at K = 1; a value of 0
# there is a CONVENTION and must be labelled as one, not printed as a measurement.

source(file.path("R", "methods.R"))
source(file.path("R", "jobs.R"))

args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir  <- if (length(args) >= 3) args[3] else file.path("results", "krecovery")
n_reps  <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
n_over  <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_

cfg <- load_config()
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

scen <- cfg$krecovery
reps <- seq_len(if (is.na(n_reps)) gint(cfg, "n_reps") else n_reps)
jobs <- expand.grid(row = seq_len(nrow(scen)), rep = reps, KEEP.OUT.ATTRS = FALSE)
jobs <- data.frame(scenario_id = scen$scenario_id[jobs$row],
                   arm = scen$arm[jobs$row],
                   d = d_of(cfg, scen$d_level[jobs$row]),
                   K_true = scen$K_true[jobs$row],
                   K_fit_max = scen$K_fit_max[jobs$row],
                   rep = jobs$rep, stringsAsFactors = FALSE)
jobs <- jobs[order(jobs$rep, jobs$scenario_id), ]
mine <- jobs[(seq_len(nrow(jobs)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]

methods <- cfg$methods[cfg$methods$method_id %in% c("M1", "M2"), , drop = FALSE]

cat(sprintf("block krecovery | task %d of %d | %d of %d cells | methods: %s\n",
            task_id, n_tasks, nrow(mine), nrow(jobs),
            paste(methods$method_id, collapse = " ")))
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO cells (the block has %d).",
                     "Pass the ORIGINAL array size as n_tasks when re-running one",
                     "task; with n_tasks = 1 every cell is reassigned."),
               task_id, n_tasks, nrow(jobs)), call. = FALSE)

empty_krow <- function(method = NA_character_, K_fit = NA_integer_,
                       error = NA_character_) {
  data.frame(method = method, K_fit = K_fit, error = error,
             seconds = NA_real_, converged = NA, iter = NA_integer_,
             stop_reason = NA_character_,
             sil_manuscript = NA_real_,   # the package's curve-space variant
             sil_pairwise   = NA_real_,   # Rousseeuw, on scores
             CH             = NA_real_,   # Calinski-Harabasz, robust whitening
             ari = NA_real_, ccr = NA_real_, n_smallest = NA_integer_,
             stringsAsFactors = FALSE)
}

rows <- vector("list", nrow(mine) * nrow(methods) * 5L)
kk <- 0L
t_start <- Sys.time()

for (i in seq_len(nrow(mine))) {
  jb <- mine[i, ]
  arm <- cfg$arms[cfg$arms$arm == jb$arm, ]
  Kmax <- as.integer(jb$K_fit_max)

  # The generator's pi comes from the config and has K_true = 2 entries. For a
  # K_true other than 2 the classes are made equiprobable, which is the neutral
  # choice: imbalance is a SEPARATE factor (the variations block owns it) and
  # confounding it with K recovery would make neither readable.
  pi_k <- rep(1 / jb$K_true, jb$K_true)

  dat <- tryCatch(
    generate_replicate(cfg, arm, d = jb$d, rep = jb$rep,
                       scenario_id = jb$scenario_id, pi = pi_k,
                       n = if (is.na(n_over)) NULL else n_over),
    error = function(e) structure(conditionMessage(e), class = "genfail"))

  if (inherits(dat, "genfail")) {
    for (m in seq_len(nrow(methods))) for (kf in seq_len(Kmax)) {
      kk <- kk + 1L
      rows[[kk]] <- cbind(jb, empty_krow(methods$method_id[m], kf,
                            paste("GENERATION:", as.character(dat))),
                          row.names = NULL)
    }
    next
  }

  for (m in seq_len(nrow(methods))) {
    id  <- methods$method_id[m]
    rob <- id == "M1"
    optns <- build_optns(cfg, rob)

    # Global scores ONCE per (cell, method): select_K() needs a common coordinate
    # system fixed ACROSS the K it compares. Whitening by each partition's own
    # pooled within-cluster scatter is degenerate and silently so -- it gives
    # trace(S^-1 W) = n*p identically for every K, so CH collapses to B alone,
    # which is monotone in K. The metric must not move with the partition.
    sc <- tryCatch(global_scores(dat, cfg, rob, optns), error = function(e) NULL)

    parts <- vector("list", Kmax)
    for (kf in seq_len(Kmax)) {
      r <- empty_krow(id, kf)
      set.seed(cell_seed(gint(cfg, "seed_master"),
                         paste0(jb$scenario_id, "|fit"), jb$rep))
      t0 <- Sys.time()
      f <- tryCatch(
        RKCFC(Ly = dat$Ly, Lt = dat$Lt, K = kf, robust = rob, optns = optns,
              init.method = cfg$global$init_method,
              maxIter = gint(cfg, "maxIter"), conv.prop = gnum(cfg, "conv_prop"),
              LOO = glog(cfg, "LOO"), parallel = FALSE, trace = 0),
        error = function(e) structure(conditionMessage(e), class = "kfail"))
      r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      if (inherits(f, "kfail")) {
        r$error <- as.character(f)
        kk <- kk + 1L; rows[[kk]] <- cbind(jb, r, row.names = NULL); next
      }
      g <- f$res_clust$g
      parts[[kf]] <- g
      r$converged <- f$conv; r$iter <- f$iter; r$stop_reason <- f$stop_reason
      r$n_smallest <- min(tabulate(g, nbins = kf))
      r$ari <- mclust::adjustedRandIndex(dat$labels, g)
      r$ccr <- if (kf == jb$K_true) ccr(dat$labels, g) else NA_real_

      # The MANUSCRIPT's silhouette: the package's curve-space variant, so the
      # simulated values are directly comparable with the submitted table. It is
      # undefined at K = 1; 0 is a convention, recorded as NA here and labelled
      # in the write-up rather than printed as if it were measured.
      if (kf >= 2L) {
        r$sil_manuscript <- tryCatch({
          df <- data.frame(id = rep(seq_along(dat$Ly), lengths(dat$Ly)),
                           g = rep(g, lengths(dat$Ly)),
                           time = unlist(dat$Lt), y = unlist(dat$Ly))
          silhouette_score(df, class = "g", id = "id", y = "y", time = "time",
                           robust = rob, option = optns)$silhouette
        }, error = function(e) NA_real_)
      }
      kk <- kk + 1L; rows[[kk]] <- cbind(jb, r, row.names = NULL)
    }

    # Score-space criteria, computed once over all K on the FIXED whitening.
    if (!is.null(sc)) {
      keep <- which(!vapply(parts, is.null, logical(1)))
      if (length(keep) >= 2L) {
        tb <- tryCatch(select_K(sc, parts[keep], distance = "mahalanobis",
                                robust = rob, silhouette_method = "pairwise"),
                       error = function(e) NULL)
        if (!is.null(tb)) for (j in seq_along(keep)) {
          idx <- which(vapply(rows[seq_len(kk)], function(z)
            identical(z$method, id) && identical(z$K_fit, keep[j]) &&
              identical(z$rep, jb$rep) &&
              identical(z$scenario_id, jb$scenario_id), logical(1)))
          if (length(idx) == 1L) {
            rows[[idx]]$CH <- tb$CH[j]
            rows[[idx]]$sil_pairwise <- tb$silhouette[j]
          }
        }
      }
    }
  }
  if (i %% 2L == 0L || i == nrow(mine))
    cat(sprintf("  %d/%d cells, %.1f min elapsed\n", i, nrow(mine),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

out <- do.call(rbind, rows[seq_len(kk)])

# PROVENANCE. Two task files in one directory can differ in their slicing or in
# the code that produced them, at exactly 100 percent completeness, and nothing
# else here can see it. The combiner refuses to merge files that disagree.
out$n_tasks   <- n_tasks
out$code_hash <- source_hash(file.path("R", "run_krecovery.R"))
f <- file.path(outdir, sprintf("krecovery_task%04d.rds", task_id))
saveRDS(out, f)
cat(sprintf("wrote %s  (%d rows, %.1f min)\n", f, nrow(out),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
