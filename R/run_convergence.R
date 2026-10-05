# THE OBJECTIVE TRACE -- block J.
#
#   Rscript R/run_convergence.R <task_id> [n_tasks] [outdir] [n_reps] [n_override]
#
# WHY THIS BLOCK EXISTS, AND WHY IT IS NOT OPTIONAL.
#
# There is no convergence argument to give, and the manuscript
# says so outright: the reassignment step is monotone BY CONSTRUCTION -- a
# subject moves only to a class with a strictly smaller residual sum -- but the
# RE-ESTIMATION step is not. For ordinary K-means it minimises the same objective
# in closed form, which is exactly what makes Lloyd's algorithm descend and
# therefore terminate. That property does not survive substituting Kendall FPCA:
# the Kendall estimator is not the minimiser of a sum of squared reconstruction
# errors, and the local-linear smoothing of the mean and covariance surfaces adds
# a bandwidth-dependent bias that does not vanish in J. So J CAN RISE across a
# refit, and no reclassification threshold repairs that.
#
# What replaces the proof is measurement, and this block is that measurement:
#
#   * the share of re-estimation steps on which J increases
#   * the DISTRIBUTION of those increases, not merely whether any occur
#   * how often a fit ends in a CYCLE rather than at a settled partition
#   * and, as a structural check, that the reassignment step descends in every
#     recorded row -- if it ever does not, the implementation contradicts the
#     one convergence property the method does have, and everything else here is
#     moot
#
# BOTH METHODS ARE TRACED, and that is deliberate. The manuscript states that the
# original K-centers functional clustering algorithm has the same gap. That is
# currently an assertion; running M2 through the same instrument turns it into a
# measurement, and a concession that comes with evidence about the comparator is
# far stronger than one that does not.
#
# WHAT THIS REPLACES. The claim in the handoff -- 37 percent of re-estimation
# steps raising J, one fit in eight ending in a cycle -- came from a SMOKE RUN of
# 12 replicates at n = 150 and d = 3, from a scratchpad script that no longer
# exists, at a separation the design later rejected as saturated. Reportable
# numbers cannot rest on that. The core grid cannot rescue it either: it stores
# `objective_final`, one scalar per fit, and discards the per-iteration trace.

source(file.path("R", "methods.R"))

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[args != "--args"]
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir  <- if (length(args) >= 3) args[3] else file.path("results", "convergence")
n_reps  <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
n_over  <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_

cfg <- load_config()
K   <- gint(cfg, "K_true")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
reps <- seq_len(if (is.na(n_reps)) gint(cfg, "n_reps") else n_reps)

# Four arms, not two. The placeholder asks for the distribution "by scenario", so
# breadth is the point here rather than depth, and a trace costs nothing beyond
# the fit that was going to happen anyway: arm 1 the Gaussian reference, arm 2
# heavy-tailed but symmetric and finite-variance, arm 4 infinite variance, arm 6
# contamination in the curves rather than the scores.
ARMS <- c(1L, 2L, 4L, 6L)

jobs <- expand.grid(arm = ARMS, rep = reps, KEEP.OUT.ATTRS = FALSE)
jobs$d <- gnum(cfg, "d_hard")
jobs$scenario_id <- sprintf("J_a%d", jobs$arm)
jobs <- jobs[order(jobs$rep, jobs$arm), ]
mine <- jobs[(seq_len(nrow(jobs)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]

METHODS <- c(M1 = "RKCFC", M2 = "KCFC")

cat(sprintf("block convergence | task %d of %d | %d of %d cells | %d methods\n",
            task_id, n_tasks, nrow(mine), nrow(jobs), length(METHODS)))
if (!is.na(n_over))
  cat(sprintf("*** n OVERRIDDEN to %d. SMOKE RUN -- not reportable. ***\n", n_over))
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO cells (this block has %d).",
                     "Pass the ORIGINAL array size as n_tasks when re-running one task."),
               task_id, n_tasks, nrow(jobs)), call. = FALSE)

# The row shape, defined once, so a failure row binds with a success row.
empty_jrow <- function(method = NA_character_, error = NA_character_) {
  data.frame(
    method = method, error = error, seconds = NA_real_,
    iter = NA_integer_, converged = NA, stop_reason = NA_character_,
    cycle_length = NA_integer_, J_final = NA_real_,
    # the trace summaries
    n_reest = NA_integer_,          # re-estimation steps actually observed
    n_reest_up = NA_integer_,       # of those, how many RAISED J
    frac_reest_up = NA_real_,
    max_rise_rel = NA_real_,        # largest single rise, as a share of J
    med_rise_rel = NA_real_,        # median rise among the steps that rose
    total_rise_rel = NA_real_,      # all rises summed, as a share of the final J
    reassign_descends_all = NA,     # the structural check
    n_reassign_rows = NA_integer_,
    ari = NA_real_,
    stringsAsFactors = FALSE)
}

# Summarise one fit's trace. `J_reestimation` is already J_before(b) minus
# J_after(b-1), so it isolates the REFIT from the reassignment that preceded it,
# and a POSITIVE value means the refit raised J.
summarise_trace <- function(f) {
  r <- empty_jrow()
  tr <- f$objective
  if (is.null(tr) || !nrow(tr)) return(r)      # K = 1 and single-iteration fits
  r$n_reassign_rows <- nrow(tr)
  # `reassignment_descends` must hold in EVERY row. It is the one property the
  # method provably has, so a FALSE here is an implementation defect and not a
  # finding about the algorithm.
  r$reassign_descends_all <- all(tr$reassignment_descends, na.rm = TRUE)

  # ALIGNMENT MATTERS HERE. J_reestimation[b] = J_before[b] - J_after[b-1], so a
  # row b that carries a value describes the refit that ran BETWEEN reassignment
  # b-1 and reassignment b, and the level it started from is J_after[b-1] -- not
  # J_after[b], and not J_before[b]. Getting this off by one would divide each
  # rise by the wrong denominator and quietly rescale every number in the block.
  idx  <- which(!is.na(tr$J_reestimation))
  d    <- tr$J_reestimation[idx]
  base <- tr$J_after[idx - 1L]
  r$n_reest <- length(d)
  if (!length(d)) return(r)

  # Relative, not absolute: J scales with n and with the noise level, so an
  # absolute rise is not comparable across arms and would make the contaminated
  # arm look worse for a reason that has nothing to do with convergence.
  rel <- d / pmax(abs(base), .Machine$double.eps)
  up  <- d > 1e-8
  r$n_reest_up    <- sum(up)
  r$frac_reest_up <- mean(up)
  r$max_rise_rel  <- if (any(up)) max(rel[up]) else 0
  r$med_rise_rel  <- if (any(up)) stats::median(rel[up]) else 0
  r$total_rise_rel <- if (any(up)) sum(d[up]) / max(abs(f$objective_final), 1e-12) else 0
  r
}

rows <- vector("list", nrow(mine) * length(METHODS))
kk <- 0L
t_start <- Sys.time()

for (i in seq_len(nrow(mine))) {
  jb  <- mine[i, ]
  arm <- cfg$arms[cfg$arms$arm == jb$arm, ]
  dat <- tryCatch(
    generate_replicate(cfg, arm, d = jb$d, rep = jb$rep,
                       scenario_id = jb$scenario_id,
                       n = if (is.na(n_over)) NULL else n_over),
    error = function(e) structure(conditionMessage(e), class = "genfail"))

  if (inherits(dat, "genfail")) {
    for (m in names(METHODS)) { kk <- kk + 1L
      rows[[kk]] <- cbind(jb, empty_jrow(m, paste("GENERATION:", as.character(dat))),
                          row.names = NULL) }
    next
  }

  for (m in names(METHODS)) {
    robust <- m == "M1"
    # ONE COMMON SEED for both methods, as everywhere else in this study: the
    # methods are being COMPARED, so a shared starting state pairs them.
    set.seed(cell_seed(gint(cfg, "seed_master"),
                       paste0(jb$scenario_id, "|fit"), jb$rep))
    optns <- build_optns(cfg, robust)
    t0 <- Sys.time()
    f <- tryCatch(
      RKCFC(Ly = dat$Ly, Lt = dat$Lt, K = K, robust = robust, optns = optns,
            init.method = cfg$global$init_method,
            maxIter = gint(cfg, "maxIter"), conv.prop = gnum(cfg, "conv_prop"),
            LOO = glog(cfg, "LOO"), parallel = FALSE, trace = 0),
      error = function(e) structure(conditionMessage(e), class = "jfail"))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (inherits(f, "jfail")) {
      r <- empty_jrow(m, as.character(f)); r$seconds <- secs
      kk <- kk + 1L; rows[[kk]] <- cbind(jb, r, row.names = NULL); next
    }

    r <- summarise_trace(f)
    r$method <- m; r$seconds <- secs
    r$iter <- f$iter; r$converged <- f$conv; r$stop_reason <- f$stop_reason
    r$cycle_length <- f$cycle_length; r$J_final <- f$objective_final
    r$ari <- mclust::adjustedRandIndex(dat$labels, f$res_clust$g)
    kk <- kk + 1L; rows[[kk]] <- cbind(jb, r[, names(empty_jrow())], row.names = NULL)
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
out$code_hash <- source_hash(file.path("R", "run_convergence.R"))

f <- file.path(outdir, sprintf("convergence_task%04d.rds", task_id))
saveRDS(out, f)
cat(sprintf("wrote %s  (%d rows, %.1f min)\n", f, nrow(out),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
