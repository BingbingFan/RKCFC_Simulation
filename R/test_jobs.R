# Checks the job list and the array slicing, which are what make the run
# restartable. Fast: no fitting, only the grid arithmetic.
#
#   Rscript R/test_jobs.R

source(file.path("config", "validate_config.R"))
source(file.path("R", "jobs.R"))

cfg <- load_config(); fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("  [%s] %-54s %s\n", if (good) "ok" else "FAIL", lab, extra))
}

cat("\n== job lists match the design ==\n")
for (b in c("pilot", "core", "variations")) {
  j <- build_jobs(cfg, b)
  nsc <- length(unique(j$scenario_id))
  ok(sprintf("%s: %d scenarios, %d cells, %d rows", b, nsc, nrow(j),
             expected_rows(cfg, b)),
     nrow(j) > 0 && expected_rows(cfg, b) > 0,
     sprintf("%d methods", expected_rows(cfg, b) / nrow(j)))
}
# COUNTS COME FROM THE CONFIG, NOT FROM A LITERAL. These were hard-coded, and on
# 2026-09-08 the variations check failed for the only bad reason a test can fail:
# a scenario was legitimately ADDED (B12, locating the observation floor) and the
# test still asserted the old total. A test that has to be edited every time the
# design grows teaches people to edit tests, which is the opposite of what it is
# for. What is worth asserting is that the job list matches the scenario tables --
# that build_jobs() and the CSVs cannot drift apart.
for (b in c("pilot", "core", "variations")) {
  tbl <- switch(b, pilot = cfg$pilot, core = cfg$core, variations = cfg$variations)
  # variations runs each config on two arms; the others are one scenario per row.
  want <- nrow(tbl) * if (b == "variations") 2L else 1L
  ok(sprintf("%s job list matches its scenario table", b),
     length(unique(build_jobs(cfg, b)$scenario_id)) == want,
     sprintf("%d rows in config -> %d scenarios", nrow(tbl), want))
}

cat("\n== separations are resolved, not left as levels ==\n")
jc <- build_jobs(cfg, "core")
# Read from the config rather than hardcoding the numbers. An earlier version
# asserted "1.5 and 3" literally and broke the moment d_easy changed, which makes
# a test that fires on a legitimate config edit -- noise, not a guard. The
# property worth asserting is that d_of() RESOLVES the level names into the
# configured values, which is the thing that can actually break.
ok("core d values are the configured hard and easy separations",
   identical(sort(unique(jc$d)), sort(c(gnum(cfg, "d_hard"), gnum(cfg, "d_easy")))),
   sprintf("%s, against configured hard %.1f and easy %.1f",
           paste(sort(unique(jc$d)), collapse = ", "),
           gnum(cfg, "d_hard"), gnum(cfg, "d_easy")))
ok("the two levels are distinct and ordered",
   gnum(cfg, "d_hard") < gnum(cfg, "d_easy"),
   "a 'hard' level at or above the easy one would invert every panel")
jv <- build_jobs(cfg, "variations")
ok("the mean-identical contrast is generated at d = 0",
   all(jv$d[jv$contrast == "identical"] == 0),
   "it is d = 0, not a change of basis")

cat("\n== array slicing: every cell run exactly once, whatever the array size ==\n")
for (nt in c(1L, 7L, 100L)) {
  j <- build_jobs(cfg, "core")
  assigned <- unlist(lapply(seq_len(nt), function(t)
    which((seq_len(nrow(j)) - 1L) %% nt + 1L == t)))
  ok(sprintf("array of %d covers every cell exactly once", nt),
     identical(sort(assigned), seq_len(nrow(j))),
     sprintf("%d cells", nrow(j)))
}

cat("\n== the array can be resized without moving any cell's data ==\n")
j <- build_jobs(cfg, "core")
s7   <- j[(seq_len(nrow(j)) - 1L) %% 7L + 1L == 3L, ]
s100 <- j[(seq_len(nrow(j)) - 1L) %% 100L + 1L == 3L, ]
ok("a cell's identity, and so its seed, does not depend on the array size",
   all(paste(s7$scenario_id, s7$rep) %in% paste(j$scenario_id, j$rep)) &&
     all(paste(s100$scenario_id, s100$rep) %in% paste(j$scenario_id, j$rep)),
   "slices differ, but cell identity and seed do not")

cat("\n== the completeness check cannot be fooled ==\n")
ok("expected is computed from config, not from observed rows",
   expected_rows(cfg, "core") == nrow(build_jobs(cfg, "core")) *
     sum(as.logical(cfg$methods$include_core)),
   sprintf("%d rows for a complete core grid", expected_rows(cfg, "core")))

cat(sprintf("\n%s  %d check%s failed\n\n",
            if (fail == 0L) "JOB LIST AND SLICING ARE CORRECT." else "JOB LIST IS BROKEN.",
            fail, if (fail == 1L) "" else "s"))
quit(status = if (fail == 0L) 0L else 1L)
