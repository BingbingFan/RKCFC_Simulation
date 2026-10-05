# RKCFC_Simulation

R code for the simulation study and the real-data analysis of the manuscript

> Fan B, et al. Robust K-Centers Functional Clustering Model for Sparsely Observed
> Non-Gaussian Longitudinal Data. Under review.

The method itself is the R package
[RKCFC](https://github.com/BingbingFan/RKCFC). The results in the manuscript were
produced with RKCFC 0.2.0 and R 4.5.1.

## Requirements

```r
install.packages(c("remotes", "fdapace", "fda", "mclust", "MASS", "lcmm", "survival", "mgcv"))
remotes::install_github("BingbingFan/RKCFC@v0.2.0")
```

Run every script from the root of this repository, because the scripts find each
other and the configuration through relative paths. `R/generate.R` stops if the
package it loads is not the univariate RKCFC.

## Layout

| Folder | Contents |
|---|---|
| `config/` | The design as tables: the six distributional arms (`arms.csv`), mean functions and eigenfunctions (`basis.csv`), score dispersion matrices (`Z.csv`), the seven methods (`methods.csv`), the five eigenspace estimators (`estimators.csv`), global parameters with the master seed (`params_global.csv`), and the scenario tables of the pilot (`scenarios_p.csv`), core grid (`scenarios_a.csv`), variations block (`scenarios_b.csv`), and class-number study (`scenarios_k.csv`). `validate_config.R` checks them. |
| `R/` | Data generation (`generate.R`), the methods (`methods.R`), the metrics (`metrics.R`), one runner and one combiner per study, the table generator (`make_tables.R`), and the real-data scripts. |

## How a study is run

Each runner divides its grid of cells into `n_tasks` slices and runs slice
`task_id`. Every cell takes its seed from the master seed in
`config/params_global.csv` and its own identity, so the result does not depend on
how the grid is divided, and any slice or cell can be rerun alone. To rerun one
slice, pass the original `n_tasks`.

For example, the core grid in one process:

```sh
Rscript R/run_cell.R core 1 1 results/core
Rscript R/combine.R results/core results/core/core
```

or in 100 slices, run `Rscript R/run_cell.R core <i> 100 results/core` for
`i` = 1, ..., 100, in any order and on any number of cores, then combine as above.

The checks to run before a full study:

```sh
Rscript config/validate_config.R
Rscript R/test_jobs.R
Rscript R/test_generate.R
Rscript R/test_methods.R
Rscript R/test_separation.R
```

When the studies have run, the tables of the manuscript are written by

```sh
Rscript R/make_tables.R results tables
```

The whole simulation is 22,510 fits. A fit of RKCFC took a median of about 150
single-core seconds and a fit of KCFC about 25, so the full study needs several
hundred core-hours and is best run in slices on several cores.

## The simulation studies

What each study varies, and where it is reported. Arms are the score or error
distributions, and `d` is the mean offset between the two classes.

| Study | What it varies | Arms | d | Methods | Replicates | Fits | Manuscript |
|---|---|---|---|---|---|---|---|
| Separation pilot | Robust Mahalanobis separation at seven values of d | All six | 0 to 4 | None, score distributions only | 20 | 840 | Table S2 |
| Pilot at full sample size | Candidate separations | Gaussian, t(1.5) | 0, 1.5, 2 | RKCFC, KCFC | 20 | 240 | Table S3 |
| Core grid | Arm and separation | All six | 1.5, 2 | All seven | 100 | 8,400 | Tables 2 and 3, S8, S10 |
| Variations block | Twelve one-factor changes to the reference design | Fleishman(0, 4), t(1.5) | 1.5, 0 in one row | RKCFC, KCFC, FPCA k-means, Kendall k-means | 50 | 4,800 | Table S4 |
| Estimation study | Five eigenspace estimators fitted per class under the true partition | All six | 1.5 | Estimators only | 50 | 3,000 | Table 4 |
| Class-number study | True K of 1, 2, and 3, with K fitted from 1 to 5 | Fleishman(0, 4) | 1.5 | RKCFC, KCFC | 50 | 1,500 | Table S9 |
| Bandwidth and mean estimator | Subject-level bandwidth halved and doubled, and an M-estimator of the mean | Gaussian, both Fleishman, t(1.5) | 1.5 | RKCFC | 50 | 800 | Table S5 |
| Estimators inside the loop | Kendall's tau, spatial sign, normal scores | Gaussian, Fleishman(0, 4), t(1.5) | 1.5 | RKCFC | 50 | 450 | Table S6 |
| Initialization | Six starting schemes | Gaussian, t(1.5) | 1.5 | RKCFC | 50 | 600 | Table S12 |
| Retained components | Number of components fixed at 1 to 4, or chosen by the fraction of variance explained | Fleishman(0, 4) | 1.5 | RKCFC | 50 | 400 | Table S13 |
| Objective trace | The objective recorded at every iteration | Gaussian, Fleishman(0, 4), t(1.5), contamination | 1.5 | RKCFC, KCFC | 50 | 400 | Table S11 |
| Design checks | Separation ladder, convergence tolerance, initialization, retained components | Gaussian, t(1.5) | 1.5, 0 to 3 in the ladder | RKCFC, KCFC | 20 | 1,080 | Settings in Table S1 |

What each study uses. Every study reads `arms.csv`, `basis.csv`, `Z.csv`, and
`params_global.csv`; the second column lists what it uses beyond those. Results
folders are the ones `make_tables.R` reads under its results root.

| Study | Also uses | Runner | Results folder | Combiner |
|---|---|---|---|---|
| Separation pilot | The grid of d in `R/run_decisions.R`, and `R/separation.R` | `Rscript R/run_decisions.R delta <task_id> <n_tasks> results/D_delta 20` | `D_delta` | `Rscript R/combine_decisions.R results/D_delta 20` |
| Pilot at full sample size | `scenarios_p.csv`, `methods.csv` | `Rscript R/run_cell.R pilot <task_id> <n_tasks> results/pilot` | `pilot` | `Rscript R/combine.R results/pilot results/pilot/pilot` |
| Core grid | `scenarios_a.csv`, `methods.csv` | `Rscript R/run_cell.R core <task_id> <n_tasks> results/core` | `core` | `Rscript R/combine.R results/core results/core/core` |
| Variations block | `scenarios_b.csv`, `methods.csv` | `Rscript R/run_cell.R variations <task_id> <n_tasks> results/variations` | `variations` | `Rscript R/combine.R results/variations results/variations/variations` |
| Estimation study | `estimators.csv` | `Rscript R/run_estim.R <task_id> <n_tasks> results/estim` | `estim` | `Rscript R/combine_estim.R results/estim` |
| Class-number study | `scenarios_k.csv`, `methods.csv`, and the third class in `basis.csv` | `Rscript R/run_krecovery.R <task_id> <n_tasks> results/krecovery` | `krecovery` | `Rscript R/combine_krecovery.R results/krecovery` |
| Bandwidth and mean estimator | The levels in `R/run_settings.R` | `Rscript R/run_settings.R smoothing <task_id> <n_tasks> results/S_smoothing` | `S_smoothing` | `Rscript R/combine_settings.R results/S_smoothing` |
| Estimators inside the loop | The levels in `R/run_settings.R` | `Rscript R/run_settings.R estimators <task_id> <n_tasks> results/S_estimators` | `S_estimators` | `Rscript R/combine_settings.R results/S_estimators` |
| Initialization | The levels in `R/run_settings.R` | `Rscript R/run_settings.R init <task_id> <n_tasks> results/S_init` | `S_init` | `Rscript R/combine_settings.R results/S_init` |
| Retained components | The levels in `R/run_settings.R` | `Rscript R/run_settings.R components <task_id> <n_tasks> results/S_components` | `S_components` | `Rscript R/combine_settings.R results/S_components` |
| Objective trace | Nothing beyond the shared tables | `Rscript R/run_convergence.R <task_id> <n_tasks> results/convergence` | `convergence` | `Rscript R/combine_convergence.R results/convergence` |
| Design checks | The levels in `R/run_decisions.R` | `Rscript R/run_decisions.R <factor> <task_id> <n_tasks> results/D_<factor> 20`, with factor `d`, `convprop`, `init`, or `nkfit` | `D_d`, `D_convprop`, `D_init`, `D_nkfit` | `Rscript R/combine_decisions.R results/D_<factor> 20` |

All studies share the data generator in `R/generate.R`, the methods in
`R/methods.R`, and the metrics in `R/metrics.R`.

`R/make_tables.R` writes one file per table. The files map to the numbering of
the manuscript as follows.

| File | Table | File | Table |
|---|---|---|---|
| `tab01.tex` | 1, design | `tabS5.tex` | S4, variations |
| `tab02.tex` | 2, adjusted Rand index on the core grid | `tabS8.tex` | S5, bandwidth and mean estimator |
| `tab05.tex` | 3, reconstruction error | `tabS9.tex` | S6, estimators inside the loop |
| `tab03.tex` | 4, estimation study | `tabStudies.tex` | S7, overview of the studies |
| `tabS1.tex` | S1, parameters | `tabCCR.tex` | S8, correct classification rate |
| `tabS2.tex` | S2, separation pilot | `tab04.tex` | S9, number of classes |
| `tabS10.tex` | S3, pilot at full sample size | `tabS3.tex` | S10, convergence |
| `tabS4.tex` | S11, objective trace | `tabS7.tex` | S12, initialization |
| `tabS6.tex` | S13, retained components | | |

## The real-data analysis

The cohort data cannot be shared publicly because of patient privacy, and they
are available from the corresponding author on reasonable request. The scripts
read the file named by the environment variable `RKCFC_DATA`. The loader in
`R/realdata_load.R` returns only the analysis set, with an integer patient ID,
and checks that none of the identifier columns of the source file is returned.
No result file of the
real-data analysis is included here, because those files hold a class label for
every patient.

`R/realdata36_settings.R` holds the settings: the 36-month clustering window, the
36-month landmark, the subject-level bandwidth of 3.6 months, the Kendall-surface
bandwidth, the reported partition, and the time origin of the reported survival
estimate. Each can be overridden by an environment variable. The defaults give
the primary analysis, RKCFC on CEA on its original scale with a Kendall-surface
bandwidth of 24 months. The log-scale sensitivity analysis uses `RKCFC36_BWK=10`,
`RKCFC36_PRIMARY=log_M1`, and `RKCFC36_RESULTS=results_real36`.

The steps, in order:

| Step | Scripts | Output |
|---|---|---|
| 1. Bandwidth sweep | `run_realdata36_f0.R`, then `combine_realdata36_f0.R` | The smallest bandwidth at which the iteration converges and starts agree |
| 2. Longer runs on the original scale | `run_realdata36_raw2s.R` at 200 and 500 iterations, then `combine_realdata36_raw2s.R` | Convergence at each candidate bandwidth |
| 3. Fits from random starts | `run_realdata36_f1.R`, one start per task, then `combine_realdata36_f1.R` | Both methods on both scales, 20 starts each |
| 4. Stability | `run_realdata36_f4.R`, 100 restarts and 500 subsample refits, then `combine_realdata_f4.R` | Agreement, clusterwise Jaccard similarity, assignment margins |
| 5. Reported partition | `realdata36_select.R <in_root> <out_root>` | The partition with the smallest objective over all starts |
| 6. Survival and diagnostics | `run_realdata36_f2.R`, `run_realdata36_f3.R`, `run_realdata36_f5.R`, with the readers `read_realdata_f1.R` and `read_realdata_f5.R`, and `check_realdata_silhouette.R` for the silhouette by number of classes | Cox models under both time origins, marker distribution, silhouettes |
| 7. Tables and figures | `make_realdata_tables.R`, `fig_realdata_traj.R`, `fig_realdata_km.R` | Table 5, Tables S14 to S16, Figures 1, S1, and S2 |

Each script states its arguments in its first lines.

## License

GPL-3, as for the RKCFC package.
