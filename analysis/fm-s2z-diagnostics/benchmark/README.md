# Reproducing the FM comparison

Run these commands from the brms repository root. They intentionally write
posterior CSVs and fitted objects outside the repository.

```sh
git clone https://github.com/adamlauretig/ny_r_talk.git /tmp/ny_r_talk
git -C /tmp/ny_r_talk checkout 9258f2b3ae0a948841f6ace266f130069ee3f66f
mkdir -p /tmp/fm-benchmark
export FM_BENCHMARK_DIR=/tmp/fm-benchmark
export NY_R_TALK_DIR=/tmp/ny_r_talk
Rscript analysis/fm-s2z-diagnostics/benchmark/calibrate_rms_prior.R
Rscript analysis/fm-s2z-diagnostics/benchmark/check_transform.R
Rscript analysis/fm-s2z-diagnostics/benchmark/fit_source.R
Rscript analysis/fm-s2z-diagnostics/benchmark/fit_semiorthogonal.R
Rscript analysis/fm-s2z-diagnostics/benchmark/analyze_fits.R
```

`fit_source.R` uses the tracked syntax-modernized source model. Its likelihood
and priors match `stan_fm_1.stan`; deprecated array syntax is updated, factor
priors are vectorized, and posterior-predictive RNG is omitted from both timed
models. Timings in `results/sampler.csv` are CmdStan's maximum reported time
over chains, not an independently measured end-to-end wall clock.

The historical native-S2Z draft in the report predates the semi-orthogonal
implementation and is retained only as context. It cannot be reproduced from
the final branch without checking out that intermediate working tree.

The 200,000-draw RMS calibration uses `set.seed(1)` and direct centered
Gaussian factors. `calibrate_rms_prior.R` reproduces the full-precision
lognormal parameters used by the benchmark.

## Expanded stress-test harness

The `expanded` directory contains a separate, configurable comparison across
rank, sparsity, degree imbalance, tied singular values, weak interactions, and
the full-rank determinant edge case. It does not reproduce the committed
results in `../results`; it is a harness for running additional experiments.
All generated executables, draws, and summaries are written outside the
repository:

```sh
export FM_EXPANDED_DIR=/tmp/fm-expanded
Rscript analysis/fm-s2z-diagnostics/benchmark/expanded/run.R
Rscript analysis/fm-s2z-diagnostics/benchmark/expanded/analyze.R
```

Use `FM_EXPANDED_SCENARIOS` and `FM_EXPANDED_MODELS` to select comma-separated
subsets. Sampling settings can be changed with `FM_EXPANDED_CHAINS`,
`FM_EXPANDED_CORES`, `FM_EXPANDED_WARMUP`, and `FM_EXPANDED_SAMPLING`.
