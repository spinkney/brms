# Experimental PNUTS engine for brms

This branch starts from [brms PR #1911](https://github.com/paul-buerkner/brms/pull/1911),
commit `36c26dbc6dedbaefb987c82469fc03fa1f42fd73`, and adds `engine = "pnuts"`.
It uses the PR's shared Stan-fit converter. The existing Stan backends remain
available. PNUTS is an experimental sampler, and model support does not imply
that every posterior will mix well.

## Installation

Install this branch into an R library of your choice:

```r
install.packages(c("remotes", "processx", "jsonlite"))
remotes::install_github("spinkney/brms@pnuts-engine-pr1911", upgrade = "never")
```

Separately build/install PNUTS and install `bridgestan==2.9.0` in a Python
environment. This fork contains the adapter, not the PNUTS sampler source.
BridgeStan needs its source distribution and a C++ toolchain to compile Stan
models. The executable and model libraries must use the same OS and CPU
architecture. This adapter targets the PNUTS CLI at commit
`81d401b1f4239afb8b1de1662728b79569ef8ddf` (and compatible versions).

```r
library(brms)
options(
  brms.pnuts.executable = "/path/to/pnuts/build/pnuts",
  brms.pnuts.python = "/path/to/venv/bin/python",
  brms.pnuts.bridgestan = "/path/to/bridgestan-2.9.0"
)
```

`PNUTS_EXECUTABLE`, `PNUTS_PYTHON` and `BRIDGESTAN` are environment-variable
alternatives. Without explicit paths, the adapter looks for `pnuts` and
`python3` on `PATH`, and lets BridgeStan locate its source distribution.
The Python environment must contain BridgeStan and NumPy.

Compiled models are cached under `tools::R_user_dir("brms", "cache")/pnuts`.
Set `options(brms.pnuts.cache_dir = "/persistent/cache")` to change this.
Compilation options can be supplied as
`stan_model_args = list(make_args = c(...), stanc_args = c(...))`.
Code, compiler flags, BridgeStan version/source, architecture and local make
configuration contribute to the cache key.

## Try a fit

```r
set.seed(481)
d <- data.frame(x = rnorm(120))
d$y <- 1.2 - 0.7 * d$x + rnorm(120, sd = 0.8)

fit <- brm(
  y ~ x, data = d, engine = "pnuts",
  chains = 4, cores = 4, iter = 2000, warmup = 1000, seed = 1921,
  control = list(target_fidelity = 0.9, gamma = 0.5)
)
summary(fit)
fixef(fit)
posterior_predict(fit, ndraws = 100)
pp_check(fit)
loo(fit)
posterior::summarise_draws(pnuts_diagnostics(fit))
```

`backend = "pnuts"` is an equivalent spelling. `iter` includes discarded
warmup, as in other brms backends. `cores` limits concurrently running chains;
each native process evaluates gradients directly through BridgeStan.
`thin` is applied when importing completed draws. `save_warmup = TRUE` enables
`pnuts_diagnostics(fit, inc_warmup = TRUE)`. `refresh` is accepted for brms API
compatibility; it does not control the native sampler's progress output.
Use `silent = 0` for compilation and per-chain completion messages.

Multilevel effects, transformed parameters and generated quantities are
imported for brms predictions. `update(fit, newdata = ...)` reuses the cached
compiled model when the Stan code is unchanged and uses a fresh output
directory. `update(stan_fit, engine = "pnuts")` switches a Stan fit to PNUTS.
Prior specifications and brms model generation are unchanged.

## Controls and diagnostics

The defaults are covariance residual adaptation (`adapt = "diag"`) in the
pilot geometry frame (`geometry_frame = "pilot"`), endpoint continuation,
`target_fidelity = 0.9`, `gamma = 0.5`, initial `step_size = 0.1`,
`max_depth = 10`, `max_memory_mb = 1024`, and `chain_timeout = 1800` seconds.
The native pilot automatically uses its bounded-rank geometry path above
its dense-pilot threshold. This adapter does not impose a 64-parameter limit.

Other available native controls are:

- `adapt`: `"none"`, `"step-size"`, `"diag"`, `"diag-fisher"`, `"pilot-only"`.
- `geometry_frame`: `"pilot"` or `"original"`.
- `method`: `"endpoint"`, `"secant"`, `"growth"`, `"soft-growth"`, `"fixed"`.
- `min_size`, `threshold`, `probability_floor`, `power`, `feature_limit`.
- `min_step_size`, `max_step_size`, `metric_regularization`, `metric_blend`.
- `geometry_pilot`, `pilot_iterations`, `pilot_geometry`, `pilot_rank`, `pilot_probes`.
- `state_dependent`, `geometry_diagnostics` (logical flags).

The native CLI validates ranges and combinations. Adaptation needs at least
20 warmup iterations; that minimum is not a recommendation for adequate
warmup. Tuning stays fixed within each tree and freezes before retained draws.
With `adapt = "none"`, supply appropriate fixed settings yourself.

PNUTS fidelity is not a Metropolis acceptance probability. Use
`target_fidelity`, not `adapt_delta`, and `max_depth`, not `max_treedepth`.
Unknown controls error rather than silently being ignored.
`pnuts_diagnostics()` reports fidelity, tree size, root selection,
depth-cap hits, numerical reflections and the other native diagnostics.
`nuts_params()` intentionally errors for PNUTS; no artificial NUTS divergence
or acceptance values are supplied. Inspect R-hat, bulk/tail ESS and the native
diagnostics before using a fit.

Random initializations are drawn independently per chain and checked for
finite target density and gradient. A numeric radius, `init = 0`, an
initialization function, or one named constrained-parameter list per chain
is supported. Partial lists fill unspecified parameters from random values.
No failed sampler run is retried or silently dropped.

## Saved runs and limitations

Pass a fresh `output_dir` to keep raw CSVs, chain logs, JSON requests and run
metadata, initial unconstrained vectors, and per-chain `.adaptation` files.
The manifest records commands, binary/library hashes, initialization work and
per-chain elapsed times. `elapsed_time(fit)` reports total process time; a
separate warmup/sampling time split is unavailable. Model compilation is not
included in those process times.

If any chain fails or times out, the entire fit fails with paths to the logs.
Partial CSV files are never converted into a successful ensemble. Native
numerical reflections are recorded; fatal configuration or pilot failures
remain errors. Raising the timeout does not repair such failures.

The fit object retains draws and native diagnostics across `saveRDS` /
`readRDS`. Raw output paths and the compiled library remain external files;
keep a persistent `output_dir` and cache if you want to preserve them.
Recompiling an absent library is possible with the configured toolchain.

This first adapter supports `algorithm = "sampling"` only. OpenCL,
`future = TRUE`, combining separate fits / `brm_multiple`, exposing Stan
functions, and Stan-function-dependent postprocessing are not supported.
Ordinary brms summaries, posterior predictions, log likelihoods and LOO are
supported. Very large fits also need enough R memory for the returned draws.

## Validation

`tests/testthat/tests.pnuts.R` covers dispatch, controls, initialization specs,
CSV conversion (including one and five chains, saved warmup and thinning),
native diagnostics, and refusal of partial or unfrozen output.

With the toolchain configured, run `tests/local/tests.pnuts.R` for real
Gaussian, multilevel Gaussian and Bernoulli fits; a Gaussian comparison with
CmdStan dense NUTS; predictions and LOO; serialization; updates; and exact
seed replay. These integration checks are opt-in and write their artifacts to
`PNUTS_TEST_OUTPUT` (or a fresh temporary directory). They establish adapter
functionality on those cases, not universal sampler convergence or speed.

Local validation on macOS arm64, R 4.5.2, BridgeStan 2.9.0 and CmdStan 2.39.0
(four chains, 1,000 warmup plus 1,000 retained draws per chain):

| Model | Maximum R-hat | Minimum bulk ESS |
| --- | ---: | ---: |
| Gaussian regression | 1.0010 | 2,399 |
| Gaussian random-intercept model | 1.0059 | 508 |
| Bernoulli regression | 1.0015 | 2,478 |

The Gaussian comparison with dense NUTS had a maximum mean difference of 1.91
combined Monte Carlo standard errors across the reported non-`lp__` variables.
Predictions, log likelihoods, LOO, RDS round trips, refits, exact seed replay,
switching from CmdStan, forced timeouts and single-chain imports were checked.

The selected existing brms regression tests have one failure in
`tests.brmsfit-methods.R:477`: a LOO comparison test uses `$` on the matrix
returned by `loo_compare()`. The same failure reproduces on the untouched PR
head with the same R libraries. Three pre-existing empty tests are skipped.
This is separate from the passing PNUTS-specific checks; a full R CMD check
has not been claimed.
