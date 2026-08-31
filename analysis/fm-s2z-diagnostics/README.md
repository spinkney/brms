# Semi-orthogonal factorization machine diagnostics

## Executive summary

The centered semi-orthogonal implementation reproduces the identifiable fit
of the unconstrained factorization machine in
[`adamlauretig/ny_r_talk`](https://github.com/adamlauretig/ny_r_talk/blob/master/stan_fm_1.stan)
while removing its continuous factor-rotation and permutation symmetries.
Across both simulated 100 by 20 data sets in that repository:

- all 2,000 cell predictors and all 2,000 double-centered interactions have
  rank-normalized split R-hat below `1.008`, with none above `1.01`;
- the ordered spectrum has maximum R-hat `1.010` and `1.004`;
- both fits have zero divergences, zero maximum-treedepth hits, and minimum
  E-BFMI above `0.66`;
- posterior-mean predictor surfaces correlate above `0.99998` with the source
  model and differ by only `0.022--0.024` RMSE; their mean absolute difference
  is `0.021--0.023` after division by the quadrature sum of posterior SDs;
- recovery of the simulated mean and interaction is effectively unchanged;
- mean leapfrog counts fall by about 82%, while CmdStan's maximum per-chain
  total time is 2.4% and 1.4% lower than the source fits.

This is a much better computational result than the earlier native
`sum_to_zero_matrix` draft: the semi-orthogonal version has about 25% lower
maximum per-chain time on both data sets and requires roughly one quarter as
many leapfrogs.

Raw reflector coordinates still have very poor cross-chain diagnostics. That
does not contradict the invariant results: the Householder representation
retains prior-only radial auxiliaries and the SVD retains paired column-sign
symmetries. Raw seeds are therefore excluded from saved draws by default.
Diagnose `sifm`, the interaction surface, and predictions instead.

## Implemented model

For field levels `i` and `j`, the full linear predictor under the default brms
intercept is

\[
\eta_{ij}=\alpha+a_i+b_j+
s\sqrt{N J}\,q_i^\mathsf{T}\operatorname{diag}(d)r_j,
\]

where

\[
\mathbf 1^\mathsf{T}Q=0,\quad Q^\mathsf{T}Q=I,
\qquad
\mathbf 1^\mathsf{T}R=0,\quad R^\mathsf{T}R=I.
\]

The Stan Householder-reflector transform is adapted from Seth Axen's
MIT-licensed
[`stan_semiorthogonal_transforms`](https://github.com/sethaxen/stan_semiorthogonal_transforms)
(revision `8f47b5760a896e363ee5e15966c4f96e017c42eb`). Its construction is based
on [Stewart](https://epubs.siam.org/doi/10.1137/0717034) and
[Nirwan and Bertschinger](https://proceedings.mlr.press/v97/nirwan19a.html).
The brms adaptation constructs Haar-Stiefel frames in the `(N - 1)` and
`(J - 1)` dimensional centered subspaces and applies an isometric Helmert lift
back to level space. A full-rank frame selected for `special = 1` instead uses
Haar measure on the relevant `SO(K)` component. The upstream copyright and
full MIT permission notice are kept in the Stan source chunk.

The normalized singular spectrum is parameterized by

\[
p\sim\operatorname{Dirichlet}(1),\qquad
e_r=\sum_{j=r}^K\frac{p_j}{j},\qquad d_r=\sqrt{e_r}.
\]

Consequently `d[1] > ... > d[K] > 0` and `sum(d^2) = 1`. This map covers the
entire ordered energy simplex without another scale redundancy. It also makes

\[
H=s\sqrt{NJ}\,Q\operatorname{diag}(d)R^\mathsf T,
\]

the double-centered interaction, with

\[
s^2=\frac{1}{NJ}\lVert H\rVert_F^2,
\]

so `sdfm` is the exact full-table RMS interaction, and `sifm` reports the
identifiable normalized singular spectrum.

Ordering removes permutations and continuous rotations almost surely.
Simultaneously flipping column `k` of both frames remains a discrete symmetry.
Near-equal singular values can also be weakly identified. For a full-rank
square centered frame, the implementation restricts at most one side to
determinant `+1`; the other side retains both signs, so no interaction surface
is excluded. If both centered margins are square at rank `K`, their remaining
relative determinant sign is a genuine disconnected component of the
full-rank interaction. It is selected by a step in the unrestricted frame's
last scalar reflector coordinate, so chains should be initialized and checked
in both determinant components when that sign is uncertain. This is especially
visible in the `2 by 2, K = 1` edge case.

## Correction applied to the upstream reflector code

The cited upstream `special = 1, N = K` branch was not special orthogonal. It
left `factors[N,N]` equal to zero, inspected `diagonal(factors)` instead of
`tau`, and assigned `+/-1` where a scalar reflector requires `tau` equal to 0
or 2. In the upstream function names, the correction is:

```stan
factors[N, N] = 1;
tau[N] = 1 - determinant_from_tau(head(tau, N - 1));
```

Fixed-parameter checks for odd and even square dimensions gave determinant
`+1` and orthogonality error below `1.8e-15`. Centered-frame tests gave column
sums below `5.0e-16` before CSV rounding.

## Models and prior calibration

The source model fits

\[
y_{ij}\sim\mathcal N(a_i+b_j+g_i^\mathsf{T}h_j,1)
\]

with independent `normal(0, 3)` main effects and independent standard-normal
factor coordinates. The comparison canonicalizes every source draw into an
intercept, two zero-sum main effects, and a double-centered interaction before
computing diagnostics.

The likelihood and complete-grid model space are the same. The prior match is
approximate for every canonical component, not only the spectrum:
canonicalization makes the intercept, main effects, and interaction mutually
dependent; it also adds Gaussian-product terms to the intercept and main
effects. The source singular values follow a product-Wishart law, for which an
exact match is not available from the implementation's elementary priors. The
benchmark therefore uses:

- fixed main-effect scales of 3, close to source marginal SDs `3.0261` and
  `2.9321`, and fixed residual SD 1;
- `normal(0, 0.736545993)` for the absorbed intercept, matching its variance
  but approximating its non-Gaussian distribution;
- `lognormal(0.7678554, 0.08103256)` for `sdfm`, calibrated to the source
  model's centered full-table RMS prior with 200,000 product-Wishart draws
  generated from direct centered Gaussian factors using `set.seed(1)`;
- a fixed `Dirichlet(1)` prior on internal `zfm_spectrum`, which maps to the
  reported invariant spectrum `sifm` but is not a configurable `sifm` prior.

These deliberate joint-prior differences are important, although 2,000
observations make the reported likelihood-dominated fits nearly identical.

## Reproducibility

- brms branch: `feature/s2z-factorization-machines`
- brms base commit: `6ae3c9b1b9bbd3f19723865da5723b84edcf5ee5`;
  the implementation is tracked on the branch above
- source repository commit: `9258f2b3ae0a948841f6ace266f130069ee3f66f`
- R 4.5.2, CmdStanR 0.9.0.9000, CmdStan 2.39.0, Apple arm64
- data: both objects in `simulated_data.rdata`
- dimensions: 2,000 complete-grid observations, 100 by 20 fields, rank 5
- chains: 4; warmup/sampling: 1,000/1,000 per chain
- `adapt_delta = 0.8`, `max_treedepth = 20`
- seeds: 123 and 216, shared with the corresponding source fit
- compilation time excluded

After sampling, the fixed internal frame, spectrum, and native sum-to-zero
priors were moved from direct model-target statements into brms's saved
`lprior` accumulator. This bookkeeping correction leaves the posterior target
unchanged, but floating-point evaluation order can prevent byte-for-byte MCMC
trajectory reproduction. Later formula-validation and prediction changes also
do not alter the fitted Stan model.

The second source data set was relabeled from 20 by 100 to 100 by 20. The
model and priors are symmetric in the two fields.

The tracked [`benchmark`](benchmark/README.md) directory contains the exact
syntax-modernized source model, fitting and analysis scripts, prior-calibration
audit, and direct transform check. Machine-readable summary outputs are in
[`results`](results/). The syntax modernization preserves the posterior but
updates deprecated Stan syntax, vectorizes factor priors, and removes generated
posterior-predictive RNG from the timed programs.

## Sampler diagnostics

| Data | Model | Divergences | Depth hits | Min E-BFMI | Mean leapfrogs | LF 90% | LF 99% | Max chain total (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Simulation 1 | Modernized source | 0 | 0 | 0.784 | 216.3 | 511 | 1023 | 51.0 |
| Simulation 1 | Native S2Z draft | 0 | 0 | 0.799 | 162.2 | 255 | 511 | 66.6 |
| Simulation 1 | Semi-orthogonal | 0 | 0 | 0.668 | 39.0 | 63 | 63 | 49.8 |
| Simulation 2 | Modernized source | 0 | 0 | 0.782 | 175.0 | 447 | 895 | 41.7 |
| Simulation 2 | Native S2Z draft | 0 | 0 | 0.806 | 127.4 | 127 | 127 | 54.7 |
| Simulation 2 | Semi-orthogonal | 0 | 0 | 0.663 | 31.0 | 31 | 31 | 41.1 |

Relative to the modernized source model, the semi-orthogonal parameterization
cuts mean leapfrog work by 82.0% and 82.3%. Maximum per-chain total time changes
by -2.4% and -1.4%, so the defensible timing conclusion is parity, not a broad
speed claim. Maximum per-chain sampling time is 5.1% higher in Simulation 1
and 23.7% lower in Simulation 2. These are CmdStan-reported chain times, not an
independently measured end-to-end wall clock.

## Diagnostics for identifiable quantities

| Data | Quantity | Model | Max R-hat | R-hat > 1.01 | Min bulk ESS | Median bulk ESS | Min tail ESS | Median bulk ESS/s |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Simulation 1 | Total surface | Modernized source | 1.0040 | 0 | 3,455 | 4,541 | 2,808 | 187.9 |
| Simulation 1 | Total surface | Semi-orthogonal | 1.0058 | 0 | 3,213 | 6,515 | 1,933 | 256.5 |
| Simulation 1 | Interaction | Modernized source | 1.0038 | 0 | 3,372 | 4,429 | 2,647 | 183.2 |
| Simulation 1 | Interaction | Semi-orthogonal | 1.0055 | 0 | 2,679 | 6,233 | 2,111 | 245.4 |
| Simulation 2 | Total surface | Modernized source | 1.0039 | 0 | 3,101 | 4,779 | 2,746 | 233.1 |
| Simulation 2 | Total surface | Semi-orthogonal | 1.0069 | 0 | 3,433 | 6,188 | 2,175 | 395.6 |
| Simulation 2 | Interaction | Modernized source | 1.0033 | 0 | 2,831 | 4,603 | 2,599 | 224.5 |
| Simulation 2 | Interaction | Semi-orthogonal | 1.0071 | 0 | 3,329 | 5,939 | 2,115 | 379.6 |

Median bulk ESS/second improves by 31--70% for invariant surfaces. Tail ESS is
lower in absolute terms for the semi-orthogonal fits, though every minimum is
above 1,900.

### Ordered spectrum

| Data | True RMS | Posterior RMS | True normalized singular values | Posterior normalized singular values | Max R-hat | Min bulk ESS |
|---|---:|---:|---|---|---:|---:|
| Simulation 1 | 1.867 | 1.880 (0.023) | .674, .492, .412, .321, .176 | .672, .498, .389, .348, .164 | 1.0096 | 431 |
| Simulation 2 | 2.245 | 2.180 (0.023) | .817, .402, .321, .226, .130 | .822, .388, .326, .215, .144 | 1.0041 | 1,837 |

The lower spectrum ESS in Simulation 1 is concentrated in its third and
fourth values, the closest pair. This is the expected weak geometry near a
spectral tie, not a rotation of the whole factor space.

## Fit and recovery

| Data | Model | Observed RMSE | True mean RMSE | True mean corr. | True mean 90% coverage | True interaction RMSE | True interaction corr. | Interaction 90% coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Simulation 1 | Modernized source | 0.8191 | 0.5576 | 0.9900 | 0.9015 | 0.5015 | 0.9633 | 0.9035 |
| Simulation 1 | Semi-orthogonal | 0.8218 | 0.5588 | 0.9900 | 0.9010 | 0.5032 | 0.9630 | 0.9035 |
| Simulation 2 | Modernized source | 0.8347 | 0.5487 | 0.9927 | 0.9190 | 0.4935 | 0.9760 | 0.9200 |
| Simulation 2 | Semi-orthogonal | 0.8390 | 0.5486 | 0.9927 | 0.9155 | 0.4933 | 0.9761 | 0.9155 |

Direct posterior comparison:

| Data | Quantity | Posterior-mean RMSE between models | Posterior-mean correlation |
|---|---|---:|---:|
| Simulation 1 | Total surface | 0.0222 | 0.999984 |
| Simulation 1 | Interaction | 0.0217 | 0.999928 |
| Simulation 2 | Total surface | 0.0237 | 0.999986 |
| Simulation 2 | Interaction | 0.0232 | 0.999940 |

The mean absolute surface difference divided cellwise by
`sqrt(sd_semi^2 + sd_source^2)` averages `0.021--0.023`. This is a quadrature
standardization, not the conventional pooled SD.
Nominal PSIS-LOO differences are +0.6 and +2.7 ELPD for the semi-orthogonal
model, but they are not usable for selection: 88--127 observations per fit
have Pareto `k > 0.7`, with maxima above 0.97. Use exact refits or K-fold
validation for a predictive comparison.

## Raw coordinates and parameter count

| Data | Model | Raw factor/reflector coordinates | Max R-hat | R-hat > 1.01 | Min bulk ESS |
|---|---:|---:|---:|---:|---:|
| Simulation 1 | Source factors | 600 | 1.036 | 403 | 159 |
| Simulation 1 | Semi reflector seeds | 570 | 2.444 | 568 | 4.9 |
| Simulation 2 | Source factors | 600 | 1.082 | 560 | 55 |
| Simulation 2 | Semi reflector seeds | 570 | 2.664 | 473 | 4.7 |

The reflector seed diagnostics are intentionally not a fit criterion. Their
directions determine the frames, their radii are likelihood-free, and paired
frame signs remain multimodal. The implementation excludes `zfm_frame*` and
`zfm_spectrum*` by default while retaining the centered frames, `sifm`, and
all quantities required for prediction. `save_pars(all = TRUE)` remains
available for bridge sampling and audits.

For this 100 by 20, rank-five model:

| Model | Total free coordinates | Interaction coordinates | Continuous excess over identifiable model |
|---|---:|---:|---:|
| Source | 720 | 600 | 36 |
| Native S2Z draft | 709 | 590 | 25 |
| Semi-orthogonal | 694 | 575 | 10 |

The semi interaction has the intrinsic 565 dimensions plus ten prior-only
reflector radii, one per reflector per field. It removes 26 continuous excess
coordinates relative to the source model.

## Geometry checks

- interaction row/column sums are below `5.2e-13` after canonicalization in
  the first 100 draws of each fit;
- the sixth singular value is below `6.5e-14` in those same checked draws;
- saved frame column sums and orthogonality errors are below `6.4e-8` and
  `2.3e-8`, respectively;
- the saved spectrum is always ordered and its squared-norm error is below
  `1.8e-8`;
- reconstructed full-table RMS differs from `sdfm` by less than `3.0e-8`.

The remaining geometry bullets use all 4,000 post-warmup draws per fit. Their
`1e-8` scale reflects eight-digit CmdStan CSV output. Direct transform tests
before serialization are accurate to roughly `1e-15`; the tracked transform
check covers rectangular frames, odd/even special-orthogonal frames, the
zero-coordinate `L = 2, K = 1` edge, and a Haar moment check.

## Limitations and recommendation

1. The comparison covers the repository's basic Gaussian model, not its more
   elaborate `stan_fm_2.stan` or `stan_fm_3.stan` variants.
2. Both data sets are complete grids. Sparse, unbalanced, and held-out-pair
   benchmarks remain necessary; the source model itself hard-codes `N * J`
   observations.
3. There is one four-chain run per data set, so timing differences this small
   should be treated as parity.
4. The independently calibrated canonical prior is only approximate: it does
   not reproduce the source model's joint dependence or exact product-Wishart
   spectrum law.
5. Raw frames still require sign alignment for visualization. Their entries
   should not be interpreted as invariant effects.
6. The reflector routine retains the upstream transform's TODO around extreme
   underflow rescaling. Standard-normal seeds did not trigger a retained-draw
   failure here, but extreme-value stress tests are still warranted.
7. Equal-size full-rank interactions have two physical determinant components.
   The benchmark rank is below both centered dimensions, so its diagnostics do
   not exercise that disconnected edge geometry.

The implementation passes the fit, invariant-diagnostic, geometry, and
performance checks on these examples. It is preferable to the native-matrix
draft: it gives an explicit RMS scale and ordered spectrum, removes continuous
rotations, uses fewer coordinates, and restores timing parity with the
unconstrained source model. Claims should remain limited to identifiable
quantities; raw reflector seeds are not suitable convergence diagnostics.
