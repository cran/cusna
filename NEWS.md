# cusna 0.1.1

## Bug fixes

* `cusna_fran()` now honours RSiena's `simstats0c` contract: `fra` and `sc`
  are periods-by-parameters matrices (one row per period). In 0.1.0 `fra`
  was the period-summed vector and `sc` held the last period only, so
  `siena07()` stopped at its first call (in `colSums(zz$fra)`). The closure
  also checks that the parameter vector passed by `siena07()` matches
  `effect_names` and `conditional`, draws its seeds from R's random number
  stream (so `sienaAlgorithmCreate(seed = )` makes a run reproducible, and the
  global seed is no longer reset on every call), and counts its calls in
  `attr(fran, "counter")$i`.
* `mom_estimate()` can now return the simulated networks: with
  `return_nets = TRUE` the fit keeps the phase-3 end networks in `fit$sims`
  (one chains-by-n-by-n integer array per period, the layout
  `cusna_gof_distribution()` takes) and, for co-evolution, the end behavior in
  `fit$sims_behavior`. In 0.1.0 the simulator was always called without its
  end networks.

## Packaging and documentation

* New dataset `s50` (three waves of friendship among 50 girls, with alcohol
  use; the RSiena teaching panel), used by the examples and vignettes.
* All examples run: no `\dontrun{}` remains, and the SAOM examples use `s50`.
  The `siena07()` vignette runs its example when RSiena is installed.
* Under `R CMD check --as-cran` (`_R_CHECK_LIMIT_CORES_`) the OpenMP engine
  uses at most two threads.
* New `stergm_simulate()`: draws STERGM transitions on the constrained spaces
  of `stergm_cmle()`.
* The description, `configure` messages and `cusna_has_cuda()` documentation
  now state that the R package builds the CPU backend of the engine; the
  CUDA backend of the same C interface is built separately.

# cusna 0.1.0

First release: a self-contained native engine (C ABI over CUDA kernels and
C++ host logic) for the simulation and estimation of network models, callable
from R without a Python runtime. The CRAN build is CPU-only from source; the
GPU path compiles when a CUDA toolkit is detected at configure time.

* **SAOM** (the RSiena model family). `saom_data()`; effect constructors
  (`cusna_effect()`, `cusna_beh_effect()`, `cusna_rate_effect()`,
  `cusna_interaction()`); the full Method-of-Moments estimator
  `mom_estimate()` / `mom_control()` returning a `cusna_fit` object with
  `summary()`, `coef()`, `vcov()`, and `as.data.frame()` methods; behavior
  co-evolution, composition change, conditional and unconditional estimation;
  multi-network co-evolution (`saom_multinet_data()`,
  `mom_estimate_multinet()`); and a `siena07()` simulation backend
  (`cusna_fran()`). Data preparation, effect preprocessing, and moment
  targets are validated bit-for-bit against the reference implementation;
  estimates agree within simulation standard errors.
* **ERGM.** Sufficient statistics (`ergm_stats()`), a TNT sampler
  (`ergm_simulate()`), pseudo-likelihood (`ergm_mple()`), and MCMC maximum
  likelihood (`ergm_mcmle()`), matching `ergm::ergm()` on benchmark models.
* **Temporal ERGM.** `tergm_mple()` (pooled MPLE with block bootstrap,
  matching `btergm`), temporal simulation (`tergm_simulate()`), and the
  separable `stergm_cmle()` (formation/persistence, matching `tergm` CMLE).
* **ALAAM.** `alaam_mple()` (exactly reproducing the corresponding `glm`),
  `alaam_mcmle()`, and a Gibbs simulator (`alaam_simulate()`).
* Low-level host statistics (`cusna_network_stats()`,
  `cusna_behavior_stats()`, `cusna_gof_distribution()`) reproduce RSiena
  targets on public datasets to machine precision.
