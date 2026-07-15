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
