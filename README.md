# cusna (native) — simulation and estimation of SAOM/RSiena, ERGM, and friends for R

`cusna` is a self-contained **native** engine — a small C interface over C++
host logic (`libcusna`) — callable from R **without a Python runtime**:

- the R package builds the engine's **CPU backend from source** (C++17,
  multi-threaded with OpenMP; `cusna_set_threads()` sets the thread count);
- a **CUDA backend** implements the same C interface and is built separately
  from the `native/` sources (it is driven by the companion Python package);
  `cusna_has_cuda()` reports `FALSE` for the R package.

## What the package provides

| Family | Functions | Validated against |
|---|---|---|
| **SAOM** (RSiena) | `saom_data()`, `cusna_effect()`, `mom_estimate()` (with `return_nets`), `cusna_fit` methods, `saom_simulate()`; behavior co-evolution, composition change, `mom_estimate_multinet()`, `cusna_fran()` for `siena07()` | data/masks/targets bit-identical to the reference; estimates within simulation SE; RSiena targets to machine zero |
| **ERGM** | `ergm_simulate()` (TNT sampler), `ergm_stats()`, `ergm_mple()`, `ergm_mcmle()` | sampler ≡ `ergm::simulate`; MLE matches `ergm::ergm()` |
| **Temporal ERGM** | `tergm_mple()` (+ block bootstrap), `tergm_simulate()`, `stergm_cmle()`, `stergm_simulate()` | matches `btergm` to machine precision; `tergm` CMLE within SE |
| **ALAAM** | `alaam_mple()`, `alaam_mcmle()`, `alaam_simulate()` | MPLE ≡ `glm`; MLE recovers observed moments |
| Low-level | `cusna_network_stats()`, `cusna_behavior_stats()`, `cusna_gof_distribution()` | RSiena Appendix B conventions, machine zero |

The underlying C interface is validated bit-for-bit in `native/test`
(see `native/VALIDATION.md`).

## Quickstart

```r
install.packages("cusna")
library(cusna)

# the bundled s50 panel: three waves of friendship among 50 girls
dat <- saom_data(s50$waves)
fit <- mom_estimate(dat,
                    effects = list(cusna_effect("density"), cusna_effect("recip"),
                                   cusna_effect("transTrip")),
                    conditional = TRUE, return_nets = TRUE)
summary(fit)

# goodness of fit on the phase-3 networks: simulated indegree distribution
cusna_gof_distribution(fit$sims[[2]], kind = 0L, levls = 0:6)

# an ERGM maximum-likelihood fit on the first wave
ergm_mcmle(s50$waves[[1]], list(ergm_term("edges"), ergm_term("mutual")),
           directed = TRUE)
```

See `vignette("cusna")` for the full tour (covariates, co-evolution,
multi-network models, TERGM/STERGM/ALAAM) and `vignette("siena07-backend")`
for driving RSiena's `siena07()` on the native simulator.

## Building from the repository

See [BUILD.md](https://github.com/artemmaltsev74-techcom/cusna/blob/main/Rpkg-native/BUILD.md).
In short: the package needs a C++17 compiler (**Rtools** on Windows); its
`configure` script vendors the native sources from `../native/`.

## License

MIT (see `LICENSE`). We compare outputs against RSiena but do not link its code.
