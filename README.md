# cusna (native) — GPU-accelerated SAOM/RSiena, ERGM, and friends for R

`cusna` is a self-contained **native** engine — a C ABI over CUDA kernels and
C++ host logic (`libcusna`) — callable from R **without a Python runtime**. It is
the native counterpart of the reticulate-based `cusna` R wrapper, modelled on the
R [`torch`](https://torch.mlverse.org/) package:

- the **CRAN build is CPU-only from source** (pure C++, needs only a C++17
  compiler);
- the **GPU path is compiled from source** when a CUDA 12.x toolkit (`nvcc`) is
  detected at configure time.

`cusna_has_cuda()` reports which build you have.

## What the package provides

| Family | Functions | Validated against |
|---|---|---|
| **SAOM** (RSiena) | `saom_data()`, `cusna_effect()`, `mom_estimate()`, `cusna_fit` methods; behavior co-evolution, composition change, `mom_estimate_multinet()`, `cusna_fran()` | data/masks/targets bit-identical to the reference; estimates within simulation SE; RSiena targets to machine zero |
| **ERGM** | `ergm_simulate()` (TNT sampler), `ergm_stats()`, `ergm_mple()`, `ergm_mcmle()` | sampler ≡ `ergm::simulate`; MLE matches `ergm::ergm()` |
| **Temporal ERGM** | `tergm_mple()` (+ block bootstrap), `tergm_simulate()`, `stergm_cmle()` | matches `btergm` to machine precision; `tergm` CMLE within SE |
| **ALAAM** | `alaam_mple()`, `alaam_mcmle()`, `alaam_simulate()` | MPLE ≡ `glm`; MLE recovers observed moments |
| Low-level | `cusna_network_stats()`, `cusna_behavior_stats()`, `cusna_gof_distribution()` | RSiena Appendix B conventions, machine zero |

The underlying C ABI is bit-for-bit validated in `native/test`
(see `native/VALIDATION.md`).

## Reviewer quickstart (CPU-only, no GPU needed)

```r
# from a checkout of the monorepo (configure vendors ../../native sources):
install.packages("cpp11")            # build-time only
# then, with a C++17 toolchain (Rtools on Windows):
#   R CMD INSTALL Rpkg-native
library(cusna)
cusna_has_cuda()                     # FALSE on the CPU-only build

# a two-wave panel and a Method-of-Moments SAOM fit, all native:
set.seed(7)
w1 <- matrix(as.integer(runif(400) < 0.12), 20, 20); diag(w1) <- 0L
w2 <- w1; flip <- sample(400, 40); w2[flip] <- 1L - w2[flip]; diag(w2) <- 0L
fit <- mom_estimate(saom_data(list(w1, w2)),
                    effects = list(cusna_effect("density"), cusna_effect("recip")))
summary(fit)

# an ERGM maximum-likelihood fit on the same data:
ergm_mcmle(w1, list(ergm_term("edges"), ergm_term("mutual")), directed = TRUE)
```

See `vignette("cusna")` for the full tour (covariates, co-evolution,
multi-network models, TERGM/STERGM/ALAAM) and `vignette("siena07-backend")`
for driving RSiena's `siena07()` on the native simulator.

## Building

See [BUILD.md](https://github.com/artemmaltsev74-techcom/cusna/blob/main/Rpkg-native/BUILD.md).
In short: the CPU-only build needs a C++17 compiler
(**Rtools** on Windows); the GPU build additionally needs a CUDA 12.x toolkit.
The native sources are vendored from `../../native/` by `configure`.

## License

MIT (see `LICENSE`). We compare outputs against RSiena but do not link its code.
