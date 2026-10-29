// cpp11 glue over the native libcusna C ABI (track A, milestone A5).
//
// The CPU-only build (CRAN default) compiles this file together with the
// vendored pure-C++ host logic (libcusna_saom_cpu.cpp) -- no CUDA needed. The
// exposed functions are the validated convention CARRIERS (host statistics):
// they reproduce the RSiena Appendix B targets to machine precision (proven in
// native/test). The full period simulator and the GPU path are compiled in when
// configure detects a CUDA toolkit (see src/Makevars.in / configure).
//
// LIBCUSNA_BUILD is defined by Makevars so CUSNA_API resolves to the local
// (exported) linkage for the statically-linked ABI functions.
#include <cpp11.hpp>
#include <vector>
#include "libcusna.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace cpp11;

[[cpp11::register]]
int cusna_abi_version_() {
  return cusna_abi_version();
}

[[cpp11::register]]
int cusna_has_cuda_() {
#ifdef CUSNA_WITH_CUDA
  return cusna_has_cuda();
#else
  return 0;  // CPU-only build: no GPU path compiled in
#endif
}

// Explicit thread-count control for the OpenMP-parallel CPU backend. An
// imperative omp_set_num_threads() call takes effect immediately and doesn't
// depend on when/whether OMP_NUM_THREADS was read by the runtime (unlike
// Sys.setenv(), whose effect on an already-initialised OpenMP runtime is
// platform-dependent). Returns the previous max thread count (0 if OpenMP is
// not compiled in, e.g. a default macOS toolchain without libomp).
[[cpp11::register]]
int cusna_set_num_threads_(int n) {
#ifdef _OPENMP
  int prev = omp_get_max_threads();
  omp_set_num_threads(n < 1 ? 1 : n);
  return prev;
#else
  return 0;
#endif
}

// Host evaluation statistics of one network (cusna_network_stats). `a` is the
// dense n*n statistic network (0/1, row-major); the descriptors are the float64
// arrays network_stats_ex sums. Returns the k effect statistics.
[[cpp11::register]]
doubles cusna_network_stats_(integers a, int n, int k, integers type_ids,
                             doubles p1, doubles p2, doubles vmat, integers vmiss,
                             integers dyn_flags, doubles vb, integers od0,
                             integers id0, doubles wtab, int symflag) {
  std::vector<uint8_t> ad(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) ad[i] = (uint8_t)(a[i] != 0);
  std::vector<uint8_t> vm(vmiss.size());
  for (R_xlen_t i = 0; i < vmiss.size(); ++i) vm[i] = (uint8_t)(vmiss[i] != 0);
  writable::doubles out(k);
  cusna_status st = cusna_network_stats(
      ad.data(), n, k,
      INTEGER(type_ids), REAL(p1), REAL(p2), REAL(vmat), vm.data(),
      INTEGER(dyn_flags), REAL(vb), INTEGER(od0), INTEGER(id0), REAL(wtab),
      symflag, REAL(out));
  if (st != CUSNA_OK) stop("cusna_network_stats failed (status %d)", (int)st);
  return out;
}

// Behavior evaluation statistics over a batch of end behavior vectors
// (cusna_behavior_stats). z is M*n, xnet the n*n predictor network. Returns
// M*kb (row-major, chain-major).
[[cpp11::register]]
doubles cusna_behavior_stats_(integers z, int n, int M, int kb, integers bty,
                              integers xnet, integers bkdist, double zbar,
                              double brange, double bsimmean) {
  std::vector<uint8_t> xn(xnet.size());
  for (R_xlen_t i = 0; i < xnet.size(); ++i) xn[i] = (uint8_t)(xnet[i] != 0);
  writable::doubles out((R_xlen_t)M * kb);
  cusna_status st = cusna_behavior_stats(
      INTEGER(z), n, M, kb, INTEGER(bty), xn.data(), INTEGER(bkdist),
      zbar, brange, bsimmean, REAL(out));
  if (st != CUSNA_OK) stop("cusna_behavior_stats failed (status %d)", (int)st);
  return out;
}

// sienaGOF auxiliary distributions over a batch of end networks
// (cusna_gof_distribution). a is M*n*n; kind 0 indegree / 1 outdegree
// (out M*nlev) / 2 geodesic (out M*(nlev+1)).
[[cpp11::register]]
doubles cusna_gof_distribution_(integers a, int n, int M, int kind,
                                integers levls, int nlev) {
  std::vector<uint8_t> ad(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) ad[i] = (uint8_t)(a[i] != 0);
  R_xlen_t cols = (kind == 2) ? (nlev + 1) : nlev;
  writable::doubles out((R_xlen_t)M * cols);
  cusna_status st = cusna_gof_distribution(
      ad.data(), n, M, kind, INTEGER(levls), nlev, REAL(out));
  if (st != CUSNA_OK) stop("cusna_gof_distribution failed (status %d)", (int)st);
  return out;
}

// ---------------------------------------------------------------------------
//  ERGM (track A-ERGM, milestone E6). All three are host functions from the
//  vendored libcusna_ergm_cpu.cpp, so they work in the CPU-only CRAN build.
//  Framing (M-3): cusna's ERGM is a sampler + MPLE demo, NOT a full MCMC-MLE.
// ---------------------------------------------------------------------------

// ERGM host statistics of one network (cusna_ergm_stats). `a` dense n*n (0/1,
// row-major); terms as type_ids (k) + parm (k); attr length n, empty for none.
[[cpp11::register]]
doubles cusna_ergm_stats_(integers a, int n, int k, int directed,
                          integers type_ids, doubles parm, doubles attr) {
  std::vector<uint8_t> ad(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) ad[i] = (uint8_t)(a[i] != 0);
  std::vector<float> pf(k), af(n, 0.0f);
  for (int e = 0; e < k; ++e) pf[e] = (float)parm[e];
  bool has_attr = attr.size() == (R_xlen_t)n;
  if (has_attr) for (int i = 0; i < n; ++i) af[i] = (float)attr[i];
  writable::doubles out(k);
  cusna_status st = cusna_ergm_stats(
      ad.data(), n, k, directed, INTEGER(type_ids), pf.data(),
      has_attr ? af.data() : nullptr, REAL(out));
  if (st != CUSNA_OK) stop("cusna_ergm_stats failed (status %d)", (int)st);
  return out;
}

// ERGM TNT simulate via the CPU mirror (the native package is CPU-only by
// default). Returns the C*samplesize*k sampled statistics, row-major
// (chain, sample, term). s0 is computed internally (passed NULL). When
// return_nets != 0 the C final networks (a_ends, [C*n*n] uint8 row-major) are
// captured too; the return value is then a list(stats = <as above>, nets =
// <integers C*n*n row-major>). With return_nets == 0 the bare stats vector is
// returned, preserving the historical (scalar) call.
[[cpp11::register]]
SEXP cusna_ergm_simulate_cpu_(integers a, int n, int k, int directed,
                              doubles theta, integers type_ids, doubles parm,
                              doubles attr, int C, double seed, int burnin,
                              int interval, int samplesize, int return_nets) {
  std::vector<uint8_t> ad(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) ad[i] = (uint8_t)(a[i] != 0);
  std::vector<float> thf(k), pf(k), af(n, 0.0f);
  for (int e = 0; e < k; ++e) { thf[e] = (float)theta[e]; pf[e] = (float)parm[e]; }
  bool has_attr = attr.size() == (R_xlen_t)n;
  if (has_attr) for (int i = 0; i < n; ++i) af[i] = (float)attr[i];

  cusna_ergm_problem p{};
  p.struct_version = cusna_abi_version();
  p.n = n; p.W = (n + 31) / 32; p.k = k; p.C = C;
  p.directed = directed; p.threads = 0;
  p.x0 = ad.data(); p.theta = thf.data(); p.type_ids = INTEGER(type_ids);
  p.parm = pf.data(); p.attr = has_attr ? af.data() : nullptr; p.s0 = nullptr;
  p.seed = (uint64_t)seed; p.burnin = burnin; p.interval = interval;
  p.samplesize = samplesize;

  writable::doubles out((R_xlen_t)C * samplesize * k);
  std::vector<uint8_t> aends;
  cusna_ergm_result r{};
  r.struct_version = cusna_abi_version();
  r.stats = REAL(out);
  if (return_nets) {
    aends.assign((size_t)C * n * n, 0);
    r.a_ends = aends.data();
  }
  cusna_status st = cusna_ergm_tnt_sample_cpu(&p, &r);
  if (st != CUSNA_OK) stop("cusna_ergm_tnt_sample_cpu failed (status %d)", (int)st);
  if (!return_nets) return (SEXP)out;
  using namespace cpp11::literals;
  writable::integers nets((R_xlen_t)C * n * n);
  for (R_xlen_t i = 0; i < (R_xlen_t)C * n * n; ++i) nets[i] = (int)aends[i];
  return (SEXP)writable::list({"stats"_nm = out, "nets"_nm = nets});
}

// ERGM MPLE change-statistic design matrix (cusna_ergm_mple_stats). Returns a
// list(X = nrow*ncol row-major doubles, y = nrow ints, nrow, ncol).
[[cpp11::register]]
list cusna_ergm_mple_(integers a, int n) {
  using namespace cpp11::literals;
  std::vector<uint8_t> ad(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) ad[i] = (uint8_t)(a[i] != 0);
  std::vector<float> X((size_t)n * n * 3);
  std::vector<uint8_t> y((size_t)n * n);
  int nrow = 0, ncol = 0;
  cusna_status st = cusna_ergm_mple_stats(ad.data(), n, X.data(), y.data(),
                                          &nrow, &ncol);
  if (st != CUSNA_OK) stop("cusna_ergm_mple_stats failed (status %d)", (int)st);
  writable::doubles Xout((R_xlen_t)nrow * ncol);
  for (int row = 0; row < nrow; ++row)
    for (int col = 0; col < ncol; ++col)
      Xout[(R_xlen_t)row * ncol + col] = (double)X[(size_t)row * 3 + col];
  writable::integers yout(nrow);
  for (int row = 0; row < nrow; ++row) yout[row] = (int)y[row];
  return writable::list({"X"_nm = Xout, "y"_nm = yout,
                         "nrow"_nm = nrow, "ncol"_nm = ncol});
}

// Temporal ERGM (TERGM) simulation via the CPU mirror: like cusna_ergm_simulate_cpu_
// but with a fixed lagged network for the temporal terms memory (7) / delrecip (8).
// Returns the C*samplesize*k sampled statistics, row-major (chain, sample, term).
// When return_nets != 0 the C final networks (a_ends, [C*n*n] uint8 row-major) are
// captured and the return value becomes list(stats = <as above>, nets = <integers
// C*n*n row-major>); with return_nets == 0 the bare stats vector is returned.
[[cpp11::register]]
SEXP cusna_ergm_tsim_cpu_(integers a, integers lag, int n, int k, int directed,
                          doubles theta, integers type_ids, doubles parm,
                          doubles attr, doubles s0, int C, double seed,
                          int burnin, int interval, int samplesize,
                          int return_nets) {
  std::vector<uint8_t> ad(a.size()), lagd(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) {
    ad[i] = (uint8_t)(a[i] != 0);
    lagd[i] = (uint8_t)(lag[i] != 0);
  }
  std::vector<float> thf(k), pf(k), af(n, 0.0f);
  for (int e = 0; e < k; ++e) { thf[e] = (float)theta[e]; pf[e] = (float)parm[e]; }
  bool has_attr = attr.size() == (R_xlen_t)n;
  if (has_attr) for (int i = 0; i < n; ++i) af[i] = (float)attr[i];

  cusna_ergm_problem p{};
  p.struct_version = cusna_abi_version();
  p.n = n; p.W = (n + 31) / 32; p.k = k; p.C = C;
  p.directed = directed; p.threads = 0;
  p.x0 = ad.data(); p.theta = thf.data(); p.type_ids = INTEGER(type_ids);
  p.parm = pf.data(); p.attr = has_attr ? af.data() : nullptr;
  p.s0 = REAL(s0);                          /* R computes s0 (incl. temporal) */
  p.seed = (uint64_t)seed; p.burnin = burnin; p.interval = interval;
  p.samplesize = samplesize;
  p.lag = lagd.data();

  writable::doubles out((R_xlen_t)C * samplesize * k);
  std::vector<uint8_t> aends;
  cusna_ergm_result r{};
  r.struct_version = cusna_abi_version();
  r.stats = REAL(out);
  if (return_nets) {
    aends.assign((size_t)C * n * n, 0);
    r.a_ends = aends.data();
  }
  cusna_status st = cusna_ergm_tnt_sample_cpu(&p, &r);
  if (st != CUSNA_OK) stop("cusna_ergm_tnt_sample_cpu failed (status %d)", (int)st);
  if (!return_nets) return (SEXP)out;
  using namespace cpp11::literals;
  writable::integers nets((R_xlen_t)C * n * n);
  for (R_xlen_t i = 0; i < (R_xlen_t)C * n * n; ++i) nets[i] = (int)aends[i];
  return (SEXP)writable::list({"stats"_nm = out, "nets"_nm = nets});
}

// ---------------------------------------------------------------------------
//  ALAAM (autologistic actor attribute model; track A-ERGM extension X1.5).
//  Host functions from the vendored libcusna_ergm_cpu.cpp -- attribute-DV.
// ---------------------------------------------------------------------------

// Per-actor ALAAM contagion change statistics (cusna_alaam_stats). `a` dense n*n
// (0/1 row-major), `y` a length-n 0/1 attribute. Returns list(con_out, con_in).
[[cpp11::register]]
list cusna_alaam_stats_(integers a, integers y, int n) {
  using namespace cpp11::literals;
  std::vector<uint8_t> ad(a.size()), yd(y.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) ad[i] = (uint8_t)(a[i] != 0);
  for (R_xlen_t i = 0; i < y.size(); ++i) yd[i] = (uint8_t)(y[i] != 0);
  writable::doubles co(n), ci(n);
  cusna_status st = cusna_alaam_stats(ad.data(), yd.data(), n, REAL(co), REAL(ci));
  if (st != CUSNA_OK) stop("cusna_alaam_stats failed (status %d)", (int)st);
  return writable::list({"con_out"_nm = co, "con_in"_nm = ci});
}

// ALAAM Gibbs sampler over the binary attribute (cusna_alaam_gibbs_cpu). X is the
// n*p row-major covariate matrix (empty for none); beta_x length p. Returns
// list(y_ends = C*n ints row-major, stats = C*2 doubles [density, contagion]).
[[cpp11::register]]
list cusna_alaam_gibbs_cpu_(integers a, int n, int directed, integers y0,
                            doubles X, int p, double beta0, double beta_c,
                            doubles beta_x, int C, double seed, int sweeps) {
  using namespace cpp11::literals;
  std::vector<uint8_t> ad(a.size()), y0d(y0.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) ad[i] = (uint8_t)(a[i] != 0);
  for (R_xlen_t i = 0; i < y0.size(); ++i) y0d[i] = (uint8_t)(y0[i] != 0);
  std::vector<float> Xf((size_t)(p > 0 ? (R_xlen_t)n * p : 0)), bxf(p > 0 ? p : 0);
  for (R_xlen_t i = 0; i < X.size(); ++i) Xf[i] = (float)X[i];
  for (int q = 0; q < p; ++q) bxf[q] = (float)beta_x[q];
  std::vector<uint8_t> yb((size_t)C * n);
  writable::doubles stats((R_xlen_t)C * 2);
  cusna_status st = cusna_alaam_gibbs_cpu(
      ad.data(), n, directed, y0d.data(), p > 0 ? Xf.data() : nullptr, p,
      beta0, beta_c, p > 0 ? bxf.data() : nullptr, C, (uint64_t)seed, sweeps,
      yb.data(), REAL(stats));
  if (st != CUSNA_OK) stop("cusna_alaam_gibbs_cpu failed (status %d)", (int)st);
  writable::integers yends((R_xlen_t)C * n);
  for (R_xlen_t i = 0; i < (R_xlen_t)C * n; ++i) yends[i] = (int)yb[i];
  return writable::list({"y_ends"_nm = yends, "stats"_nm = stats});
}

// ---------------------------------------------------------------------------
//  STERGM (separable temporal ERGM; track A-ERGM extension X2). The formation
//  and dissolution sub-processes are constrained ERGMs: the TNT sampler toggles
//  only the permitted (free) dyads. Reuses the CPU sampler with a permitted mask.
// ---------------------------------------------------------------------------

// ERGM TNT simulate with a toggle-permitted (free) dyad mask (cusna_ergm_problem
// ->permitted): like cusna_ergm_simulate_cpu_ but the chain toggles only the free
// dyads (`permitted` is a dense n*n 0/1 row-major mask). s0 is computed internally
// (standard terms). Returns the C*samplesize*k sampled statistics, row-major. When
// return_nets != 0 the C final (masked) networks (a_ends, [C*n*n] uint8 row-major)
// are captured too and the return value becomes list(stats = <as above>, nets =
// <integers C*n*n row-major>); with return_nets == 0 the bare stats vector is
// returned, preserving the historical (scalar) call.
[[cpp11::register]]
SEXP cusna_ergm_sim_masked_cpu_(integers a, integers permitted, int n, int k,
                                int directed, doubles theta, integers type_ids,
                                doubles parm, doubles attr, int C, double seed,
                                int burnin, int interval, int samplesize,
                                int return_nets) {
  std::vector<uint8_t> ad(a.size()), pm(a.size());
  for (R_xlen_t i = 0; i < a.size(); ++i) {
    ad[i] = (uint8_t)(a[i] != 0);
    pm[i] = (uint8_t)(permitted[i] != 0);
  }
  std::vector<float> thf(k), pf(k), af(n, 0.0f);
  for (int e = 0; e < k; ++e) { thf[e] = (float)theta[e]; pf[e] = (float)parm[e]; }
  bool has_attr = attr.size() == (R_xlen_t)n;
  if (has_attr) for (int i = 0; i < n; ++i) af[i] = (float)attr[i];

  cusna_ergm_problem p{};
  p.struct_version = cusna_abi_version();
  p.n = n; p.W = (n + 31) / 32; p.k = k; p.C = C;
  p.directed = directed; p.threads = 0;
  p.x0 = ad.data(); p.theta = thf.data(); p.type_ids = INTEGER(type_ids);
  p.parm = pf.data(); p.attr = has_attr ? af.data() : nullptr; p.s0 = nullptr;
  p.seed = (uint64_t)seed; p.burnin = burnin; p.interval = interval;
  p.samplesize = samplesize;
  p.permitted = pm.data();                  /* free-dyad mask (STERGM) */

  writable::doubles out((R_xlen_t)C * samplesize * k);
  std::vector<uint8_t> aends;
  cusna_ergm_result r{};
  r.struct_version = cusna_abi_version();
  r.stats = REAL(out);
  if (return_nets) {
    aends.assign((size_t)C * n * n, 0);
    r.a_ends = aends.data();
  }
  cusna_status st = cusna_ergm_tnt_sample_cpu(&p, &r);
  if (st != CUSNA_OK) stop("cusna_ergm_tnt_sample_cpu (masked) failed (status %d)", (int)st);
  if (!return_nets) return (SEXP)out;
  using namespace cpp11::literals;
  writable::integers nets((R_xlen_t)C * n * n);
  for (R_xlen_t i = 0; i < (R_xlen_t)C * n * n; ++i) nets[i] = (int)aends[i];
  return (SEXP)writable::list({"stats"_nm = out, "nets"_nm = nets});
}
