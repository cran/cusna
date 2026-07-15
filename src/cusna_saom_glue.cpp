// cpp11 glue over the native SAOM period simulator + typed/cross statistic
// carriers (track A-RAPI, milestones R4 + R3b). The Robbins-Monro loop in
// R/mom.R calls cusna_simulate_period_cpu_ once per period/batch; the targets
// assembly (R/targets.R) uses the typed/cross wrappers. All of these drive the
// EXISTING C ABI (native/libcusna.h) -- no ABI extension. The single-network
// simulator marshals cusna_saom_problem from R arrays (float64 host descriptors,
// A4.3), mirroring native/test/test_native_estimate.py. Masks default to the
// identity (NULL) case (s50 / complete data); the estimator bar (baseline /
// covariate / coevolution / multinet) needs no explicit masks.
#include <cpp11.hpp>
#include <vector>
#include <cstring>
#include "libcusna.h"

using namespace cpp11;

namespace {

// Pack a dense n*n row-major uint8 adjacency into out/in bit rows (n*W words).
// out[i] bit j = a[i,j] (out-ties of i); in[i] bit j = a[j,i] (in-ties of i),
// matching the unpack() convention in libcusna_saom_cpu.cpp.
void pack_out_in(const uint8_t *dense, int n, int W,
                 std::vector<uint32_t> &out, std::vector<uint32_t> &in) {
  out.assign((size_t)n * W, 0u);
  in.assign((size_t)n * W, 0u);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (dense[(size_t)i * n + j]) {
        out[(size_t)i * W + (j >> 5)] |= (1u << (j & 31));
        in[(size_t)j * W + (i >> 5)] |= (1u << (i & 31));
      }
}

// dense uint8 from an R integer vector (row-major, values != 0 -> 1)
std::vector<uint8_t> dense_u8(const integers &m) {
  std::vector<uint8_t> d(m.size());
  for (R_xlen_t i = 0; i < m.size(); ++i) d[i] = (uint8_t)(m[i] != 0);
  return d;
}

// Pack one dense n*n mask into out-oriented bit rows (out[i] bit j = m[i,j]).
void pack_single(const uint8_t *dense, int n, int W, uint32_t *out) {
  std::memset(out, 0, (size_t)n * W * sizeof(uint32_t));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (dense[(size_t)i * n + j]) out[(size_t)i * W + (j >> 5)] |= (1u << (j & 31));
}

}  // namespace

// One SAOM period, M chains, on the CPU backend. `x0` is the dense n*n start
// network (row-major). `theta` (length k) and `lam` (scalar) are broadcast to
// the M chains. Descriptors are the float64 host arrays (vmat/wtab/vb/p1/p2);
// vmiss is k*n. With kb>0 the behavior sub-block is marshalled and the behavior
// scores/end values are returned (behavior eval statistics are computed in R via
// cusna_behavior_stats). Returns a named list of the result buffers.
[[cpp11::register]]
list cusna_simulate_period_cpu_(
    int n, int k, int M, double seed,
    integers x0, double lam, doubles theta,
    integers type_ids, integers etypes, integers dyn_flags, integers dyc_idx,
    doubles vmat, doubles wtab, doubles vb, doubles p1, doubles p2,
    integers vmiss, int maxdegree, int conditional, int target_dist, int symflag,
    int kb, double lamb, doubles thetab, integers bty, integers z0,
    integers bkdist, int zmin, int zmax, double zbar, double brange,
    double bsimmean, int want_ends) {
  const int W = (n + 31) / 32;

  // start network -> packed bit rows
  std::vector<uint8_t> x0d = dense_u8(x0);
  std::vector<uint32_t> x0_out, x0_in;
  pack_out_in(x0d.data(), n, W, x0_out, x0_in);

  // per-chain rate + eval params (float32, broadcast over M)
  std::vector<float> lamv((size_t)M, (float)lam);
  std::vector<float> thv((size_t)M * k);
  for (int c = 0; c < M; ++c)
    for (int e = 0; e < k; ++e) thv[(size_t)c * k + e] = (float)theta[e];

  // float64 host descriptors (A4.3) + int descriptors
  std::vector<double> vmat64(vmat.begin(), vmat.end());
  std::vector<double> wtab64(wtab.begin(), wtab.end());
  std::vector<double> vb64(vb.begin(), vb.end());
  std::vector<double> p1_64(p1.begin(), p1.end());
  std::vector<double> p2_64(p2.begin(), p2.end());
  std::vector<int32_t> tids(type_ids.begin(), type_ids.end());
  std::vector<int32_t> etys(etypes.begin(), etypes.end());
  std::vector<int32_t> dynf(dyn_flags.begin(), dyn_flags.end());
  std::vector<int32_t> dyci(dyc_idx.begin(), dyc_idx.end());
  std::vector<uint8_t> vm(vmiss.size());
  for (R_xlen_t i = 0; i < vmiss.size(); ++i) vm[i] = (uint8_t)(vmiss[i] != 0);

  cusna_saom_problem p;
  std::memset(&p, 0, sizeof(p));
  p.struct_version = cusna_abi_version();
  p.n = n; p.W = W; p.k = k; p.M = M; p.threads = 0;
  p.x0_out = x0_out.data(); p.x0_in = x0_in.data();
  p.lam = lamv.data(); p.theta = thv.data();
  p.type_ids = tids.data(); p.etypes = etys.data();
  p.dyn_flags = dynf.data(); p.dyc_idx = dyci.data();
  p.vmiss = vm.data();
  p.vmat64 = vmat64.data(); p.wtab64 = wtab64.data(); p.vb64 = vb64.data();
  p.p1_64 = p1_64.data(); p.p2_64 = p2_64.data();
  p.maxdegree = maxdegree; p.conditional = conditional;
  p.target_dist = target_dist; p.symflag = symflag;
  p.seed = (uint64_t)seed; p.nact = n; p.brange = 1.0f;

  // outputs
  writable::doubles stats((R_xlen_t)M * k);
  writable::doubles scores_out((R_xlen_t)M * (k + 1));
  writable::doubles ntim_out((R_xlen_t)M);
  writable::integers dist_out((R_xlen_t)M);
  writable::integers dist_stat_out((R_xlen_t)M);
  writable::integers steps_out((R_xlen_t)M);
  std::vector<float> scores((size_t)M * (k + 1));
  std::vector<int32_t> dist(M), dist_stat(M), steps(M);
  std::vector<double> ntim(M);

  cusna_saom_result r;
  std::memset(&r, 0, sizeof(r));
  r.struct_version = cusna_abi_version();
  r.stats = REAL(stats); r.scores = scores.data();
  r.ntim = ntim.data(); r.dist = dist.data();
  r.dist_stat = dist_stat.data(); r.steps = steps.data();

  // behavior sub-block
  std::vector<float> lambv, thbv;
  std::vector<int32_t> btyv, z0v, bkv;
  std::vector<float> bscores;
  std::vector<int32_t> bdist, z_end;
  if (kb > 0) {
    lambv.assign((size_t)M, (float)lamb);
    thbv.resize((size_t)M * kb);
    for (int c = 0; c < M; ++c)
      for (int e = 0; e < kb; ++e) thbv[(size_t)c * kb + e] = (float)thetab[e];
    btyv.assign(bty.begin(), bty.end());
    z0v.assign(z0.begin(), z0.end());
    bkv.assign(bkdist.begin(), bkdist.end());
    p.kb = kb; p.lamb = lambv.data(); p.thetab = thbv.data();
    p.bty = btyv.data(); p.z0 = z0v.data(); p.bkdist = bkv.data();
    p.vb64 = vb64.data();
    p.zmin = zmin; p.zmax = zmax;
    p.zbar = (float)zbar; p.brange = (float)brange; p.bsimmean = (float)bsimmean;
    bscores.assign((size_t)M * (kb + 1), 0.f);
    bdist.assign(M, 0); z_end.assign((size_t)M * n, 0);
    r.bdist = bdist.data(); r.bscores = bscores.data(); r.z_end = z_end.data();
  }

  std::vector<uint8_t> a_ends;
  if (want_ends) {
    a_ends.assign((size_t)M * n * n, 0);
    r.want_end_nets = 1; r.a_ends = a_ends.data();
  }

  cusna_status st = cusna_saom_simulate_period_cpu(&p, &r);
  if (st != CUSNA_OK) stop("cusna_saom_simulate_period_cpu failed (status %d)", (int)st);

  // copy native buffers into the R outputs
  for (int i = 0; i < M * (k + 1); ++i) scores_out[i] = (double)scores[i];
  for (int i = 0; i < M; ++i) {
    ntim_out[i] = ntim[i]; dist_out[i] = dist[i];
    dist_stat_out[i] = dist_stat[i]; steps_out[i] = steps[i];
  }

  writable::list out;
  out.push_back({"stats"_nm = stats});
  out.push_back({"scores"_nm = scores_out});
  out.push_back({"ntim"_nm = ntim_out});
  out.push_back({"dist"_nm = dist_out});
  out.push_back({"dist_stat"_nm = dist_stat_out});
  out.push_back({"steps"_nm = steps_out});
  if (kb > 0) {
    writable::doubles bsc((R_xlen_t)M * (kb + 1));
    writable::integers bd((R_xlen_t)M);
    writable::integers ze((R_xlen_t)M * n);
    for (int i = 0; i < M * (kb + 1); ++i) bsc[i] = (double)bscores[i];
    for (int i = 0; i < M; ++i) bd[i] = bdist[i];
    for (int i = 0; i < M * n; ++i) ze[i] = z_end[i];
    out.push_back({"bscores"_nm = bsc});
    out.push_back({"bdist"_nm = bd});
    out.push_back({"z_end"_nm = ze});
  }
  if (want_ends) {
    writable::integers ae((R_xlen_t)M * n * n);
    for (R_xlen_t i = 0; i < (R_xlen_t)M * n * n; ++i) ae[i] = a_ends[i];
    out.push_back({"a_ends"_nm = ae});
  }
  return out;
}

// Typed (creation/endowment) statistic over a batch of end networks
// (cusna_typed_stats; R3b). a_batch is M*n*n (row-major), F/x0/masks dense n*n.
[[cpp11::register]]
doubles cusna_typed_stats_(integers a_batch, int n, int M, int etype,
                           doubles F, integers x0, integers keep, integers force,
                           integers keepS, integers forceS, integers miss_m,
                           integers miss_m1) {
  std::vector<uint8_t> ab = dense_u8(a_batch), x0d = dense_u8(x0);
  std::vector<uint8_t> kp = dense_u8(keep), fo = dense_u8(force);
  std::vector<uint8_t> kpS = dense_u8(keepS), foS = dense_u8(forceS);
  std::vector<uint8_t> mm = dense_u8(miss_m), mm1 = dense_u8(miss_m1);
  writable::doubles out((R_xlen_t)M);
  cusna_status st = cusna_typed_stats(ab.data(), n, M, etype, REAL(F), x0d.data(),
      kp.data(), fo.data(), kpS.data(), foS.data(), mm.data(), mm1.data(),
      REAL(out));
  if (st != CUSNA_OK) stop("cusna_typed_stats failed (status %d)", (int)st);
  return out;
}

// Cross-network statistic over a batch of end networks (cusna_cross_stats; R3b).
[[cpp11::register]]
doubles cusna_cross_stats_(integers a_batch, integers wp, int n, int M, int t,
                           double parm) {
  std::vector<uint8_t> ab = dense_u8(a_batch), wpd = dense_u8(wp);
  writable::doubles out((R_xlen_t)M);
  cusna_status st = cusna_cross_stats(ab.data(), wpd.data(), n, M, t, parm, REAL(out));
  if (st != CUSNA_OK) stop("cusna_cross_stats failed (status %d)", (int)st);
  return out;
}

// One multi-network SAOM period, M chains, CPU backend (cusna_saom_simulate_
// multi_cpu; R4b). `x0`/`perm` are D*n*n dense (row-major); lam (length D) and
// theta (length k) broadcast to the M chains; descriptors are float32 (the
// kernel path). Returns the per-net rate scores, effect scores, elapsed time,
// ministeps and the D*M end networks -- the eval / cross statistics are computed
// in R from a_ends (cusna_network_stats + cusna_cross_stats), matching the host
// split in cpu.py::simulate_batch_cpu_multi.
[[cpp11::register]]
list cusna_simulate_multi_cpu_(
    int n, int k, int M, int D, double seed,
    integers x0, doubles lam, doubles theta,
    integers type_ids, integers etypes, integers netof, integers netref,
    doubles wtab, doubles p1, doubles p2, doubles vmat,
    integers perm, integers symf) {
  const int W = (n + 31) / 32;
  std::vector<uint8_t> x0d = dense_u8(x0), permd = dense_u8(perm);

  std::vector<uint32_t> x0_out((size_t)D * n * W), x0_in((size_t)D * n * W),
                        permb((size_t)D * n * W);
  for (int d = 0; d < D; ++d) {
    std::vector<uint32_t> o, in;
    pack_out_in(x0d.data() + (size_t)d * n * n, n, W, o, in);
    std::copy(o.begin(), o.end(), x0_out.begin() + (size_t)d * n * W);
    std::copy(in.begin(), in.end(), x0_in.begin() + (size_t)d * n * W);
    pack_single(permd.data() + (size_t)d * n * n, n, W,
                permb.data() + (size_t)d * n * W);
  }

  std::vector<float> lamv((size_t)M * D), thv((size_t)M * k);
  for (int c = 0; c < M; ++c) {
    for (int d = 0; d < D; ++d) lamv[(size_t)c * D + d] = (float)lam[d];
    for (int e = 0; e < k; ++e) thv[(size_t)c * k + e] = (float)theta[e];
  }
  std::vector<float> wtabf(wtab.begin(), wtab.end()), p1f(p1.begin(), p1.end()),
                     p2f(p2.begin(), p2.end()), vmatf(vmat.begin(), vmat.end());
  std::vector<int32_t> tids(type_ids.begin(), type_ids.end()),
    etys(etypes.begin(), etypes.end()), nof(netof.begin(), netof.end()),
    nref(netref.begin(), netref.end()), symv(symf.begin(), symf.end());
  std::vector<int32_t> dyci((size_t)k, -1);            // no dyadic covariates here
  std::vector<int32_t> actv(n);
  for (int i = 0; i < n; ++i) actv[i] = i;

  cusna_saom_multi_problem p;
  std::memset(&p, 0, sizeof(p));
  p.struct_version = cusna_abi_version();
  p.n = n; p.W = W; p.k = k; p.M = M; p.D = D; p.threads = 0;
  p.x0_out = x0_out.data(); p.x0_in = x0_in.data();
  p.lam = lamv.data(); p.theta = thv.data();
  p.type_ids = tids.data(); p.etypes = etys.data();
  p.netof = nof.data(); p.netref = nref.data();
  p.wtab = wtabf.data(); p.p1 = p1f.data(); p.p2 = p2f.data(); p.vmat = vmatf.data();
  p.dyc_idx = dyci.data();                              // all -1 (dyc = NULL)
  p.perm = permb.data(); p.symf = symv.data();
  p.seed = (uint64_t)seed; p.act = actv.data(); p.nact = n;

  std::vector<float> rscores((size_t)M * D), scores((size_t)M * k);
  std::vector<double> ntim(M);
  std::vector<int32_t> steps(M);
  std::vector<uint8_t> a_ends((size_t)M * D * n * n);
  cusna_saom_multi_result r;
  std::memset(&r, 0, sizeof(r));
  r.struct_version = cusna_abi_version();
  r.rscores = rscores.data(); r.scores = scores.data();
  r.ntim = ntim.data(); r.steps = steps.data(); r.a_ends = a_ends.data();

  cusna_status st = cusna_saom_simulate_multi_cpu(&p, &r);
  if (st != CUSNA_OK) stop("cusna_saom_simulate_multi_cpu failed (status %d)", (int)st);

  writable::doubles rsc_out((R_xlen_t)M * D), sc_out((R_xlen_t)M * k);
  writable::doubles ntim_out((R_xlen_t)M);
  writable::integers steps_out((R_xlen_t)M);
  writable::integers ae_out((R_xlen_t)M * D * n * n);
  for (R_xlen_t i = 0; i < (R_xlen_t)M * D; ++i) rsc_out[i] = (double)rscores[i];
  for (R_xlen_t i = 0; i < (R_xlen_t)M * k; ++i) sc_out[i] = (double)scores[i];
  for (int i = 0; i < M; ++i) { ntim_out[i] = ntim[i]; steps_out[i] = steps[i]; }
  for (R_xlen_t i = 0; i < (R_xlen_t)M * D * n * n; ++i) ae_out[i] = a_ends[i];

  writable::list out;
  out.push_back({"rscores"_nm = rsc_out});
  out.push_back({"scores"_nm = sc_out});
  out.push_back({"ntim"_nm = ntim_out});
  out.push_back({"steps"_nm = steps_out});
  out.push_back({"a_ends"_nm = ae_out});
  return out;
}
