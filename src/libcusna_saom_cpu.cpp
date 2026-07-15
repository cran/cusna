/* ==========================================================================
 *  libcusna_saom_cpu.cpp -- CPU mirror of the SAOM period simulator (A3).
 *
 *  Pure C++ (no CUDA): compiles for the CRAN / macOS CPU-only build. A faithful
 *  float64 port of prototype/cusna/saom/cpu.py::_simulate_one (the reference the
 *  23-check harness certifies), driven behind the SAME C ABI struct as the GPU
 *  path (cusna_saom_simulate_period_cpu). Chains are independent -> OpenMP.
 *
 *  Scope (A3): the simulation core -- trajectory, MoM scores, distances, the
 *  end network, behavior and rate outputs. The masked EVALUATION statistics are
 *  computed on the host by batch_stats (stats.py) for BOTH backends; that host
 *  port is milestone A4, so result->stats is left untouched here.
 *
 *  Bit-identity target: cpu.py float64. Everything is bit-reproducible except
 *  where the C++ and numba runtime math libraries round exp/log differently;
 *  sqrt is IEEE-correctly-rounded so degree effects are exact regardless.
 * ========================================================================== */
#include <cmath>
#include <cstdint>
#include <cstdlib>   // std::labs (gcc/libstdc++ is stricter than MSVC here)
#include <cstring>
#include <vector>

#include "libcusna.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

inline void philox4(uint32_t c0, uint32_t c1, uint32_t c2, uint32_t c3,
                    uint32_t k0, uint32_t k1, uint32_t out[4]) {
    for (int i = 0; i < 10; ++i) {
        uint64_t p0 = (uint64_t)0xD2511F53u * c0;
        uint64_t p1 = (uint64_t)0xCD9E8D57u * c2;
        uint32_t n0 = (uint32_t)((p1 >> 32) ^ c1 ^ k0);
        uint32_t n1 = (uint32_t)(p1 & 0xFFFFFFFFu);
        uint32_t n2 = (uint32_t)((p0 >> 32) ^ c3 ^ k1);
        uint32_t n3 = (uint32_t)(p0 & 0xFFFFFFFFu);
        c0 = n0; c1 = n1; c2 = n2; c3 = n3;
        k0 += 0x9E3779B9u; k1 += 0xBB67AE85u;
    }
    out[0] = c0; out[1] = c1; out[2] = c2; out[3] = c3;
}
inline double u01(uint32_t x) { return ((double)x + 0.5) * (1.0 / 4294967296.0); }

// float64 mirror of cpu.py::_ego_ds. a is dense n*n row-major uint8.
double ego_ds(int t, const uint8_t *a, int n, const long *indeg,
              const long *outdeg, long rdeg_i, const int *tp,
              double vi, double vj, double dyw, double p1, double p2,
              int i, int j, double delta, const double *wrow) {
    const uint8_t *ai = a + (size_t)i * n;
    const uint8_t *aj = a + (size_t)j * n;
    switch (t) {
        case 22: return delta * dyw;
        case 0:  return delta;
        case 1:  return delta * (double)aj[i];
        case 2: { long tt = 0; for (int h = 0; h < n; ++h) if (ai[h]) { if (a[(size_t)h*n+j]) ++tt; if (aj[h]) ++tt; } return delta * (double)tt; }
        case 3: { long c3 = 0; for (int h = 0; h < n; ++h) if (aj[h] && a[(size_t)h*n+i]) ++c3; return delta * (double)c3; }
        case 4:  return delta * (double)(indeg[j] + (1 - ai[j]));
        case 5: { int target = ai[j] ? 1 : 0; long cnt = 0;
                  for (int h = 0; h < n; ++h) if (ai[h] && aj[h] && tp[h] == target) ++cnt;
                  if (tp[j] > 0) ++cnt; return delta * (double)cnt; }
        case 6: { long cnt = 0; if (aj[i]) for (int h = 0; h < n; ++h) if (ai[h] && a[(size_t)h*n+j]) ++cnt;
                  for (int h = 0; h < n; ++h) if (ai[h] && a[(size_t)h*n+i] && aj[h]) ++cnt; return delta * (double)cnt; }
        case 7:  return delta * std::sqrt((double)(indeg[j] + (1 - ai[j])));
        case 8:  return delta * (double)outdeg[j];
        case 9:  return delta * std::sqrt((double)outdeg[j]);
        case 10: return delta * vi;
        case 11: return delta * vj;
        case 12: return delta * (1.0 - std::fabs(vi - vj) * p1 - p2);
        case 13: return delta * (std::fabs(vi - vj) < 1e-9 ? 1.0 : 0.0);
        case 14: return delta * (double)indeg[i];
        case 15: return delta * std::sqrt((double)indeg[i]);
        case 16: return delta * (2.0 * (double)outdeg[i] + delta);
        case 17: { double d = (double)outdeg[i], dn = d + delta; return dn * std::sqrt(dn) - d * std::sqrt(d); }
        case 18: { double c = p1, d = (double)outdeg[i]; return std::fmin(d + delta, c) - std::fmin(d, c); }
        case 19: return delta * ((double)rdeg_i + (double)aj[i] * ((double)outdeg[i] + delta));
        case 20: if (vi > vj + 1e-9) return delta; if (std::fabs(vi - vj) <= 1e-9) return delta * 0.5; return 0.0;
        case 21: return delta * vi * vj;
        case 23: { long c3 = 0; for (int h = 0; h < n; ++h) if (aj[h] && a[(size_t)h*n+i]) ++c3; return delta * (double)c3; }
        case 24: { int tpf = 0; for (int h = 0; h < n; ++h) if (ai[h] && a[(size_t)h*n+j]) ++tpf; return delta * wrow[tpf]; }
        case 25: return delta * (vi + vj);
        case 26: { int c = 0; for (int h = 0; h < n; ++h) if (ai[h] && aj[h]) ++c; return delta * wrow[c]; }
        case 27: { int c = 0; for (int h = 0; h < n; ++h) if (a[(size_t)h*n+i] && a[(size_t)h*n+j]) ++c; return delta * wrow[c]; }
        case 28: { int c = 0; for (int h = 0; h < n; ++h) if (a[(size_t)h*n+i] && aj[h]) ++c; return delta * wrow[c]; }
        case 29: { int xij = ai[j]; double acc = 0.0;
                   for (int h = 0; h < n; ++h) if (aj[h] && h != i) {
                       int p = 0; for (int m = 0; m < n; ++m) if (ai[m] && a[(size_t)m*n+h]) ++p;
                       acc += xij ? (wrow[p] - wrow[p-1]) : (wrow[p+1] - wrow[p]); }
                   return delta * acc; }
        case 30: { int xij = ai[j]; double acc = 0.0;
                   for (int h = 0; h < n; ++h) if (a[(size_t)h*n+j] && h != i) {
                       int p = 0; for (int m = 0; m < n; ++m) if (ai[m] && a[(size_t)h*n+m]) ++p;
                       acc += xij ? (wrow[p] - wrow[p-1]) : (wrow[p+1] - wrow[p]); }
                   return delta * acc; }
        case 31: return delta * (vj >= p1 - 1e-12 ? 1.0 : 0.0);
    }
    return 0.0;
}

double beh_ds(int t, double dlt, double zi_c, long nb, double zsum,
              long nH, long nE, long nL, double brange, double bsimmean) {
    switch (t) {
        case 0: return dlt;
        case 1: return (2.0 * zi_c + dlt) * dlt;
        case 2: return nb > 0 ? dlt * zsum / (double)nb : 0.0;
        case 3: { if (nb == 0) return 0.0; double tc = dlt > 0 ? (double)(nH - nE - nL) : (double)(nL - nE - nH); return tc / (brange * (double)nb); }
        case 4: { if (nb == 0) return 0.0; double tc = dlt > 0 ? (double)(nH - nE - nL) : (double)(nL - nE - nH); return tc / brange - (double)nb * bsimmean; }
    }
    return 0.0;
}

const int MAX_STEPS = 4000000;

// All problem-level (chain-invariant) inputs, unpacked to dense form once.
struct Ctx {
    int n, k, kb, nrate;
    const uint8_t *a0;      // data start network (dense n*n) -- for distance
    const uint8_t *a0sim;   // simulation start network (dense n*n)
    const uint8_t *kdist;   // dyads counted in distance (dense n*n)
    const uint8_t *perm;    // toggleable dyads (dense n*n)
    double lam, lamb;
    const double *th, *thb;
    const int *type_ids, *etys, *dynf, *dyci, *cint0, *cint1, *cint2, *bty, *rkind;
    const double *wtab;     // k*(n+1) promoted to double
    const double *p1, *p2, *vmat, *dymat, *ralpha, *ru;
    const int *z0;
    int zmin, zmax; double zbar, brange, bsimmean;
    const int *bkdist;
    const int *act; int nact; const uint8_t *inact; int has_comp, symflag, maxdeg;
    int cond, target_dist; uint32_t k0, k1;
    int need_tp, need_rdeg, nd;
};

struct ChainOut {
    long dist, dist_stat, steps; double t, rate_score, brate_score, bdist;
};

// float64 mirror of cpu.py::_simulate_one for one chain.
void simulate_one(const Ctx &C, uint32_t chain, uint8_t *a, int *z,
                  double *scores, double *bscores, double *rsc,
                  long *adist, ChainOut &o) {
    const int n = C.n, k = C.k, kb = C.kb, nrate = C.nrate;
    std::vector<long> indeg(n, 0), outdeg(n, 0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            a[(size_t)i*n+j] = C.a0sim[(size_t)i*n+j];
            if (a[(size_t)i*n+j]) { indeg[j]++; outdeg[i]++; }
        }
    for (int i = 0; i < n; ++i) z[i] = C.z0[i];
    for (int e = 0; e < k; ++e) scores[e] = 0.0;
    for (int e = 0; e < kb; ++e) bscores[e] = 0.0;
    for (int q = 0; q < nrate; ++q) rsc[q] = 0.0;

    std::vector<double> df(n), wsum(k), crate(n, 1.0), cup(kb ? kb : 1), cdn(kb ? kb : 1);
    std::vector<double> dsbuf((size_t)n * k), ubuf(nrate ? nrate : 1), rsum(nrate ? nrate : 1);
    std::vector<uint8_t> okbuf(n);
    std::vector<int> tp(n);

    long dist = 0, bdist = 0, msteps = 0, msteps_beh = 0;
    double t = 0.0, rho_tau = 0.0;
    uint32_t rr[4];

    while (true) {
        philox4((uint32_t)msteps, chain, 0u, 0u, C.k0, C.k1, rr);
        int i;
        bool is_net;
        if (nrate > 0) {
            double S0 = 0.0;
            for (int q = 0; q < nrate; ++q) rsum[q] = 0.0;
            for (int idx = 0; idx < C.nact; ++idx) {
                int i2 = C.act[idx]; double s = 0.0;
                for (int q = 0; q < nrate; ++q) {
                    double u = C.rkind[q] == 0 ? C.ru[(size_t)q*n+i2]
                             : C.rkind[q] == 1 ? (double)outdeg[i2] : (double)indeg[i2];
                    ubuf[q] = u; s += C.ralpha[q] * u;
                }
                if (s > 30.0) s = 30.0; else if (s < -30.0) s = -30.0;
                double c = std::exp(s); crate[i2] = c; S0 += c;
                for (int q = 0; q < nrate; ++q) rsum[q] += c * ubuf[q];
            }
            double tau = -std::log(u01(rr[0])) / (C.lam * S0);
            double tau_eff = tau;
            if (!C.cond && t + tau > 1.0) tau_eff = 1.0 - t;
            rho_tau += tau_eff * S0;
            for (int q = 0; q < nrate; ++q) rsc[q] -= tau_eff * C.lam * rsum[q];
            t += tau;
            if (!C.cond && t > 1.0) break;
            is_net = true;
            double tgt = u01(rr[2]) * S0, acc = 0.0; i = C.act[C.nact - 1];
            for (int idx = 0; idx < C.nact; ++idx) { acc += crate[C.act[idx]]; if (acc >= tgt) { i = C.act[idx]; break; } }
        } else {
            t += -std::log(u01(rr[0])) / ((double)C.nact * (C.lam + C.lamb));
            if (!C.cond && t > 1.0) break;
            is_net = kb == 0 || u01(rr[1]) < C.lam / (C.lam + C.lamb);
            long idx = (long)(u01(rr[2]) * C.nact); if (idx >= C.nact) idx = C.nact - 1;
            i = C.act[idx];
        }

        if (!is_net) {  // behavior ministep
            int zi = z[i]; double zi_c = (double)zi - C.zbar;
            long nb = 0, nH = 0, nE = 0, nL = 0; double zsum = 0.0;
            for (int j = 0; j < n; ++j) if (a[(size_t)i*n+j]) {
                nb++; int zj = z[j]; zsum += (double)zj - C.zbar;
                if (zj > zi) nH++; else if (zj == zi) nE++; else nL++;
            }
            bool up_ok = zi < C.zmax, dn_ok = zi > C.zmin;
            double f_up = 0.0, f_dn = 0.0;
            for (int e = 0; e < kb; ++e) {
                cup[e] = beh_ds(C.bty[e], 1.0, zi_c, nb, zsum, nH, nE, nL, C.brange, C.bsimmean);
                cdn[e] = beh_ds(C.bty[e], -1.0, zi_c, nb, zsum, nH, nE, nL, C.brange, C.bsimmean);
                f_up += C.thb[e] * cup[e]; f_dn += C.thb[e] * cdn[e];
            }
            double mx = 0.0; if (up_ok && f_up > mx) mx = f_up; if (dn_ok && f_dn > mx) mx = f_dn;
            double pdn = dn_ok ? std::exp(f_dn - mx) : 0.0, p0 = std::exp(-mx), pup = up_ok ? std::exp(f_up - mx) : 0.0;
            double psum = pdn + p0 + pup; pdn /= psum; p0 /= psum; pup /= psum;
            double uc = u01(rr[3]); int d = uc < pdn ? -1 : (uc < pdn + p0 ? 0 : 1);
            for (int e = 0; e < kb; ++e) { double s = d == 1 ? cup[e] : (d == -1 ? cdn[e] : 0.0);
                if (up_ok) s -= pup * cup[e]; if (dn_ok) s -= pdn * cdn[e]; bscores[e] += s; }
            if (d != 0) { if (C.bkdist[i]) bdist += std::labs(zi + d - C.z0[i]) - std::labs((long)zi - C.z0[i]); z[i] = zi + d; }
            msteps_beh++; msteps++;
            if (C.cond && dist >= C.target_dist) break;
            if (msteps >= MAX_STEPS) break;
            continue;
        }

        if (C.need_tp) { for (int h = 0; h < n; ++h) tp[h] = 0;
            for (int x = 0; x < n; ++x) if (a[(size_t)i*n+x]) for (int h = 0; h < n; ++h) tp[h] += a[(size_t)x*n+h]; }
        long rdeg_i = 0;
        if (C.need_rdeg) for (int h = 0; h < n; ++h) if (a[(size_t)i*n+h] && a[(size_t)h*n+i]) rdeg_i++;
        double zvi = (double)z[i] - C.zbar;

        double gbest = -1e30, lmax = -1e30; int jstar = i;
        for (int j = 0; j < n; ++j) {
            double f = 0.0; bool ok = true;
            if (j == i) { for (int e = 0; e < k; ++e) dsbuf[(size_t)j*k+e] = 0.0; }
            else {
                if (!C.perm[(size_t)i*n+j]) ok = false;
                else if (C.maxdeg > 0 && a[(size_t)i*n+j] == 0 && outdeg[i] >= C.maxdeg) ok = false;
                if (ok) {
                    double delta = 1.0 - 2.0 * a[(size_t)i*n+j];
                    double zvj = (double)z[j] - C.zbar;
                    for (int e = 0; e < k; ++e) {
                        if (C.type_ids[e] == 32) {  // interaction
                            double d = C.cint2[e] < 0 ? delta : 1.0;
                            d *= dsbuf[(size_t)j*k+C.cint0[e]] * dsbuf[(size_t)j*k+C.cint1[e]];
                            if (C.cint2[e] >= 0) d *= dsbuf[(size_t)j*k+C.cint2[e]];
                            dsbuf[(size_t)j*k+e] = d; f += C.th[e] * d; continue;
                        }
                        double vi = C.dynf[e] ? zvi : C.vmat[(size_t)e*n+i];
                        double vj = C.dynf[e] ? zvj : C.vmat[(size_t)e*n+j];
                        double dyw = C.dyci[e] >= 0 ? C.dymat[((size_t)C.dyci[e]*n+i)*n+j] : 0.0;
                        double d = ego_ds(C.type_ids[e], a, n, indeg.data(), outdeg.data(),
                                          rdeg_i, tp.data(), vi, vj, dyw, C.p1[e], C.p2[e],
                                          i, j, delta, C.wtab + (size_t)e*(n+1));
                        if ((C.etys[e] == 1 && a[(size_t)i*n+j] == 1) ||
                            (C.etys[e] == 2 && a[(size_t)i*n+j] == 0)) d = 0.0;
                        dsbuf[(size_t)j*k+e] = d; f += C.th[e] * d;
                    }
                }
            }
            okbuf[j] = ok ? 1 : 0; df[j] = f;
            if (ok && f > lmax) lmax = f;
            philox4((uint32_t)msteps, chain, 1u, (uint32_t)j, C.k0, C.k1, rr);
            if (ok) { double g = f - std::log(-std::log(u01(rr[0])));
                if (g > gbest || (g == gbest && j < jstar)) { gbest = g; jstar = j; } }
        }

        double sume = 0.0; for (int e = 0; e < k; ++e) wsum[e] = 0.0;
        for (int j = 0; j < n; ++j) if (okbuf[j]) { double e_ = std::exp(df[j] - lmax); sume += e_;
            for (int e = 0; e < k; ++e) wsum[e] += e_ * dsbuf[(size_t)j*k+e]; }
        for (int e = 0; e < k; ++e) scores[e] += dsbuf[(size_t)jstar*k+e] - wsum[e] / sume;
        msteps++;

        if (jstar != i) {
            int j = jstar; long delta = 1 - 2 * a[(size_t)i*n+j]; long chg = C.symflag ? 2 : 1;
            if (C.kdist[(size_t)i*n+j]) dist += (a[(size_t)i*n+j] == C.a0[(size_t)i*n+j]) ? chg : -chg;
            indeg[j] += delta; outdeg[i] += delta; a[(size_t)i*n+j] = 1 - a[(size_t)i*n+j];
            if (C.symflag) { indeg[i] += delta; outdeg[j] += delta; a[(size_t)j*n+i] = 1 - a[(size_t)j*n+i]; }
        }
        if (nrate > 0) for (int q = 0; q < nrate; ++q)
            rsc[q] += C.rkind[q] == 0 ? C.ru[(size_t)q*n+i] : C.rkind[q] == 1 ? (double)outdeg[i] : (double)indeg[i];

        if (C.cond && dist >= C.target_dist) break;
        if (msteps >= MAX_STEPS) break;
    }

    if (C.has_comp) for (int i = 0; i < n; ++i) if (C.inact[i])
        for (int j = 0; j < n; ++j) { a[(size_t)i*n+j] = C.a0[(size_t)i*n+j]; a[(size_t)j*n+i] = C.a0[(size_t)j*n+i]; }

    long dist_stat = 0;
    for (int i = 0; i < n; ++i) { long cnt = 0;
        for (int j = 0; j < n; ++j) if (C.kdist[(size_t)i*n+j] && a[(size_t)i*n+j] != C.a0[(size_t)i*n+j]) cnt++;
        adist[i] = cnt; dist_stat += cnt; }

    o.dist = dist; o.dist_stat = dist_stat; o.steps = msteps; o.t = t;
    o.rate_score = nrate > 0 ? (double)(msteps - msteps_beh) / C.lam - rho_tau
                             : (double)(msteps - msteps_beh) / C.lam - (double)C.nact;
    o.brate_score = kb > 0 ? (double)msteps_beh / C.lamb - (double)C.nact : 0.0;
    o.bdist = (double)bdist;
}

// unpack a packed bit-row array (n*W uint32) to a dense n*n uint8 (row-major)
void unpack(const uint32_t *rows, int n, int W, uint8_t *dense) {
    std::memset(dense, 0, (size_t)n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            dense[(size_t)i*n+j] = (rows[(size_t)i*W + (j >> 5)] >> (j & 31)) & 1u;
}

}  // namespace

extern "C" CUSNA_API cusna_status
cusna_saom_simulate_period_cpu(const cusna_saom_problem *pb, cusna_saom_result *r) {
    if (!pb || !r) return CUSNA_ERR_NULL_ARG;
    if (pb->struct_version != CUSNA_ABI_VERSION) return CUSNA_ERR_SHAPE;
    const int n = pb->n, W = pb->W, k = pb->k, M = pb->M, kb = pb->kb, nrate = pb->nrate;
    if (n < 1 || n > CUSNA_MAX_N || W != (n + 31) / 32 || k < 1 || M < 1) return CUSNA_ERR_SHAPE;

    // unpack networks + masks to dense (cpu.py works on dense a[i,j])
    std::vector<uint8_t> a0(n * n), a0sim(n * n), kdist(n * n), perm(n * n);
    unpack(pb->x0_out, n, W, a0.data());
    unpack(pb->sim_out ? pb->sim_out : pb->x0_out, n, W, a0sim.data());
    if (pb->kdist) unpack(pb->kdist, n, W, kdist.data()); else { for (int i=0;i<n;++i) for (int j=0;j<n;++j) kdist[(size_t)i*n+j] = (i!=j); }
    if (pb->perm) unpack(pb->perm, n, W, perm.data()); else { for (int i=0;i<n;++i) for (int j=0;j<n;++j) perm[(size_t)i*n+j] = (i!=j); }
    // statistic-network masks (A & keep) | force for the eval statistics (A4.2)
    std::vector<uint8_t> keepd((size_t)n*n, 1), forced((size_t)n*n, 0);
    if (pb->keep_out) unpack(pb->keep_out, n, W, keepd.data());
    else for (int i=0;i<n;++i) for (int j=0;j<n;++j) keepd[(size_t)i*n+j] = (i!=j);
    if (pb->force_out) unpack(pb->force_out, n, W, forced.data());

    // theta / thetab are per-chain (M*k) but constant across chains here; take row 0
    std::vector<double> th(k); for (int e = 0; e < k; ++e) th[e] = pb->theta[e];
    std::vector<double> thb(kb ? kb : 1, 0.0); for (int e = 0; e < kb; ++e) thb[e] = pb->thetab ? pb->thetab[e] : 0.0;
    // descriptors -> double; prefer the exact float64 fields (A4.3) when given,
    // else promote the float32 kernel descriptors
    std::vector<double> wtab((size_t)k * (n + 1), 0.0);
    for (size_t q = 0; q < (size_t)k * (n + 1); ++q)
        wtab[q] = pb->wtab64 ? pb->wtab64[q] : (pb->wtab ? pb->wtab[q] : 0.0);
    std::vector<double> p1(k), p2(k), vmat((size_t)k * n);
    for (int e = 0; e < k; ++e) {
        p1[e] = pb->p1_64 ? pb->p1_64[e] : (pb->p1 ? pb->p1[e] : 0.0);
        p2[e] = pb->p2_64 ? pb->p2_64[e] : (pb->p2 ? pb->p2[e] : 0.0);
    }
    for (size_t q = 0; q < (size_t)k * n; ++q)
        vmat[q] = pb->vmat64 ? pb->vmat64[q] : (pb->vmat ? pb->vmat[q] : 0.0);
    std::vector<double> vb(n, 0.0);
    for (int i = 0; i < n; ++i) vb[i] = pb->vb64 ? pb->vb64[i] : (pb->vb ? pb->vb[i] : 0.0);
    std::vector<double> ralpha(nrate ? nrate : 1, 0.0), ru((size_t)(nrate ? nrate : 1) * n, 0.0);
    for (int q = 0; q < nrate; ++q) ralpha[q] = pb->ralpha[q];
    for (size_t q = 0; nrate && q < (size_t)nrate * n; ++q) ru[q] = pb->ru[q];
    std::vector<double> dymat(1 * n * n, 0.0);  // A2/A3 configs carry no dyadic covariate

    std::vector<int> z0(n, 0); for (int i = 0; i < n && pb->z0; ++i) z0[i] = pb->z0[i];
    std::vector<int> bkdist(n, 1); for (int i = 0; i < n && pb->bkdist; ++i) bkdist[i] = pb->bkdist[i];
    std::vector<int> act(n); int nact = pb->act ? pb->nact : n;
    if (pb->act) for (int i = 0; i < n; ++i) act[i] = pb->act[i]; else for (int i = 0; i < n; ++i) act[i] = i;
    std::vector<uint8_t> inact(n, 0);
    if (pb->has_comp && pb->inactw) for (int i = 0; i < n; ++i) inact[i] = (pb->inactw[i >> 5] >> (i & 31)) & 1u;

    bool need_tp = false, need_rdeg = false;
    for (int e = 0; e < k; ++e) { if (pb->type_ids[e] == 5) need_tp = true; if (pb->type_ids[e] == 19) need_rdeg = true; }

    Ctx C;
    C.n = n; C.k = k; C.kb = kb; C.nrate = nrate;
    C.a0 = a0.data(); C.a0sim = a0sim.data(); C.kdist = kdist.data(); C.perm = perm.data();
    C.lam = pb->lam ? pb->lam[0] : 1.0; C.lamb = pb->lamb ? pb->lamb[0] : 0.0;
    C.th = th.data(); C.thb = thb.data();
    C.type_ids = pb->type_ids; C.etys = pb->etypes; C.dynf = pb->dyn_flags;
    C.dyci = pb->dyc_idx; C.cint0 = pb->int0; C.cint1 = pb->int1; C.cint2 = pb->int2;
    C.bty = pb->bty; C.rkind = pb->rkind;
    C.wtab = wtab.data(); C.p1 = p1.data(); C.p2 = p2.data(); C.vmat = vmat.data();
    C.dymat = dymat.data(); C.ralpha = ralpha.data(); C.ru = ru.data();
    C.z0 = z0.data(); C.zmin = pb->zmin; C.zmax = pb->zmax;
    C.zbar = pb->zbar; C.brange = pb->brange; C.bsimmean = pb->bsimmean;
    C.bkdist = bkdist.data(); C.act = act.data(); C.nact = nact; C.inact = inact.data();
    C.has_comp = pb->has_comp; C.symflag = pb->symflag; C.maxdeg = pb->maxdegree;
    C.cond = pb->conditional; C.target_dist = pb->target_dist;
    C.k0 = (uint32_t)(pb->seed & 0xFFFFFFFFu); C.k1 = (uint32_t)(pb->seed >> 32);
    C.need_tp = need_tp; C.need_rdeg = need_rdeg; C.nd = 1;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int c = 0; c < M; ++c) {
        std::vector<uint8_t> a((size_t)n * n);
        std::vector<int> z(n);
        std::vector<double> sc(k), bsc(kb ? kb : 1), rs(nrate ? nrate : 1);
        std::vector<long> ad(n);
        ChainOut o;
        simulate_one(C, (uint32_t)c, a.data(), z.data(), sc.data(), bsc.data(), rs.data(), ad.data(), o);
        if (r->dist) r->dist[c] = (int32_t)o.dist;
        if (r->dist_stat) r->dist_stat[c] = (int32_t)o.dist_stat;
        if (r->steps) r->steps[c] = (int32_t)o.steps;
        if (r->ntim) r->ntim[c] = o.t;
        if (r->scores) { r->scores[(size_t)c*(k+1)] = (float)o.rate_score;
            for (int e = 0; e < k; ++e) r->scores[(size_t)c*(k+1)+1+e] = (float)sc[e]; }
        if (kb > 0) {
            if (r->bdist) r->bdist[c] = (int32_t)(long)o.bdist;
            if (r->z_end) for (int i = 0; i < n; ++i) r->z_end[(size_t)c*n+i] = z[i];
            if (r->bscores) { r->bscores[(size_t)c*(kb+1)] = (float)o.brate_score;
                for (int e = 0; e < kb; ++e) r->bscores[(size_t)c*(kb+1)+1+e] = (float)bsc[e]; }
        }
        if (nrate > 0) {
            if (r->rscores) for (int q = 0; q < nrate; ++q) r->rscores[(size_t)c*nrate+q] = (float)rs[q];
            if (r->actor_dist) for (int i = 0; i < n; ++i) r->actor_dist[(size_t)c*n+i] = (int32_t)ad[i];
        }
        if (r->want_end_nets && r->a_ends)
            for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
                r->a_ends[((size_t)c*n+i)*n+j] = a[(size_t)i*n+j];
        // eval statistics on the masked end network (A4.2): (A & keep) | force,
        // then the shared host stat function (float64 descriptors promoted above)
        if (r->stats) {
            std::vector<uint8_t> masked((size_t)n*n);
            for (size_t q = 0; q < (size_t)n*n; ++q) masked[q] = (a[q] & keepd[q]) | forced[q];
            cusna_network_stats(masked.data(), n, k, pb->type_ids, p1.data(), p2.data(),
                                vmat.data(), pb->vmiss, pb->dyn_flags, vb.data(),
                                pb->od0, pb->id0, wtab.data(), pb->symflag,
                                r->stats + (size_t)c * k);
        }
    }
    return CUSNA_OK;
}

/* ---- A4.5: multi-network CPU mirror (cpu.py::_simulate_one_multi) --------- */
namespace {
// cross-network contribution: toggle own-net tie i->j reading source net aw's
// CURRENT state (float64 port of cpu.py::_cross_ds). ax is the own net (for `to`).
double cross_ds(int t, const uint8_t *aw, const uint8_t *ax, int n,
                int i, int j, double delta, double parm) {
    switch (t) {
        case 33: return delta * (double)aw[(size_t)i*n+j];
        case 34: return delta * (double)aw[(size_t)j*n+i];
        case 35: return delta * (double)(aw[(size_t)i*n+j] & aw[(size_t)j*n+i]);
        case 36: { long c = 0; for (int h = 0; h < n; ++h) if (aw[(size_t)i*n+h] && aw[(size_t)j*n+h]) ++c;
                   return delta * (parm >= 2.0 ? std::sqrt((double)c) : (double)c); }
        case 37: { long c = 0; for (int h = 0; h < n; ++h) if (aw[(size_t)i*n+h] && ax[(size_t)h*n+j]) ++c;
                   return delta * (double)c; }
    }
    return 0.0;
}

struct MultiCtx {
    int n, k, D;
    const uint8_t *a0;        // D*n*n data start (== sim start) for restore/dist
    const double *lam;        // per-chain per-net rates (M*D) -- read row c
    const double *th;         // per-chain params (M*k)
    const int *type_ids, *netof, *netref, *etys, *dyci;
    const double *wtab, *p1, *p2, *vmat, *dymat;
    const uint8_t *perm;      // D*n*n toggleable dyads
    const int *symf;          // D
    const int *act; int nact; const uint8_t *inact; int has_comp;
    int need_tp; uint32_t k0, k1;
};

// one unconditional multi-network chain -> end nets a (D*n*n) + scores/rscores.
void simulate_one_multi(const MultiCtx &C, uint32_t chain, uint8_t *a,
                        double *scores, double *rscores, double &ntim, long &steps) {
    const int n = C.n, k = C.k, D = C.D;
    const double *lam_d = C.lam + (size_t)chain * D;
    const double *th = C.th + (size_t)chain * k;
    std::vector<long> indeg((size_t)D * n, 0), outdeg((size_t)D * n, 0);
    for (int d = 0; d < D; ++d)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                uint8_t v = C.a0[((size_t)d*n+i)*n+j];
                a[((size_t)d*n+i)*n+j] = v;
                if (v) { indeg[(size_t)d*n+j]++; outdeg[(size_t)d*n+i]++; }
            }
    double lam_tot = 0.0; for (int d = 0; d < D; ++d) lam_tot += lam_d[d];
    std::vector<long> msteps_d(D, 0);
    for (int e = 0; e < k; ++e) scores[e] = 0.0;
    std::vector<double> df(n), wsum(k), dsbuf((size_t)n * k);
    std::vector<uint8_t> okbuf(n);
    std::vector<int> tp(n);
    long msteps = 0; double t = 0.0; uint32_t rr[4];

    while (true) {
        philox4((uint32_t)msteps, chain, 0u, 0u, C.k0, C.k1, rr);
        t += -std::log(u01(rr[0])) / ((double)C.nact * lam_tot);
        if (t > 1.0) break;
        double uv = u01(rr[1]) * lam_tot; int d = 0; double acc = lam_d[0];
        while (d < D - 1 && uv >= acc) { d += 1; acc += lam_d[d]; }
        long idx = (long)(u01(rr[2]) * C.nact); if (idx >= C.nact) idx = C.nact - 1;
        int i = C.act[idx];
        uint8_t *ad = a + (size_t)d * n * n;          // chosen net dense
        const long *ind = &indeg[(size_t)d*n]; const long *outd = &outdeg[(size_t)d*n];

        if (C.need_tp) { for (int h = 0; h < n; ++h) tp[h] = 0;
            for (int x = 0; x < n; ++x) if (ad[(size_t)i*n+x])
                for (int h = 0; h < n; ++h) tp[h] += ad[(size_t)x*n+h]; }

        double gbest = -1e30, lmax = -1e30; int jstar = i;
        for (int j = 0; j < n; ++j) {
            double f = 0.0; bool ok = true;
            if (j == i) { for (int e = 0; e < k; ++e) dsbuf[(size_t)j*k+e] = 0.0; }
            else if (!C.perm[((size_t)d*n+i)*n+j]) ok = false;
            else {
                double delta = 1.0 - 2.0 * ad[(size_t)i*n+j];
                for (int e = 0; e < k; ++e) {
                    if (C.netof[e] != d) { dsbuf[(size_t)j*k+e] = 0.0; continue; }
                    double dd;
                    if (C.type_ids[e] >= 33)
                        dd = cross_ds(C.type_ids[e], a + (size_t)C.netref[e]*n*n, ad,
                                      n, i, j, delta, C.p1[e]);
                    else {
                        double dyw = C.dyci[e] >= 0 ? C.dymat[((size_t)C.dyci[e]*n+i)*n+j] : 0.0;
                        dd = ego_ds(C.type_ids[e], ad, n, ind, outd, 0, tp.data(),
                                    C.vmat[(size_t)e*n+i], C.vmat[(size_t)e*n+j], dyw,
                                    C.p1[e], C.p2[e], i, j, delta, C.wtab + (size_t)e*(n+1));
                    }
                    if ((C.etys[e] == 1 && ad[(size_t)i*n+j] == 1) ||
                        (C.etys[e] == 2 && ad[(size_t)i*n+j] == 0)) dd = 0.0;
                    dsbuf[(size_t)j*k+e] = dd; f += th[e] * dd;
                }
            }
            okbuf[j] = ok ? 1 : 0; df[j] = f;
            if (ok && f > lmax) lmax = f;
            philox4((uint32_t)msteps, chain, 1u, (uint32_t)j, C.k0, C.k1, rr);
            if (ok) { double g = f - std::log(-std::log(u01(rr[0])));
                if (g > gbest || (g == gbest && j < jstar)) { gbest = g; jstar = j; } }
        }
        double sume = 0.0; for (int e = 0; e < k; ++e) wsum[e] = 0.0;
        for (int j = 0; j < n; ++j) if (okbuf[j]) { double e_ = std::exp(df[j] - lmax); sume += e_;
            for (int e = 0; e < k; ++e) wsum[e] += e_ * dsbuf[(size_t)j*k+e]; }
        for (int e = 0; e < k; ++e) if (C.netof[e] == d) scores[e] += dsbuf[(size_t)jstar*k+e] - wsum[e] / sume;
        msteps++; msteps_d[d]++;
        if (jstar != i) {
            int j = jstar; long delta = 1 - 2 * ad[(size_t)i*n+j];
            indeg[(size_t)d*n+j] += delta; outdeg[(size_t)d*n+i] += delta;
            ad[(size_t)i*n+j] = 1 - ad[(size_t)i*n+j];
            if (C.symf[d]) { indeg[(size_t)d*n+i] += delta; outdeg[(size_t)d*n+j] += delta;
                ad[(size_t)j*n+i] = 1 - ad[(size_t)j*n+i]; }
        }
        if (msteps >= MAX_STEPS) break;
    }
    if (C.has_comp) for (int d = 0; d < D; ++d) for (int i = 0; i < n; ++i) if (C.inact[i])
        for (int j = 0; j < n; ++j) {
            a[((size_t)d*n+i)*n+j] = C.a0[((size_t)d*n+i)*n+j];
            a[((size_t)d*n+j)*n+i] = C.a0[((size_t)d*n+j)*n+i]; }
    for (int d = 0; d < D; ++d) rscores[d] = (double)msteps_d[d] / lam_d[d] - (double)C.nact;
    ntim = t; steps = msteps;
}
}  // namespace

extern "C" CUSNA_API cusna_status
cusna_saom_simulate_multi_cpu(const cusna_saom_multi_problem *pb,
                              cusna_saom_multi_result *r) {
    if (!pb || !r) return CUSNA_ERR_NULL_ARG;
    if (pb->struct_version != CUSNA_ABI_VERSION) return CUSNA_ERR_SHAPE;
    const int n = pb->n, W = pb->W, k = pb->k, M = pb->M, D = pb->D;
    if (n < 1 || n > CUSNA_MAX_N || W != (n + 31) / 32 || k < 1 || k > CUSNA_MAXEFF
        || M < 1 || D < 2 || D > CUSNA_MAXNET) return CUSNA_ERR_SHAPE;

    // unpack the D start networks + per-net permitted masks to dense
    std::vector<uint8_t> a0((size_t)D*n*n), perm((size_t)D*n*n, 1);
    for (int d = 0; d < D; ++d) {
        unpack(pb->x0_out + (size_t)d*n*W, n, W, a0.data() + (size_t)d*n*n);
        if (pb->perm) unpack(pb->perm + (size_t)d*n*W, n, W, perm.data() + (size_t)d*n*n);
        else for (int i=0;i<n;++i) for (int j=0;j<n;++j) perm[((size_t)d*n+i)*n+j] = (i!=j);
    }
    // descriptors float32 -> double
    std::vector<double> lam((size_t)M*D), th((size_t)M*k), wtab((size_t)k*(n+1), 0.0),
        p1(k), p2(k), vmat((size_t)k*n);
    for (size_t q = 0; q < (size_t)M*D; ++q) lam[q] = pb->lam[q];
    for (size_t q = 0; q < (size_t)M*k; ++q) th[q] = pb->theta[q];
    for (size_t q = 0; q < (size_t)k*(n+1); ++q) wtab[q] = pb->wtab ? pb->wtab[q] : 0.0;
    for (int e = 0; e < k; ++e) { p1[e] = pb->p1 ? pb->p1[e] : 0.0; p2[e] = pb->p2 ? pb->p2[e] : 0.0; }
    for (size_t q = 0; q < (size_t)k*n; ++q) vmat[q] = pb->vmat ? pb->vmat[q] : 0.0;
    std::vector<double> dymat(1, 0.0);
    if (pb->dyc) { dymat.assign((size_t)n*n, 0.0); for (size_t q = 0; q < (size_t)n*n; ++q) dymat[q] = pb->dyc[q]; }

    std::vector<int> act(n); int nact = pb->act ? pb->nact : n;
    if (pb->act) for (int i = 0; i < n; ++i) act[i] = pb->act[i]; else for (int i = 0; i < n; ++i) act[i] = i;
    std::vector<uint8_t> inact(n, 0);  // multi golden has no comp restore path here
    bool need_tp = false; for (int e = 0; e < k; ++e) if (pb->type_ids[e] == 5) need_tp = true;

    MultiCtx C;
    C.n = n; C.k = k; C.D = D; C.a0 = a0.data();
    C.lam = lam.data(); C.th = th.data();
    C.type_ids = pb->type_ids; C.netof = pb->netof; C.netref = pb->netref; C.etys = pb->etypes;
    C.dyci = pb->dyc_idx; C.wtab = wtab.data(); C.p1 = p1.data(); C.p2 = p2.data();
    C.vmat = vmat.data(); C.dymat = dymat.data(); C.perm = perm.data(); C.symf = pb->symf;
    C.act = act.data(); C.nact = nact; C.inact = inact.data(); C.has_comp = 0;
    C.need_tp = need_tp;
    C.k0 = (uint32_t)(pb->seed & 0xFFFFFFFFu); C.k1 = (uint32_t)(pb->seed >> 32);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int c = 0; c < M; ++c) {
        std::vector<uint8_t> a((size_t)D*n*n);
        std::vector<double> sc(k), rs(D); double nt; long st;
        simulate_one_multi(C, (uint32_t)c, a.data(), sc.data(), rs.data(), nt, st);
        if (r->ntim) r->ntim[c] = nt;
        if (r->steps) r->steps[c] = (int32_t)st;
        if (r->rscores) for (int d = 0; d < D; ++d) r->rscores[(size_t)c*D+d] = (float)rs[d];
        if (r->scores) for (int e = 0; e < k; ++e) r->scores[(size_t)c*k+e] = (float)sc[e];
        if (r->a_ends) for (size_t q = 0; q < (size_t)D*n*n; ++q) r->a_ends[(size_t)c*D*n*n+q] = a[q];
    }
    return CUSNA_OK;
}

/* ---- A4.1: host evaluation statistics of one network (network_stats_ex) --- */
namespace {
inline double deg_stat(int t, double pm, double cur, double start) {
    if (t == 8 || t == 9) {  // outPop family
        if (pm >= 1.0) return t == 9 ? std::sqrt(cur) : cur;
        return pm == 0.0 ? start : cur + start;
    }                        // inAct family (14/15)
    if (pm >= 1.0) return t == 15 ? std::sqrt(cur) : cur;
    double s = pm == 0.0 ? start : cur + start;
    return t == 15 ? std::sqrt(s) : s;
}
}  // namespace

extern "C" CUSNA_API cusna_status
cusna_network_stats(const uint8_t *a, int32_t n, int32_t k,
                    const int32_t *type_ids, const double *p1, const double *p2,
                    const double *vmat, const uint8_t *vmiss, const int32_t *dyn_flags,
                    const double *vb, const int32_t *od0, const int32_t *id0,
                    const double *wtab, int32_t symflag, double *out) {
    if (!a || !type_ids || !out || n < 1 || n > CUSNA_MAX_N || k < 1) return CUSNA_ERR_NULL_ARG;
    std::vector<double> indeg(n, 0.0), outdeg(n, 0.0);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
        if (a[(size_t)i*n+j]) { outdeg[i] += 1.0; indeg[j] += 1.0; }

    bool need_aa = false, need_fb = false, need_bf = false, need_bb = false;
    for (int e = 0; e < k; ++e) { int t = type_ids[e];
        if (t==2||t==3||t==5||t==6||t==23||t==24||t==29) need_aa = true;
        if (t==26||t==30) need_fb = true; if (t==27) need_bf = true; if (t==28) need_bb = true; }
    std::vector<int> aa, fb, bf, bb;  // two-path / shared-partner counts
    auto mm = [&](std::vector<int> &M, bool rowT, bool colT) {
        M.assign((size_t)n*n, 0);
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { int s = 0;
            for (int h = 0; h < n; ++h) {
                int l = rowT ? a[(size_t)h*n+i] : a[(size_t)i*n+h];
                int r2 = colT ? a[(size_t)j*n+h] : a[(size_t)h*n+j];
                s += l & r2; }
            M[(size_t)i*n+j] = s; } };
    if (need_aa) mm(aa, false, false);  // a@a
    if (need_fb) mm(fb, false, true);   // a@a^T
    if (need_bf) mm(bf, true, false);   // a^T@a
    if (need_bb) mm(bb, true, true);    // a^T@a^T

    for (int e = 0; e < k; ++e) {
        int t = type_ids[e]; double acc = 0.0;
        const double *wr = wtab + (size_t)e * (n + 1);
        auto V = [&](int i){ return dyn_flags && dyn_flags[e] ? (vb?vb[i]:0.0) : (vmat?vmat[(size_t)e*n+i]:0.0); };
        auto miss = [&](int i){ return vmiss && vmiss[(size_t)e*n+i]; };
        double pm = p1 ? (double)p1[e] : 0.0, q2 = p2 ? (double)p2[e] : 0.0;
        switch (t) {
            case 0: for (int i=0;i<n;++i) acc += outdeg[i]; if (symflag) acc *= 0.5; break;
            case 1: for (int i=0;i<n;++i) for (int j=0;j<n;++j) acc += a[(size_t)i*n+j] && a[(size_t)j*n+i]; break;
            case 2: for (size_t q=0;q<(size_t)n*n;++q) if (a[q]) acc += aa[q]; break;
            case 3: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)j*n+i]) acc += aa[(size_t)i*n+j]; acc/=3.0; break;
            case 23:for (size_t q=0;q<(size_t)n*n;++q) if (a[q]) acc += aa[q]; acc/=6.0; break;
            case 5: for (size_t q=0;q<(size_t)n*n;++q) if (a[q] && aa[q]>0) acc += 1.0; break;
            case 6: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j] && a[(size_t)j*n+i]) acc += aa[(size_t)i*n+j]; break;
            case 4: for (int j=0;j<n;++j) acc += indeg[j]*indeg[j]; break;
            case 7: for (int j=0;j<n;++j) acc += std::pow(indeg[j],1.5); break;
            case 8: case 9: for (int i=0;i<n;++i) acc += indeg[i]*deg_stat(t,pm,outdeg[i],od0?od0[i]:0.0); break;
            case 14: case 15: for (int i=0;i<n;++i) acc += outdeg[i]*deg_stat(t,pm,indeg[i],id0?id0[i]:0.0); break;
            case 16: for (int i=0;i<n;++i) acc += outdeg[i]*outdeg[i]; break;
            case 17: for (int i=0;i<n;++i) acc += std::pow(outdeg[i],1.5); break;
            case 18: for (int i=0;i<n;++i) acc += std::fmin(outdeg[i],pm); break;
            case 19: for (int i=0;i<n;++i){ double rd=0; for (int j=0;j<n;++j) rd += a[(size_t)i*n+j] && a[(size_t)j*n+i]; acc += outdeg[i]*rd; } break;
            case 10: for (int i=0;i<n;++i) if (!miss(i)) acc += outdeg[i]*V(i); break;
            case 11: for (int i=0;i<n;++i) if (!miss(i)) acc += indeg[i]*V(i); break;
            case 12: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j] && !miss(i) && !miss(j)) acc += 1.0 - std::fabs(V(i)-V(j))*pm - q2; break;
            case 13: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j] && !miss(i) && !miss(j)) acc += std::fabs(V(i)-V(j))<1e-9?1.0:0.0; break;
            case 20: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j]) { double val; if (miss(i)||miss(j)) val=0.5; else { double d=V(i)-V(j); val = (d>1e-9?1.0:0.0)+0.5*(std::fabs(d)<=1e-9?1.0:0.0); } acc += val; } break;
            case 21: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j] && !miss(i) && !miss(j)) acc += V(i)*V(j); break;
            case 25: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j] && !miss(i) && !miss(j)) acc += V(i)+V(j); break;
            case 31: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j] && !miss(j) && V(j) >= pm - 1e-12) acc += 1.0; break;
            case 24: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j]) acc += wr[aa[(size_t)i*n+j]]; break;
            case 26: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j]) acc += wr[fb[(size_t)i*n+j]]; break;
            case 27: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j]) acc += wr[bf[(size_t)i*n+j]]; break;
            case 28: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (a[(size_t)i*n+j]) acc += wr[bb[(size_t)i*n+j]]; break;
            case 29: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (i!=j) acc += wr[aa[(size_t)i*n+j]]; break;
            case 30: for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (i!=j) acc += wr[fb[(size_t)i*n+j]]; break;
            default: break;  // X(22)/interaction(32) -> A4 later
        }
        out[e] = acc;
    }
    return CUSNA_OK;
}

extern "C" CUSNA_API cusna_status
cusna_typed_stats(const uint8_t *a_batch, int32_t n, int32_t M, int32_t etype,
                  const double *F, const uint8_t *x0,
                  const uint8_t *keep, const uint8_t *force,
                  const uint8_t *keepS, const uint8_t *forceS,
                  const uint8_t *miss_m, const uint8_t *miss_m1, double *out) {
    if (!a_batch || !F || !x0 || !out || n < 1 || M < 1) return CUSNA_ERR_NULL_ARG;
    const size_t nn = (size_t)n * n;  // F has a zero diagonal, so diagonal terms drop
    for (int c = 0; c < M; ++c) {
        const uint8_t *a = a_batch + (size_t)c * nn;
        double s = 0.0;
        for (size_t q = 0; q < nn; ++q) {
            if (etype == 1) {  // creation: gained = ((a&keep)|force) & ~x0 & ~miss_m
                int masked = (a[q] & (keep ? keep[q] : 1)) | (force ? force[q] : 0);
                if (masked && !x0[q] && !(miss_m && miss_m[q])) s += F[q];
            } else {           // endowment: lost = x0 & ~((a&keepS)|forceS) & ~miss_m1
                int cur = (a[q] & (keepS ? keepS[q] : 1)) | (forceS ? forceS[q] : 0);
                if (x0[q] && !cur && !(miss_m1 && miss_m1[q])) s += F[q];
            }
        }
        out[c] = (etype == 1) ? s : -s;   // network endowment reported negated
    }
    return CUSNA_OK;
}

extern "C" CUSNA_API cusna_status
cusna_gof_distribution(const uint8_t *a_batch, int32_t n, int32_t M, int32_t kind,
                       const int32_t *levls, int32_t nlev, double *out) {
    if (!a_batch || !levls || !out || n < 1 || M < 1 || nlev < 1) return CUSNA_ERR_NULL_ARG;
    const size_t nn = (size_t)n * n;
    if (kind == 0 || kind == 1) {  // in/out degree: cumulative actor counts
        std::vector<int> deg(n);
        for (int c = 0; c < M; ++c) {
            const uint8_t *a = a_batch + (size_t)c * nn;
            for (int v = 0; v < n; ++v) deg[v] = 0;
            for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
                if (a[(size_t)i*n+j]) deg[kind == 0 ? j : i]++;   // indeg / outdeg
            for (int L = 0; L < nlev; ++L) {
                long cnt = 0; for (int v = 0; v < n; ++v) if (deg[v] <= levls[L]) cnt++;
                out[(size_t)c * nlev + L] = (double)cnt;
            }
        }
        return CUSNA_OK;
    }
    if (kind == 2) {  // geodesic: ordered pairs at distance <= level (+ unreachable)
        int maxk = 0; for (int L = 0; L < nlev; ++L) if (levls[L] > maxk) maxk = levls[L];
        std::vector<uint8_t> cur(nn), nxt(nn);
        std::vector<long> counts(maxk + 1, 0);
        for (int c = 0; c < M; ++c) {
            const uint8_t *reach = a_batch + (size_t)c * nn;
            for (size_t q = 0; q < nn; ++q) cur[q] = reach[q] ? 1 : 0;
            for (int k = 1; k <= maxk; ++k) {
                long cn = 0;
                for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
                    if (i != j && cur[(size_t)i*n+j]) cn++;
                counts[k] = cn;
                if (k < maxk) {
                    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
                        uint8_t r = cur[(size_t)i*n+j];
                        if (!r) for (int m = 0; m < n; ++m)
                            if (cur[(size_t)i*n+m] && reach[(size_t)m*n+j]) { r = 1; break; }
                        nxt[(size_t)i*n+j] = r;
                    }
                    for (int i = 0; i < n; ++i) nxt[(size_t)i*n+i] = 0;
                    cur.swap(nxt);
                }
            }
            for (int L = 0; L < nlev; ++L)
                out[(size_t)c * (nlev+1) + L] = (double)counts[levls[L]];
            out[(size_t)c * (nlev+1) + nlev] = (double)((long)n * (n-1) - counts[maxk]);
        }
        return CUSNA_OK;
    }
    return CUSNA_ERR_UNSUPPORTED;
}

/* ---- A4.5: behavior evaluation statistics (behavior_stats, convention 9) --- */
extern "C" CUSNA_API cusna_status
cusna_behavior_stats(const int32_t *z_batch, int32_t n, int32_t M, int32_t kb,
                     const int32_t *bty, const uint8_t *xnet,
                     const int32_t *bkdist, double zbar, double brange,
                     double bsimmean, double *out) {
    if (!z_batch || !bty || !xnet || !out || n < 1 || n > CUSNA_MAX_N || M < 1 || kb < 1)
        return CUSNA_ERR_NULL_ARG;
    // ok[i]: actor i is observed at both period ends (kmask_dist); default all-1.
    // It filters egos AND alters, and (via stat_values) zeroes missing egos' u.
    std::vector<double> ok(n, 1.0);
    if (bkdist) for (int i = 0; i < n; ++i) ok[i] = bkdist[i] ? 1.0 : 0.0;
    // chain-invariant neighbour counts from the predictor network xnet (n*n):
    // outdeg = all out-alters (avAlt); nb = non-missing out-alters (avSim/totSim).
    std::vector<double> outdeg(n, 0.0), nb(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double od = 0.0, nbi = 0.0;
        for (int j = 0; j < n; ++j) if (xnet[(size_t)i*n+j]) { od += 1.0; nbi += ok[j]; }
        outdeg[i] = od; nb[i] = nbi;
    }
    std::vector<double> u(n);
    for (int c = 0; c < M; ++c) {
        const int32_t *z = z_batch + (size_t)c * n;
        for (int i = 0; i < n; ++i) u[i] = ok[i] ? ((double)z[i] - zbar) : 0.0;
        for (int e = 0; e < kb; ++e) {
            int t = bty[e]; double acc = 0.0;
            switch (t) {
                case 0: for (int i = 0; i < n; ++i) acc += u[i]; break;         // linear
                case 1: for (int i = 0; i < n; ++i) acc += u[i] * u[i]; break;  // quad
                case 2: for (int i = 0; i < n; ++i) {                           // avAlt
                    if (outdeg[i] <= 0.0) continue;   // avg 0; u_i=0 for missing egos
                    double as = 0.0;
                    for (int j = 0; j < n; ++j) if (xnet[(size_t)i*n+j]) as += u[j];
                    acc += u[i] * (as / outdeg[i]);
                } break;
                case 3: case 4: for (int i = 0; i < n; ++i) {                   // avSim/totSim
                    if (ok[i] == 0.0) continue;       // ego filter (missing egos drop)
                    double ss = 0.0;
                    for (int j = 0; j < n; ++j)
                        if (xnet[(size_t)i*n+j] && ok[j] != 0.0)
                            ss += 1.0 - std::fabs(u[i] - u[j]) / brange;
                    ss -= nb[i] * bsimmean;           // totSim -nb*simMean term (conv. 9)
                    if (t == 3) ss = nb[i] > 0.0 ? ss / nb[i] : 0.0;   // avSim divides
                    acc += ss;
                } break;
            }
            out[(size_t)c * kb + e] = acc;
        }
    }
    return CUSNA_OK;
}

/* ---- A4.5: cross-network evaluation statistic (cross_stat, multi-net) ------ */
extern "C" CUSNA_API cusna_status
cusna_cross_stats(const uint8_t *a_batch, const uint8_t *wp, int32_t n, int32_t M,
                  int32_t t, double parm, double *out) {
    if (!a_batch || !wp || !out || n < 1 || n > CUSNA_MAX_N || M < 1)
        return CUSNA_ERR_NULL_ARG;
    const size_t nn = (size_t)n * n;
    // `from`: star = wp @ wp^T (chain-invariant), optional sqrt for parm >= 2
    std::vector<double> star;
    if (t == 36) {
        star.assign(nn, 0.0);
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
            long s = 0;
            for (int h = 0; h < n; ++h) s += wp[(size_t)i*n+h] & wp[(size_t)j*n+h];
            star[(size_t)i*n+j] = (parm >= 2.0) ? std::sqrt((double)s) : (double)s;
        }
    }
    for (int c = 0; c < M; ++c) {
        const uint8_t *a = a_batch + (size_t)c * nn;
        double acc = 0.0;
        switch (t) {
            case 33: for (size_t q = 0; q < nn; ++q) acc += (a[q] & wp[q]); break;
            case 34: for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
                        acc += (a[(size_t)i*n+j] & wp[(size_t)j*n+i]); break;
            case 35: for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
                        acc += (a[(size_t)i*n+j] & wp[(size_t)i*n+j] & wp[(size_t)j*n+i]); break;
            case 36: for (size_t q = 0; q < nn; ++q) if (a[q]) acc += star[q]; break;
            case 37: for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
                        if (!a[(size_t)i*n+j]) continue;
                        long s = 0;  // (wp @ a)_ij = sum_h wp_ih * a_hj
                        for (int h = 0; h < n; ++h) s += wp[(size_t)i*n+h] & a[(size_t)h*n+j];
                        acc += (double)s;
                     } break;
            default: return CUSNA_ERR_UNSUPPORTED;
        }
        out[c] = acc;
    }
    return CUSNA_OK;
}

/* ---- CPU-only ABI utilities (A5) ------------------------------------------
 * In the GPU build these live in libcusna_saom.cu; a CPU-only build (the CRAN R
 * package Makevars defines CUSNA_CPU_ONLY) does not compile the .cu, so provide
 * the pure-C++ ones here. cusna_has_cuda is intentionally absent -- it needs the
 * CUDA runtime; the CPU-only glue reports no GPU without calling it. */
#ifdef CUSNA_CPU_ONLY
extern "C" CUSNA_API int32_t cusna_abi_version(void) { return CUSNA_ABI_VERSION; }
extern "C" CUSNA_API const char *cusna_version_string(void) {
    return "cusna-native 0.1 (track A, CPU-only build)";
}
extern "C" CUSNA_API const char *cusna_last_error(void) { return ""; }
#endif
