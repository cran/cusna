/* ==========================================================================
 *  libcusna_ergm_cpu.cpp -- host C++ for the native ERGM port (no CUDA).
 *
 *  Pure host code, so it compiles both in the native CMake build and in the CRAN
 *  CPU-only R package build (no nvcc). Milestones:
 *    E3: cusna_ergm_stats -- float64 port of tnt_gpu.py::network_stats_ergm (the
 *        7 ERGM terms), the host statistics carrier and the s0 provider for both
 *        samplers when problem->s0 is NULL.
 *    E4: cusna_ergm_tnt_sample_cpu -- the C++/OpenMP TNT-MH mirror (added later).
 *
 *  Summation order matches numpy's triu / row-major traversal where it is cheap
 *  to do so, so covariate/gwesp/gwdsp targets land within machine zero (~1e-13)
 *  of the Python reference; the integer terms (edges/mutual/nodematch) are exact.
 * ========================================================================== */
#include <cmath>
#include <cstring>
#include <vector>

#include "libcusna.h"

/* Term type ids -- must match _TERM_IDS in prototype/cusna/ergm/tnt_gpu.py and
 * the switch in ergm_kernels.cu. */
enum {
    ETERM_EDGES = 0, ETERM_MUTUAL = 1, ETERM_NODECOV = 2, ETERM_NODEMATCH = 3,
    ETERM_ABSDIFF = 4, ETERM_GWESP = 5, ETERM_GWDSP = 6
};

/* Host evaluation statistics of ONE network for the ERGM term list (E3).
 * Faithful float64 port of network_stats_ergm: a is dense n*n uint8 (row-major),
 * attr widened to double (NULL => zeros), parm[e] is the gwesp/gwdsp decay. */
extern "C" CUSNA_API cusna_status cusna_ergm_stats(
    const uint8_t *a, int32_t n, int32_t k, int32_t directed,
    const int32_t *type_ids, const float *parm, const float *attr, double *out) {
    if (!a || !type_ids || !out) return CUSNA_ERR_NULL_ARG;
    if (n < 1 || k < 1 || k > CUSNA_ERGM_MAXTERMS) return CUSNA_ERR_SHAPE;

    const double half = directed ? 1.0 : 0.5;
    std::vector<double> v((size_t)n, 0.0);
    if (attr) for (int i = 0; i < n; ++i) v[i] = (double)attr[i];

    for (int e = 0; e < k; ++e) {
        const int t = type_ids[e];
        double s = 0.0;
        switch (t) {
        case ETERM_EDGES: {                       /* a.sum() * half            */
            long long cnt = 0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) cnt += a[(size_t)i * n + j];
            s = (double)cnt * half;
            break;
        }
        case ETERM_MUTUAL: {                       /* (a * a.T).sum() / 2       */
            long long cnt = 0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    cnt += (long long)a[(size_t)i * n + j] * a[(size_t)j * n + i];
            s = (double)cnt / 2.0;
            break;
        }
        case ETERM_NODECOV: {                      /* sum a_ij*(v_i+v_j) * half */
            double acc = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (a[(size_t)i * n + j]) acc += v[i] + v[j];
            s = acc * half;
            break;
        }
        case ETERM_NODEMATCH: {                    /* sum a_ij*[v_i==v_j] * half*/
            long long cnt = 0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (a[(size_t)i * n + j] && std::fabs(v[i] - v[j]) < 1e-9) ++cnt;
            s = (double)cnt * half;
            break;
        }
        case ETERM_ABSDIFF: {                      /* sum a_ij*|v_i-v_j| * half */
            double acc = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (a[(size_t)i * n + j]) acc += std::fabs(v[i] - v[j]);
            s = acc * half;
            break;
        }
        case ETERM_GWESP:
        case ETERM_GWDSP: {                        /* undirected; w(s) weights  */
            if (directed) return CUSNA_ERR_UNSUPPORTED;  /* gwesp/gwdsp undirected */
            const double alpha = (double)parm[e];
            const double ea = std::exp(alpha);
            const double oma = 1.0 - std::exp(-alpha);
            double acc = 0.0;
            /* upper triangle i<j in row-major order (== numpy triu_indices(n,1)) */
            for (int i = 0; i < n; ++i) {
                const uint8_t *ai = a + (size_t)i * n;
                for (int j = i + 1; j < n; ++j) {
                    int cn = 0;                    /* common neighbours = (a@a)_ij */
                    for (int m = 0; m < n; ++m)
                        cn += ai[m] & a[(size_t)m * n + j];
                    if (t == ETERM_GWESP) {        /* edges only (mask a_ij>0)  */
                        if (ai[j]) acc += ea * (1.0 - std::pow(oma, (double)cn));
                    } else {                       /* gwdsp: all dyads, w(0)=0  */
                        acc += ea * (1.0 - std::pow(oma, (double)cn));
                    }
                }
            }
            s = acc;
            break;
        }
        default:
            return CUSNA_ERR_UNSUPPORTED;
        }
        out[e] = s;
    }
    return CUSNA_OK;
}

/* ==========================================================================
 *  E4: cusna_ergm_tnt_sample_cpu -- pure C++/OpenMP TNT mirror (no CUDA).
 *
 *  A faithful float64 transcription of the ergm_tnt kernel: host Philox4x32-10
 *  (same counter layout), the exact TNT proposal, the 7 change statistics, and
 *  the edge-list swap-remove -- one chain per OpenMP iteration (chains are
 *  independent, no races). This is the FIRST CPU ERGM sampler: there is no
 *  bit-exact reference (the Python sampler is GPU-only), so it is validated
 *  DISTRIBUTIONALLY against the GPU path and ergm::simulate. The integer change
 *  statistics (edges/mutual) can stay bit-close to the GPU; gwesp/gwdsp diverge
 *  where the CPU std::pow differs from the device pow by ~1 ULP near an accept
 *  boundary -- exactly the documented GPU-float vs CPU-float split.
 *
 *  Attribute arithmetic mirrors the kernel's float32 (fabsf / float add) so the
 *  covariate terms match to the bit for the integer-valued attributes used here.
 * ========================================================================== */
#include <cstdint>

static inline void ephilox4(uint32_t c0, uint32_t c1, uint32_t c2, uint32_t c3,
                            uint32_t k0, uint32_t k1, uint32_t r[4]) {
    for (int i = 0; i < 10; ++i) {
        uint64_t p0 = (uint64_t)0xD2511F53u * c0;
        uint64_t p1 = (uint64_t)0xCD9E8D57u * c2;
        uint32_t hi0 = (uint32_t)(p0 >> 32), lo0 = (uint32_t)p0;
        uint32_t hi1 = (uint32_t)(p1 >> 32), lo1 = (uint32_t)p1;
        c0 = hi1 ^ c1 ^ k0; c1 = lo1; c2 = hi0 ^ c3 ^ k1; c3 = lo0;
        k0 += 0x9E3779B9u; k1 += 0xBB67AE85u;
    }
    r[0] = c0; r[1] = c1; r[2] = c2; r[3] = c3;
}
static inline double eu01d(uint32_t x) {
    return ((double)x + 0.5) * 2.3283064365386963e-10;
}
static inline double egw_w(int s, double ea, double oma) {
    return ea * (1.0 - std::pow(oma, (double)s));
}

extern "C" CUSNA_API cusna_status
cusna_ergm_tnt_sample_cpu(const cusna_ergm_problem *p, cusna_ergm_result *r) {
    if (!p || !r) return CUSNA_ERR_NULL_ARG;
    if (p->struct_version != CUSNA_ABI_VERSION) return CUSNA_ERR_SHAPE;
    const int n = p->n, k = p->k, C = p->C, directed = p->directed ? 1 : 0;
    const int S = p->samplesize;
    if (n < 1 || n >= 65536 || k < 1 || k > CUSNA_ERGM_MAXTERMS || C < 1 || S < 1
        || p->interval < 1 || p->burnin < 0)
        return CUSNA_ERR_SHAPE;
    if (!p->x0 || !p->theta || !p->type_ids || !p->parm) return CUSNA_ERR_NULL_ARG;

    /* start statistics: caller-provided, or the host port */
    std::vector<double> s0buf;
    const double *s0 = p->s0;
    if (!s0) {
        s0buf.resize((size_t)k);
        cusna_status st = cusna_ergm_stats(p->x0, n, k, directed, p->type_ids,
                                           p->parm, p->attr, s0buf.data());
        if (st != CUSNA_OK) return st;
        s0 = s0buf.data();
    }

    const double ndyads = directed ? (double)n * (n - 1)
                                   : (double)n * (n - 1) / 2.0;
    const size_t maxe = directed ? (size_t)n * (n - 1) : (size_t)n * (n - 1) / 2;

    /* STERGM (X2): when p->permitted is set, restrict the TNT proposal to the
     * permitted (free) dyads. Enumerate them once (shared, read-only) and use
     * their count as the effective dyad universe; change statistics are still
     * evaluated on the full network. NULL => all dyads free (base path). */
    const bool masked = (p->permitted != nullptr);
    std::vector<uint32_t> fdyads;
    if (masked) {
        fdyads.reserve(maxe);
        if (directed) {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (i != j && p->permitted[(size_t)i * n + j])
                        fdyads.push_back(((uint32_t)i << 16) | (uint32_t)j);
        } else {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < i; ++j)
                    if (p->permitted[(size_t)i * n + j])
                        fdyads.push_back(((uint32_t)i << 16) | (uint32_t)j);
        }
        if (fdyads.empty()) {   /* no free dyads => constrained network fixed at x0 */
            int nx = 0;
            for (size_t q = 0; q < (size_t)n * n; ++q) nx += p->x0[q];
            if (!directed) nx /= 2;
            for (int c = 0; c < C; ++c) {
                if (r->stats)
                    for (int is = 0; is < S; ++is)
                        for (int e = 0; e < k; ++e)
                            r->stats[((size_t)c * S + is) * k + e] = s0[e];
                if (r->nacc) r->nacc[c] = 0;
                if (r->nedges) r->nedges[c] = nx;
                if (r->a_ends)
                    std::memcpy(r->a_ends + (size_t)c * n * n, p->x0, (size_t)n * n);
            }
            return CUSNA_OK;
        }
    }
    const double ndyads_eff = masked ? (double)fdyads.size() : ndyads;
    const double P = 0.5, Q = 0.5, DP = P * ndyads_eff, DO = DP / Q;

    std::vector<double> ea((size_t)k), oma((size_t)k);
    for (int e = 0; e < k; ++e) {
        ea[e] = std::exp((double)p->parm[e]);
        oma[e] = 1.0 - std::exp(-(double)p->parm[e]);
    }
    std::vector<double> v((size_t)n, 0.0);
    if (p->attr) for (int i = 0; i < n; ++i) v[i] = (double)p->attr[i];

    const uint32_t k0 = (uint32_t)(p->seed & 0xFFFFFFFFull);
    const uint32_t k1 = (uint32_t)(p->seed >> 32);
    const long long nsteps = (long long)p->burnin
        + (long long)p->interval * (S - 1);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int c = 0; c < C; ++c) {
        std::vector<uint8_t> a(p->x0, p->x0 + (size_t)n * n);   /* dense start */
        std::vector<uint32_t> elist(maxe ? maxe : 1, 0u);
        int nedges = 0;
        if (directed) {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (i != j && a[(size_t)i * n + j]
                        && (!masked || p->permitted[(size_t)i * n + j]))
                        elist[nedges++] = ((uint32_t)i << 16) | (uint32_t)j;
        } else {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < i; ++j)
                    if (a[(size_t)i * n + j]
                        && (!masked || p->permitted[(size_t)i * n + j]))
                        elist[nedges++] = ((uint32_t)i << 16) | (uint32_t)j;
        }

        double stats[CUSNA_ERGM_MAXTERMS];
        for (int e = 0; e < k; ++e) stats[e] = s0[e];
        int nacc = 0, isample = 0;

        for (long long step = 0; step <= nsteps; ++step) {
            if (step >= p->burnin && isample < S
                && (step - p->burnin) % p->interval == 0) {
                if (r->stats)
                    for (int e = 0; e < k; ++e)
                        r->stats[((size_t)c * S + isample) * k + e] = stats[e];
                ++isample;
            }
            if (step == nsteps) break;

            uint32_t rr[4];
            ephilox4((uint32_t)(step & 0xFFFFFFFF), (uint32_t)c,
                     (uint32_t)(step >> 32), 0u, k0, k1, rr);
            int t, h;
            bool edgestate;
            double logratio;
            if (eu01d(rr[0]) < P && nedges > 0) {
                int ei = (int)(eu01d(rr[1]) * nedges);
                if (ei >= nedges) ei = nedges - 1;
                const uint32_t pk = elist[ei];
                t = (int)(pk >> 16);
                h = (int)(pk & 0xFFFF);
                edgestate = true;
                logratio = (nedges == 1) ? std::log(1.0 / (DP + Q))
                                         : std::log((double)nedges / (DO + nedges));
            } else {
                const double u = eu01d(rr[1]);
                if (masked) {                        /* pick a uniform free dyad */
                    long long idx = (long long)(u * ndyads_eff);
                    if (idx >= (long long)ndyads_eff) idx = (long long)ndyads_eff - 1;
                    const uint32_t pk = fdyads[(size_t)idx];
                    t = (int)(pk >> 16);
                    h = (int)(pk & 0xFFFF);
                } else if (directed) {
                    long long idx = (long long)(u * ndyads);
                    if (idx >= (long long)ndyads) idx = (long long)ndyads - 1;
                    t = (int)(idx / (n - 1));
                    int rem = (int)(idx % (n - 1));
                    h = rem + (rem >= t ? 1 : 0);
                } else {
                    long long idx = (long long)(u * ndyads);
                    if (idx >= (long long)ndyads) idx = (long long)ndyads - 1;
                    int lo = 1, hi = n - 1;
                    while (lo < hi) {
                        int mid = (lo + hi + 1) >> 1;
                        if ((long long)mid * (mid - 1) / 2 <= idx) lo = mid;
                        else hi = mid - 1;
                    }
                    t = lo;
                    h = (int)(idx - (long long)t * (t - 1) / 2);
                }
                edgestate = a[(size_t)t * n + h] != 0;
                if (edgestate)
                    logratio = (nedges == 1) ? std::log(1.0 / (DP + Q))
                                             : std::log((double)nedges / (DO + nedges));
                else
                    logratio = (nedges == 0) ? std::log(DP + Q)
                                             : std::log(1.0 + DO / (nedges + 1.0));
            }

            const double delta = edgestate ? -1.0 : 1.0;
            double ds[CUSNA_ERGM_MAXTERMS];
            double ip = 0.0;
            for (int e = 0; e < k; ++e) {
                double d = 0.0;
                switch (p->type_ids[e]) {
                case 0: d = delta; break;                        /* edges       */
                case 1:                                          /* mutual      */
                    d = delta * (double)(a[(size_t)h * n + t] != 0);
                    break;
                case 2: {                                        /* nodecov     */
                    float sum = (float)v[t] + (float)v[h];       /* kernel float */
                    d = delta * (double)sum;
                    break;
                }
                case 3: {                                        /* nodematch   */
                    float diff = std::fabs((float)v[t] - (float)v[h]);
                    d = delta * ((diff < 1e-9f) ? 1.0 : 0.0);
                    break;
                }
                case 4: {                                        /* absdiff     */
                    float diff = std::fabs((float)v[t] - (float)v[h]);
                    d = delta * (double)diff;
                    break;
                }
                case 5: {                                        /* gwesp       */
                    int cn = 0;
                    double dd = 0.0;
                    for (int kk = 0; kk < n; ++kk)
                        if (a[(size_t)t * n + kk] && a[(size_t)h * n + kk]) {
                            ++cn;
                            int e1 = 0, e2 = 0;
                            for (int m = 0; m < n; ++m) {
                                e1 += a[(size_t)t * n + m] & a[(size_t)kk * n + m];
                                e2 += a[(size_t)h * n + m] & a[(size_t)kk * n + m];
                            }
                            dd += egw_w(e1 + (int)delta, ea[e], oma[e]) - egw_w(e1, ea[e], oma[e]);
                            dd += egw_w(e2 + (int)delta, ea[e], oma[e]) - egw_w(e2, ea[e], oma[e]);
                        }
                    d = delta * egw_w(cn, ea[e], oma[e]) + dd;
                    break;
                }
                case 6: {                                        /* gwdsp       */
                    double dd = 0.0;
                    for (int kk = 0; kk < n; ++kk)               /* k in N(h): (t,k) */
                        if (a[(size_t)h * n + kk] && kk != t) {
                            int dsp = 0;
                            for (int m = 0; m < n; ++m)
                                dsp += a[(size_t)t * n + m] & a[(size_t)kk * n + m];
                            dd += egw_w(dsp + (int)delta, ea[e], oma[e]) - egw_w(dsp, ea[e], oma[e]);
                        }
                    for (int kk = 0; kk < n; ++kk)               /* k in N(t): (h,k) */
                        if (a[(size_t)t * n + kk] && kk != h) {
                            int dsp = 0;
                            for (int m = 0; m < n; ++m)
                                dsp += a[(size_t)h * n + m] & a[(size_t)kk * n + m];
                            dd += egw_w(dsp + (int)delta, ea[e], oma[e]) - egw_w(dsp, ea[e], oma[e]);
                        }
                    d = dd;
                    break;
                }
                case 7:                                          /* memory (TERGM) */
                    d = delta * (double)(p->lag &&
                        p->lag[(size_t)t * n + h]);              /* lagged tie t->h */
                    break;
                case 8:                                          /* delrecip (TERGM) */
                    d = delta * (double)(p->lag &&
                        p->lag[(size_t)h * n + t]);              /* lagged tie h->t */
                    break;
                default: break;
                }
                ds[e] = d;
                ip += (double)p->theta[e] * d;
            }

            if (std::log(eu01d(rr[3])) < ip + logratio) {
                a[(size_t)t * n + h] ^= 1u;
                if (!directed) a[(size_t)h * n + t] ^= 1u;
                if (edgestate) {
                    const uint32_t pk = ((uint32_t)t << 16) | (uint32_t)h;
                    const uint32_t pk2 = ((uint32_t)h << 16) | (uint32_t)t;
                    for (int ei = 0; ei < nedges; ++ei)
                        if (elist[ei] == pk || (!directed && elist[ei] == pk2)) {
                            elist[ei] = elist[nedges - 1];
                            break;
                        }
                    --nedges;
                } else {
                    elist[nedges++] = ((uint32_t)t << 16) | (uint32_t)h;
                }
                for (int e = 0; e < k; ++e) stats[e] += ds[e];
                ++nacc;
            }
        }
        if (r->nacc) r->nacc[c] = nacc;
        if (r->nedges) r->nedges[c] = nedges;
        if (r->a_ends)
            std::memcpy(r->a_ends + (size_t)c * n * n, a.data(), (size_t)n * n);
    }
    return CUSNA_OK;
}

/* ==========================================================================
 *  E5: cusna_ergm_mple_stats -- MPLE change-statistic design matrix (host).
 *
 *  MPLE reduces ERGM estimation to a logistic regression of y_ij on the 0->1
 *  change statistics delta_ij (Strauss & Ikeda 1990); building that design matrix
 *  is the parallel part, the GLM fit itself is standard (R glm() / numpy IRLS,
 *  left to the caller). This is the DEMO from prototype/cusna/ergm/mple_gpu.py --
 *  three DIRECTED terms (edges, mutual, ttriple) -- honestly the M-3 boundary:
 *  cusna's ERGM is a sampler + this MPLE demo, NOT a general MCMC-MLE estimator.
 *  Not extended to the full 7-term catalog (no reference to validate against; the
 *  sampler already carries those terms).
 *
 *  For every off-diagonal dyad (i,j) in row-major order it writes one design row
 *      [ edges=1, mutual=x_ji, ttriple=|out(i)&in(j)|+|out(i)&out(j)|+|in(i)&in(j)| ]
 *  into X (nrow x 3, row-major) and the response x_ij into y. Integer change
 *  statistics -> exact vs mple_change_stats_numpy (< 2^24, exact in float32).
 *  Caller allocates X[n*n*3] and y[n*n]; the diagonal is skipped, so *nrow =
 *  n*(n-1), *ncol = 3.
 * ========================================================================== */
extern "C" CUSNA_API cusna_status cusna_ergm_mple_stats(
    const uint8_t *a, int32_t n, float *X, uint8_t *y,
    int32_t *nrow, int32_t *ncol) {
    if (!a || !X || !y) return CUSNA_ERR_NULL_ARG;
    if (n < 2 || n >= 65536) return CUSNA_ERR_SHAPE;

    size_t row = 0;
    for (int i = 0; i < n; ++i) {
        const uint8_t *oi = a + (size_t)i * n;                 /* out-ties of i */
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            const uint8_t *oj = a + (size_t)j * n;             /* out-ties of j */
            int tp = 0;   /* two-paths / transitive-triple roles closed by i->j */
            for (int m = 0; m < n; ++m) {
                tp += oi[m] & a[(size_t)m * n + j];            /* out(i) & in(j) */
                tp += oi[m] & oj[m];                           /* out(i) & out(j) */
                tp += a[(size_t)m * n + i] & a[(size_t)m * n + j]; /* in(i)&in(j) */
            }
            X[row * 3 + 0] = 1.0f;
            X[row * 3 + 1] = (float)a[(size_t)j * n + i];      /* mutual = x_ji  */
            X[row * 3 + 2] = (float)tp;                        /* ttriple        */
            y[row] = a[(size_t)i * n + j];                     /* response x_ij  */
            ++row;
        }
    }
    if (nrow) *nrow = (int32_t)row;
    if (ncol) *ncol = 3;
    return CUSNA_OK;
}

/* ==========================================================================
 *  X1.5: ALAAM (autologistic actor attribute model) -- attribute-DV extension.
 *
 *  ALAAM models a BINARY actor attribute y at a FIXED (exogenous) network A: the
 *  dependent variable is the attribute, not a tie. It is the SAME bit-kernel
 *  primitive applied to (network-row AND attribute-vector) instead of (row AND
 *  row) -- the change statistic of the contagion configuration when actor i flips
 *  its attribute is the masked popcount of A's row i against the attribute vector.
 *  Reference-free on CRAN (ALAAMEE is a GitHub/Python package); validated by the
 *  exact equivalence of the MPLE to a plain logistic regression on the native
 *  change statistics (see native/test/compare_alaam.R).
 *
 *  cusna_alaam_stats -- per-actor contagion change statistics:
 *      con_out[i] = sum_j A[i,j] * y[j]   (out-neighbours of i with y=1)
 *      con_in[i]  = sum_j A[j,i] * y[j]   (in-neighbours of i with y=1)
 *  For an undirected (symmetric) A the two coincide and give the neighbour count
 *  with y=1 -- exactly the ALAAM contagion change statistic. Integer-exact host.
 * ========================================================================== */
extern "C" CUSNA_API cusna_status cusna_alaam_stats(
    const uint8_t *a, const uint8_t *y, int32_t n,
    double *con_out, double *con_in) {
    if (!a || !y || !con_out || !con_in) return CUSNA_ERR_NULL_ARG;
    if (n < 1 || n >= 65536) return CUSNA_ERR_SHAPE;
    for (int i = 0; i < n; ++i) {
        const uint8_t *ai = a + (size_t)i * n;          /* out-ties of i        */
        long long co = 0, ci = 0;
        for (int j = 0; j < n; ++j) {
            co += (long long)ai[j] * y[j];              /* out-neighbour j, y_j */
            ci += (long long)a[(size_t)j * n + i] * y[j]; /* in-neighbour j, y_j */
        }
        con_out[i] = (double)co;
        con_in[i]  = (double)ci;
    }
    return CUSNA_OK;
}

/* ALAAM Gibbs sampler over the binary attribute at a FIXED network (X1.5). The
 * full conditional of the autologistic model is
 *     logit P(y_i = 1 | y_{-i}) = beta0 + beta_c * c_i + sum_p beta_x[p]*X[i,p],
 * with the contagion c_i = sum_j A[i,j]*y_j (+ sum_j A[j,i]*y_j when directed) --
 * the same masked popcount the MPLE regresses on, so estimator and simulator use
 * one statistic. One chain per OpenMP iteration; each runs `sweeps` full Gibbs
 * sweeps (actors in index order) from the start attribute y0, host Philox RNG
 * (independent per chain). Writes each chain's end attribute into y_ends [C*n]
 * and, if stats != NULL, [density = sum_i y_i, contagion = sum over edges y_i*y_j]
 * per chain into stats [C*2]. X is [n*p] row-major (NULL/p=0 for none). */
extern "C" CUSNA_API cusna_status cusna_alaam_gibbs_cpu(
    const uint8_t *a, int32_t n, int32_t directed,
    const uint8_t *y0, const float *X, int32_t p,
    double beta0, double beta_c, const float *beta_x,
    int32_t C, uint64_t seed, int32_t sweeps,
    uint8_t *y_ends, double *stats) {
    if (!a || !y0 || !y_ends) return CUSNA_ERR_NULL_ARG;
    if (n < 1 || n >= 65536 || C < 1 || sweeps < 0) return CUSNA_ERR_SHAPE;
    if (p < 0 || (p > 0 && (!X || !beta_x))) return CUSNA_ERR_NULL_ARG;

    const uint32_t k0 = (uint32_t)(seed & 0xFFFFFFFFull);
    const uint32_t k1 = (uint32_t)(seed >> 32);
    const double halfc = directed ? 1.0 : 0.5;   /* undirected edge counted once */

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int c = 0; c < C; ++c) {
        std::vector<uint8_t> y(y0, y0 + (size_t)n);
        /* fixed per-actor covariate contribution beta_x . X[i,:] */
        std::vector<double> xb((size_t)n, 0.0);
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int q = 0; q < p; ++q) s += (double)beta_x[q] * (double)X[(size_t)i * p + q];
            xb[i] = s;
        }
        for (int sw = 0; sw < sweeps; ++sw) {
            for (int i = 0; i < n; ++i) {
                const uint8_t *ai = a + (size_t)i * n;
                long long ci = 0;
                for (int j = 0; j < n; ++j) ci += (long long)ai[j] * y[j];
                if (directed)
                    for (int j = 0; j < n; ++j) ci += (long long)a[(size_t)j * n + i] * y[j];
                const double lo = beta0 + beta_c * (double)ci + xb[i];
                const double prob = 1.0 / (1.0 + std::exp(-lo));
                uint32_t rr[4];
                long long ctr = (long long)sw * n + i;
                ephilox4((uint32_t)(ctr & 0xFFFFFFFF), (uint32_t)c,
                         (uint32_t)(ctr >> 32), 0u, k0, k1, rr);
                y[i] = (eu01d(rr[0]) < prob) ? 1u : 0u;
            }
        }
        std::memcpy(y_ends + (size_t)c * n, y.data(), (size_t)n);
        if (stats) {
            long long dens = 0;
            double con = 0.0;
            for (int i = 0; i < n; ++i) {
                dens += y[i];
                if (y[i]) {
                    const uint8_t *ai = a + (size_t)i * n;
                    long long ci = 0;
                    for (int j = 0; j < n; ++j) ci += (long long)ai[j] * y[j];
                    con += (double)ci;
                }
            }
            stats[(size_t)c * 2 + 0] = (double)dens;
            stats[(size_t)c * 2 + 1] = con * halfc;   /* sum over edges y_i*y_j */
        }
    }
    return CUSNA_OK;
}
