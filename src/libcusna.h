/* ==========================================================================
 *  libcusna.h -- C ABI for the native cusna SAOM/ERGM engine (track A).
 *
 *  Goal: a self-contained native library (CUDA kernels + C++ host logic) that
 *  R (via cpp11/Rcpp) and other callers can drive WITHOUT a Python runtime,
 *  modelled on the R `torch` package (CRAN build is CPU-only; the GPU path is
 *  compiled from source when a CUDA toolkit is present).
 *
 *  The ABI is intentionally split into two layers so the port can be proven in
 *  stages (see native/DESIGN.md):
 *
 *    Layer 1 -- kernel ABI (milestones A1, A2, A3).
 *      cusna_saom_simulate_period() runs M chains of ONE SAOM period from
 *      ALREADY-marshalled inputs (packed bit rows, effect descriptor arrays,
 *      RSiena masks). It is a thin wrapper over the CUDA kernel saom_forward
 *      (native/saom_kernels.cu, extracted verbatim from the validated NVRTC
 *      source). A1 proves bit-identity of the baseline slice (density/recip/
 *      transTrip) against native/testdata/s50_baseline_golden.npz; A2 exercises
 *      the rest of the same kernel; A3 adds the CPU mirror behind the SAME
 *      struct.
 *
 *    Layer 2 -- estimator ABI (milestones A4, A5). [DECLARED, NOT YET DEFINED]
 *      cusna_mom_estimate() will own the RSiena data conventions (Appendix B),
 *      the Method-of-Moments Robbins-Monro loop and GoF -- the ~2000 host lines
 *      that are the dominant cost and risk of the port. Sketched at the bottom
 *      of this header; its concrete shape is fixed in A4.
 *
 *  Conventions:
 *    - extern "C", fixed-width integer types, plain-old-data structs, caller
 *      allocates all output buffers (NULL to skip an optional output).
 *    - structs carry a `struct_version` first field so fields can be appended
 *      without breaking old callers (check cusna_abi_version()).
 *    - all matrices are row-major; "packed bit rows" are uint32 words, W words
 *      per actor, bit (j & 31) of word (j >> 5) = tie to actor j (little-endian
 *      within the word, matching numpy packbits(bitorder="little")).
 *
 *  This header is hand-written and stable; saom_kernels.cu is generated. Keep
 *  MAXEFF/MAXBEFF/MAXRATE in sync with the kernel #defines.
 * ========================================================================== */
#ifndef LIBCUSNA_H
#define LIBCUSNA_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#  ifdef LIBCUSNA_BUILD
#    define CUSNA_API __declspec(dllexport)
#  else
#    define CUSNA_API __declspec(dllimport)
#  endif
#else
#  define CUSNA_API __attribute__((visibility("default")))
#endif

/* Kernel capacity limits -- MUST match the #defines in saom_kernels.cu. */
#define CUSNA_MAXEFF  24  /* network evaluation effects per model            */
#define CUSNA_MAXBEFF 6   /* behavior effects per model                      */
#define CUSNA_MAXRATE 4   /* non-basic rate effects per network              */
#define CUSNA_MAXNET  4   /* co-evolving networks (multi-net kernel)         */
#define CUSNA_MAX_N   8192

/* ABI/library version. cusna_abi_version() returns CUSNA_ABI_VERSION at build
 * time; bump on any breaking struct change. */
#define CUSNA_ABI_VERSION 1

/* Status codes returned by every entry point. */
typedef enum cusna_status {
    CUSNA_OK              = 0,
    CUSNA_ERR_NULL_ARG    = 1,  /* required pointer was NULL                  */
    CUSNA_ERR_SHAPE       = 2,  /* n/W/k out of range or inconsistent         */
    CUSNA_ERR_NO_CUDA     = 3,  /* GPU path requested, no CUDA device         */
    CUSNA_ERR_CUDA        = 4,  /* a CUDA runtime call failed                 */
    CUSNA_ERR_UNSUPPORTED = 5,  /* option not implemented at this milestone   */
    CUSNA_ERR_INTERNAL    = 6
} cusna_status;

/* --------------------------------------------------------------------------
 *  Layer 1: SAOM period-simulation kernel ABI (A1/A2/A3).
 *
 *  cusna_saom_problem mirrors the argument list of the saom_forward kernel.
 *  Shapes below use: n actors, W = ceil(n/32) words, k effects, M chains.
 *  Pointers marked [opt] may be NULL when the corresponding feature is off
 *  (kb==0, has_comp==0, nrate==0); the A1 baseline leaves all of them NULL.
 * -------------------------------------------------------------------------- */
typedef struct cusna_saom_problem {
    int32_t struct_version;   /* set to CUSNA_ABI_VERSION                     */

    /* dimensions */
    int32_t n;                /* actors (<= CUSNA_MAX_N)                      */
    int32_t W;                /* uint32 words per actor = (n + 31) / 32       */
    int32_t k;                /* network effects (<= CUSNA_MAXEFF)            */
    int32_t M;                /* number of chains (CUDA blocks)              */
    int32_t threads;          /* block size; 0 = auto (64/128/256 by n)      */

    /* start network, packed bit rows, length n*W each (row-major over actors)*/
    const uint32_t *x0_out;   /* out-ties of the data start network          */
    const uint32_t *x0_in;    /* in-ties  of the data start network          */
    /* simulation start network (== x0 unless composition change clears ties) */
    const uint32_t *sim_out;  /* [opt] NULL => use x0_out                     */
    const uint32_t *sim_in;   /* [opt] NULL => use x0_in                      */

    /* rates and evaluation parameters (per chain) */
    const float *lam;         /* length M: basic rate per chain              */
    const float *theta;       /* length M*k, row-major: eval params per chain*/

    /* effect descriptors, length k unless noted */
    const int32_t *type_ids;  /* effect type id (see saom_kernels.cu ego_ds) */
    const int32_t *etypes;    /* 0 eval, 1 creation, 2 endowment             */
    const float   *wtab;      /* k*(n+1): gwesp/gwdsp weight tables           */
    const float   *p1;        /* length k: effect parameter 1                */
    const float   *p2;        /* length k: effect parameter 2                */
    const float   *vmat;      /* k*n: per-effect covariate values            */
    const uint8_t *vmiss;     /* k*n: covariate missing flags (statistics)   */
    const int32_t *dyc_idx;   /* length k: index into dyc, or -1             */
    const float   *dyc;       /* [opt] nd*n*n dyadic covariates              */
    const int32_t *int0;      /* length k: interaction component 0 or -1     */
    const int32_t *int1;      /* length k: interaction component 1 or -1     */
    const int32_t *int2;      /* length k: interaction component 2 or -1     */
    const int32_t *dyn_flags; /* length k: effect covariate = behavior (0/1) */

    /* RSiena masks, packed bit rows length n*W each */
    const uint32_t *kdist;    /* dyads counted in the period distance        */
    const uint32_t *keep_out; /* statistic network: (A & keep) | force       */
    const uint32_t *keep_in;
    const uint32_t *force_out;
    const uint32_t *force_in;
    const uint32_t *perm;     /* toggleable (permitted) dyads                */

    /* start degrees, length n (for luseStart degree statistics) */
    const int32_t *od0;
    const int32_t *id0;

    /* scalar controls */
    int32_t maxdegree;        /* > 0 forbids introductions at the bound      */
    int32_t need_mask;        /* precomputed count-need bits (choice loop)   */
    int32_t conditional;      /* 1 = stop at target_dist, internal rate 1    */
    int32_t target_dist;      /* conditional stopping distance               */
    int32_t symflag;          /* symmetric network: paired cell flips        */
    uint64_t seed;            /* Philox key (64-bit)                         */

    /* composition change [opt; has_comp==0 => everybody active] */
    int32_t has_comp;         /* run the setLeaversBack restore pass         */
    int32_t nact;             /* number of active actors                     */
    const int32_t  *act;      /* [opt] length n: active actor ids first nact */
    const uint32_t *inactw;   /* [opt] length W: inactive-actor bitmask      */

    /* dependent behavior (co-evolution) [opt; kb==0 => pure network model] */
    int32_t kb;               /* behavior effects (<= CUSNA_MAXBEFF)         */
    const float   *lamb;      /* [opt] length M: behavior rate per chain     */
    const float   *thetab;    /* [opt] length M*kb: behavior params          */
    const int32_t *bty;       /* [opt] length kb: behavior effect type ids   */
    const int32_t *z0;        /* [opt] length n: start behavior values       */
    const float   *vb;        /* [opt] length n: centered start behavior     */
    const int32_t *bkdist;    /* [opt] length n: actor in behavior distance  */
    int32_t zmin, zmax;
    float   zbar, brange, bsimmean;

    /* non-basic rate effects [opt; nrate==0 => constant rate] */
    int32_t nrate;            /* outRate/inRate/RateX effects (<= MAXRATE)    */
    const int32_t *rkind;     /* [opt] length nrate: 0 RateX,1 outRate,2 inRate*/
    const float   *ralpha;    /* [opt] length nrate: rate-effect parameters  */
    const float   *ru;        /* [opt] nrate*n: covariate values (kind 0)    */

    /* float64 host-statistic descriptors (A4.3) [opt]. Host statistics live in
     * float64; when these are non-NULL the CPU backend uses them (exact) for
     * both the choice dynamics and the eval statistics, instead of promoting the
     * float32 kernel descriptors -- so covariate/gwesp results match the
     * float64 reference to machine zero. Same shapes as the float32 fields;
     * ignored by the GPU path. */
    const double *vmat64;     /* [opt] k*n covariate values                   */
    const double *wtab64;     /* [opt] k*(n+1) gwesp/gwdsp weight tables      */
    const double *vb64;       /* [opt] n centered start behavior              */
    const double *p1_64;      /* [opt] k effect parameter 1                   */
    const double *p2_64;      /* [opt] k effect parameter 2                   */
} cusna_saom_problem;

/* Output buffers, caller-allocated. NULL to skip an optional output. Types
 * match the kernel's native output types (no widening), so a bit-identity
 * test compares these directly against the frozen golden. */
typedef struct cusna_saom_result {
    int32_t struct_version;

    double  *stats;      /* [M*k]      final masked statistics (double)      */
    double  *ntim;       /* [M]        elapsed simulated time                */
    float   *scores;     /* [M*(k+1)]  col 0 = rate score, then k effects    */
    int32_t *dist;       /* [M]        signed simulated distance vs x0       */
    int32_t *dist_stat;  /* [M]        |x_end xor x0| over kdist             */
    int32_t *steps;      /* [M]        ministeps taken                       */

    /* behavior outputs [opt; only when problem->kb > 0] */
    int32_t *bdist;      /* [opt] [M]                                        */
    float   *bscores;    /* [opt] [M*(kb+1)]                                 */
    int32_t *z_end;      /* [opt] [M*n]                                      */

    /* rate-effect outputs [opt; only when problem->nrate > 0] */
    float   *rscores;    /* [opt] [M*nrate]                                  */
    int32_t *actor_dist; /* [opt] [M*n] per-ego masked tie differences       */

    /* end networks [opt] -- dense uint8 n*n per chain, row-major, for GoF/
     * typed-effect statistics. Set want_end_nets and provide a_ends[M*n*n]. */
    int32_t  want_end_nets;
    uint8_t *a_ends;     /* [opt] [M*n*n]                                    */
} cusna_saom_result;

/* Run M chains of one SAOM period on the GPU. Returns CUSNA_ERR_NO_CUDA if no
 * device is present (use cusna_has_cuda() to dispatch, or the _cpu variant). */
CUSNA_API cusna_status cusna_saom_simulate_period(
    const cusna_saom_problem *problem, cusna_saom_result *result);

/* CPU mirror (numba path port, milestone A3). Same struct; bit-identical to
 * the GPU path where the kernel guarantees it. Until A3: CUSNA_ERR_UNSUPPORTED.*/
CUSNA_API cusna_status cusna_saom_simulate_period_cpu(
    const cusna_saom_problem *problem, cusna_saom_result *result);

/* --------------------------------------------------------------------------
 *  Layer 1b: multi-network co-evolution kernel ABI (A4.5).
 *
 *  cusna_saom_multi_problem mirrors the saom_multi kernel (native/saom_kernels.cu
 *  extracted verbatim; the "hot loop only" launcher, transcribed from
 *  gpu.py::SaomGpuSimulator.simulate_multi). D co-evolving one-mode networks
 *  (2..CUSNA_MAXNET) share one actor set; the kernel picks a network per
 *  ministep in proportion to its rate, toggles one tie, and accumulates the
 *  per-net rate scores + per-effect scores. dist_stat and the masked EVALUATION
 *  statistics (including the cross-network effects) are computed on the HOST from
 *  the returned end networks -- the same split as the single-net path.
 *
 *  Effect blocks are per network: `netof[e]` = which network effect e belongs
 *  to; cross-network effects (type id >= 33: crprod/crprodRecip/crprodMutual/
 *  from/to) read the CURRENT state of source network `netref[e]` in the dynamics
 *  (and the period-START predictor in the host statistics). Behavior,
 *  interactions, and transTies/transRecTrip/reciAct are not supported here.
 * -------------------------------------------------------------------------- */
typedef struct cusna_saom_multi_problem {
    int32_t struct_version;   /* set to CUSNA_ABI_VERSION                     */
    int32_t n;                /* actors (shared across the D networks)        */
    int32_t W;                /* uint32 words per actor = (n + 31) / 32       */
    int32_t k;                /* total network effects across nets (<=MAXEFF) */
    int32_t M;                /* number of chains                             */
    int32_t D;                /* co-evolving networks (2..CUSNA_MAXNET)       */
    int32_t threads;          /* block size; 0 = auto (64/128/256 by n)       */

    const uint32_t *x0_out;   /* D*n*W packed out-ties of the D start nets    */
    const uint32_t *x0_in;    /* D*n*W packed in-ties                         */
    const float   *lam;       /* M*D  per-chain per-net basic rate            */
    const float   *theta;     /* M*k  eval params per chain (per-net blocks)  */

    const int32_t *type_ids;  /* k: effect type id                           */
    const int32_t *etypes;    /* k: 0 eval, 1 creation, 2 endowment          */
    const int32_t *netof;     /* k: network each effect belongs to (0..D-1)  */
    const int32_t *netref;    /* k: source network for cross effects, else 0 */
    const float   *wtab;      /* k*(n+1): gwesp/gwdsp weight tables           */
    const float   *p1;        /* k                                           */
    const float   *p2;        /* k                                           */
    const float   *vmat;      /* k*n: per-effect covariate values            */
    const float   *dyc;       /* [opt] nd*n*n dyadic covariates              */
    const int32_t *dyc_idx;   /* k: index into dyc, or -1                    */
    const uint32_t *perm;     /* D*n*W: per-net toggleable (permitted) dyads */
    const int32_t *symf;      /* D: per-net symmetric flag (paired flips)    */
    uint64_t seed;            /* Philox key (64-bit)                         */
    const int32_t *act;       /* [opt] n: active actor ids first nact        */
    int32_t nact;             /* active actor count (== n without comp)      */
} cusna_saom_multi_problem;

/* Outputs, caller-allocated, NULL to skip. The kernel writes rate/effect scores
 * with the MAXNET/MAXEFF stride internally; the launcher slices them to D/k. */
typedef struct cusna_saom_multi_result {
    int32_t struct_version;
    float   *rscores;    /* [opt] [M*D]     per-net basic-rate scores        */
    float   *scores;     /* [opt] [M*k]     effect scores                    */
    double  *ntim;       /* [opt] [M]       elapsed simulated time           */
    int32_t *steps;      /* [opt] [M]       ministeps taken                  */
    uint8_t *a_ends;     /* [opt] [M*D*n*n] dense end networks (row-major)   */
} cusna_saom_multi_result;

/* Run M chains of one multi-network SAOM period on the GPU (saom_multi). Returns
 * CUSNA_ERR_NO_CUDA if no device is present. dist_stat and the masked eval /
 * cross statistics are the caller's job (cusna_network_stats + cusna_cross_stats
 * over a_ends), matching the Python simulate_multi host split. */
CUSNA_API cusna_status cusna_saom_simulate_multi(
    const cusna_saom_multi_problem *problem, cusna_saom_multi_result *result);

/* CPU mirror of the multi-network simulator (A4.5). Same structs; pure C++ (no
 * CUDA) for the CRAN / macOS CPU-only build. Bit-identical trajectory + integer
 * outputs vs cpu.py::_simulate_one_multi; scores to float32 (ABI). */
CUSNA_API cusna_status cusna_saom_simulate_multi_cpu(
    const cusna_saom_multi_problem *problem, cusna_saom_multi_result *result);

/* --------------------------------------------------------------------------
 *  Utilities
 * -------------------------------------------------------------------------- */

/* Pack a dense row-major uint8 (n*n) adjacency into out/in bit rows (n*W each).
 * Caller allocates out_rows and in_rows with n*W uint32 words; *W_out = W. */
CUSNA_API cusna_status cusna_pack_bits(
    const uint8_t *dense, int32_t n,
    uint32_t *out_rows, uint32_t *in_rows, int32_t *W_out);

/* Philox4x32-10 one-block known-answer test (RNG parity check for A1). */
CUSNA_API cusna_status cusna_philox4x32(
    const uint32_t counter[4], const uint32_t key[2], uint32_t out[4]);

/* Host evaluation statistics of ONE network (float64 port of stats.py
 * network_stats_ex; milestone A4.1). `a` is the pre-transformed statistic
 * network (dense n*n uint8, row-major); descriptors as in cusna_saom_problem
 * (v = vb where dyn_flags[e], else vmat[e*n+i]); od0/id0 are start degrees for
 * parm<=0 degree effects. Fills out[k]. This is the shared host code both
 * backends use for eval statistics and targets. */
CUSNA_API cusna_status cusna_network_stats(
    const uint8_t *a, int32_t n, int32_t k,
    const int32_t *type_ids, const double *p1, const double *p2,
    const double *vmat, const uint8_t *vmiss, const int32_t *dyn_flags,
    const double *vb, const int32_t *od0, const int32_t *id0,
    const double *wtab, int32_t symflag, double *out);

/* Creation/endowment statistic of one typed effect over a batch of end networks
 * (float64 port of stats.py typed_stats_batch; convention 16, egoPlusAltX).
 * etype: 1 creation (gained ties), 2 endowment (lost ties, reported NEGATED).
 * Dense n*n uint8 masks (0/1): keep/force (statistic network), keepS/forceS
 * (structural-only replacement), miss_m/miss_m1 (missing at period start/end);
 * a_batch is M*n*n, F the (n*n) per-pair effect matrix, x0 the start network.
 * Fills out[M]. */
CUSNA_API cusna_status cusna_typed_stats(
    const uint8_t *a_batch, int32_t n, int32_t M, int32_t etype,
    const double *F, const uint8_t *x0,
    const uint8_t *keep, const uint8_t *force,
    const uint8_t *keepS, const uint8_t *forceS,
    const uint8_t *miss_m, const uint8_t *miss_m1, double *out);

/* sienaGOF auxiliary distributions over a batch of end networks (float64 port
 * of gof.py; milestone A4.5). kind: 0 indegree, 1 outdegree -- cumulative counts
 * of actors with degree <= levls[.], out is M*nlev. kind 2 geodesic -- counts of
 * ordered pairs (i!=j) at geodesic distance <= levls[.] plus the unreachable
 * count, out is M*(nlev+1). a_batch is M*n*n (the masked statistic network);
 * levls holds nlev cut levels. The Mahalanobis test is left to the caller (a
 * small pinv/quadratic form, natural in R/numpy). */
CUSNA_API cusna_status cusna_gof_distribution(
    const uint8_t *a_batch, int32_t n, int32_t M, int32_t kind,
    const int32_t *levls, int32_t nlev, double *out);

/* Behavior evaluation statistics over a batch of end behavior vectors (float64
 * port of stats.py behavior_stats; milestone A4.5, convention 9). Carries the
 * RSiena StatisticCalculator::calculateBehaviorStatistics conventions:
 *   - values centered by the overall mean (zbar); actors missing at either
 *     period end (bkdist[i]==0) contribute nothing (their centered value is
 *     forced to 0, and they are dropped as ego and as alter);
 *   - the predictor network is the START-of-period observed network less ties
 *     missing at either end (xnet = network_less_missing_start, dense n*n uint8,
 *     a cross-lagged moment);
 *   - avSim/totSim exclude alters missing at either end from BOTH the similarity
 *     sum and the neighbour count nb; totSim subtracts nb*simMean (the
 *     -outdeg*simMean term); avAlt uses the full out-degree and does not filter
 *     alters. bty holds kb effect type ids (0 linear, 1 quad, 2 avAlt, 3 avSim,
 *     4 totSim). z_batch is M*n int32; fills out[M*kb]. */
CUSNA_API cusna_status cusna_behavior_stats(
    const int32_t *z_batch, int32_t n, int32_t M, int32_t kb,
    const int32_t *bty, const uint8_t *xnet, const int32_t *bkdist,
    double zbar, double brange, double bsimmean, double *out);

/* Cross-network evaluation statistic of ONE cross effect over a batch of end
 * networks (float64 port of stats.py cross_stat; milestone A4.5, multi-net).
 * `a_batch` is M*n*n dense uint8 -- the OWN masked end network per chain; `wp`
 * is the n*n dense uint8 predictor = the SOURCE network at its period START
 * (network_less_missing_start). t: 33 crprod (sum a_ij*wp_ij), 34 crprodRecip
 * (a_ij*wp_ji), 35 crprodMutual (a_ij*wp_ij*wp_ji), 36 from (a_ij*(wp@wp^T)_ij,
 * sqrt when parm>=2), 37 to (a_ij*(wp@a)_ij, second step over the OWN net).
 * Fills out[M]. */
CUSNA_API cusna_status cusna_cross_stats(
    const uint8_t *a_batch, const uint8_t *wp, int32_t n, int32_t M,
    int32_t t, double parm, double *out);

/* --------------------------------------------------------------------------
 *  Layer 1c: ERGM TNT sampler ABI (track A-ERGM, milestones E2-E5).
 *
 *  cusna_ergm_tnt_sample runs C independent TNT Metropolis-Hastings chains of an
 *  exponential-family random graph model, ONE CUDA thread per chain (chains are
 *  lighter than SAOM periods) -- a faithful transcription of prototype/cusna/ergm/
 *  tnt_gpu.py::ErgmTntSampler.sample over the ergm_tnt kernel (native/
 *  ergm_kernels.cu, extracted verbatim). The kernel reproduces the ergm package's
 *  TNT proposal exactly (Mp_TNT, P=0.5); the Python sampler already matches
 *  ergm::simulate distributionally (comparison_ergm.md, max|z|=2.26 over 4 models).
 *
 *  Terms (type id): 0 edges | 1 mutual (directed) | 2 nodecov | 3 nodematch |
 *  4 absdiff | 5 gwesp | 6 gwdsp (undirected). At most CUSNA_ERGM_MAXTERMS.
 *  Edge packing uses 16-bit node ids, so n < 65536. This is a SAMPLER + MPLE, not
 *  a full MCMC-MLE estimator (framing M-3, kept uniform in the docs/paper).
 * -------------------------------------------------------------------------- */
#define CUSNA_ERGM_MAXTERMS 8   /* must match MAXTERMS in ergm_kernels.cu       */

typedef struct cusna_ergm_problem {
    int32_t struct_version;   /* set to CUSNA_ABI_VERSION                       */
    int32_t n;                /* actors (< 65536; edge list uses 16-bit ids)   */
    int32_t W;                /* uint32 words per actor = (n + 31) / 32        */
    int32_t k;                /* terms (<= CUSNA_ERGM_MAXTERMS)                */
    int32_t C;                /* independent chains (one CUDA thread each)     */
    int32_t directed;         /* 1 = directed (enables the mutual term)        */
    int32_t threads;          /* block size; 0 => 64 (TPB in sample())         */

    /* start network, dense n*n uint8 row-major (a_ij = tie i->j). The launcher
     * packs the out-adjacency bit rows and builds the per-chain edge list exactly
     * as sample() does (row-major nonzero; strictly-lower triangle if undirected,
     * so the RNG's edge-index draws select the same edges as the Python path). */
    const uint8_t *x0;

    /* term descriptors, length k */
    const float   *theta;     /* natural parameters (float32, as the kernel)   */
    const int32_t *type_ids;  /* term type ids (see the list above)            */
    const float   *parm;      /* per-term parameter (gwesp/gwdsp decay); else 0*/
    const float   *attr;      /* [opt] length n node attribute; NULL => zeros  */

    /* start-network statistics, length k: the kernel initialises each chain's
     * running statistics to these, then accumulates change statistics. [opt]:
     * NULL asks the launcher to compute them via cusna_ergm_stats (E3); until E3
     * s0 is REQUIRED (else CUSNA_ERR_UNSUPPORTED). Passing the reference-computed
     * s0 makes the GPU path bit-identical to the Python-GPU golden in isolation
     * from the host stats port. */
    const double  *s0;

    /* sampling schedule: samples are recorded at burnin + s*interval MH steps
     * (s = 0..samplesize-1); the final chain state equals the last sample. */
    uint64_t seed;
    int32_t burnin;
    int32_t interval;
    int32_t samplesize;

    /* [opt] temporal ERGM (TERGM): a fixed lagged network Y^{t-1}, dense n*n
     * uint8 row-major, enabling the temporal terms memory (id 7, change stat =
     * lag_ij) and delrecip (id 8, change stat = lag_ji). NULL when no temporal
     * terms are used (base behaviour unchanged). Currently the CPU sampler
     * honours these; the GPU sampler rejects temporal terms (X1b follow-on). */
    const uint8_t *lag;

    /* [opt] toggle-permitted (free) dyad mask, dense n*n uint8 row-major: 1 = the
     * dyad may toggle, 0 = it is fixed at its x0 value. NULL => all dyads free
     * (base behaviour unchanged). This drives the separable STERGM sub-samplers
     * (X2): the formation TNT toggles only the dyads absent at t-1 (x0 = union),
     * the dissolution TNT only those present at t-1 (x0 = intersection). The TNT
     * proposal is restricted to the free dyad set; change statistics are still
     * evaluated on the full network. Only the CPU sampler honours it. */
    const uint8_t *permitted;
} cusna_ergm_problem;

typedef struct cusna_ergm_result {
    int32_t struct_version;
    double  *stats;    /* [C*samplesize*k] absolute statistics at sample points */
    int32_t *nacc;     /* [opt] [C] accepted proposals per chain                */
    int32_t *nedges;   /* [opt] [C] final edge count per chain                  */
    uint8_t *a_ends;   /* [opt] [C*n*n] dense final networks (row-major)        */
} cusna_ergm_result;

/* Run C independent ERGM TNT chains on the GPU (ergm_tnt kernel). Returns
 * CUSNA_ERR_NO_CUDA if no device is present (use the _cpu mirror). Bit-identical
 * to the Python-GPU sampler on identical inputs: integer-valued terms exactly,
 * and gwesp/gwdsp too (the kernel computes their weights in float64 and there is
 * no cross-thread reduction, unlike SAOM's atomicAdd -- so GPU-vs-GPU is exact). */
CUSNA_API cusna_status cusna_ergm_tnt_sample(
    const cusna_ergm_problem *problem, cusna_ergm_result *result);

/* CPU mirror (milestone E4): a pure C++/OpenMP TNT-MH loop, no CUDA, for the CRAN
 * CPU-only / macOS build. There is NO bit-exact reference (the Python sampler is
 * GPU-only), so it is validated DISTRIBUTIONALLY against the GPU path and
 * ergm::simulate. Until E4: CUSNA_ERR_UNSUPPORTED. */
CUSNA_API cusna_status cusna_ergm_tnt_sample_cpu(
    const cusna_ergm_problem *problem, cusna_ergm_result *result);

/* Host statistics of ONE network for the ERGM term list (float64 port of
 * tnt_gpu.py network_stats_ergm; milestone E3). Dense n*n uint8 `a`; term
 * descriptors as in cusna_ergm_problem (attr widened to double, NULL => zeros).
 * Fills out[k]. Machine-zero match to the Python reference; also the s0 provider
 * for the samplers when problem->s0 is NULL. Until E3: CUSNA_ERR_UNSUPPORTED. */
CUSNA_API cusna_status cusna_ergm_stats(
    const uint8_t *a, int32_t n, int32_t k, int32_t directed,
    const int32_t *type_ids, const float *parm, const float *attr, double *out);

/* ERGM MPLE change-statistic design matrix (milestone E5; framing M-3: cusna's
 * ERGM is a sampler + this MPLE DEMO, not a general MCMC-MLE estimator). Ports the
 * mple_gpu.py demo -- three DIRECTED terms (edges, mutual, ttriple). For every
 * off-diagonal dyad (i,j) in row-major order it writes one design row
 *   [ 1, x_ji, |out(i)&in(j)|+|out(i)&out(j)|+|in(i)&in(j)| ]
 * into X (row-major) and the response x_ij into y. Integer change statistics ->
 * exact vs mple_change_stats_numpy. Caller allocates X[n*n*3] and y[n*n]; the
 * diagonal is skipped, so *nrow = n*(n-1), *ncol = 3. Host-only (works CPU-only
 * for CRAN); the pseudo-likelihood GLM fit is standard and left to the caller. */
CUSNA_API cusna_status cusna_ergm_mple_stats(
    const uint8_t *a, int32_t n, float *X, uint8_t *y,
    int32_t *nrow, int32_t *ncol);

/* --------------------------------------------------------------------------
 *  ALAAM (autologistic actor attribute model; track A-ERGM extension X1.5).
 *
 *  ALAAM models a BINARY actor attribute y at a FIXED (exogenous) network A -- the
 *  dependent variable is the attribute, not a tie. It is the shared bit-kernel
 *  primitive on (network-row AND attribute-vector): the change statistic of the
 *  contagion configuration when actor i flips its attribute is the masked popcount
 *  of A's row i against y. No CRAN reference (ALAAMEE is GitHub/Python); validated
 *  by the exact equivalence of the MPLE to a logistic regression on these change
 *  statistics (native/test/compare_alaam.R). Host-only (works CPU-only for CRAN).
 * -------------------------------------------------------------------------- */

/* Per-actor contagion change statistics for the ALAAM MPLE:
 *   con_out[i] = sum_j A[i,j]*y[j]   (out-neighbours of i with y=1)
 *   con_in[i]  = sum_j A[j,i]*y[j]   (in-neighbours of i with y=1)
 * Undirected (symmetric) A: the two coincide = neighbour count with y=1. `a` is
 * dense n*n uint8 row-major, `y` a length-n 0/1 attribute. Integer-exact. */
CUSNA_API cusna_status cusna_alaam_stats(
    const uint8_t *a, const uint8_t *y, int32_t n,
    double *con_out, double *con_in);

/* ALAAM Gibbs sampler over the binary attribute at a fixed network. The full
 * conditional is logit P(y_i=1|y_{-i}) = beta0 + beta_c*c_i + sum_p beta_x[p]*X[i,p],
 * c_i = sum_j A[i,j]*y_j (+ sum_j A[j,i]*y_j when directed) -- the same contagion
 * the MPLE regresses on. One chain per OpenMP iteration, `sweeps` Gibbs sweeps
 * from y0, host Philox RNG. Writes end attributes into y_ends[C*n]; if stats!=NULL
 * writes [density, contagion] per chain into stats[C*2]. X is [n*p] row-major
 * (NULL/p=0 for none), beta_x length p. */
CUSNA_API cusna_status cusna_alaam_gibbs_cpu(
    const uint8_t *a, int32_t n, int32_t directed,
    const uint8_t *y0, const float *X, int32_t p,
    double beta0, double beta_c, const float *beta_x,
    int32_t C, uint64_t seed, int32_t sweeps,
    uint8_t *y_ends, double *stats);

CUSNA_API int32_t     cusna_abi_version(void);      /* CUSNA_ABI_VERSION      */
CUSNA_API const char *cusna_version_string(void);   /* e.g. "cusna-native 0.1"*/
CUSNA_API int32_t     cusna_has_cuda(void);         /* 1 if a device is usable*/
CUSNA_API const char *cusna_last_error(void);       /* human-readable detail  */

/* --------------------------------------------------------------------------
 *  Layer 2: estimator ABI (A4/A5) -- DECLARED FOR DIRECTION, DEFINED IN A4.
 *
 *  The high-level entry that R will actually call. It will own the full RSiena
 *  marshalling (Appendix B conventions), target computation, the Robbins-Monro
 *  MoM loop (mom.py) and sienaGOF (gof.py). Left opaque here so callers and the
 *  R glue can be sketched; the concrete struct is frozen once A4 re-validates
 *  the 23-check harness on the native stack.
 * -------------------------------------------------------------------------- */
typedef struct cusna_mom_options cusna_mom_options;   /* opaque, defined in A4 */
typedef struct cusna_mom_result  cusna_mom_result;    /* opaque, defined in A4 */

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif  /* LIBCUSNA_H */
