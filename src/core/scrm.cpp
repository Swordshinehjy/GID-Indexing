// ============================================================================
//  scrm.cpp — SCRM w-free two-step method (Specular-Constrained Reciprocal Metric)
//
//  Step one (enumerate only (u,v), without w):
//      t_i = q_z,i / q_s                       (pure data)
//      H   = h − u·t_i ,  K = k − v·t_i
//      q_xy,i² = (X·H² + Y·K² − 2Zc·H·K) / (1 − g_uv/q_s²)      (3 seed peaks → 3×3 solve for X,Y,Zc)
//      Λ   = (XY − Zc²) / (q_s² − g_uv)
//      z_a² = X + v²Λ ,  z_b² = Y + u²Λ ,  C = Zc + u·v·Λ       (Λ transform → in-plane metric)
//  Criteria: z_a²>0, z_b²>0, z_a²z_b²−C²>0, q_s² > g_uv (⇔ ζ²>0), unit-cell range, 2D Niggli.
//  Note: q_xy² cannot be written directly as A·H²+B·K²−2C·H·K — that holds only for the (001) contact plane; on a tilted contact plane
//      the omitted normal component (κH+λK+ρL)² can reach 95% of q_xy² (measured 0.94 for DIP uvw=(−1,2,1)).
//      The expression above is an identity for the solved (X,Y,Zc) (including the (AB−C²)(hv−ku)²/q_s² term), see wiki §5.3.
//
//  Step two (fully analytic w, l):
//      ζ   = sqrt(q_s² − g_uv)
//      s3  = A h u + B k v − C(h v + k u)
//      g_p = (q_z·q_s − s3)/ζ                  (synonymous with gz in improved)
//      solve [h k l]·(κ,λ,ρ) = g_p             (3 seed peaks)
//      w   = round((ζ − uκ − vλ)/ρ)            ← w derived from the specular constraint nᵀt = q_spec
//      l_j = round((g_p,j − h_jκ − k_jλ)/ρ)
//      G*  = {A+κ², B+λ², ρ², −C+κλ, κρ, λρ}   → cellFromFrame → refinement/reduction
//
//  Derivation in docs/SCRM_specular_two_step_design_wiki.md §5 (step one) and §6 (step two).
// ============================================================================
#include "scrm.h"

#include "engine.h" // alignSolution (completes specular rod points as m·uvw)
#include "linalg.h"
#include "parallel.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gid {
namespace scrm {
namespace {

constexpr double PI2 = 6.28318530717958647692;
constexpr double DEG_PER_RAD = 180.0 / M_PI;

inline double clamp1(double x) { return x > 1.0 ? 1.0 : (x < -1.0 ? -1.0 : x); }

// ---------------------------------------------------------------------------
//  Funnel diagnostics (enabled by environment variable SCRM_TRACE=1)
//  gFunnel[i] records the number of candidates that "pass stage i"; the difference between adjacent stages is the number eliminated there.
//  Parallel runNormal calls accumulate atomically; runSpecIndex prints everything at the end.
// ---------------------------------------------------------------------------
enum FunnelStage {
    kStageCall = 0, // runNormal invoked
    kStageSeed,     // seed peaks complete and argument valid
    kStageRho,      // trySol entry: ρ>0 and (κ,λ) finite
    kStageW,        // w integer consistency (including fixedUvw)
    kStageCell,     // cellFromFrame + restrictCell(c/α/β)
    kStageSpecPre,  // specular peak self-consistency |qSpec−qs| ≤ dqspec (pre-check, coarse filter)
    kStageRms3,     // per-peak l recovered and rmsQxyz ≤ dqxyzCutoff
    kStageRmsZ,     // rmsQz ≤ 3·dqzCutoff
    kStageRaw,      // entering the refinement queue (raw-key dedup and refineTop truncation already done)
    kStageReduce,   // refinement + reduceCell succeeded
    kStageRmsQ,     // rmsQxyz ≤ dqxyzCutoff after evaluate
    kStageLattice,  // latticeIndex ≤ 1
    kStageRmsXY,    // rmsQxy/rmsQz ≤ 3×threshold
    kStageSpecFin,  // |qSpec−qs| ≤ dqspec after reduction
    kStageRange,    // restrictCell(a/b/c/vol)
    kStageUnique,   // solution key first occurrence (dedup passed)
    kStageSolution, // finally stored
    kStageCount
};

inline const char* funnelName(int stage) {
    static const char* const names[kStageCount] = {
        "runNormal call",  "seed peaks complete",       "ρ>0 and (κ,λ) finite", "w integer consistent",
        "cell construction+out-of-plane constraint", "specular peak self-consistency (pre)", "rmsQxyz (pre)",     "rmsQz (pre)",
        "entering refinement queue",     "reduceCell succeeded",   "rmsQxyz (post-reduction)",   "latticeIndex≤1",
        "rmsQxy/rmsQz (post-reduction)", "specular peak self-consistency (post-reduction)", "unit-cell range",     "dedup passed",
        "finally stored"};
    return (stage >= 0 && stage < kStageCount) ? names[stage] : "?";
}

std::atomic<long long> gFunnel[kStageCount];
std::atomic<int> gTraceFlag{0}; // 0 not probed / 1 on / -1 off

inline bool envFlag(const char* name) {
#if defined(_MSC_VER)
    char* buf = nullptr;
    size_t sz = 0;
    const bool has = (_dupenv_s(&buf, &sz, name) == 0 && buf != nullptr);
    std::free(buf);
    return has;
#else
    return std::getenv(name) != nullptr;
#endif
}

inline bool traceOn() {
    int t = gTraceFlag.load(std::memory_order_relaxed);
    if (t == 0) {
        t = envFlag("SCRM_TRACE") ? 1 : -1;
        gTraceFlag.store(t, std::memory_order_relaxed);
    }
    return t > 0;
}

inline void bump(int stage) { gFunnel[stage].fetch_add(1, std::memory_order_relaxed); }







inline void funnelDump() {
    if (!traceOn()) return;
    std::fprintf(stderr, "---- SCRM funnel statistics (in-process cumulative; difference between adjacent stages = eliminated at that stage) ----\n");
    for (int i = 0; i < kStageCount; ++i)
        std::fprintf(stderr, "  [%2d] %-22s %12lld\n", i, funnelName(i),
                     gFunnel[i].load(std::memory_order_relaxed));
    std::fprintf(stderr, "----------------------------------------------\n");
}

inline uint64_t mixHash(uint64_t h, long long x) {
    h ^= (uint64_t)(x + 0x9e3779b97f4a7c15LL);
    h *= 1099511628211ull;
    return h;
}

// Dedup key for in-plane solutions (0.15 Å / 0.3°, same tolerance as improved's abgKey)
inline uint64_t abgKey(int u, int v, double a, double b, double ga) {
    auto q = [](double x, double t) -> long long { return (long long)std::llround(x / t); };
    uint64_t h = 1469598103934665603ull;
    h = mixHash(h, u);
    h = mixHash(h, v);
    h = mixHash(h, q(a, 0.15));
    h = mixHash(h, q(b, 0.15));
    h = mixHash(h, q(ga, 0.3));
    return h;
}

inline uint64_t cellKeyOf(const Cell& c) {
    auto q = [](double x, double t) -> long long { return (long long)std::llround(x / t); };
    uint64_t h = 1469598103934665603ull;
    h = mixHash(h, q(c.a, 0.02));
    h = mixHash(h, q(c.b, 0.02));
    h = mixHash(h, q(c.c, 0.02));
    h = mixHash(h, q(c.al, 0.1));
    h = mixHash(h, q(c.be, 0.1));
    h = mixHash(h, q(c.ga, 0.1));
    return h;
}

// Coarse dedup key before refinement (0.25 Å / 2°): raws in the same basin of attraction are refined only once
inline uint64_t rawKeyOf(const Cell& c, int w) {
    auto q = [](double x, double t) -> long long { return (long long)std::llround(x / t); };
    uint64_t h = 1469598103934665603ull;
    h = mixHash(h, q(c.a, 0.25));
    h = mixHash(h, q(c.b, 0.25));
    h = mixHash(h, q(c.c, 0.25));
    h = mixHash(h, q(c.al, 2.0));
    h = mixHash(h, q(c.be, 2.0));
    h = mixHash(h, q(c.ga, 2.0));
    h = mixHash(h, w);
    return h;
}

struct RawIP {
    int u = 0, v = 0;
    double A = 0.0, B = 0.0, C = 0.0;
    double a = 0.0, b = 0.0, ga = 90.0;
    double guv = 0.0, zeta = 0.0;
};

struct RawSol {
    Cell cell;
    int w = 0;
    double rms = 0.0;
    std::vector<Hkl> hkl;
};

// ---------------------------------------------------------------------------
//  Per-peak assignment / refinement of an in-plane candidate (tail of step one; a separate function avoids overly deep lambda nesting)
// ---------------------------------------------------------------------------
struct IpCtx {
    const std::vector<Peak>& peaks;
    const std::vector<double>& tOf;
    int N = 0;
    double qs2 = 0.0;
    int R = 4;
    int L2 = 6;
};

// With the current (A,B,C), pick the R nearest (h,k) groups for each peak and return RMSD(q_xy).
//   In-plane prediction: q_xy² = [A H² + B K² − 2C·HK − (AB−C²)(h·v−k·u)²/q_s²] / (1 − g_uv/q_s²)
static double assignHK(const IpCtx& ctx, InPlaneCand& c) {
    const double du = (double)c.u, dv = (double)c.v;
    const double det2D = c.A * c.B - c.C * c.C;
    const double denIp = 1.0 - c.guv / ctx.qs2;
    c.hk.assign((size_t)ctx.N, {});
    if (!(denIp > 1e-9) || !(det2D > 0.0)) return 1e30;
    double sq = 0.0;
    int used = 0; // number of peaks that actually have (h,k) candidates: skipped peaks are not counted in the denominator, otherwise RMSD is systematically underestimated
    const int R = std::min(ctx.R, 8);
    for (int i = 0; i < ctx.N; ++i) {
        const double Q = ctx.peaks[i].qxy;
        const double uz = du * ctx.tOf[i], vz = dv * ctx.tOf[i];
        // Keep only the R nearest (h,k) groups: use "ordered min-heap-style insertion" instead of fully sorting 169 elements;
        // the resulting set and order (ascending distance) match a full sort, but the per-peak per-round vector allocation is removed.
        double bd[8];
        Hkl bh[8];
        int nb = 0;
        for (int h = -ctx.L2; h <= ctx.L2; ++h)
            for (int k = -ctx.L2; k <= ctx.L2; ++k) {
                const double H = h - uz, K = k - vz;
                const double rod = (double)h * dv - (double)k * du;
                const double num =
                    c.A * H * H + c.B * K * K - 2.0 * c.C * H * K - det2D * rod * rod / ctx.qs2;
                if (!(num >= 0.0) || !std::isfinite(num)) continue;
                const double qp = std::sqrt(num / denIp);
                const double d = std::fabs(qp - Q);
                if (nb == R && d >= bd[R - 1]) continue;
                int pos = (nb < R) ? nb : R - 1;
                while (pos > 0 && d < bd[pos - 1]) {
                    bd[pos] = bd[pos - 1];
                    bh[pos] = bh[pos - 1];
                    --pos;
                }
                bd[pos] = d;
                bh[pos] = Hkl{h, k, 0};
                if (nb < R) ++nb;
            }
        for (int j = 0; j < nb; ++j) c.hk[(size_t)i].push_back(bh[j]);
        if (nb > 0) {
            sq += bd[0] * bd[0];
            ++used;
        }
    }
    // If no peak has a feasible (h,k), RMSD is undefined; return an unacceptable value (the caller rejects by threshold)
    if (used == 0) return 1e30;
    return std::sqrt(sq / (double)used);
}

// In-plane metric refinement: weighted LSQ over all peaks, two rounds of "solve → reassign (h,k)".
//   The three-step seed solve has only 3 equations, so measurement error biases (A,B,C) and RMSD(q_xy) often sits near the threshold,
//   while step two's analytic inversion (s3, g_p, w) is extremely sensitive to (A,B,C) — one refinement round usually reduces rms by an order of magnitude.
//   The model is the same as step one: Q² = X·(H²+Q²u²/q_s²) + Y·(K²+Q²v²/q_s²) + Zc·(−2HK−2Q²uv/q_s²),
//   weight 1/Q² (error propagation of Q²).
static void refineIp(const IpCtx& ctx, const Params& prm, InPlaneCand& c) {
    const double du = (double)c.u, dv = (double)c.v;
    const double du2 = du * du, dv2 = dv * dv, duv = du * dv;
    const double qs2 = ctx.qs2;
    for (int it = 0; it < 2; ++it) {
        Mat Ap{3, 3, 0.0}, bp{3, 1, 0.0};
        int used = 0;
        for (int i = 0; i < ctx.N; ++i) {
            if (c.hk[(size_t)i].empty()) continue;
            const double q2 = ctx.peaks[i].qxy * ctx.peaks[i].qxy;
            const double H = (double)c.hk[(size_t)i][0].h - du * ctx.tOf[i];
            const double K = (double)c.hk[(size_t)i][0].k - dv * ctx.tOf[i];
            const double row[3] = {H * H + du2 * q2 / qs2, K * K + dv2 * q2 / qs2,
                                   -2.0 * H * K - 2.0 * duv * q2 / qs2};
            const double wgt = 1.0 / q2;
            for (int rr = 0; rr < 3; ++rr) {
                for (int ss = 0; ss < 3; ++ss) Ap(rr, ss) += wgt * row[rr] * row[ss];
                bp(rr, 0) += wgt * row[rr] * q2;
            }
            ++used;
        }
        if (used < 3) return;
        Mat x;
        if (!matSolve(Ap, bp, x)) return;
        const double X = x(0, 0), Y = x(1, 0), Zc = x(2, 0);
        const double den = qs2 - du2 * X - dv2 * Y + 2.0 * duv * Zc;
        if (!(den > 1e-12)) return;
        const double Lam = (X * Y - Zc * Zc) / den;
        const double za2 = X + dv2 * Lam, zb2 = Y + du2 * Lam;
        if (!(za2 > 1e-12) || !(zb2 > 1e-12)) return;
        const double A2 = za2, B2 = zb2, C2 = Zc + duv * Lam;
        if (!(A2 * B2 - C2 * C2 > 0.0)) return;
        const double za = std::sqrt(A2), zb = std::sqrt(B2);
        const double cg = clamp1(C2 / (za * zb));
        const double sg = std::sqrt(std::max(1e-12, 1.0 - cg * cg));
        const double a2 = PI2 / (za * sg), b2 = PI2 / (zb * sg);
        const double ga2 = std::acos(cg) * DEG_PER_RAD;
        if (!(a2 > 0.0) || !(b2 > 0.0) || !std::isfinite(ga2)) return;
        if (a2 < prm.aLo || b2 < prm.bLo) return;
        if (prm.restrictCell && (a2 > prm.aHi || b2 > prm.bHi || ga2 < prm.gaLo || ga2 > prm.gaHi))
            return;
        if (prm.niggli && !niggli2D(a2, b2, ga2)) return;
        const double guv = du2 * X + dv2 * Y - 2.0 * duv * Zc;
        if (!(qs2 > guv)) return;

        InPlaneCand nc = c;
        nc.A = A2;
        nc.B = B2;
        nc.C = C2;
        nc.a = a2;
        nc.b = b2;
        nc.ga = ga2;
        nc.guv = guv;
        nc.zeta = std::sqrt(qs2 - guv);
        nc.rmsQxy = assignHK(ctx, nc);
        if (!(nc.rmsQxy <= c.rmsQxy)) return; // stop if it did not improve
        c = std::move(nc);
    }
}

} // namespace

// ---------------------------------------------------------------------------
//  Step one: enumerate only (u,v)
//    gateMul  : candidate RMSD(q_xy) threshold (× dqxyCutoff)
//    scanMul  : scan budget (× maxPart1) — walk in ascending area order and stop when the budget is exhausted
// ---------------------------------------------------------------------------
static std::vector<InPlaneCand> runInPlaneImpl(const std::vector<Peak>& peaksIn,
                                               const Params& prm, const Progress& prog,
                                               double gateMul, int scanMul,
                                               Stats* st = nullptr) {
    // Points on the specular rod (q_xy≈0) do not participate in the in-plane fit (their rows are identically zero)
    std::vector<Peak> peaks;
    for (const Peak& p : peaksIn)
        if (p.qxy > 1e-6 && std::isfinite(p.qxy) && std::isfinite(p.qz)) peaks.push_back(p);
    const int N = (int)peaks.size();
    if (st) st->peaksUsed = std::max(st->peaksUsed, N);
    if (N < 3 || !(prm.qspec > 0.0)) return {};

    const double qs = prm.qspec, qs2 = qs * qs;
    const int L1 = std::max(1, prm.hkLimit1);
    const int W1 = 2 * L1 + 1;
    const int W2 = W1 * W1;
    const int nTri = std::min(N, std::max(3, prm.linesToPermute));

    // t = q_z/q_s: the only entry point for out-of-plane information (pure data, w does not appear)
    std::vector<double> tOf(N), qlen(N);
    for (int i = 0; i < N; ++i) {
        tOf[i] = peaks[i].qz / qs;
        qlen[i] = std::sqrt(peaks[i].qxy * peaks[i].qxy + peaks[i].qz * peaks[i].qz);
    }
    const double qTol = 5e-4;

    // Initial three reflections: two seed sets (|q| ascending / q_xy ascending), same as the original and improved versions
    auto pickSeeds = [&](bool byQlen) {
        std::vector<int> ord(N);
        for (int i = 0; i < N; ++i) ord[i] = i;
        if (byQlen)
            std::sort(ord.begin(), ord.end(), [&](int i, int j) {
                if (qlen[i] != qlen[j]) return qlen[i] < qlen[j];
                return i < j;
            });
        else
            std::sort(ord.begin(), ord.end(), [&](int i, int j) {
                if (peaks[i].qxy != peaks[j].qxy) return peaks[i].qxy < peaks[j].qxy;
                return i < j;
            });
        std::vector<int> out;
        double rep = 0.0;
        bool first = true;
        for (int i : ord) {
            const double key = byQlen ? qlen[i] : peaks[i].qxy;
            const double tol = byQlen ? qTol : 1e-3;
            if (first || std::fabs(key - rep) > tol) {
                rep = key;
                first = false;
                out.push_back(i);
                if ((int)out.size() >= nTri) break;
            }
        }
        if ((int)out.size() < 3)
            for (int i : ord) {
                if (std::find(out.begin(), out.end(), i) == out.end()) out.push_back(i);
                if ((int)out.size() >= 3) break;
            }
        return out;
    };

    std::vector<std::array<int, 3>> triples;
    std::vector<std::array<int, 3>> triplesSorted;
    for (bool byQlen : {true, false}) {
        const std::vector<int> sel = pickSeeds(byQlen);
        const int nT = std::min((int)sel.size(), N);
        for (int i = 0; i < nT; ++i)
            for (int j = i + 1; j < nT; ++j)
                for (int k = j + 1; k < nT; ++k) {
                    const std::array<int, 3> t = {sel[i], sel[j], sel[k]};
                    std::array<int, 3> key = t;
                    std::sort(key.begin(), key.end());
                    if (std::find(triplesSorted.begin(), triplesSorted.end(), key) ==
                        triplesSorted.end()) {
                        triplesSorted.push_back(key);
                        triples.push_back(t);
                    }
                }
    }
    if (triples.empty() && N >= 3) triples.push_back({0, 1, 2});

    // (u,v) half-plane: (u,v) and (−u,−v) give the same (A,B,C) (with h,k negated simultaneously)
    std::vector<std::pair<int, int>> uvs;
    const int UL = std::max(0, prm.uvLimit);
    if (prm.fixedUvw) {
        uvs.push_back({prm.fu, prm.fv});
    } else {
        for (int u = -UL; u <= UL; ++u)
            for (int v = -UL; v <= UL; ++v) {
                if (u < 0) continue;
                if (u == 0 && v < 0) continue;
                uvs.push_back({u, v});
            }
    }
    if (uvs.empty()) uvs.push_back({0, 0});
    if (st) st->uvDirs += (int)uvs.size();

    // Row table: one per (u,v). Row = (H² + Q²u²/q_s², K² + Q²v²/q_s², −2HK − 2Q²uv/q_s²),
    // right-hand side = Q² (Q = peak q_xy). The added constant terms come from the Λ reparameterization after
    // "folding the normal component (κH+λK+ρL)² into the in-plane model": q_xy² = (X H² + Y K² − 2Zc·HK)/(1 − g_uv/q_s²), where H = h − u·t,
    // K = k − v·t, t = q_z/q_s. Without these three terms (old implementation), the error on a tilted contact plane can reach 95% of q_xy².
    std::vector<std::vector<double>> tabs(uvs.size());
    for (size_t ui = 0; ui < uvs.size(); ++ui) {
        const double du = (double)uvs[ui].first, dv = (double)uvs[ui].second;
        const double du2 = du * du, dv2 = dv * dv, duv = du * dv;
        std::vector<double>& tab = tabs[ui];
        tab.assign((size_t)N * W2 * 3, 0.0);
        for (int i = 0; i < N; ++i) {
            const double uz = du * tOf[i], vz = dv * tOf[i];
            const double q2 = peaks[i].qxy * peaks[i].qxy;
            const double add0 = du2 * q2 / qs2;
            const double add1 = dv2 * q2 / qs2;
            const double add2 = -2.0 * duv * q2 / qs2;
            double* row = tab.data() + (size_t)i * W2 * 3;
            for (int h = -L1; h <= L1; ++h)
                for (int k = -L1; k <= L1; ++k) {
                    const double H = h - uz, K = k - vz;
                    double* o = row + ((size_t)((h + L1) * W1 + (k + L1))) * 3;
                    o[0] = H * H + add0;
                    o[1] = K * K + add1;
                    o[2] = -2.0 * H * K + add2;
                }
        }
    }

    struct Task {
        int ui;
        int tri;
    };
    std::vector<Task> tasks;
    tasks.reserve(uvs.size() * triples.size());
    for (size_t ui = 0; ui < uvs.size(); ++ui)
        for (size_t t = 0; t < triples.size(); ++t) tasks.push_back(Task{(int)ui, (int)t});

    std::mutex mu;
    std::unordered_map<uint64_t, RawIP> best;
    std::atomic<int> done{0};

    parallelForEach((int)tasks.size(), [&](int ti) {
        if (prog.stop()) return;
        const Task tk = tasks[(size_t)ti];
        const int u = uvs[(size_t)tk.ui].first, v = uvs[(size_t)tk.ui].second;
        const double du = (double)u, dv = (double)v;
        const double du2 = du * du, dv2 = dv * dv, duv = du * dv;
        const std::array<int, 3>& tr = triples[(size_t)tk.tri];
        const double* const tabBase = tabs[(size_t)tk.ui].data();
        const double* const rowA = tabBase + (size_t)tr[0] * W2 * 3;
        const double* const rowB = tabBase + (size_t)tr[1] * W2 * 3;
        const double* const rowC = tabBase + (size_t)tr[2] * W2 * 3;
        const double q1 = peaks[tr[0]].qxy * peaks[tr[0]].qxy;
        const double q2 = peaks[tr[1]].qxy * peaks[tr[1]].qxy;
        const double q3 = peaks[tr[2]].qxy * peaks[tr[2]].qxy;

        std::vector<RawIP> local;
        std::unordered_set<uint64_t> localSeen;
        const int LA = (prm.seed3Half > 0) ? std::min(L1, prm.seed3Half) : L1;
        for (int h1 = -LA; h1 <= LA; ++h1)
            for (int k1 = -LA; k1 <= LA; ++k1) {
                const double* m0 = rowA + ((size_t)((h1 + L1) * W1 + (k1 + L1))) * 3;
                const double m0x = m0[0], m0y = m0[1], m0z = m0[2];
                const double s0 = std::fabs(m0x) + std::fabs(m0y) + std::fabs(m0z);
                for (int h2 = -L1; h2 <= L1; ++h2)
                    for (int k2 = -L1; k2 <= L1; ++k2) {
                        const double* m1 = rowB + ((size_t)((h2 + L1) * W1 + (k2 + L1))) * 3;
                        const double m1x = m1[0], m1y = m1[1], m1z = m1[2];
                        const double s01 = s0 * (std::fabs(m1x) + std::fabs(m1y) + std::fabs(m1z));
                        for (int h3 = -L1; h3 <= L1; ++h3)
                            for (int k3 = -L1; k3 <= L1; ++k3) {
                                const double* m2 = rowC + ((size_t)((h3 + L1) * W1 + (k3 + L1))) * 3;
                                const double m2x = m2[0], m2y = m2[1], m2z = m2[2];
                                const double c0x = m1y * m2z - m1z * m2y;
                                const double c0y = m1z * m2x - m1x * m2z;
                                const double c0z = m1x * m2y - m1y * m2x;
                                const double det = m0x * c0x + m0y * c0y + m0z * c0z;
                                const double s2 = std::fabs(m2x) + std::fabs(m2y) + std::fabs(m2z);
                                if (std::fabs(det) < 1e-9 * std::max(1e-12, s01 * s2)) continue;
                                const double inv = 1.0 / det;
                                // Cramer: replace the three columns with (q1,q2,q3) → (X, Y, Zc)
                                const double X =
                                    (q1 * c0x - m0y * (q2 * m2z - m1z * q3) +
                                     m0z * (q2 * m2y - m1y * q3)) * inv;
                                const double Y =
                                    (m0x * (q2 * m2z - m1z * q3) - q1 * (m1x * m2z - m1z * m2x) +
                                     m0z * (m1x * q3 - q2 * m2x)) * inv;
                                const double Zc =
                                    (m0x * (m1y * q3 - q2 * m2y) -
                                     m0y * (m1x * q3 - q2 * m2x) + q1 * c0z) * inv;
                                // Λ transform: (X,Y,Zc) → in-plane metric (z_a², z_b², C=z_a z_b cosγ)
                                const double den = qs2 - du2 * X - dv2 * Y + 2.0 * duv * Zc;
                                if (!(den > 1e-12)) continue; // ⇔ q_s² > g_uv ⇔ ζ² > 0
                                const double Lam = (X * Y - Zc * Zc) / den;
                                const double za2 = X + dv2 * Lam, zb2 = Y + du2 * Lam;
                                if (!(za2 > 1e-12) || !(zb2 > 1e-12)) continue;
                                const double A = za2, B = zb2, C = Zc + duv * Lam;
                                const double det2 = A * B - C * C; // = (z_a z_b sinγ)²
                                if (!(det2 > 0.0)) continue;       // in-plane metric is positive definite
                                // ---- squared-domain range filter (no more sqrt) ----
                                //   a² = PI2²·B/det2, b² = PI2²·A/det2 (derived from a=PI2/(z_a sinγ)),
                                //   so a>=aLo ⇔ PI2²·B >= aLo²·det2, no need to compute a/b.
                                const double a2n = PI2 * PI2 * B; // = a²·det2
                                const double b2n = PI2 * PI2 * A; // = b²·det2
                                if (prm.restrictCell) {
                                    if (a2n < prm.aLo * prm.aLo * det2 ||
                                        a2n > prm.aHi * prm.aHi * det2)
                                        continue;
                                    if (b2n < prm.bLo * prm.bLo * det2 ||
                                        b2n > prm.bHi * prm.bHi * det2)
                                        continue;
                                } else {
                                    // Hard lower bound (default 3 Å): degenerate "sliver cells" like 0.2 Å can fit
                                    // q_xy very well (small area, ranked first in the candidate list), but step two
                                    // can never succeed. The upper bound is applied only under restrictCell.
                                    if (a2n < prm.aLo * prm.aLo * det2) continue;
                                    if (b2n < prm.bLo * prm.bLo * det2) continue;
                                }
                                // ---- squared-domain 2D Niggli (equivalent to niggli2D, saving acos/cos) ----
                                //   a<=b ⇔ z_a>=z_b ⇔ A>=B; the 1e-9 Å tolerance becomes, in the squared domain,
                                //   PI2(z_b−z_a) <= 1e-9·sinγ·z_a z_b.
                                //   |b cosγ|<=a/2 ⇔ 2|C|<=B; the 1e-9 Å tolerance ⇔ + (1e-9/π)·sinγ·z_a·B.
                                const double za = std::sqrt(A), zb = std::sqrt(B);
                                const double cg = clamp1(C / (za * zb));
                                const double sg = std::sqrt(std::max(1e-12, 1.0 - cg * cg));
                                if (prm.niggli) {
                                    if (A < B && PI2 * (zb - za) > 1e-9 * sg * za * zb) continue;
                                    if (2.0 * std::fabs(C) > B + (1e-9 / M_PI) * sg * za * B)
                                        continue;
                                }
                                // Compute a,b,γ only for survivors (acos/division are deferred to here)
                                const double a = PI2 / (za * sg), b = PI2 / (zb * sg);
                                if (!std::isfinite(a) || !std::isfinite(b) || a <= 0.0 || b <= 0.0)
                                    continue;
                                const double ga = std::acos(cg) * DEG_PER_RAD;
                                if (!std::isfinite(ga)) continue;
                                if (prm.restrictCell) {
                                    if (ga < prm.gaLo || ga > prm.gaHi) continue;
                                }
                                const double guv = du2 * X + dv2 * Y - 2.0 * duv * Zc;
                                const uint64_t key = abgKey(u, v, a, b, ga);
                                if (localSeen.insert(key).second)
                                    local.push_back(RawIP{u,          v,    A, B, C, a,
                                                          b,          ga,   guv,
                                                          std::sqrt(qs2 - guv)});
                            }
                    }
            }
        if (prog.report)
            prog((double)(++done) / (double)tasks.size() * 0.9, "SCRM step one: enumerate (u,v)");
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(mu);
            // Deterministic merge: for the same abgKey keep only the representative with the lexicographically smallest (A,B,C).
            // "First come wins" is not allowed — thread merge order is nondeterministic, which would let the same data produce
            // different candidate sets (and hence different solutions) across runs, making benchmarks irreproducible.
            for (const RawIP& r : local) {
                const uint64_t key = abgKey(r.u, r.v, r.a, r.b, r.ga);
                auto it = best.find(key);
                if (it == best.end()) {
                    best.emplace(key, r);
                    continue;
                }
                const RawIP& o = it->second;
                const bool smaller = (r.A != o.A) ? (r.A < o.A) : ((r.B != o.B) ? (r.B < o.B)
                                                                              : (r.C < o.C));
                if (smaller) it->second = r;
            }
        }
    }, prm.threads);

    if (prog.stop()) return {};

    std::vector<RawIP> uniq;
    uniq.reserve(best.size());
    for (const auto& kv : best) uniq.push_back(kv.second);
    prog(0.92, "SCRM step one: in-plane indexing of all peaks");

    // Candidate order: in-plane parallelogram area ascending (same as the original/improved versions, preventing overfitted large cells from ranking first)
    // Area = 4π²/√det2D (det2D = A·B−C²), so area ascending ⇔ det2D descending, no sin/sqrt needed.
    auto areaDet2 = [](const RawIP& x) { return x.A * x.B - x.C * x.C; };
    std::sort(uniq.begin(), uniq.end(), [&](const RawIP& x, const RawIP& y) {
        const double dx = areaDet2(x), dy = areaDet2(y);
        if (dx != dy) return dx > dy;
        if (x.a != y.a) return x.a < y.a;
        if (x.b != y.b) return x.b < y.b;
        return x.ga < y.ga;
    });

    // Take the R nearest (h,k) per peak and compute RMSD(q_xy); in area order, stop when maxPart1 is filled or
    // the scan budget (scanMul × maxPart1) is exhausted.
    // Note: in area-ascending order the true candidate often ranks in the hundreds (large overfitted cells give smaller rms), so
    // filter with the RMSD(q_xy) threshold (not just maxPart1) — under a tight threshold the true value ranks only in the teens.
    const int R = std::max(1, prm.redundant);
    const int L2 = std::max(1, prm.hkLimit2);
    const double rmsGate = prm.dqxyCutoff * gateMul + 1e-9;
    const int scanBudget = std::max(1, prm.maxPart1) * std::max(1, scanMul);

    const IpCtx ctx{peaks, tOf, N, qs2, R, L2};

    std::vector<InPlaneCand> cands;
    cands.reserve(std::min<size_t>(uniq.size(), (size_t)std::max(1, prm.maxPart1) + 1));
    int scanned = 0;
    for (const RawIP& r : uniq) {
        if ((int)cands.size() >= prm.maxPart1) break;
        if (scanned >= scanBudget) break;
        ++scanned;
        if (!(1.0 - r.guv / qs2 > 1e-9)) continue;
        InPlaneCand c;
        c.u = r.u;
        c.v = r.v;
        c.A = r.A;
        c.B = r.B;
        c.C = r.C;
        c.a = r.a;
        c.b = r.b;
        c.ga = r.ga;
        c.guv = r.guv;
        c.zeta = r.zeta;
        c.rmsQxy = assignHK(ctx, c);
        if (!(c.rmsQxy <= rmsGate)) continue;
        refineIp(ctx, prm, c);
        cands.push_back(std::move(c));
    }
    prog(1.0, "SCRM step one complete");
    if (st) st->ipCands += (long long)cands.size();
    return cands;
}

// Public entry point: tight threshold (1.2 × dqxy) + 8× scan budget — on DIP data the true in-plane solution ranks only in the teens
// (in area order there are thousands of candidates ahead of it, but only a dozen or so pass the q_xy threshold).
// When nothing is found, runSpecIndex falls back with a wide threshold / more candidates.
std::vector<InPlaneCand> runInPlane(const std::vector<Peak>& peaks, const Params& prm,
                                    const Progress& prog, Stats* st) {
    return runInPlaneImpl(peaks, prm, prog, 1.2, 8, st);
}

// ---------------------------------------------------------------------------
//  Step two: analytic (κ,λ,ρ,w,l)
// ---------------------------------------------------------------------------
std::vector<Solution> runNormal(const std::vector<Peak>& peaksIn, const InPlaneCand& cand,
                                const Params& prm, Stats* st) {
    std::vector<Peak> peaks;
    for (const Peak& p : peaksIn)
        if (p.qxy > 1e-6 && std::isfinite(p.qxy) && std::isfinite(p.qz)) peaks.push_back(p);
    const int N = (int)peaks.size();
    if (st) {
        st->normalCalls += 1;
        st->peaksUsed = std::max(st->peaksUsed, N);
    }
    if (N < 3 || !(prm.qspec > 0.0)) return {};
    if (traceOn()) bump(kStageCall);

    const double qs = prm.qspec, qs2 = qs * qs;
    const int u = cand.u, v = cand.v;
    const double du = (double)u, dv = (double)v;
    const double A = cand.A, B = cand.B, C = cand.C;
    const double za = std::sqrt(A), zb = std::sqrt(B);
    const double det2d = A * B - C * C;
    if (!(det2d > 0.0)) return {};
    const double guv = A * du * du + B * dv * dv - 2.0 * C * du * dv;
    if (!(qs2 > guv)) return {};
    const double zeta = std::sqrt(qs2 - guv);
    if (!(zeta > 1e-9)) return {};

    const int R = std::max(1, prm.redundant);
    const int L = std::max(1, prm.lLimit);
    const int Ls = std::max(0, std::min(prm.seedL, L));
    if (st) st->seedLUsed = std::max(st->seedLUsed, Ls);

    // R groups of (h,k) per peak → s3 and g_p
    std::vector<std::vector<double>> s3(N), gp(N);
    std::vector<int> sub3(N, 0);
    for (int i = 0; i < N; ++i) {
        const int nr = std::min(R, (int)cand.hk[i].size());
        sub3[i] = nr;
        s3[i].resize(nr);
        gp[i].resize(nr);
        for (int j = 0; j < nr; ++j) {
            const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
            const double s = A * (double)h * du + B * (double)k * dv -
                             C * ((double)h * dv + (double)k * du);
            s3[i][j] = s;
            gp[i][j] = (peaks[i].qz * qs - s) / zeta;
        }
    }

    // Seed pool: peaks in ascending q_z with distinct q_z
    std::vector<int> ord(N);
    for (int i = 0; i < N; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int i, int j) { return peaks[i].qz < peaks[j].qz; });
    std::vector<int> byQz;
    for (int i : ord) {
        bool dup = false;
        for (int p : byQz)
            if (std::fabs(peaks[i].qz - peaks[p].qz) < 0.01) {
                dup = true;
                break;
            }
        if (!dup) byQz.push_back(i);
    }
    if (byQz.size() < 3) return {};

    // Pool = the 6 smallest q_z plus 4 quantiles over the full q_z distribution.
    // Taking only "the 6 smallest q_z" is insufficient: the three seed peaks form the 3×3 system
    //   κh + λk + ρl = gp ,
    // and ρ can only be determined via "some peak with l≠0". On a near-(001) contact plane, reflections with l=0 have q_z determined only by
    // κh+λk, and the peaks with the smallest q_z all have l=0; no combination of them can solve for ρ,
    // so step two yields only trivial solutions (measured on synthetic (001) data: allPeak=0).
    // The pool must contain peaks with clearly larger q_z (l≠0 or larger indices).
    std::vector<int> pick;
    auto pushPool = [&](int i) {
        if (std::find(pick.begin(), pick.end(), i) == pick.end()) pick.push_back(i);
    };
    for (size_t i = 0; i < byQz.size() && i < 6; ++i) pushPool(byQz[i]);
    for (double q : {0.35, 0.55, 0.75, 0.95}) {
        const size_t idx =
            (size_t)std::llround(q * (double)(byQz.size() - 1));
        pushPool(byQz[std::min(idx, byQz.size() - 1)]);
    }
    std::sort(pick.begin(), pick.end(),
              [&](int a, int b) { return peaks[a].qz < peaks[b].qz; });

    // Seed triples:
    //   Group 1: "the three smallest q_z" (historical behavior; tilted contact planes such as DIP rely on it).
    //   Group 2: "greedy full rank" — the (h,k) of the three seeds are pairwise non-collinear. The first two columns of the 3×3 system are the seeds'
    //         (h,k); when the (h,k) of three peaks are collinear the first two columns are linearly dependent and det is always 0.
    //   Group 3+: added in descending q_z span — a large span necessarily includes high q_z peaks, making ρ solvable.
    // Stop only when a group yields a "clearly good" raw (rms ≤ half the threshold): normal data runs only group 1, adding no overhead.
    std::vector<std::array<int, 3>> seedTris;
    seedTris.push_back({byQz[0], byQz[1], byQz[2]});
    {
        std::vector<int> full; // each step requires a nonzero cross product with the (h,k) of the selected seeds (increasing rank step by step)
        for (int i : ord) {
            if (full.size() >= 3) break;
            if (cand.hk[i].empty()) continue;
            const Hkl a = cand.hk[i][0];
            bool indep = true;
            for (int s : full) {
                const Hkl b = cand.hk[s][0];
                if (a.h * b.k - b.h * a.k == 0) {
                    indep = false;
                    break;
                }
            }
            if (indep) full.push_back(i);
        }
        if (full.size() == 3) {
            const std::array<int, 3> t{full[0], full[1], full[2]};
            if (!(t == seedTris[0])) seedTris.push_back(t);
        }
    }
    {
        std::vector<std::array<int, 3>> rest;
        const int np = (int)pick.size();
        for (int i = 0; i < np; ++i)
            for (int j = i + 1; j < np; ++j)
                for (int k = j + 1; k < np; ++k) {
                    const std::array<int, 3> t{pick[i], pick[j], pick[k]};
                    if (std::find(seedTris.begin(), seedTris.end(), t) != seedTris.end()) continue;
                    rest.push_back(t);
                }
        std::sort(rest.begin(), rest.end(), [&](const std::array<int, 3>& a,
                                                const std::array<int, 3>& b) {
            const double sa = peaks[a[2]].qz - peaks[a[0]].qz;
            const double sb = peaks[b[2]].qz - peaks[b[0]].qz;
            if (sa != sb) return sa > sb;
            return a < b;
        });
        const size_t kMaxTri = 4; // total group limit is 6
        for (size_t i = 0; i < rest.size() && i < kMaxTri; ++i) seedTris.push_back(rest[i]);
    }
    if (traceOn()) bump(kStageSeed);

    std::vector<RawSol> raws;
    const bool constrained = (prm.system != CrystalSystem::Triclinic);
    const std::array<bool, 6> freeMask = systemFreeMask(prm.system);

    // Synthesize the unit cell from (κ,λ,ρ) and S3 = sgn·ζ, and analytically recover l per peak:
    //   q_z·q_s = σ(h,k) + (κh+λk+ρl)·S3 ,  S3 = κu+λv+ρw = sgn·ζ
    // The physical parameterization requires ρ = √G33 > 0; first normalize globally to ρ>0 (negating (κ,λ,ρ) together with sgn when necessary).
    // m and −m give the same G*=G0+mmᵀ, and (w,l) are identical after negating sgn as well, so of the original "± two branches"
    // one is necessarily rejected by ρ>0 in the first row and the other yields a solution that exactly duplicates the normalized one; each Cramer solution is tried only once.
    // After normalization sgn = sign(ρ_original); data with S3=−ζ are covered by negating the whole ρ<0 Cramer solution.
    auto trySol = [&](double kap, double lam, double rho) {
        double sgn = 1.0;
        if (rho < 0.0) {
            kap = -kap;
            lam = -lam;
            rho = -rho;
            sgn = -1.0;
        }
        if (!(rho > 1e-9) || !std::isfinite(kap) || !std::isfinite(lam)) return;
        if (traceOn()) bump(kStageRho);
        // Integer consistency of w: the strongest cheap criterion (w is data-derived, with no range limit)
        const double sS3 = sgn * zeta;
        const double wf = (sS3 - du * kap - dv * lam) / rho;
        const int w = (int)std::llround(wf);
        if (!(std::fabs(wf - (double)w) <= 0.35)) return;
        if (prm.fixedUvw && w != prm.fw) return;
        if (!prm.fixedUvw && std::abs(w) > prm.wLimit) return; // w range (default 3)
        if (traceOn()) bump(kStageW);

        // Build the reciprocal metric directly: G* = {A+κ², B+λ², ρ², −C+κλ, κρ, λρ}
        // (mathematically equivalent to recipMetric(cellFromFrame(...)); avoids the redundant and expensive G/nu recomputation in
        //  cellFromFrame + reciprocalMatrix + nuFromUvw).
        const double G0 = A + kap * kap, G1 = B + lam * lam, G2 = rho * rho;
        const double G3 = -C + kap * lam, G4 = kap * rho, G5 = lam * rho;
        // Specular peak self-consistency: cheapest and most selective, moved before cellFromFrame;
        // most spurious solutions are rejected here, avoiding the cost of cellFromFrame (2 acos) for every combination.
        const double wd = (double)w;
        const double qs2p = G0 * du * du + G1 * dv * dv + G2 * wd * wd +
                            2.0 * (G3 * du * dv + G4 * du * wd + G5 * dv * wd);
        if (!(qs2p > 0.0)) return;
        const double qsPred = std::sqrt(qs2p);
        if (std::fabs(qsPred - qs) > prm.dqspec) return;
        if (traceOn()) bump(kStageSpecPre);

        // Construct a Cell only for survivors that pass specular peak self-consistency (needed by restrictCell and the final raw.cell)
        Cell cell;
        if (!cellFromFrame(za, zb, C, kap, lam, rho, cell)) return;
        if (prm.restrictCell) {
            if (cell.c < prm.cLo || cell.c > prm.cHi) return;
            if (cell.al < prm.alLo || cell.al > prm.alHi) return;
            if (cell.be < prm.beLo || cell.be > prm.beHi) return;
        }
        if (traceOn()) bump(kStageCell);

        // nu = G* n / q_spec, computed directly (replacing the double recomputation in reciprocalMatrix + nuFromUvw)
        const std::array<double, 6> G{G0, G1, G2, G3, G4, G5};
        const std::array<double, 3> nu{(G0 * du + G3 * dv + G4 * wd) / qsPred,
                                       (G3 * du + G1 * dv + G5 * wd) / qsPred,
                                       (G4 * du + G5 * dv + G2 * wd) / qsPred};
        std::vector<Hkl> ah((size_t)N);
        double sq = 0.0, sqz = 0.0;
        const double rmsThr2 = prm.dqxyzCutoff * prm.dqxyzCutoff;
        int skipped = 0;
        int used = 0; // number of peaks actually accumulated into the residual (skipped peaks neither add to sq nor enter the denominator)
        for (int i = 0; i < N; ++i) {
            double best = 1e300, bz = 0.0;
            Hkl bh{0, 0, 0};
            for (int j = 0; j < sub3[i]; ++j) {
                const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
                const double lc = (sgn * gp[i][j] - (double)h * kap - (double)k * lam) / rho;
                const int l = (int)std::llround(lc);
                if (std::abs(l) > L) continue;
                double qxp = 0.0, qzp = 0.0;
                forwardNu(G, nu, h, k, l, qxp, qzp);
                const double dx = qxp - peaks[i].qxy;
                const double dz = qzp - peaks[i].qz;
                const double d = dx * dx + dz * dz;
                if (d < best) {
                    best = d;
                    bz = dz * dz;
                    bh = Hkl{h, k, l};
                }
            }
            if (best > 1e299) {
                // All (h,k) candidates of this peak yield |l| > lLimit (lLimit defaults to 6; high-order reflections
                // with high q_z easily exceed it). Discarding the whole candidate would let one or two out-of-range peaks ruin all solutions,
                // so skip this peak instead; still fail when more than 1/4 are skipped, so RMSD does not become meaningless.
                // Leaving (0,0,0) is fine: assignHkl later reassigns all peaks anyway.
                if (++skipped * 4 > N) return;
                ah[(size_t)i] = Hkl{0, 0, 0};
                continue;
            }
            sq += best;
            sqz += bz;
            ++used;
            ah[(size_t)i] = bh;
            // Early residual exit: the denominator still uses the total peak count N (more lenient than used; only for abandoning clearly out-of-spec candidates early,
            // ensuring peaks already skipped cannot wrongly kill solutions that may still pass)
            if (sq > (double)N * rmsThr2) return;
        }
        if (used == 0) return; // no valid peaks → RMSD undefined, fail immediately
        const double rms = std::sqrt(sq / (double)used);
        if (!(rms <= prm.dqxyzCutoff)) return;
        if (traceOn()) bump(kStageRms3);
        if (std::sqrt(sqz / (double)used) > prm.dqzCutoff * 3.0) return;
        if (traceOn()) bump(kStageRmsZ);
        RawSol r;
        r.cell = cell;
        r.w = w;
        r.rms = rms;
        r.hkl = std::move(ah);
        raws.push_back(std::move(r));
    };

    long long trials = 0; // number of step-two 3×3 combinations (maxTrials budget)
    for (size_t ti = 0; ti < seedTris.size(); ++ti) {
        const int pa = seedTris[ti][0], pb = seedTris[ti][1], pc = seedTris[ti][2];
        if (sub3[pa] == 0 || sub3[pb] == 0 || sub3[pc] == 0) continue;
        const size_t rawBase = raws.size();
        for (int i0 = 0; i0 < sub3[pa]; ++i0) {
            const int h1 = cand.hk[pa][i0].h, k1 = cand.hk[pa][i0].k;
            const double g1 = gp[pa][i0];
            for (int i1 = 0; i1 < sub3[pb]; ++i1) {
                const int h2 = cand.hk[pb][i1].h, k2 = cand.hk[pb][i1].k;
                const double g2 = gp[pb][i1];
                for (int i2 = 0; i2 < sub3[pc]; ++i2) {
                    const int h3 = cand.hk[pc][i2].h, k3 = cand.hk[pc][i2].k;
                    const double g3 = gp[pc][i2];
                    for (int l1 = -Ls; l1 <= Ls; ++l1)
                        for (int l2 = -Ls; l2 <= Ls; ++l2)
                            for (int l3 = -Ls; l3 <= Ls; ++l3) {
                                if (prm.maxTrials > 0 && trials >= (long long)prm.maxTrials)
                                    goto trials_done;
                                ++trials;
                                const double c0x = (double)k2 * l3 - (double)l2 * k3;
                                const double c0y = (double)l2 * h3 - (double)h2 * l3;
                                const double c0z = (double)h2 * k3 - (double)k2 * h3;
                                const double det = (double)h1 * c0x + (double)k1 * c0y +
                                                   (double)l1 * c0z;
                                const double sc = (std::fabs((double)h1) + std::fabs((double)k1) +
                                                   std::fabs((double)l1)) *
                                                  (std::fabs((double)h2) + std::fabs((double)k2) +
                                                   std::fabs((double)l2)) *
                                                  (std::fabs((double)h3) + std::fabs((double)k3) +
                                                   std::fabs((double)l3));
                                if (std::fabs(det) < 1e-9 * std::max(1e-12, sc)) continue;
                                const double inv = 1.0 / det;
                                const double kap =
                                    (g1 * c0x - (double)k1 * (g2 * l3 - (double)l2 * g3) +
                                     (double)l1 * (g2 * k3 - (double)k2 * g3)) * inv;
                                const double lam =
                                    ((double)h1 * (g2 * l3 - (double)l2 * g3) -
                                     g1 * ((double)h2 * l3 - (double)l2 * h3) +
                                     (double)l1 * ((double)h2 * g3 - g2 * h3)) * inv;
                                const double rho =
                                    ((double)h1 * ((double)k2 * g3 - g2 * k3) -
                                     (double)k1 * ((double)h2 * g3 - g2 * h3) + g1 * c0z) * inv;
                                // ρ normalized to >0, which determines sgn (see trySol comment)
                                trySol(kap, lam, rho);
                            }
                }
            }
        }
        // Group-stop criterion: stop only when this group's raws are already "clearly good" (≤ half the Δq_xyz threshold).
        // Testing just "raws non-empty" is insufficient — degenerate seed groups can also produce raws whose rms sits right at the threshold,
        // and breaking early would shut out the group that can actually produce good solutions (measured on synthetic (001) data).
        if (raws.size() > rawBase) {
            double best = 1e30;
            for (size_t q = rawBase; q < raws.size(); ++q) best = std::min(best, raws[q].rms);
            if (best <= 0.5 * prm.dqxyzCutoff) break;
        }
    }
trials_done:

    if (st) {
        st->secondTrials += trials; // one Cramer solve per 3×3 seed combination (ρ already normalized, no ± duplication)
        st->rawSol += (long long)raws.size();
    }

    if (raws.empty()) return {};

    // Refinement budget: sort by rms + refine each basin of attraction only once (same as the improved tail)
    std::sort(raws.begin(), raws.end(),
              [](const RawSol& x, const RawSol& y) { return x.rms < y.rms; });
    const int kRefineLimit = std::max(1, std::min(prm.refineTop, 200));
    {
        std::vector<RawSol> picked;
        picked.reserve((size_t)std::min((int)raws.size(), kRefineLimit));
        std::unordered_set<uint64_t> qseen;
        for (const RawSol& r : raws) {
            if (!qseen.insert(rawKeyOf(r.cell, r.w)).second) continue;
            picked.push_back(r);
            if ((int)picked.size() >= kRefineLimit) break;
        }
        raws.swap(picked);
    }

    // Final assignment reuses step one's (h,k) range (default hkLimit2=6, matching the previous hard-coded value):
    // otherwise, after the user increases hkLimit2, step one has candidates but final assignment beyond that range cannot be solved.
    const int L2 = std::max(1, prm.hkLimit2);
    const std::vector<Hkl> cands = hklGrid(L2, L2, L2);
    std::vector<Solution> sols;
    sols.reserve(raws.size());
    std::unordered_set<uint64_t> seen;
    for (const RawSol& r : raws) {
        if (traceOn()) bump(kStageRaw);
        Cell cell = r.cell;
        Hkl uvw0{u, v, r.w};
        Orient orient = Orient::fromUvw(uvw0);
        std::vector<Hkl> ah = r.hkl;
        // |G| upper bound: shrink the assignment grid from 2197 to tens–hundreds
        double maxG = 1e30;
        if (constrained) {
            // Compare axis settings before projection: tetragonal/hexagonal/trigonal have 3 unique axes, orthorhombic 6 axis permutations,
            // monoclinic 3 unique axes (hexagonal also includes the equivalent γ=60° b->b-a setting). Take the setting with the smallest cost (post-projection
            // joint assignment residual); only the winner proceeds to refinement; otherwise, once the correct solution lands in a
            // non-default setting, projection would destroy it.
            // Sort by projection cost first (pure metric computation, very cheap); usually the first is the setting containing the correct solution;
            // run full assignment on only a few settings, avoiding the cost of 6 settings × per-candidate assignment.
            // Two-stage: first score and sort each setting using evaluate (O(peak count), no hkl grid search),
            std::vector<std::pair<double, std::array<int, 9>>> settings;
            for (const std::array<int, 9>& T0 : axisSettings(prm.system, r.cell.ga)) {
                Cell cv;
                double score = 1e30;
                if (transformCell(r.cell, T0, cv)) {
                    applySystemConstraints(prm.system, cv);
                    const Hkl uv = applyT(uvw0, T0);
                    std::vector<Hkl> hk;
                    hk.reserve(ah.size());
                    for (const Hkl& x : ah) hk.push_back(applyT(x, T0));
                    const Quality qv = evaluate(cv, Orient::fromUvw(uv), peaks, hk);
                    if (std::isfinite(qv.rmsQxyz)) score = qv.rmsQxyz + 1e-3 * axisSettingCost(prm.system, r.cell, T0);
                }
                settings.push_back({score, T0});
            }
            std::sort(settings.begin(), settings.end(),
                      [](const std::pair<double, std::array<int, 9>>& x,
                         const std::pair<double, std::array<int, 9>>& y) {
                          return x.first < y.first;
                      });
            double bestRms = 1e30;
            Cell bestCell = cell;
            Hkl bestUvw = uvw0;
            std::vector<Hkl> bestHkl;
            double bestMaxG = 1e30;
            for (const auto& sp : settings) {
                const std::array<int, 9>& T = sp.second;
                if (sp.first >= 1e30) break; // all settings whose transform failed are sorted to the end
                Cell cv;
                if (!transformCell(r.cell, T, cv)) continue;
                applySystemConstraints(prm.system, cv);
                const Hkl uv = applyT(uvw0, T);
                const Orient ov = Orient::fromUvw(uv);
                std::vector<Hkl> hk;
                hk.reserve(ah.size());
                for (const Hkl& x : ah) hk.push_back(applyT(x, T));
                const double mg = assignMaxG(cv, ov, peaks, hk);
                const double rms = assignHkl(cv, ov, peaks, cands, hk, nullptr, 1e30, 1.0, mg);
                if (!std::isfinite(rms)) continue;
                if (rms < bestRms) {
                    bestRms = rms;
                    bestCell = cv;
                    bestUvw = uv;
                    bestHkl = std::move(hk); // move only after a better one is found, avoiding loss of the saved result
                    bestMaxG = mg;
                }
                // Already clearly within tolerance (≤ half tolerance): the other settings cannot overturn this; stop early
                if (bestRms <= 0.5 * prm.dqxyzCutoff) break;
            }
            if (!(bestRms < 1e29)) continue;
            cell = bestCell;
            uvw0 = bestUvw;
            ah = std::move(bestHkl);
            orient = Orient::fromUvw(uvw0);
            maxG = bestMaxG;
        } else {
            const std::array<double, 6> G0 = recipMetric(r.cell);
            const std::array<double, 3> nu0 = orient.get(G0);
            double qmax = 0.0, dmax = 0.0;
            for (size_t i = 0; i < peaks.size(); ++i) {
                double qx = 0.0, qz = 0.0;
                forwardNu(G0, nu0, ah[i].h, ah[i].k, ah[i].l, qx, qz);
                const double dx = qx - peaks[i].qxy, dz = qz - peaks[i].qz;
                dmax = std::max(dmax, std::sqrt(dx * dx + dz * dz));
                qmax = std::max(qmax, std::sqrt(peaks[i].qxy * peaks[i].qxy +
                                                peaks[i].qz * peaks[i].qz));
            }
            maxG = qmax + 3.0 * dmax + 0.05;
        }
        refineCell(peaks, ah, orient, cell, freeMask, 60, prm.system);
        if (constrained) applySystemConstraints(prm.system, cell);
        // Recompute the |G| upper bound from the projected cell and the assigned indices (no longer the pre-projection estimate/full box)
        if (constrained) maxG = assignMaxG(cell, orient, peaks, ah);
        assignHkl(cell, orient, peaks, cands, ah, nullptr, 1e30, 1.0, maxG);
        refineCell(peaks, ah, orient, cell, freeMask, 60, prm.system);
        if (constrained) applySystemConstraints(prm.system, cell);
        if (constrained) maxG = assignMaxG(cell, orient, peaks, ah);
        assignHkl(cell, orient, peaks, cands, ah, nullptr, 1e30, 1.0, maxG);

        // After projection, index aliasing may lock LM in a local minimum: run one progressive reassignment rescue within the constrained subspace.
        // Triggered only when rmsQxyz after refinement is still > 0.3×Δq tolerance (no residual can be recovered once within tolerance).
        if (constrained) {
            maxG = assignMaxG(cell, orient, peaks, ah);
            refineConstrained(peaks, orient, cands, prm.system, cell, ah, maxG,
                              0.3 * prm.dqxyzCutoff, 2.0 * prm.dqxyzCutoff);
            // Monoclinic: canonical β branch (β>90 uses the equivalent c->-c setting, a purely descriptive transform);
            // when all l are 0 and the specular direction w=0, c and β are inseparable, so normalize to β=90 to avoid a β family.
            if (prm.system == CrystalSystem::Monoclinic) {
                monoclinicBranch(cell, uvw0, ah);
                orient = Orient::fromUvw(uvw0);
                bool allL0 = true, hasW = (uvw0.l != 0);
                for (const Hkl& x : ah)
                    if (x.l != 0) { allL0 = false; break; }
                if (allL0 && !hasW) {
                    const Cell keepCell = cell;
                    const Hkl keepUvw = uvw0;
                    monoclinicSnapBeta90(cell, ah);
                    if (!(evaluate(cell, orient, peaks, ah).rmsQxyz <= prm.dqxyzCutoff)) {
                        cell = keepCell;
                        uvw0 = keepUvw;
                    }
                    orient = Orient::fromUvw(uvw0);
                }
            }
        }

        // Crystal-system-constrained mode skips whole-cell reduction: reduceCell changes basis vectors and would rewrite a γ=120° hexagonal cell
        // or an α=γ=90° monoclinic cell into an equivalent form that does not follow the crystal-system standard axis setting.
        Cell rc = cell;
        Hkl ruvw = uvw0;
        std::vector<Hkl> rhkl;
        if (!constrained) {
            const Reduction red = reduceCell(cell);
            if (!red.ok) continue;
            if (traceOn()) bump(kStageReduce);
            rc = red.cell;
            ruvw = applyT(uvw0, red.T);
            rhkl.reserve(ah.size());
            for (const Hkl& x : ah) rhkl.push_back(applyT(x, red.T));
        } else {
            rhkl = ah;
        }

        const Quality q = evaluate(rc, Orient::fromUvw(ruvw), peaks, rhkl);
        if (!(q.rmsQxyz <= prm.dqxyzCutoff)) continue;
        if (traceOn()) bump(kStageRmsQ);
        if (!acceptLatticeIndex(prm.system, rhkl)) continue;
        if (traceOn()) bump(kStageLattice);
        if (q.rmsQxy > prm.dqxyCutoff * 3.0 || q.rmsQz > prm.dqzCutoff * 3.0) continue;
        if (traceOn()) bump(kStageRmsXY);
        if (std::fabs(qSpecOf(reciprocalMatrix(rc), ruvw) - qs) > prm.dqspec) continue;
        if (traceOn()) bump(kStageSpecFin);
        const double vol = cellVolume(rc);
        if (prm.restrictCell) {
            if (rc.a < prm.aLo || rc.a > prm.aHi || rc.b < prm.bLo || rc.b > prm.bHi ||
                rc.c < prm.cLo || rc.c > prm.cHi)
                continue;
            if (vol < prm.volLo || vol > prm.volHi) continue;
        }
        if (traceOn()) bump(kStageRange);
        Solution s;
        s.uvw = ruvw;
        s.nu = Orient::fromUvw(ruvw).get(recipMetric(rc));
        s.cell = rc;
        s.vol = vol;
        s.qspec = qSpecOf(reciprocalMatrix(rc), ruvw);
        s.rmsQxyz = q.rmsQxyz;
        s.rmsQxy = q.rmsQxy;
        s.rmsQz = q.rmsQz;
        s.dqspec = std::fabs(s.qspec - qs);
        s.hasPlane = true; // allow alignSolution to complete specular rod points as m·uvw
        s.hkl = std::move(rhkl);
        // Sign representative normalization: (u,v,w) and (−u,−v,−w) are the same set of lattice planes; unify to the branch with the larger u+v+w
        // (depending on the reduction T, the same physical solution may land in either branch). fixedUvw keeps the user-specified direction.
        if (!prm.fixedUvw) canonicalizeUvw(s);
        s.dist.resize(s.hkl.size());
        for (size_t i = 0; i < s.hkl.size(); ++i) {
            double qx = 0.0, qz = 0.0;
            forwardNu(recipMetric(rc), s.nu, s.hkl[i].h, s.hkl[i].k, s.hkl[i].l, qx, qz);
            s.dist[i] = std::hypot(qx - peaks[i].qxy, qz - peaks[i].qz);
        }
        // Constrained mode deduplicates by a crystal-system-invariant cell fingerprint + face normal (invariants under basis change): equivalent solutions
        // produced by axis permutations/axis settings (e.g. the 6 axis settings of one orthorhombic lattice) are merged; the triclinic path keeps the original key.
        uint64_t key = constrained ? constrainedCellKey(rc, s.uvw) : cellKeyOf(rc);
        if (!constrained) {
            key = mixHash(key, (long long)s.uvw.h * 10000 + (long long)s.uvw.k * 100 + s.uvw.l);
        }
        if (!seen.insert(key).second) continue;
        if (traceOn()) bump(kStageUnique);
        sols.push_back(std::move(s));
        if (traceOn()) bump(kStageSolution);
    }
    std::sort(sols.begin(), sols.end(),
              [](const Solution& x, const Solution& y) { return x.rmsQxyz < y.rmsQxyz; });
    if ((int)sols.size() > prm.maxSolutions) sols.resize((size_t)prm.maxSolutions);
    if (st) st->finalSol += (long long)sols.size();
    return sols;
}

// ---------------------------------------------------------------------------
//  Driver
// ---------------------------------------------------------------------------
std::vector<Solution> runSpecIndex(const std::vector<Peak>& peaksIn, const Params& prm,
                                   const Progress& prog, Stats* st) {
    if (!(prm.qspec > 0.0)) return {};
    const auto t0 = std::chrono::steady_clock::now();
    Stats acc;
    acc.peaksIn = (int)peaksIn.size();
    // Computation input uses |q_xy| uniformly (the sign only reflects the mirror handedness of the detection geometry).
    std::vector<Peak> peaks = peaksIn;
    normalizeQxySigns(peaks);

    // "Clearly good" criterion: RMSD ≤ 8% of the acceptance threshold (leaving a 20% margin below the 10% threshold).
    // If judged directly at ≤10%, a pseudo-supercell just at the line (e.g. aspirin 12 peaks, rms≈0.0048 ≤ 0.005)
    // would be treated as "clearly good", skipping the entire seedL=3 fallback round and losing the true solution (rms≈0.0036,
    // seed |l| up to 3). With the margin, such solutions merely count as "passing" and no longer block the fallback.
    const double clearlyGoodTh = 0.08 * prm.dqxyzCutoff;

    // Step two: run once over the given in-plane candidate list (candidate-level parallelism + early stop)
    auto pass = [&](const std::vector<InPlaneCand>& cands, const Params& p) {
        std::vector<Solution> all;
        const int np = (int)cands.size();
        if (np == 0) return all;
        int nt = (p.threads > 0) ? p.threads : (int)hardwareThreads();
        nt = std::max(1, std::min(nt, np));
        const int batch = std::max(1, nt);
        double bestRms = 1e30;
        int sinceImprove = 0;
        for (int base = 0; base < np; base += batch) {
            if (prog.stop()) break;
            const int end = std::min(np, base + batch);
            {
                char buf[176];
                std::snprintf(buf, sizeof(buf),
                              "SCRM step two: candidate %d/%d   a=%.3f b=%.3f gamma=%.2f  rms=%.4f",
                              base + 1, np, cands[(size_t)base].a, cands[(size_t)base].b,
                              cands[(size_t)base].ga, cands[(size_t)base].rmsQxy);
                prog((double)base / (double)np, buf);
            }
            std::vector<std::vector<Solution>> part((size_t)(end - base));
            std::vector<Stats> partStats((size_t)(end - base));
            parallelForEach(end - base, [&](int k) {
                Params p1 = p;
                p1.threads = 1; // inner level serial, avoiding nested-parallel oversubscription
                part[(size_t)k] = runNormal(peaks, cands[(size_t)(base + k)], p1,
                                            &partStats[(size_t)k]);
            }, nt);

            bool advanced = false;
            for (auto& vec : part)
                for (Solution& s : vec) {
                    if (s.rmsQxyz < bestRms) {
                        bestRms = s.rmsQxyz;
                        advanced = true;
                    }
                    all.push_back(std::move(s));
                }
            for (const Stats& ps : partStats) {
                acc.normalCalls += ps.normalCalls;
                acc.peaksUsed = std::max(acc.peaksUsed, ps.peaksUsed);
                acc.secondTrials += ps.secondTrials;
                acc.rawSol += ps.rawSol;
                acc.seedLUsed = std::max(acc.seedLUsed, ps.seedLUsed);
            }
            sinceImprove = advanced ? 0 : sinceImprove + (end - base);
            if (p.earlyStop) {
                // Stop early only for solutions "clearly better than the acceptance threshold" (same criterion as clearlyGood).
                if (bestRms <= clearlyGoodTh && end >= std::min(np, 160)) break;
                if (!all.empty() && sinceImprove >= 160) break;
            }
        }
        return all;
    };

    // "Clearly good" solution: same criterion as the early stop
    auto clearlyGood = [&](const std::vector<Solution>& v) {
        for (const Solution& s : v)
            if (s.rmsQxyz <= clearlyGoodTh) return true;
        return false;
    };

    // Round 1: tight-threshold in-plane candidates + default seed l radius (fast)
    const std::vector<InPlaneCand> cands = runInPlane(peaks, prm, prog, &acc);
    if (prog.stop()) return {};
    std::vector<Solution> all = pass(cands, prm);
    // Fallback 1: when there is no "clearly good" solution, widen the seed l radius to 3 (±3³=343 combinations).
    // The seeds are the three peaks with the smallest q_z, whose true |l| may reach 2~3 (e.g. (0,0,2) for DIP);
    // the default ±1 would miss solutions for the whole round; but widening makes step two 12× more expensive, so do it only when needed.
    if (!clearlyGood(all) && prm.seedL < 3) {
        Params p2 = prm;
        p2.seedL = 3;
        p2.earlyStop = false;
        std::vector<Solution> wide = pass(cands, p2);
        all.insert(all.end(), std::make_move_iterator(wide.begin()),
                   std::make_move_iterator(wide.end()));
    }
    // Fallback 2: still no "clearly good" solution → relax the in-plane threshold (3 × dqxy), the scan budget and the candidate cap,
    // and rerun step one. This happens only when "nothing was found" (cases such as synthetic/exact data where the candidate list is extremely long
    // and the true value is pushed beyond the hundreds by many small-area degenerate solutions), trading time for recall.
    if (!clearlyGood(all)) {
        Params pw = prm;
        pw.maxPart1 = std::max(prm.maxPart1, 500);
        const std::vector<InPlaneCand> wideC = runInPlaneImpl(peaks, pw, prog, 3.0, 2, &acc);
        acc.wideFallback = true;
        if (!wideC.empty()) {
            Params p3 = pw;
            p3.seedL = 3;
            p3.earlyStop = false;
            std::vector<Solution> wide = pass(wideC, p3);
            all.insert(all.end(), std::make_move_iterator(wide.begin()),
                       std::make_move_iterator(wide.end()));
        }
    }

    std::sort(all.begin(), all.end(),
              [](const Solution& a, const Solution& b) { return a.rmsQxyz < b.rmsQxyz; });
    std::vector<Solution> out;
    std::unordered_set<uint64_t> seen;
    const bool constrained = (prm.system != CrystalSystem::Triclinic);
    for (Solution& s : all) {
        if (!prm.fixedUvw) canonicalizeUvw(s); // idempotent: ensures a unique representative after merging multiple rounds
        uint64_t key = constrained ? constrainedCellKey(s.cell, s.uvw) : cellKeyOf(s.cell);
        if (!constrained)
            key = mixHash(key, (long long)s.uvw.h * 10000 + (long long)s.uvw.k * 100 + s.uvw.l);
        if (!seen.insert(key).second) continue;
        alignSolution(peaks, s); // complete specular rod points as m·uvw; expand hkl/dist to the input peak count
        out.push_back(std::move(s));
        if ((int)out.size() >= prm.maxSolutions) break;
    }
    funnelDump();
    if (st) {
        acc.finalSol = (long long)out.size();
        acc.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        *st = acc;
    }
    return out;
}

Params fromSpecParams(const SpecParams& sp) {
    Params p;
    p.qspec = sp.qspec;
    p.system = sp.system;
    p.uvLimit = sp.uvLimit;
    p.fixedUvw = sp.fixedUvw;
    p.fu = sp.fu;
    p.fv = sp.fv;
    p.fw = sp.fw;
    p.hkLimit1 = sp.hkLimit1;
    p.hkLimit2 = sp.hkLimit2;
    p.seed3Half = sp.seed3Half;
    p.linesToPermute = sp.linesToPermute;
    p.redundant = sp.redundant;
    p.maxPart1 = sp.maxPart1;
    p.lLimit = sp.lLimit;
    p.wLimit = sp.wLimit;
    p.refineTop = sp.refineTop;
    p.maxSolutions = sp.maxSolutions;
    p.earlyStop = sp.earlyStop;
    p.dqxyCutoff = sp.dqxyCutoff;
    p.dqzCutoff = sp.dqzCutoff;
    p.dqxyzCutoff = sp.dqxyzCutoff;
    p.dqspec = sp.dqspec;
    p.restrictCell = sp.restrictCell;
    p.aLo = sp.aLo;
    p.aHi = sp.aHi;
    p.bLo = sp.bLo;
    p.bHi = sp.bHi;
    p.cLo = sp.cLo;
    p.cHi = sp.cHi;
    p.alLo = sp.alLo;
    p.alHi = sp.alHi;
    p.beLo = sp.beLo;
    p.beHi = sp.beHi;
    p.gaLo = sp.gaLo;
    p.gaHi = sp.gaHi;
    p.volLo = sp.volLo;
    p.volHi = sp.volHi;
    p.niggli = sp.niggli;
    p.threads = sp.threads;
    return p;
}

} // namespace scrm
} // namespace gid
