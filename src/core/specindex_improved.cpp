// ============================================================================
//  specindex_improved.cpp — improved two-step method (Part ONE + Part TWO + refinement)
//
//  Code comes from index2_other/src/core/specindex.cpp, ported wholesale into gid::improved
//  namespace (algorithm unchanged), coexisting with the original port (part1/part2) in specindex.cpp.
// ============================================================================
#include "specindex.h"
#include "linalg.h"
#include "parallel.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gid {
namespace improved {
namespace {

constexpr double PI2 = 6.28318530717958647692;

inline double clamp1(double x) { return x > 1.0 ? 1.0 : (x < -1.0 ? -1.0 : x); }

// ---------------------------------------------------------------------------
//  Part ONE
// ---------------------------------------------------------------------------
struct AbgRaw {
    int u = 0, v = 0;
    double a = 0, b = 0, ga = 0;
    double area = 0; // in-plane parallelogram area a b sin γ (computed once on insert; the sort comparator no longer calls sin)
};

// Quantized dedup key (tolerances 0.15 / 0.15 / 0.3, matching MATLAB uniquetol)
inline uint64_t abgKey(int u, int v, double a, double b, double ga) {
    auto q = [](double x, double t) -> long long { return (long long)std::llround(x / t); };
    const long long ka = q(a, 0.15), kb = q(b, 0.15), kg = q(ga, 0.3);
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](long long x) {
        h ^= (uint64_t)(x + 0x9e3779b97f4a7c15LL);
        h *= 1099511628211ull;
    };
    mix(u); mix(v); mix(ka); mix(kb); mix(kg);
    return h;
}

// Full unit-cell dedup key (0.02 Å / 0.1°)
inline uint64_t cellKey(const Cell& c) {
    auto q = [](double x, double t) -> long long { return (long long)std::llround(x / t); };
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](long long x) {
        h ^= (uint64_t)(x + 0x9e3779b97f4a7c15LL);
        h *= 1099511628211ull;
    };
    mix(q(c.a, 0.02)); mix(q(c.b, 0.02)); mix(q(c.c, 0.02));
    mix(q(c.al, 0.1)); mix(q(c.be, 0.1)); mix(q(c.ga, 0.1));
    return h;
}

// Coarse dedup key before refinement (0.05 Å / 0.5°) + w: raw solutions in the same basin of attraction are refined only once
inline uint64_t rawKey(const Cell& c, int w) {
    auto q = [](double x, double t) -> long long { return (long long)std::llround(x / t); };
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](long long x) {
        h ^= (uint64_t)(x + 0x9e3779b97f4a7c15LL);
        h *= 1099511628211ull;
    };
    mix(q(c.a, 0.25)); mix(q(c.b, 0.25)); mix(q(c.c, 0.25));
    mix(q(c.al, 2.0)); mix(q(c.be, 2.0)); mix(q(c.ga, 2.0));
    mix(w);
    return h;
}

// Fixed-size 4x3 Householder least squares (hot spot of Part TWO: ~7e5 calls per candidate).
// Numerically equivalent to lstsqQR in linalg.cpp, but uses only stack arrays, avoiding 3 dynamic matrix allocations per call.
inline bool lstsq4x3(const double A[4][3], const double b[4], double x[3]) {
    double R[4][3];
    double B[4];
    for (int i = 0; i < 4; ++i) {
        B[i] = b[i];
        R[i][0] = A[i][0];
        R[i][1] = A[i][1];
        R[i][2] = A[i][2];
    }
    for (int k = 0; k < 3; ++k) {
        double nrm2 = 0.0;
        for (int i = k; i < 4; ++i) nrm2 += R[i][k] * R[i][k];
        const double nrm = std::sqrt(nrm2);
        if (nrm < 1e-300) return false;
        const double alpha = (R[k][k] >= 0.0) ? -nrm : nrm;
        double v[4] = {0.0, 0.0, 0.0, 0.0};
        for (int i = k; i < 4; ++i) v[i] = R[i][k];
        v[k] -= alpha;
        double vtv = 0.0;
        for (int i = k; i < 4; ++i) vtv += v[i] * v[i];
        if (vtv < 1e-300) return false;
        for (int j = k; j < 3; ++j) {
            double dot = 0.0;
            for (int i = k; i < 4; ++i) dot += v[i] * R[i][j];
            const double f = 2.0 * dot / vtv;
            for (int i = k; i < 4; ++i) R[i][j] -= f * v[i];
        }
        double dotb = 0.0;
        for (int i = k; i < 4; ++i) dotb += v[i] * B[i];
        const double fb = 2.0 * dotb / vtv;
        for (int i = k; i < 4; ++i) B[i] -= fb * v[i];
    }
    for (int i = 2; i >= 0; --i) {
        double s = B[i];
        for (int j = i + 1; j < 3; ++j) s -= R[i][j] * x[j];
        if (std::fabs(R[i][i]) < 1e-300) return false;
        x[i] = s / R[i][i];
    }
    return true;
}

} // namespace

std::vector<AbgCand> runPartOne(const std::vector<Peak>& peaksIn, const SpecParams& prm,
                                const Progress& prog) {
    // Points on the specular rod (q_xy ≈ 0) take no part in the in-plane fit (they are described by q_spec)
    std::vector<Peak> peaks;
    for (const Peak& p : peaksIn)
        if (p.qxy > 1e-6 && std::isfinite(p.qxy) && std::isfinite(p.qz)) peaks.push_back(p);
    const int N = (int)peaks.size();
    if (N < 3 || !(prm.qspec > 0.0)) return {};

    const double qs = prm.qspec, qs2 = qs * qs;
    const int L1 = std::max(1, prm.hkLimit1);
    const int W1 = 2 * L1 + 1;
    const int nTri = std::min(N, std::max(3, prm.linesToPermute));

    // Starting three reflections: both seed sets are taken, each expanded as "any 3 of the first
    //   nTri peaks": (A) ascending |q| = sqrt(qxy^2+qz^2) (matches the original
    //       function_INITIALIZE_GUI, deduplicated with uniquetol, absolute tolerance 5e-4);
    //   (B) ascending q_xy (the approach the improved version originally used).
    // Using only one ordering loses solutions: the true DIP solution lies in the |q|-sorted seeds
    // (the four peaks with the smallest q_xy are almost coplanar), while synthetic data with many small-angle peaks conversely needs the q_xy-sorted seeds.
    std::vector<double> qlen(N);
    for (int i = 0; i < N; ++i)
        qlen[i] = std::sqrt(peaks[i].qxy * peaks[i].qxy + peaks[i].qz * peaks[i].qz);
    const double qTol = 5e-4; // equivalent to the original tol_diff_qxym/qmax * max|q|

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

    // Note: triples are **not sorted by peak index**; element order keeps the "seed order" (sel[i] is
    // the one with the smallest |q| or q_xy among the three), because the "first-seed nearest-neighbour box"
    // of Part ONE tightens only the (h,k) box of tr[0]; dedup compares only sorted copies, so the two seed sets cannot give the same triple.
    std::vector<std::array<int, 3>> triples;
    std::vector<std::array<int, 3>> triplesSorted; // used only for deduplication
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

    // (u,v) ∈ [-L,L]^2; (0,0) must be kept -- it corresponds to the most common (001) contact plane.
    // (u,v) and (-u,-v) give exactly the same (a,b,gamma) solution set (negate h,k simultaneously; see
    // the 3x3 coefficient matrix of Part ONE, which depends on H,K only through products), so only the half-plane is enumerated, halving the candidates and the subsequent Part TWO work.
    std::vector<std::pair<int, int>> uvs;
    const int UL = std::max(0, prm.uvLimit);
    if (prm.fixedUvw) {
        // Fixed contact plane: search only the given (u,v)
        uvs.push_back({prm.fu, prm.fv});
    } else {
        for (int u = -UL; u <= UL; ++u)
            for (int v = -UL; v <= UL; ++v) {
                if (u < 0) continue;          // half-plane
                if (u == 0 && v < 0) continue; // boundary: (0,v) and (0,-v) are equivalent
                uvs.push_back({u, v});
            }
    }
    if (uvs.empty()) uvs.push_back({0, 0});

    // ------------------------------------------------------------------
    //  Row table: row[peak][(h,k)] = (m0, m1, m2), one table per (u,v), all built up front.
    //  Previously built on demand inside the parallel block; now every (u,v)×triple task uses them directly.
    // ------------------------------------------------------------------
    const int W2 = W1 * W1;
    std::vector<std::vector<double>> tabs(uvs.size());
    for (size_t ui = 0; ui < uvs.size(); ++ui) {
        const double du = (double)uvs[ui].first, dv = (double)uvs[ui].second;
        std::vector<double>& tab = tabs[ui];
        tab.assign((size_t)N * W2 * 3, 0.0);
        for (int i = 0; i < N; ++i) {
            const double Q = peaks[i].qxy, Z = peaks[i].qz;
            if (!(Q > 1e-9)) continue; // row stays all zero -> the det criterion rejects it naturally
            const double iQ2 = 1.0 / (Q * Q);
            const double uz = du * Z / qs, vz = dv * Z / qs;
            double* row = tab.data() + (size_t)i * W2 * 3;
            for (int h = -L1; h <= L1; ++h)
                for (int k = -L1; k <= L1; ++k) {
                    const double H = h - uz, K = k - vz;
                    double* o = row + ((size_t)((h + L1) * W1 + (k + L1))) * 3;
                    o[0] = H * H * iQ2 + du * du / qs2;
                    o[1] = K * K * iQ2 + dv * dv / qs2;
                    o[2] = -2.0 * (H * K * iQ2 + du * dv / qs2);
                }
        }
    }

    // Task = (u,v) × triple. Previously chunked only by (u,v), so with a small `uvLimit` (default 1 -> 6 chunks)
    // the cores were underused; distributing by task keeps every core busy.
    struct Task {
        int ui;
        int tri;
    };
    std::vector<Task> tasks;
    tasks.reserve(uvs.size() * triples.size());
    for (size_t ui = 0; ui < uvs.size(); ++ui)
        for (size_t t = 0; t < triples.size(); ++t) tasks.push_back(Task{(int)ui, (int)t});

    std::mutex mu;
    std::vector<AbgRaw> out;
    std::unordered_set<uint64_t> seen;
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

        std::vector<AbgRaw> local;
        std::unordered_set<uint64_t> localSeen; // in-place dedup to limit memory
        // A 6-fold loop replaces the old "one flat loop + 6 integer divmods each iteration" (integer division took
        // nearly half of the inner-loop time) and hoists quantities independent of the third peak outward.
        // (h,k)=(0,0) is allowed: (00l) reflections on the specular rod are common starting reflections and the
        // original 6-tuple also includes 0; when u=v=0 that row is always zero, so the det criterion rejects it naturally.
        // First-seed nearest-neighbour box: the seed peak with the shortest indices in the triple (tr[0], seeds chosen
        // by ascending |q| / q_xy) needs (h,k) enumerated only within [-LA, LA]². With seed3Half=0, LA=L1, matching the original
        // implementation exactly; with LA=1 the inner iteration count drops 5.4x (7²→3²).
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
                        const double s1 = std::fabs(m1x) + std::fabs(m1y) + std::fabs(m1z);
                        // c2 = m0 x m1 depends only on (h1,k1,h2,k2), hoisted out of the third-peak loop
                        const double c2x = m0y * m1z - m0z * m1y;
                        const double c2y = m0z * m1x - m0x * m1z;
                        const double c2z = m0x * m1y - m0y * m1x;
                        const double s01 = s0 * s1;
                        for (int h3 = -L1; h3 <= L1; ++h3)
                            for (int k3 = -L1; k3 <= L1; ++k3) {
                                const double* m2 = rowC + ((size_t)((h3 + L1) * W1 + (k3 + L1))) * 3;
                                const double m2x = m2[0], m2y = m2[1], m2z = m2[2];
                                // det = m0 . (m1 x m2)
                                const double c0x = m1y * m2z - m1z * m2y;
                                const double c0y = m1z * m2x - m1x * m2z;
                                const double c0z = m1x * m2y - m1y * m2x;
                                const double det = m0x * c0x + m0y * c0y + m0z * c0z;
                                const double s2 = std::fabs(m2x) + std::fabs(m2y) + std::fabs(m2z);
                                if (std::fabs(det) < 1e-9 * std::max(1e-12, s01 * s2)) continue;
                                // x = M^-1 * [1,1,1]; the three columns of adj = m1xm2, m2xm0, m0xm1
                                const double c1x = m2y * m0z - m2z * m0y;
                                const double c1y = m2z * m0x - m2x * m0z;
                                const double c1z = m2x * m0y - m2y * m0x;
                                const double inv = 1.0 / det;
                                const double X = (c0x + c1x + c2x) * inv;
                                const double Y = (c0y + c1y + c2y) * inv;
                                const double Zc = (c0z + c1z + c2z) * inv;
                                // den = q_s^2 - g_uv is an algebraic identity (the Lam terms cancel
                                // when g_uv is expanded), so the original "SUB2: q_s^2 > g_uv" is equivalent to
                                // den > 0, so g_uv / Cz need not be computed separately.
                                const double den = qs2 - du2 * X - dv2 * Y + 2.0 * duv * Zc;
                                if (!(den > 1e-12)) continue;
                                const double Lam = (X * Y - Zc * Zc) / den;
                                const double za2 = X + dv2 * Lam;
                                const double zb2 = Y + du2 * Lam;
                                if (!(za2 > 1e-12) || !(zb2 > 1e-12)) continue;
                                // Cheap equivalent form of 2D Niggli reduction (in the squared domain, avoiding sqrt/acos):
                                //   a <= b          <=> za >= zb  <=> za2 >= zb2
                                //   |b cos g| <= a/2 <=> 2|K| <= zb2, where K = za zb cos g
                                // A little slack is allowed here for a "superset" filter, and survivors are still judged
                                // exactly by niggli2D(a,b,ga), so the candidate set is unchanged.
                                const double K = Zc + duv * Lam;
                                if (prm.niggli) {
                                    if (za2 < zb2 * (1.0 - 1e-9)) continue;
                                    if (2.0 * std::fabs(K) > zb2 * (1.0 + 1e-6)) continue;
                                }
                                const double za = std::sqrt(za2), zb = std::sqrt(zb2);
                                const double cg = clamp1(K / (za * zb));
                                const double sg = std::sqrt(std::max(1e-12, 1.0 - cg * cg));
                                const double a = PI2 / (za * sg), b = PI2 / (zb * sg);
                                if (!std::isfinite(a) || !std::isfinite(b) || a <= 0.0 || b <= 0.0)
                                    continue;
                                const double ga = std::acos(cg) * 180.0 / M_PI;
                                if (!std::isfinite(ga)) continue;
                                if (prm.restrictCell) {
                                    if (a < prm.aLo || a > prm.aHi) continue;
                                    if (b < prm.bLo || b > prm.bHi) continue;
                                    if (ga < prm.gaLo || ga > prm.gaHi) continue;
                                }
                                if (prm.niggli && !niggli2D(a, b, ga)) continue;
                                const uint64_t key = abgKey(u, v, a, b, ga);
                                const double area = a * b * std::sin(ga * M_PI / 180.0);
                                if (localSeen.insert(key).second)
                                    local.push_back(AbgRaw{u, v, a, b, ga, area});
                            }
                    }
            }
        if (prog.report)
            prog((double)(++done) / (double)tasks.size() * 0.9, "Part 1: enumerating (u,v) candidates");
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(mu);
            if ((int)out.size() < 400000) out.insert(out.end(), local.begin(), local.end());
        }
    }, prm.threads);

    if (prog.stop()) return {};

    // Deduplication (post-pass, A9)
    std::vector<AbgRaw> uniq;
    uniq.reserve(std::min<size_t>(out.size(), 20000));
    for (const AbgRaw& r : out) {
        const uint64_t key = abgKey(r.u, r.v, r.a, r.b, r.ga);
        if (seen.insert(key).second) uniq.push_back(r);
    }
    prog(0.92, "Part 1: indexing all peaks");

    // area = a b sin γ is precomputed on insert (once per candidate). The old comparator called sin
    // once per element on every comparison, i.e. hundreds of thousands of candidates × ~2×10⁷ comparisons = ~10⁷ sin calls, the largest
    // single item in Part ONE (~38% of Part ONE at 265k candidates in measurements). The key values are exactly as computed before, and
    // the sort order is bit-identical (in the squared domain det2 = A·B−C² is order-equivalent to the area, but sqrt/acos round-trips
    // would reorder nearly tied candidates and thereby shift the maxPart1 truncation boundary, so the original area key is kept).
    std::sort(uniq.begin(), uniq.end(), [&](const AbgRaw& x, const AbgRaw& y) {
        if (x.area != y.area) return x.area < y.area;
        if (x.a != y.a) return x.a < y.a;
        if (x.b != y.b) return x.b < y.b;
        return x.ga < y.ga;
    });

    // For each candidate take the R nearest (h,k) per peak and compute RMSD(q_xy).
    // Walk in area order and stop once maxPart1 are collected -- the result matches "compute all, then sort and truncate",
    // but avoids a full (h,k) table scan for tens of thousands of candidates (most of which are truncated right after).
    const int R = std::max(1, prm.redundant);
    const int L2 = std::max(1, prm.hkLimit2);
    std::vector<AbgCand> cands;
    cands.reserve(std::min<size_t>(uniq.size(), (size_t)std::max(1, prm.maxPart1) + 1));
    for (const AbgRaw& r : uniq) {
        if ((int)cands.size() >= prm.maxPart1) break;
        double za = 0, zb = 0, C = 0;
        abgToZ(r.a, r.b, r.ga, za, zb, C);
        const double za2 = za * za, zb2 = zb * zb;
        const double du = (double)r.u, dv = (double)r.v;
        const double guv = du * du * za2 + dv * dv * zb2 - 2.0 * du * dv * C;
        const double den = 1.0 - guv / qs2;
        if (!(den > 1e-9) || !(qs2 > guv)) continue;
        const double zeta = std::sqrt(qs2 - guv);

        AbgCand c;
        c.u = r.u; c.v = r.v; c.a = r.a; c.b = r.b; c.ga = r.ga;
        c.za = za; c.zb = zb; c.C = C; c.guv = guv; c.zeta = zeta;
        c.hk.assign(N, {});
        double sq = 0.0;
        for (int i = 0; i < N; ++i) {
            const double Q = peaks[i].qxy, Z = peaks[i].qz;
            const double uz = du * Z / qs, vz = dv * Z / qs;
            struct Item { double d; Hkl hk; };
            std::vector<Item> items;
            items.reserve((2 * L2 + 1) * (2 * L2 + 1));
            for (int h = -L2; h <= L2; ++h)
                for (int k = -L2; k <= L2; ++k) {
                    // (h,k)=(0,0) is the (00l) reflection on the specular rod and must take part in indexing;
                    // otherwise these peaks get only poor alternative (h,k) for dq_xy and the whole candidate is rejected by RMSD.
                    const double H = h - uz, K = k - vz;
                    const double rod = (h * dv - k * du);
                    const double num = za2 * H * H + zb2 * K * K - 2.0 * C * H * K -
                                       (za2 * zb2 - C * C) * rod * rod / qs2;
                    if (num < 0.0) continue;
                    const double qp = std::sqrt(num / den);
                    items.push_back({std::fabs(qp - Q), Hkl{h, k, 0}});
                }
            std::sort(items.begin(), items.end(),
                      [](const Item& x, const Item& y) { return x.d < y.d; });
            const int nr = std::min(R, (int)items.size());
            for (int j = 0; j < nr; ++j) c.hk[i].push_back(items[j].hk);
            if (nr > 0) sq += items[0].d * items[0].d;
        }
        c.rmsQxy = std::sqrt(sq / (double)N);
        if (c.rmsQxy <= prm.dqxyCutoff * 3.0 + 1e-9) cands.push_back(std::move(c));
    }

    prog(1.0, "Part 1 complete");
    return cands;
}

// ---------------------------------------------------------------------------
//  Part TWO
// ---------------------------------------------------------------------------
std::vector<Solution> runPartTwo(const std::vector<Peak>& peaksIn, const AbgCand& cand,
                                 const SpecParams& prm) {
    std::vector<Peak> peaks;
    for (size_t i = 0; i < peaksIn.size(); ++i) {
        const Peak& p = peaksIn[i];
        if (p.qxy > 1e-6 && std::isfinite(p.qxy) && std::isfinite(p.qz)) peaks.push_back(p);
    }
    const int N = (int)peaks.size();
    if (N < 3 || !(prm.qspec > 0.0)) return {};

    const double qs = prm.qspec, qs2 = qs * qs;
    // The crystal-system constraint applies only during refinement: Part ONE / the second-step solve still use a triclinic parametrization,
    // and the output unit cell is projected and locked by applySystemConstraints + systemFreeMask.
    const bool constrained = (prm.system != CrystalSystem::Triclinic);
    const double za = cand.za, zb = cand.zb, C = cand.C, zeta = cand.zeta;
    const double za2 = za * za, zb2 = zb * zb;
    const double denom = std::sqrt(std::max(1e-18, za2 * zb2 - C * C)); // za zb sin gamma
    const int R = std::max(1, prm.redundant);
    const int u = cand.u, v = cand.v;

    // Per-peak, per-candidate sub3 / g_z / P (squared in-plane modulus)
    std::vector<std::vector<double>> sub3(N), gz(N), pla(N);
    for (int i = 0; i < N; ++i) {
        const int nr = std::min(R, (int)cand.hk[i].size());
        sub3[i].resize(nr);
        gz[i].resize(nr);
        pla[i].resize(nr);
        for (int j = 0; j < nr; ++j) {
            const double h = cand.hk[i][j].h, k = cand.hk[i][j].k;
            const double s3 = h * u * za2 + k * v * zb2 - (h * v + k * u) * C;
            sub3[i][j] = s3;
            gz[i][j] = (peaks[i].qz * qs - s3) / zeta;
            pla[i][j] = h * h * za2 + k * k * zb2 - 2.0 * h * k * C;
        }
    }

    // Pick 3 peaks with the smallest q_z that are mutually distinct (0.01)
    std::vector<int> ord(N);
    for (int i = 0; i < N; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(),
              [&](int i, int j) { return peaks[i].qz < peaks[j].qz; });
    std::vector<int> pick;
    for (int i : ord) {
        bool dup = false;
        for (int p : pick)
            if (std::fabs(peaks[i].qz - peaks[p].qz) < 0.01) { dup = true; break; }
        if (!dup) pick.push_back(i);
        if ((int)pick.size() == 3) break;
    }
    if ((int)pick.size() < 3) return {};

    const int L = std::max(1, prm.lLimit);
    const int NL = 2 * L + 1;
    const int W = std::max(0, prm.wLimit);
    // Fixed contact plane: w takes only the given value
    const int wLo = prm.fixedUvw ? prm.fw : -W;
    const int wHi = prm.fixedUvw ? prm.fw : W;
    const int nrCombo = (int)std::pow((double)R, 3.0);

    // Scoring: returns rms_qxyz; optionally outputs the best (h,k,l) per peak
    auto score = [&](double kap, double lam, double rho, int nUse, std::vector<Hkl>* outHkl,
                     double* rmsQz, double* rmsQxy, double cutoff) -> double {
        double sq = 0.0, sqz = 0.0, sqx = 0.0;
        const int n = std::min(nUse, N);
        // Saving hkl is unnecessary when only the first 5 peaks are used for the "validity pre-filter" (1e3~1e5 calls per candidate),
        // so when outHkl == nullptr no buffer is allocated.
        if (outHkl) outHkl->assign((size_t)n, Hkl{0, 0, 0});
        const double lim = (cutoff < 1e29) ? cutoff * cutoff * (double)N : 1e300;
        for (int i = 0; i < n; ++i) {
            double best = 1e300, bz = 0.0, bx = 0.0;
            Hkl bh{0, 0, 0};
            const double gm = std::sqrt(peaks[i].qxy * peaks[i].qxy + peaks[i].qz * peaks[i].qz);
            const int nr = (int)sub3[i].size();
            for (int j = 0; j < nr; ++j) {
                const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
                const double lc = (gz[i][j] - kap * h - lam * k) / rho;
                const double lr = std::round(lc);
                // l is rounded to the nearest integer and judged by its residual; |lc-lr| < 1e-4 is no longer required:
                // measured peak positions (≈3 decimals) and the ρ error of the 3-peak fit both push lc away from an integer
                // by ~1e-2, so a hard gate would mark every peak as unindexed and no solution could pass the scoring.
                if (std::fabs(lr) > (double)L) continue;
                const double gp = kap * h + lam * k + rho * lr;
                const double qzp = (gp * zeta + sub3[i][j]) / qs;
                const double gpred = std::sqrt(std::max(0.0, pla[i][j] + gp * gp));
                const double dqz = peaks[i].qz - qzp;
                const double dqg = gpred - gm;
                const double d = dqg * dqg;
                if (d < best) {
                    best = d;
                    bz = dqz * dqz;
                    const double qxyp = std::sqrt(std::max(0.0, gpred * gpred - qzp * qzp));
                    const double dqx = qxyp - peaks[i].qxy;
                    bx = dqx * dqx;
                    bh = Hkl{h, k, (int)lr};
                }
            }
            if (best > 1e299) { best = cutoff * cutoff * 100.0; bz = best; bx = best; }
            sq += best; sqz += bz; sqx += bx;
            if (outHkl) (*outHkl)[(size_t)i] = bh;
            // Early exit: the residual already exceeds the caller-supplied threshold, so the remaining peaks need not be computed
            if (sq > lim) return std::sqrt(sq / (double)n) + 1.0;
        }
        if (rmsQz) *rmsQz = std::sqrt(sqz / (double)n);
        if (rmsQxy) *rmsQxy = std::sqrt(sqx / (double)n);
        return std::sqrt(sq / (double)n);
    };

    struct Raw {
        Cell cell;
        int w = 0;
        double rms = 0.0;
        std::vector<Hkl> hkl;
    };
    std::mutex mu;
    std::vector<Raw> raws;
    const int nCheck = std::min(N, 5);
    // Weak-validation threshold: the sum of squared residuals over the first nCheck peaks is a lower bound on the full sum, so if the full RMS is
    // <= dqxyzCutoff, the RMS over the first nCheck peaks is necessarily <= sqrt(N/nCheck)*dqxyzCutoff.
    // max(4, sqrt(N/5)+0.5) guarantees a "superset of the necessary condition" for any peak count, so nothing is wrongly killed.
    const double loose =
        std::max(4.0, std::sqrt((double)N / (double)std::max(1, nCheck)) + 0.5) * prm.dqxyzCutoff;

    // Generate one raw solution from (κ,λ,ρ,w): the l-triple enumeration and difference-vector voting paths share the same criteria
    auto emit = [&](double kap, double lam, double rho, int w, std::vector<Raw>& local) {
        if (!(rho > 1e-9)) return;
        // Closed-form inversion (wiki §3.4)
        const double D2 = kap * kap * zb2 + 2.0 * kap * lam * C + lam * lam * za2 + za2 * zb2 - C * C;
        if (!(D2 > 0.0)) return;
        const double D = std::sqrt(D2);
        const double c = PI2 * D / (rho * denom);
        if (!(c > 0.0) || !std::isfinite(c)) return;
        // Cheap checks first (A5: c first)
        if (prm.restrictCell && (c < prm.cLo || c > prm.cHi)) return;
        // Specular-peak self-consistency: q_spec(solution) = sqrt(g_uv + ((u,v,w)·(κ,λ,ρ))²). This check is the cheapest and the most selective
        // (dqspec is usually 0.01), so it is placed before acos/Niggli/volume.
        const double gpar = (double)u * kap + (double)v * lam + (double)w * rho;
        const double gsp = std::sqrt(std::max(0.0, cand.guv + gpar * gpar));
        if (std::fabs(gsp - qs) > prm.dqspec) return;
        // Then weak validation (A6): residuals of the first nCheck peaks. Full scoring can pass only if this passes,
        // so it is hoisted before Niggli / volume / acos, where the vast majority of spurious solutions are eliminated.
        const double rms0 = score(kap, lam, rho, nCheck, nullptr, nullptr, nullptr, loose);
        if (rms0 > loose) return;
        const double cosAl = clamp1(-(kap * C / za + lam * za) / D);
        const double cosBe = clamp1(-(kap * zb + lam * C / zb) / D);
        const double al = std::acos(cosAl) * 180.0 / M_PI;
        const double be = std::acos(cosBe) * 180.0 / M_PI;
        if (!std::isfinite(al) || !std::isfinite(be)) return;
        if (prm.restrictCell) {
            if (al < prm.alLo || al > prm.alHi) return;
            if (be < prm.beLo || be > prm.beHi) return;
        }
        Cell cell;
        zToAbg(za, zb, C, cell.a, cell.b, cell.ga);
        cell.c = c; cell.al = al; cell.be = be;
        const double vol = cellVolume(cell);
        if (!(vol > 0.0)) return;
        if (prm.restrictCell && (vol < prm.volLo || vol > prm.volHi)) return;
        // 3D Niggli is a "reduced cell" criterion: high-symmetry cells (e.g. a hexagonal cell with γ=120°) are not necessarily
        // reduced, and reduceCell is not run in constrained mode, so this step is skipped to avoid killing true solutions.
        if (prm.niggli && !constrained) {
            if (cell.b * cell.b > cell.c * cell.c) return;
            if (!niggli3D(cell)) return;
        }
        std::vector<Hkl> ah;
        double rz = 0, rx = 0;
        const double rms = score(kap, lam, rho, N, &ah, &rz, &rx, 1e30);
        if (rms > prm.dqxyzCutoff) return;
        if (rz > prm.dqzCutoff * 3.0) return;
        Raw r;
        r.cell = cell;
        r.w = w;
        r.rms = rms;
        r.hkl = ah;
        local.push_back(std::move(r));
    };

    // ---------------- Path A: l-triple enumeration (original to the improved version) ----------------
    auto runEnum = [&]() {
    parallelFor(0, NL * NL * NL, [&](int b, int e) {
        std::vector<Raw> local;
        for (int t = b; t < e; ++t) {
            int x = t;
            const int l3 = x % NL - L; x /= NL;
            const int l2 = x % NL - L; x /= NL;
            const int l1 = x % NL - L;
            const int ls[3] = {l1, l2, l3};
            for (int ci = 0; ci < nrCombo; ++ci) {
                int y = ci;
                const int r3 = y % R; y /= R;
                const int r2 = y % R; y /= R;
                const int r1 = y % R;
                const int rs[3] = {r1, r2, r3};
                if (rs[0] >= (int)sub3[pick[0]].size() || rs[1] >= (int)sub3[pick[1]].size() ||
                    rs[2] >= (int)sub3[pick[2]].size())
                    continue;
                const int hs[3] = {cand.hk[pick[0]][rs[0]].h, cand.hk[pick[1]][rs[1]].h,
                                   cand.hk[pick[2]][rs[2]].h};
                const int ks[3] = {cand.hk[pick[0]][rs[0]].k, cand.hk[pick[1]][rs[1]].k,
                                   cand.hk[pick[2]][rs[2]].k};
                const double gzs[3] = {gz[pick[0]][rs[0]], gz[pick[1]][rs[1]], gz[pick[2]][rs[2]]};
                for (int w = wLo; w <= wHi; ++w) {
                    if (u == 0 && v == 0 && w == 0) continue;
                    // 4x3 least squares: row 1 is the specular constraint (u,v,w)·(κ,λ,ρ) = ζ,
                    // and the remaining 3 rows are g_z = κh + λk + ρl for the three starting peaks.
                    // The constraint row cannot be eliminated exactly down to 2 unknowns: when two starting peaks share
                    // the same (h,k) (as DIP's (1,1,-1) and (1,1,0) do), the 2x2 normal equations are singular,
                    // whereas the 4x3 system still determines (κ,λ,ρ) because l differs.
                    double A4[4][3], rhs4[4], klr[3];
                    A4[0][0] = u; A4[0][1] = v; A4[0][2] = w; rhs4[0] = zeta;
                    for (int j = 0; j < 3; ++j) {
                        A4[j + 1][0] = hs[j];
                        A4[j + 1][1] = ks[j];
                        A4[j + 1][2] = ls[j];
                        rhs4[j + 1] = gzs[j];
                    }
                    if (!lstsq4x3(A4, rhs4, klr)) continue;
                    emit(klr[0], klr[1], klr[2], w, local);
                }
            }
        }
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(mu);
            raws.insert(raws.end(), local.begin(), local.end());
        }
    }, prm.threads);
    };

    // ---------------- Path C: 2-peak minimal sample + per-peak early exit (deterministic RANSAC) ----------------
    // There are only 3 unknowns (κ,λ,ρ) and the specular row is itself one equation, so **two peaks form a minimal solution**:
    //     [ u   v   w  ] [κ]   [ ζ  ]
    //     [ h1  k1  l1 ] [λ] = [ gz1]
    //     [ h2  k2  l2 ] [ρ]   [ gz2]
    // Number of hypotheses ((2L+1)·R)² = 52² = 2704 (per w, per pair), 52x fewer than 140608 for the 3-peak 4×3 case.
    // Each hypothesis first gets an early-exit test on the "minimum g_z-channel residual"; only very few reach the full criteria chain in emit.
    // The seed pairs are 2 pairs out of "the 3 peaks with the smallest q_z"; because q_xy>0 has already filtered out specular-rod points,
    // a seed cannot be parallel to the specular row (q_xy would be 0 in that case).
    auto runRansac = [&]() {
        std::vector<Raw> local;
        // Per-peak minimum g_z residual: l is back-solved directly from (κ,λ,ρ) and rounded to the nearest integer
        auto gzRes = [&](int i, double kap, double lam, double rho) {
            double best = 1e300;
            const int nr = (int)sub3[i].size();
            for (int j = 0; j < nr; ++j) {
                const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
                const double lr = std::round((gz[i][j] - kap * h - lam * k) / rho);
                if (std::fabs(lr) > (double)L) continue;
                const double r = std::fabs(gz[i][j] - kap * h - lam * k - rho * lr);
                if (r < best) best = r;
            }
            return best;
        };
        // g_z residual -> q_z residual: q_z = (gp·ζ + sub3)/q_s ⇒ δq_z = δgp·ζ/q_s.
        // 3·dqzCutoff is taken as the upper bound, matching the per-peak maxAbs gate in emit (a superset of it, so nothing is wrongly killed).
        const double tolGz = 3.0 * prm.dqzCutoff * qs / std::max(1e-12, zeta);
        // Order of w: 0 first ((u,v,0) planes are the most common), then ±1, ±2; affects only early exits, not the result set
        int ws[32], nws = 0;
        for (int d = 0; d <= W; ++d) {
            if (wLo <= d && d <= wHi) ws[nws++] = d;
            if (d != 0 && -d >= wLo && -d <= wHi) ws[nws++] = -d;
        }
        const int pairSet[2][2] = {{0, 1}, {0, 2}};
        for (int iw = 0; iw < nws; ++iw) {
            const int w = ws[iw];
            if (u == 0 && v == 0 && w == 0) continue;
            // Specular row r1 = (u,v,w); its cross-product term c3 = r1 × r2 does not depend on the third row
            for (int ip = 0; ip < 2; ++ip) {
                const int ia = pick[pairSet[ip][0]], ib = pick[pairSet[ip][1]];
                const int na = (int)sub3[ia].size(), nb = (int)sub3[ib].size();
                for (int ja = 0; ja < na; ++ja) {
                    const int h1 = cand.hk[ia][ja].h, k1 = cand.hk[ia][ja].k;
                    const double g1 = gz[ia][ja];
                    const double s1 = std::fabs((double)u) + std::fabs((double)v) + std::fabs((double)w);
                    const double s2 = std::fabs((double)h1) + std::fabs((double)k1);
                    for (int jb = 0; jb < nb; ++jb) {
                        const int h2 = cand.hk[ib][jb].h, k2 = cand.hk[ib][jb].k;
                        const double g2 = gz[ib][jb];
                        for (int l1 = -L; l1 <= L; ++l1) {
                            // c3 = r1 × r2 (r2 is fixed: (h1,k1,l1))
                            const double c3x = (double)v * (double)l1 - (double)w * (double)k1;
                            const double c3y = (double)w * (double)h1 - (double)u * (double)l1;
                            const double c3z = (double)u * (double)k1 - (double)v * (double)h1;
                            const double s3 = s2 + std::fabs((double)l1);
                            for (int l2 = -L; l2 <= L; ++l2) {
                                // c1 = r2 × r3 , c2 = r3 × r1
                                const double c1x = (double)k1 * (double)l2 - (double)l1 * (double)k2;
                                const double c1y = (double)l1 * (double)h2 - (double)h1 * (double)l2;
                                const double c1z = (double)h1 * (double)k2 - (double)k1 * (double)h2;
                                const double c2x = (double)k2 * (double)w - (double)l2 * (double)v;
                                const double c2y = (double)l2 * (double)u - (double)h2 * (double)w;
                                const double c2z = (double)h2 * (double)v - (double)k2 * (double)u;
                                const double det = (double)u * c1x + (double)v * c1y + (double)w * c1z;
                                const double sc = s1 * s3 * (std::fabs((double)h2) +
                                                             std::fabs((double)k2) +
                                                             std::fabs((double)l2));
                                if (std::fabs(det) < 1e-9 * std::max(1e-12, sc)) continue;
                                const double inv = 1.0 / det;
                                const double kap = (zeta * c1x + g1 * c2x + g2 * c3x) * inv;
                                const double lam = (zeta * c1y + g1 * c2y + g2 * c3y) * inv;
                                const double rho = (zeta * c1z + g1 * c2z + g2 * c3z) * inv;
                                if (!(rho > 1e-9)) continue;
                                // Early exit: check 2 non-seed peaks first; only if they pass is the full chain handed to emit
                                bool ok = true;
                                int nChk = 0;
                                for (int i = 0; i < N && ok && nChk < 2; ++i) {
                                    if (i == ia || i == ib) continue;
                                    if (gzRes(i, kap, lam, rho) > tolGz) ok = false;
                                    ++nChk;
                                }
                                if (!ok) continue;
                                emit(kap, lam, rho, w, local);
                            }
                        }
                    }
                }
            }
        }
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(mu);
            raws.insert(raws.end(), local.begin(), local.end());
        }
    };

    // ---------------- Path D: integer least squares (difference-vector lattice + interval sweep, no bins) ----------------
    //  1) ρ candidates: |Δgz|/n of peak pairs sharing (h,k) (n=1..2L) -- a "1D difference-vector lattice"; no histogram,
    //     instead direct relative-tolerance clustering + sorting by support count (removes bin error and harmonic heuristics).
    //  2) For each (ρ,w): write the feasible set of t as "one interval per peak × per (h,k) × per l"
    //     (|gz-A-tB-ρl| ≤ tol is linear in t), then sweep endpoints to obtain the **t covered by the most peaks**.
    //     This is the **exact global optimum** of t for the given (ρ,w), not a histogram peak.
    //  3) Take the midpoint of each candidate cell -> one round of (κ,λ) least squares -> emit.
    auto runIlls = [&](bool& votable) {
        votable = false;
        if (u == 0 && v == 0) return; // the (t,ρ) parametrization degenerates; leave it to the u=v=0 enumeration
        const double uv2 = (double)u * u + (double)v * v;
        const double tolGz = 3.0 * prm.dqzCutoff * qs / std::max(1e-12, zeta);
        const double rhoMin = 0.05, rhoMax = 2.2;

        std::unordered_map<long long, std::vector<double>> groups;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < (int)sub3[i].size(); ++j)
                groups[(long long)cand.hk[i][j].h * 1000 + cand.hk[i][j].k].push_back(gz[i][j]);
        for (const auto& kv : groups)
            if (kv.second.size() >= 2) {
                votable = true;
                break;
            }
        if (!votable) return;

        // 1) Difference-vector lattice -> ρ candidates (clustered at 0.2% relative tolerance, cluster mean as representative)
        std::vector<double> raw;
        for (const auto& kv : groups) {
            const std::vector<double>& g = kv.second;
            for (size_t a = 0; a < g.size(); ++a)
                for (size_t b = a + 1; b < g.size(); ++b) {
                    const double d = std::fabs(g[a] - g[b]);
                    if (!(d > rhoMin)) continue;
                    const int nmax = std::min(2 * L, (int)std::floor(d / rhoMin));
                    for (int nn = 1; nn <= nmax; ++nn) {
                        const double r = d / (double)nn;
                        if (r < rhoMin || r > rhoMax) continue;
                        raw.push_back(r);
                    }
                }
        }
        if (raw.empty()) return;
        std::sort(raw.begin(), raw.end());
        struct RCl {
            double v = 0.0;
            int n = 0;
        };
        std::vector<RCl> cl;
        for (double r : raw) {
            if (!cl.empty() && std::fabs(r - cl.back().v) < 0.002 * cl.back().v + 2e-5) {
                RCl& c = cl.back();
                c.v += (r - c.v) / (double)(c.n + 1);
                ++c.n;
            } else {
                cl.push_back(RCl{r, 1});
            }
        }
        std::sort(cl.begin(), cl.end(), [](const RCl& x, const RCl& y) {
            if (x.n != y.n) return x.n > y.n;
            return x.v > y.v; // for equal support prefer the larger ρ (harmonic suppression)
        });
        const int maxRho = 24;
        if ((int)cl.size() > maxRho) cl.resize(maxRho);

        int ws[32], nws = 0;
        for (int d = 0; d <= W; ++d) {
            if (wLo <= d && d <= wHi) ws[nws++] = d;
            if (d != 0 && -d >= wLo && -d <= wHi) ws[nws++] = -d;
        }

        std::vector<Raw> local;
        std::vector<std::pair<double, int>> ev; // (endpoint, ±1)
        ev.reserve(2 * (size_t)N * (size_t)R * (size_t)(2 * L + 1));
        std::vector<std::pair<int, double>> cells; // (coverage count, t)
        cells.reserve(2048);

        for (int iw = 0; iw < nws; ++iw) {
            const int w = ws[iw];
            if (u == 0 && v == 0 && w == 0) continue;
            for (const RCl& rc : cl) {
                if (rc.n < 3) continue; // minimum support: 2 votes collide too easily
                const double rho = rc.v;
                const double crhs = zeta - (double)w * rho;
                const double k0 = (double)u * crhs / uv2;
                const double l0 = (double)v * crhs / uv2;
                // 2) Endpoint set
                ev.clear();
                for (int i = 0; i < N; ++i)
                    for (int j = 0; j < (int)sub3[i].size(); ++j) {
                        const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
                        const double Bv = (double)v * (double)h - (double)u * (double)k;
                        if (std::fabs(Bv) < 1e-9) continue;
                        const double base = gz[i][j] - (k0 * (double)h + l0 * (double)k);
                        for (int l = -L; l <= L; ++l) {
                            const double c = base - rho * (double)l;
                            double lo = (c - tolGz) / Bv, hi = (c + tolGz) / Bv;
                            if (lo > hi) std::swap(lo, hi);
                            ev.push_back({lo, +1});
                            ev.push_back({hi, -1});
                        }
                    }
                if (ev.empty()) continue;
                std::sort(ev.begin(), ev.end(),
                          [](const std::pair<double, int>& x, const std::pair<double, int>& y) {
                              if (x.first != y.first) return x.first < y.first;
                              return x.second > y.second; // at the same position +1 before -1 (closed intervals)
                          });
                // 3) Sweep: each "covered cell" records (coverage count, midpoint)
                cells.clear();
                int cur = 0;
                for (size_t e = 0; e < ev.size();) {
                    size_t f = e;
                    while (f + 1 < ev.size() && ev[f + 1].first == ev[e].first) ++f;
                    for (size_t k2 = e; k2 <= f; ++k2) cur += ev[k2].second;
                    const double tA = ev[e].first;
                    const double tB = (f + 1 < ev.size()) ? ev[f + 1].first : tA;
                    if (cur >= 2 && tB > tA) cells.push_back({cur, 0.5 * (tA + tB)});
                    e = f + 1;
                }
                if (cells.empty()) continue;
                const int kWant = 6;
                if ((int)cells.size() > kWant)
                    std::partial_sort(cells.begin(), cells.begin() + kWant, cells.end(),
                                      [](const std::pair<int, double>& x,
                                         const std::pair<int, double>& y) {
                                          if (x.first != y.first) return x.first > y.first;
                                          return x.second < y.second;
                                      });
                const int nc = std::min(kWant, (int)cells.size());
                for (int ci = 0; ci < nc; ++ci) {
                    const double t = cells[(size_t)ci].second;
                    double kap = k0 + (double)v * t;
                    double lam = l0 - (double)u * t;
                    // One round of (κ,λ) least squares (ρ, w, and the specular row fixed)
                    {
                        double Saa = (double)u * u, Sab = (double)u * v, Sbb = (double)v * v;
                        double Sa = (double)u * crhs, Sb = (double)v * crhs;
                        for (int i = 0; i < N; ++i) {
                            double best = 1e300, bh = 0, bk = 0, bq = 0;
                            for (int j = 0; j < (int)sub3[i].size(); ++j) {
                                const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
                                const double lr = std::round((gz[i][j] - kap * h - lam * k) / rho);
                                if (std::fabs(lr) > (double)L) continue;
                                const double r = std::fabs(gz[i][j] - kap * h - lam * k - rho * lr);
                                if (r < best) {
                                    best = r;
                                    bh = h;
                                    bk = k;
                                    bq = gz[i][j] - rho * lr;
                                }
                            }
                            if (best > 1e299) continue;
                            Saa += bh * bh;
                            Sab += bh * bk;
                            Sbb += bk * bk;
                            Sa += bh * bq;
                            Sb += bk * bq;
                        }
                        const double dd = Saa * Sbb - Sab * Sab;
                        if (std::fabs(dd) > 1e-12 * std::max(1.0, Saa * Sbb)) {
                            kap = (Sbb * Sa - Sab * Sb) / dd;
                            lam = (Saa * Sb - Sab * Sa) / dd;
                        }
                    }
                    emit(kap, lam, rho, w, local);
                }
            }
        }
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(mu);
            raws.insert(raws.end(), local.begin(), local.end());
        }
    };

    // ---------------- Path B: PeriodVote/DDL (difference vectors + periodic voting) ----------------
    // Two peaks sharing (h,k) satisfy g_z difference = ρ·Δl, independent of (κ,λ). Folding these difference vectors onto the
    // ρ axis by Δl and voting (equivalent to finding the dominant period of the difference-vector set via a periodogram) yields ρ candidates;
    // then along the specular constraint uκ+vλ = ζ-wρ, (κ,λ) collapse to a single parameter t, and a second-level vote is run over t.
    // Its cost is independent of the NL³×R³×(2W+1) brute force and approximately linear in the number of peaks.
    auto runVote = [&](bool& votable) {
        votable = false;
        const double uv2 = (double)u * u + (double)v * v;
        if (!(uv2 > 1e-12)) return; // with (u,v)=(0,0) the constraint row degenerates; use the enumeration path instead

        const double rhoMin = 0.05, rhoMax = 2.2;
        const int nBins = 4096;
        const double binW = (rhoMax - rhoMin) / (double)nBins;
        std::vector<double> hist((size_t)nBins, 0.0);

        // 1) Group by (h,k); pairwise differences within a group give ρ = Δg_z / Δl
        std::unordered_map<long long, std::vector<double>> groups;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < (int)sub3[i].size(); ++j) {
                const long long key = (long long)cand.hk[i][j].h * 1000 + cand.hk[i][j].k;
                groups[key].push_back(gz[i][j]);
            }
        // With no "peak pair sharing (h,k)" the difference-vector method carries no information; let the caller fall back to enumeration
        for (const auto& kv : groups)
            if (kv.second.size() >= 2) { votable = true; break; }
        if (!votable) return;
        for (const auto& kv : groups) {
            const std::vector<double>& g = kv.second;
            for (size_t a = 0; a < g.size(); ++a)
                for (size_t b2 = a + 1; b2 < g.size(); ++b2) {
                    const double d = std::fabs(g[a] - g[b2]);
                    if (!(d > rhoMin)) continue;
                    const int nmax = std::min(2 * L, (int)std::floor(d / rhoMin));
                    for (int nn = 1; nn <= nmax; ++nn) {
                        const double r = d / (double)nn;
                        if (r < rhoMin || r >= rhoMax) continue;
                        hist[(size_t)((r - rhoMin) / binW)] += 1.0;
                    }
                }
        }
        if (groups.empty()) return;

        // 2) Take the ρ values with the highest vote counts (merging adjacent bins)
        std::vector<std::pair<double, double>> pk; // (vote count, ρ)
        for (int i = 0; i < nBins; ++i) {
            if (!(hist[(size_t)i] > 0.0)) continue;
            double v = hist[(size_t)i];
            if (i > 0) v += 0.5 * hist[(size_t)(i - 1)];
            if (i + 1 < nBins) v += 0.5 * hist[(size_t)(i + 1)];
            pk.push_back({v, rhoMin + ((double)i + 0.5) * binW});
        }
        if (pk.empty()) return;
        // Subharmonic suppression: ρ/2, ρ/3… of the true period ρ also explain the same difference vectors as integer multiples,
        // so their vote counts often tie or even win. Hence within a tie prefer the **larger** ρ (the fundamental),
        // and apply non-maximum suppression over neighbouring bins.
        std::sort(pk.begin(), pk.end(),
                  [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                      if (std::fabs(a.first - b.first) > 0.9) return a.first > b.first;
                      return a.second > b.second;
                  });
        std::vector<double> rhos;
        for (const auto& p : pk) {
            bool dup = false;
            for (double r : rhos)
                if (std::fabs(r - p.second) < 0.02 * p.second + 4.0 * binW) { dup = true; break; }
            if (!dup) rhos.push_back(p.second);
            if ((int)rhos.size() >= 32) break;
        }

        // 3) For each (ρ,w), fix t by a second-level vote -> (κ,λ)
        const double tMin = -6.0, tMax = 6.0;
        const int nTB = 2400;
        const double tW = (tMax - tMin) / (double)nTB;
        std::vector<double> hT((size_t)nTB, 0.0);
        std::vector<Raw> local;

        for (double rho : rhos) {
            for (int w = wLo; w <= wHi; ++w) {
                const double crhs = zeta - (double)w * rho; // uκ+vλ
                const double kap0 = (double)u * crhs / uv2;
                const double lam0 = (double)v * crhs / uv2;
                std::fill(hT.begin(), hT.end(), 0.0);
                bool any = false;
                for (int i = 0; i < N; ++i)
                    for (int j = 0; j < (int)sub3[i].size(); ++j) {
                        const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
                        const double Bv = (double)v * h - (double)u * k;
                        if (std::fabs(Bv) < 1e-9) continue;
                        const double A = kap0 * h + lam0 * k;
                        for (int l = -L; l <= L; ++l) {
                            const double t = (gz[i][j] - A - rho * (double)l) / Bv;
                            if (t < tMin || t >= tMax) continue;
                            hT[(size_t)((t - tMin) / tW)] += 1.0;
                            any = true;
                        }
                    }
                if (!any) continue;
                std::vector<std::pair<double, double>> tp;
                for (int i = 0; i < nTB; ++i) {
                    if (!(hT[(size_t)i] > 0.0)) continue;
                    double v2 = hT[(size_t)i];
                    if (i > 0) v2 += 0.5 * hT[(size_t)(i - 1)];
                    if (i + 1 < nTB) v2 += 0.5 * hT[(size_t)(i + 1)];
                    tp.push_back({v2, tMin + ((double)i + 0.5) * tW});
                }
                std::sort(tp.begin(), tp.end(),
                          [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                              return a.first > b.first;
                          });
                int taken = 0;
                for (const auto& p : tp) {
                    double kap = kap0 + (double)v * p.second;
                    double lam = lam0 - (double)u * p.second;
                    // Refine (κ,λ) for two more rounds using all peaks (l taken as the nearest integer)
                    for (int it = 0; it < 2; ++it) {
                        double Saa = (double)u * u, Sab = (double)u * v, Sbb = (double)v * v;
                        double Sa = (double)u * crhs, Sb = (double)v * crhs;
                        for (int i = 0; i < N; ++i) {
                            double best = 1e300, bh = 0, bk = 0, bq = 0;
                            for (int j = 0; j < (int)sub3[i].size(); ++j) {
                                const int h = cand.hk[i][j].h, k = cand.hk[i][j].k;
                                const double lc = (gz[i][j] - kap * h - lam * k) / rho;
                                const double lr = std::round(lc);
                                if (std::fabs(lr) > (double)L) continue;
                                const double res =
                                    std::fabs(gz[i][j] - kap * h - lam * k - rho * lr);
                                if (res < best) {
                                    best = res;
                                    bh = h;
                                    bk = k;
                                    bq = gz[i][j] - rho * lr;
                                }
                            }
                            if (best > 1e299) continue;
                            Saa += bh * bh;
                            Sab += bh * bk;
                            Sbb += bk * bk;
                            Sa += bh * bq;
                            Sb += bk * bq;
                        }
                        const double dd = Saa * Sbb - Sab * Sab;
                        if (!(std::fabs(dd) > 1e-12 * std::max(1.0, Saa * Sbb))) break;
                        kap = (Sbb * Sa - Sab * Sb) / dd;
                        lam = (Saa * Sb - Sab * Sa) / dd;
                    }
                    emit(kap, lam, rho, w, local);
                    if (++taken >= 10) break;
                }
            }
        }
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(mu);
            raws.insert(raws.end(), local.begin(), local.end());
        }
    };

    if (prm.part1Alg == Part1Algorithm::PeriodVoteDdl) {
        bool votable = false;
        runVote(votable);
        // Fall back within this candidate only when voting is "structurally impossible" (no peak pair sharing (h,k));
        // if votes exist but yield no solution, the (a,b,gamma) candidate itself is invalid and is skipped.
        if (!votable) runEnum();
    } else if (prm.part1Alg == Part1Algorithm::Ransac) {
        runRansac();
    } else if (prm.part1Alg == Part1Algorithm::Ills) {
        bool votable = false;
        runIlls(votable);
        // With no same-family peak pairs / with (u,v)=(0,0) the (t,ρ) parametrization degenerates; fall back to enumeration
        if (!votable) runEnum();
    } else {
        runEnum();
    }

    if (raws.empty()) return {};
    std::sort(raws.begin(), raws.end(), [](const Raw& a, const Raw& b) { return a.rms < b.rms; });
    // Refinement budget: at most refineTop raw solutions per candidate (capped at 200),
    // and solutions in the same basin of attraction (identical after 0.25 Å / 2° quantization) are refined only once.
    // One refinement = 2 rounds of LM(60) + 2 passes of 13³ grid assignment, the main cost of Part TWO.
    // Note: too small a cap lets "multiple starting points in the same basin" exhaust the budget and the solution count collapses
    // (on DIP, cap=32 yields only 1 solution while cap=200 yields 12).
    const int kRefineLimit = std::max(1, std::min(prm.refineTop, 200));
    {
        std::vector<Raw> picked;
        picked.reserve((size_t)std::min((int)raws.size(), kRefineLimit));
        std::unordered_set<uint64_t> qseen;
        for (const Raw& r : raws) {
            if (!qseen.insert(rawKey(r.cell, r.w)).second) continue;
            picked.push_back(r);
            if ((int)picked.size() >= kRefineLimit) break;
        }
        raws.swap(picked);
    }

    // Refine + reduce (A10: Top-K only). Under a crystal-system constraint only that system's free parameters are fitted
    // (Cubic: a only; Tetragonal/Hexagonal: a,c; Orthorhombic: a,b,c …).
    const std::array<bool, 6> freeMask = systemFreeMask(prm.system);
    std::vector<Solution> sols;
    sols.reserve(raws.size());
    std::unordered_set<uint64_t> seen;
    // The final assignment uses the first-step (h,k) range (default hkLimit2=6, matching the previous hard-coded value)
    const int L2 = std::max(1, prm.hkLimit2);
    const std::vector<Hkl> cands = hklGrid(L2, L2, L2);
    for (const Raw& r : raws) {
        Cell cell = r.cell;
        Hkl uvw0{u, v, r.w};
        std::vector<Hkl> ah = r.hkl;
        Orient orient = Orient::fromUvw(uvw0);
        // Upper bound on |G| for assignment candidates: if the current assignment's maximum deviation is dmax, indices with |G| > qmax + 3·dmax
        // cannot be closer than the current assignment (the cell barely moves during refinement, so a 3x margin suffices), which prunes the
        // candidate table of 13³=2197 hkl down to the tens or hundreds actually relevant.
        double maxG = 1e30;
        if (constrained) {
            // Compare axis settings before projection: 3 unique axes for tetragonal/hexagonal/trigonal, 6 axis permutations for orthorhombic,
            // 3 unique axes for monoclinic (hexagonal also includes the equivalent b->b-a setting with γ=60°). Take the setting with the lowest
            // cost (joint assignment residual after projection); only the winner enters subsequent refinement, otherwise
            // the correct solution would be destroyed by the projection as soon as it falls on a non-default setting.
            // Sort by projection cost first (pure metric computation, very cheap); usually the first one is the setting of the correct solution,
            // and full assignment is run for only a few settings, avoiding the assignment cost of 6 settings × per candidate.
            // Two-stage: first score and sort every setting with evaluate (O(number of peaks), no hkl grid search),
            // then run full assignHkl only for the best setting (grid assignment is about an order of magnitude more expensive than evaluate);
            // if its residual is already within tolerance, stop; otherwise try the next settings in score order.
            std::vector<std::pair<double, std::array<int, 9>>> settings;
            for (const std::array<int, 9>& T0 : axisSettings(prm.system, r.cell.ga))
                settings.push_back({axisSettingCost(prm.system, r.cell, T0), T0});
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
                if (sp.first >= 1e30) break; // settings that failed the transform are all sorted last
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
                    bestHkl = std::move(hk); // move only after a better result is found, so saved results are not lost
                    bestMaxG = mg;
                }
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
        // Reassign all peaks ((h,k,l) may change after refinement) and refine one more round
        // After projection, recompute the |G| upper bound from the cell and the assigned indices (no longer the pre-projection estimate / full box)
        if (constrained) maxG = assignMaxG(cell, orient, peaks, ah);
        assignHkl(cell, orient, peaks, cands, ah, nullptr, 1e30, 1.0, maxG);
        refineCell(peaks, ah, orient, cell, freeMask, 60, prm.system);
        if (constrained) applySystemConstraints(prm.system, cell);
        if (constrained) maxG = assignMaxG(cell, orient, peaks, ah);
        assignHkl(cell, orient, peaks, cands, ah, nullptr, 1e30, 1.0, maxG);
        // After projection, index aliasing can lock LM in a local minimum: run one progressive reassignment rescue inside the constrained subspace.
        // Triggered only when rmsQxyz after refinement is still > 0.3× the Δq tolerance (nothing to recover once within tolerance).
        if (constrained) {
            maxG = assignMaxG(cell, orient, peaks, ah);
            const double rmsGate = 0.3 * prm.dqxyzCutoff;
            const double rmsSkip = 2.0 * prm.dqxyzCutoff;
            refineConstrained(peaks, orient, cands, prm.system, cell, ah, maxG, rmsGate,
                              rmsSkip);
            // Monoclinic: canonical β branch (for β>90 use the equivalent c->-c setting, a purely descriptive transform);
            // when all l are 0 and the specular direction has w=0, c and β are inseparable, so canonicalize to β=90 to avoid a family of β.
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

        // Fixed contact-plane mode: the user-specified (u,v,w) must be preserved exactly, so no whole-cell reduction is performed.
        // Reduction would transform directions such as (0 0 1) into the reduced basis (e.g. into (-1 0 -1)),
        // causing out-of-plane crystal planes to leak into the "fixed plane" search results.
        // Crystal-system constrained mode likewise skips reduction: reduceCell changes the basis vectors and would rewrite a γ=120° hexagonal cell
        // or an α=γ=90° monoclinic cell into an equivalent form that does not follow the standard axis setting of the crystal system.
        Cell rc = cell;
        Hkl ruvw = uvw0;
        std::vector<Hkl> rhkl = ah;
        if (!prm.fixedUvw && !constrained) {
            const Reduction red = reduceCell(cell);
            if (!red.ok) continue;
            rc = red.cell;
            ruvw = applyT(Hkl{u, v, r.w}, red.T);
            rhkl.clear();
            rhkl.reserve(ah.size());
            for (const Hkl& x : ah) rhkl.push_back(applyT(x, red.T));
        }

        const Quality q = evaluate(rc, Orient::fromUvw(ruvw), peaks, rhkl);
        if (!(q.rmsQxyz <= prm.dqxyzCutoff)) continue;
        // Supercell check: the observed reciprocal points must generate that reciprocal lattice (index == 1), otherwise it is a d-fold supercell
        if (!acceptLatticeIndex(prm.system, rhkl)) continue;
        if (q.rmsQxy > prm.dqxyCutoff * 3.0 || q.rmsQz > prm.dqzCutoff * 3.0) continue;
        // Specular-peak criterion: the refined cell must still reproduce the measured q_spec, otherwise it is overfitting in a purely
        // geometric sense (a large cell + high-index (hkl) can fit any peak position)
        {
            const double qsp = qSpecOf(reciprocalMatrix(rc), ruvw);
            if (std::fabs(qsp - qs) > prm.dqspec) continue;
        }
        const double vol = cellVolume(rc);
        if (prm.restrictCell) {
            if (rc.a < prm.aLo || rc.a > prm.aHi || rc.b < prm.bLo || rc.b > prm.bHi ||
                rc.c < prm.cLo || rc.c > prm.cHi)
                continue;
            if (vol < prm.volLo || vol > prm.volHi) continue;
        }
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
        s.hkl = rhkl;
        // Sign-representative normalization: (u,v,w) and (−u,−v,−w) are the same set of crystal planes; unify to the branch with the larger u+v+w
        // (depending on the reduction T, the same physical solution may fall into either branch). fixedUvw keeps the user-specified direction.
        if (!prm.fixedUvw) canonicalizeUvw(s);
        // Constrained mode deduplicates by a crystal-system-invariant cell fingerprint + plane normal (basis-change invariants): equivalent solutions
        // from axis permutations/settings (e.g. the 6 axis settings of the same orthorhombic lattice) merge; the triclinic path keeps the original key.
        uint64_t key = constrained ? constrainedCellKey(rc, s.uvw) : cellKey(rc);
        if (!constrained) {
            key ^= (uint64_t)(s.uvw.h * 10000 + s.uvw.k * 100 + s.uvw.l) * 1099511628211ull;
        }
        if (!seen.insert(key).second) continue;
        sols.push_back(std::move(s));
    }
    std::sort(sols.begin(), sols.end(),
              [](const Solution& a, const Solution& b) { return a.rmsQxyz < b.rmsQxyz; });
    if ((int)sols.size() > prm.maxSolutions) sols.resize(prm.maxSolutions);
    return sols;
}

// ---------------------------------------------------------------------------
//  Driver
// ---------------------------------------------------------------------------
std::vector<Solution> runSpecIndex(const std::vector<Peak>& peaksIn, const SpecParams& prm,
                                   const Progress& prog) {
    // Computation input is uniformly in |q_xy| (the sign only reflects the mirror handedness of the detection geometry).
    std::vector<Peak> peaks = peaksIn;
    normalizeQxySigns(peaks);
    std::vector<AbgCand> abg = improved::runPartOne(peaks, prm, prog);
    if (prog.stop()) return {};
    if (abg.empty()) return {};

    // A2: candidate-level parallelism. Previously parallel inside (chunked by l-triple); tasks were too fine-grained and uneven in cost,
    // giving only ~3.7x on 8 threads; now distributed over candidates at the outer level (dynamic queue) with Part TWO serial inside.
    // D1: process candidates in ascending area order in batches, with adaptive early stopping at the end of each batch.
    const int np = (int)abg.size();
    int nt = (prm.threads > 0) ? prm.threads : (int)hardwareThreads();
    nt = std::max(1, std::min(nt, np));
    const int batch = std::max(1, nt);

    // Early-stop criteria: a "very good" solution has appeared (rmsQxyz <= 20% of the threshold) and enough candidates have been processed,
    // or a large batch of consecutive candidates failed to improve the best solution. Candidates are in ascending in-plane area, later ones overfit more easily.
    const double kGoodRms = 0.2 * prm.dqxyzCutoff;
    // The cost of one candidate decides "how many more candidates to scan after a good solution is found": PeriodVote/DDL takes only a few ms per candidate,
    // so scanning more is needed to list the solutions converged from different (a,b,gamma); enumeration takes ~40 ms per candidate,
    // so the lower bound is used to avoid ending too early (which would yield only one solution).
    const bool cheapCand = (prm.part1Alg != Part1Algorithm::Improved &&
                            prm.part1Alg != Part1Algorithm::Original);
    const int kGoodMinCands = cheapCand ? 160 : 64;
    const int kStallCands = cheapCand ? 160 : 64;

    std::vector<Solution> all;
    double bestRms = 1e30;
    int sinceImprove = 0;
    for (int base = 0; base < np; base += batch) {
        if (prog.stop()) break;
        const int end = std::min(np, base + batch);
        {
            char buf[176];
            std::snprintf(buf, sizeof(buf),
                          "Part 2: candidate %d/%d   a=%.3f b=%.3f gamma=%.2f  rms=%.4f", base + 1, np,
                          abg[base].a, abg[base].b, abg[base].ga, abg[base].rmsQxy);
            prog((double)base / (double)np, buf);
        }
        std::vector<std::vector<Solution>> part((size_t)(end - base));
        parallelForEach(end - base, [&](int k) {
            SpecParams p1 = prm;
            p1.threads = 1; // serial inner loop to avoid nested-parallel oversubscription
            part[(size_t)k] = improved::runPartTwo(peaks, abg[base + k], p1);
        }, nt);

        bool advanced = false;
        for (auto& v : part)
            for (auto& s : v) {
                if (s.rmsQxyz < bestRms) {
                    bestRms = s.rmsQxyz;
                    advanced = true;
                }
                all.push_back(std::move(s));
            }
        sinceImprove = advanced ? 0 : sinceImprove + (end - base);

        if (prm.earlyStop) {
            if (bestRms <= kGoodRms && end >= std::min(np, kGoodMinCands)) break;
            if (!all.empty() && sinceImprove >= kStallCands) break;
        }
    }
    std::sort(all.begin(), all.end(),
              [](const Solution& a, const Solution& b) { return a.rmsQxyz < b.rmsQxyz; });
    // Global deduplication (constrained mode uses the crystal-system invariant key, merging equivalent solutions from axis settings/permutations)
    std::vector<Solution> out;
    std::unordered_set<uint64_t> seen;
    const bool constrained = (prm.system != CrystalSystem::Triclinic);
    for (Solution& s : all) {
        uint64_t key =
            constrained ? constrainedCellKey(s.cell, s.uvw) : cellKey(s.cell);
        if (!constrained) {
            key ^= (uint64_t)(s.uvw.h * 10000 + s.uvw.k * 100 + s.uvw.l) * 1099511628211ull;
        }
        if (!seen.insert(key).second) continue;
        out.push_back(std::move(s));
        if ((int)out.size() >= prm.maxSolutions) break;
    }
    // Fallback (happens at most once):
    //   - cheap per-candidate implementations (PeriodVote/DDL, Ransac, Ills) find no solution in the whole round -> rerun with l-triple enumeration;
    //   - Part ONE used the "first-seed nearest-neighbour box" and the whole round finds nothing -> rerun with the full box.
    const bool cheapAlg = (prm.part1Alg != Part1Algorithm::Improved &&
                           prm.part1Alg != Part1Algorithm::Original);
    if (out.empty() && (cheapAlg || prm.seed3Half > 0)) {
        SpecParams p2 = prm;
        if (cheapAlg) p2.part1Alg = Part1Algorithm::Improved;
        p2.seed3Half = 0; // the fallback always uses the full box, so the recursion necessarily terminates
        return improved::runSpecIndex(peaks, p2, prog);
    }
    return out;
}

} // namespace improved
} // namespace gid
