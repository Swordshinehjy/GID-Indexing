// ============================================================================
//  fixedplane.cpp — see fixedplane.h
//
//  Three paths; returns as soon as any one succeeds (fastest first):
//    0) Initial-cell fast path: refineAssign (alternating LM iterations) — the fast lane when an approximate cell is already available
//    1) Two-step method when specular rod peaks exist (fixedUvw + rod fundamental p)
//    2) FP-DL deterministic search (Fixed-Plane Direct Linear)
//
//  Mathematical basis of FP-DL (G = G*, 6 components [G11,G22,G33,G12,G13,G23]; t = primitive (u,v,w)):
//    (A)  q_i²  = x_i G x_iᵀ                             linear in G
//    (B)  q_z,i = ν·x_i,   ν = G tᵀ/√(t G tᵀ)             linear in ν
//    (C)  G tᵀ  = (t·ν)·ν                                 derived from q_s = t·ν, linear in G
//  Derivation of (C): q_s = √(t G tᵀ) = √(t · G tᵀ) = √(q_s (t·ν)) ⇒ q_s = t·ν ⇒ G tᵀ = (t·ν)ν.
//
//  The whole search is therefore purely algebraic and needs no LM in the loop:
//    S1 three peaks + integer Cramer -> ν (exact 3×3 solve, det≠0 as the exact rank filter)
//       + a "crystal-system compatibility" hard screen on the (C) rows only (cubic/tetragonal/hexagonal pin the ν direction here)
//       + sort score = (C) residual + q_z residual, each first divided by its own tolerance (dimensionless)
//    S2 3 rows of (C) + m rows of (A) -> G (linear LS; triclinic is full rank at m ≥ 3)
//    S3 3 rounds of alternating "nearest-neighbour assignment + linear re-solve"
//    S4 only surviving candidates (usually single digits) go through refineCell/refineAssign polish
//
//  Why random RANSAC was removed: the assignment search space is |cands|^nSeed (728⁵ ≈ 2e14 by default),
//  so 2000 random trials always had a hit probability of 0; it had in effect degenerated into "local optimization with random restarts",
//  while each trial still paid for 60 LM rounds (numerical Jacobian).
//
//  Contact-plane semantics: the model depends only on the direction of (u,v,w) (ν is invariant under n → m·n); the entry normalizes to a primitive direction;
//  rod-peak q_z is an integer multiple of the fundamental, and q_spec is used only for posterior validation, never as a solving constraint.
// ============================================================================
#include "fixedplane.h"
#include "engine.h"
#include "parallel.h"
#include "specindex.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <numeric>
#include <set>
#include <vector>

namespace gid {
namespace {

// ============================================================================
//  Numerics shared with nospec.cpp (deliberately copied: nospec.cpp is bound by the
//  cross-repository sync rules of tools/verify_core_sync.ps1, this file is not.
//  When changing anything here, keep behaviour consistent with the corresponding functions in nospec.cpp.)
// ============================================================================

// 3x3 integer determinant (exact prefilter)
inline int intDet3(const int a[3], const int b[3], const int c[3]) {
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
           a[2] * (b[0] * c[1] - b[1] * c[0]);
}

// Cramer numerators for the minimal 3-peak sample of the (B) equations: the rows of adj(M) are the cross products of the other two rows, all integers.
// Returns det (also serving as the rank prefilter), w = adj(M)·rhs; the true solution is ν = w/det.
//
// 4.3M calls in S1, so **no division is done here**: the later screens (fundamental sign, (C) crystal-system compatibility)
// are all independent of the scale of ν and can be judged directly on w; only the few candidates that pass need w/det.
// When det < 0 the caller negates w (ν and −ν are equivalent; normalizing to the det > 0 branch is enough).
inline int cramerNumerators(const int a[3], const int b[3], const int c[3], const double rhs[3],
                            double w[3]) {
    const int wx[3] = {b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2],
                       b[0] * c[1] - b[1] * c[0]};
    const int wy[3] = {c[1] * a[2] - c[2] * a[1], c[2] * a[0] - c[0] * a[2],
                       c[0] * a[1] - c[1] * a[0]};
    const int wz[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                       a[0] * b[1] - a[1] * b[0]};
    const int det = a[0] * wx[0] + a[1] * wx[1] + a[2] * wx[2];
    w[0] = (double)(wx[0]) * rhs[0] + (double)(wx[1]) * rhs[1] + (double)(wx[2]) * rhs[2];
    w[1] = (double)(wy[0]) * rhs[0] + (double)(wy[1]) * rhs[1] + (double)(wy[2]) * rhs[2];
    w[2] = (double)(wz[0]) * rhs[0] + (double)(wz[1]) * rhs[1] + (double)(wz[2]) * rhs[2];
    return det;
}

// Number of free parameters of the (A) equations (matches the x length of fillMetricArr exactly)
constexpr int metricFreeCount(CrystalSystem sys) {
    switch (sys) {
        case CrystalSystem::Cubic: return 1;
        case CrystalSystem::Tetragonal:
        case CrystalSystem::Hexagonal:
        case CrystalSystem::Trigonal: return 2;
        case CrystalSystem::Orthorhombic: return 3;
        case CrystalSystem::Monoclinic: return 4;
        default: return 6;
    }
}

// Fill the full G* from the free components (constrained components completed per crystal system)
std::array<double, 6> fillMetricArr(CrystalSystem sys, const double* x) {
    std::array<double, 6> G{0, 0, 0, 0, 0, 0};
    switch (sys) {
        case CrystalSystem::Cubic:
            G[0] = x[0]; G[1] = G[0]; G[2] = G[0];
            break;
        case CrystalSystem::Tetragonal:
            G[0] = x[0]; G[2] = x[1]; G[1] = G[0];
            break;
        case CrystalSystem::Hexagonal:
        case CrystalSystem::Trigonal:
            G[0] = x[0]; G[2] = x[1]; G[1] = G[0]; G[3] = 0.5 * G[0];
            break;
        case CrystalSystem::Orthorhombic:
            G[0] = x[0]; G[1] = x[1]; G[2] = x[2];
            break;
        case CrystalSystem::Monoclinic:
            G[0] = x[0]; G[1] = x[1]; G[2] = x[2]; G[4] = x[3];
            break;
        default:
            for (size_t i = 0; i < 6; ++i) G[i] = x[i];
            break;
    }
    return G;
}

// Design row of the (A) equations for the free metric components.
// Constrained components (G22=G11, G33=G11, G12=G11/2 …) must be **folded into** the free parameters,
// never "pick only the free columns" — otherwise q² = G11·h² is completely wrong for cubic.
inline void metricRow(CrystalSystem sys, double h, double k, double l, double* a) {
    switch (sys) {
        case CrystalSystem::Cubic:
            a[0] = h * h + k * k + l * l;
            return;
        case CrystalSystem::Tetragonal:
            a[0] = h * h + k * k;
            a[1] = l * l;
            return;
        case CrystalSystem::Hexagonal:
        case CrystalSystem::Trigonal:
            a[0] = h * h + k * k + h * k;
            a[1] = l * l;
            return;
        case CrystalSystem::Orthorhombic:
            a[0] = h * h; a[1] = k * k; a[2] = l * l;
            return;
        case CrystalSystem::Monoclinic:
            a[0] = h * h; a[1] = k * k; a[2] = l * l; a[3] = 2.0 * h * l;
            return;
        default:
            a[0] = h * h; a[1] = k * k; a[2] = l * l;
            a[3] = 2.0 * h * k; a[4] = 2.0 * h * l; a[5] = 2.0 * k * l;
            return;
    }
}

// Solve min||A x - b|| (normal equations + Gauss elimination), A: m*n row-major
bool solveLS(const std::vector<double>& A, const std::vector<double>& b, int m, int n,
             std::vector<double>& x) {
    std::vector<double> M((size_t)n * n, 0.0), v((size_t)n, 0.0);
    for (int i = 0; i < m; ++i)
        for (int p = 0; p < n; ++p) {
            const double ap = A[(size_t)i * n + p];
            v[(size_t)p] += ap * b[(size_t)i];
            for (int q = p; q < n; ++q) M[(size_t)p * n + q] += ap * A[(size_t)i * n + q];
        }
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < p; ++q) M[(size_t)p * n + q] = M[(size_t)q * n + p];
    double dmax = 0.0;
    for (int p = 0; p < n; ++p) dmax = std::max(dmax, std::fabs(M[(size_t)p * n + p]));
    if (!(dmax > 0.0)) return false;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(M[(size_t)i * n + k]) > std::fabs(M[(size_t)p * n + k])) p = i;
        if (std::fabs(M[(size_t)p * n + k]) < 1e-11 * dmax) return false;  // rank deficient
        if (p != k) {
            for (int j = 0; j < n; ++j) std::swap(M[(size_t)p * n + j], M[(size_t)k * n + j]);
            std::swap(v[(size_t)p], v[(size_t)k]);
        }
        for (int i = k + 1; i < n; ++i) {
            const double f = M[(size_t)i * n + k] / M[(size_t)k * n + k];
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) M[(size_t)i * n + j] -= f * M[(size_t)k * n + j];
            v[(size_t)i] -= f * v[(size_t)k];
        }
    }
    x.assign((size_t)n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double s = v[(size_t)i];
        for (int j = i + 1; j < n; ++j) s -= M[(size_t)i * n + j] * x[(size_t)j];
        x[(size_t)i] = s / M[(size_t)i * n + i];
    }
    for (double d : x)
        if (!std::isfinite(d)) return false;
    return true;
}

bool positiveDefinite(const std::array<double, 6>& g) {
    if (!(g[0] > 0) || !(g[1] > 0) || !(g[2] > 0)) return false;
    const double m2 = g[0] * g[1] - g[3] * g[3];
    if (!(m2 > 0)) return false;
    const double d = g[0] * (g[1] * g[2] - g[5] * g[5]) - g[3] * (g[3] * g[2] - g[5] * g[4]) +
                     g[4] * (g[3] * g[5] - g[1] * g[4]);
    return d > 0;
}

// ============================================================================
//  FP-DL specific
// ============================================================================

inline int gcd3(int a, int b, int c) {
    a = std::abs(a);
    b = std::abs(b);
    c = std::abs(c);
    return std::gcd(std::gcd(a, b), c);
}

inline double q2Of(const Peak& p) { return p.qxy * p.qxy + p.qz * p.qz; }

// Compatibility of the fundamental p with the measured rod fundamental pMeas: p must be pMeas / n (n = observed order).
// Under single-order ambiguity (rod.ambiguous) the caller passes pMeas/n instead of pMeas.
inline bool rodConsistent(double p, double pMeas, double tol, int nmax) {
    if (!(pMeas > 0.0)) return true;
    for (int n = 1; n <= std::max(1, nmax); ++n)
        if (std::fabs(p - pMeas / (double)n) <= tol) return true;
    return false;
}

// Given ν, enumerate indices inside the index box that are compatible with q_z: use the largest |ν| component as the solving axis
// (one division, and the least sensitive), enumerate the other two axes. Complexity drops from |box| to the product of the other two ranges.
// Returns whether one was found; bestD outputs the minimum |Δq_z|, bestX outputs that index.
bool bestQzCand(const int lim[3], const std::array<double, 3>& nu, double qz, double tol,
                double& bestD, Hkl* bestX) {
    int ax = 0;
    double amax = std::fabs(nu[0]);
    for (int t = 1; t < 3; ++t)
        if (std::fabs(nu[t]) > amax) { amax = std::fabs(nu[t]); ax = t; }
    if (!(amax > 1e-9)) return false;
    const int a1 = (ax == 0) ? 1 : 0;
    const int a2 = (ax == 2) ? 1 : 2;
    bool found = false;
    bestD = 1e30;
    const double inv = 1.0 / nu[ax];
    for (int i = -lim[a1]; i <= lim[a1]; ++i)
        for (int j = -lim[a2]; j <= lim[a2]; ++j) {
            int x[3] = {0, 0, 0};
            x[a1] = i;
            x[a2] = j;
            const double s = nu[0] * x[0] + nu[1] * x[1] + nu[2] * x[2];
            const double v = (qz - s) * inv;
            const int iv = (int)std::lround(v);
            if (iv < -lim[ax] || iv > lim[ax]) continue;
            const double d = std::fabs((iv - v) * nu[ax]);
            if (d > tol) continue;
            x[ax] = iv;
            if (x[0] == 0 && x[1] == 0 && x[2] == 0) continue;
            if (d < bestD) {
                bestD = d;
                found = true;
                if (bestX) *bestX = Hkl{x[0], x[1], x[2]};
            }
        }
    return found;
}

// Enumerate **all** indices compatible with q_z (not just the nearest one): q_z constrains only the normal component,
// the in-plane (h,k) are completely free, and taking only "the nearest one" would exclude the true solution.
//
// Sort key = (|x|², |Δq_z|): q_z acts only as a **gate** (tol is 3× the acceptance tolerance, very loose); sorting looks at index length.
// The reverse (comparing Δq_z first) pushes the true indices out of the top ranks: a long vector such as (4,-3,2) can beat the
// ~1e-3 measurement noise of the true index with an exact 1e-4 hit — the total failure of cubic-012 was caused by exactly this.
// Returns the count (<= maxOut).
int qzCands(const int lim[3], const std::array<double, 3>& nu, double qz, double tol,
            Hkl* out, int maxOut) {
    int ax = 0;
    double amax = std::fabs(nu[0]);
    for (int t = 1; t < 3; ++t)
        if (std::fabs(nu[t]) > amax) { amax = std::fabs(nu[t]); ax = t; }
    if (!(amax > 1e-9)) return 0;
    const int a1 = (ax == 0) ? 1 : 0;
    const int a2 = (ax == 2) ? 1 : 2;
    const double inv = 1.0 / nu[ax];
    struct Item {
        double d;
        int n2;
        Hkl x;
    };
    std::vector<Item> items;
    for (int i = -lim[a1]; i <= lim[a1]; ++i)
        for (int j = -lim[a2]; j <= lim[a2]; ++j) {
            int x[3] = {0, 0, 0};
            x[a1] = i;
            x[a2] = j;
            const double s = nu[0] * x[0] + nu[1] * x[1] + nu[2] * x[2];
            const double v = (qz - s) * inv;
            const int iv = (int)std::lround(v);
            if (iv < -lim[ax] || iv > lim[ax]) continue;
            const double d = std::fabs((iv - v) * nu[ax]);
            if (d > tol) continue;
            x[ax] = iv;
            if (x[0] == 0 && x[1] == 0 && x[2] == 0) continue;
            items.push_back({d, x[0] * x[0] + x[1] * x[1] + x[2] * x[2],
                             Hkl{x[0], x[1], x[2]}});
        }
    const int n = (int)items.size();
    if (n == 0) return 0;
    const int keep = std::min(n, std::max(1, maxOut));
    std::partial_sort(items.begin(), items.begin() + keep, items.end(),
                      [](const Item& a, const Item& b) {
                          if (a.n2 != b.n2) return a.n2 < b.n2;
                          return a.d < b.d;
                      });
    for (int i = 0; i < keep; ++i) out[i] = items[(size_t)i].x;
    return keep;
}

// Solve min||A x − b|| via normal equations (A: m*n row-major, n <= 6, no heap allocation)
bool solveLSArr(const double* A, const double* b, int m, int n, double* x) {
    double M[6][6] = {{0}};
    double v[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < m; ++i)
        for (int p = 0; p < n; ++p) {
            const double ap = A[(size_t)i * n + p];
            v[p] += ap * b[i];
            for (int q = p; q < n; ++q) M[p][q] += ap * A[(size_t)i * n + q];
        }
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < p; ++q) M[p][q] = M[q][p];
    double dmax = 0.0;
    for (int p = 0; p < n; ++p) dmax = std::max(dmax, std::fabs(M[p][p]));
    if (!(dmax > 0.0)) return false;
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(M[i][k]) > std::fabs(M[piv][k])) piv = i;
        if (std::fabs(M[piv][k]) < 1e-11 * dmax) return false;   // rank deficient
        if (piv != k) {
            for (int j = 0; j < n; ++j) std::swap(M[piv][j], M[k][j]);
            std::swap(v[piv], v[k]);
        }
        for (int i = k + 1; i < n; ++i) {
            const double f = M[i][k] / M[k][k];
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) M[i][j] -= f * M[k][j];
            v[i] -= f * v[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = v[i];
        for (int j = i + 1; j < n; ++j) s -= M[i][j] * x[j];
        x[i] = s / M[i][i];
        if (!std::isfinite(x[i])) return false;
    }
    return true;
}

// Design matrix of the (C) constraint rows G tᵀ = (t·ν)ν (linear in the free metric components).
// A is written row-major as 3×nf, the row stride is given by the caller (solveMetricN passes nf for compact storage),
// b is filled by the caller as tn·ν. **A depends only on (t, crystal system), not on ν**, so it must be precomputed on the hot path.
void fillConstraintCols(const Hkl& t0, CrystalSystem sys, double* A, int stride) {
    const int nf = metricFreeCount(sys);
    const double tt[3] = {(double)t0.h, (double)t0.k, (double)t0.l};
    // Coefficients = the metric for a unit free parameter, multiplied by t.
    for (int p = 0; p < nf; ++p) {
        double e[6] = {0, 0, 0, 0, 0, 0};
        e[p] = 1.0;
        const std::array<double, 6> g = fillMetricArr(sys, e);
        const double mm[3][3] = {{g[0], g[3], g[4]}, {g[3], g[1], g[5]}, {g[4], g[5], g[2]}};
        for (int i = 0; i < 3; ++i)
            A[(size_t)i * stride + p] = mm[i][0] * tt[0] + mm[i][1] * tt[1] + mm[i][2] * tt[2];
    }
}

// Single-row fill version of (C) (used by solveMetricN): also fills b = tn·ν and returns tn
double fillConstraintRows(const Hkl& t0, const std::array<double, 3>& nu, CrystalSystem sys,
                          double* A, double* b, int stride) {
    fillConstraintCols(t0, sys, A, stride);
    const int nf = metricFreeCount(sys);
    const double tn = t0.h * nu[0] + t0.k * nu[1] + t0.l * nu[2];
    for (int i = 0; i < 3; ++i) b[i] = tn * nu[i];
    (void)nf;
    return tn;
}

// S1's "crystal-system compatibility" prescreener: precomputes the whole least-squares projection of the (C) rows as a linear function of ν.
//
// This step **must be finished during pool filtering**: for cubic nf=1, (C) forces ν ∥ (u,v,w)
// (G = kI ⇒ G tᵀ = k t must equal (t·ν)ν), while tetragonal/hexagonal nf=2 leave 1 degree of freedom.
// If this step were left only inside solveMetricN (after the beam truncation), the pool ordering would be blind to the crystal-system constraint,
// and the S1 cost Σ(Δq_z)² precisely **favours** spurious ν that satisfy the q_z equations exactly (measurement noise
// makes the true ν residual larger instead), so the true ν would be squeezed out wholesale at the beam cut.
//
// The residual judged here is the least-squares residual of the (C) rows, which is **not larger** than the joint residual in solveMetricN,
// so a candidate eliminated early would certainly be eliminated in solveMetricN too: no solution is lost.
// When nf >= 3 or the system is rank deficient (with this t, (C) no longer constrains the ν direction, e.g. the tetragonal (110) plane), no judgement is made.
//
// Implementation folds the "least-squares residual" into a quadratic form independent of b: for b = tn·ν,
//     |b − A x*|² = tn² · νᵀQν,      Q = I − A(AᵀA)⁻¹Aᵀ
// the criterion res/tn⁴ ≤ tol² becomes νᵀQν ≤ (tol·tn)² — no sqrt and no division on the hot path,
// only 6 multiply-adds (4.3M calls; the arithmetic here directly determines the S1 runtime).
struct ConstraintScreen {
    bool active = false;
    double tol = 0.05;
    Hkl t{0, 0, 1};
    double Q[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    void build(const Hkl& t0, CrystalSystem sys, double consistTol) {
        active = false;
        tol = consistTol;
        t = t0;
        const int nf = metricFreeCount(sys);
        if (nf < 1 || nf >= 3) return;   // 3 rows / >=3 unknowns: almost always satisfiable exactly, no screening power
        double A[3 * 6] = {0};
        fillConstraintCols(t0, sys, A, nf);
        double M[6][6] = {{0}};
        for (int p = 0; p < nf; ++p)
            for (int q = 0; q < nf; ++q) {
                double s = 0.0;
                for (int i = 0; i < 3; ++i) s += A[i * nf + p] * A[i * nf + q];
                M[p][q] = s;
            }
        // Compute M^-1 (nf <= 2, direct Gauss-Jordan)
        double W[6][6] = {{0}};
        for (int p = 0; p < nf; ++p) W[p][p] = 1.0;
        double dmax = 0.0;
        for (int p = 0; p < nf; ++p) dmax = std::max(dmax, std::fabs(M[p][p]));
        if (!(dmax > 0.0)) return;
        for (int k = 0; k < nf; ++k) {
            int piv = k;
            for (int i = k + 1; i < nf; ++i)
                if (std::fabs(M[i][k]) > std::fabs(M[piv][k])) piv = i;
            if (std::fabs(M[piv][k]) < 1e-11 * dmax) return;   // rank deficient: make no judgement
            if (piv != k)
                for (int j = 0; j < nf; ++j) {
                    std::swap(M[piv][j], M[k][j]);
                    std::swap(W[piv][j], W[k][j]);
                }
            const double d = M[k][k];
            for (int j = 0; j < nf; ++j) {
                M[k][j] /= d;
                W[k][j] /= d;
            }
            for (int i = 0; i < nf; ++i) {
                if (i == k) continue;
                const double f = M[i][k];
                if (f == 0.0) continue;
                for (int j = 0; j < nf; ++j) {
                    M[i][j] -= f * M[k][j];
                    W[i][j] -= f * W[k][j];
                }
            }
        }
        // Q = I − A W Aᵀ
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int p = 0; p < nf; ++p)
                    for (int q = 0; q < nf; ++q) s += A[i * nf + p] * W[p][q] * A[j * nf + q];
                Q[i * 3 + j] = (i == j ? 1.0 : 0.0) - s;
            }
        active = true;
    }

    // nu may be any vector proportional to the true ν (the criterion is scale-invariant).
    // Returns whether it is compatible; ratioOut returns the relative residual of the (C) rows (always 0 when inactive).
    bool pass(const double* nu, double& ratioOut) const {
        ratioOut = 0.0;
        const double tn = t.h * nu[0] + t.k * nu[1] + t.l * nu[2];
        if (!(tn > 0.0)) return false;
        if (!active) return true;
        const double tq = nu[0] * (Q[0] * nu[0] + Q[1] * nu[1] + Q[2] * nu[2]) +
                          nu[1] * (Q[3] * nu[0] + Q[4] * nu[1] + Q[5] * nu[2]) +
                          nu[2] * (Q[6] * nu[0] + Q[7] * nu[1] + Q[8] * nu[2]);
        const double lim = tol * tn;
        if (tq > lim * lim) return false;
        ratioOut = std::sqrt(std::max(0.0, tq)) / tn;
        return true;
    }
};

// 3 constraint rows of (C) + n rows of (A) -> linear LS for G* -> unit cell
// Fixed-array version: the combinatorial enumeration builds a 6×(3+n) design matrix on every call, so no per-call heap allocation.
bool solveMetricN(const Peak* pk, const Hkl* hkl, int n, const std::array<double, 3>& nu,
                  const Hkl& t0, CrystalSystem sys, double consistTol, Cell& cellOut) {
    const int nf = metricFreeCount(sys);
    if (nf < 1 || n < 1) return false;
    const int nRow = std::min(n, 128);
    const int m = 3 + nRow;
    double A[(3 + 128) * 6] = {0};
    double b[3 + 128] = {0};

    const double tn = fillConstraintRows(t0, nu, sys, A, b, nf);
    const double tt[3] = {(double)t0.h, (double)t0.k, (double)t0.l};

    // (A) rows: q² = x G xᵀ
    for (int i = 0; i < nRow; ++i) {
        double row[6] = {0, 0, 0, 0, 0, 0};
        metricRow(sys, (double)hkl[i].h, (double)hkl[i].k, (double)hkl[i].l, row);
        for (int p = 0; p < nf; ++p) A[(size_t)(3 + i) * nf + p] = row[p];
        b[3 + i] = q2Of(pk[i]);
    }

    double x[6] = {0, 0, 0, 0, 0, 0};
    if (!solveLSArr(A, b, m, nf, x)) return false;
    const std::array<double, 6> G = fillMetricArr(sys, x);

    // Relative residual of the constraint rows: under high-symmetry systems, ν contradicting the symmetry is eliminated here
    const double mm[3][3] = {{G[0], G[3], G[4]}, {G[3], G[1], G[5]}, {G[4], G[5], G[2]}};
    double res = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double lhs = mm[i][0] * tt[0] + mm[i][1] * tt[1] + mm[i][2] * tt[2];
        const double d = lhs - tn * nu[i];
        res += d * d;
    }
    const double scale = std::max(1e-12, tn * tn);
    if (!(std::sqrt(res) / scale <= consistTol)) return false;
    if (!positiveDefinite(G)) return false;
    if (!cellFromRecipMetric(G, cellOut)) return false;
    return cellVolume(cellOut) > 1.0;
}

// Unified posterior: uvw normalization + sub/supercell + rod-peak order + fundamental self-consistency
std::vector<Solution> postFilter(std::vector<Solution> in, const std::vector<Peak>& pk,
                                 const RodPeriod& rod, const FixedPlaneParams& prm,
                                 const Hkl& t0) {
    std::vector<Solution> out;
    out.reserve(in.size());
    for (Solution& s : in) {
        if (!(s.cell.a > 0.0) || !(s.cell.b > 0.0) || !(s.cell.c > 0.0)) continue;
        // Normalize (u,v,w) to a primitive direction: ν and all (q_xy,q_z) predictions are unchanged (hkl needs no transform),
        // only qspec becomes the fundamental — so the order derivation in alignSolution is an integer.
        const int g = gcd3(s.uvw.h, s.uvw.k, s.uvw.l);
        if (g == 0) continue;
        if (g > 1) s.uvw = Hkl{s.uvw.h / g, s.uvw.k / g, s.uvw.l / g};
        // Fixed contact plane: the user-specified (u,v,w) is a hard constraint. Any other direction (e.g. the equivalent
        // direction (0 0 1) -> (-1 0 -1) under a reduced cell basis) is not a solution that "fixes this crystal plane". A negated
        // normal describes the same plane and is normalized to t0.
        if (s.uvw.h == -t0.h && s.uvw.k == -t0.k && s.uvw.l == -t0.l)
            s.uvw = t0;
        else if (s.uvw.h != t0.h || s.uvw.k != t0.k || s.uvw.l != t0.l)
            continue;
        s.qspec = qSpecOf(reciprocalMatrix(s.cell), s.uvw);
        if (!(s.qspec > 1e-9)) continue;
        // Subcell / supercell: the observed reciprocal points must generate this reciprocal lattice (latticeIndex == 1)
        if (latticeIndex(s.hkl) > 1) continue;
        // Rod-peak order must be a positive integer — the only external check of the absolute normal scale
        if (rod.ok) {
            bool ok = true;
            for (int i : rod.peakIdx) {
                const double n = pk[(size_t)i].qz / s.qspec;
                const double nr = std::round(n);
                if (!(nr >= 1.0) || std::fabs(n - nr) > 0.15) { ok = false; break; }
            }
            if (!ok) continue;
            s.dqspec = std::fabs(s.qspec - rod.p);
            if (!(s.dqspec <= prm.dqspec)) continue;
        } else {
            s.dqspec = 0.0;
        }
        out.push_back(std::move(s));
    }
    return out;
}

// ----------------------------------------------------------------------------
//  FP-DL: deterministic search
// ----------------------------------------------------------------------------
std::vector<Solution> fpDeterministic(const std::vector<Peak>& pk,
                                      const std::vector<int>& metricIdx, const Hkl& t0,
                                      const Orient& orient, const RodPeriod& rod,
                                      const FixedPlaneParams& prm, const std::vector<Hkl>& box,
                                      const Progress& prog) {
    const int nM = (int)metricIdx.size();
    if (nM < 6) {
        prog(0.0, "Fixed contact plane: fewer than 6 ordinary peaks, skipping deterministic search");
        return {};
    }
    const int lim[3] = {std::max(1, prm.maxH), std::max(1, prm.maxK), std::max(1, prm.maxL)};

    // Sort ordinary peaks ascending by |q| (the smallest-|q| peaks have the shortest indices and the best conditioning)
    std::vector<int> ord = metricIdx;
    std::sort(ord.begin(), ord.end(), [&pk](int a, int b) { return q2Of(pk[a]) < q2Of(pk[b]); });
    std::vector<Peak> pkM((size_t)nM);
    for (int i = 0; i < nM; ++i) pkM[(size_t)i] = pk[(size_t)ord[(size_t)i]];

    // Seed index shell: |hkl|² <= seedShell; ascending by |x|² (the true solution is usually a short vector, so try short ones first)
    std::vector<Hkl> shell;
    shell.reserve(box.size());
    const int shellMax = std::max(1, prm.seedShell);
    for (const Hkl& x : box) {
        const int n2 = x.h * x.h + x.k * x.k + x.l * x.l;
        if (n2 == 0 || n2 > shellMax) continue;
        shell.push_back(x);
    }
    if ((int)shell.size() < 4) return {};
    std::sort(shell.begin(), shell.end(), [](const Hkl& a, const Hkl& b) {
        const int na = a.h * a.h + a.k * a.k + a.l * a.l;
        const int nb = b.h * b.h + b.k * b.k + b.l * b.l;
        if (na != nb) return na < nb;
        if (a.h != b.h) return a.h < b.h;
        if (a.k != b.k) return a.k < b.k;
        return a.l < b.l;
    });
    // The first seed uses only a half-space: ν → −ν with x → −x gives the same solution (cell and |hkl| unchanged)
    std::vector<Hkl> shellFirst;
    for (const Hkl& x : shell)
        if (x.h > 0 || (x.h == 0 && (x.k > 0 || (x.k == 0 && x.l > 0)))) shellFirst.push_back(x);
    if (shellFirst.empty()) return {};
    const int kFirst = std::min<int>((int)shellFirst.size(), 12);

    // Seed triples
    const int nSeedPk = std::min(std::max(3, prm.seedPeaks), nM);
    std::vector<std::array<int, 3>> triples;
    for (int i = 0; i < nSeedPk; ++i)
        for (int j = i + 1; j < nSeedPk; ++j)
            for (int k = j + 1; k < nSeedPk; ++k) triples.push_back({i, j, k});
    if (triples.empty()) return {};

    const int S = (int)shell.size();
    const long long ops = (long long)triples.size() * (long long)kFirst * (long long)S * (long long)S;
    if (ops > (long long)std::max(1000, prm.maxCands) * 100LL) {
        prog(0.0, "Fixed contact plane: enumeration too large, reduce the hkl range / seed peak count / seed shell");
        return {};
    }

    const double qzTol = std::max(prm.qzTolMul * prm.acceptTol, 0.02);
    const double loose = 3.0 * prm.acceptTol;
    const int probe = std::min(nM, std::max(3, prm.minRowsMetric) + 4);

    // The pool stores sort scores **normalized to the tolerance** (dimensionless, see the end of S1), not raw q_z residuals:
    // raw residuals would always rank the "spurious ν that satisfy the q_z equations exactly" first.
    struct NuCand {
        std::array<double, 3> nu{0, 0, 0};
        double cost = 0.0;
    };

    // ---- S1: three peaks -> ν (3×3 integer Cramer) -----------------------------------------
    std::vector<NuCand> pool;
    std::mutex mu;
    std::atomic<int> done{0};
    const int nTasks = (int)triples.size() * kFirst;
    const int beamKeep = std::max(16, prm.beam);
    const auto tS1 = std::chrono::steady_clock::now();
    std::atomic<long long> nGatePass{0}, nProbePk{0};
    ConstraintScreen screen;
    screen.build(t0, prm.system, prm.consistTol);

    // ---- Cubic fast path: bucket the "index enumeration of the 2nd and 3rd seeds" by m = t·x ----
    // For cubic, G = kI ⇒ (C) becomes k·t = (t·ν)ν ⇒ **ν ∥ t**, hence q_z,i = λ·(t·x_i),
    // and the (C) compatibility of the three peaks depends only on m_i = t·x_i, independent of the in-plane component of the indices:
    //     q_z,1/m_1 = q_z,2/m_2 = q_z,3/m_3 = λ
    // So given (three peaks, first index x1), m2/m3 are pinned by q_z, and only indices in those buckets need enumeration.
    // The inner enumeration drops from S² (80²=6400) to a few m buckets (on the order of 9 indices in this case).
    // Note the first seed's q_z may be ≈ 0 (in-plane peak); λ is then a small quantity rather than a fixed value,
    // so we use a "feasible interval of λ" instead of a single λ: a solution exists only if the three intervals intersect, covering both cases.
    // Other crystal systems (nf=2 leaves 1 degree of freedom, nf>=3 is unconstrained) lack the "m only" structure and use generic enumeration.
    const bool cubicFast = (prm.system == CrystalSystem::Cubic) && screen.active;
    std::vector<int> shellM, bucketOff, bucketItem, allIdx;
    int mMax = 0;
    if (cubicFast) {
        shellM.resize((size_t)S);
        for (int i = 0; i < S; ++i) {
            shellM[(size_t)i] = t0.h * shell[(size_t)i].h + t0.k * shell[(size_t)i].k +
                                t0.l * shell[(size_t)i].l;
            mMax = std::max(mMax, std::abs(shellM[(size_t)i]));
        }
        bucketOff.assign((size_t)(2 * mMax + 2), 0);
        for (int i = 0; i < S; ++i) ++bucketOff[(size_t)(shellM[(size_t)i] + mMax + 1)];
        for (int m = 0; m <= 2 * mMax; ++m) bucketOff[(size_t)(m + 1)] += bucketOff[(size_t)m];
        bucketItem.resize((size_t)S);
        std::vector<int> cursor(bucketOff.begin(), bucketOff.end() - 1);
        for (int i = 0; i < S; ++i)
            bucketItem[(size_t)cursor[(size_t)(shellM[(size_t)i] + mMax)]++] = i;
    }
    allIdx.resize((size_t)S);
    for (int i = 0; i < S; ++i) allIdx[(size_t)i] = i;

    parallelForEach(nTasks, [&](int task) {
        if (prog.stop()) return;
        const std::array<int, 3>& tr = triples[(size_t)(task / kFirst)];
        const Hkl& x1 = shellFirst[(size_t)(task % kFirst)];
        const double rhs[3] = {pkM[(size_t)tr[0]].qz, pkM[(size_t)tr[1]].qz,
                               pkM[(size_t)tr[2]].qz};
        const int a[3] = {x1.h, x1.k, x1.l};
        std::vector<NuCand> local;

        // Index candidates for the 2nd/3rd seeds: the whole shell in the generic case; only m buckets on the cubic fast path
        const int* p2 = nullptr;
        const int* p3 = nullptr;
        int n2c = 0, n3c = 0;
        std::vector<int> sel2, sel3;
        bool fastUse = false;
        if (cubicFast) {
            // Each peak yields a feasible interval for λ: λ_i = q_z,i/m_i ± qzTol/|m_i|; the three intervals must intersect,
            // otherwise this (x1, b, c) cannot satisfy (C) and the whole combination needs no enumeration.
            // The intersection upper bound is folded back into an integer range for m3 → only a few buckets are visited instead of the whole shell.
            const int m1 = t0.h * a[0] + t0.k * a[1] + t0.l * a[2];
            if (m1 != 0) {
                double llo = (rhs[0] - qzTol) / (double)m1;
                double lhi = (rhs[0] + qzTol) / (double)m1;
                if (llo > lhi) std::swap(llo, lhi);
                // Interval division num/den (den spanning 0 ⇒ unbounded)
                auto divRange = [](double nlo, double nhi, double dlo, double dhi, double& lo,
                                   double& hi) {
                    if (dlo <= 0.0 && dhi >= 0.0) {
                        lo = -1e30;
                        hi = 1e30;
                        return;
                    }
                    const double c0 = nlo / dlo, c1 = nlo / dhi;
                    const double c2 = nhi / dlo, c3 = nhi / dhi;
                    lo = std::min(std::min(c0, c1), std::min(c2, c3));
                    hi = std::max(std::max(c0, c1), std::max(c2, c3));
                };
                auto pushBucket = [&](int ml, int mh, std::vector<int>& out) {
                    ml = std::max(ml, -mMax);
                    mh = std::min(mh, mMax);
                    for (int m = ml; m <= mh; ++m)
                        for (int k = bucketOff[(size_t)(m + mMax)];
                             k < bucketOff[(size_t)(m + mMax + 1)]; ++k)
                            out.push_back(bucketItem[(size_t)k]);
                };
                sel2.clear();
                sel3.clear();
                for (int m2 = -mMax; m2 <= mMax; ++m2) {
                    double alo = llo, ahi = lhi;
                    if (m2 == 0) {
                        // λ·0 must explain this peak's q_z (i.e. it must fall within |q_z| <= qzTol)
                        if (std::fabs(rhs[1]) > qzTol) continue;
                    } else {
                        const double c = rhs[1] / (double)m2;
                        const double hw = qzTol / std::fabs((double)m2);
                        alo = std::max(alo, c - hw);
                        ahi = std::min(ahi, c + hw);
                    }
                    if (!(alo <= ahi)) continue;
                    double m3lo = 0.0, m3hi = 0.0;
                    divRange(rhs[2] - qzTol, rhs[2] + qzTol, alo, ahi, m3lo, m3hi);
                    if (m3hi < (double)(-mMax - 1) || m3lo > (double)(mMax + 1)) continue;
                    m3lo = std::max(m3lo, (double)(-mMax));
                    m3hi = std::min(m3hi, (double)(mMax));
                    const int i3lo = (int)std::ceil(m3lo - 1e-9);
                    const int i3hi = (int)std::floor(m3hi + 1e-9);
                    if (i3lo > i3hi) continue;
                    pushBucket(m2, m2, sel2);
                    pushBucket(i3lo, i3hi, sel3);
                }
                std::sort(sel3.begin(), sel3.end());
                sel3.erase(std::unique(sel3.begin(), sel3.end()), sel3.end());
                // The three λ intervals do not intersect pairwise ⇒ this (x1, three peaks) yields no combination satisfying (C),
                // and not even "scan the whole shell" is needed (exactly the usual state of useless triples with in-plane seed peaks).
                if (sel2.empty() || sel3.empty()) return;
                fastUse = (long long)sel2.size() * (long long)sel3.size() <
                          (long long)S * (long long)S;
            }
        }
        if (fastUse) {
            p2 = sel2.data();
            n2c = (int)sel2.size();
            p3 = sel3.data();
            n3c = (int)sel3.size();
        } else {
            p2 = allIdx.data();
            n2c = (int)allIdx.size();
            p3 = p2;
            n3c = n2c;
        }
        for (int i2 = 0; i2 < n2c; ++i2) {
            const int b[3] = {shell[(size_t)p2[i2]].h, shell[(size_t)p2[i2]].k,
                              shell[(size_t)p2[i2]].l};
            for (int i3 = 0; i3 < n3c; ++i3) {
                const int c[3] = {shell[(size_t)p3[i3]].h, shell[(size_t)p3[i3]].k,
                                  shell[(size_t)p3[i3]].l};
                // ν = w/det. Screening starts from w (scale-invariant); division is reserved for the few candidates that pass.
                double w[3];
                int det = cramerNumerators(a, b, c, rhs, w);
                if (det == 0) continue;                 // integer rank prefilter
                if (det < 0) {                          // ν ~ −w are equivalent; normalize to the det > 0 branch
                    det = -det;
                    w[0] = -w[0];
                    w[1] = -w[1];
                    w[2] = -w[2];
                }
                // The fundamental p = t·ν must be positive (ν and −ν differ only by an overall negation; p > 0 fixes a unique branch)
                const double pw = t0.h * w[0] + t0.k * w[1] + t0.l * w[2];
                if (!(pw > 0.0)) continue;
                // The crystal-system constraint must act during pool filtering: for cubic/tetragonal/hexagonal, (C) has already pinned the ν direction,
                // so this step rejects many spurious ν and guarantees that candidates passing here also pass in solveMetricN.
                double consRatio = 0.0;
                if (!screen.pass(w, consRatio)) continue;
                const double inv = 1.0 / (double)det;
                const double nu[3] = {w[0] * inv, w[1] * inv, w[2] * inv};
                // Compatible with the measured rod fundamental: p ≈ pMeas / n (this is exactly where the "order" belongs)
                const double p = pw * inv;
                if (rod.ok && !rodConsistent(p, rod.p, prm.dqspec, prm.rodNmax)) continue;
                ++nGatePass;
                // The few smallest peaks must all have a q_z-compatible index (cheap, and it blocks the long tail)
                double cost = 0.0;
                bool ok = true;
                for (int i = 0; i < probe && ok; ++i) {
                    double d = 0.0;
                    ok = bestQzCand(lim, {nu[0], nu[1], nu[2]}, pkM[(size_t)i].qz, qzTol, d,
                                    nullptr);
                    if (ok) cost += d * d;
                }
                if (!ok) continue;
                // Sort score = sum of squares of two "residual/tolerance" ratios:
                //   1) relative residual of the (C) rows — agreement with the fixed plane + crystal system (identically 0 for low-symmetry systems,
                //      automatically degrading to pure q_z sorting);
                //   2) q_z residual — agreement with the actual observations.
                // Sorting by 2) alone would be flooded by spurious ν that "satisfy the q_z equations exactly" (their residual is 0,
                // while the true ν carries measurement noise), and the beam cut would squeeze out the true ν wholesale, so 1) must be in the ordering.
                const double qzScore = cost / (std::max(1, probe) * qzTol * qzTol);
                const double cScore = consRatio / prm.consistTol;
                if (hasEnv("GID_NOEMIT")) { ++nProbePk; continue; }   // temporary: measure the effect of push_back
                local.push_back(NuCand{{nu[0], nu[1], nu[2]}, qzScore + cScore * cScore});
            }
        }
        if (!local.empty()) {
            const size_t keep = 256;
            if (local.size() > keep) {
                std::nth_element(local.begin(), local.begin() + (ptrdiff_t)keep, local.end(),
                                 [](const NuCand& x, const NuCand& y) { return x.cost < y.cost; });
                local.resize(keep);
            }
            std::lock_guard<std::mutex> lk(mu);
            pool.insert(pool.end(), local.begin(), local.end());
        }
        const int cur = ++done;
        if (prog.report && (cur % 128) == 0)
            prog(0.4 * (double)cur / (double)nTasks, "Fixed contact plane: enumerating seed solutions");
    }, prm.threads);

    if (prog.stop() || pool.empty()) return {};
    std::sort(pool.begin(), pool.end(),
              [](const NuCand& x, const NuCand& y) { return x.cost < y.cost; });
    if (hasEnv("GID_DEBUG")) {
        std::fprintf(stderr, "[fixedplane] poolRaw=%zu score@0=%.3e score@N=%.3e\n", pool.size(),
                     pool.front().cost,
                     (int)pool.size() > beamKeep ? pool[(size_t)beamKeep - 1].cost
                                                 : pool.back().cost);
        for (int i = 0; i < 5 && i < (int)pool.size(); ++i)
            std::fprintf(stderr, "    nu[%d]=(%.4f %.4f %.4f) score=%.3e\n", i, pool[i].nu[0],
                         pool[i].nu[1], pool[i].nu[2], pool[i].cost);
        double bestDev = 1e30;
        size_t bestIdx = 0;
        const double tt2 = (double)(t0.h * t0.h + t0.k * t0.k + t0.l * t0.l);
        for (size_t i = 0; i < pool.size(); ++i) {
            const std::array<double, 3>& v = pool[i].nu;
            const double tn = t0.h * v[0] + t0.k * v[1] + t0.l * v[2];
            const double perp2 = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) - tn * tn / tt2;
            if (tn > 0.0 && perp2 < bestDev) { bestDev = perp2; bestIdx = i; }
        }
        std::fprintf(stderr, "    closest-to-t: idx=%zu/%zu dev=%.5f score=%.3e nu=(%.4f %.4f %.4f)\n",
                     bestIdx, pool.size(), std::sqrt(std::max(0.0, bestDev)), pool[bestIdx].cost,
                     pool[bestIdx].nu[0], pool[bestIdx].nu[1], pool[bestIdx].nu[2]);
    }
    if ((int)pool.size() > beamKeep) pool.resize((size_t)beamKeep);

    const bool dbg = hasEnv("GID_DEBUG");
    if (dbg) {
        const double msS1 = std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - tS1).count();
        std::fprintf(stderr,
                     "[fixedplane] nM=%d shell=%d first=%d ops=%lld pool=%zu gatePass=%lld "
                     "S1=%.0fms\n",
                     nM, S, kFirst, ops, pool.size(), nGatePass.load(), msS1);
    }

    // ---- S2/S3: solve G + alternate "nearest-neighbour assignment + linear re-solve" -------
    //
    // Runs in two stages: the **seed cells** solved from the thousands of ν in the S1 beam are highly redundant (after quantizing to
    // 1e-3 Å / 1e-2° usually only single-digit distinct solutions remain), while polishing each seed costs
    // 3 rounds of assignHkl (each rebuilding and sorting a table of 729 theoretical points) + refineAssign;
    // measured: 6000 seeds take 1.0 s, 95% of which is duplicated work. Hence:
    //   2a) enumerate "ν × 5-peak index combinations" in parallel -> seed cell + quick residual gate, with global quantized deduplication
    //   2b) alternate assignment + LM polish only on the deduplicated seeds (usually a few dozen)
    constexpr int kCandMax = 3;
    constexpr int kRowSeed = 5;
    constexpr int kQuickCand = 32;   // candidate table length used by the quick residual gate (see below)
    constexpr int kSeedCap = 4096;   // upper bound on deduplicated seeds: keeps 2b bounded
    auto cellKey = [](const Cell& c) {
        return std::array<long long, 6>{std::llround(c.a * 1000.0),  std::llround(c.b * 1000.0),
                                        std::llround(c.c * 1000.0),  std::llround(c.al * 100.0),
                                        std::llround(c.be * 100.0), std::llround(c.ga * 100.0)};
    };

    std::vector<Cell> seedCells;
    seedCells.reserve(256);
    std::atomic<int> nCombo{0}, nSolve{0}, nSeed{0}, nNoCand{0};
    const auto tS2 = std::chrono::steady_clock::now();
    const int nPool = (int)pool.size();
    {
        std::mutex muS;
        std::set<std::array<long long, 6>> seen;
        parallelFor(0, nPool, [&](int b, int e) {
            std::vector<Cell> local;
            for (int ci = b; ci < e; ++ci) {
                if (prog.stop()) return;
                const std::array<double, 3> nu = pool[(size_t)ci].nu;

                // ---- Minimal sample: the first kRowSeed smallest peaks × their q_z-compatible indices ----
                // Key 1: q_z = ν·x constrains only the normal component; in-plane (h,k) are completely free. Taking only "the nearest
                //   candidate" would point the (A) rows at wrong in-plane indices and send the metric flying.
                // Key 2: **3 peaks are not enough**. Triclinic has nf = 6; beyond the 3 constraint rows there are only 3 rows of (A),
                //   and peaks whose q² equations are proportional, like (0,0,1) and (0,0,2), immediately make the normal equations rank
                //   deficient, so solveLSArr's rank criterion rejects the true solution (this is why the measured DIP cases failed entirely).
                //   Taking 5 peaks -> 3 + 5 = 8 rows for 6 unknowns: full rank plus 2 degrees of freedom for discrimination.
                if (nM < kRowSeed) {
                    ++nNoCand;
                    continue;
                }
                Hkl cand[kRowSeed][kCandMax];
                int nc[kRowSeed] = {0, 0, 0, 0, 0};
                bool candOk = true;
                for (int s = 0; s < kRowSeed; ++s) {
                    nc[s] = qzCands(lim, nu, pkM[(size_t)s].qz, qzTol, cand[s], kCandMax);
                    if (nc[s] == 0) {
                        candOk = false;
                        break;
                    }
                }
                if (!candOk) {
                    ++nNoCand;
                    continue;
                }
                const bool trace = dbg && std::getenv("GID_TRACE") != nullptr && ci < 40;
                if (trace) {
                    std::fprintf(stderr, "  ci=%d nu=(%.4f %.4f %.4f) score=%.3e\n", ci, nu[0],
                                 nu[1], nu[2], pool[(size_t)ci].cost);
                    for (int s = 0; s < kRowSeed; ++s) {
                        std::fprintf(stderr, "    pk%d q=(%.4f,%.4f) nc=%d:", s,
                                     pkM[(size_t)s].qxy, pkM[(size_t)s].qz, nc[s]);
                        for (int t = 0; t < nc[s]; ++t)
                            std::fprintf(stderr, " (%d %d %d)", cand[s][t].h, cand[s][t].k,
                                         cand[s][t].l);
                        std::fprintf(stderr, "\n");
                    }
                }

                int nLocal = 0;
                int idx[kRowSeed] = {0, 0, 0, 0, 0};
                for (;;) {
                    Hkl x[kRowSeed];
                    for (int s = 0; s < kRowSeed; ++s) x[s] = cand[s][idx[s]];
                    if (dbg) ++nCombo;
                    Cell c3;
                    if (solveMetricN(pkM.data(), x, kRowSeed, nu, t0, prm.system, prm.consistTol,
                                     c3)) {
                        if (dbg) ++nSolve;
                        // Quick residual: for the next few peaks take the smallest (q_xy,q_z) residual among their q_z-compatible candidates.
                        // The candidate table here **must not be truncated too short**: the criterion is a min, so more candidates only make it
                        // looser and can never kill the true solution; truncating to 6 would exclude long indices (e.g. (2,1,0), n2=5, when 6
                        // shorter indices with the same q_z already precede it) and this gate would throw away the correct cell
                        // (this is why cubic-012 failed entirely).
                        const std::array<double, 6> G3 = recipMetric(c3);
                        const std::array<double, 3> nu3 = orient.get(G3);
                        bool ok = true;
                        for (int i = kRowSeed; i < std::min(nM, 10) && ok; ++i) {
                            Hkl cc[kQuickCand];
                            const int ncc = qzCands(lim, nu, pkM[(size_t)i].qz, qzTol, cc,
                                                    kQuickCand);
                            double best = 1e30;
                            for (int t = 0; t < ncc; ++t) {
                                double qx = 0.0, qz = 0.0;
                                forwardNu(G3, nu3, cc[t].h, cc[t].k, cc[t].l, qx, qz);
                                const double dx = qx - pkM[(size_t)i].qxy;
                                const double dz = qz - pkM[(size_t)i].qz;
                                best = std::min(best, std::sqrt(dx * dx + dz * dz));
                            }
                            if (!(best <= qzTol)) ok = false;
                        }
                        if (trace)
                            std::fprintf(stderr,
                                         "      x=(%d %d %d)(%d %d %d)(%d %d %d)(%d %d %d)(%d %d "
                                         "%d) cell a=%.3f b=%.3f c=%.3f al=%.1f be=%.1f ga=%.1f "
                                         "quick=%d\n",
                                         x[0].h, x[0].k, x[0].l, x[1].h, x[1].k, x[1].l, x[2].h,
                                         x[2].k, x[2].l, x[3].h, x[3].k, x[3].l, x[4].h, x[4].k,
                                         x[4].l, c3.a, c3.b, c3.c, c3.al, c3.be, c3.ga, (int)ok);
                        if (ok) {
                            ++nSeed;
                            local.push_back(c3);
                            if (++nLocal >= 4) break;   // enough seeds for this ν (the rest are near-duplicates)
                        }
                    }
                    // mixed-radix carry
                    int s = 0;
                    for (; s < kRowSeed; ++s) {
                        if (++idx[s] < nc[s]) break;
                        idx[s] = 0;
                    }
                    if (s == kRowSeed) break;   // combination enumeration finished
                }
            }
            if (!local.empty()) {
                std::lock_guard<std::mutex> lk(muS);
                for (const Cell& c : local) {
                    if (seedCells.size() >= (size_t)kSeedCap) break;
                    if (seen.insert(cellKey(c)).second) seedCells.push_back(c);
                }
            }
        }, prm.threads);
    }

    // ---- S3: alternate "nearest-neighbour assignment + linear re-solve"; only survivors reach LM ----
    //
    // Before LM, do one **all-peak residual sort and truncation**: a single assignHkl (all peaks) already suffices to pick out
    // the slightly biased seeds that are "nearly duplicates w.r.t. the factor set", and is more than an order of magnitude cheaper
    // than running assignHkl+LM for every seed. The seeds have already passed the 10-peak quick residual gate, so this only trims the long tail.
    std::vector<Solution> out;
    {
        struct Ranked {
            double rms;
            int idx;
        };
        std::vector<Ranked> ranked((size_t)seedCells.size());
        parallelFor(0, (int)seedCells.size(), [&](int b, int e) {
            std::vector<Hkl> ah;
            for (int si = b; si < e; ++si)
                ranked[(size_t)si] = {assignHkl(seedCells[(size_t)si], orient, pkM, box, ah), si};
        }, prm.threads);
        std::sort(ranked.begin(), ranked.end(),
                  [](const Ranked& a, const Ranked& b) { return a.rms < b.rms; });
        const int kPolishMax = 64;
        int nPolish = std::min<int>((int)ranked.size(), kPolishMax);
        while (nPolish > 0 && !(ranked[(size_t)nPolish - 1].rms < loose)) --nPolish;
        std::vector<Cell> polishSeeds;
        polishSeeds.reserve((size_t)nPolish);
        for (int i = 0; i < nPolish; ++i) polishSeeds.push_back(seedCells[(size_t)ranked[(size_t)i].idx]);

        std::mutex mu2;
        const int nPolishSeeds = (int)polishSeeds.size();
        parallelFor(0, nPolishSeeds, [&](int b, int e) {
            std::vector<Solution> local;
            std::vector<Hkl> ah((size_t)nM, Hkl{0, 0, 0});
            for (int si = b; si < e; ++si) {
                if (prog.stop()) return;
                Cell cell = polishSeeds[(size_t)si];
                double rms = 1e30;
                bool alive = true;
                for (int it = 0; it < 3 && alive; ++it) {
                    rms = assignHkl(cell, orient, pkM, box, ah);
                    if (!(rms < loose)) { alive = false; break; }
                    const std::array<double, 3> nuNow = orient.get(recipMetric(cell));
                    if (!solveMetricN(pkM.data(), ah.data(), nM, nuNow, t0, prm.system,
                                      prm.consistTol, cell))
                        alive = false;
                }
                if (!alive || !(rms < loose)) continue;

                double r2 = 0.0;
                const Cell cell2 = refineAssign(pkM, orient, box, cell, prm.system, 5, r2);
                const bool use2 = r2 < rms;
                const Cell best = use2 ? cell2 : cell;
                if (!((use2 ? r2 : rms) < prm.acceptTol)) continue;

                // All-peak (including rod peaks) assignment: rod peaks are pinned to n·t0 here, guarding the normal scale
                std::vector<Hkl> ahAll;
                std::vector<double> distAll;
                assignHkl(best, orient, pk, box, ahAll, &distAll);
                const Quality q = evaluate(best, orient, pk, ahAll);
                if (!(q.rmsQxyz < prm.acceptTol)) continue;

                Solution s;
                s.uvw = t0;
                s.nu = orient.get(recipMetric(best));
                s.cell = best;
                s.vol = cellVolume(best);
                s.qspec = qSpecOf(reciprocalMatrix(best), t0);
                s.hkl = std::move(ahAll);
                s.dist = std::move(distAll);
                s.rmsQxy = q.rmsQxy;
                s.rmsQz = q.rmsQz;
                s.rmsQxyz = q.rmsQxyz;
                local.push_back(std::move(s));
            }
            if (!local.empty()) {
                std::lock_guard<std::mutex> lk(mu2);
                out.insert(out.end(), local.begin(), local.end());
            }
        }, prm.threads);
    }

    if (dbg) {
        const double msS2 = std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - tS2).count();
        std::fprintf(stderr,
                     "[fixedplane] pool=%d noCand=%d combo=%d solve=%d seed=%d seeds=%zu out=%zu "
                     "nM=%d qzTol=%.5f loose=%.5f acceptTol=%.5f S2=%.0fms\n",
                     nPool, nNoCand.load(), nCombo.load(), nSolve.load(), nSeed.load(),
                     seedCells.size(), out.size(), nM, qzTol, loose, prm.acceptTol, msS2);
    }

    return out;
}

} // namespace

std::vector<Solution> runFixedPlane(const std::vector<Peak>& peaksIn, const FixedPlaneParams& prm,
                                    const Progress& prog) {
    std::vector<Peak> pk = peaksIn;
    // The sign of q_xy reflects only the mirror handedness of the detection geometry; forwardNu on the model side always yields non-negative q_xy.
    // Data is forced to |q_xy| (normalizeQxySigns) before solving, so it no longer depends on the autoSign switch:
    // that switch is kept only for compatibility with existing UI/parameter interfaces; on and off give the same result.
    normalizeQxySigns(pk);
    if (pk.empty()) return {};

    // The contact plane keeps only its direction: the model is exactly invariant under integer multiples of uvw (ν = G nᵀ/√(n G nᵀ)),
    // and after normalization qspec is the rod fundamental, so the order derivation in alignSolution is an integer.
    const int g0 = gcd3(prm.u, prm.v, prm.w);
    if (g0 == 0) return {};
    const Hkl t0{prm.u / g0, prm.v / g0, prm.w / g0};
    const Orient orient = Orient::fromUvw(t0);

    const std::vector<Hkl> box =
        hklGrid(std::max(1, prm.maxH), std::max(1, prm.maxK), std::max(1, prm.maxL));
    if (box.empty()) return {};

    // Specular rod: fundamental + order (not min q_z)
    const RodPeriod rod = detectRodPeriod(pk, prm.rodTol, prm.rodNmax);
    std::vector<int> metricIdx;                     // ordinary peaks (q_xy > rodTol)
    for (int i = 0; i < (int)pk.size(); ++i)
        if (!(std::fabs(pk[i].qxy) <= prm.rodTol)) metricIdx.push_back(i);

    // ---------------------------------------------------------------- Path 0
    std::vector<Solution> accepted;
    {
        Cell cell = prm.start;
        applySystemConstraints(prm.system, cell);
        double rms = 0.0;
        cell = refineAssign(pk, orient, box, cell, prm.system, prm.refineIter, rms);
        std::vector<Hkl> ah;
        std::vector<double> dist;
        assignHkl(cell, orient, pk, box, ah, &dist);
        const Quality q = evaluate(cell, orient, pk, ah);
        Solution s;
        s.uvw = t0;
        s.nu = orient.get(recipMetric(cell));
        s.cell = cell;
        s.vol = cellVolume(cell);
        s.qspec = qSpecOf(reciprocalMatrix(cell), t0);
        s.hkl = ah;
        s.dist = dist;
        s.rmsQxy = q.rmsQxy;
        s.rmsQz = q.rmsQz;
        s.rmsQxyz = q.rmsQxyz;
        std::vector<Solution> ok = postFilter({std::move(s)}, pk, rod, prm, t0);
        if (!ok.empty()) {
            if (ok.front().rmsQxyz <= 0.3 * prm.acceptTol) {
                prog(1.0, "Fixed contact plane: initial-cell refinement complete");
                return ok;
            }
            // If refinement did not converge, it **must not** enter the results: an initial LM that converges to a local minimum (e.g. cubic
            // a=8.25 with rms=0.13) would impersonate "Top-1" in the final ranking and squeeze out the true solution
            // found by FP-DL below. The same acceptance threshold as fpDeterministic is used here.
            if (ok.front().rmsQxyz <= prm.acceptTol)
                accepted.insert(accepted.end(), ok.begin(), ok.end());
        }
    }

    // ---------------------------------------------------------------- Path 1
    // With (u,v,w) known, the two-step method is a strong constraint (every 3 peaks + the specular condition solve 3 unknowns),
    // but its q_spec must be the rod **fundamental**: when the first-order peak is missing, try p/n as well.
    if (rod.ok && (int)pk.size() >= 4) {
        std::vector<double> qsTries{rod.p};
        if (rod.ambiguous)
            for (int n = 2; n <= std::max(2, prm.rodNmax); ++n) qsTries.push_back(rod.p / (double)n);
        for (double qs : qsTries) {
            if (prog.stop()) break;
            SpecParams sp;
            sp.qspec = qs;
            sp.part1Alg = Part1Algorithm::Improved;
            sp.fixedUvw = true;
            sp.fu = t0.h;
            sp.fv = t0.k;
            sp.fw = t0.l;
            const int hkAbs = std::max(std::abs(prm.maxH), std::abs(prm.maxK));
            sp.hkLimit1 = std::clamp(hkAbs, 1, 4);
            sp.hkLimit2 = std::clamp(hkAbs, 2, 8);
            sp.lLimit = std::max(1, std::abs(prm.maxL));
            sp.linesToPermute = std::max(3, std::min(5, (int)pk.size()));
            sp.redundant = 4;
            sp.maxPart1 = 40;
            sp.refineTop = 200;
            sp.maxSolutions = std::max(1, prm.topN);
            sp.dqxyzCutoff = std::max(1e-6, prm.acceptTol);
            sp.dqspec = std::max(1e-6, prm.dqspec);
            sp.threads = prm.threads;
            std::vector<Solution> sols = improved::runSpecIndex(pk, sp, prog);
            sols = postFilter(std::move(sols), pk, rod, prm, t0);
            if (!sols.empty()) {
                prog(1.0, "Fixed contact plane: two-step method complete");
                return sols;
            }
        }
    }

    // ---------------------------------------------------------------- Path 2
    {
        std::vector<Solution> sols =
            fpDeterministic(pk, metricIdx, t0, orient, rod, prm, box, prog);
        sols = postFilter(std::move(sols), pk, rod, prm, t0);
        accepted.insert(accepted.end(), sols.begin(), sols.end());
    }

    // ------------------------------------------------- sort + quantized dedup + Top-N
    std::sort(accepted.begin(), accepted.end(), [](const Solution& a, const Solution& b) {
        if (a.rmsQxyz != b.rmsQxyz) return a.rmsQxyz < b.rmsQxyz;
        const Cell& x = a.cell;
        const Cell& y = b.cell;
        if (x.a != y.a) return x.a < y.a;
        if (x.b != y.b) return x.b < y.b;
        if (x.c != y.c) return x.c < y.c;
        if (x.al != y.al) return x.al < y.al;
        if (x.be != y.be) return x.be < y.be;
        return x.ga < y.ga;
    });
    std::vector<Solution> out;
    std::set<std::array<long long, 6>> seen;
    for (Solution& s : accepted) {
        const double sc = 1000.0;
        std::array<long long, 6> key{std::llround(s.cell.a * sc), std::llround(s.cell.b * sc),
                                     std::llround(s.cell.c * sc), std::llround(s.cell.al * sc),
                                     std::llround(s.cell.be * sc), std::llround(s.cell.ga * sc)};
        if (!seen.insert(key).second) continue;
        out.push_back(std::move(s));
        if ((int)out.size() >= std::max(1, prm.topN)) break;
    }
    prog(1.0, "Fixed contact plane: complete");
    return out;
}

} // namespace gid
