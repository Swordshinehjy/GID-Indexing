// ============================================================================
//  nospec.cpp — see nospec.h
// ============================================================================
#include "nospec.h"
#include "parallel.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <unordered_set>
#include <vector>

namespace gid {
namespace {

constexpr int kMaxRows = kNoSpecMaxRows;

// States are copied in the millions, so fields use the narrowest integers (indices < 128) to reduce memory-bandwidth usage
struct State {
    double nu[3] = {0, 0, 0};
    double cost = 0.0;
    // LS row residual of the (A) equations (q² = hkl·G*·hkl) in q units: meaningful only once the rows are full rank; <0 means not yet solved
    float residA = -1.0f;
    // Physical consistency |νᵀ(G*)⁻¹ν − 1|: <0 means not yet computed (computed only after full rank)
    float residC = -1.0f;
    int16_t n = 0;
    int16_t pk[kMaxRows] = {0};       // assigned peak indices
    int8_t hkl[3 * kMaxRows] = {0};   // corresponding (h,k,l)
};

// State sort key (smaller is better):
//   * not yet full rank (residA < 0): can only be sorted by the q_z residual (weak criterion), and ranks after all full-rank states;
//   * after full rank: use the "normalized worst violation" max(rms/pruneRms, |s1−1|/consistTol), putting the two criteria with different
//     dimensions on the same ruler (for sorting only, no hard elimination).
// Note: at the first full-rank row (n = nf+1, only 1 degree of freedom) rms is almost always 0, so the ranking at this level carries
// very little information — that is why the beam width must be several thousand; the cause is here, not a poorly written key (measured: changing it does not help).
inline double stateKey(const State& s, double pruneRms, double consistTol) {
    if (s.residA < 0.0f) return 1e3 + s.cost;
    const double a = (pruneRms > 0.0) ? (double)s.residA / pruneRms : (double)s.residA;
    const double c = (s.residC >= 0.0f && consistTol > 0.0) ? (double)s.residC / consistTol : 0.0;
    return std::max(a, c);
}
// 3x3 integer determinant (exact prefilter, A1)
inline int intDet3(const int a[3], const int b[3], const int c[3]) {
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
           a[2] * (b[0] * c[1] - b[1] * c[0]);
}

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

// Solve min||A x - b|| (normal equations + Gauss elimination), A: m*n row-major
bool solveLS(const std::vector<double>& A, const std::vector<double>& b, int m, int n,
             std::vector<double>& x) {
    std::vector<double> M((size_t)n * n, 0.0), v(n, 0.0);
    for (int i = 0; i < m; ++i)
        for (int p = 0; p < n; ++p) {
            const double ap = A[(size_t)i * n + p];
            v[p] += ap * b[i];
            for (int q2 = p; q2 < n; ++q2) M[(size_t)p * n + q2] += ap * A[(size_t)i * n + q2];
        }
    for (int p = 0; p < n; ++p)
        for (int q2 = 0; q2 < p; ++q2) M[(size_t)p * n + q2] = M[(size_t)q2 * n + p];
    // Relative pivot threshold: no unique G* can be solved when the normal equations are rank deficient (the peak combination does not span nf dimensions)
    double dmax = 0.0;
    for (int p = 0; p < n; ++p) dmax = std::max(dmax, std::fabs(M[(size_t)p * n + p]));
    if (!(dmax > 0.0)) return false;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(M[(size_t)i * n + k]) > std::fabs(M[(size_t)p * n + k])) p = i;
        if (std::fabs(M[(size_t)p * n + k]) < 1e-11 * dmax) return false;
        if (p != k) {
            for (int j = 0; j < n; ++j) std::swap(M[(size_t)p * n + j], M[(size_t)k * n + j]);
            std::swap(v[p], v[k]);
        }
        for (int i = k + 1; i < n; ++i) {
            const double f = M[(size_t)i * n + k] / M[(size_t)k * n + k];
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) M[(size_t)i * n + j] -= f * M[(size_t)k * n + j];
            v[i] -= f * v[k];
        }
    }
    x.assign(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double s = v[i];
        for (int j = i + 1; j < n; ++j) s -= M[(size_t)i * n + j] * x[j];
        x[i] = s / M[(size_t)i * n + i];
    }
    for (double d : x)
        if (!std::isfinite(d)) return false;
    return true;
}

// Fill the full G* from the free components (constrained components completed per crystal system). x points to the nf free components.
// Fixed-size pointer version: called for every state in the pruning inner loop, avoiding construction of a std::vector.
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

std::array<double, 6> fillMetric(CrystalSystem sys, const std::vector<double>& x) {
    return fillMetricArr(sys, x.data());
}

// Number of free parameters of the (A) equations q² = hkl·G*·hkl (matches the x length of fillMetric exactly)
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

// Design row of the (A) equations for the free metric components.
// Key: constrained components (G22=G11, G33=G11, G12=G11/2 …) must be **folded into** the free parameters,
// not dropped as in "pick only the free columns" — otherwise q² = G11·h² is completely wrong for cubic (peaks such as (0,1,0)
// can never be explained), and for hexagonal it treats 2 degrees of freedom as 3 and the matrix becomes rank deficient.
inline void metricRow(CrystalSystem sys, double h, double k, double l, double* a) {
    switch (sys) {
        case CrystalSystem::Cubic:
            a[0] = h * h + k * k + l * l;
            return;
        case CrystalSystem::Tetragonal: // G11 = G22
            a[0] = h * h + k * k;
            a[1] = l * l;
            return;
        case CrystalSystem::Hexagonal: // G11 = G22, G12 = G11/2
        case CrystalSystem::Trigonal:
            a[0] = h * h + k * k + h * k;
            a[1] = l * l;
            return;
        case CrystalSystem::Orthorhombic:
            a[0] = h * h; a[1] = k * k; a[2] = l * l;
            return;
        case CrystalSystem::Monoclinic: // the only free off-diagonal component G13
            a[0] = h * h; a[1] = k * k; a[2] = l * l; a[3] = 2.0 * h * l;
            return;
        default:
            a[0] = h * h; a[1] = k * k; a[2] = l * l;
            a[3] = 2.0 * h * k; a[4] = 2.0 * h * l; a[5] = 2.0 * k * l;
            return;
    }
}

// Exact Cramer solution for the minimal 3-peak sample of the (B) equations: the three rows of M = the hkl of the three seed peaks, the rows of adj(M) =
// the cross products of the other two rows; cross products and determinant are all integers (magnitude << 2^31), exact; the three 3×3 solves need only one division.
// det is passed in by the caller from intDet3 (also serving as the rank prefilter).
void cramerNu(const int* xa, const int* xb, const int* xc, const double rhs[3], int det,
              double out[3]) {
    const int* rows[3] = {xa, xb, xc};
    double D[3][3];
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) D[r][c] = (double)rows[r][c];
    for (int col = 0; col < 3; ++col) {
        for (int r = 0; r < 3; ++r) D[r][col] = rhs[r];
        out[col] = (D[0][0] * (D[1][1] * D[2][2] - D[1][2] * D[2][1]) -
                    D[0][1] * (D[1][0] * D[2][2] - D[1][2] * D[2][0]) +
                    D[0][2] * (D[1][0] * D[2][1] - D[1][1] * D[2][0])) / (double)det;
        for (int r = 0; r < 3; ++r) D[r][col] = (double)rows[r][col];
    }
}

// Symmetric 3x3 inverse (G* -> G*^-1), input order [G11,G22,G33,G12,G13,G23]
bool invMetric(const std::array<double, 6>& g, std::array<double, 6>& inv) {    const double d = g[0] * (g[1] * g[2] - g[5] * g[5]) - g[3] * (g[3] * g[2] - g[5] * g[4]) +
                     g[4] * (g[3] * g[5] - g[1] * g[4]);
    if (!(std::fabs(d) > 1e-30)) return false;
    const double id = 1.0 / d;
    inv[0] = (g[1] * g[2] - g[5] * g[5]) * id;
    inv[1] = (g[0] * g[2] - g[4] * g[4]) * id;
    inv[2] = (g[0] * g[1] - g[3] * g[3]) * id;
    inv[3] = (g[4] * g[5] - g[3] * g[2]) * id;
    inv[4] = (g[3] * g[5] - g[1] * g[4]) * id;
    inv[5] = (g[3] * g[4] - g[0] * g[5]) * id;
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

// ----------------------------------------------------------------------------
// B2: LLL reduction + sphere-decoding assignment
// ----------------------------------------------------------------------------

// 3x3 LLL reduction (δ=0.75). rows = 3 basis vectors (Cartesian); outputs the unimodular integer matrix U
// (recording the row operations: reduced = U · input) and the reduced rows. The sphere-decoding enumeration radius on an oblique basis is set by the longest
// basis vector; a long c axis (DIP data c≈25 Å, c* very short) is exactly the ill-conditioned case, and LLL swaps the basis for nearly orthogonal
// short vectors, tightening the bound and sharply cutting the enumeration — this is the premise of B2.
//   Convention: lattice point q = xᵀ · rows (x integer), ||q||² = xᵀ(rows·rowsᵀ)x = xᵀG*x.
//   The lattice is unchanged after reduction: x = Uᵀ x' (U unimodular ⇒ integer one-to-one correspondence).
void lll3(double rows[3][3], int U[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) U[i][j] = (i == j) ? 1 : 0;
    double bs[3][3] = {};  // Gram–Schmidt
    double mu[3][3] = {};
    const auto gso = [&]() {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) bs[i][j] = rows[i][j];
            for (int j = 0; j < i; ++j) {
                double d = 0.0;
                for (int t = 0; t < 3; ++t) d += rows[i][t] * bs[j][t];
                const double n2 = bs[j][0] * bs[j][0] + bs[j][1] * bs[j][1] + bs[j][2] * bs[j][2];
                mu[i][j] = (n2 > 1e-30) ? d / n2 : 0.0;
                for (int t = 0; t < 3; ++t) bs[i][t] -= mu[i][j] * bs[j][t];
            }
        }
    };
    gso();
    const auto norm2 = [](const double v[3]) {
        return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    };
    int k = 1;
    int guard = 0;
    while (k < 3 && ++guard < 100) {
        // size reduction: b_k -= round(μ_kj)·b_j
        for (int j = k - 1; j >= 0; --j) {
            const double r = std::floor(mu[k][j] + 0.5);
            if (r == 0.0) continue;
            for (int t = 0; t < 3; ++t) rows[k][t] -= r * rows[j][t];
            U[k][0] -= (int)r * U[j][0];
            U[k][1] -= (int)r * U[j][1];
            U[k][2] -= (int)r * U[j][2];
            gso();
        }
        // Lovász condition: if violated, swap b_k and b_{k-1} and step back
        const double lhs = norm2(bs[k]);
        const double rhs = (0.75 - mu[k][k - 1] * mu[k][k - 1]) * norm2(bs[k - 1]);
        if (lhs >= rhs) {
            ++k;
        } else {
            for (int t = 0; t < 3; ++t) std::swap(rows[k][t], rows[k - 1][t]);
            std::swap(U[k][0], U[k - 1][0]);
            std::swap(U[k][1], U[k - 1][1]);
            std::swap(U[k][2], U[k - 1][2]);
            gso();
            k = std::max(k - 1, 1);
        }
    }
}

// Cholesky of G* ([G11,G22,G33,G12,G13,G23]): G* = RᵀR, with R upper triangular.
// Returns false = not positive definite (the caller already has a positiveDefinite gate; this is only a fallback).
bool cholMetric(const std::array<double, 6>& g, double R[3][3]) {
    const double m[3][3] = {{g[0], g[3], g[4]}, {g[3], g[1], g[5]}, {g[4], g[5], g[2]}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) R[i][j] = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i; j < 3; ++j) {
            double s = m[i][j];
            for (int t = 0; t < i; ++t) s -= R[t][i] * R[t][j];
            if (i == j) {
                if (!(s > 1e-30)) return false;
                R[i][i] = std::sqrt(s);
            } else {
                R[i][j] = s / R[i][i];
            }
        }
    }
    return true;
}

// Single-peak assignment: find hkl in the **whole integer space** minimizing (Δq_z/σz)² + (Δ|q|/σq)².
//   (A) ellipsoid: ||R x||² = xᵀG*x ∈ [q²−εA, q²+εA] (the tolerance on observed |q| converted to a bandwidth on q²);
//   (B) slab: |nu·x − q_z| ≤ tolZ.
// Fincke–Pohst enumerates downward from x3 (R upper triangular ⇒ y3 contains only x3); the bandwidth uses the sphere upper bound, each candidate is re-tested.
// Returns false if the candidate count exceeds the budget (the caller falls back to the grid nearest neighbour). Returns true and writes out on success.
bool sphereAssign(const std::array<double, 6>& G, const std::array<double, 3>& nu,
                  const std::vector<Peak>& pk, int i, double tolZ, double tolQ,
                  std::array<int, 3>& out) {
    // Flow: chol(G*) gives the lattice basis rows → LLL reduction (rows mixed by U) → the transformed metric G'' = U·G*·Uᵀ
    // → Cholesky back to upper-triangular R — LLL keeps R well-conditioned, so the triangular-enumeration assumption of Fincke–Pohst holds.
    double Rb[3][3];
    if (!cholMetric(G, Rb)) return false;
    int U[3][3];
    lll3(Rb, U);
    const double gm[3][3] = {{G[0], G[3], G[4]}, {G[3], G[1], G[5]}, {G[4], G[5], G[2]}};
    double R[3][3];
    {
        double g2[3][3];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                double s = 0.0;
                for (int t = 0; t < 3; ++t)
                    for (int u = 0; u < 3; ++u) s += (double)U[a][t] * gm[t][u] * (double)U[b][u];
                g2[a][b] = s;
            }
        bool ok = true;
        for (int a = 0; a < 3 && ok; ++a) {
            for (int b = 0; b < 3; ++b) R[a][b] = 0.0;
            for (int b = a; b < 3; ++b) {
                double s = g2[a][b];
                for (int t = 0; t < a; ++t) s -= R[t][a] * R[t][b];
                if (a == b) {
                    if (!(s > 1e-30)) { ok = false; break; }
                    R[a][a] = std::sqrt(s);
                } else {
                    R[a][b] = s / R[a][a];
                }
            }
        }
        if (!ok) return false;
    }
    const double qz = pk[i].qz;
    const double qxy = pk[i].qxy;
    const double q2 = qz * qz + qxy * qxy;
    const double q = std::sqrt(std::max(q2, 1e-12));
    const double epsA = 2.0 * q * tolQ;          // Δ|q| → Δq²
    const double rHi2 = q2 + epsA;
    const double rHi = std::sqrt(std::max(rHi2, 1e-12));
    const double wz = tolZ * tolZ, wq = tolQ * tolQ;
    // x = Uᵀx' ⇒ nu·x = (U nu)·x'
    const double nuR[3] = {U[0][0] * nu[0] + U[0][1] * nu[1] + U[0][2] * nu[2],
                           U[1][0] * nu[0] + U[1][1] * nu[1] + U[1][2] * nu[2],
                           U[2][0] * nu[0] + U[2][1] * nu[1] + U[2][2] * nu[2]};
    double best = 1e300;
    std::array<int, 3> bx{0, 0, 0};
    bool found = false;
    long long budget = 20000;  // per-peak enumeration budget: on overrun, fall back to the grid (only happens when G* deviates badly)
    const int x3hi = (int)std::floor(rHi / R[2][2]);
    for (int x3 = -x3hi; x3 <= x3hi && budget > 0; ++x3) {
        const double y3 = R[2][2] * x3;
        const double r2 = rHi2 - y3 * y3;
        if (r2 < 0.0) continue;
        const double r2s = std::sqrt(r2);
        const int x2lo = (int)std::ceil((-R[1][2] * x3 - r2s) / R[1][1]);
        const int x2hi = (int)std::floor((-R[1][2] * x3 + r2s) / R[1][1]);
        for (int x2 = x2lo; x2 <= x2hi && budget > 0; ++x2) {
            const double y2 = R[1][1] * x2 + R[1][2] * x3;
            const double r1 = r2 - y2 * y2;
            if (r1 < 0.0) continue;
            const double r1s = std::sqrt(r1);
            const int x1lo = (int)std::ceil((-R[0][1] * x2 - R[0][2] * x3 - r1s) / R[0][0]);
            const int x1hi = (int)std::floor((-R[0][1] * x2 - R[0][2] * x3 + r1s) / R[0][0]);
            for (int x1 = x1lo; x1 <= x1hi && budget > 0; ++x1) {
                --budget;
                const double y1 = R[0][0] * x1 + R[0][1] * x2 + R[0][2] * x3;
                const double n2 = y1 * y1 + y2 * y2 + y3 * y3;
                const double dq = (n2 - q2) / (2.0 * q);
                if (dq > tolQ || dq < -tolQ) continue;  // bandwidth re-test (the sphere enumeration used only the upper bound)
                const double dz = nuR[0] * x1 + nuR[1] * x2 + nuR[2] * x3 - qz;
                if (dz > tolZ || dz < -tolZ) continue;
                const double score = dz * dz / wz + dq * dq / wq;
                if (score < best) {
                    best = score;
                    // x = Uᵀ x' (inverse mapping of the U row operations, see the lll3 comment)
                    bx = {U[0][0] * x1 + U[1][0] * x2 + U[2][0] * x3,
                          U[0][1] * x1 + U[1][1] * x2 + U[2][1] * x3,
                          U[0][2] * x1 + U[1][2] * x2 + U[2][2] * x3};
                    found = true;
                }
            }
        }
    }
    if (found) out = bx;
    return found;
}

} // namespace

std::vector<Solution> runNoSpec(const std::vector<Peak>& peaksIn, const NoSpecParams& prmIn,
                                const Progress& prog) {
    // P1-1: sanitize the user parameters first. The UI allows arbitrary integers; out-of-range values cause int8_t index truncation
    // or an explosion of the enumeration size, so all later logic uses only the sanitized copy.
    const NoSpecParams prm = prmIn.sanitized();
    // Computation uniformly uses |q_xy|: the sign reflects only the mirror handedness of the detection geometry, not structural information.
    std::vector<Peak> inp = peaksIn;
    normalizeQxySigns(inp);
    // Consistent with the criterion in engine::alignSolution: only peaks with q_xy>0 are processed (q_xy=0 specular rod peaks
    // do not belong to the no-specular model). Exact duplicate points (Friedel pairs) are also deduplicated:
    // duplicate peaks leave the (A) rows with very few independent q² values, making the metric equations underdetermined and producing spurious solutions.
    std::vector<int> srcPos;      // peaks of the input table participating in the search (original order)
    std::vector<Peak> uniq;       // deduplicated peaks
    std::vector<int> srcToUniq;   // srcPos position -> uniq index
    for (int i = 0; i < (int)inp.size(); ++i) {
        const Peak& p = inp[i];
        if (!(p.qxy > 1e-6) || !std::isfinite(p.qxy) || !std::isfinite(p.qz)) continue;
        int u = -1;
        for (int j = 0; j < (int)uniq.size(); ++j)
            if (std::fabs(uniq[j].qxy - p.qxy) < 1e-6 && std::fabs(uniq[j].qz - p.qz) < 1e-6) {
                u = j;
                break;
            }
        if (u < 0) {
            u = (int)uniq.size();
            uniq.push_back(p);
        }
        srcPos.push_back(i);
        srcToUniq.push_back(u);
    }
    const int M = (int)srcPos.size();
    // ascending by |g|
    std::vector<int> ord(uniq.size());
    for (int i = 0; i < (int)ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&uniq](int a, int b) {
        const double qa = uniq[a].qxy * uniq[a].qxy + uniq[a].qz * uniq[a].qz;
        const double qb = uniq[b].qxy * uniq[b].qxy + uniq[b].qz * uniq[b].qz;
        return qa < qb;
    });
    std::vector<int> posOf(uniq.size());
    std::vector<Peak> pk(uniq.size());
    for (int i = 0; i < (int)ord.size(); ++i) {
        posOf[ord[i]] = i;
        pk[i] = uniq[ord[i]];
    }
    for (int& u : srcToUniq) u = posOf[u];
    const int N = (int)pk.size();
    // nf = number of free metric components, must match the x length of metricRow / fillMetric.
    // Note: do not compute it as "how many length parameters the crystal system has": for hexagonal a=b, but G12 is bound to G11/2,
    // so only a and c are truly free; hence metricFreeCount is used uniformly.
    const int nf = metricFreeCount(prm.system);

    // P1-2: std::clamp requires lo <= hi. With few peaks the upper bound can fall below the lower bound (UB); guard the minimum peak count first.
    if (N < 3) return {};
    const int minRows = std::clamp(prm.minRows, std::max(3, nf + 2), kMaxRows);
    if (N < minRows) return {};
    const int nSeed = std::clamp(prm.seedPeaks, 3, std::max(3, std::min(N, kMaxRows)));

    // ---- index box / seed shell / expansion enumeration table
    // Box: independent range per axis (a long c axis needs l much larger than h,k, otherwise the true solution is outside the search space)
    const int lim[3] = {prm.seedH, prm.seedK, prm.seedL};
    // Seed shell: seed peaks have the smallest |g| and their indices must be short vectors, so only |hkl|² <= seedShell is taken
    std::vector<std::array<int, 3>> shell;
    for (int h = -lim[0]; h <= lim[0]; ++h)
        for (int k = -lim[1]; k <= lim[1]; ++k)
            for (int l = -lim[2]; l <= lim[2]; ++l) {
                if (h == 0 && k == 0 && l == 0) continue;
                if (h * h + k * k + l * l > prm.seedShell) continue;
                shell.push_back({h, k, l});
            }
    if ((int)shell.size() < 4) return {};
    // Global negation x -> -x gives the same solution (nu -> -nu; cell and hkl differ only in sign), so seed 1 takes only half
    std::vector<std::array<int, 3>> shellHalf;
    for (const auto& x : shell)
        if (x[0] > 0 || (x[0] == 0 && (x[1] > 0 || (x[1] == 0 && x[2] > 0))))
            shellHalf.push_back(x);
    // During expansion, find "indices in the box compatible with the current nu" per peak: use the largest |nu| component as the solving axis (one division,
    // and the least sensitive), enumerate the other two. Complexity drops from |box| to the product of the other two ranges; the long axis is no longer scanned point by point.
    std::vector<std::array<int, 3>> pairTab[3];
    for (int ax = 0; ax < 3; ++ax) {
        const int a1 = (ax == 0) ? 1 : 0;
        const int a2 = (ax == 2) ? 1 : 2;
        pairTab[ax].reserve((size_t)(2 * lim[a1] + 1) * (size_t)(2 * lim[a2] + 1));
        for (int i = -lim[a1]; i <= lim[a1]; ++i)
            for (int j = -lim[a2]; j <= lim[a2]; ++j) {
                std::array<int, 3> x{0, 0, 0};
                x[a1] = i;
                x[a2] = j;
                pairTab[ax].push_back(x);
            }
    }
    const int S = (int)shell.size();
    const int SH = (int)shellHalf.size();

    // P2-1: criterion scale factors moved to NoSpecParams (they can be adjusted uniformly when data noise changes)
    const double tol = std::max(prm.qzTolMul * prm.dqzCutoff, prm.qzTolFloor);
    // (A) pruning threshold: once the rows suffice to solve the metric, spurious states are eliminated by the least-squares row residual of the q² equations
    const double pruneRms =
        prm.metricPrune > 0.0 ? prm.metricPrune : prm.metricPruneMul * prm.dqxyzCutoff;

    // Physical-consistency tolerance |νᵀ(G*)⁻¹ν − 1|: a soft criterion during search, a hard gate in S3 (the same number)
    const double consistTol = 0.12;

    // LS solution of the (A) equations + two criteria:
    //   * metric positive definiteness: non-PD = physically impossible cell. Acts as a hard gate **only on the last row** (where it is
    //     bit-for-bit equivalent to the S3 rejection criterion); with few rows it is only a soft component of the sort key.
    //   * |νᵀG*⁻¹ν − 1|: the true solution necessarily satisfies it (|n̂|=1), spurious states generally deviate; enters the sort key.
    // solved=false means the normal equations are rank deficient (no judgement possible; the state is kept).
    struct MetricEval {
        bool solved = false;
        bool pd = false;
        double rms = 0.0;
        double consist = -1.0;
    };
    const auto evalMetric = [&](const State& st) -> MetricEval {
        MetricEval ev;
        double M[6][6] = {}, v[6] = {}, x[6] = {};
        const int m = st.n;
        // Row cache: metricRow would be computed twice (building the normal equations + computing residuals after back-substitution); with m ≤ kMaxRows, nf ≤ 6,
        // precomputing once on the stack halves the calls to this hot function.
        double rows[kMaxRows][6];
        double rhs[kMaxRows];
        for (int i = 0; i < m; ++i) {
            const double h = st.hkl[3 * i], k = st.hkl[3 * i + 1], l = st.hkl[3 * i + 2];
            double* a = rows[i];
            metricRow(prm.system, h, k, l, a);
            const double bi =
                pk[st.pk[i]].qxy * pk[st.pk[i]].qxy + pk[st.pk[i]].qz * pk[st.pk[i]].qz;
            rhs[i] = bi;
            for (int p = 0; p < nf; ++p) {
                v[p] += a[p] * bi;
                for (int q2 = p; q2 < nf; ++q2) M[p][q2] += a[p] * a[q2];
            }
        }
        for (int p = 0; p < nf; ++p)
            for (int q2 = 0; q2 < p; ++q2) M[p][q2] = M[q2][p];
        double dmax = 0.0;
        for (int p = 0; p < nf; ++p) dmax = std::max(dmax, std::fabs(M[p][p]));
        if (!(dmax > 0.0)) return ev;
        for (int c0 = 0; c0 < nf; ++c0) {
            int p = c0;
            for (int i = c0 + 1; i < nf; ++i)
                if (std::fabs(M[i][c0]) > std::fabs(M[p][c0])) p = i;
            if (std::fabs(M[p][c0]) < 1e-11 * dmax) return ev; // rank deficient: draw no conclusion
            if (p != c0) {
                for (int j = 0; j < nf; ++j) std::swap(M[p][j], M[c0][j]);
                std::swap(v[p], v[c0]);
            }
            for (int i = c0 + 1; i < nf; ++i) {
                const double f = M[i][c0] / M[c0][c0];
                if (f == 0.0) continue;
                for (int j = c0; j < nf; ++j) M[i][j] -= f * M[c0][j];
                v[i] -= f * v[c0];
            }
        }
        for (int i = nf - 1; i >= 0; --i) {
            double s = v[i];
            for (int j = i + 1; j < nf; ++j) s -= M[i][j] * x[j];
            x[i] = s / M[i][i];
            if (!std::isfinite(x[i])) return ev;
        }
        ev.solved = true;
        double res = 0.0;
        for (int i = 0; i < m; ++i) {
            const double* a = rows[i];
            const double b = rhs[i];
            double s = 0.0;
            for (int j = 0; j < nf; ++j) s += a[j] * x[j];
            const double dq = (s - b) / (2.0 * std::sqrt(std::max(b, 1e-12)));
            res += dq * dq;
        }
        ev.rms = std::sqrt(res / (double)m);
        const std::array<double, 6> G = fillMetricArr(prm.system, x);
        ev.pd = positiveDefinite(G);
        if (!ev.pd) return ev;
        std::array<double, 6> Gi{};
        if (invMetric(G, Gi)) {
            const double s1 = st.nu[0] * (Gi[0] * st.nu[0] + Gi[3] * st.nu[1] + Gi[4] * st.nu[2]) +
                              st.nu[1] * (Gi[3] * st.nu[0] + Gi[1] * st.nu[1] + Gi[5] * st.nu[2]) +
                              st.nu[2] * (Gi[4] * st.nu[0] + Gi[5] * st.nu[1] + Gi[2] * st.nu[2]);
            ev.consist = std::fabs(s1 - 1.0);
        }
        return ev;
    };

    // Seed peak combinations: any 3 of the first nSeed peaks (the 3 smallest-|g| peaks may be linearly dependent,
    // in which case the 3x3 system is necessarily singular, so other combinations must be allowed)
    std::vector<std::array<int, 3>> seedTriples;
    for (int i = 0; i < nSeed; ++i)
        for (int j = i + 1; j < nSeed; ++j)
            for (int k = j + 1; k < nSeed; ++k) seedTriples.push_back({i, j, k});
    if (seedTriples.empty()) seedTriples.push_back({0, 1, 2});
    const int nTri = (int)seedTriples.size();

    // Expansion order: |g| descending (seeds are the smallest peaks and naturally come last)
    std::vector<int> expOrder(N);
    for (int i = 0; i < N; ++i) expOrder[i] = i;
    std::sort(expOrder.begin(), expOrder.end(), [&pk](int a, int b) {
        const double qa = pk[a].qxy * pk[a].qxy + pk[a].qz * pk[a].qz;
        const double qb = pk[b].qxy * pk[b].qxy + pk[b].qz * pk[b].qz;
        return qa > qb;
    });
    // One extra row (>minRows): used for (A) pruning and discrimination. More rows make pruning effective sooner, but each extra row
    // multiplies the state count by a branching factor; here only 1 extra row is added, the rest is left to the joint reassignment in S5.
    const int targetRows = std::min(std::min(N, kMaxRows), minRows + 1);

    // P1-1: individual parameters are already clamped by sanitize, but the product "box × shell × beam width" can still explode
    // (e.g. seedShell=64 with lim=32). Estimate the actual enumeration size; if over the limit, give up with a message.
    long long pairMax = 0;
    for (int ax = 0; ax < 3; ++ax) pairMax = std::max(pairMax, (long long)pairTab[ax].size());
    const long long seedOps = (long long)nTri * (long long)SH * (long long)S * (long long)S;
    const long long expandOps = (long long)prm.beam * (long long)targetRows * pairMax;
    const long long kSeedOpsCap = 200000000LL;   // upper bound on the number of 3x3 Cramer solves in S1
    const long long kExpandOpsCap = 50000000LL;  // upper bound on S2 expansion nodes within a single chunk
    if (seedOps > kSeedOpsCap || expandOps > kExpandOpsCap) {
        prog(0.0, "No specular: enumeration too large, reduce the index box / seed shell / beam width");
        return {};
    }

    // Candidates (G*, nu, assigned rows): candidate generation (beam search / B1 minimal-sample verification) all converge here;
    // the subsequent S5+S6 (joint reassignment + refinement + reduction) are fully shared by the three methods.
    struct Cand {
        std::array<double, 6> G;
        std::array<double, 3> nu;
        double resid = 0.0;
        std::vector<std::array<int, 3>> hkl; // assigned peaks
    };
    std::vector<Cand> cands;

    const bool dbgOn = hasEnv("GID_DEBUG");
    auto nowMs = [] { return std::chrono::duration<double, std::milli>(
                             std::chrono::steady_clock::now().time_since_epoch()).count(); };
    const int totalChunks = nTri * SH;

    // S1+S2: chunk by (seed combination, hkl of seed 1); beam search within a chunk (chunks run in parallel)
    std::mutex mu;
    std::vector<State> states;
    std::atomic<int> done{0};
    const double tS1S2 = nowMs();

    const size_t keepPerThread = std::max<size_t>(prm.beam * 2, 1024);
    std::atomic<long long> usExpand{0}, usMetric{0}, usSort{0};
    const auto byKey = [pruneRms, consistTol](const State& x, const State& y) {
        return stateKey(x, pruneRms, consistTol) < stateKey(y, pruneRms, consistTol);
    };
    // Two-stage pruning: the next layer's candidate count can reach beam×|pairTab| (tens of thousands), and solving a 6×6 normal equation for each is too expensive.
    // First truncate to keepMetric by the cheap q_z cost, then compute the (A) residual / positive definiteness / consistency for that subset.
    const size_t keepMetric = std::max<size_t>(keepPerThread, (size_t)prm.beam * 2);
    // Dynamic scheduling (task = one chunk): chunk workloads vary enormously (some chunks have thousands of states alive at
    // the end, others only a few); with static partitioning the total time is set by the slowest thread — this long tail is why
    // the wall clock barely changed across 8/12/16 threads. With a task queue every core keeps working.
    parallelForEach(totalChunks, [&](int ch) {
        {
            if (prog.stop()) return;
            const int ti = ch / SH;
            const int i0 = ch % SH;
            const auto& tr = seedTriples[ti];
            std::vector<State> cur;
            // Expansion peak sequence: the 3 S1 seed peaks are **identical for all states in this chunk**, so "the next unassigned
            // peak" is state-independent (previously every state scanned all N peaks); it is computed here outside the chunk.
            std::vector<int> expList;
            expList.reserve((size_t)N);
            for (int oi = 0; oi < N; ++oi) {
                const int i = expOrder[oi];
                if (i == tr[0] || i == tr[1] || i == tr[2]) continue;
                expList.push_back(i);
            }
            // S1: exact nu from the three seed peaks (Cramer + integer determinant prefilter)
            const double rhs[3] = {pk[tr[0]].qz, pk[tr[1]].qz, pk[tr[2]].qz};
            const int* xa = shellHalf[i0].data();
            for (int i1 = 0; i1 < S; ++i1) {
                for (int i2 = 0; i2 < S; ++i2) {
                    const int* xb = shell[i1].data();
                    const int* xc = shell[i2].data();
                    const int det = intDet3(xa, xb, xc);
                    if (det == 0) continue; // integer determinant prefilter (A1)
                    double nuV[3];
                    cramerNu(xa, xb, xc, rhs, det, nuV);
                    const double nu0 = nuV[0], nu1 = nuV[1], nu2 = nuV[2];
                    if (!std::isfinite(nu0) || !std::isfinite(nu1) || !std::isfinite(nu2))
                        continue;
                    State st;
                    st.nu[0] = nu0; st.nu[1] = nu1; st.nu[2] = nu2;
                    st.cost = 0.0;
                    st.n = 3;
                    const int* src[3] = {xa, xb, xc};
                    for (int j = 0; j < 3; ++j) {
                        st.pk[j] = (int16_t)tr[j];
                        for (int t = 0; t < 3; ++t) st.hkl[3 * j + t] = (int8_t)src[j][t];
                    }
                    cur.push_back(st);
                }
            }

            // S2: per-peak expansion (|g| descending; large-q peaks have more "independent" index rows, so (A) reaches full rank sooner).
            // next is reused across layers (swapped with cur): creating/destroying a tens-of-MB vector per layer would let the
            // first-touch page cost dominate the algorithm itself (measured: this was one of the main costs).
            std::vector<State> next;
            while (!cur.empty() && cur.front().n < targetRows) {
                const int level = cur.front().n;
                const int li = level - 3;
                if (li < 0 || li >= (int)expList.size()) break;
                const int p = expList[(size_t)li];
                next.clear();
                next.reserve(std::min<size_t>(cur.size() * 4 + 64, keepMetric));
                const auto tExp0 = std::chrono::steady_clock::now();
                for (const State& st : cur) {
                    const int cn = st.n;
                    const double qz = pk[p].qz;
                    // solving axis = largest |nu| component: xi = (qz - contribution of the other two axes)/nu_ax, rounded to the nearest integer
                    int ax = 0;
                    double amax = std::fabs(st.nu[0]);
                    for (int t = 1; t < 3; ++t) {
                        const double a = std::fabs(st.nu[t]);
                        if (a > amax) { amax = a; ax = t; }
                    }
                    if (!(amax > 1e-9)) continue;
                    const double nax = st.nu[ax];
                    const double inv = 1.0 / nax;
                    const double vmax = (double)lim[ax] + 0.5;
                    for (const auto& pr2 : pairTab[ax]) {
                        const double s =
                            pr2[0] * st.nu[0] + pr2[1] * st.nu[1] + pr2[2] * st.nu[2];
                        const double v = (qz - s) * inv;
                        if (v > vmax || v < -vmax) continue;
                        const int iv = (int)(v + (v >= 0.0 ? 0.5 : -0.5));
                        if (iv < -lim[ax] || iv > lim[ax]) continue;
                        // residual = (iv - v) * nu_ax (faster than another dot product, and unaffected by rounding)
                        const double d = (iv - v) * nax;
                        if (d > tol || d < -tol) continue;
                        std::array<int, 3> x{pr2[0], pr2[1], pr2[2]};
                        x[ax] = iv;
                        State ns = st;
                        ns.pk[cn] = (int16_t)p;
                        ns.hkl[3 * cn] = (int8_t)x[0];
                        ns.hkl[3 * cn + 1] = (int8_t)x[1];
                        ns.hkl[3 * cn + 2] = (int8_t)x[2];
                        ns.n = (int16_t)(cn + 1);
                        ns.cost = st.cost + d * d;
                        next.push_back(ns);
                    }
                }
                if (dbgOn)
                    usExpand += std::chrono::duration_cast<std::chrono::microseconds>(
                                    std::chrono::steady_clock::now() - tExp0)
                                    .count();
                const auto tMet0 = std::chrono::steady_clock::now();
                // (A) pruning (two-stage): once the row count exceeds the metric's degrees of freedom
                //   1) first truncate to keepMetric by the cheap q_z cost (this layer's candidate count can reach tens of thousands);
                //   2) then solve the metric: drop rows whose residual exceeds the threshold; add the positive-definiteness gate **only on the last row**
                //      (there S3 uses exactly these rows, so the criterion is bit-for-bit equivalent; with few rows the metric is not yet pinned and
                //      an early kill would prune the true solution branch — measured); consistency |νᵀG*⁻¹ν−1| only enters the sort key, never a hard kill.
                if (!next.empty() && next.front().n > nf && pruneRms > 0.0) {
                    if (next.size() > keepMetric) {
                        std::nth_element(next.begin(), next.begin() + keepMetric, next.end(), byKey);
                        next.resize(keepMetric);
                    }
                    // Positive definiteness / consistency are judged **only on the last row**: there S3 uses exactly these rows, so the criterion is
                    // bit-for-bit equivalent; with few rows the metric is not yet pinned and an early kill harms the true solution (measured: it prunes the whole true branch).
                    const bool lastRow = (next.front().n >= targetRows);
                    size_t w = 0;
                    for (size_t i = 0; i < next.size(); ++i) {
                        const MetricEval ev = evalMetric(next[i]);
                        if (!ev.solved) { next[i].residA = -1.0f; next[w++] = next[i]; continue; }
                        if (!(ev.rms <= pruneRms)) continue;
                        if (lastRow && !ev.pd) continue;
                        next[i].residA = (float)ev.rms;
                        next[i].residC = ev.pd ? (float)ev.consist : -1.0f;
                        next[w++] = next[i];
                    }
                    next.resize(w);
                }
                if (dbgOn)
                    usMetric += std::chrono::duration_cast<std::chrono::microseconds>(
                                    std::chrono::steady_clock::now() - tMet0)
                                    .count();
                const auto tSort0 = std::chrono::steady_clock::now();
                if (next.size() > (size_t)prm.beam) {
                    std::nth_element(next.begin(), next.begin() + prm.beam, next.end(), byKey);
                    next.resize(prm.beam);
                }
                if (dbgOn)
                    usSort += std::chrono::duration_cast<std::chrono::microseconds>(
                                  std::chrono::steady_clock::now() - tSort0)
                                  .count();
                cur.swap(next);
            }
            if (!cur.empty()) {
                // Each chunk submits only its best perChunkKeep states: the global pool is capped at 600k, and pushing all of cur would
                // truncate later chunks wholesale (the true solution may be among them). 600 chunks × 512 ≈ 300k, safely within the gate.
                const size_t perChunkKeep = 512;
                if (cur.size() > perChunkKeep) {
                    std::nth_element(cur.begin(), cur.begin() + perChunkKeep, cur.end(), byKey);
                    cur.resize(perChunkKeep);
                }
                std::lock_guard<std::mutex> lk(mu);
                if (states.size() < 600000) states.insert(states.end(), cur.begin(), cur.end());
            }
            if (prog.report && (++done) % 64 == 0)
                prog((double)done / (double)totalChunks * 0.6, "No specular: beam search progressive expansion");
        }
    }, prm.threads);

    if (prog.stop()) return {};
    if (dbgOn) {
        std::fprintf(stderr, "[nospec] S1+S2 = %.1f ms  states=%zu  (expand %.0f/metric %.0f/sort %.0f ms)\n",
                     nowMs() - tS1S2, states.size(), usExpand / 1000.0, usMetric / 1000.0,
                     usSort / 1000.0);
        size_t nGood = 0;
        for (const State& st : states)
            if (st.residA >= 0.0f) ++nGood;
        std::fprintf(stderr, "[nospec] fullRank=%zu / pool=%zu\n", nGood, states.size());
    }
    if (states.empty()) return {};
    std::sort(states.begin(), states.end(), byKey);
    if ((int)states.size() > prm.beam * 4) states.resize(prm.beam * 4);
    if (dbgOn)
        std::fprintf(stderr, "[nospec] sort+trim = %.1f ms  kept=%zu\n", nowMs() - tS1S2,
                     states.size());
    prog(0.65, "No specular: solving the reciprocal metric");

    // S3+S4: solve (A) -> G* -> unit cell -> filter (final candidate generation for the beam search / B2 paths)
    long long dbg[6] = {0, 0, 0, 0, 0, 0};
    cands.reserve(states.size());
    for (const State& st : states) {
        const int m = st.n;
        if (m < nf) { ++dbg[0]; continue; }
        std::vector<double> A((size_t)m * nf, 0.0), b(m, 0.0);
        for (int i = 0; i < m; ++i) {
            const double h = st.hkl[3 * i], k = st.hkl[3 * i + 1], l = st.hkl[3 * i + 2];
            double a[6];
            metricRow(prm.system, h, k, l, a);
            for (int j = 0; j < nf; ++j) A[(size_t)i * nf + j] = a[j];
            const int pi = st.pk[i];
            b[i] = pk[pi].qxy * pk[pi].qxy + pk[pi].qz * pk[pi].qz;
        }
        std::vector<double> x;
        if (!solveLS(A, b, m, nf, x)) { ++dbg[1]; continue; }
        // residual
        double res = 0.0;
        for (int i = 0; i < m; ++i) {
            const double h = st.hkl[3 * i], k = st.hkl[3 * i + 1], l = st.hkl[3 * i + 2];
            double a[6];
            metricRow(prm.system, h, k, l, a);
            double s = 0.0;
            for (int j = 0; j < nf; ++j) s += a[j] * x[j];
            // Convert to q units: d(q^2) ~ 2|q| dq
            const double qi = std::sqrt(std::max(b[i], 1e-12));
            const double dq = (s - b[i]) / (2.0 * qi);
            res += dq * dq;
        }
        const std::array<double, 6> G = fillMetric(prm.system, x);
        if (!positiveDefinite(G)) { ++dbg[2]; continue; }
        Cell cell;
        if (!cellFromRecipMetric(G, cell)) { ++dbg[3]; continue; }
        const double vol = cellVolume(cell);
        if (!(vol > 1.0)) { ++dbg[3]; continue; }
        if (prm.restrictCell) {
            if (cell.a < prm.aLo || cell.a > prm.aHi || cell.b < prm.bLo || cell.b > prm.bHi ||
                cell.c < prm.cLo || cell.c > prm.cHi)
                continue;
            if (vol < prm.volLo || vol > prm.volHi) continue;
        }
        // consistency nu^T G*^-1 nu = 1
        std::array<double, 6> Gi{};
        if (!invMetric(G, Gi)) { ++dbg[4]; continue; }
        const std::array<double, 3> nu{st.nu[0], st.nu[1], st.nu[2]};
        const double s1 = nu[0] * (Gi[0] * nu[0] + Gi[3] * nu[1] + Gi[4] * nu[2]) +
                          nu[1] * (Gi[3] * nu[0] + Gi[1] * nu[1] + Gi[5] * nu[2]) +
                          nu[2] * (Gi[4] * nu[0] + Gi[5] * nu[1] + Gi[2] * nu[2]);
        if (!(std::fabs(s1 - 1.0) < consistTol)) { ++dbg[5]; continue; }
        Cand c;
        c.G = G;
        c.nu = nu;
        c.resid = std::sqrt(res / (double)m);
        for (int i = 0; i < m; ++i) c.hkl.push_back({st.hkl[3 * i], st.hkl[3 * i + 1], st.hkl[3 * i + 2]});
        cands.push_back(std::move(c));
    }
    if (hasEnv("GID_DEBUG"))
        std::fprintf(stderr,
                     "[nospec] states=%zu cands=%zu  reject(rows=%lld, lsq=%lld, pd=%lld, "
                     "cell=%lld, inv=%lld, consist=%lld)\n",
                     states.size(), cands.size(), dbg[0], dbg[1], dbg[2], dbg[3], dbg[4], dbg[5]);
    if (cands.empty()) return {};
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
        return a.resid < b.resid;
    });
    if (hasEnv("GID_DEBUG")) {
        std::unordered_set<uint64_t> keys;
        int exactN = 0, uniq = 0;
        for (const Cand& c : cands) {
            if (c.resid > 1e-9) break;
            ++exactN;
            Cell cc;
            if (!cellFromRecipMetric(c.G, cc)) continue;
            const uint64_t k = cellKey(cc);
            if (!keys.insert(k).second) continue;
            ++uniq;
            if (uniq <= 20)
                std::fprintf(stderr, "[nospec] exact cand %s nu=(%.3f,%.3f,%.3f) V=%.1f\n",
                             fmtCell(cc).c_str(), c.nu[0], c.nu[1], c.nu[2], cellVolume(cc));
        }
        std::fprintf(stderr, "[nospec] exact-resid cands=%d unique=%d\n", exactN, uniq);
    }
    if ((int)cands.size() > prm.refineTop * 3) cands.resize(prm.refineTop * 3);

    // S5+S6: joint reassignment + re-solve + reduction
    const std::vector<Hkl> grid = hklGrid(prm.gridH, prm.gridK, prm.gridL);
    // The refinement freeMask is a **cell-parameter** mask (a,b,c,α,β,γ), not the same thing as a reciprocal-metric mask:
    // passing a metric mask would, for cubic/tetragonal, "free only a" while b,c stay at their initial values, a≠b=c, breaking the symmetry.
    const std::array<bool, 6> freeAll = systemFreeMask(prm.system);
    const bool useSphere = (prm.method == NoSpecMethod::Sphere);
    // |q| bandwidth for sphere decoding: the same measure as the final per-peak maximum-deviation gate
    const double sphereTolQ = prm.maxAbsMul * prm.dqxyzCutoff;
    prog(0.7, "No specular: joint reassignment and refinement");

    std::vector<Solution> out;
    std::mutex mu2;
    std::atomic<double> minRms{1e30};
    std::atomic<int> cntRms{0}, cntRed{0}, cntQ{0}, cntLi{0};
    parallelFor(0, (int)cands.size(), [&](int cb, int ce) {
        std::vector<Solution> local;
        for (int ci = cb; ci < ce; ++ci) {
            if (prog.stop()) return;
            Cand c = cands[ci];
            std::vector<std::array<int, 3>> assign(N, {0, 0, 0});
            for (int i = 0; i < (int)c.hkl.size() && i < N; ++i) assign[i] = c.hkl[i];
            double rms = 1e30;
            Cell cell;
            std::array<double, 3> nu = c.nu;
            std::array<double, 6> G = c.G;
            for (int it = 0; it < 4; ++it) {
                // Assignment: weighted (q_xy, q_z) nearest neighbour.
                // B2 (Sphere): do not use the reassignment grid box; run LLL + sphere decoding directly in the **whole integer space** —
                // (A) ellipsoid (|q| bandwidth) ∩ (B) slab (q_z tolerance); fall back to the grid nearest neighbour when a single-peak
                // enumeration exceeds the budget or the current G* is ill-conditioned, so the iteration is never interrupted.
                double tot = 0.0;
                for (int i = 0; i < N; ++i) {
                    bool done1 = false;
                    if (useSphere) {
                        std::array<int, 3> x{0, 0, 0};
                        if (sphereAssign(G, nu, pk, i, tol, sphereTolQ, x)) {
                            assign[i] = {x[0], x[1], x[2]};
                            done1 = true;
                        }
                    }
                    if (!done1) {
                        double best = 1e300;
                        int bh = 0, bk = 0, bl = 0;
                        for (const Hkl& g : grid) {
                            double qx = 0.0, qz = 0.0;
                            forwardNu(G, nu, g.h, g.k, g.l, qx, qz);
                            const double dx = qx - pk[i].qxy, dz = qz - pk[i].qz;
                            const double d2 = dx * dx + dz * dz;
                            if (d2 < best) { best = d2; bh = g.h; bk = g.k; bl = g.l; }
                        }
                        assign[i] = {bh, bk, bl};
                    }
                    double qx = 0.0, qz = 0.0;
                    forwardNu(G, nu, assign[i][0], assign[i][1], assign[i][2], qx, qz);
                    const double dx = qx - pk[i].qxy, dz = qz - pk[i].qz;
                    tot += dx * dx + dz * dz;
                }
                rms = std::sqrt(tot / (double)N);
                // re-solve (A)
                {
                    std::vector<double> A((size_t)N * nf, 0.0), b(N, 0.0);
                    for (int i = 0; i < N; ++i) {
                        const double h = assign[i][0], k = assign[i][1], l = assign[i][2];
                        double a[6];
                        metricRow(prm.system, h, k, l, a);
                        for (int j = 0; j < nf; ++j) A[(size_t)i * nf + j] = a[j];
                        b[i] = pk[i].qxy * pk[i].qxy + pk[i].qz * pk[i].qz;
                    }
                    std::vector<double> x;
                    if (!solveLS(A, b, N, nf, x)) break;
                    const std::array<double, 6> G2 = fillMetric(prm.system, x);
                    if (!positiveDefinite(G2)) break;
                    G = G2;
                }
                // re-solve (B)
                {
                    std::vector<double> A((size_t)N * 3, 0.0), b(N, 0.0);
                    for (int i = 0; i < N; ++i) {
                        A[3 * i] = assign[i][0];
                        A[3 * i + 1] = assign[i][1];
                        A[3 * i + 2] = assign[i][2];
                        b[i] = pk[i].qz;
                    }
                    std::vector<double> x;
                    if (!solveLS(A, b, N, 3, x)) break;
                    nu = {x[0], x[1], x[2]};
                }
                if (!cellFromRecipMetric(G, cell)) break;
            }
            if (!cellFromRecipMetric(G, cell)) continue;
            {
                double cur = minRms.load();
                while (rms < cur && !minRms.compare_exchange_weak(cur, rms)) {}
            }
            if (rms > prm.dqxyzCutoff * 1.5) { ++cntRms; continue; }

            // deterministic LM refinement (nu fixed)
            std::vector<Hkl> hkls;
            hkls.reserve(N);
            for (const auto& a : assign) hkls.push_back(Hkl{a[0], a[1], a[2]});
            const Orient orient = Orient::fromNu(nu);
            refineCell(pk, hkls, orient, cell, freeAll, 60);
            applySystemConstraints(prm.system, cell); // re-apply constraints such as a=b=c after refinement
            std::vector<Hkl> ah;
            assignHkl(cell, orient, pk, grid, ah);
            refineCell(pk, ah, orient, cell, freeAll, 60);
            applySystemConstraints(prm.system, cell);
            assignHkl(cell, orient, pk, grid, ah);
            hkls = ah;

            const Reduction red = reduceCell(cell);
            if (!red.ok) { ++cntRed; continue; }
            Cell rc = red.cell;
            // nu' = T^-1 nu
            std::array<int, 9> T = red.T;
            const int det = T[0] * (T[4] * T[8] - T[5] * T[7]) - T[1] * (T[3] * T[8] - T[5] * T[6]) +
                            T[2] * (T[3] * T[7] - T[4] * T[6]);
            if (std::abs(det) != 1) continue;
            // S = T^-1 = adj(T)/det
            std::array<int, 9> S{};
            S[0] = (T[4] * T[8] - T[5] * T[7]) * det;
            S[1] = (T[2] * T[7] - T[1] * T[8]) * det;
            S[2] = (T[1] * T[5] - T[2] * T[4]) * det;
            S[3] = (T[5] * T[6] - T[3] * T[8]) * det;
            S[4] = (T[0] * T[8] - T[2] * T[6]) * det;
            S[5] = (T[2] * T[3] - T[0] * T[5]) * det;
            S[6] = (T[3] * T[7] - T[4] * T[6]) * det;
            S[7] = (T[1] * T[6] - T[0] * T[7]) * det;
            S[8] = (T[0] * T[4] - T[1] * T[3]) * det;
            std::array<double, 3> nuR{};
            for (int i = 0; i < 3; ++i)
                nuR[i] = S[3 * i] * nu[0] + S[3 * i + 1] * nu[1] + S[3 * i + 2] * nu[2];
            std::vector<Hkl> rhkl;
            rhkl.reserve(hkls.size());
            for (const Hkl& x : hkls) rhkl.push_back(applyT(x, T));

            // Negating the normal and all hkl = the same solution (Friedel pair); normalize so the first non-zero component > 0,
            // otherwise the ± solutions would appear twice in the results as two records
            const int fz = (std::fabs(nuR[0]) > 1e-9) ? 0 : ((std::fabs(nuR[1]) > 1e-9) ? 1 : 2);
            if (nuR[fz] < 0.0) {
                nuR[0] = -nuR[0];
                nuR[1] = -nuR[1];
                nuR[2] = -nuR[2];
                for (Hkl& x : rhkl) { x.h = -x.h; x.k = -x.k; x.l = -x.l; }
            }

            // Contact-plane consistency: the free nu is only a continuous estimate of the "out-of-plane normal"; if only the nearest integer
            // direction were written into the result as a label, the model would still carry a lateral tilt, and the (00l) out-of-plane peaks along (u,v,w)
            // would yield a q_xy that does not exist (especially obvious for orthorhombic 001). Here nu is snapped to the exact normal of that integer
            // plane, so the reported (u,v,w) truly participates as a constraint in the error computation.
            const Mat3 R = reciprocalMatrix(rc);
            Hkl plane{0, 0, 0};
            const bool hasPlane = uvwFromNu(R, nuR, plane, 6);
            if (hasPlane) nuR = nuFromUvw(R, plane);

            const Orient ro = Orient::fromNu(nuR);
            const Quality q = evaluate(rc, ro, pk, rhkl);
            if (!(q.rmsQxyz <= prm.dqxyzCutoff)) { ++cntQ; continue; }
            // Per-peak maximum deviation: RMS alone would let through spurious solutions where individual peaks deviate far (hkl does not match q)
            if (q.maxAbs > prm.maxAbsMul * prm.dqxyzCutoff) { ++cntQ; continue; }
            // Supercell check: the observed reciprocal points must generate this reciprocal lattice (index == 1)
            if (latticeIndex(rhkl) > 1) { ++cntLi; continue; }
            if (q.rmsQxy > prm.dqxyCutoff * 3.0 || q.rmsQz > prm.dqzCutoff * 3.0) continue;
            const double vol = cellVolume(rc);
            if (prm.restrictCell) {
                if (rc.a < prm.aLo || rc.a > prm.aHi || rc.b < prm.bLo || rc.b > prm.bHi ||
                    rc.c < prm.cLo || rc.c > prm.cHi)
                    continue;
                if (vol < prm.volLo || vol > prm.volHi) continue;
            }
            // Contact plane: the integer direction was snapped above and used for error computation; when there is no integer contact plane, keep
            // the "—" semantics per P2-4 and do not invent indices.
            const Hkl uvw = hasPlane ? plane : Hkl{0, 0, 0};
            std::array<double, 6> Gi{};
            const std::array<double, 6> GR = recipMetric(rc);
            invMetric(GR, Gi);
            std::array<double, 3> nuT{Gi[0] * nuR[0] + Gi[3] * nuR[1] + Gi[4] * nuR[2],
                                      Gi[3] * nuR[0] + Gi[1] * nuR[1] + Gi[5] * nuR[2],
                                      Gi[4] * nuR[0] + Gi[5] * nuR[1] + Gi[2] * nuR[2]};
            const double nn = std::sqrt(nuT[0] * nuT[0] + nuT[1] * nuT[1] + nuT[2] * nuT[2]);
            double qspec = 0.0;
            if (hasPlane && nn > 1e-12) {
                qspec = std::sqrt((double)(uvw.h * uvw.h + uvw.k * uvw.k + uvw.l * uvw.l)) / nn;
            }
            Solution s;
            s.hasPlane = hasPlane;
            s.uvw = uvw;
            s.nu = nuR;
            s.cell = rc;
            s.vol = vol;
            s.qspec = qspec;
            s.rmsQxyz = q.rmsQxyz;
            s.rmsQxy = q.rmsQxy;
            s.rmsQz = q.rmsQz;
            s.dqspec = 0.0;
            // Expand the deduplicated hkl back to the original peak-table order (duplicates share the same indices), aligning with alignSolution
            s.hkl.resize(M);
            for (int t = 0; t < M; ++t) s.hkl[t] = rhkl[srcToUniq[t]];
            // Per-peak deviation: must use the solution's own forward model (nu), otherwise hkl does not match the measured (q_xy,q_z)
            s.dist.assign(M, 0.0);
            for (int t = 0; t < M; ++t) {
                const Peak& pt = pk[srcToUniq[t]];
                double qx = 0.0, qz = 0.0;
                forwardNu(GR, nuR, s.hkl[t].h, s.hkl[t].k, s.hkl[t].l, qx, qz);
                const double dx = qx - pt.qxy, dz = qz - pt.qz;
                s.dist[t] = std::sqrt(dx * dx + dz * dz);
            }
            local.push_back(std::move(s));
        }
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(mu2);
            out.insert(out.end(), local.begin(), local.end());
        }
    }, prm.threads);

    if (hasEnv("GID_DEBUG"))
        std::fprintf(stderr,
                     "[nospec] refined=%zu  minRms=%.5f  reject(rms=%d, reduce=%d, quality=%d, "
                     "supercell=%d)\n",
                     out.size(), minRms.load(), cntRms.load(), cntRed.load(), cntQ.load(),
                     cntLi.load());
    // Sorting: (1) residual — those below the "tolerance floor tier" are treated as equally good; (2) within a tier, smaller volume wins.
    // A larger cell with more degrees of freedom necessarily has a smaller residual; sorting purely by residual would let overfitted large cells monopolize the top.
    const double floorRms = prm.rmsFloorMul * prm.dqxyzCutoff;
    const auto scoreOf = [floorRms](const Solution& s) { return std::max(s.rmsQxyz, floorRms); };
    std::sort(out.begin(), out.end(), [&scoreOf](const Solution& a, const Solution& b) {
        const double sa = scoreOf(a), sb = scoreOf(b);
        if (sa < sb) return true;
        if (sb < sa) return false;
        return a.vol < b.vol;
    });
    std::vector<Solution> res;
    std::unordered_set<uint64_t> seen;
    for (Solution& s : out) {
        uint64_t key = cellKey(s.cell);
        key ^= (uint64_t)(s.uvw.h * 10000 + s.uvw.k * 100 + s.uvw.l) * 1099511628211ull;
        if (!seen.insert(key).second) continue;
        res.push_back(std::move(s));
        if ((int)res.size() >= prm.maxSolutions) break;
    }
    prog(1.0, "No specular: complete");
    if (hasEnv("GID_DEBUG")) {
        for (int i = 0; i < (int)res.size() && i < 5; ++i)
            std::fprintf(stderr, "[nospec] final#%d rms=%.5g %s uvw=(%d,%d,%d) V=%.2f\n", i,
                         res[i].rmsQxyz, fmtCell(res[i].cell).c_str(), res[i].uvw.h, res[i].uvw.k,
                         res[i].uvw.l, res[i].vol);
    }
    return res;
}

} // namespace gid
