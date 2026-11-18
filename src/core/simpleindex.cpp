// ============================================================================
//  simpleindex.cpp — simple indexing (out-of-plane 001) automatic search
//
//  Formulas (analyzer simple indexing, out-of-plane fixed to 001):
//      q_xy² = h²A + k²B − 2hkC
//      q_z   = hP + kQ + lR
//  The solution is carried by an ordinary reciprocal metric:
//      G* = [ A+P²,  P Q−C, P R
//             P Q−C, B+Q²,  Q R
//             P R,    Q R,   R² ]
//  Thus nu=(P,Q,R), uvw=(0,0,1), fully compatible with predictQ / UI / export.
// ============================================================================
#include "simpleindex.h"

#include "engine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace gid {
namespace {

struct X6 {
    double A = 0.0, B = 0.0, C = 0.0; // in-plane
    double P = 0.0, Q = 0.0, R = 0.0; // out-of-plane
};

struct PlaneHyp {
    X6 x;
    int inl = 0;
    double rms = 0.0;
    double cost = 0.0; // mean |h|+|k| (complexity, for MDL ordering)
};

struct Scored {
    Solution s;
    int inl = 0;
    double cost = 0.0; // mean |h|+|k|+|l|
};

// Residual cost per extra mean index unit (Å⁻¹).
// After raising the (h,k) limit, high-index / supercell planes can always push rms lower (more
// degrees of freedom), so sorting purely by rms systematically favors supercells. A linear penalty
// far below the tolerance (tol=0.05) prefers smaller, more compact index interpretations when rms is nearly equal.
constexpr double kIndexPenalty = 0.002;

double meanCost(const std::vector<Hkl>& hkl) {
    int n = 0;
    double sum = 0.0;
    for (const Hkl& x : hkl) {
        if (x.h == 0 && x.k == 0 && x.l == 0) continue;
        sum += std::abs(x.h) + std::abs(x.k) + std::abs(x.l);
        ++n;
    }
    return n ? sum / n : 0.0;
}

// Ranking key for in-plane hypotheses: inlier count first, then "rms + complexity penalty"
inline double planeScore(const PlaneHyp& p) { return p.rms + kIndexPenalty * p.cost; }

double qxyOf(const X6& x, int h, int k) {
    const double u = x.A * h * h + x.B * k * k - 2.0 * x.C * h * k;
    return std::sqrt(std::max(0.0, u));
}

double qzOf(const X6& x, int h, int k, int l) { return x.P * h + x.Q * k + x.R * l; }

// 3x3 linear system, Gauss-Jordan with partial pivoting.
bool solve3(const double m0[3][3], const double b0[3], double x[3]) {
    double a[3][4];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) a[i][j] = m0[i][j];
        a[i][3] = b0[i];
    }
    for (int col = 0; col < 3; ++col) {
        int piv = col;
        for (int r = col + 1; r < 3; ++r)
            if (std::fabs(a[r][col]) > std::fabs(a[piv][col])) piv = r;
        if (std::fabs(a[piv][col]) < 1e-12) return false;
        if (piv != col)
            for (int j = 0; j < 4; ++j) std::swap(a[piv][j], a[col][j]);
        const double d = a[col][col];
        for (int j = col; j < 4; ++j) a[col][j] /= d;
        for (int r = 0; r < 3; ++r) {
            if (r == col) continue;
            const double f = a[r][col];
            if (f == 0.0) continue;
            for (int j = col; j < 4; ++j) a[r][j] -= f * a[col][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        x[i] = a[i][3];
        if (!std::isfinite(x[i])) return false;
    }
    return true;
}

// Canonical representative of an in-plane row: row = [h², k², −2hk]; keep only one (h,k) per row.
// The row set depends only on the index box H yet is called millions of times during the search:
// cache it and precompute the three coefficients, avoiding rebuilding std::set and repeated multiplication each time.
struct Row {
    int h = 0, k = 0;
    double r0 = 0.0, r1 = 0.0, r2 = 0.0; // h², k², −2hk
};

const std::vector<Row>& rowsFor(int H) {
    static thread_local std::map<int, std::vector<Row>> cache;
    auto it = cache.find(H);
    if (it != cache.end()) return it->second;
    std::vector<Row> rows;
    std::set<std::array<int, 3>> seen;
    for (int h = -H; h <= H; ++h) {
        for (int k = -H; k <= H; ++k) {
            if (h == 0 && k == 0) continue;
            if (!seen.insert({h * h, k * k, h * k}).second) continue;
            rows.push_back(Row{h, k, (double)(h * h), (double)(k * k), -2.0 * h * k});
        }
    }
    return cache.emplace(H, std::move(rows)).first->second;
}

// Theoretical q_xy list (deduplicated by q, keeping the low-index representative)
std::vector<std::pair<double, Hkl>> theoryXy(const X6& x, const std::vector<Row>& rows) {
    std::vector<std::pair<double, Hkl>> out;
    out.reserve(rows.size());
    for (const Row& r : rows) {
        const double q = std::sqrt(std::max(0.0, x.A * r.r0 + x.B * r.r1 + x.C * r.r2));
        if (!(q > 1e-9) || !std::isfinite(q)) continue;
        bool dup = false;
        for (const auto& o : out) {
            if (std::fabs(o.first - q) < 1e-6) {
                dup = true;
                break;
            }
        }
        if (!dup) out.push_back({q, Hkl{r.h, r.k, 0}});
    }
    std::sort(out.begin(), out.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    return out;
}

bool validPlane(double A, double B, double C) {
    if (!(std::isfinite(A) && std::isfinite(B) && std::isfinite(C))) return false;
    if (A <= 1e-9 || B <= 1e-9) return false;
    // q_xy² must be non-negative for all (h,k)
    return A * B - C * C >= -1e-6 * (1.0 + A * B);
}

// Each peak independently finds the nearest q_xy; q_xy is a projected quantity, so several peaks (different l) sharing one (h,k) is normal.
int assignPlane(const std::vector<Peak>& pk, const X6& x, const SimpleParams& prm,
                const std::vector<Row>& rows, std::vector<Hkl>& hkl, double& rms) {
    const auto th = theoryXy(x, rows);
    hkl.assign(pk.size(), Hkl{});
    double ss = 0.0;
    int n = 0;
    for (size_t i = 0; i < pk.size(); ++i) {
        if (!(pk[i].qxy > prm.qMin)) continue; // rod peaks / origin do not constrain the in-plane part
        double best = 1e30;
        Hkl bh{};
        for (const auto& t : th) {
            const double d = std::fabs(t.first - pk[i].qxy);
            if (d < best) {
                best = d;
                bh = t.second;
            }
        }
        if (best <= prm.tolXy * 1.1) {
            hkl[i] = bh;
            ss += best * best;
            ++n;
        }
    }
    rms = n > 0 ? std::sqrt(ss / n) : 0.0;
    return n;
}

bool lsPlane(const std::vector<Peak>& pk, const std::vector<Hkl>& hkl, double& A, double& B,
             double& C) {
    double m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double b[3] = {0, 0, 0};
    for (size_t i = 0; i < pk.size(); ++i) {
        const Hkl& x = hkl[i];
        if (x.h == 0 && x.k == 0) continue;
        if (x.h == 0 && x.k == 0 && x.l == 0) continue;
        const double row[3] = {(double)(x.h * x.h), (double)(x.k * x.k), -2.0 * x.h * x.k};
        const double u = pk[i].qxy * pk[i].qxy;
        for (int a = 0; a < 3; ++a) {
            for (int c = 0; c < 3; ++c) m[a][c] += row[a] * row[c];
            b[a] += row[a] * u;
        }
    }
    double sol[3];
    if (!solve3(m, b, sol)) return false;
    A = sol[0];
    B = sol[1];
    C = sol[2];
    return validPlane(A, B, C);
}

bool lsOut(const std::vector<Peak>& pk, const std::vector<Hkl>& hkl, double& P, double& Q,
           double& R) {
    double m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double b[3] = {0, 0, 0};
    int n = 0;
    for (size_t i = 0; i < pk.size(); ++i) {
        const Hkl& x = hkl[i];
        if (x.h == 0 && x.k == 0 && x.l == 0) continue;
        const double row[3] = {(double)x.h, (double)x.k, (double)x.l};
        for (int a = 0; a < 3; ++a) {
            for (int c = 0; c < 3; ++c) m[a][c] += row[a] * row[c];
            b[a] += row[a] * pk[i].qz;
        }
        ++n;
    }
    if (n < 4) return false;
    double sol[3];
    if (!solve3(m, b, sol)) return false;
    P = sol[0];
    Q = sol[1];
    R = sol[2];
    return std::isfinite(P) && std::isfinite(Q) && std::isfinite(R) && R > 1e-6;
}

bool gstarPD(const X6& x) {
    const double G11 = x.A + x.P * x.P;
    const double G22 = x.B + x.Q * x.Q;
    const double G33 = x.R * x.R;
    const double G12 = x.P * x.Q - x.C;
    const double G13 = x.P * x.R;
    const double G23 = x.Q * x.R;
    if (!(G11 > 0.0 && G22 > 0.0 && G33 > 0.0)) return false;
    if (G11 * G22 - G12 * G12 <= 0.0) return false;
    const double det = G11 * (G22 * G33 - G23 * G23) - G12 * (G12 * G33 - G23 * G13) +
                       G13 * (G12 * G23 - G22 * G13);
    return det > 0.0 && std::isfinite(det);
}

bool buildCell(const X6& x, Cell& cell) {
    if (!gstarPD(x)) return false;
    const std::array<double, 6> gs{x.A + x.P * x.P, x.B + x.Q * x.Q, x.R * x.R,
                                   x.P * x.Q - x.C, x.P * x.R,        x.Q * x.R};
    if (!cellFromRecipMetric(gs, cell)) return false;
    const std::array<double, 6> g2 = recipMetric(cell);
    for (int i = 0; i < 6; ++i) {
        if (std::fabs(g2[i] - gs[i]) > 2e-4 * (1.0 + std::fabs(gs[i]))) return false;
    }
    if (!(cell.a > 0.5 && cell.b > 0.5 && cell.c > 0.5)) return false;
    if (cell.a > 1000.0 || cell.b > 1000.0 || cell.c > 1000.0) return false;
    if (!(cell.al > 0.1 && cell.al < 179.9 && cell.be > 0.1 && cell.be < 179.9 &&
          cell.ga > 0.1 && cell.ga < 179.9))
        return false;
    return true;
}

PlaneHyp refinePlane(const std::vector<Peak>& pk, X6 x, const SimpleParams& prm,
                     const std::vector<Row>& rows) {
    std::vector<Hkl> hkl;
    double rms = 0.0;
    int n = assignPlane(pk, x, prm, rows, hkl, rms);
    for (int it = 0; it < 3; ++it) {
        double A, B, C;
        if (lsPlane(pk, hkl, A, B, C)) {
            x.A = A;
            x.B = B;
            x.C = C;
        }
        const int n2 = assignPlane(pk, x, prm, rows, hkl, rms);
        if (n2 == n) break;
        n = n2;
    }
    PlaneHyp h;
    h.x = x;
    h.inl = n;
    h.rms = rms;
    h.cost = meanCost(hkl);
    return h;
}

// In-plane hypothesis generation: anchor on small-q_xy peaks, three-row exact solve + all-peak consensus + linear refinement.
//
// Hanchor is both the anchor box and the index box for consensus/refinement. This matters:
//   * With the consensus box at the full hkLimit, any hypothesis can reassign peaks to high-index rows, systematically
//     lowering rms, and the top ranks all become supercell / high-index planes (worse after raising hk);
//   * Evaluating hypotheses inside the anchor box naturally favors low-index interpretations; data that truly
//     need high indices get a deeper Hanchor from the staged search in runSimpleIndex, and refineFull finally fixes hkl at the full hkLimit.
std::vector<PlaneHyp> generatePlanes(const std::vector<Peak>& pk, const SimpleParams& prm,
                                     const Progress& prog, int Hanchor, int maxAnchors) {
    std::vector<int> sorted;
    for (size_t i = 0; i < pk.size(); ++i)
        if (pk[i].qxy > prm.qMin && std::isfinite(pk[i].qxy) && std::isfinite(pk[i].qz))
            sorted.push_back((int)i);
    std::sort(sorted.begin(), sorted.end(),
              [&](int a, int b) { return pk[a].qxy < pk[b].qxy; });

    // Anchor selection: ascending q_xy, at most 2 peaks per identical q_xy.
    //   * Low q_xy => low indices, guaranteeing completeness at small Hanchor;
    //   * With two peaks at the same q_xy, the 3x3 system forces two different rows to satisfy q(r1)=q(r2),
    //     which is itself a valid constraint (this is how aspirin's correct plane was found), so do not delete them all as redundant;
    //   * A third peak at the same q_xy only makes the triple nearly singular, so skip it.
    std::vector<int> idx;
    for (int i : sorted) {
        int same = 0;
        for (int a : idx)
            if (std::fabs(pk[a].qxy - pk[i].qxy) < 1e-6) ++same;
        if (same >= 2) continue;
        idx.push_back(i);
        if ((int)idx.size() >= maxAnchors) break;
    }
    if (idx.size() < 3) { // if peaks are too dense, fall back to unfiltered so at least one triple exists
        idx.assign(sorted.begin(), sorted.begin() + std::min<size_t>(sorted.size(), 6));
    }

    int nonRod = (int)sorted.size();
    const int minCount = std::max(3, (int)std::ceil(0.35 * nonRod));

    const std::vector<Row>& rows = rowsFor(Hanchor);
    std::map<std::array<long long, 3>, PlaneHyp> best;

    auto prune = [&]() {
        if (best.size() < 20000) return;
        std::vector<PlaneHyp> v;
        v.reserve(best.size());
        for (auto& kv : best) v.push_back(kv.second);
        std::sort(v.begin(), v.end(), [](const PlaneHyp& a, const PlaneHyp& b) {
            if (a.inl != b.inl) return a.inl > b.inl;
            return planeScore(a) < planeScore(b);
        });
        best.clear();
        for (size_t i = 0; i < v.size() && i < 6000; ++i) {
            const X6& x = v[i].x;
            best[{llround(x.A * 1000), llround(x.B * 1000), llround(x.C * 1000)}] = v[i];
        }
    };

    const int nA = (int)idx.size();
    for (int ia = 0; ia < nA; ++ia) {
        if (prog.stop()) return {};
        for (int ib = ia + 1; ib < nA; ++ib) {
            for (int ic = ib + 1; ic < nA; ++ic) {
                const int pi = idx[ia], pj = idx[ib], p1 = idx[ic];
                const double u1 = pk[pi].qxy * pk[pi].qxy;
                const double u2 = pk[pj].qxy * pk[pj].qxy;
                const double u3 = pk[p1].qxy * pk[p1].qxy;
                for (const Row& r1 : rows) {
                    for (const Row& r2 : rows) {
                        for (const Row& r3 : rows) {
                            double m[3][3] = {{r1.r0, r1.r1, r1.r2},
                                              {r2.r0, r2.r1, r2.r2},
                                              {r3.r0, r3.r1, r3.r2}};
                            const double rhs[3] = {u1, u2, u3};
                            double sol[3];
                            if (!solve3(m, rhs, sol)) continue;
                            if (!validPlane(sol[0], sol[1], sol[2])) continue;

                            X6 x;
                            x.A = sol[0];
                            x.B = sol[1];
                            x.C = sol[2];
                            // Lattice-density gate: an overly dense lattice can match almost any peak
                            const auto th = theoryXy(x, rows);
                            if (th.size() < 3) continue;
                            double minGap = 1e30;
                            for (size_t t = 1; t < th.size(); ++t)
                                minGap = std::min(minGap, th[t].first - th[t - 1].first);
                            if (minGap < prm.tolXy * 0.5) continue;

                            std::vector<Hkl> hkl;
                            double rms = 0.0;
                            const int n = assignPlane(pk, x, prm, rows, hkl, rms);
                            if (n < minCount) continue;
                            PlaneHyp h = refinePlane(pk, x, prm, rows);
                            if (h.inl < minCount) continue;
                            const auto key = std::array<long long, 3>{
                                llround(h.x.A * 1000), llround(h.x.B * 1000),
                                llround(h.x.C * 1000)};
                            auto it = best.find(key);
                            if (it == best.end() || it->second.inl < h.inl) best[key] = h;
                        }
                    }
                }
            }
        }
        prune();
    }

    std::vector<PlaneHyp> out;
    out.reserve(best.size());
    for (auto& kv : best) out.push_back(kv.second);
    std::sort(out.begin(), out.end(), [](const PlaneHyp& a, const PlaneHyp& b) {
        if (a.inl != b.inl) return a.inl > b.inl;
        return a.rms < b.rms;
    });
    // Deduplicate in parameter space so near-identical solutions of one plane do not take slots
    std::vector<PlaneHyp> uniq;
    for (const PlaneHyp& h : out) {
        bool dup = false;
        for (const PlaneHyp& o : uniq) {
            if (std::fabs(h.x.A - o.x.A) < 1e-3 * (1.0 + std::fabs(o.x.A)) &&
                std::fabs(h.x.B - o.x.B) < 1e-3 * (1.0 + std::fabs(o.x.B)) &&
                std::fabs(h.x.C - o.x.C) < 1e-3 * (1.0 + std::fabs(o.x.C))) {
                dup = true;
                break;
            }
        }
        if (!dup) uniq.push_back(h);
        if ((int)uniq.size() >= prm.maxPlanes) break;
    }
    return uniq;
}

// c* candidates: specular rod fundamental/orders + Δq_z of peak clusters sharing the same q_xy.
std::vector<double> estimateR(const std::vector<Peak>& pk, const SimpleParams& prm,
                              const RodPeriod& rp) {
    std::vector<double> raw;
    auto add = [&](double R) {
        if (!(R > 0.03 && R < 5.0) || !std::isfinite(R)) return;
        for (double v : raw)
            if (std::fabs(v - R) <= 1e-3 * (1.0 + std::fabs(R))) return;
        raw.push_back(R);
    };

    if (rp.ok && rp.p > 0.0) {
        const double base = rp.ambiguous ? (rp.qz.empty() ? rp.p : rp.qz.front()) : rp.p;
        if (base > 0.0)
            for (int n = 1; n <= prm.lLimit; ++n) add(base / n);
    }
    // Secondary: q_z differences within a same-q_xy peak cluster (shared h,k) = Δl·R
    for (size_t i = 0; i < pk.size(); ++i) {
        if (!(pk[i].qxy > prm.qMin) || !std::isfinite(pk[i].qz)) continue;
        for (size_t j = i + 1; j < pk.size(); ++j) {
            if (!(pk[j].qxy > prm.qMin) || !std::isfinite(pk[j].qz)) continue;
            if (std::fabs(pk[i].qxy - pk[j].qxy) > prm.tolXy) continue;
            const double dz = std::fabs(pk[i].qz - pk[j].qz);
            if (!(dz > prm.qMin)) continue;
            for (int n = 1; n <= prm.lLimit; ++n) add(dz / n);
        }
    }
    // Fallback: coarse scan
    for (double R = 0.15; R <= 1.55; R += 0.05) add(R);

    std::vector<double> out;
    auto push = [&](double R) {
        if (!(R > 0.03 && R < 5.0) || !std::isfinite(R)) return;
        for (double v : out)
            if (std::fabs(v - R) <= 1e-3 * (1.0 + std::fabs(R))) return;
        out.push_back(R);
    };
    // Constrain R by the rod period: q_rod must equal n·R. Returns false if the candidate contradicts the measured period.
    auto withRod = [&](double R0, double& R) {
        R = R0;
        if (!rp.ok || rp.qz.empty()) return true;
        double num = 0.0, den = 0.0;
        for (double q : rp.qz) {
            const int n = (int)std::llround(q / R0);
            if (n < 1 || std::fabs(q - n * R0) > std::max(prm.tolZ, 0.15 * R0)) return false;
            num += n * q;
            den += (double)n * n;
        }
        if (den > 0.0) R = num / den;
        for (double q : rp.qz) {
            const int n = (int)std::llround(q / R);
            if (n < 1 || std::fabs(q - n * R) > std::max(prm.tolZ, 0.15 * R)) return false;
        }
        return true;
    };

    // Rod peaks are a hard period constraint: coarse-scan candidates violating n·R are discarded outright (most of the 0.05 steps).
    for (double R0 : raw) {
        double R = 0.0;
        if (withRod(R0, R)) push(R);
    }
    // If the rod peaks contradict every candidate (e.g. a spurious q_xy≈0 peak), fall back to no discarding to avoid failing the whole path.
    if (out.empty())
        for (double R0 : raw) push(R0);
    if (out.size() > 16) out.resize(16);
    return out;
}

std::vector<std::vector<Hkl>> buildCands(const std::vector<Peak>& pk, const X6& x,
                                         const SimpleParams& prm) {
    std::vector<std::vector<Hkl>> out(pk.size());
    const int H = prm.hkLimit;
    for (size_t i = 0; i < pk.size(); ++i) {
        if (!std::isfinite(pk[i].qxy) || !std::isfinite(pk[i].qz)) continue;
        for (int h = -H; h <= H; ++h) {
            for (int k = -H; k <= H; ++k) {
                const double q = qxyOf(x, h, k);
                if (std::fabs(q - pk[i].qxy) <= prm.tolXy * 1.25) out[i].push_back(Hkl{h, k, 0});
            }
        }
        std::sort(out[i].begin(), out[i].end(), [&](const Hkl& a, const Hkl& b) {
            const double da = std::fabs(qxyOf(x, a.h, a.k) - pk[i].qxy);
            const double db = std::fabs(qxyOf(x, b.h, b.k) - pk[i].qxy);
            if (std::fabs(da - db) > 1e-9) return da < db;
            return std::abs(a.h) + std::abs(a.k) < std::abs(b.h) + std::abs(b.k);
        });
        if (out[i].size() > 16) out[i].resize(16);
    }
    return out;
}

double tiltScore(double p, double q, const std::vector<Peak>& pk,
                 const std::vector<std::vector<Hkl>>& cands, double R, double tolZ) {
    double s = 0.0;
    const double inv = 1.0 / R;
    for (size_t i = 0; i < pk.size(); ++i) {
        if (cands[i].empty()) continue;
        const double zr = pk[i].qz * inv;
        double best = -1.0;
        for (const Hkl& x : cands[i]) {
            const double d = zr - (p * x.h + q * x.k);
            const double fr = d - std::round(d);
            const double t = fr * R / tolZ;
            const double rho = 1.0 - t * t;
            if (rho > best) best = rho;
        }
        if (best > 0.0) s += best;
    }
    return s;
}

double wrapHalf(double v) {
    while (v >= 0.5) v -= 1.0;
    while (v < -0.5) v += 1.0;
    return v;
}

std::vector<std::pair<double, double>> houghTilt(const std::vector<Peak>& pk,
                                                 const std::vector<std::vector<Hkl>>& cands,
                                                 double R, const SimpleParams& prm) {
    const int G = std::max(9, prm.grid);
    const double step = 1.0 / G;
    struct Node {
        double s, p, q;
    };
    std::vector<Node> grid;
    grid.reserve(G * G);
    for (int gi = 0; gi < G; ++gi) {
        const double p = -0.5 + gi * step;
        for (int gj = 0; gj < G; ++gj) {
            const double q = -0.5 + gj * step;
            grid.push_back(Node{tiltScore(p, q, pk, cands, R, prm.tolZ), p, q});
        }
    }
    std::sort(grid.begin(), grid.end(), [](const Node& a, const Node& b) { return a.s > b.s; });

    std::vector<std::pair<double, double>> seeds;
    const int want = std::max(1, prm.maxTilts) * 2;
    for (const Node& n : grid) {
        bool close = false;
        for (const auto& s : seeds) {
            double dp = std::fabs(n.p - s.first);
            double dq = std::fabs(n.q - s.second);
            dp = std::min(dp, 1.0 - dp);
            dq = std::min(dq, 1.0 - dq);
            if (std::max(dp, dq) < 2.0 * step) {
                close = true;
                break;
            }
        }
        if (close) continue;
        seeds.push_back({n.p, n.q});
        if ((int)seeds.size() >= want) break;
    }

    std::vector<std::pair<double, double>> out;
    for (const auto& s0 : seeds) {
        double p = s0.first, q = s0.second;
        double best = tiltScore(p, q, pk, cands, R, prm.tolZ);
        double st = step;
        for (int it = 0; it < 16; ++it) {
            bool moved = false;
            const double dp[4] = {st, -st, 0.0, 0.0};
            const double dq[4] = {0.0, 0.0, st, -st};
            for (int k = 0; k < 4; ++k) {
                const double np = wrapHalf(p + dp[k]);
                const double nq = wrapHalf(q + dq[k]);
                const double sc = tiltScore(np, nq, pk, cands, R, prm.tolZ);
                if (sc > best + 1e-12) {
                    best = sc;
                    p = np;
                    q = nq;
                    moved = true;
                }
            }
            if (!moved) st *= 0.5;
        }
        bool dup = false;
        for (const auto& o : out) {
            double dp = std::fabs(wrapHalf(p - o.first));
            double dq = std::fabs(wrapHalf(q - o.second));
            if (dp < 1e-3 && dq < 1e-3) {
                dup = true;
                break;
            }
        }
        if (!dup) out.push_back({p, q});
    }
    if ((int)out.size() > prm.maxTilts) out.resize(prm.maxTilts);
    return out;
}

void assignTilt(const std::vector<Peak>& pk, const std::vector<std::vector<Hkl>>& cands, const X6& x,
                const SimpleParams& prm, std::vector<Hkl>& hkl) {
    hkl.assign(pk.size(), Hkl{});
    for (size_t i = 0; i < pk.size(); ++i) {
        if (cands[i].empty()) continue;
        double best = 1e30;
        Hkl bh{};
        for (const Hkl& hk : cands[i]) {
            const double lf = (pk[i].qz - x.P * hk.h - x.Q * hk.k) / x.R;
            const int l = (int)std::llround(lf);
            if (std::abs(l) > prm.lLimit) continue;
            const double d = std::fabs(lf - l) * x.R;
            if (d < best) {
                best = d;
                bh = Hkl{hk.h, hk.k, l};
            }
        }
        if (best <= prm.tolZ * 1.25) hkl[i] = bh;
    }
}

struct TheoryPoint {
    double qxy, qz;
    Hkl hkl;
    int cost;
};

std::vector<TheoryPoint> buildTheory(const X6& x, const SimpleParams& prm) {
    std::vector<TheoryPoint> out;
    std::set<std::array<long long, 2>> seen;
    const int H = prm.hkLimit, L = prm.lLimit;
    for (int h = -H; h <= H; ++h) {
        for (int k = -H; k <= H; ++k) {
            for (int l = -L; l <= L; ++l) {
                if (h == 0 && k == 0 && l == 0) continue;
                const double qx = qxyOf(x, h, k);
                const double qz = qzOf(x, h, k, l);
                if (!std::isfinite(qx) || !std::isfinite(qz)) continue;
                const long long kx = llround(qx * 1e5);
                const long long kz = llround(qz * 1e5);
                if (!seen.insert({kx, kz}).second) continue;
                const int cost = std::abs(h) + std::abs(k) + std::abs(l);
                out.push_back(TheoryPoint{qx, qz, Hkl{h, k, l}, cost});
            }
        }
    }
    return out;
}

void assignFull(const std::vector<Peak>& pk, const X6& x, const SimpleParams& prm, double tol,
                std::vector<Hkl>& hkl, int& inl, double& rmsXy, double& rmsZ, double& rms3,
                int* maxIdx = nullptr, bool wantDist = false, std::vector<double>* distOut = nullptr) {
    const auto th = buildTheory(x, prm);
    hkl.assign(pk.size(), Hkl{});
    if (wantDist && distOut) distOut->assign(pk.size(), std::numeric_limits<double>::quiet_NaN());
    double sx = 0.0, sz = 0.0, s3 = 0.0;
    int n = 0, mi = 0;
    for (size_t i = 0; i < pk.size(); ++i) {
        if (!std::isfinite(pk[i].qxy) || !std::isfinite(pk[i].qz)) continue;
        double best = 1e30;
        const TheoryPoint* bp = nullptr;
        for (const auto& t : th) {
            const double dx = t.qxy - pk[i].qxy;
            const double dz = t.qz - pk[i].qz;
            const double d2 = dx * dx + dz * dz;
            if (d2 < best) {
                best = d2;
                bp = &t;
            }
        }
        if (!bp) continue;
        const double d = std::sqrt(best);
        if (d <= tol) {
            hkl[i] = bp->hkl;
            const double dx = bp->qxy - pk[i].qxy;
            const double dz = bp->qz - pk[i].qz;
            sx += dx * dx;
            sz += dz * dz;
            s3 += dx * dx + dz * dz;
            mi = std::max(mi, bp->cost);
            ++n;
            if (wantDist && distOut) (*distOut)[i] = d;
        }
    }
    inl = n;
    rmsXy = n > 0 ? std::sqrt(sx / n) : 0.0;
    rmsZ = n > 0 ? std::sqrt(sz / n) : 0.0;
    rms3 = n > 0 ? std::sqrt(s3 / n) : 0.0;
    if (maxIdx) *maxIdx = mi;
}

bool refineFull(const std::vector<Peak>& pk, const SimpleParams& prm, X6& x, std::vector<Hkl>& hkl,
                int& inl, double& rmsXy, double& rmsZ, double& rms3) {
    std::vector<Hkl> cur = hkl;
    X6 bx = x;
    std::vector<Hkl> bestHkl = hkl;
    int bestInl = -1;
    double best3 = 1e30, bestXy = 0.0, bestZ = 0.0;

    for (int it = 0; it <= prm.icpIters; ++it) {
        if (it > 0) {
            double A, B, C;
            if (lsPlane(pk, cur, A, B, C)) {
                X6 trial = x;
                trial.A = A;
                trial.B = B;
                trial.C = C;
                if (gstarPD(trial)) x = trial;
            }
            double P, Q, R;
            if (lsOut(pk, cur, P, Q, R)) {
                X6 trial = x;
                trial.P = P;
                trial.Q = Q;
                trial.R = R;
                if (gstarPD(trial)) x = trial;
            }
            assignFull(pk, x, prm, prm.tol * 1.25, cur, inl, rmsXy, rmsZ, rms3);
        } else {
            assignFull(pk, x, prm, prm.tol * 1.25, cur, inl, rmsXy, rmsZ, rms3);
        }
        if (inl > bestInl || (inl == bestInl && rms3 < best3)) {
            bestInl = inl;
            best3 = rms3;
            bestXy = rmsXy;
            bestZ = rmsZ;
            bestHkl = cur;
            bx = x;
        }
        if (it > 0) {
            bool same = cur.size() == hkl.size();
            if (same)
                for (size_t i = 0; i < cur.size(); ++i)
                    if (cur[i].h != bestHkl[i].h || cur[i].k != bestHkl[i].k ||
                        cur[i].l != bestHkl[i].l) {
                        same = false;
                        break;
                    }
            if (same) break;
        }
        hkl = cur;
    }
    x = bx;
    hkl = bestHkl;
    inl = bestInl;
    rmsXy = bestXy;
    rmsZ = bestZ;
    rms3 = best3;
    return inl >= 3;
}

bool makeSolution(const std::vector<Peak>& pk, const X6& x, const SimpleParams& prm,
                  const RodPeriod& rp, Solution& s, int& inl, int& maxIdx) {
    Cell cell;
    if (!buildCell(x, cell)) return false;
    std::vector<Hkl> hkl;
    std::vector<double> dist;
    double rmsXy = 0.0, rmsZ = 0.0, rms3 = 0.0;
    assignFull(pk, x, prm, prm.tol, hkl, inl, rmsXy, rmsZ, rms3, &maxIdx, true, &dist);
    if (inl < 3) return false;

    s = Solution{};
    s.uvw = Hkl{0, 0, 1};
    s.hasPlane = true;
    s.nu = {x.P, x.Q, x.R};
    s.cell = cell;
    s.vol = cellVolume(cell);
    s.qspec = x.R;
    s.dqspec = rp.ok ? std::fabs(x.R - rp.p) : 0.0;
    s.rmsQxy = rmsXy;
    s.rmsQz = rmsZ;
    s.rmsQxyz = rms3;
    s.hkl = std::move(hkl);
    s.dist = std::move(dist);
    return true;
}

// Whether two solutions predict the measured peaks identically (= the data cannot distinguish them). Only peaks
// assigned in both are compared; if one assigns a peak and the other does not, inlier coverage differs and they are not equivalent.
bool samePrediction(const Solution& a, const Solution& b) {
    int common = 0;
    const size_t n = std::min(a.hkl.size(), b.hkl.size());
    for (size_t i = 0; i < n; ++i) {
        const Hkl& x = a.hkl[i];
        const Hkl& y = b.hkl[i];
        const bool xa = !(x.h == 0 && x.k == 0 && x.l == 0);
        const bool ya = !(y.h == 0 && y.k == 0 && y.l == 0);
        if (!xa && !ya) continue;
        if (xa != ya) return false;
        const QPoint p = predictQ(a, x.h, x.k, x.l);
        const QPoint q = predictQ(b, y.h, y.k, y.l);
        if (std::fabs(p.qxy - q.qxy) > 1e-4 || std::fabs(p.qz - q.qz) > 1e-4) return false;
        ++common;
    }
    return common > 0;
}

} // namespace

std::vector<Solution> runSimpleIndex(const std::vector<Peak>& peaksIn, const SimpleParams& prm,
                                     const Progress& prog) {
    std::vector<Peak> pk = peaksIn;
    normalizeQxySigns(pk);

    int nonRod = 0;
    for (const Peak& p : pk)
        if (p.qxy > prm.qMin && std::isfinite(p.qxy) && std::isfinite(p.qz)) ++nonRod;
    if (nonRod < 3) return {};

    // Staged in-plane search: start the anchor/consensus box at 2 and stop when the result is good
    // enough; deepen step by step to hkLimit only if needed. The candidate count is proportional to
    // rows(H)³ (~1.7k at H=2, ~27k at H=4, ~64k at H=6), so deepening on demand is one to two orders of magnitude faster than using hkLimit directly; most data satisfy inl >= 90% of non-rod peaks at H=2.
    // Deepening only adds candidates and never degrades quality; the final hkl is reassigned by refineFull at the full hkLimit.
    prog(0.02, "Simple indexing: in-plane hypotheses");
    std::vector<PlaneHyp> planes;
    const int needInl = (int)std::ceil(0.9 * nonRod);
    for (int ha = std::min(2, prm.hkLimit); ha <= prm.hkLimit; ++ha) {
        if (prog.stop()) break;
        prog(std::min(0.30, 0.02 + 0.05 * (ha - 2)),
             ha == 2 ? "Simple indexing: in-plane hypotheses" : "Simple indexing: deepening in-plane search box");
        std::vector<PlaneHyp> more = generatePlanes(pk, prm, prog, ha, ha <= 3 ? 10 : 12);
        planes.insert(planes.end(), more.begin(), more.end());
        std::sort(planes.begin(), planes.end(), [](const PlaneHyp& a, const PlaneHyp& b) {
            if (a.inl != b.inl) return a.inl > b.inl;
            return planeScore(a) < planeScore(b);
        });
        // Different anchor boxes may find the same plane: deduplicate in parameter space
        std::vector<PlaneHyp> uniq;
        for (const PlaneHyp& h : planes) {
            bool dup = false;
            for (const PlaneHyp& o : uniq)
                if (std::fabs(h.x.A - o.x.A) < 1e-3 * (1.0 + std::fabs(o.x.A)) &&
                    std::fabs(h.x.B - o.x.B) < 1e-3 * (1.0 + std::fabs(o.x.B)) &&
                    std::fabs(h.x.C - o.x.C) < 1e-3 * (1.0 + std::fabs(o.x.C))) {
                    dup = true;
                    break;
                }
            if (!dup) uniq.push_back(h);
        }
        planes.swap(uniq);
        if (prm.maxPlanes > 0 && (int)planes.size() > prm.maxPlanes) planes.resize(prm.maxPlanes);
        if (!planes.empty() && planes.front().inl >= needInl) break;
    }
    if (planes.empty()) return {};

    prog(0.35, "Simple indexing: out-of-plane period and tilt");
    RodPeriod rp = detectRodPeriod(pk, prm.rodTol, prm.lLimit);
    std::vector<double> Rlist = estimateR(pk, prm, rp);
    if (Rlist.empty()) return {};

    std::vector<Scored> out;
    const int total = (int)planes.size() * (int)Rlist.size();
    int done = 0;
    for (const PlaneHyp& ph : planes) {
        if (prog.stop()) break;
        const auto cands = buildCands(pk, ph.x, prm);
        for (double R : Rlist) {
            if (prog.stop()) break;
            const auto tilts = houghTilt(pk, cands, R, prm);
            for (const auto& tq : tilts) {
                X6 x = ph.x;
                x.R = R;
                x.P = tq.first * R;
                x.Q = tq.second * R;
                std::vector<Hkl> hkl;
                assignTilt(pk, cands, x, prm, hkl);
                int inl = 0, maxIdx = 0;
                double rxy = 0.0, rz = 0.0, r3 = 0.0;
                if (!refineFull(pk, prm, x, hkl, inl, rxy, rz, r3)) continue;
                Solution s;
                if (!makeSolution(pk, x, prm, rp, s, inl, maxIdx)) continue;
                out.push_back(Scored{s, inl, meanCost(s.hkl)});
            }
            ++done;
            if (total > 0)
                prog(0.35 + 0.6 * (double)done / (double)total, "Simple indexing: refining candidates");
        }
    }

    std::sort(out.begin(), out.end(), [](const Scored& a, const Scored& b) {
        if (a.inl != b.inl) return a.inl > b.inl;
        // The complexity penalty ranks supercells after the primitive cell under equal fit (especially important when raising the hk limit)
        return a.s.rmsQxyz + kIndexPenalty * a.cost < b.s.rmsQxyz + kIndexPenalty * b.cost;
    });

    std::vector<Solution> res;
    for (const Scored& c : out) {
        bool dup = false;
        for (const Solution& o : res) {
            const double ca[6] = {c.s.cell.a, c.s.cell.b, c.s.cell.c,
                                  c.s.cell.al, c.s.cell.be, c.s.cell.ga};
            const double oa[6] = {o.cell.a, o.cell.b, o.cell.c,
                                  o.cell.al, o.cell.be, o.cell.ga};
            bool same = true;
            for (int i = 0; i < 6; ++i) {
                if (std::fabs(ca[i] - oa[i]) > 1e-3 * (1.0 + std::fabs(oa[i]))) {
                    same = false;
                    break;
                }
            }
            if (same)
                for (int i = 0; i < 3; ++i)
                    if (std::fabs(c.s.nu[i] - o.nu[i]) > 1e-3) {
                        same = false;
                        break;
                    }
            // Equivalent solutions indistinguishable at the data level (c doubled, global sign flip, ...): if
            // per-peak predictions agree they are duplicates; keep only the higher-ranked branch (smaller rms or smaller indices).
            if (!same && samePrediction(c.s, o)) same = true;
            if (same) {
                dup = true;
                break;
            }
        }
        if (!dup) res.push_back(c.s);
        if ((int)res.size() >= prm.topK) break;
    }
    prog(0.99, "Simple indexing: done");
    return res;
}

} // namespace gid
