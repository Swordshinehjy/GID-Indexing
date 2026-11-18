// ============================================================================
//  refine.cpp — see refine.h
// ============================================================================
#include "refine.h"

#include <algorithm>
#include <cmath>

namespace gid {
namespace {

std::array<double, 6> recipMetricOf(const Mat3& r) {
    std::array<double, 6> g{};
    g[0] = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    g[1] = r[3] * r[3] + r[4] * r[4] + r[5] * r[5];
    g[2] = r[6] * r[6] + r[7] * r[7] + r[8] * r[8];
    g[3] = r[0] * r[3] + r[1] * r[4] + r[2] * r[5];
    g[4] = r[0] * r[6] + r[1] * r[7] + r[2] * r[8];
    g[5] = r[3] * r[6] + r[4] * r[7] + r[5] * r[8];
    return g;
}

struct Theo {
    double qxy = 0.0, qz = 0.0;
    int h = 0, k = 0, l = 0;
};

std::vector<Theo> theoTable(const Cell& cell, const Orient& orient, const std::vector<Hkl>& cands,
                           double maxG) {
    const std::array<double, 6> G = recipMetricOf(reciprocalMatrix(cell));
    const std::array<double, 3> nu = orient.get(G);
    const double lim2 = (maxG < 1e29) ? maxG * maxG : 1e300;
    std::vector<Theo> pts;
    pts.reserve(cands.size());
    for (const Hkl& c : cands) {
        // |G|^2 = (h,k,l) G* (h,k,l)^T; indices exceeding maxG cannot be the nearest neighbour of any peak
        const double g2 = G[0] * c.h * c.h + G[1] * c.k * c.k + G[2] * c.l * c.l +
                          2.0 * (G[3] * c.h * c.k + G[4] * c.h * c.l + G[5] * c.k * c.l);
        if (g2 > lim2) continue;
        Theo t;
        forwardNu(G, nu, c.h, c.k, c.l, t.qxy, t.qz);
        t.h = c.h; t.k = c.k; t.l = c.l;
        pts.push_back(t);
    }
    // No more sorting/deduplication: coincident theoretical points (qxy,qz) are exactly identical and do not affect the nearest-neighbour distance;
    // ties in the nearest-neighbour query are ordered by (h,k,l) lexicographic order, making the result deterministic.
    return pts;
}

// Uniform grid on the (qxy,qz) plane: exact nearest-neighbour query, reducing the O(peaks×points) linear scan
// to "start cell + ring expansion", stopping when a ring's geometric lower bound exceeds the current best (equivalent to linear scan).
class NearestGrid {
  public:
    NearestGrid(const std::vector<Theo>& pts, double wz) : pts_(pts) {
        const int n = (int)pts.size();
        double x0 = 1e300, x1 = -1e300, z0 = 1e300, z1 = -1e300;
        for (const Theo& p : pts) {
            x0 = std::min(x0, p.qxy);
            x1 = std::max(x1, p.qxy);
            z0 = std::min(z0, p.qz);
            z1 = std::max(z1, p.qz);
        }
        x0_ = x0;
        z0_ = z0;
        const double sx = std::max(x1 - x0, 1e-9);
        const double swz = std::max(std::fabs(wz) * (z1 - z0), 1e-9);
        // Goal: total cell count ≈ n/2 (two points per cell on average), distributed by weighted span, approximately square
        const double cells = std::max(4.0, 0.5 * (double)n);
        const double tot = sx + swz;
        nx_ = (int)std::lround(std::sqrt(cells * sx / tot * 2.0));
        nz_ = (int)std::lround(std::sqrt(cells * swz / tot * 2.0));
        nx_ = std::max(1, std::min(nx_, 128));
        nz_ = std::max(1, std::min(nz_, 128));
        cx_ = sx / (double)nx_;
        cz_ = std::max(z1 - z0, 0.0) / (double)nz_;
        head_.assign((size_t)nx_ * nz_, -1);
        next_.assign((size_t)n, -1);
        for (int i = 0; i < n; ++i) {
            const size_t c = (size_t)cellZ(pts[i].qz) * nx_ + (size_t)cellX(pts[i].qxy);
            next_[i] = head_[c];
            head_[c] = i;
        }
        // Geometric lower bound for ring expansion: only use axes with "more than 1 cell" (an axis with 1 cell provides no bound)
        const double bx = (nx_ > 1) ? cx_ : 1e300;
        const double bz = (nz_ > 1) ? std::fabs(wz) * cz_ : 1e300;
        minCell_ = std::min(bx, bz);
    }

    // Return the nearest point index; bestD2 is the squared weighted distance. Equivalent to linear scan; when d² is exactly equal, take the lexicographically smallest (h,k,l).
    int nearest(double px, double pz, double wz, double& bestD2) {
        int bi = -1;
        double bd = 1e300;
        const int ix = cellX(px), iz = cellZ(pz);
        const int maxRing = std::max(nx_, nz_);
        for (int r = 0; r <= maxRing; ++r) {
            if (r > 0) {
                const double lb = (double)(r - 1) * minCell_;
                if (bd <= lb * lb) break;
            }
            if (r == 0) {
                scanCell((size_t)iz * nx_ + (size_t)ix, px, pz, wz, bi, bd);
                continue;
            }
            for (int dz = -r; dz <= r; ++dz) {
                const int zz = iz + dz;
                if (zz < 0 || zz >= nz_) continue;
                const int stepx = (dz == -r || dz == r) ? 1 : 2 * r;
                for (int dx = -r; dx <= r; dx += stepx) {
                    const int xx = ix + dx;
                    if (xx < 0 || xx >= nx_) continue;
                    scanCell((size_t)zz * nx_ + (size_t)xx, px, pz, wz, bi, bd);
                }
            }
        }
        bestD2 = bd;
        return bi;
    }

  private:
    int cellX(double x) const {
        if (!(x > x0_)) return 0;
        const int i = (int)((x - x0_) / cx_);
        return i >= nx_ ? nx_ - 1 : i;
    }
    int cellZ(double z) const {
        if (!(z > z0_)) return 0;
        const int i = (int)((z - z0_) / cz_);
        return i >= nz_ ? nz_ - 1 : i;
    }
    void scanCell(size_t c, double px, double pz, double wz, int& bi, double& bd) const {
        for (int idx = head_[c]; idx >= 0; idx = next_[idx]) {
            const Theo& p = pts_[idx];
            const double dx = px - p.qxy;
            const double dz = wz * (pz - p.qz);
            const double d2 = dx * dx + dz * dz;
            if (d2 < bd || (d2 == bd && (bi < 0 || hklLess(p, pts_[bi]))) ) {
                bd = d2;
                bi = idx;
            }
        }
    }
    static bool hklLess(const Theo& a, const Theo& b) {
        if (a.h != b.h) return a.h < b.h;
        if (a.k != b.k) return a.k < b.k;
        return a.l < b.l;
    }

    const std::vector<Theo>& pts_;
    double x0_ = 0.0, z0_ = 0.0, cx_ = 1.0, cz_ = 1.0, minCell_ = 1e300;
    int nx_ = 1, nz_ = 1;
    std::vector<int> head_, next_;
};

// Small symmetric positive-definite system (Gauss + partial pivoting), n <= 9
bool solveSPD(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(A[(size_t)i * n + k]) > std::fabs(A[(size_t)p * n + k])) p = i;
        if (std::fabs(A[(size_t)p * n + k]) < 1e-14) return false;
        if (p != k) {
            for (int j = 0; j < n; ++j) std::swap(A[(size_t)p * n + j], A[(size_t)k * n + j]);
            std::swap(b[p], b[k]);
        }
        for (int i = k + 1; i < n; ++i) {
            const double f = A[(size_t)i * n + k] / A[(size_t)k * n + k];
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) A[(size_t)i * n + j] -= f * A[(size_t)k * n + j];
            b[i] -= f * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[(size_t)i * n + j] * b[j];
        b[i] = s / A[(size_t)i * n + i];
    }
    return true;
}

void cellToVec(const Cell& c, double* x) {
    x[0] = c.a; x[1] = c.b; x[2] = c.c; x[3] = c.al; x[4] = c.be; x[5] = c.ga;
}
Cell vecToCell(const double* x) { return Cell{x[0], x[1], x[2], x[3], x[4], x[5]}; }

} // namespace

std::array<double, 3> Orient::get(const std::array<double, 6>& G) const {
    if (useNu) return nu;
    const double n[3] = {(double)uvw.h, (double)uvw.k, (double)uvw.l};
    const double q2 = G[0] * n[0] * n[0] + G[1] * n[1] * n[1] + G[2] * n[2] * n[2] +
                      2.0 * (G[3] * n[0] * n[1] + G[4] * n[0] * n[2] + G[5] * n[1] * n[2]);
    if (!(q2 > 1e-24)) return {0, 0, 0};
    const double qs = std::sqrt(q2);
    return {(G[0] * n[0] + G[3] * n[1] + G[4] * n[2]) / qs,
            (G[3] * n[0] + G[1] * n[1] + G[5] * n[2]) / qs,
            (G[4] * n[0] + G[5] * n[1] + G[2] * n[2]) / qs};
}

void forwardNu(const std::array<double, 6>& G, const std::array<double, 3>& nu, int h, int k,
               int l, double& qxy, double& qz) {
    const double g2 = G[0] * h * h + G[1] * k * k + G[2] * l * l +
                      2.0 * (G[3] * h * k + G[4] * h * l + G[5] * k * l);
    const double z = nu[0] * h + nu[1] * k + nu[2] * l;
    qz = z;
    qxy = std::sqrt(std::max(0.0, g2 - z * z));
}

std::vector<Hkl> hklGrid(int maxH, int maxK, int maxL) {
    std::vector<Hkl> out;
    out.reserve((size_t)(2 * maxH + 1) * (2 * maxK + 1) * (2 * maxL + 1));
    for (int h = -maxH; h <= maxH; ++h)
        for (int k = -maxK; k <= maxK; ++k)
            for (int l = -maxL; l <= maxL; ++l) {
                if (h == 0 && k == 0 && l == 0) continue;
                out.push_back(Hkl{h, k, l});
            }
    return out;
}

double assignHkl(const Cell& cell, const Orient& orient, const std::vector<Peak>& peaks,
                 const std::vector<Hkl>& cands, std::vector<Hkl>& outHkl,
                 std::vector<double>* outDist, double cutoff, double weightQz, double maxG) {
    outHkl.assign(peaks.size(), Hkl{0, 0, 0});
    if (outDist) outDist->assign(peaks.size(), 0.0);
    if (peaks.empty()) return 0.0;
    const std::vector<Theo> pts = theoTable(cell, orient, cands, maxG);
    if (pts.empty()) return 1e30;

    const double limit = cutoff < 1e29 ? cutoff * cutoff * (double)peaks.size() : 1e300;
    double total = 0.0;
    NearestGrid grid(pts, weightQz);
    for (size_t i = 0; i < peaks.size(); ++i) {
        double best = 1e300;
        const int bi = grid.nearest(peaks[i].qxy, peaks[i].qz, weightQz, best);
        if (bi < 0) return 1e30;
        outHkl[i] = Hkl{pts[bi].h, pts[bi].k, pts[bi].l};
        if (outDist) (*outDist)[i] = std::sqrt(best);
        total += best;
        if (total > limit) return std::sqrt(total / (double)peaks.size());
    }
    return std::sqrt(total / (double)peaks.size());
}

bool refineCell(const std::vector<Peak>& peaks, const std::vector<Hkl>& hkl, const Orient& orient,
                Cell& cell, const std::array<bool, 6>& freeMask, int maxIter,
                CrystalSystem sys) {
    const size_t n = peaks.size();
    if (n == 0 || hkl.size() != n) return false;
    int idx[6];
    int nf = 0;
    for (int i = 0; i < 6; ++i)
        if (freeMask[i]) idx[nf++] = i;
    if (nf == 0) return true;

    double x[6];
    cellToVec(cell, x);
    // Crystal-system-constrained "reparameterisation": one free variable drives a group of dependent cell parameters, so the residual model matches the constraints.
    // The old implementation only froze dependent parameters and projected after the call; LM was actually solving an inconsistent model
    // (e.g. in tetragonal, a changes while b does not), which systematically biased the step direction.
    // Triclinic goes through vecToCell, bit-for-bit equivalent to the original implementation.
    auto expand = [&](const double* p) -> Cell {
        switch (sys) {
            case CrystalSystem::Cubic:
                return Cell{p[0], p[0], p[0], 90.0, 90.0, 90.0};
            case CrystalSystem::Tetragonal:
                return Cell{p[0], p[0], p[2], 90.0, 90.0, 90.0};
            case CrystalSystem::Hexagonal:
            case CrystalSystem::Trigonal:
                return Cell{p[0], p[0], p[2], 90.0, 90.0, 120.0};
            case CrystalSystem::Orthorhombic:
                return Cell{p[0], p[1], p[2], 90.0, 90.0, 90.0};
            case CrystalSystem::Monoclinic:
                return Cell{p[0], p[1], p[2], 90.0, p[4], 90.0};
            default:
                return vecToCell(p);
        }
    };
    const double lo[6] = {2.0, 2.0, 2.0, 40.0, 40.0, 40.0};
    const double hi[6] = {60.0, 60.0, 60.0, 140.0, 140.0, 140.0};

    auto residuals = [&](const double* p, std::vector<double>& r) {
        const Cell cc = expand(p);
        const std::array<double, 6> G = recipMetricOf(reciprocalMatrix(cc));
        const std::array<double, 3> nu = orient.get(G);
        for (size_t i = 0; i < n; ++i) {
            double qx = 0.0, qz = 0.0;
            forwardNu(G, nu, hkl[i].h, hkl[i].k, hkl[i].l, qx, qz);
            r[2 * i] = qx - peaks[i].qxy;
            r[2 * i + 1] = qz - peaks[i].qz;
        }
    };
    // All workspaces preallocated: this function runs 2 rounds on every raw solution of Part TWO, each up to 60
    // iterations × up to 10 trials; the old style created 5~7 vectors per round (up to hundreds of heap
    // allocations per refinement), one of the main costs across thousands of refinement calls.
    std::vector<double> r(2 * n), rt(2 * n), rc(2 * n);
    std::vector<double> J((size_t)2 * n * nf, 0.0);
    std::vector<double> A((size_t)nf * nf, 0.0), g((size_t)nf, 0.0);
    std::vector<double> Ad((size_t)nf * nf, 0.0), gd((size_t)nf, 0.0);
    auto cost = [&](const double* p) {
        residuals(p, rc);
        double s = 0.0;
        for (double v : rc) s += v * v;
        return s;
    };

    double f = cost(x);
    if (!std::isfinite(f)) return false;
    double lambda = 1e-3;
    for (int iter = 0; iter < maxIter; ++iter) {
        residuals(x, r);
        for (int j = 0; j < nf; ++j) {
            const int pj = idx[j];
            const double hstep = 1e-5 * std::max(1.0, std::fabs(x[pj]));
            const double save = x[pj];
            x[pj] = save + hstep;
            residuals(x, rt);
            x[pj] = save;
            for (size_t i = 0; i < 2 * n; ++i) J[i * nf + j] = (rt[i] - r[i]) / hstep;
        }
        std::fill(A.begin(), A.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);
        for (int a = 0; a < nf; ++a) {
            for (size_t i = 0; i < 2 * n; ++i) g[(size_t)a] += J[i * nf + a] * r[i];
            for (int b = a; b < nf; ++b) {
                double s = 0.0;
                for (size_t i = 0; i < 2 * n; ++i) s += J[i * nf + a] * J[i * nf + b];
                A[(size_t)a * nf + b] = s;
                A[(size_t)b * nf + a] = s;
            }
        }
        bool okStep = false;
        for (int tr = 0; tr < 10 && !okStep; ++tr) {
            std::copy(A.begin(), A.end(), Ad.begin());
            std::copy(g.begin(), g.end(), gd.begin());
            for (int a = 0; a < nf; ++a) Ad[(size_t)a * nf + a] += lambda * (Ad[(size_t)a * nf + a] + 1e-6);
            for (int a = 0; a < nf; ++a) gd[a] = -gd[a];
            if (!solveSPD(Ad, gd, nf)) { lambda *= 10.0; continue; }
            double xt[6];
            for (int i = 0; i < 6; ++i) xt[i] = x[i];
            for (int a = 0; a < nf; ++a) {
                const int pj = idx[a];
                xt[pj] = std::min(hi[pj], std::max(lo[pj], xt[pj] + gd[a]));
            }
            const double ft = cost(xt);
            if (std::isfinite(ft) && ft < f) {
                const double rel = (f - ft) / std::max(1e-30, f);
                // Convergence criterion: stop directly when both the relative cost decrease and the parameter step are negligible.
                // The old version only checked rel < 1e-13, so in practice it almost always ran the full maxIter iterations (nf+1
                // residual evaluations per iteration), the main cost of the Part TWO refinement stage.
                double stepRel = 0.0;
                for (int a = 0; a < nf; ++a) {
                    const int pj = idx[a];
                    stepRel = std::max(stepRel,
                                       std::fabs(gd[(size_t)a]) / std::max(1.0, std::fabs(x[pj])));
                }
                for (int i = 0; i < 6; ++i) x[i] = xt[i];
                f = ft;
                lambda = std::max(1e-9, lambda / 5.0);
                okStep = true;
                if (rel < 1e-12 || stepRel < 1e-10) iter = maxIter;
            } else {
                lambda *= 8.0;
            }
        }
        if (!okStep) break;
    }
    cell = expand(x);
    return true;
}

Cell refineAssign(const std::vector<Peak>& peaks, const Orient& orient, const std::vector<Hkl>& cands,
                  Cell cell, CrystalSystem sys, int iters, double& rmsOut, double maxG) {
    const std::array<bool, 6> freeMask = systemFreeMask(sys);
    applySystemConstraints(sys, cell);
    const size_t n = peaks.size();
    if (n == 0) {
        rmsOut = 0.0;
        return cell;
    }

    // Graduated assignment: sort by |q| ascending, first use only the shortest few peaks to establish the cell,
    // then add long-vector peaks in stages.
    // Why this is necessary: alternating "nearest-neighbour assignment + LM" gets locked onto index aliases of long-vector peaks
    // when the initial value deviates by 1~2%. Example: for an orthorhombic cell a=10,b=12,c=17, (2,2,0) and (1,3,0) differ in |q| by only
    // 0.03; with a 2% initial offset the alias is actually closer, and once assigned, LM pushes a,b toward accommodating the alias
    // side, after which every round reproduces the same alias (a self-consistent wrong minimum); no amount of better initialisation recovers it.
    // Short-vector peak indices are unique under the same initial offset; using them to first tighten the cell makes the alias disappear.
    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&peaks](size_t a, size_t b) {
        const double qa = peaks[a].qxy * peaks[a].qxy + peaks[a].qz * peaks[a].qz;
        const double qb = peaks[b].qxy * peaks[b].qxy + peaks[b].qz * peaks[b].qz;
        return qa < qb;
    });

    constexpr int kStage = 6;
    std::vector<size_t> stages;   // peaks included at each stage (strictly increasing, merged automatically when few peaks)
    for (int st = 1; st <= kStage; ++st) {
        const size_t m = (st == kStage)
                             ? n
                             : std::min(n, std::max<size_t>(5, (n * (size_t)st) / kStage + 1));
        if (stages.empty() || m > stages.back()) stages.push_back(m);
    }
    std::vector<Peak> sub;
    sub.reserve(n);
    std::vector<Hkl> ah;
    for (size_t si = 0; si < stages.size(); ++si) {
        const size_t m = stages[si];
        sub.clear();
        for (size_t i = 0; i < m; ++i) sub.push_back(peaks[order[i]]);
        // The first few stages do only 2 rounds (small subsets, cheap iterations); the last stage uses the caller's iteration budget
        const bool last = (si + 1 == stages.size());
        const int nit = last ? std::max(1, iters) : 2;
        for (int it = 0; it < nit; ++it) {
            const double rms = assignHkl(cell, orient, sub, cands, ah, nullptr, 1e30, 1.0, maxG);
            if (!std::isfinite(rms) || rms > 1e29) break;
            Cell next = cell;
            if (!refineCell(sub, ah, orient, next, freeMask, 80, sys)) break;
            applySystemConstraints(sys, next);
            const double d = std::fabs(next.a - cell.a) + std::fabs(next.b - cell.b) +
                             std::fabs(next.c - cell.c) + std::fabs(next.al - cell.al) +
                             std::fabs(next.be - cell.be) + std::fabs(next.ga - cell.ga);
            cell = next;
            if (d < 1e-9) break;
        }
    }
    rmsOut = assignHkl(cell, orient, peaks, cands, ah, nullptr, 1e30, 1.0, maxG);
    return cell;
}

Quality evaluate(const Cell& cell, const Orient& orient, const std::vector<Peak>& peaks,
                 const std::vector<Hkl>& hkl) {
    Quality q;
    const size_t n = peaks.size();
    if (n == 0 || hkl.size() != n) return q;
    const std::array<double, 6> G = recipMetricOf(reciprocalMatrix(cell));
    const std::array<double, 3> nu = orient.get(G);
    double sx = 0.0, sz = 0.0, sg = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double qx = 0.0, qz = 0.0;
        forwardNu(G, nu, hkl[i].h, hkl[i].k, hkl[i].l, qx, qz);
        const double g2 = G[0] * hkl[i].h * hkl[i].h + G[1] * hkl[i].k * hkl[i].k +
                          G[2] * hkl[i].l * hkl[i].l +
                          2.0 * (G[3] * hkl[i].h * hkl[i].k + G[4] * hkl[i].h * hkl[i].l +
                                 G[5] * hkl[i].k * hkl[i].l);
        const double dx = qx - peaks[i].qxy;
        const double dz = qz - peaks[i].qz;
        const double dg = std::sqrt(std::max(0.0, g2)) -
                          std::sqrt(peaks[i].qxy * peaks[i].qxy + peaks[i].qz * peaks[i].qz);
        sx += dx * dx;
        sz += dz * dz;
        sg += dg * dg;
        q.maxAbs = std::max(q.maxAbs, std::max(std::fabs(dx), std::fabs(dz)));
    }
    q.rmsQxy = std::sqrt(sx / (double)n);
    q.rmsQz = std::sqrt(sz / (double)n);
    q.rmsQxyz = std::sqrt(sg / (double)n);
    q.used = (int)n;
    return q;
}

double assignMaxG(const Cell& cell, const Orient& orient, const std::vector<Peak>& peaks,
                  const std::vector<Hkl>& hkl) {
    if (hkl.size() != peaks.size() || peaks.empty()) return 1e30;
    const std::array<double, 6> G = recipMetricOf(reciprocalMatrix(cell));
    const std::array<double, 3> nu = orient.get(G);
    double qmax = 0.0, dmax = 0.0;
    int used = 0;
    for (size_t i = 0; i < peaks.size(); ++i) {
        qmax = std::max(qmax, std::hypot(peaks[i].qxy, peaks[i].qz));
        const Hkl& x = hkl[i];
        if (x.h == 0 && x.k == 0 && x.l == 0) continue; // unassigned peaks do not contribute to dmax
        double qx = 0.0, qz = 0.0;
        forwardNu(G, nu, x.h, x.k, x.l, qx, qz);
        dmax = std::max(dmax, std::hypot(qx - peaks[i].qxy, qz - peaks[i].qz));
        ++used;
    }
    // Without a valid assignment the deviation cannot be estimated; fall back to the full box (no filtering).
    if (used == 0 || !std::isfinite(dmax)) return 1e30;
    return qmax + 3.0 * dmax + 0.05;
}

bool refineConstrained(const std::vector<Peak>& peaks, const Orient& orient,
                       const std::vector<Hkl>& cands, CrystalSystem sys, Cell& cell,
                       std::vector<Hkl>& hkl, double maxG, double rmsGood, double rmsWorst) {
    if (sys == CrystalSystem::Triclinic || peaks.empty()) return false;
    applySystemConstraints(sys, cell);
    // Trigger condition: only solutions still "not good enough" after refinement get rescued. The threshold is 0.3× the user Δq_xyz tolerance:
    // for solutions already far below tolerance (≤ 30%), alias reassignment has no residual to recover, and skipping saves
    // most of the rescue-refinement cost (the original implementation ran unconditionally, taking ~78% of the second stage's time).
    if (std::isfinite(rmsGood) || std::isfinite(rmsWorst)) {
        const Quality q0 = evaluate(cell, orient, peaks, hkl);
        if (q0.rmsQxyz <= rmsGood) return false;
        // Solutions far outside tolerance: 6 LM iterations cannot bring the residual back into the convergence basin, so skip the rescue.
        if (q0.rmsQxyz > rmsWorst) return false;
    }
    // Baseline: joint residual of reassigning the current cell within the (maxG-limited) candidate box
    std::vector<Hkl> base;
    const double rmsBefore = assignHkl(cell, orient, peaks, cands, base, nullptr, 1e30, 1.0, maxG);
    if (!std::isfinite(rmsBefore)) return false;

    double rmsTmp = 0.0;
    Cell alt = refineAssign(peaks, orient, cands, cell, sys, 6, rmsTmp, maxG);
    applySystemConstraints(sys, alt);
    std::vector<Hkl> altHkl;
    const double rmsAfter = assignHkl(alt, orient, peaks, cands, altHkl, nullptr, 1e30, 1.0, maxG);
    if (!std::isfinite(rmsAfter) || !(rmsAfter < rmsBefore - 1e-9)) return false;
    cell = alt;
    hkl = std::move(altHkl);
    return true;
}

} // namespace gid
