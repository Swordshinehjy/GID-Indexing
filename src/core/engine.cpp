// ============================================================================
//  engine.cpp — see engine.h
// ============================================================================
#include "engine.h"

#include <algorithm>
#include <cmath>

namespace gid {

const char* modeName(IndexMode m) {
    switch (m) {
        case IndexMode::WithSpec: return "With specular (two-step)";
        case IndexMode::NoSpec: return "No q_spec (one-step)";
        case IndexMode::FixedPlane: return "Fixed contact plane (uvw)";
        case IndexMode::NoSpecTwoSteps: return "No q_spec (two-step)";
        case IndexMode::Simple: return "Simple indexing (out-of-plane 001)";
    }
    return "";
}

RodPeriod detectRodPeriod(const std::vector<Peak>& peaks, double rodTol, int nmax) {
    RodPeriod rp;
    for (int i = 0; i < (int)peaks.size(); ++i) {
        const Peak& q = peaks[i];
        if (!(std::fabs(q.qxy) <= rodTol)) continue;
        if (!(q.qz > 1e-6) || !std::isfinite(q.qz)) continue;
        rp.peakIdx.push_back(i);
    }
    rp.nRod = (int)rp.peakIdx.size();
    if (rp.nRod == 0) return rp;

    std::vector<double> z;
    z.reserve(rp.peakIdx.size());
    for (int i : rp.peakIdx) z.push_back(peaks[i].qz);
    std::sort(z.begin(), z.end());
    z.erase(std::unique(z.begin(), z.end(),
                        [](double a, double b) { return std::fabs(a - b) < 1e-6; }),
            z.end());
    // With only one independent q_z, n cannot be determined from the data (p could be z/1, z/2 …)
    rp.ambiguous = (z.size() == 1);

    // Take "the smallest n that explains all z_i": z_i / (z_min/n) must all be close to positive integers.
    // Taking the smallest n is equivalent to "prefer p = the smallest observed value"; the fundamental frequency is compressed only when multi-order structure forces it.
    const double relTol = 0.15;
    for (int n = 1; n <= std::max(1, nmax); ++n) {
        const double cand = z.front() / (double)n;
        if (!(cand > 1e-9)) break;
        bool ok = true;
        for (double zi : z) {
            const double r = zi / cand;
            const double rn = std::round(r);
            if (!(rn >= 1.0) || std::fabs(r - rn) > relTol) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;
        rp.ok = true;
        rp.p = cand;
        rp.qz = z;
        rp.order.resize(z.size());
        for (size_t i = 0; i < z.size(); ++i) rp.order[i] = (int)std::lround(z[i] / rp.p);
        break;
    }
    return rp;
}

bool detectQSpec(const std::vector<Peak>& peaks, double& qspec, int& rodCount) {
    // Keep the strict-zero criterion of the UI / tool side (they feed this number directly as q_spec)
    const RodPeriod rp = detectRodPeriod(peaks, 1e-6, 6);
    qspec = rp.p;
    rodCount = rp.nRod;
    return rp.ok;
}

void alignSolution(const std::vector<Peak>& peaks, Solution& s) {
    const size_t n = peaks.size();
    // Before alignment: hkl / dist cover only peaks with q_xy > 0, in per-index correspondence
    const std::vector<Hkl> hklIn = s.hkl;
    const std::vector<double> distIn = s.dist;
    if (hklIn.size() == n) {
        // P1-3: dist must be the same length as hkl, otherwise the later "dist[i] paired with peaks[i]" silently misaligns
        if (s.dist.size() != n) s.dist.assign(n, std::nan(""));
        return;
    }
    std::vector<Hkl> full(n, Hkl{0, 0, 0});
    std::vector<double> dist(n, std::nan(""));
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        // Inputs already had |q_xy| taken uniformly at the solve entry; using |q_xy| here only protects direct callers.
        if (std::fabs(peaks[i].qxy) > 1e-6) {
            if (j < hklIn.size()) {
                full[i] = hklIn[j];
                if (j < distIn.size()) dist[i] = distIn[j];
            }
            ++j;
        } else if (s.hasPlane && s.qspec > 1e-9) {
            // Point on the specular rod: index is m*(u,v,w)
            const double m = peaks[i].qz / s.qspec;
            const int mi = (int)std::round(m);
            if (std::fabs(m - mi) < 0.15 && mi != 0) {
                full[i] = Hkl{s.uvw.h * mi, s.uvw.k * mi, s.uvw.l * mi};
                const QPoint q = predictQ(s, full[i].h, full[i].k, full[i].l);
                dist[i] = std::hypot(q.qxy - peaks[i].qxy, q.qz - peaks[i].qz);
            }
        }
        // hasPlane == false (P2-4): rod points have no trustworthy contact-plane model, so no indices are fabricated:
        // hkl stays (0,0,0), dist stays NaN.
    }
    s.hkl = std::move(full);
    s.dist = std::move(dist);
}

EngineResult runIndexing(const std::vector<Peak>& peaks, const EngineParams& prm,
                         const Progress& prog) {
    EngineResult res;
    if (peaks.empty()) {
        res.message = "No peak data";
        return res;
    }
    // Computations use |q_xy| uniformly: the sign only reflects the mirror handedness of the detection geometry, not structural information.
    // Normalise once here so the solvers and alignSolution see the same peak data.
    std::vector<Peak> pk = peaks;
    normalizeQxySigns(pk);
    std::vector<Solution> sols;
    switch (prm.mode) {
        case IndexMode::WithSpec:
            if (!(prm.spec.qspec > 0.0)) {
                res.message = "q_spec not specified: enter it in the settings, or keep q_xy=0 specular rod peaks in the data";
                return res;
            }
            sols = runSpecIndex(pk, prm.spec, prog);
            break;
        case IndexMode::NoSpec:
            sols = runNoSpec(pk, prm.nospec, prog);
            break;
        case IndexMode::FixedPlane:
            sols = runFixedPlane(pk, prm.fixed, prog);
            break;
        case IndexMode::NoSpecTwoSteps:
            sols = nospec2::runNoSpecTwoSteps(pk, prm.nospec2, prog);
            break;
        case IndexMode::Simple:
            sols = runSimpleIndex(pk, prm.simple, prog);
            break;
    }
    if (prog.stop()) {
        res.stopped = true;
        res.message = "Stopped";
        return res;
    }
    // The two-step branch of the fixed-contact-plane mode (including specular rod peaks) produces hkl covering only
    // peaks with q_xy>0, returned in that subset's order, not per-index aligned with the input peaks (rod peaks are dropped inside the two-step method);
    // the same alignment must be done, otherwise predictQ maps hkl to the wrong peaks, showing as "solution does not match the panel qxy/qz".
    for (Solution& s : sols) alignSolution(peaks, s);
    res.solutions = std::move(sols);
    res.ok = !res.solutions.empty();
    if (res.solutions.empty())
        res.message = "No solution meeting the criteria: relax the error thresholds / widen the index range / adjust the unit-cell range";
    else
        res.message = "Found " + std::to_string(res.solutions.size()) + " candidate solutions";
    return res;
}

} // namespace gid
