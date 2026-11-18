// ============================================================================
//  types.h — Input/output structures shared by the three indexing modes
// ============================================================================
#pragma once

#include "cell.h"
#include "refine.h"

#include <cmath>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

namespace gid {

// Normalisation of the q_xy sign for computation inputs.
//
// The sign of q_xy depends only on the mirror handedness of the detection geometry and is not structural information: the model side (forwardNu /
// forwardHkl) always gives a non-negative q_xy. To avoid inconsistent conventions from the GUI, file reading, batch processing, and reproducibility tools
// each handling the sign differently, take |q_xy| uniformly here, and call it only once before data enters the solve.
inline void normalizeQxySigns(std::vector<Peak>& peaks) {
    for (Peak& p : peaks) p.qxy = std::fabs(p.qxy);
}

// Helper that only checks whether an environment variable is set: MSVC reports C4996 for a bare getenv (_CRT_SECURE_NO_WARNINGS
// is not necessarily enabled in every project), so wrap it uniformly here. Inside the kernel it only toggles debug output (GID_DEBUG).
inline bool hasEnv(const char* name) {
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
    return std::getenv(name) != nullptr;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
}

// One indexing solution: contact plane (uvw) + unit cell + (hkl) per peak + deviation statistics
//
// Contract (P1-3 / P2-4):
//   * After alignSolution, hkl and dist are both the same length as the input peaks and correspond one to one by index;
//   * For peaks no model can explain (rod points + hasPlane == false), hkl is (0,0,0) and dist is NaN;
//   * hasPlane == false means nu does not correspond to any integer contact plane (a common result of the no-q_spec mode),
//     in which case the values of uvw / qspec / dqspec are meaningless, and the UI should display "—".
struct Solution {
    Hkl uvw;                          // contact plane indices (set to (0,0,0) when hasPlane == false)
    bool hasPlane = true;             // false = no integer contact plane exists
    std::array<double, 3> nu{0, 0, 0}; // out-of-plane component coefficients (for exact forward modelling)
    Cell cell;
    double vol = 0.0;
    double qspec = 0.0;      // |G_spec| computed from the solution (0 when hasPlane == false)
    double rmsQxyz = 0.0;
    double rmsQxy = 0.0;
    double rmsQz = 0.0;
    double dqspec = 0.0;     // |q_spec(solution) - q_spec(measured)|
    std::vector<Hkl> hkl;    // one entry per peak
    std::vector<double> dist;// one entry per peak, same length as hkl
};

// Contact-plane sign representative (unified convention across Part1 methods)
//
// (u,v,w) and (−u,−v,−w) are opposite normals of the same crystal planes, and the model is strictly symmetric to them:
// q_spec = sqrt(n G* nᵀ) is an even function; and after negating hkl ((n,hkl) → (−n,−hkl)),
// each peak's q_xy and q_z are unchanged, so the unit cell, rms, and peak assignment all stay the same; only the sign representative differs.
// Unified rule (deterministic, independent of enumeration order / thread scheduling):
//   1) take the branch with larger u+v+w (i.e. u+v+w > 0);
//   2) when u+v+w == 0 (e.g. (1,−1,0), (0,1,−1)), take the branch with "first non-zero component > 0".
// Call sites: each method calls it before deduplication/output; fixedUvw (user-specified direction) does not call it.
inline bool uvwIsCanonical(const Hkl& n) {
    const int s = n.h + n.k + n.l;
    if (s != 0) return s > 0;
    if (n.h != 0) return n.h > 0;
    if (n.k != 0) return n.k > 0;
    return n.l > 0;
}

inline void canonicalizeUvw(Solution& s) {
    if (!s.hasPlane) return;
    if (s.uvw.h == 0 && s.uvw.k == 0 && s.uvw.l == 0) return;
    if (uvwIsCanonical(s.uvw)) return;
    s.uvw = Hkl{-s.uvw.h, -s.uvw.k, -s.uvw.l};
    s.nu = {-s.nu[0], -s.nu[1], -s.nu[2]};
    for (Hkl& x : s.hkl) x = Hkl{-x.h, -x.k, -x.l};
}

// Unified progress/cancellation interface (the kernel must not touch the UI)
struct Progress {
    std::function<void(double, const std::string&)> report; // value 0..1
    std::function<bool()> cancelled;
    void operator()(double v, const std::string& msg) const {
        if (report) report(v, msg);
    }
    bool stop() const { return cancelled && cancelled(); }
};

// Forward model of a solution: prefer nu (a no-q_spec solution has only the out-of-plane normal n̂, with no trustworthy integer contact plane (u,v,w);
// using uvw for forward modelling would make hkl mismatch the measured (q_xy,q_z)), otherwise use the integer contact plane.
// Both modes give identical results when an integer contact plane exists, so this single path suffices.
inline QPoint predictQ(const Solution& s, int h, int k, int l) {
    const double n2 = s.nu[0] * s.nu[0] + s.nu[1] * s.nu[1] + s.nu[2] * s.nu[2];
    if (n2 > 1e-12) {
        double qxy = 0.0, qz = 0.0;
        forwardNu(recipMetric(s.cell), s.nu, h, k, l, qxy, qz);
        return QPoint{qxy, qz};
    }
    return forwardHkl(s.cell, s.uvw, h, k, l);
}

inline std::string fmtCell(const Cell& c) {
    char buf[160];
    std::snprintf(buf, sizeof(buf), "a=%.4f b=%.4f c=%.4f  al=%.3f be=%.3f ga=%.3f", c.a, c.b, c.c,
                  c.al, c.be, c.ga);
    return std::string(buf);
}

} // namespace gid
