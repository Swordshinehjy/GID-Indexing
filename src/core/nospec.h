// ============================================================================
//  nospec.h — one-step indexing without a specular peak (unknown q_spec)
//
//  Corresponds to GIDInd paper sections 2.1.2 / 3.2 and Appendix C (docs/gidind_nospec.py):
//    9 unknowns: the 6 independent components of the reciprocal metric G* (determine a,b,c,α,β,γ)
//                + the out-of-plane coefficient nu = (a*.n, b*.n, c*.n) (determines orientation)
//    The two equation sets are linear in and share no unknowns:
//      (A) q_xyz^2 = q_xy^2 + q_z^2 = h^2 G11 + k^2 G22 + l^2 G33
//                    + 2hk G12 + 2hl G13 + 2kl G23
//      (B) q_z     = h*nu_0 + k*nu_1 + l*nu_2
//    Coupling is only through the physical consistency nu^T (G*)^-1 nu = 1 (|n_hat| = 1)
//
//  Flow (progressive expansion, keeping the 74^6 brute force to a manageable size):
//    S0 drop zero peaks, sort ascending by q_xyz (smallest q -> smallest Laue indices)
//    S1 3 seed peaks: enumerate (h,k,l) in the seed box, solve nu exactly by 3x3, integer |det|>=1 prefilter
//    S2 beam search progressive expansion: add one peak at a time, predict q_z from the current nu, keep candidates within tolerance
//    S3 once 6 rows (4 for monoclinic) are reached, solve (A) -> G* -> full unit cell
//    S4 physical/crystallographic filters: positive definiteness, ranges, consistency nu^T G*^-1 nu = 1, Niggli
//    S5 joint reassignment: choose hkl over (-6..6)^3 by (q_xy,q_z) weighted residual, then re-solve (A)(B)
//    S6 reduction + deviation statistics + deduplication + sorting
// ============================================================================
#pragma once

#include "types.h"

#include <algorithm>
#include <cmath>

namespace gid {

// Maximum rows per state (= maximum number of peaks participating in the search). State in nospec.cpp allocates fixed arrays by this.
constexpr int kNoSpecMaxRows = 16;

// Solution methods without spec (see the "no q_spec method" section of docs/algorithms-spec.md):
//   * Beam   — beam search (default): S1 seeding + S2 per-peak expansion + S3 metric solve. Discriminating power comes from search
//              breadth; slow but stable, reliable for all crystal systems (including high-symmetry cubic/hexagonal).
//   * Sphere — assignment replaced by LLL + sphere decoding: S5's "nearest neighbour in the grid box" becomes (A) ellipsoid
//              (Cholesky + LLL reduction, then Fincke–Pohst enumeration) ∩ (B) slab integer search,
//              free of the reassignment grid box, and naturally searches the whole integer space.
enum class NoSpecMethod { Beam = 0, Sphere = 1 };

struct NoSpecParams {
    CrystalSystem system = CrystalSystem::Triclinic;

    // Solution method (see above). Default = beam search (most stable: recovers the true cell even for highly
    // symmetric systems; all measured DIP + synthetic regression cases pass).
    NoSpecMethod method = NoSpecMethod::Beam;

    // Index enumeration box (each axis takes an independent range). A long axis (alkyl chain / long c) usually needs l much
    // larger than h,k; e.g. the manual solution for the DIP data contains (-1,0,6): with an l limit of only 2 the true solution is not in the search space at all.
    int seedH = 2, seedK = 2, seedL = 6;
    // Index shell of the seed peaks: |hkl|² <= seedShell. Seeds are the peaks with smallest |g|, whose indices must be short
    // vectors; using a shell instead of the whole box keeps the B³ seed enumeration manageable (B is the box volume; the shell is only a small part).
    int seedShell = 6;
    int seedPeaks = 6;  // take any 3 of the smallest |g| peaks as seeds
    int minRows = 8;    // minimum rows needed to solve (A) (triclinic >=8; the row surplus provides discrimination)
    // During expansion, prune by the least-squares row residual of (A) (the q² equation): spurious states are dropped as soon as the rows suffice to solve the metric.
    // <=0 means metricPruneMul * dqxyzCutoff is used.
    double metricPrune = 0.0;
    double dqzCutoff = 0.01;
    double dqxyzCutoff = 0.05;
    double dqxyCutoff = 0.05;

    // ---- scale factors of the four criteria (P2-1) -----------------------------------------
    // Previously hard-coded inside nospec.cpp, so they could not be adjusted uniformly when data noise changed.
    double qzTolMul = 3.0;       // q_z compatibility tolerance during search = max(qzTolMul*dqzCutoff, qzTolFloor)
    double qzTolFloor = 0.02;    // absolute floor of that tolerance (q units)
    double metricPruneMul = 0.4; // (A) row-residual pruning threshold (used when metricPrune <= 0)
    double maxAbsMul = 3.0;      // per-peak maximum deviation gate = maxAbsMul * dqxyzCutoff
    double rmsFloorMul = 0.25;   // "same tier" residual floor when sorting = rmsFloorMul * dqxyzCutoff

    int gridH = 6, gridK = 6, gridL = 6; // joint reassignment grid
    int beam = 8000;                     // beam width
    int refineTop = 300;
    int maxSolutions = 50;

    bool restrictCell = false;
    double aLo = 3, aHi = 30, bLo = 3, bHi = 30, cLo = 3, cHi = 30;
    double alLo = 60, alHi = 120, beLo = 60, beHi = 120, gaLo = 60, gaHi = 120;
    double volLo = 50, volHi = 100000;

    int threads = -1;

    // Entry parameter sanitization: the UI lets users type arbitrary integers; out-of-range values cause two kinds of accidents
    //   * |h|,|k|,|l| > 127 — indices in State are stored as int8_t, so writing truncates;
    //   * enumeration box / grid / beam width unbounded — the search node count explodes and the UI freezes for a long time.
    // Here all size parameters are clamped to safe ranges, and the downstream clamp precondition lo <= hi is guaranteed.
    NoSpecParams sanitized() const {
        auto pos = [](double v, double dflt) { return (std::isfinite(v) && v > 0.0) ? v : dflt; };
        auto nonneg = [](double v, double dflt) { return (std::isfinite(v) && v >= 0.0) ? v : dflt; };
        NoSpecParams p = *this;
        p.method = static_cast<NoSpecMethod>(
            std::clamp(static_cast<int>(method), 0, 2));
        p.seedH = std::clamp(seedH, 0, 32);
        p.seedK = std::clamp(seedK, 0, 32);
        p.seedL = std::clamp(seedL, 0, 32);
        p.seedShell = std::clamp(seedShell, 1, 64);
        p.seedPeaks = std::clamp(seedPeaks, 3, 8);
        p.minRows = std::clamp(minRows, 3, kNoSpecMaxRows);
        p.gridH = std::clamp(gridH, 0, 8);
        p.gridK = std::clamp(gridK, 0, 8);
        p.gridL = std::clamp(gridL, 0, 8);
        p.beam = std::clamp(beam, 1, 100000);
        p.refineTop = std::clamp(refineTop, 1, 2000);
        p.maxSolutions = std::clamp(maxSolutions, 1, 500);
        p.dqzCutoff = pos(dqzCutoff, 0.01);
        p.dqxyzCutoff = pos(dqxyzCutoff, 0.05);
        p.dqxyCutoff = pos(dqxyCutoff, 0.05);
        p.qzTolMul = std::clamp(pos(qzTolMul, 3.0), 0.0, 20.0);
        p.qzTolFloor = std::clamp(nonneg(qzTolFloor, 0.02), 0.0, 1.0);
        p.metricPruneMul = std::clamp(pos(metricPruneMul, 0.4), 0.0, 20.0);
        p.maxAbsMul = std::clamp(pos(maxAbsMul, 3.0), 0.0, 50.0);
        p.rmsFloorMul = std::clamp(pos(rmsFloorMul, 0.25), 0.0, 20.0);
        if (!(std::isfinite(metricPrune) && metricPrune > 0.0)) p.metricPrune = 0.0;
        return p;
    }
};

std::vector<Solution> runNoSpec(const std::vector<Peak>& peaks, const NoSpecParams& prm,
                                const Progress& prog);

} // namespace gid
