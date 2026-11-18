// ============================================================================
//  fixedplane.h — indexing with a user-fixed contact plane (u v w)
//  Uses no known cell parameters: searches cell a,b,c,α,β,γ together with each peak's (hkl),
//  least-squares fitting (q_xy, q_z) over all peaks, and cuts degrees of freedom per crystal system.
//
//  Three paths; returns as soon as any one succeeds (see the top comment in fixedplane.cpp):
//    0) Initial-cell fast path: refineAssign (alternating LM iterations)
//    1) Two-step method when specular rod peaks exist (fixedUvw + rod fundamental p)
//    2) FP-DL deterministic search (Fixed-Plane Direct Linear):
//         three peaks -> 3×3 integer Cramer solve for ν -> "G tᵀ = (t·ν)ν" + "q² = x G xᵀ" linear LS
//         -> alternate "nearest-neighbour assignment + linear re-solve" -> LM polish on survivors only
//  The old "random-seed RANSAC" was removed: the assignment search space is |cands|^nSeed (728⁵ ≈ 2e14 by default),
//  so 2000 random trials always had a hit probability of 0; it had degenerated into "local optimization with random restarts".
//
//  Semantics of the contact plane (u,v,w) (important):
//    The model depends only on its **direction** — ν = G nᵀ/√(n G nᵀ) is invariant under n → m·n,
//    so the entry point normalizes by gcd to a primitive direction; when reporting, qspec is the rod fundamental.
//    The observed q_z on a rod is n·fundamental (n unknown); anywhere q_spec is treated as a hard constraint one must
//    use the fundamental p (detectRodPeriod), never the minimum observed q_z.
// ============================================================================
#pragma once

#include "types.h"

namespace gid {

struct FixedPlaneParams {
    int u = 0, v = 0, w = 1;              // fixed contact plane (entry normalizes to a primitive direction)
    int maxH = 4, maxK = 4, maxL = 4;     // hkl search range (also the FP-DL index box)
    CrystalSystem system = CrystalSystem::Triclinic;
    Cell start{10.0, 10.0, 15.0, 90.0, 90.0, 90.0};

    // ---- FP-DL deterministic search (main path when no specular rod peak) ------------------
    int seedPeaks = 8;         // seed pool = the N ordinary peaks with smallest |q| (q_xy > rodTol)
    int seedShell = 6;         // seed index shell |hkl|² <= seedShell (seed 1 uses only the 12 shortest)
    int minRowsMetric = 6;     // minimum row count for the first (A) system (6 for triclinic, else clamped to nf+3)
    int beam = 1500;           // number of ν candidates kept
    int maxCands = 400000;     // enumeration size gate (×100 bounds the actual candidate count)
    double qzTolMul = 3.0;     // q_z compatibility tolerance = max(qzTolMul*acceptTol, 0.02)
    double consistTol = 0.05;  // constraint-row relative residual gate |G tᵀ − (t·ν)ν| / (t·ν)²

    // ---- specular rod peaks (order / fundamental) ------------------------------------------
    double rodTol = 0.02;      // |q_xy| <= rodTol is treated as a rod peak
    int rodNmax = 6;           // highest order on the rod (also used for the p/n scan under single-order ambiguity)
    double dqspec = 0.01;      // tolerance on |q_spec(solution) − fundamental|

    double acceptTol = 0.05;
    int topN = 10;
    int refineIter = 6;
    // Current state: q_xy is already forced to |q_xy| at the data entry (normalizeQxySigns), so this switch no longer
    // affects the result; the field is kept only for compatibility with existing UI/parameter interfaces (on == off).
    bool autoSign = true;                 // legacy switch: q_xy sign compensation (now equivalent)
    int threads = -1;
};

std::vector<Solution> runFixedPlane(const std::vector<Peak>& peaks, const FixedPlaneParams& prm,
                                    const Progress& prog);

} // namespace gid
