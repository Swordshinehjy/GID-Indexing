// ============================================================================
//  simpleindex.h — fifth indexing mode: simple indexing (out-of-plane fixed to 001, automatic search)
//
//  Model (consistent with the analyzer's "simple indexing" formulas; the out-of-plane contact plane defaults to (0 0 1)):
//      q_xy² = h²·A + k²·B − 2·h·k·C
//      q_z   = h·P + k·Q + l·R
//  where A=a*∥², B=b*∥², C=cosγ∥, P=a*⊥, Q=b*⊥, R=c*.
//
//  Given (hkl) the model is linear in the 6 parameters: one 3x3 in-plane and one 3x3 out-of-plane system, both closed-form solvable.
//  The search has three stages:
//      S1 generate (A,B,C) hypotheses from q_xy only (three-peak minimal sample + all-peak consensus);
//      S2 estimate R from q_z only and run a Hough search for tilt on the torus (P/R,Q/R);
//      S3 nearest-neighbor reassignment + alternating refinement of the two linear least-squares systems.
//
//  The solution is still represented by an ordinary reciprocal metric:
//      G* = [ A+P²,  P·Q−C, P·R
//             P·Q−C, B+Q²,  Q·R
//             P·R,    Q·R,   R²  ]
//  Thus nu=(P,Q,R), uvw=(0,0,1), and the existing predictQ / UI / export paths can be used directly.
// ============================================================================
#pragma once

#include "types.h"

namespace gid {

struct SimpleParams {
    double tol = 0.05;    // full-model 2D match tolerance (Å⁻¹)
    double tolXy = 0.05;  // in-plane q_xy tolerance
    double tolZ = 0.05;   // out-of-plane q_z tolerance
    double rodTol = 0.04; // q_xy <= rodTol is treated as a specular rod peak
    double qMin = 0.05;   // ignore peaks near the origin

    int hkLimit = 2;   // in-plane index box |h|,|k| <= hkLimit
    int lLimit = 6;    // out-of-plane order box |l| <= lLimit
    int maxPlanes = 12; // number of in-plane hypotheses kept
    int maxTilts = 4;   // tilt candidates kept per (in-plane, R)
    int topK = 10;      // number of solutions output
    int icpIters = 4;   // reassignment iterations
    int grid = 33;      // tilt Hough grid
};

// Automatic search. Returns candidate solutions sorted by descending inlier count -> ascending RMS -> ascending index complexity.
std::vector<Solution> runSimpleIndex(const std::vector<Peak>& peaks, const SimpleParams& prm,
                                     const Progress& prog);

} // namespace gid
