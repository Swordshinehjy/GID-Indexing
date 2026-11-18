// ============================================================================
//  refine.h — Numerical refinement: hkl nearest-neighbour assignment + Levenberg-Marquardt unit-cell fitting + quality statistics
//  Shared by the three indexing modes:
//    * With-q_spec two-step method: LM refinement after Part TWO (review A10 only for Top-K)
//    * No-q_spec one-step method: joint reassignment and refitting in S5/S6
//    * Fixed contact plane (uvw): LM fitting from random seed initial values with alternating iteration
//
//  Unified forward model (needs only the out-of-plane component coefficients nu and the reciprocal metric G*):
//      q_z  = h*nu_0 + k*nu_1 + l*nu_2
//      |g|^2 = (h,k,l) G* (h,k,l)^T
//      q_xy = sqrt(|g|^2 - q_z^2)
//  For an integer contact plane (uvw), nu = G* (u,v,w)^T / q_spec (varies with the unit cell);
//  without specular peaks, nu is fitted directly from the q_z equations (Orient::useNu).
// ============================================================================
#pragma once

#include "cell.h"

#include <vector>

namespace gid {

struct Peak {
    double qxy = 0.0;
    double qz = 0.0;
};

// Out-of-plane direction: either an integer contact plane (uvw) or a direct nu
struct Orient {
    Hkl uvw{0, 0, 1};
    bool useNu = false;
    std::array<double, 3> nu{0, 0, 0};

    // Solve nu from the current reciprocal metric
    std::array<double, 3> get(const std::array<double, 6>& G) const;
    static Orient fromUvw(const Hkl& uvw) { return Orient{uvw, false, {0, 0, 0}}; }
    static Orient fromNu(const std::array<double, 3>& nu) { return Orient{{0, 0, 1}, true, nu}; }
};

struct Quality {
    double rmsQxy = 0.0;
    double rmsQz = 0.0;
    double rmsQxyz = 0.0;
    double maxAbs = 0.0;
    int used = 0;
};

// hkl candidates (excluding 000)
std::vector<Hkl> hklGrid(int maxH, int maxK, int maxL);

// Forward model: given the reciprocal metric and nu
void forwardNu(const std::array<double, 6>& G, const std::array<double, 3>& nu, int h, int k,
               int l, double& qxy, double& qz);

// Nearest-neighbour assignment; cutoff is used for early pruning in random trials
// maxG: consider only candidate indices with |G| <= maxG. If the caller knows "the maximum deviation dmax of the current assignment",
// pass qmax + k*dmax (points with larger |G| cannot be closer than the existing assignment), which shrinks the candidate table from
// the whole hkl grid (e.g. 13³=2197) to a small fraction of it; the default 1e30 means no filtering.
double assignHkl(const Cell& cell, const Orient& orient, const std::vector<Peak>& peaks,
                 const std::vector<Hkl>& cands, std::vector<Hkl>& outHkl,
                 std::vector<double>* outDist = nullptr, double cutoff = 1e30,
                 double weightQz = 1.0, double maxG = 1e30);

// LM refinement of the 6 unit-cell parameters (parameters with freeMask=false stay unchanged)
// When sys != Triclinic, enable true reparameterisation: one free variable drives the dependent cell parameters of that
// crystal system (Cubic: a->(a,b,c); Tetragonal/Hex: a->(a,b), etc.), so the objective matches the constraints.
// Triclinic (default) is line-for-line equivalent to the old implementation.
bool refineCell(const std::vector<Peak>& peaks, const std::vector<Hkl>& hkl, const Orient& orient,
                Cell& cell, const std::array<bool, 6>& freeMask, int maxIter = 80,
                CrystalSystem sys = CrystalSystem::Triclinic);

// Alternate "assign + refine" until convergence
Cell refineAssign(const std::vector<Peak>& peaks, const Orient& orient, const std::vector<Hkl>& cands,
                  Cell cell, CrystalSystem sys, int iters, double& rmsOut,
                  double maxG = 1e30);

// Estimate the |G| upper bound for assignment from the current (projected) unit cell and the assigned hkl:
// indices with |G| > qmax + 3·dmax + 0.05 cannot become candidates closer than the existing assignment.
// Returns 1e30 (no filtering) when there is no valid assignment.
double assignMaxG(const Cell& cell, const Orient& orient, const std::vector<Peak>& peaks,
                  const std::vector<Hkl>& hkl);

// Crystal-system-constrained "rescue refinement" (used by WithSpec): the candidate cell is obtained analytically from a triclinic
// parameterisation, and after projection to the target crystal system the existing (h,k,l) assignment may lock LM into an alias minimum.
// Re-solve within the constrained subspace using refineAssign's graduated assignment; replace cell / hkl only if the joint RMS strictly decreases.
// Return false directly when rmsNow <= rmsGood (solutions already within tolerance need no rescue, avoiding 5–9× overhead);
// likewise return directly when rmsNow > rmsWorst (solutions far outside tolerance cannot be recovered by 6 LM steps,
// so skipping saves the vast majority of futile rescues). Do nothing when sys == Triclinic (return false).
bool refineConstrained(const std::vector<Peak>& peaks, const Orient& orient,
                       const std::vector<Hkl>& cands, CrystalSystem sys, Cell& cell,
                       std::vector<Hkl>& hkl, double maxG, double rmsGood,
                       double rmsWorst = 1e30);

// Deviation statistics
Quality evaluate(const Cell& cell, const Orient& orient, const std::vector<Peak>& peaks,
                 const std::vector<Hkl>& hkl);

} // namespace gid
