// ============================================================================
//  engine.h — Unified entry: five indexing modes + parameter container + q_spec auto-detection
// ============================================================================
#pragma once

#include "fixedplane.h"
#include "nospec.h"
#include "nospec_twosteps.h"
#include "simpleindex.h"
#include "specindex.h"

namespace gid {

enum class IndexMode {
    WithSpec,        // two-step method with specular peaks (known q_spec)
    NoSpec,          // one-step method without q_spec
    FixedPlane,      // user-specified fixed contact plane (u v w)
    NoSpecTwoSteps,  // two-step method without q_spec (q_s as the 4th unknown in step one)
    Simple           // simple indexing (out-of-plane fixed to (001), automatic search)
};

const char* modeName(IndexMode m);

struct EngineParams {
    IndexMode mode = IndexMode::WithSpec;
    SpecParams spec;
    NoSpecParams nospec;
    FixedPlaneParams fixed;
    nospec2::Params nospec2; // fourth mode: no q_spec (two-step method)
    SimpleParams simple;     // fifth mode: simple indexing (out-of-plane 001)
};

struct EngineResult {
    bool ok = false;
    bool stopped = false;
    std::string message;
    std::vector<Solution> solutions;
};

// ----------------------------------------------------------------------------
//  Specular rod peak "fundamental frequency + order"
//
//  A q_z observed on the rod must equal n·p (p = fundamental frequency = |G(u,v,w)^T|, n ∈ Z⁺ unknown):
//  when the first-order peak is missing (very thin film, systematic extinction, weak (00l)…), the smallest observed value is 2p, 3p …
//  The old implementation returned min(q_z) directly and treated it as the modulus of (u,v,w), so anywhere q_spec is a hard
//  constraint (|q_spec(solution) − q_s| ≤ dqspec in the two-step Part ONE/TWO) implicitly asserted n = 1,
//  and the result was either total failure of the path or a subcell/supercell integer-multiplied along the normal.
//
//  Here the fundamental frequency is solved inversely as "the smallest n that explains all rod peaks": if n = 1 explains them, leave it alone;
//  only on multi-order data such as {2p, 3p} is the fundamental frequency compressed to p.
// ----------------------------------------------------------------------------
struct RodPeriod {
    bool ok = false;
    double p = 0.0;             // fundamental frequency |G(u,v,w)^T|
    int nRod = 0;               // number of rod peaks
    bool ambiguous = false;     // only one independent q_z: order undecidable, downstream must scan p/n
    std::vector<double> qz;     // deduplicated rod peak q_z (ascending)
    std::vector<int> order;     // order aligned with qz = round(qz/p)
    std::vector<int> peakIdx;   // indices of rod peaks in the input peak table
};
// rodTol: |q_xy| <= rodTol is treated as a rod peak; nmax: highest order on the rod
RodPeriod detectRodPeriod(const std::vector<Peak>& peaks, double rodTol, int nmax = 6);

// Legacy-compatible interface (UI "auto-detect" / tool probes): qspec now returns the **fundamental frequency** rather than the smallest observed q_z;
// rodCount is still the number of rod peaks. When the first-order peak exists, the return value is exactly the same as the old implementation.
bool detectQSpec(const std::vector<Peak>& peaks, double& qspec, int& rodCount);

// Align the solution's hkl to the full input peak table (points on the specular rod are filled back as m*(u,v,w))
void alignSolution(const std::vector<Peak>& peaks, Solution& s);

EngineResult runIndexing(const std::vector<Peak>& peaks, const EngineParams& prm,
                         const Progress& prog);

} // namespace gid
