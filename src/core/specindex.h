// ============================================================================
//  specindex.h — two-step indexing with specular peaks (known q_spec)
//
//  Implementations:
//    * improved (specindex_improved.cpp, from index2_other): Part ONE solves (a,b,γ),
//      Part TWO solves (κ,λ,ρ) by l-triple enumeration / difference-vector voting / RANSAC / interval scan.
//  The original GIDInd part1/part2 port (Original path) was removed; `Original` is simply forwarded to the improved version.
//  The current scrm (scrm.h/.cpp) is a standalone implementation: step one enumerates only in-plane (u,v), without w;
//  step two analytically recovers (κ,λ,ρ,w,l) from the specular constraint; it does not follow this file's Part ONE/Part TWO structure,
//  and plugs into the dispatcher below as Part1Algorithm::Scrm (gid::scrm::runSpecIndex can also be called explicitly).
//  The default algorithm remains PeriodVoteDdl; Scrm is only a newly added option.
//
//  Outputs Solution; both implementations include full unit-cell synthesis + LM refinement + Niggli/Buerger reduction.
// ============================================================================
#pragma once

#include "types.h"

namespace gid {

// Which implementation to use
enum class Part1Algorithm {
    Original = 0, // original GIDInd port (deprecated: path removed, now equivalent to forwarding to the improved version)
    Improved = 1, // improved version (two-step method from index2_other: l-triple enumeration)
    PeriodVoteDdl = 2,   // PeriodVote/DDL (Δ-Direct Lattice period voting): difference vectors (Δg_z of peak pairs with the same (h,k)) + period voting to fix ρ directly, then solve (κ,λ) (default)
    Ransac = 3,   // 2-peak minimal sample (specular row + two peaks = 3 equations) + per-peak early exit = deterministic RANSAC
    Ills = 4,     // integer least squares: ρ from difference-vector lattice clustering + exact t by interval scan (no bins, no top-K truncation)
    Scrm = 5      // SCRM (step one without w): in-plane (u,v) solves (a,b,γ), step two analytically recovers (κ,λ,ρ,w,l) from the specular constraint
};

// A candidate from the original Part ONE (used to bridge the original part2.cpp)
struct AbgCand {
    int u = 0, v = 0;
    double a = 0, b = 0, ga = 0;
    double za = 0, zb = 0, C = 0;   // C = z_a z_b cos(gamma)
    double guv = 0, zeta = 0;
    double rmsQxy = 0;
    std::vector<std::vector<Hkl>> hk; // [peak][R], sorted by |dqxy| ascending
};

struct SpecParams {
    double qspec = 0.0;
    // Default PeriodVote/DDL (Δ-Direct Lattice period voting): essentially independent of the l/w range; per candidate it is two orders of magnitude cheaper than enumeration
    // (DIP 12 peaks: 0.4 s with default parameters vs 1.4 s for the improved version).
    Part1Algorithm part1Alg = Part1Algorithm::PeriodVoteDdl;

    // Crystal-system constraint: default Triclinic = unconstrained (exactly the old behavior). With another crystal system,
    // every cell synthesized by Part TWO is first projected onto that system, LM refinement fits only the free
    // parameters of that system (systemFreeMask), and the 3D Niggli criterion and whole-cell reduction are skipped — reduction changes basis vectors and
    // would rewrite a γ=120° hexagonal cell or an α=γ=90° monoclinic cell into an equivalent form that does not follow the standard axis setting.
    // Note: the search phase (Part ONE / step-two analytic recovery) is still triclinic parameterization; the crystal-system constraint takes effect at the refinement stage,
    // and all final output solutions strictly satisfy the chosen crystal system.
    CrystalSystem system = CrystalSystem::Triclinic;

    // Search scale
    int redundant = 4;      // R: number of (h,k) candidates kept per peak
    int uvLimit = 2;        // (u,v) ∈ [-L, L]
    int hkLimit1 = 3;       // Part ONE sextuple range
    int hkLimit2 = 6;       // Part ONE subsequent indexing (h,k) range
    int lLimit = 6;         // l range for Part TWO
    int wLimit = 2;         // w ∈ [-L, L]
    int linesToPermute = 4; // Part ONE initial three reflections taken from the first n peaks

    // Part ONE speed-up: for the **first** of the three initial reflections (the one with the shortest indices), (h,k) enumerate only [-n,n]²
    // ("nearest-neighbor small box"); if the whole round finds nothing, runSpecIndex automatically reruns once with the full box.
    // 0 = off (full box, exactly as before); with 1 the inner iteration count drops 5.4× (7²→3²).
    int seed3Half = 0;

    int maxPart1 = 200;      // number of (a,b,gamma) candidates passed to Part TWO (used for truncation by the improved version)
    int refineTop = 200;    // number of candidates to refine (used by the improved version; each candidate is also quantized/deduplicated internally and capped at 16)
    int maxSolutions = 50;
    bool earlyStop = true;  // improved version: end the Part TWO candidate loop early once a sufficiently good solution exists

    // Fixed contact plane (known uvw): search only this orientation, avoiding enumeration over (u,v) and w
    bool fixedUvw = false;
    int fu = 0, fv = 0, fw = 1;

    // Criteria/thresholds
    bool niggli = true;
    bool restrictCell = false;
    double dqxyCutoff = 0.05;
    double dqzCutoff = 0.05;
    double dqxyzCutoff = 0.05;
    double dqspec = 0.01;

    double aLo = 3, aHi = 30, bLo = 3, bHi = 30, cLo = 3, cHi = 30;
    double alLo = 60, alHi = 120, beLo = 60, beHi = 120, gaLo = 60, gaHi = 120;
    double volLo = 50, volHi = 100000;

    int threads = -1;
};

// Dispatcher: Original/Improved/PeriodVoteDdl/Ransac/Ills forward to the improved implementation;
// Scrm forwards to the standalone implementation gid::scrm::runSpecIndex (parameters mapped via scrm::fromSpecParams).
// (The original GIDInd part1/part2 port path and the old SCRM path have both been removed, and gid::runPartTwo was removed as well).
std::vector<Solution> runSpecIndex(const std::vector<Peak>& peaks, const SpecParams& prm,
                                   const Progress& prog);

// Improved two-step method (specindex_improved.cpp, from index2_other/src/core)
namespace improved {
std::vector<AbgCand> runPartOne(const std::vector<Peak>& peaks, const SpecParams& prm,
                                const Progress& prog);
std::vector<Solution> runPartTwo(const std::vector<Peak>& peaks, const AbgCand& cand,
                                 const SpecParams& prm);
std::vector<Solution> runSpecIndex(const std::vector<Peak>& peaks, const SpecParams& prm,
                                   const Progress& prog);
} // namespace improved

} // namespace gid
