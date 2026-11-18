// ============================================================================
//  nospec_twosteps.h — Fourth indexing mode: no q_spec (two-step method)
//
//  See docs/NoSpec_no_specular_peak_two_step_design_algorithm_wiki.md for the derivation.
//  The only model difference from the specular-peak two-step method (SCRM): q_s = sqrt(n^T G n) changes from "observed" to
//  "quantity to estimate", so step one's continuous unknowns go from 3 (X,Y,Zc) to 4 (X,Y,Zc,q_s),
//  and the minimum peak count goes 3 → 4; step-two equations (28)–(34) are unchanged, only the measured q_s is replaced by the estimate τ.
//
//  Implementation strategy (standalone file, no existing algorithm implementation modified):
//    * Step one: outer 1-D search over τ = q_s. For a given τ, the row vector (22) depends only on τ,
//      so call the public interface of scrm::runInPlane directly, scoring with the in-plane residual R(τ);
//      first a coarse scan, then golden-section refinement of each local minimum, yielding a few τ candidates.
//    * Step two: call the public interface of scrm::runNormal on the surviving τ ((28)–(34) reused verbatim).
//      Setting scrm::Params::dqspec to a huge value disables the two hard self-consistency checks for "q_s is observed"
//      (scrm.cpp:834 / :1143), without adding a switch to SCRM or modifying its source.
//    * Output: q_spec is computed from the solution itself (Solution::qspec), dqspec is set to 0 (no measured value to compare).
//
//  Parameters: the search ranges, criteria, and reduction switches of steps one/two are all declared by this file's own Params,
//  visible and tunable item by item in the UI parameter page and probes; the .cpp performs a single call adaptation to scrm::Params,
//  embedding no scrm::Params and reusing none of its default configuration.
//
//  Relation to the one-step method nospec.cpp: mutually independent, no mutual modification; the one-step method solves 9 continuous unknowns +
//  beam search, while this method fixes the normal to integers (u,v,w), replacing the continuous dimensions with discrete enumeration.
// ============================================================================
#pragma once

// Depends only on the common types and the crystal-system enum; scrm::Params does not appear in this header (see the mapping function in the .cpp).
#include "specindex.h"
#include "types.h"

#include <vector>

namespace gid {
namespace nospec2 {

// ---------------------------------------------------------------------------
//  Parameters
// ---------------------------------------------------------------------------
// This struct is the fourth mode's **own** parameter set: it neither embeds nor inherits scrm::Params.
// Every search range, criterion, and reduction switch of steps one/two is explicitly declared here; the UI parameter page and probes
// correspond one to one; the .cpp performs a single call adaptation to scrm::Params (mapping function) and reuses none of its defaults.
struct Params {
    // ---- τ = q_spec 1-D search ----
    // <= 0 means infer automatically from the peak table:
    //   lo = max(0.5·min(q_z>0), 0.15·max|q|), hi = max(2.5·max(q_z), 1.2·max|q|).
    double qsLo = 0.0;
    double qsHi = 0.0;
    int qsSteps = 12;      // coarse scan grid points
    int qsRefine = 6;      // golden-section iterations per local minimum (0 = grid points only)
    int maxTau = 8;        // number of τ entering the full step two (§4.5.2 multiple roots need integer-consistency adjudication)

    // Rod peak soft prior (design draft §8 A5 / §7.5.4): peaks with q_xy≈0 do not participate in fitting, but q_z,rod/m are
    // natural candidate positions for q_s; they serve only as τ candidates + score discounts, not as a measured q_spec, nor as hard constraints.
    // rodQxyTol <= 0 disables this completely.
    double rodQxyTol = 0.02;
    int rodOrders = 3;
    double rodPrior = 0.5;   // <1 = prioritise entering step two

    // Lightweight budget for the τ coarse scan (also visible in the UI; only used for R(τ) scoring of τ, not for the final solve).
    // The coarse scan evaluates dozens of τ; expanding the full search budget would slow it down significantly (measured 3.4 s → 8.6 s,
    // and all five items expanded can reach 25 s).
    // To make the coarse scan fully consistent with the final stage, set these 5 items equal to the step-one main parameters below.
    int coarseHk1 = 2;
    int coarseHk2 = 4;
    int coarseSeedPeaks = 3;
    int coarseRedundant = 2;
    int coarseMaxPart1 = 40;

    // ---- Crystal system and step one: in-plane direction / index enumeration ----
    CrystalSystem system = CrystalSystem::Triclinic;
    int uvLimit = 2;         // representative half-plane of (u,v) ∈ [-uvLimit, uvLimit]²
    int hk1 = 3;             // step-one seed (h,k) radius
    int hk2 = 6;             // all-peak (h,k) assignment radius
    int seed3Half = 0;       // >0 shrinks the first seed's (h,k) radius to this value
    int seedPeaks = 4;       // initial seed peak pool size
    int redundant = 4;       // number of (h,k) sets kept per peak
    int maxPart1 = 200;      // step-one candidate cap
    double dqxyCutoff = 0.05;

    // ---- Step two: out-of-plane integer enumeration ----
    int seedL = 1;           // seed l radius (±seedL)
    int lLimit = 6;          // per-peak acceptance |l| ≤ lLimit
    int wLimit = 3;          // |w| ≤ wLimit derived from the data
    int maxTrials = 0;       // step-two 3×3 trial cap; 0 = unlimited
    int refineTop = 200;     // number of candidates entering final refinement
    int seedLFallback = 3;   // widen seedL to this value when the first τ yields no good solution; 0 = disable fallback
    double dqzCutoff = 0.05;

    // ---- Unit-cell criteria and reduction ----
    bool restrictCell = false;
    double aLo = 3, aHi = 30, bLo = 3, bHi = 30, cLo = 3, cHi = 30;
    double alLo = 60, alHi = 120, beLo = 60, beHi = 120, gaLo = 60, gaHi = 120;
    double volLo = 50, volHi = 100000;
    bool niggli = true;

    // ---- Solve scale / output ----
    int maxIpCands = 16;     // in-plane candidate cap per τ sent to step two (best by rmsQxy)
    int maxSolutions = 50;
    double dqxyzCutoff = 0.05;
    int threads = -1;        // -1 = automatic; coarse scan is fixed to single thread

    Params sanitized() const;
};

// ---------------------------------------------------------------------------
//  Diagnostics
// ---------------------------------------------------------------------------
struct Stats {
    double seconds = 0.0;
    int peaksIn = 0;          // input peak count
    int peaksUsed = 0;        // peaks participating in fitting (q_xy > 0)
    double qsLoActual = 0.0;  // actually used τ search interval (0 means automatic inference)
    double qsHiActual = 0.0;
    int tauTried = 0;         // number of τ evaluated in step one (coarse scan + refinement)
    int tauKept = 0;          // number of τ entering step two
    std::vector<double> tauUsed;    // list of τ entering step two
    std::vector<double> tauBestRms; // minimum rmsQxyz obtained by each τ in step two (kBig = no solution)
    long long ipCands = 0;    // cumulative in-plane candidates produced by step one
    long long secondTrials = 0;
    long long rawSol = 0;
    int finalSol = 0;
    double bestTau = 0.0;       // best-scoring τ (diagnostic only, not used for standalone adjudication)
    double bestTauScore = 0.0;  // corresponding τ score (≈ best rmsQxyz + 0.05·in-plane residual)
};

// Full pipeline. The returned Solution::hkl covers only the peaks participating in the search (same as SCRM),
// and is expanded to the full input peak table by engine::runIndexing or the caller using alignSolution.
std::vector<Solution> runNoSpecTwoSteps(const std::vector<Peak>& peaks, const Params& prm,
                                        const Progress& prog, Stats* st = nullptr);

} // namespace nospec2
} // namespace gid
