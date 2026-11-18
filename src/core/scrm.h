// ============================================================================
//  scrm.h — SCRM two-step method (Specular-Constrained Reciprocal Metric, step one without w)
//
//  Corresponds to docs/SCRM_specular_two_step_design_wiki.md:
//    Step one: enumerate only in-plane directions (u,v), solve the 3 in-plane unknowns (X,Y,Zc) and recover
//            (A,B,C) = (z_a², z_b², z_a z_b cosγ); w does not appear, and l is eliminated by the measured
//            q_z (out-of-plane information enters only as t = q_z/q_s).
//    Step two: solve (κ,λ,ρ) from the specular constraint nᵀt = q_spec (⇔ S3 = ±ζ) and three seed peaks,
//            then analytically recover w and the per-peak l.
//
//  This implementation is standalone: it does not call specindex_improved's seed enumeration, six-fold (h,k) loops,
//  voting/DDL or maxPart1 truncation, nor does it share data structures with it; it reuses only low-level utilities
//  (linear solver, hklGrid, refineCell/assignHkl, Image2D, peak-list types).
// ============================================================================
#pragma once

#include "specindex.h"

namespace gid {
namespace scrm {

// ---------------------------------------------------------------------------
//  Parameters (self-contained; does not change the default behavior of the shared SpecParams)
// ---------------------------------------------------------------------------
struct Params {
    double qspec = 0.0;          // specular q (required; supplied by detectQSpec or the caller)
    CrystalSystem system = CrystalSystem::Triclinic; // crystal-system constraint (see SpecParams::system)

    // Step one: in-plane direction enumeration and index range
    int uvLimit = 2;             // (u,v) ∈ [-uvLimit, uvLimit]², only half-plane representatives are used:
                                 // u ≥ 0; when u == 0 also require v ≥ 0 ((u,v) and (−u,−v) give
                                 // the same (A,B,C) with (h,k) negated as a whole). Deterministic rule, independent of traversal
                                 // order; uvLimit=2 gives 13 directions per round. Note this is deduplication at the enumeration level,
                                 // not the same as negating (κ,λ,ρ) as a whole in step two when ρ<0;
                                 // (u,v,w) in the output is normalized to a single sign representative by canonicalizeUvw (see wiki §4).
    bool fixedUvw = false;       // fixed contact plane (fu,fv,fw): solve only this direction
    int fu = 0, fv = 0, fw = 1;
    int hkLimit1 = 3;            // L1: step-one seed (h,k) radius
    int hkLimit2 = 6;            // L2: (h,k) assignment radius for all peaks
    int seed3Half = 0;           // if >0, shrink the h,k radius of the first seed to this value
    int linesToPermute = 4;      // initial seed peak pool size
    int redundant = 4;           // R: number of (h,k) groups kept per peak
    int maxPart1 = 200;          // upper bound on step-one candidates (entering step two)
    int seedL = 1;               // step-two seed l radius (±L_seed)

    // Step two: search budget and ranges for the analytic inversion
    int lLimit = 6;              // accepted range for per-peak l: |l| ≤ lLimit
    int wLimit = 3;              // the data-derived w must satisfy |w| ≤ wLimit
                                 // (note: SpecParams/UI default is 2, dual track; see wiki §5.2 for the measured comparison)
    int maxTrials = 0;           // step-two 3×3 trial limit; 0 = unlimited
    int refineTop = 200;         // number of candidates entering final refinement
    int maxSolutions = 50;       // output solution limit
    bool earlyStop = true;       // stop early once a sufficiently good solution exists

    // Thresholds / feasible region
    double dqxyCutoff = 0.05;    // step-one RMSD(q_xy) threshold multiplier
    double dqzCutoff = 0.05;     // q_z residual threshold (×3 used for early exit)
    double dqxyzCutoff = 0.05;   // full-q residual threshold
    double dqspec = 0.01;        // specular-q self-consistency threshold
    bool restrictCell = false;
    double aLo = 3, aHi = 30, bLo = 3, bHi = 30, cLo = 3, cHi = 30;
    double alLo = 60, alHi = 120, beLo = 60, beHi = 120;
    double gaLo = 60, gaHi = 120;
    double volLo = 50, volHi = 100000;
    bool niggli = true;
    int threads = -1;
};

// ---------------------------------------------------------------------------
//  Diagnostics (runtime information beyond the solution set)
// ---------------------------------------------------------------------------
// Counting convention: the integer counters below are all **multi-round accumulations** of runSpecIndex (the seedL fallback round
// and the wide-threshold fallback round each run once), not per-round values. For example, with uvLimit=2 each round has 13 (u,v) directions (half-plane representatives,
// duplicates with the same solution as (−u,−v) already removed); uvDirs=26 in test.csv = 13×2 (the wide-threshold fallback ran step one twice); normalCalls counts candidate **invocations**,
// and the same candidate is counted again in the seedL 1→3 fallback round (DIP 200 candidates → 400).
struct Stats {
    double seconds = 0.0;        // total elapsed time (filled by runSpecIndex)
    int peaksIn = 0;             // number of input peaks
    int peaksUsed = 0;           // number of peaks used in the in-plane fit (excluding q_xy≈0)
    int uvDirs = 0;              // sum of (u,v) direction counts over all step-one enumeration runs (multi-round); each round is
                                 // the number of half-plane representatives (uvLimit=2 → 13), excluding opposite-normal duplicates
    long long ipCands = 0;       // cumulative in-plane candidates passing the gate over all step-one rounds
    long long normalCalls = 0;   // cumulative number of candidates entering step two (fallback rounds double-count)
    long long secondTrials = 0;  // cumulative step-two 3×3 Cramer solves (no longer ×2 after merging ±S3)
    long long rawSol = 0;        // cumulative raw step-two solutions (including duplicates across rounds/candidates, not deduplicated)
    long long finalSol = 0;      // number of final output solutions
    int seedLUsed = 0;           // step-two seedL actually used (maximum over rounds)
    bool wideFallback = false;   // whether the wide-threshold fallback was triggered
};

// Step-one candidate: w-free in-plane solution
struct InPlaneCand {
    int u = 0, v = 0;
    double A = 0.0, B = 0.0, C = 0.0; // A = z_a², B = z_b², C = z_a z_b cosγ
    double a = 0.0, b = 0.0, ga = 90.0;
    double guv = 0.0, zeta = 0.0;     // guv = u²A+v²B−2uvC, zeta = sqrt(q_s²−guv)
    double rmsQxy = 0.0;
    std::vector<std::vector<Hkl>> hk; // [peak][R]; l is obtained analytically in step two
};

// Step one: enumerate only (u,v); returns candidates sorted by in-plane parallelogram area ascending
std::vector<InPlaneCand> runInPlane(const std::vector<Peak>& peaks, const Params& prm,
                                    const Progress& prog, Stats* st = nullptr);

// Step two: given one in-plane candidate, analytically recover (κ,λ,ρ,w,l) → unit cell → refinement/reduction
std::vector<Solution> runNormal(const std::vector<Peak>& peaks, const InPlaneCand& cand,
                                const Params& prm, Stats* st = nullptr);

// Full pipeline (step one + step two + candidate-level parallelism/early stop/fallback); standalone entry point
std::vector<Solution> runSpecIndex(const std::vector<Peak>& peaks, const Params& prm,
                                   const Progress& prog, Stats* st = nullptr);

// Maps the shared SpecParams (UI/dispatcher parameters) to this implementation's Params.
// Fields absent from SpecParams (e.g. the step-two seed radius seedL) keep this file's defaults.
Params fromSpecParams(const SpecParams& sp);

} // namespace scrm
} // namespace gid
