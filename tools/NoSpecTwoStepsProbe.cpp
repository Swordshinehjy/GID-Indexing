// ============================================================================
//  NoSpecTwoStepsProbe.cpp — independent probe for the fourth method "no q_spec (two-step)"
//
//  Default data: DIP 12 peaks (including the {0.000, 1.776} rod peak). The rod peak does not join the
//  fit and is not treated as a measured q_spec; its q_z only serves as a soft prior candidate for τ = q_s (design note §8 A5 / §7.5.4).
//  --no-rod removes the rod peak to test the global τ search on pure 11 peaks with no prior at all.
//  Criterion: whether the manual cell a=7.13 b=8.48 c=16.67 / 89.4 87.8 89.7 is reproduced.
//  Also the self-test entry of this method: returns 0 when the manual cell is hit, otherwise 1.
//
//  Usage:
//    GIDNoSpec2.exe [--peaks FILE] [--lo X] [--hi X] [--steps N] [--maxTau N]
//                   [--ip N] [--maxSol N] [--dxyz X] [--threads N]
//                   [--uv N] [--hk1 N] [--hk2 N] [--redundant N] [--max1 N]
//                   [--seedL N] [--l N] [--w N] [--refineTop N] ...
//                   [--engine] [--no-rod] [--no-ref] [--help]
//    --peaks FILE  each line "q_xy q_z", '#' starts a comment; built-in DIP 12 peaks if omitted
//    --engine      go through engine::runIndexing (verify the UI/engine integration path)
//    --no-rod      drop rod peaks with q_xy≈0 (pure data without prior)
//    --no-ref      disable the manual-cell criterion for custom peak tables; only require solutions
//    Both step-one and step-two search ranges are fully overridable (one-to-one with nospec2::Params /
//    the UI parameter page); see --help for the full list.
// ============================================================================
#include "../src/PathIO.hpp"
#include "engine.h"
#include "nospec_twosteps.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace gid;

static const std::vector<Peak> kDipAll = {
    {0.000, 1.776}, {0.368, 0.091}, {0.396, 1.160}, {0.428, 1.069}, {0.689, 0.984},
    {0.734, 0.182}, {0.806, 1.232}, {1.140, 0.171}, {1.223, 0.083}, {1.529, 0.790},
    {1.625, 0.697}, {2.274, 0.953},
};
static const Cell kRef{7.13, 8.48, 16.67, 89.4, 87.8, 89.7};
static const Hkl kRefUvw{-1, 2, 1};

static std::vector<Peak> dipNoSpec() {
    std::vector<Peak> out;
    for (const Peak& p : kDipAll)
        if (p.qxy > 1e-6) out.push_back(p);
    return out;
}

static std::vector<Peak> loadPeaks(const char* path) {
    std::vector<Peak> out;
    auto f = pathio::input(path);
    if (!f) {
        std::printf("cannot open %s\n", path);
        return out;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line[0] == '#') continue;
        for (char& c : line)
            if (c == ',' || c == ';' || c == '\t') c = ' ';
        std::istringstream is(line);
        double qxy = 0.0, qz = 0.0;
        if (is >> qxy >> qz) out.push_back(Peak{qxy, qz});
    }
    return out;
}

static void usage() {
    std::printf(
        "GIDNoSpec2 -- no-q_spec two-step indexing probe\n"
        "  --peaks FILE   per line: q_xy q_z (default: built-in DIP 11 peaks)\n"
        "  --lo X --hi X  q_spec search range (0 = auto)\n"
        "  --steps N      coarse tau grid points\n"
        "  --maxTau N     tau candidates entering step two\n"
        "  --ip N         in-plane candidates per tau\n"
        "  --maxSol N     max solutions\n"
        "  --dxyz X       rmsQxyz cutoff\n"
        "  --threads N    0 = auto\n"
        "  --refine N     golden-section iterations per local minimum\n"
        "  --coarseHk1 N --coarseHk2 N --coarseSeedPeaks N --coarseRedundant N\n"
        "  --coarseMax1 N  lightweight budget used only for tau scoring\n"
        "  --rodTol X --rodOrders N --rodPrior X   rod soft prior\n"
        "  --uv N --hk1 N --hk2 N --seed3 N --seedPeaks N --redundant N --max1 N\n"
        "                 step-1 in-plane enumeration ranges\n"
        "  --dqxy X       step-1 rmsQxy cutoff\n"
        "  --seedL N --seedLFallback N --l N --w N --trials N --refineTop N\n"
        "                 step-2 out-of-plane enumeration ranges\n"
        "  --dqz X        step-2 q_z cutoff\n"
        "  --restrict N --niggli N   cell range filter / Niggli reduction (0|1)\n"
        "  --engine       go through engine::runIndexing\n"
        "  --no-rod       drop q_xy~0 rod peaks (no prior)\n"
        "  --no-ref       only require at least one solution\n");
}

int main(int argc, char** argv) {
    std::vector<Peak> pk = kDipAll; // includes the rod peak: soft τ prior only, not a measured q_spec
    nospec2::Params prm;
    bool useEngine = false, refMode = true;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto nextDouble = [&](double& x) {
            if (i + 1 < argc) x = std::atof(argv[++i]);
        };
        auto nextInt = [&](int& x) {
            if (i + 1 < argc) x = std::atoi(argv[++i]);
        };
        if (a == "--help" || a == "-h") { usage(); return 0; }
        else if (a == "--peaks" && i + 1 < argc) {
            pk = loadPeaks(argv[++i]);
            if (pk.empty()) {
                std::printf("no peaks loaded\n");
                return 1;
            }
        }
        else if (a == "--lo") nextDouble(prm.qsLo);
        else if (a == "--hi") nextDouble(prm.qsHi);
        else if (a == "--steps") nextInt(prm.qsSteps);
        else if (a == "--maxTau") nextInt(prm.maxTau);
        else if (a == "--ip") nextInt(prm.maxIpCands);
        else if (a == "--maxSol") nextInt(prm.maxSolutions);
        else if (a == "--dxyz") nextDouble(prm.dqxyzCutoff);
        else if (a == "--threads") nextInt(prm.threads);
        else if (a == "--refine") nextInt(prm.qsRefine);
        else if (a == "--coarseHk1") nextInt(prm.coarseHk1);
        else if (a == "--coarseHk2") nextInt(prm.coarseHk2);
        else if (a == "--coarseSeedPeaks") nextInt(prm.coarseSeedPeaks);
        else if (a == "--coarseRedundant") nextInt(prm.coarseRedundant);
        else if (a == "--coarseMax1") nextInt(prm.coarseMaxPart1);
        else if (a == "--rodTol") nextDouble(prm.rodQxyTol);
        else if (a == "--rodOrders") nextInt(prm.rodOrders);
        else if (a == "--rodPrior") nextDouble(prm.rodPrior);
        else if (a == "--uv") nextInt(prm.uvLimit);
        else if (a == "--hk1") nextInt(prm.hk1);
        else if (a == "--hk2") nextInt(prm.hk2);
        else if (a == "--seed3") nextInt(prm.seed3Half);
        else if (a == "--seedPeaks") nextInt(prm.seedPeaks);
        else if (a == "--redundant") nextInt(prm.redundant);
        else if (a == "--max1") nextInt(prm.maxPart1);
        else if (a == "--dqxy") nextDouble(prm.dqxyCutoff);
        else if (a == "--seedL") nextInt(prm.seedL);
        else if (a == "--seedLFallback") nextInt(prm.seedLFallback);
        else if (a == "--l") nextInt(prm.lLimit);
        else if (a == "--w") nextInt(prm.wLimit);
        else if (a == "--trials") nextInt(prm.maxTrials);
        else if (a == "--refineTop") nextInt(prm.refineTop);
        else if (a == "--dqz") nextDouble(prm.dqzCutoff);
        else if (a == "--restrict") { int t = 0; nextInt(t); prm.restrictCell = (t != 0); }
        else if (a == "--niggli") { int t = 1; nextInt(t); prm.niggli = (t != 0); }
        else if (a == "--engine") useEngine = true;
        else if (a == "--no-rod") { pk = dipNoSpec(); prm.rodQxyTol = 0.0; }
        else if (a == "--no-ref") refMode = false;
    }
    if (prm.threads == 0) prm.threads = -1;

    const Mat3 Rref = reciprocalMatrix(kRef);
    const double qsRef = qSpecOf(Rref, kRefUvw);
    const double vRef = cellVolume(kRef);

    std::printf("=== no-q_spec TWO-STEP | %zu peaks | qs range=[%.4f, %.4f] (0=auto) | qspec(ref)=%.6f ===\n",
                pk.size(), prm.qsLo, prm.qsHi, qsRef);
    std::printf("params: steps=%d refine=%d maxTau=%d ip=%d maxSol=%d dxyz=%.4f threads=%d %s\n",
                prm.qsSteps, prm.qsRefine, prm.maxTau, prm.maxIpCands, prm.maxSolutions,
                prm.dqxyzCutoff, prm.threads, useEngine ? "via engine" : "direct");

    Progress prog;
    nospec2::Stats st;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<Solution> sols;
    if (useEngine) {
        EngineParams ep;
        ep.mode = IndexMode::NoSpecTwoSteps;
        ep.nospec2 = prm;
        sols = runIndexing(pk, ep, prog).solutions;
    } else {
        sols = nospec2::runNoSpecTwoSteps(pk, prm, prog, &st);
    }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // Align to the full input peak table: the hkl returned by the core covers only searching peaks (as in SCRM).
    for (Solution& s : sols) alignSolution(pk, s);

    if (!useEngine) {
        std::printf("\nstats: peaks %d/%d | tau tried=%d kept=%d | ipCands=%lld"
                    " secondTrials=%lld rawSol=%lld final=%d"
                    " | bestTau=%.6f (score=%.5f) | %.2f s\n",
                    st.peaksUsed, st.peaksIn, st.tauTried, st.tauKept, st.ipCands,
                    st.secondTrials, st.rawSol, st.finalSol, st.bestTau, st.bestTauScore, sec);
        std::printf("tau range=[%.4f, %.4f]\n", st.qsLoActual, st.qsHiActual);
        for (size_t i = 0; i < st.tauUsed.size(); ++i) {
            const double b = (i < st.tauBestRms.size()) ? st.tauBestRms[i] : 0.0;
            std::printf("  tau[%zu]=%.5f  bestRms=%.5f\n", i, st.tauUsed[i], b);
        }
    } else {
        std::printf("\n(via engine, %.2f s)\n", sec);
    }

    std::printf("\n=== %zu solutions ===\n", sols.size());
    int firstHit = -1;
    for (size_t i = 0; i < sols.size() && i < 20; ++i) {
        const Solution& s = sols[i];
        const double dv = std::fabs(s.vol - vRef) / vRef;
        const bool hit = refMode && s.rmsQxyz < 0.01 && dv < 0.05 &&
                         std::fabs(s.qspec - qsRef) < 0.05;
        if (hit && firstHit < 0) firstHit = (int)i;
        std::printf("#%-3zu uvw=(%2d %2d %2d) %s V=%.1f qspec=%.4f rms=%.4f"
                    " xy=%.4f z=%.4f dqspec=%.3f%s\n",
                    i, s.uvw.h, s.uvw.k, s.uvw.l, fmtCell(s.cell).c_str(), s.vol, s.qspec,
                    s.rmsQxyz, s.rmsQxy, s.rmsQz, s.dqspec, hit ? "  <== HIT" : "");
    }
    if (refMode) std::printf("\nfirstHit=%d  qspecRef=%.6f  Vref=%.1f\n", firstHit, qsRef, vRef);

    const bool pass = refMode ? (firstHit >= 0) : (!sols.empty());
    std::printf("[%s] no-q_spec two-step: %s\n", pass ? "PASS" : "FAIL",
                pass ? "manual cell / solution found" : "no matching solution");
    return pass ? 0 : 1;
}
