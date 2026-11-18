// ============================================================================
//  NoSpecProbe.cpp — single-step method test probe without q_spec (unknown q_spec)
//  Data: the 12 DIP peaks from DipProbe, with the {0.000, 1.776} specular rod peak removed (q_xy=0)
//  Criterion: whether the manual cell is reproduced (a=7.13 b=8.48 c=16.67, al=89.4 be=87.8 ga=89.7)
// ============================================================================
#include "../src/PathIO.hpp"
#include "engine.h"
#include "nospec.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

// Manual solution
static const Cell kRef{7.13, 8.48, 16.67, 89.4, 87.8, 89.7};
static const Hkl kRefUvw{-1, 2, 1};

static std::vector<Peak> dipNoSpec() {
    std::vector<Peak> out;
    for (const Peak& p : kDipAll)
        if (p.qxy > 1e-6) out.push_back(p);
    return out;
}

static double msSince(const std::chrono::steady_clock::time_point& t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// Each line "q_xy q_z" (comma/semicolon/tab also accepted), '#' starts a comment.
// iostream is used instead of fopen/sscanf: MSVC reports C4996 for the latter (tool projects may
// not enable _CRT_SECURE_NO_WARNINGS), with identical semantics.
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
        double a = 0.0, b = 0.0;
        if (is >> a >> b) out.push_back(Peak{a, b});
    }
    // The sign of q_xy only reflects the mirror chirality of the detector geometry: after reading, uniformly take |q_xy|.
    normalizeQxySigns(out);
    return out;
}

// Consistent with the "participates in search" criterion in nospec.cpp: q_xy>1e-6 and finite.
// Specular rod points with q_xy≈0 do not belong to the model without specular, are culled by the solver, and take no hkl slot.
static bool inSearch(const Peak& p) {
    return p.qxy > 1e-6 && std::isfinite(p.qxy) && std::isfinite(p.qz);
}

// Consistent with the deduplication criterion in nospec.cpp (a position difference <1e-6 counts as the same reciprocal point, e.g. a Friedel pair).
// Returns the index of an earlier peak at the same position, or -1; used to tag "merged into which entry" in the output.
static int duplicateOf(const std::vector<Peak>& pk, size_t i) {
    for (size_t j = 0; j < i; ++j)
        if (inSearch(pk[j]) && std::fabs(pk[j].qxy - pk[i].qxy) < 1e-6 &&
            std::fabs(pk[j].qz - pk[i].qz) < 1e-6)
            return (int)j;
    return -1;
}

int main(int argc, char** argv) {
    std::vector<Peak> pk = dipNoSpec();
    bool refMode = true; // print the manual-cell comparison only when using the built-in DIP data
    bool useEngine = false;
    NoSpecParams prm;
    prm.threads = -1;
    prm.maxSolutions = 50;
    // Allow command-line overrides (convenient for parameter tuning and comparison)
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--withspec")) {
            pk = kDipAll; // keep the {0,1.776} specular peak
            continue;
        }
        if (!std::strcmp(argv[i], "--engine")) {
            useEngine = true; // go through runIndexing (same path as the UI, including alignSolution)
            continue;
        }
        if (i + 1 >= argc) break;
        const char* key = argv[i];
        const char* val = argv[i + 1];
        const double v = std::atof(val);
        ++i;
        if (!std::strcmp(key, "--peaks")) {
            pk = loadPeaks(val);
            refMode = false;
        } else if (!std::strcmp(key, "--system")) {
            prm.system = static_cast<CrystalSystem>((int)v);
        } else if (!std::strcmp(key, "--method")) prm.method = static_cast<NoSpecMethod>((int)v);
        else if (!std::strcmp(key, "--beam")) prm.beam = (int)v;
        else if (!std::strcmp(key, "--seed")) {
            prm.seedH = prm.seedK = prm.seedL = (int)v;
        }
        else if (!std::strcmp(key, "--seedPeaks")) prm.seedPeaks = (int)v;
        else if (!std::strcmp(key, "--shell")) prm.seedShell = (int)v;
        else if (!std::strcmp(key, "--prune")) prm.metricPrune = v;
        else if (!std::strcmp(key, "--seedH")) prm.seedH = (int)v;
        else if (!std::strcmp(key, "--seedK")) prm.seedK = (int)v;
        else if (!std::strcmp(key, "--seedL")) prm.seedL = (int)v;
        else if (!std::strcmp(key, "--minRows")) prm.minRows = (int)v;
        else if (!std::strcmp(key, "--refineTop")) prm.refineTop = (int)v;
        else if (!std::strcmp(key, "--maxSol")) prm.maxSolutions = (int)v;
        else if (!std::strcmp(key, "--dzq")) prm.dqzCutoff = v;
        else if (!std::strcmp(key, "--dxyz")) prm.dqxyzCutoff = v;
        else if (!std::strcmp(key, "--grid"))
            prm.gridH = prm.gridK = prm.gridL = (int)v;
        else if (!std::strcmp(key, "--threads")) prm.threads = (int)v;
        else std::printf("unknown option %s\n", key);
    }

    double qs = 0.0;
    int rods = 0;
    detectQSpec(pk, qs, rods);
    const Mat3 Rref = reciprocalMatrix(kRef);
    const double qsRef = qSpecOf(Rref, kRefUvw);
    const double vRef = cellVolume(kRef);
    std::printf("peaks=%zu  detectQSpec=%.4f rods=%d  system=%d\n", pk.size(), qs, rods,
                (int)prm.system);
    if (refMode)
        std::printf("ref cell: %s  V=%.2f  qspec(uvw=-1,2,1)=%.4f\n", fmtCell(kRef).c_str(), vRef,
                    qsRef);
    std::printf("prm: beam=%d seed=(%d,%d,%d) seedPeaks=%d minRows=%d refineTop=%d grid=%d "
                "dqz=%.4f dqxyz=%.4f\n",
                prm.beam, prm.seedH, prm.seedK, prm.seedL, prm.seedPeaks, prm.minRows,
                prm.refineTop, prm.gridH, prm.dqzCutoff, prm.dqxyzCutoff);
    std::fflush(stdout);

    const auto t0 = std::chrono::steady_clock::now();
    Progress prog;
    prog.report = [t0](double v, const std::string& msg) {
        std::printf("  [%3.0f%% %7.0f ms] %s\n", v * 100.0, msSince(t0), msg.c_str());
        std::fflush(stdout);
    };

    const std::vector<Solution> sols = useEngine
                                          ? runIndexing(pk, EngineParams{IndexMode::NoSpec, {}, prm}, prog).solutions
                                          : runNoSpec(pk, prm, prog);
    const double ms = msSince(t0);

    std::printf("\n=== %zu solutions in %.0f ms ===\n", sols.size(), ms);
    int firstHit = -1;
    for (size_t i = 0; i < sols.size(); ++i) {
        const Solution& s = sols[i];
        const double dv = std::fabs(s.vol - vRef) / vRef;
        const bool hit = refMode && s.rmsQxyz < 0.01 && dv < 0.05 &&
                         std::fabs(s.qspec - qsRef) < 0.05;
        if (hit && firstHit < 0) firstHit = (int)i;
        std::printf("#%-3zu uvw=(%2d %2d %2d) %s V=%.1f qspec=%.4f rms=%.4f xy=%.4f z=%.4f nu=(%.3f,%.3f,%.3f)%s\n",
                    i, s.uvw.h, s.uvw.k, s.uvw.l, fmtCell(s.cell).c_str(), s.vol, s.qspec,
                    s.rmsQxyz, s.rmsQxy, s.rmsQz, s.nu[0], s.nu[1], s.nu[2], hit ? "  <== HIT" : "");
    }
    if (refMode) std::printf("\nfirstHit=%d  qspecRef=%.4f\n", firstHit, qsRef);

    // nu of the reference solution itself (to check whether the solved orientation is equivalent)
    const std::array<double, 3> nuRef = nuFromUvw(Rref, kRefUvw);

    // Residual of the manual cell under this model (theoretical lower bound without specular data); only meaningful for the built-in DIP data
    if (refMode) {
        std::printf("ref nu=(%.4f,%.4f,%.4f)\n", nuRef[0], nuRef[1], nuRef[2]);
        std::vector<Hkl> ah;
        const std::vector<Hkl> grid = hklGrid(6, 6, 6);
        const Orient o = Orient::fromNu(nuRef);
        assignHkl(kRef, o, pk, grid, ah);
        const Quality q = evaluate(kRef, o, pk, ah);
        std::printf("ref-model rmsQxyz=%.5f xy=%.5f z=%.5f\n", q.rmsQxyz, q.rmsQxy, q.rmsQz);
        std::printf("ref assignment:\n");
        for (size_t i = 0; i < pk.size(); ++i) {
            double qxy = 0.0, qz = 0.0;
            forwardNu(recipMetric(kRef), nuRef, ah[i].h, ah[i].k, ah[i].l, qxy, qz);
            std::printf("   peak(%.3f,%.3f) -> (%2d %2d %2d)  model(%.3f,%.3f) d=(%.4f,%.4f)\n",
                        pk[i].qxy, pk[i].qz, ah[i].h, ah[i].k, ah[i].l, qxy, qz, qxy - pk[i].qxy,
                        qz - pk[i].qz);
        }
    }

    // Forward comparison: a solution without q_spec has only nu; the UI/export always goes through nu (predictQ)
    if (!sols.empty()) {
        // types.h contract: Solution::hkl / dist match the input peaks in length and line up one-to-one
        // by index only after alignSolution. runNoSpec returns the pre-alignment form, its hkl covering
        // only peaks in the search (rod points with q_xy≈0 and invalid points culled; duplicates still
        // take one entry each). Reading by pk.size() index would go out of bounds and shift peaks after
        // the rod, so expand per the contract first (rod points filled in by alignSolution or left as (0,0,0)).
        Solution s0 = sols[0];
        alignSolution(pk, s0);

        // Error comparison counts only searching peaks, consistent with Solution::rms*; after alignment
        // the rod points' hkl are filled-in/placeholder values, not model predictions, and are excluded.
        std::vector<Peak> pkSearch;
        std::vector<Hkl> ahSearch;
        for (size_t i = 0; i < pk.size(); ++i) {
            if (!inSearch(pk[i])) continue;
            pkSearch.push_back(pk[i]);
            ahSearch.push_back(s0.hkl[i]);
        }
        const Quality qn = evaluate(s0.cell, Orient::fromNu(s0.nu), pkSearch, ahSearch);
        const Quality qu = evaluate(s0.cell, Orient::fromUvw(s0.uvw), pkSearch, ahSearch);
        std::printf("\npredict check on #0: uvw=(%d %d %d)\n", s0.uvw.h, s0.uvw.k, s0.uvw.l);
        std::printf("  via nu  (current UI path): rmsQxyz=%.4f (xy=%.4f z=%.4f)\n", qn.rmsQxyz,
                    qn.rmsQxy, qn.rmsQz);
        std::printf("  via uvw (old path)       : rmsQxyz=%.4f (xy=%.4f z=%.4f)\n", qu.rmsQxyz,
                    qu.rmsQxy, qu.rmsQz);
        std::printf("  model covers %zu/%zu peaks; excluded peaks are tagged below\n",
                    pkSearch.size(), pk.size());
        const std::array<double, 6> G0 = recipMetric(s0.cell);
        std::printf("  per-peak (mea -> nu-model -> uvw-model):\n");
        for (size_t i = 0; i < pk.size(); ++i) {
            double nx = 0.0, nz = 0.0;
            forwardNu(G0, s0.nu, s0.hkl[i].h, s0.hkl[i].k, s0.hkl[i].l, nx, nz);
            const QPoint u = forwardHkl(s0.cell, s0.uvw, s0.hkl[i].h, s0.hkl[i].k, s0.hkl[i].l);
            std::string tag;
            if (!inSearch(pk[i])) {
                tag = (!std::isfinite(pk[i].qxy) || !std::isfinite(pk[i].qz))
                          ? "  <invalid: excluded from search>"
                          : "  <rod q_xy≈0: excluded from search>";
            } else {
                const int d = duplicateOf(pk, i);
                if (d >= 0) tag = "  <dup of #" + std::to_string(d) + ">";
            }
            std::printf(
                "    #%2zu (%6.3f,%6.3f) hkl(%2d %2d %2d) nu(%6.3f,%6.3f) uvw(%6.3f,%6.3f)%s\n",
                i, pk[i].qxy, pk[i].qz, s0.hkl[i].h, s0.hkl[i].k, s0.hkl[i].l, nx, nz, u.qxy,
                u.qz, tag.c_str());
        }
    }
    return 0;
}
