// ============================================================================
//  SimpleProbe.cpp — test probe for the fifth mode "simple indexing (out-of-plane 001)"
//
//  Usage:
//    GIDSimpleProbe <peaks.csv>                 # run simple indexing only
//    GIDSimpleProbe <peaks.csv> --compare-spec  # also run the q_spec two-step method for comparison
//  Optional arguments:
//    --qspec <v>     q_spec used by the q_spec mode (default: auto-detect)
//    --tol <v>       full-model 2D match tolerance            (default 0.05)
//    --tolXy <v>     in-plane q_xy tolerance                  (default 0.05)
//    --tolZ <v>      out-of-plane q_z tolerance               (default 0.05)
//    --rod <v>       specular rod peak q_xy threshold         (default 0.04)
//    --qmin <v>      ignore peaks near the origin             (default 0.05)
//    --hk <n>        in-plane index box |h|,|k| <= n          (default 2)
//    --l <n>         out-of-plane order box |l| <= n          (default 6)
//    --top <n>       number of solutions to output            (default 10)
//    --grid <n>      tilt Hough grid side length              (default 33)
//    --plane <n>     number of in-plane hypotheses kept       (default 12)
//    --tilt <n>      tilt candidates per (in-plane,R)         (default 4)
//    --icp <n>       reassignment rounds                      (default 4)
//    --show <n>      number of solutions with per-peak assignment printed (default 1)
//
//  Peak file: each line "q_xy q_z" (comma/semicolon/tab also accepted), '#' starts a comment.
//   Default path: D:\code\giwaxs-examples\peak_finder\aspirin.csv
// ============================================================================
#include "../src/PathIO.hpp"
#include "engine.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

using namespace gid;

static double msSince(const std::chrono::steady_clock::time_point& t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
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
        // UTF-8 BOM (common in CSV exported from Excel): feeding it straight to >> double makes the first line fail to parse
        if (line.size() >= 3 && (unsigned char)line[0] == 0xEF && (unsigned char)line[1] == 0xBB &&
            (unsigned char)line[2] == 0xBF)
            line.erase(0, 3);
        if (!line.empty() && line[0] == '#') continue;
        for (char& c : line)
            if (c == ',' || c == ';' || c == '\t') c = ' ';
        std::istringstream is(line);
        double a = 0.0, b = 0.0;
        if (is >> a >> b) out.push_back(Peak{a, b});
    }
    normalizeQxySigns(out);
    return out;
}

static int nInlier(const Solution& s) {
    int n = 0;
    for (const Hkl& x : s.hkl)
        if (x.h != 0 || x.k != 0 || x.l != 0) ++n;
    return n;
}

static void printSolution(const Solution& s, size_t i, int nPeak) {
    std::printf("#%-2zu uvw=(%2d %2d %2d) %s V=%.1f qspec=%.4f nu=(%.3f,%.3f,%.3f)\n"
                "      rms=%.4f xy=%.4f z=%.4f  inliers=%d/%d\n",
                i, s.uvw.h, s.uvw.k, s.uvw.l, fmtCell(s.cell).c_str(), s.vol, s.qspec, s.nu[0],
                s.nu[1], s.nu[2], s.rmsQxyz, s.rmsQxy, s.rmsQz, nInlier(s), nPeak);
}

static void printPeaks(const std::vector<Peak>& pk, const Solution& s) {
    std::printf("      %-22s %-16s %-22s %s\n", "peak(q_xy,q_z)", "hkl", "model(q_xy,q_z)", "d");
    for (size_t i = 0; i < pk.size(); ++i) {
        const Hkl& x = s.hkl[i];
        const QPoint q = predictQ(s, x.h, x.k, x.l);
        const double d = std::hypot(q.qxy - pk[i].qxy, q.qz - pk[i].qz);
        char hb[32];
        std::snprintf(hb, sizeof(hb), "(%2d %2d %2d)", x.h, x.k, x.l);
        std::printf("      (%7.4f,%7.4f)  %-16s (%7.4f,%8.4f)  %7.4f%s\n", pk[i].qxy, pk[i].qz, hb,
                    q.qxy, q.qz, d, (x.h == 0 && x.k == 0 && x.l == 0) ? "  [unassigned]" : "");
    }
}

int main(int argc, char** argv) {
    std::string path = "D:\\code\\giwaxs-examples\\peak_finder\\aspirin.csv";
    bool compareSpec = false;
    bool haveQspec = false;
    double qspecIn = 0.0;
    int show = 1;
    int specShow = 2;

    SimpleParams sp;
    SpecParams specPrm;

    for (int i = 1; i < argc; ++i) {
        const char* key = argv[i];
        if (!std::strcmp(key, "--compare-spec")) {
            compareSpec = true;
            continue;
        }
        if (key[0] != '-') {
            path = key;
            continue;
        }
        if (i + 1 >= argc) break;
        const char* val = argv[i + 1];
        const double v = std::atof(val);
        ++i;
        if (!std::strcmp(key, "--qspec")) {
            qspecIn = v;
            haveQspec = v > 0.0;
        } else if (!std::strcmp(key, "--tol")) sp.tol = v;
        else if (!std::strcmp(key, "--tolXy")) sp.tolXy = v;
        else if (!std::strcmp(key, "--tolZ")) sp.tolZ = v;
        else if (!std::strcmp(key, "--rod")) sp.rodTol = v;
        else if (!std::strcmp(key, "--qmin")) sp.qMin = v;
        else if (!std::strcmp(key, "--hk")) sp.hkLimit = (int)v;
        else if (!std::strcmp(key, "--l")) sp.lLimit = (int)v;
        else if (!std::strcmp(key, "--top")) sp.topK = (int)v;
        else if (!std::strcmp(key, "--grid")) sp.grid = (int)v;
        else if (!std::strcmp(key, "--plane")) sp.maxPlanes = (int)v;
        else if (!std::strcmp(key, "--tilt")) sp.maxTilts = (int)v;
        else if (!std::strcmp(key, "--icp")) sp.icpIters = (int)v;
        else if (!std::strcmp(key, "--show")) show = (int)v;
        else std::printf("unknown option %s\n", key);
    }

    const std::vector<Peak> pk = loadPeaks(path.c_str());
    if (pk.empty()) return 1;

    double qsDet = 0.0;
    int rods = 0;
    detectQSpec(pk, qsDet, rods);
    const RodPeriod rp = detectRodPeriod(pk, sp.rodTol, sp.lLimit);
    std::printf("file=%s  peaks=%zu  detectQSpec=%.4f rods=%d\n", path.c_str(), pk.size(), qsDet,
                rods);
    std::printf("rod period: ok=%d p=%.4f n=%d ambiguous=%d qz=[", rp.ok ? 1 : 0, rp.p, rp.nRod,
                rp.ambiguous ? 1 : 0);
    for (size_t i = 0; i < rp.qz.size(); ++i)
        std::printf("%s%.4f", i ? "," : "", rp.qz[i]);
    std::printf("]\n");
    std::printf("simple prm: tol=%.4f tolXy=%.4f tolZ=%.4f rod=%.4f qmin=%.3f hk=%d l=%d "
                "top=%d grid=%d plane=%d tilt=%d icp=%d\n",
                sp.tol, sp.tolXy, sp.tolZ, sp.rodTol, sp.qMin, sp.hkLimit, sp.lLimit, sp.topK,
                sp.grid, sp.maxPlanes, sp.maxTilts, sp.icpIters);
    std::fflush(stdout);

    EngineParams ep;
    ep.mode = IndexMode::Simple;
    ep.simple = sp;

    Progress prog;
    const auto t0 = std::chrono::steady_clock::now();
    prog.report = [t0](double v, const std::string& msg) {
        std::printf("  [%3.0f%% %6.0f ms] %s\n", v * 100.0, msSince(t0), msg.c_str());
        std::fflush(stdout);
    };
    const EngineResult res = runIndexing(pk, ep, prog);
    const double ms = msSince(t0);

    std::printf("\n=== Simple: %zu solutions in %.0f ms (%s) ===\n", res.solutions.size(), ms,
                res.message.c_str());
    for (size_t i = 0; i < res.solutions.size(); ++i) printSolution(res.solutions[i], i, (int)pk.size());

    if (!res.solutions.empty() && show > 0) {
        std::printf("\n-- per-peak assignment of the best solution --\n");
        printPeaks(pk, res.solutions.front());
    }

    if (compareSpec) {
        SpecParams sPrm;
        sPrm.qspec = haveQspec ? qspecIn : qsDet;
        sPrm.threads = -1;
        sPrm.maxSolutions = 20;
        if (!(sPrm.qspec > 0.0)) {
            std::printf("\n=== WithSpec: no q_spec, skipped (specify one with --qspec) ===\n");
        } else {
            EngineParams sep;
            sep.mode = IndexMode::WithSpec;
            sep.spec = sPrm;
            const auto t1 = std::chrono::steady_clock::now();
            const EngineResult r2 = runIndexing(pk, sep, prog);
            const double ms2 = msSince(t1);
            std::printf("\n=== WithSpec: %zu solutions in %.0f ms (qspec=%.4f) ===\n",
                        r2.solutions.size(), ms2, sPrm.qspec);
            const int lim = std::min<int>((int)r2.solutions.size(), specShow);
            for (int i = 0; i < lim; ++i) printSolution(r2.solutions[i], i, (int)pk.size());
            if (!r2.solutions.empty()) {
                std::printf("\n-- per-peak assignment of the WithSpec best solution --\n");
                printPeaks(pk, r2.solutions.front());
            }
        }
    }
    return 0;
}
