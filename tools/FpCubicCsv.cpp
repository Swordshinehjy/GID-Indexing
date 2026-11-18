// Temporary probe: read peaks from a CSV -> fixed contact plane search (cubic, fixed plane u v w)
#include "../src/PathIO.hpp"
#include "fixedplane.h"
#include "engine.h"
#include "refine.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace gid;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: FpCubicCsv <csv> [u v w] [maxHKL] [system 0..6]\n");
        return 1;
    }
    std::vector<Peak> pk;
    {
        auto in = pathio::input(argv[1]);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            double x = 0.0, z = 0.0;
            if (std::sscanf(line.c_str(), "%lf,%lf", &x, &z) == 2) pk.push_back(Peak{x, z});
        }
    }
    // The sign of q_xy only reflects the mirror chirality of the detector geometry: after reading, uniformly take |q_xy| (same convention as the GUI/core).
    normalizeQxySigns(pk);
    FixedPlaneParams fp;
    if (argc >= 5) {
        fp.u = std::atoi(argv[2]);
        fp.v = std::atoi(argv[3]);
        fp.w = std::atoi(argv[4]);
    } else {
        fp.u = 0; fp.v = 1; fp.w = 4;
    }
    if (argc >= 6) {
        fp.maxH = fp.maxK = fp.maxL = std::atoi(argv[5]);
    }
    if (argc >= 7) fp.system = (CrystalSystem)std::atoi(argv[6]);
    if (argc >= 8) fp.beam = std::atoi(argv[7]);
    if (argc >= 9) fp.start.a = fp.start.b = fp.start.c = std::atof(argv[8]);
    // Optional: specify the k / l upper bounds separately (useful for reproducing asymmetric boxes such as 4/4/6)
    if (argc >= 10) fp.maxK = std::atoi(argv[9]);
    if (argc >= 11) fp.maxL = std::atoi(argv[10]);
    // Optional: argv[11] == "engine" uses runIndexing (same path as the UI, including alignSolution)
    const bool useEngine = argc >= 12 && (std::strcmp(argv[11], "engine") == 0 ||
                                          std::strcmp(argv[11], "engine-nosign") == 0);
    if (argc >= 12 && (std::strcmp(argv[11], "nosign") == 0 ||
                       std::strcmp(argv[11], "engine-nosign") == 0))
        fp.autoSign = false;

    std::printf("peaks=%zu uvw=(%d %d %d) max=%d,%d,%d system=%s\n", pk.size(), fp.u, fp.v, fp.w,
                fp.maxH, fp.maxK, fp.maxL, crystalSystemName(fp.system).c_str());
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<Solution> sols;
    if (useEngine) {
        EngineParams ep;
        ep.mode = IndexMode::FixedPlane;
        ep.fixed = fp;
        sols = runIndexing(pk, ep, Progress{}).solutions;
    } else {
        sols = runFixedPlane(pk, fp, Progress{});
    }
    const double ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::printf("sols=%zu  total=%.0f ms%s\n", sols.size(), ms, useEngine ? "  (engine)" : "");
    for (size_t i = 0; i < sols.size(); ++i) {
        const Solution& s = sols[i];
        std::printf("[%zu] uvw=(%d %d %d) %s V=%.1f qspec=%.4f rms=%.5f (xy=%.5f z=%.5f)\n", i,
                    s.uvw.h, s.uvw.k, s.uvw.l, fmtCell(s.cell).c_str(), s.vol, s.qspec, s.rmsQxyz,
                    s.rmsQxy, s.rmsQz);
        for (size_t p = 0; p < pk.size() && p < s.hkl.size(); ++p) {
            const QPoint q = predictQ(s, s.hkl[p].h, s.hkl[p].k, s.hkl[p].l);
            std::printf("      pk%2zu obs(%7.4f,%7.4f) -> (%2d %2d %2d) pred(%7.4f,%7.4f) "
                        "dqxy=%7.4f dqz=%7.4f\n",
                        p, pk[p].qxy, pk[p].qz, s.hkl[p].h, s.hkl[p].k, s.hkl[p].l, q.qxy, q.qz,
                        q.qxy - pk[p].qxy, q.qz - pk[p].qz);
        }
    }
    return 0;
}
