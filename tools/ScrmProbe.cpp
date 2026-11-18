// ============================================================================
//  ScrmProbe.cpp — scrm self-test probe (not compiled by default, see GID_BUILD_SCRM_PROBE in CMakeLists)
//
//  Usage: GIDScrmProbe [--w N] [--top N] [--csv PATH]
//    --w N     set scrm's wLimit (default from scrm::Params, i.e. 3; UI/SpecParams default 2)
//    --top N   print the top N solutions of each implementation (default 5)
//    --csv PATH specify the CSV path (each line "q_xy,q_z")
//    1) first run the built-in DIP 12 peaks, print them next to the existing gid::runSpecIndex and check whether the manual solution is hit;
//    2) then run the CSV; if no parameters are given, try D:/code/giwaxs-examples/peak_finder/test.csv.
//
//  Discipline: read input and print only; never write any file.
// ============================================================================
#include "engine.h"
#include "refine.h"
#include "scrm.h"
#include "specindex.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace gid;

static const std::vector<Peak> kDip = {
    {0.000, 1.776}, {0.368, 0.091}, {0.396, 1.160}, {0.428, 1.069}, {0.689, 0.984},
    {0.734, 0.182}, {0.806, 1.232}, {1.140, 0.171}, {1.223, 0.083}, {1.529, 0.790},
    {1.625, 0.697}, {2.274, 0.953},
};

static int matchedRings(const Solution& s, double tol) {
    if (s.dist.empty()) return -1; // existing implementation does not fill dist -> not applicable
    int n = 0;
    for (size_t i = 0; i < s.hkl.size() && i < s.dist.size(); ++i)
        if ((s.hkl[i].h || s.hkl[i].k || s.hkl[i].l) && std::isfinite(s.dist[i]) && s.dist[i] <= tol)
            ++n;
    return n;
}

static void printSols(const char* tag, const std::vector<Solution>& sols, double seconds, int nTop = 5) {
    std::printf("  [%s] %zu sols, %.3f s\n", tag, sols.size(), seconds);
    for (int i = 0; i < (int)sols.size() && i < nTop; ++i) {
        const Solution& s = sols[(size_t)i];
        const int mr = matchedRings(s, 0.05);
        std::printf("    #%d uvw=(%d %d %d) a=%.3f b=%.3f c=%.3f al=%.2f be=%.2f ga=%.2f "
                    "V=%.1f rmsQxyz=%.4f rmsQxy=%.4f rmsQz=%.4f qspec=%.4f matched=%s\n",
                    i + 1, s.uvw.h, s.uvw.k, s.uvw.l, s.cell.a, s.cell.b, s.cell.c,
                    s.cell.al, s.cell.be, s.cell.ga, s.vol, s.rmsQxyz, s.rmsQxy, s.rmsQz,
                    s.qspec, mr < 0 ? "n/a" : std::to_string(mr).c_str());
    }
}

static bool manualHit(const Solution& s) {
    const double mc[3] = {7.13, 8.48, 16.67};
    const double mAng[3] = {89.4, 87.8, 89.7};
    const double sc[3] = {s.cell.a, s.cell.b, s.cell.c};
    const double sAng[3] = {s.cell.al, s.cell.be, s.cell.ga};
    for (int i = 0; i < 3; ++i)
        if (std::fabs(sc[i] - mc[i]) > 0.5 || std::fabs(sAng[i] - mAng[i]) > 3.0) return false;
    const int u0 = s.uvw.h, v0 = s.uvw.k, w0 = s.uvw.l;
    // Parallel to (-1,2,1) (an overall sign flip is allowed)
    const bool same = (u0 == -1 && v0 == 2 && w0 == 1);
    const bool neg = (u0 == 1 && v0 == -2 && w0 == -1);
    return same || neg;
}

static bool similarCell(const Solution& a, const Solution& b) {
    return std::fabs(a.cell.a - b.cell.a) < 0.5 && std::fabs(a.cell.b - b.cell.b) < 0.5 &&
           std::fabs(a.cell.c - b.cell.c) < 0.5 && std::fabs(a.cell.al - b.cell.al) < 3.0 &&
           std::fabs(a.cell.be - b.cell.be) < 3.0 && std::fabs(a.cell.ga - b.cell.ga) < 3.0;
}

static std::vector<Peak> loadCsv(const std::string& path, bool& ok) {
    std::vector<Peak> out;
    std::ifstream in(path);
    ok = (bool)in;
    if (!ok) return out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        double a = 0, b = 0;
        if (std::sscanf(line.c_str(), "%lf , %lf", &a, &b) != 2) continue;
        out.push_back(Peak{a, b});
    }
    return out;
}

static const char* algName(int a) {
    switch (a) {
        case 0: return "Original";
        case 1: return "Improved";
        case 2: return "PeriodVoteDdl";
        case 3: return "Ransac";
        case 4: return "Ills";
        case 5: return "Scrm";
    }
    return "?";
}

// Run only one Part1 method, for cross-method comparison (uvw sign representatives, time, solution count, first-solution rms)
static void runAlg(const char* name, const std::vector<Peak>& peaks, double qspec, int alg,
                   int nTop) {
    SpecParams p;
    p.qspec = qspec;
    p.part1Alg = (Part1Algorithm)alg;
    Progress prog;
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<Solution> sols = runSpecIndex(peaks, p, prog);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const std::string tag = std::string("[") + name + " " + algName(alg) + "]";
    printSols(tag.c_str(), sols, sec, nTop);
}

static void runBoth(const char* name, const std::vector<Peak>& peaks, double qspec,
                    const scrm::Params& sprm, int nTop) {
    std::printf("\n=== %s : %zu peaks, qspec=%.4f ===\n", name, peaks.size(), qspec);

    SpecParams eprm;
    eprm.qspec = qspec;
    Progress prog;
    auto t0 = std::chrono::steady_clock::now();
    auto es = runSpecIndex(peaks, eprm, prog);
    const double et = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printSols("existing runSpecIndex", es, et, nTop);

    scrm::Stats st;
    auto t1 = std::chrono::steady_clock::now();
    auto ss = scrm::runSpecIndex(peaks, sprm, prog, &st);
    const double stime = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    std::printf("  [scrm] stats: uvDirs=%d ipCands=%lld normalCalls=%lld secondTrials=%lld "
                "rawSol=%lld finalSol=%lld seedL=%d fallback=%d\n",
                st.uvDirs, st.ipCands, st.normalCalls, st.secondTrials, st.rawSol,
                st.finalSol, st.seedLUsed, (int)st.wideFallback);
    printSols("scrm", ss, stime, nTop);

    if (!es.empty()) {
        int mir = -1;
        for (size_t i = 0; i < ss.size(); ++i)
            if (similarCell(ss[i], es[0])) { mir = (int)i; break; }
        std::printf("  does scrm contain the existing#1 cell: %s\n",
                    mir >= 0 ? ("YES #" + std::to_string(mir + 1)).c_str() : "NO");
    }

    if (std::string(name) == "DIP") {
        int hit = -1;
        for (size_t i = 0; i < ss.size(); ++i)
            if (manualHit(ss[i])) { hit = (int)i; break; }
        if (hit >= 0)
            std::printf("  manual hit: YES at #%d\n", hit + 1);
        else
            std::printf("  manual hit: NO\n");
    }
}

int main(int argc, char** argv) {
    scrm::Params sprm;
    std::string csv = "D:/code/giwaxs-examples/peak_finder/test.csv";
    // Optional: --w N overrides wLimit; --csv PATH specifies the CSV; fixed contact plane direction u v w (2nd-4th positional args)
    int pos = 0;
    int nTop = 5;
    int algSel = 0; // 0=off; 1..4=single method; -1=run all four methods
    int fixed[3] = {0, 0, 1};
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--alg") == 0 && i + 1 < argc) {
            const char* a = argv[++i];
            algSel = (std::strcmp(a, "all") == 0) ? -1 : std::atoi(a);
            continue;
        }
        if (std::strcmp(argv[i], "--w") == 0 && i + 1 < argc) {
            sprm.wLimit = std::atoi(argv[++i]);
            continue;
        }
        if (std::strcmp(argv[i], "--csv") == 0 && i + 1 < argc) {
            csv = argv[++i];
            continue;
        }
        if (std::strcmp(argv[i], "--top") == 0 && i + 1 < argc) {
            nTop = std::atoi(argv[++i]);
            continue;
        }
        if (std::strcmp(argv[i], "--refine-top") == 0 && i + 1 < argc) {
            sprm.refineTop = std::atoi(argv[++i]);
            continue;
        }
        if (pos < 3) fixed[pos++] = std::atoi(argv[i]);
    }
    if (pos == 3) {
        sprm.fixedUvw = true;
        sprm.fu = fixed[0];
        sprm.fv = fixed[1];
        sprm.fw = fixed[2];
        std::printf("[cli] fixed uvw = (%d %d %d)\n", sprm.fu, sprm.fv, sprm.fw);
    }
    std::printf("[cli] wLimit = %d, csv = %s\n", sprm.wLimit, csv.c_str());

    double qsDip = 0.0;
    int rods = 0;
    detectQSpec(kDip, qsDip, rods);
    sprm.qspec = qsDip;
    runBoth("DIP", kDip, qsDip, sprm, nTop);
    if (algSel != 0) {
        const int a0 = (algSel < 0) ? 1 : algSel;
        const int a1 = (algSel < 0) ? 4 : algSel;
        for (int a = a0; a <= a1; ++a) runAlg("DIP", kDip, qsDip, a, nTop);
    }

    bool ok = false;
    const std::vector<Peak> test = loadCsv(csv, ok);
    if (!ok || test.empty()) {
        std::printf("\n[csv] skipped (%s not found or empty)\n", csv.c_str());
        return 0;
    }
    double qsCsv = 0.0;
    int rodsCsv = 0;
    detectQSpec(test, qsCsv, rodsCsv);
    sprm.qspec = qsCsv;
    runBoth("CSV", test, qsCsv, sprm, nTop);
    if (algSel != 0) {
        const int a0 = (algSel < 0) ? 1 : algSel;
        const int a1 = (algSel < 0) ? 4 : algSel;
        for (int a = a0; a <= a1; ++a) runAlg("CSV", test, qsCsv, a, nTop);
    }
    return 0;
}
