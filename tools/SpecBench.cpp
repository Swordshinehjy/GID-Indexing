// ============================================================================
//  SpecBench.cpp — benchmark of the q_spec implementations (improved version / PeriodVote+DDL /
//  RANSAC / integer least squares / SCRM)
//
//  SCRM is an independent implementation (src/core/scrm.cpp); step one is timed with
//  scrm::runInPlane; the other four share improved::runPartOne.
//  Usage:
//    GIDSpecBench              default parameters, each case run 3 times taking the best
//    GIDSpecBench --l 8 --w 4 --R 6 --uv 2 --hk1 3 --rows 5
// ============================================================================
#include "scrm.h"
#include "specindex.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace gid;

// DIP measured peak table (12 peaks, including the specular rod peak {0.000, 1.776})
static const std::vector<Peak> kDip = {
    {0.000, 1.776}, {0.368, 0.091}, {0.396, 1.160}, {0.428, 1.069}, {0.689, 0.984},
    {0.734, 0.182}, {0.806, 1.232}, {1.140, 0.171}, {1.223, 0.083}, {1.529, 0.790},
    {1.625, 0.697}, {2.274, 0.953},
};

// Synthetic triclinic cell (same ground truth as tools/SelfTest.cpp), peak table from forward modeling
struct Case {
    const char* name;
    std::vector<Peak> peaks;
    double qspec = 0.0;
};

static double msSince(const std::chrono::steady_clock::time_point& t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// DIP manual solution: uvw=(-1,2,1) or the equivalent (1,-2,-1), in-plane ~7.15/8.49, c ~16.7
static bool isManualDip(const Solution& s) {
    const double a = std::min(s.cell.a, s.cell.b), b = std::max(s.cell.a, s.cell.b);
    return std::fabs(a - 7.15) < 0.25 && std::fabs(b - 8.49) < 0.25 && std::fabs(s.cell.c - 16.70) < 0.6;
}

static const char* algTag(Part1Algorithm a) {
    switch (a) {
        case Part1Algorithm::Original: return "Original";
        case Part1Algorithm::PeriodVoteDdl: return "PeriodVote/DDL";
        case Part1Algorithm::Ransac: return "RANSAC ";
        case Part1Algorithm::Ills: return "ILS/int ";
        case Part1Algorithm::Scrm: return "SCRM    ";
        case Part1Algorithm::Improved: break;
    }
    return "Improved";
}

static void runOne(const Case& cs, Part1Algorithm alg, const SpecParams& base, int rows) {
    SpecParams sp = base;
    sp.qspec = cs.qspec;
    sp.part1Alg = alg;

    // Time step one separately (the improved family uses improved::runPartOne; SCRM uses scrm::runInPlane)
    size_t ncand = 0;
    double p1Best = 1e30;
    for (int r = 0; r < rows; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        if (alg == Part1Algorithm::Scrm) {
            ncand = scrm::runInPlane(cs.peaks, scrm::fromSpecParams(sp), Progress{}, nullptr).size();
        } else {
            ncand = improved::runPartOne(cs.peaks, sp, Progress{}).size();
        }
        p1Best = std::min(p1Best, msSince(t0));
    }

    // Full pipeline (Part ONE + Part TWO)
    std::vector<Solution> sols;
    double best = 1e30;
    for (int r = 0; r < rows; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        sols = runSpecIndex(cs.peaks, sp, Progress{});
        best = std::min(best, msSince(t0));
    }

    const char* tag = algTag(alg);
    int hit = -1, nhit = 0;
    for (size_t i = 0; i < sols.size(); ++i)
        if (isManualDip(sols[i])) {
            if (hit < 0) hit = (int)i;
            ++nhit;
        }
    std::printf("  %s  P1 %6.1f ms (%7zu cands)  total %7.1f ms  %2zu sols  manual sol #%d (%d hits)",
                tag, p1Best, ncand, best, sols.size(), hit, nhit);
    if (!sols.empty()) {
        const Solution& s = sols.front();
        std::printf("  top1 uvw=(%d %d %d) a=%.3f b=%.3f c=%.3f rms=%.4f", s.uvw.h, s.uvw.k,
                    s.uvw.l, s.cell.a, s.cell.b, s.cell.c, s.rmsQxyz);
    }
    std::printf("\n");
    std::fflush(stdout);
}

// Part TWO cost for a single candidate (single-threaded, excluding Part ONE / candidate-level parallelism)
static void benchPart2(const Case& cs, const SpecParams& base, int rows) {
    SpecParams sp = base;
    sp.qspec = cs.qspec;
    sp.part1Alg = Part1Algorithm::Improved;
    const std::vector<AbgCand> abg = improved::runPartOne(cs.peaks, sp, Progress{});
    if (abg.empty()) {
        std::printf("  [Part2] no Part ONE candidates\n");
        return;
    }
    const AbgCand& cand = abg.front();
    std::printf("\n[%s] single-candidate Part TWO (u=%d v=%d a=%.3f b=%.3f ga=%.2f, single-threaded)\n", cs.name,
                cand.u, cand.v, cand.a, cand.b, cand.ga);
    for (Part1Algorithm alg : {Part1Algorithm::Improved, Part1Algorithm::PeriodVoteDdl,
                               Part1Algorithm::Ransac, Part1Algorithm::Ills}) {
        // SCRM's second step is an independent analytic implementation (scrm::runNormal) and does not join this table's Part TWO comparison
        SpecParams p2 = sp;
        p2.part1Alg = alg;
        double best = 1e30;
        size_t nsol = 0;
        for (int r = 0; r < rows; ++r) {
            const auto t0 = std::chrono::steady_clock::now();
            const std::vector<Solution> sols = improved::runPartTwo(cs.peaks, cand, p2);
            best = std::min(best, msSince(t0));
            nsol = sols.size();
        }
        std::printf("  %s %8.2f ms  (%zu sols)\n", algTag(alg), best, nsol);
    }
    std::fflush(stdout);
}

// Per-call cost of each atom in the refinement stage (to see where Part TWO time goes)
static void micro() {
    Cell c{9.1, 11.3, 15.7, 88.0, 95.0, 99.0};
    const Hkl uvw{0, 0, 1};
    const Orient o = Orient::fromUvw(uvw);
    std::vector<Peak> peaks;
    std::vector<Hkl> ah;
    for (int h = -2; h <= 2; ++h)
        for (int k = -2; k <= 2; ++k)
            for (int l = 0; l <= 2; ++l) {
                if (!h && !k && !l) continue;
                const QPoint q = forwardHkl(reciprocalMatrix(c), uvw, h, k, l);
                if (q.qxy < 0.05) continue;
                peaks.push_back(Peak{q.qxy, q.qz});
                ah.push_back(Hkl{h, k, l});
            }
    while (peaks.size() > 21) {
        peaks.pop_back();
        ah.pop_back();
    }
    const std::vector<Hkl> cands = hklGrid(6, 6, 6);
    double qmax = 0.0;
    for (const Peak& p : peaks) qmax = std::max(qmax, std::hypot(p.qxy, p.qz));
    const std::array<bool, 6> freeAll{true, true, true, true, true, true};
    const int nRep = 3000;
    std::vector<Hkl> out;

    auto bench = [&](const char* name, auto&& fn) {
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < nRep; ++i) fn();
        std::printf("  %-34s %8.2f us\n", name, msSince(t0) * 1000.0 / nRep);
    };

    std::printf("\n[micro benchmark] %zu peaks, candidate grid %zu\n", peaks.size(), cands.size());
    bench("assignHkl(no |G| filter)", [&] {
        assignHkl(c, o, peaks, cands, out);
    });
    bench("assignHkl(maxG=qmax+0.2)", [&] {
        assignHkl(c, o, peaks, cands, out, nullptr, 1e30, 1.0, qmax + 0.2);
    });
    bench("refineCell(60 iter)", [&] {
        Cell cc = c;
        refineCell(peaks, ah, o, cc, freeAll, 60);
    });
    bench("reduceCell", [&] { reduceCell(c); });
    bench("evaluate", [&] { evaluate(c, o, peaks, ah); });
    std::fflush(stdout);
}

int main(int argc, char** argv) {
    SpecParams sp; // the defaults in specindex.h are the UI default parameters
    int rows = 3;
    bool wantMicro = false;
    bool withOrig = false; // the original GIDInd is very slow, so it is off by default
    bool wantPart2 = false;
    for (int i = 1; i < argc; ++i) {
        auto next = [&](double dflt) {
            if (i + 1 < argc) return std::atof(argv[++i]);
            return dflt;
        };
        if (!std::strcmp(argv[i], "--l")) sp.lLimit = (int)next(sp.lLimit);
        else if (!std::strcmp(argv[i], "--w")) sp.wLimit = (int)next(sp.wLimit);
        else if (!std::strcmp(argv[i], "--R")) sp.redundant = (int)next(sp.redundant);
        else if (!std::strcmp(argv[i], "--uv")) sp.uvLimit = (int)next(sp.uvLimit);
        else if (!std::strcmp(argv[i], "--hk1")) sp.hkLimit1 = (int)next(sp.hkLimit1);
        else if (!std::strcmp(argv[i], "--hk2")) sp.hkLimit2 = (int)next(sp.hkLimit2);
        else if (!std::strcmp(argv[i], "--rows")) rows = (int)next(rows);
        else if (!std::strcmp(argv[i], "--max1")) sp.maxPart1 = (int)next(sp.maxPart1);
        else if (!std::strcmp(argv[i], "--seed3")) sp.seed3Half = (int)next(sp.seed3Half);
        else if (!std::strcmp(argv[i], "--refine")) sp.refineTop = (int)next(sp.refineTop);
        else if (!std::strcmp(argv[i], "--no-early")) sp.earlyStop = false;
        else if (!std::strcmp(argv[i], "--micro")) wantMicro = true;
        else if (!std::strcmp(argv[i], "--orig")) withOrig = true;
        else if (!std::strcmp(argv[i], "--part2")) wantPart2 = true;
    }
    std::printf("params: uv=%d hk1=%d hk2=%d R=%d l=%d w=%d maxPart1=%d seed3Half=%d refineTop=%d "
                "earlyStop=%d rows=%d\n",
                sp.uvLimit, sp.hkLimit1, sp.hkLimit2, sp.redundant, sp.lLimit, sp.wLimit,
                sp.maxPart1, sp.seed3Half, sp.refineTop, (int)sp.earlyStop, rows);

    std::vector<Case> cases;
    cases.push_back({"DIP 12 peaks", kDip, 1.776});
    // Synthetic triclinic: forward model a=9.1 b=11.3 c=15.7 al=88 be=95 ga=99, uvw=(0,0,1)
    {
        Cell c{9.1, 11.3, 15.7, 88.0, 95.0, 99.0};
        const Hkl uvw{0, 0, 1};
        std::vector<Peak> pk;
        for (int h = -2; h <= 2; ++h)
            for (int k = -2; k <= 2; ++k)
                for (int l = 0; l <= 2; ++l) {
                    if (!h && !k && !l) continue;
                    const QPoint q = forwardHkl(reciprocalMatrix(c), uvw, h, k, l);
                    if (q.qxy < 0.05) continue;
                    pk.push_back(Peak{q.qxy, q.qz});
                }
        std::sort(pk.begin(), pk.end(), [](const Peak& x, const Peak& y) { return x.qxy < y.qxy; });
        if (pk.size() > 20) pk.resize(20);
        // Specular rod peak: q_z = |(0,0,1) direction| = qSpecOf
        const double qs = qSpecOf(reciprocalMatrix(c), uvw);
        pk.insert(pk.begin(), Peak{0.0, qs});
        cases.push_back({"synthetic triclinic 20 peaks", pk, qs});
    }

    for (const Case& cs : cases) {
        std::printf("\n[%s] %zu peaks  q_spec=%.4f\n", cs.name, cs.peaks.size(), cs.qspec);
        if (withOrig) runOne(cs, Part1Algorithm::Original, sp, 1);
        runOne(cs, Part1Algorithm::Improved, sp, rows);
        runOne(cs, Part1Algorithm::PeriodVoteDdl, sp, rows);
        runOne(cs, Part1Algorithm::Ransac, sp, rows);
        runOne(cs, Part1Algorithm::Ills, sp, rows);
        runOne(cs, Part1Algorithm::Scrm, sp, rows);
        if (wantPart2) benchPart2(cs, sp, 2);
    }
    if (wantMicro) micro();
    return 0;
}
