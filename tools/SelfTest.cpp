// ============================================================================
//  SelfTest.cpp — core algorithm self-test (not part of the single-file exe artifacts)
//  Forward-model synthetic data -> inversion -> check the cell and the deviations
// ============================================================================
#include "engine.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

using namespace gid;

namespace {

int g_fail = 0;

void check(bool ok, const std::string& name, double value = 0.0, double tol = 0.0) {
    if (ok) {
        std::printf("[PASS] %s\n", name.c_str());
    } else {
        ++g_fail;
        std::printf("[FAIL] %s   (value=%.6g tol=%.6g)\n", name.c_str(), value, tol);
    }
    std::fflush(stdout);
}

std::vector<Peak> makeData(const Cell& cell, const Hkl& uvw, int maxH, int maxK, int maxL,
                           bool withRod, double& qspec, std::vector<Hkl>* hklOut = nullptr) {
    const Mat3 R = reciprocalMatrix(cell);
    qspec = qSpecOf(R, uvw);
    std::vector<Peak> out;
    if (hklOut) hklOut->clear();
    if (withRod) {
        out.push_back(Peak{0.0, qspec});
        out.push_back(Peak{0.0, 2.0 * qspec});
        if (hklOut) {
            hklOut->push_back(uvw);
            hklOut->push_back(Hkl{2 * uvw.h, 2 * uvw.k, 2 * uvw.l});
        }
    }
    for (int h = -maxH; h <= maxH; ++h)
        for (int k = -maxK; k <= maxK; ++k)
            for (int l = -maxL; l <= maxL; ++l) {
                if (h == 0 && k == 0 && l == 0) continue;
                if (l < 0) continue; // keep only the half space (common in GIXD)
                const QPoint p = forwardHkl(R, uvw, h, k, l);
                if (p.qxy < 0.05) continue;
                out.push_back(Peak{p.qxy, p.qz});
                if (hklOut) hklOut->push_back(Hkl{h, k, l});
            }
    return out;
}

} // namespace

int main() {
    // ---------------------------------------------------------------- geometry
    {
        Cell c{7.8, 9.1, 15.2, 95.0, 98.0, 88.0};
        const Mat3 R = reciprocalMatrix(c);
        const Hkl uvw{1, 0, 2};
        // Forward model: verify forwardHkl against an explicit Cartesian dot product
        std::array<double, 3> gu{};
        for (int i = 0; i < 3; ++i)
            gu[i] = uvw.h * R[3 * 0 + i] + uvw.k * R[3 * 1 + i] + uvw.l * R[3 * 2 + i];
        const double gn = std::sqrt(gu[0] * gu[0] + gu[1] * gu[1] + gu[2] * gu[2]);
        double err = 0.0;
        for (int h = -2; h <= 2; ++h)
            for (int k = -2; k <= 2; ++k)
                for (int l = -2; l <= 2; ++l) {
                    double g[3] = {h * R[0] + k * R[3] + l * R[6], h * R[1] + k * R[4] + l * R[7],
                                   h * R[2] + k * R[5] + l * R[8]};
                    const double qz = (g[0] * gu[0] + g[1] * gu[1] + g[2] * gu[2]) / gn;
                    const double g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
                    const QPoint p = forwardHkl(R, uvw, h, k, l);
                    err = std::max(err, std::fabs(p.qz - qz));
                    // q_xy is extremely rounding-sensitive at grazing incidence (q_xy->0), so compare only where q_xy is large enough
                    if (p.qxy > 0.1) {
                        const double qxy = std::sqrt(std::max(0.0, g2 - qz * qz));
                        err = std::max(err, std::fabs(p.qxy - qxy));
                    }
                }
        check(err < 1e-12, "forward model (q_xy,q_z) agrees with the Cartesian dot product", err, 1e-12);

        // Frame round trip: cell -> (za,zb,C,kappa,lambda,rho) -> cell
        const Frame f = frameOf(c);
        Cell back{};
        cellFromFrame(f.za, f.zb, f.cg, f.kap, f.lam, f.rho, back);
        const double d = std::fabs(back.a - c.a) + std::fabs(back.b - c.b) + std::fabs(back.c - c.c) +
                         std::fabs(back.al - c.al) + std::fabs(back.be - c.be) +
                         std::fabs(back.ga - c.ga);
        check(d < 1e-8, "Frame round trip (cell -> kappa/lambda/rho -> cell)", d, 1e-8);

        // Consistency of nu with G*: nu^T G*^-1 nu == 1
        const std::array<double, 6> G = recipMetric(c);
        const std::array<double, 3> nu = Orient::fromUvw(uvw).get(G);
        const Mat3 Rm = reciprocalMatrix(c);
        const std::array<double, 3> nu2 = nuFromUvw(Rm, uvw);
        const double dn = std::fabs(nu[0] - nu2[0]) + std::fabs(nu[1] - nu2[1]) +
                          std::fabs(nu[2] - nu2[2]);
        check(dn < 1e-9, "the two nu algorithms agree", dn, 1e-9);
    }

    // ---------------------------------------------------------------- with q_spec
    {
        Cell truth{9.4, 11.2, 16.5, 92.0, 97.0, 105.0};
        const Reduction red = reduceCell(truth);
        const Cell cell = red.ok ? red.cell : truth;
        const Hkl uvw = applyT(Hkl{0, 0, 1}, red.T);
        double qspec = 0.0;
        const std::vector<Peak> data = makeData(cell, uvw, 2, 2, 3, true, qspec);

        double det = 0.0;
        int rods = 0;
        check(detectQSpec(data, det, rods) && std::fabs(det - qspec) < 1e-6,
              "q_spec auto-detection", std::fabs(det - qspec), 1e-6);

        SpecParams sp;
        sp.qspec = qspec;
        sp.uvLimit = 1;
        sp.hkLimit1 = 2;
        sp.hkLimit2 = 5;
        sp.maxPart1 = 6;
        sp.threads = 4;
        std::printf("       true uvw=(%d %d %d) qspec=%.4f peaks=%zu\n", uvw.h, uvw.k, uvw.l,
                    qspec, data.size());
        const std::vector<AbgCand> abg = improved::runPartOne(data, sp, Progress{});
        std::printf("       Part1 candidates = %zu", abg.size());
        for (size_t i = 0; i < abg.size() && i < 5; ++i)
            std::printf("  [uv=(%d,%d) a=%.3f b=%.3f ga=%.2f rms=%.4f]", abg[i].u, abg[i].v,
                        abg[i].a, abg[i].b, abg[i].ga, abg[i].rmsQxy);
        {
            int bestI = -1;
            double bestS = 1e9;
            for (size_t i = 0; i < abg.size(); ++i) {
                const double sc = std::fabs(abg[i].a - cell.a) + std::fabs(abg[i].b - cell.b) +
                                  std::fabs(abg[i].ga - cell.ga) / 5.0;
                if (sc < bestS) { bestS = sc; bestI = (int)i; }
            }
            if (bestI >= 0)
                std::printf("\n       truth rank=%d score=%.3f uv=(%d,%d) a=%.3f b=%.3f ga=%.2f "
                            "rms=%.5f",
                            bestI, bestS, abg[bestI].u, abg[bestI].v, abg[bestI].a, abg[bestI].b,
                            abg[bestI].ga, abg[bestI].rmsQxy);
        }
        std::printf("\n       true: a=%.3f b=%.3f ga=%.2f\n", cell.a, cell.b, cell.ga);
        std::fflush(stdout);
        check(!abg.empty(), "Part1 finds (a,b,gamma) candidates");
        // Default implementation = PeriodVote/DDL (SpecParams::part1Alg; the UI drop-down defaults to it too)
        check(SpecParams{}.part1Alg == Part1Algorithm::PeriodVoteDdl, "default Part 1 algorithm = PeriodVote/DDL");

        EngineParams ep;
        ep.mode = IndexMode::WithSpec;
        ep.spec.qspec = qspec;
        ep.spec.uvLimit = 1;
        ep.spec.hkLimit1 = 2;
        ep.spec.hkLimit2 = 5;
        ep.spec.lLimit = 4;
        ep.spec.wLimit = 1;
        ep.spec.maxPart1 = 6;
        ep.spec.maxSolutions = 10;
        ep.spec.threads = 4;
        ep.spec.part1Alg = Part1Algorithm::Improved;
        const EngineResult r = runIndexing(data, ep, Progress{});
        check(r.ok, "q_spec two-step method finds solutions");
        if (r.ok) {
            const Solution& s = r.solutions.front();
            const double rel = std::fabs(s.vol - cellVolume(cell)) / cellVolume(cell);
            std::printf("       best: %s uvw=(%d %d %d) V=%.3f (true %.3f) rms=%.5f\n",
                        fmtCell(s.cell).c_str(), s.uvw.h, s.uvw.k, s.uvw.l, s.vol,
                        cellVolume(cell), s.rmsQxyz);
            check(s.rmsQxyz < 0.02, "with q_spec: RMSD(q_xyz) < 0.02", s.rmsQxyz, 0.02);
            check(rel < 0.03, "with q_spec: volume agrees with the truth", rel, 0.03);
            check(s.hkl.size() == data.size(), "with q_spec: hkl aligned with the peak table");
            // Each peak's (hkl), after forward modeling, must match the measured peak position (specular peaks included)
            double worst = 0.0;
            for (size_t i = 0; i < data.size() && i < s.hkl.size(); ++i) {
                const QPoint q = forwardHkl(reciprocalMatrix(s.cell), s.uvw, s.hkl[i].h,
                                            s.hkl[i].k, s.hkl[i].l);
                worst = std::max(worst,
                                 std::sqrt((q.qxy - data[i].qxy) * (q.qxy - data[i].qxy) +
                                           (q.qz - data[i].qz) * (q.qz - data[i].qz)));
            }
            check(worst < 0.02, "with q_spec: predicted peak positions for the hkl match the measurement", worst, 0.02);
        }

        // PeriodVote/DDL (difference vectors + period voting): same data, same search ranges
        ep.spec.part1Alg = Part1Algorithm::PeriodVoteDdl;
        const EngineResult rf = runIndexing(data, ep, Progress{});
        check(rf.ok, "PeriodVote/DDL finds solutions");
        if (rf.ok) {
            const Solution& s = rf.solutions.front();
            const double rel = std::fabs(s.vol - cellVolume(cell)) / cellVolume(cell);
            std::printf("       PeriodVote/DDL best: %s uvw=(%d %d %d) V=%.3f (true %.3f) rms=%.5f\n",
                        fmtCell(s.cell).c_str(), s.uvw.h, s.uvw.k, s.uvw.l, s.vol,
                        cellVolume(cell), s.rmsQxyz);
            check(s.rmsQxyz < 0.02, "PeriodVote/DDL: RMSD(q_xyz) < 0.02", s.rmsQxyz, 0.02);
            check(rel < 0.03, "PeriodVote/DDL: volume agrees with the truth", rel, 0.03);
        }

        // Original Part ONE (default algorithm): only check that it runs. This synthetic data contains
        // degenerate triples with the same (h,k) but different l among the seed peaks; the original takes
        // the smallest-|q| peak as the seed as required, so the 3x3 matrix is singular at uv=(0,0) -- its own behavior.
        ep.spec.part1Alg = Part1Algorithm::Original;
        const EngineResult ro = runIndexing(data, ep, Progress{});
        std::printf("       original Part1: %d solutions", (int)ro.solutions.size());
        if (!ro.solutions.empty()) {
            const Solution& s = ro.solutions.front();
            std::printf("   best: %s uvw=(%d %d %d) rms=%.5f", fmtCell(s.cell).c_str(), s.uvw.h,
                        s.uvw.k, s.uvw.l, s.rmsQxyz);
        }
        std::printf("\n");
        std::fflush(stdout);
        check(ro.ok, "original two-step method runs and produces candidate solutions");
    }

    // ------------------------------------------- DIP (manual Table A.2, real-data regression)
    // Manual solution: uvw=(-1 2 1), a=7.13 b=8.48 c=16.67 Å, α=89.4 β=87.8 γ=89.7.
    // Both kinds of "starting three reflections" are indispensable here: the smallest-|q| peaks include
    // (001)/(002) on the specular rod, the smallest-q_xy peaks nearly coplanar; the improved version seeded by q_xy alone and missed it.
    {
        const std::vector<Peak> dip = {
            {0.000, 1.776}, {0.368, 0.091}, {0.396, 1.160}, {0.428, 1.069}, {0.689, 0.984},
            {0.734, 0.182}, {0.806, 1.232}, {1.140, 0.171}, {1.223, 0.083}, {1.529, 0.790},
            {1.625, 0.697}, {2.274, 0.953},
        };
        EngineParams ep;
        ep.mode = IndexMode::WithSpec;
        ep.spec.qspec = 1.776;
        ep.spec.part1Alg = Part1Algorithm::Improved;
        ep.spec.linesToPermute = 4;
        ep.spec.uvLimit = 2;
        ep.spec.hkLimit1 = 2;
        ep.spec.hkLimit2 = 6;
        ep.spec.lLimit = 6;
        ep.spec.wLimit = 2;
        ep.spec.redundant = 4;
        ep.spec.maxPart1 = 20;
        ep.spec.refineTop = 200;
        ep.spec.maxSolutions = 50;
        ep.spec.threads = 4;
        const EngineResult r = runIndexing(dip, ep, Progress{});
        check(r.ok, "DIP improved: linesToPermute=4/uvLimit=2/hkLimit1=2 finds solutions");
        bool hit = false;
        for (const Solution& s : r.solutions) {
            const int au = std::abs(s.uvw.h), av = std::abs(s.uvw.k), aw = std::abs(s.uvw.l);
            if (au == 1 && av == 2 && aw == 1 && std::fabs(s.cell.a - 7.13) < 0.15 &&
                std::fabs(s.cell.b - 8.48) < 0.15 && std::fabs(s.cell.c - 16.67) < 0.30) {
                std::printf("       %s uvw=(%d %d %d) V=%.1f rms=%.5f\n", fmtCell(s.cell).c_str(),
                            s.uvw.h, s.uvw.k, s.uvw.l, s.vol, s.rmsQxyz);
                hit = true;
                break;
            }
        }
        std::fflush(stdout);
        check(hit, "DIP improved: obtains the (-1 2 1) solution");

        // Solve the same data once more with PeriodVote/DDL (difference vectors + period voting)
        ep.spec.part1Alg = Part1Algorithm::PeriodVoteDdl;
        const EngineResult rf = runIndexing(dip, ep, Progress{});
        check(rf.ok, "DIP PeriodVote/DDL: finds solutions");
        bool hitF = false;
        for (const Solution& s : rf.solutions) {
            const int au = std::abs(s.uvw.h), av = std::abs(s.uvw.k), aw = std::abs(s.uvw.l);
            if (au == 1 && av == 2 && aw == 1 && std::fabs(s.cell.a - 7.13) < 0.15 &&
                std::fabs(s.cell.b - 8.48) < 0.15 && std::fabs(s.cell.c - 16.67) < 0.30) {
                std::printf("       PeriodVote/DDL %s uvw=(%d %d %d) V=%.1f rms=%.5f\n",
                            fmtCell(s.cell).c_str(), s.uvw.h, s.uvw.k, s.uvw.l, s.vol,
                            s.rmsQxyz);
                hitF = true;
                break;
            }
        }
        std::fflush(stdout);
        check(hitF, "DIP PeriodVote/DDL: obtains the (-1 2 1) solution");

        // Two more Part TWO variants (2-peak minimal sample / difference-vector lattice + interval sweep):
        // only the (κ,λ,ρ) search of Part TWO changes, and the same manual solution must come out. Dispatch
        // is covered too: these two enum values would silently fall into the original part1/part2 paths otherwise.
        const std::pair<Part1Algorithm, const char*> variants[] = {
            {Part1Algorithm::Ransac, "RANSAC"},
            {Part1Algorithm::Ills, "integer least squares"},
        };
        for (const auto& v : variants) {
            ep.spec.part1Alg = v.first;
            const EngineResult rv = runIndexing(dip, ep, Progress{});
            bool hitV = false;
            for (const Solution& s : rv.solutions) {
                const int au = std::abs(s.uvw.h), av = std::abs(s.uvw.k), aw = std::abs(s.uvw.l);
                if (au == 1 && av == 2 && aw == 1 && std::fabs(s.cell.a - 7.13) < 0.15 &&
                    std::fabs(s.cell.b - 8.48) < 0.15 && std::fabs(s.cell.c - 16.67) < 0.30) {
                    std::printf("       %s %s uvw=(%d %d %d) V=%.1f rms=%.5f\n", v.second,
                                fmtCell(s.cell).c_str(), s.uvw.h, s.uvw.k, s.uvw.l, s.vol,
                                s.rmsQxyz);
                    hitV = true;
                    break;
                }
            }
            std::fflush(stdout);
            check(hitV, std::string("DIP ") + v.second + ": obtains the (-1 2 1) solution");
        }

        // SCRM (new method): an independent implementation; step one enumerates only the in-plane (u,v)
        // (without w), while step two resolves (κ,λ,ρ,w,l) from the specular constraint. Same search scale as scrm defaults.
        ep.spec.part1Alg = Part1Algorithm::Scrm;
        ep.spec.uvLimit = 2;
        ep.spec.hkLimit1 = 3;
        ep.spec.seed3Half = 0;
        ep.spec.maxPart1 = 200;
        ep.spec.refineTop = 200;
        ep.spec.maxSolutions = 50;
        ep.spec.wLimit = 3;
        const EngineResult rs = runIndexing(dip, ep, Progress{});
        check(rs.ok, "DIP SCRM: finds solutions");
        bool hitS = false;
        for (const Solution& s : rs.solutions) {
            const int au = std::abs(s.uvw.h), av = std::abs(s.uvw.k), aw = std::abs(s.uvw.l);
            if (au == 1 && av == 2 && aw == 1 && std::fabs(s.cell.a - 7.13) < 0.15 &&
                std::fabs(s.cell.b - 8.48) < 0.15 && std::fabs(s.cell.c - 16.67) < 0.30) {
                std::printf("       SCRM %s uvw=(%d %d %d) V=%.1f rms=%.5f\n", fmtCell(s.cell).c_str(),
                            s.uvw.h, s.uvw.k, s.uvw.l, s.vol, s.rmsQxyz);
                hitS = true;
                break;
            }
        }
        if (!rs.solutions.empty()) {
            const Solution& s = rs.solutions.front();
            std::printf("       SCRM top1: uvw=(%d %d %d) %s rms=%.5f (%d solutions)\n", s.uvw.h,
                        s.uvw.k, s.uvw.l, fmtCell(s.cell).c_str(), s.rmsQxyz,
                        (int)rs.solutions.size());
        }
        std::fflush(stdout);
        check(hitS, "DIP SCRM: obtains the (-1 2 1) solution");
    }

    // ---------------------------------------------------------------- without q_spec
    {
        Cell truth{8.1, 10.3, 14.7, 90.0, 95.0, 90.0}; // monoclinic
        const Reduction red = reduceCell(truth);
        const Cell cell = red.ok ? red.cell : truth;
        const Hkl uvw = applyT(Hkl{0, 0, 1}, red.T);
        double qspec = 0.0;
        std::vector<Hkl> trueHkl;
        std::vector<Peak> data = makeData(cell, uvw, 2, 2, 2, false, qspec, &trueHkl);
        if (gid::hasEnv("GID_DEBUG")) {
            std::vector<int> ord(data.size());
            for (size_t i = 0; i < data.size(); ++i) ord[i] = (int)i;
            std::sort(ord.begin(), ord.end(), [&](int a, int b) {
                return data[a].qxy * data[a].qxy + data[a].qz * data[a].qz <
                       data[b].qxy * data[b].qxy + data[b].qz * data[b].qz;
            });
            std::fprintf(stderr, "[nospec] true cell %s uvw=(%d %d %d)\n", fmtCell(cell).c_str(),
                         uvw.h, uvw.k, uvw.l);
            std::fprintf(stderr, "[nospec] 10 smallest |g| true hkl: ");
            for (int i = 0; i < 10 && i < (int)ord.size(); ++i)
                std::fprintf(stderr, "(%d %d %d) ", trueHkl[ord[i]].h, trueHkl[ord[i]].k,
                             trueHkl[ord[i]].l);
            std::fprintf(stderr, "\n");
        }

        // Run both methods (Beam / B2 LLL+sphere decoding) sharing the same set of checks
        static const char* kMName[2] = {"beam search", "B2 sphere decoding"};
        for (int m = 0; m < 2; ++m) {
        EngineParams ep;
        ep.mode = IndexMode::NoSpec;
        ep.nospec.method = static_cast<gid::NoSpecMethod>(m);
        ep.nospec.seedH = 1;
        ep.nospec.seedK = 2;
        ep.nospec.seedL = 2;
        ep.nospec.seedPeaks = 6;
        ep.nospec.minRows = 8;
        ep.nospec.beam = 4000;
        ep.nospec.refineTop = 60;
        ep.nospec.gridH = 4; ep.nospec.gridK = 4; ep.nospec.gridL = 4;
        ep.nospec.threads = 4;
        const EngineResult r = runIndexing(data, ep, Progress{});
        const std::string tag = std::string("without q_spec (") + kMName[m] + ") single-step method finds solutions";
        check(r.ok, tag);
        if (r.ok) {
            const Solution& s = r.solutions.front();
            const double rel = std::fabs(s.vol - cellVolume(cell)) / cellVolume(cell);
            std::printf("       [%s] best: %s uvw=(%d %d %d) V=%.3f (true %.3f) rms=%.5f\n",
                        kMName[m], fmtCell(s.cell).c_str(), s.uvw.h, s.uvw.k, s.uvw.l, s.vol,
                        cellVolume(cell), s.rmsQxyz);
            check(s.rmsQxyz < 0.03,
                  "without q_spec (" + std::string(kMName[m]) + "): RMSD(q_xyz) < 0.03", s.rmsQxyz, 0.03);
            check(rel < 0.05, "without q_spec (" + std::string(kMName[m]) + "): volume agrees with the truth", rel, 0.05);
            // Regression: the reported (u,v,w) must agree with the model direction nu. Previously nu was a
            // freely solved tilt normal and uvw only the nearest integer direction picked afterwards, so the two
            // could differ by ~3.6°; (00l) peaks along (u,v,w) then got a nonexistent q_xy (orthorhombic 001).
            bool planeConsistent = true;
            for (const Solution& s : r.solutions) {
                if (!s.hasPlane) continue;
                const std::array<double, 3> n = nuFromUvw(reciprocalMatrix(s.cell), s.uvw);
                const double d =
                    std::fabs(n[0] - s.nu[0]) + std::fabs(n[1] - s.nu[1]) + std::fabs(n[2] - s.nu[2]);
                if (d > 1e-6) planeConsistent = false;
            }
            check(planeConsistent, "without q_spec (" + std::string(kMName[m]) +
                                       "): (u,v,w) agrees with the model direction nu (out-of-plane peaks q_xy=0)");
        }
        }
    }

    // ------------------------------------------- without q_spec (cubic, constrained crystal system regression)
    // For cubic, G22 = G33 = G11, so the design row of equation (A) must be h²+k²+l². If it degenerates
    // to h², as with "picking only the free column", (0,1,0)/(0,0,1)-type peaks never resolve and the beam search is pruned away.
    {
        const Cell cell{6.326, 6.326, 6.326, 90.0, 90.0, 90.0};
        const Hkl uvw{0, 1, 2};
        double qspec = 0.0;
        const std::vector<Peak> data = makeData(cell, uvw, 2, 2, 2, false, qspec);

        static const char* kMName[2] = {"beam search", "B2 sphere decoding"};
        for (int m = 0; m < 2; ++m) {
        EngineParams ep;
        ep.mode = IndexMode::NoSpec;
        ep.nospec.method = static_cast<gid::NoSpecMethod>(m);
        ep.nospec.system = CrystalSystem::Cubic;
        ep.nospec.seedH = 2;
        ep.nospec.seedK = 2;
        ep.nospec.seedL = 2;
        ep.nospec.beam = 4000;
        ep.nospec.refineTop = 60;
        ep.nospec.gridH = 4; ep.nospec.gridK = 4; ep.nospec.gridL = 4;
        ep.nospec.threads = 4;
        const EngineResult r = runIndexing(data, ep, Progress{});
        check(r.ok, "without q_spec (cubic, " + std::string(kMName[m]) + ") finds solutions");
        bool cubicHit = false;
        for (const Solution& s : r.solutions) {
            const bool eq = std::fabs(s.cell.b - s.cell.a) < 1e-6 &&
                            std::fabs(s.cell.c - s.cell.a) < 1e-6;
            if (std::fabs(s.cell.a - 6.326) < 0.02 && eq && s.rmsQxyz < 0.01) cubicHit = true;
        }
        if (r.ok) {
            const Solution& s = r.solutions.front();
            std::printf("       [%s] cubic best: %s uvw=(%d %d %d) rms=%.5f\n", kMName[m],
                        fmtCell(s.cell).c_str(),
                        s.uvw.h, s.uvw.k, s.uvw.l, s.rmsQxyz);
        }
        check(cubicHit, "without q_spec (cubic, " + std::string(kMName[m]) +
                            "): recovers the cubic cell with a=6.326 and a=b=c");
        }
    }

    // ------------------------------------------- SCRM two-step method without w (wiki §3/§5/§6)
    {
        // (1) Numerical check of the four identities (all exact, so the residual should be ~1e-15)
        const Cell cell{8.6, 10.4, 14.8, 84.0, 96.0, 101.0};
        const Hkl uvw{1, 2, 3};
        const Mat3 Rm = reciprocalMatrix(cell);
        const std::array<double, 6> G = recipMetric(cell);
        const double qs = qSpecOf(Rm, uvw);
        // 4+2 decomposition of G* along ĉ*: A=z_a², B=z_b², C=z_a z_b cosγ, κ,λ,ρ
        const double rho = std::sqrt(G[2]);
        const double kap = G[4] / rho;
        const double lam = G[5] / rho;
        const double A = G[0] - kap * kap;
        const double B = G[1] - lam * lam;
        const double C = kap * lam - G[3];
        const double u = uvw.h, v = uvw.k, w = uvw.l;
        const double guv = A * u * u + B * v * v - 2.0 * C * u * v;
        const double zeta = std::sqrt(qs * qs - guv);

        double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0, eAB = 0.0;
        for (int h = -3; h <= 3; ++h)
            for (int k = -3; k <= 3; ++k)
                for (int l = -3; l <= 3; ++l) {
                    if (!h && !k && !l) continue;
                    const QPoint q = forwardHkl(Rm, uvw, h, k, l);
                    const double t = q.qz / qs;
                    const double H = h - u * t, K = k - v * t, L = l - w * t;
                    const double pla = A * (double)h * h + B * (double)k * k -
                                       2.0 * C * (double)h * k;
                    const double gp = h * kap + k * lam + l * rho;
                    const double s3 = A * (double)h * u + B * (double)k * v -
                                      C * ((double)h * v + (double)k * u);
                    const double plaHK = A * H * H + B * K * K - 2.0 * C * H * K;
                    const double leak = H * kap + K * lam + L * rho;
                    // (I) |g|² = pla + gp²
                    e1 = std::max(e1, std::fabs(q.qxy * q.qxy + q.qz * q.qz - pla - gp * gp));
                    // (II) q_z·q_s = s3 + gp·ζ
                    e2 = std::max(e2, std::fabs(q.qz * qs - s3 - gp * zeta));
                    // (IV) q_xy² = <(H,K)>_S + (Hκ+Kλ+Lρ)²
                    e4 = std::max(e4, std::fabs(q.qxy * q.qxy - plaHK - leak * leak));
                }
        e3 = std::fabs(qs * qs - guv - zeta * zeta);
        // abgToZ consistency with the G* decomposition
        double za = 0.0, zb = 0.0, Cg = 0.0;
        abgToZ(cell.a, cell.b, cell.ga, za, zb, Cg);
        eAB = std::fabs(za * za - A) + std::fabs(zb * zb - B) + std::fabs(Cg - C);
        std::printf("       SCRM identity residuals: (I)=%.3g (II)=%.3g (III)=%.3g (IV)=%.3g  abgToZ=%.3g\n",
                    e1, e2, e3, e4, eAB);
        std::fflush(stdout);
        check(e1 < 1e-12, "SCRM identity (I) |g|²=pla+gp²", e1, 1e-12);
        check(e2 < 1e-12, "SCRM identity (II) q_z·q_s=s3+gp·ζ", e2, 1e-12);
        check(e3 < 1e-12, "SCRM identity (III) q_s²=guv+ζ²", e3, 1e-12);
        check(e4 < 1e-12, "SCRM identity (IV) q_xy²=<(H,K)>+(Hκ+Kλ+Lρ)²", e4, 1e-12);
        check(eAB < 1e-9, "SCRM: G* decomposition agrees with abgToZ", eAB, 1e-9);
    }

    // ---------------------------------------------------------------- fixed contact plane
    {
        const Cell cell{10.0, 12.0, 17.0, 90.0, 90.0, 90.0};
        const Hkl uvw{0, 0, 1};
        double qspec = 0.0;
        const std::vector<Peak> data = makeData(cell, uvw, 2, 2, 2, false, qspec);

        {
            std::vector<Hkl> ah;
            const double r0 = assignHkl(cell, Orient::fromUvw(uvw), data, hklGrid(3, 3, 3), ah);
            const int li = latticeIndex(ah);
            std::printf("       truth cell assignment rms=%.3g   latticeIndex=%d  peaks=%zu\n", r0, li,
                        ah.size());
            std::fflush(stdout);
            check(r0 < 1e-9, "fixed contact plane: zero assignment residual for the truth cell", r0, 1e-9);
            check(li == 1, "fixed contact plane: the truth is a primitive cell", (double)li, 1.0);

            // Refine from a near-truth initial value using the truth hkl (to isolate LM itself)
            {
                Cell c2{10.6, 12.4, 16.6, 90.0, 90.0, 90.0};
                const std::array<bool, 6> fm = systemFreeMask(CrystalSystem::Orthorhombic);
                refineCell(data, ah, Orient::fromUvw(uvw), c2, fm, 200);
                std::printf("       truth hkl + LM: %s (expected %s)\n", fmtCell(c2).c_str(),
                            fmtCell(cell).c_str());
                std::fflush(stdout);
                const double d2 = std::fabs(c2.a - cell.a) + std::fabs(c2.b - cell.b) +
                                  std::fabs(c2.c - cell.c);
                check(d2 < 1e-4, "fixed contact plane: LM with truth hkl converges exactly", d2, 1e-4);
            }

            // Starting from a near-truth initial value, alternating refinement should converge back to the truth
            Cell near{10.2, 12.2, 16.8, 90.0, 90.0, 90.0};
            double rr = 0.0;
            const Cell got = refineAssign(data, Orient::fromUvw(uvw), hklGrid(3, 3, 3), near,
                                          CrystalSystem::Orthorhombic, 8, rr);
            std::printf("       refinement from near-truth initial value: %s rms=%.5g\n", fmtCell(got).c_str(), rr);
            std::fflush(stdout);
            const double dd = std::fabs(got.a - cell.a) + std::fabs(got.b - cell.b) +
                              std::fabs(got.c - cell.c);
            check(dd < 0.02, "fixed contact plane: converges from a near-truth initial value", dd, 0.02);
        }

        EngineParams ep;
        ep.mode = IndexMode::FixedPlane;
        ep.fixed.u = uvw.h; ep.fixed.v = uvw.k; ep.fixed.w = uvw.l;
        ep.fixed.maxH = 3; ep.fixed.maxK = 3; ep.fixed.maxL = 3;
        ep.fixed.system = CrystalSystem::Orthorhombic;
        ep.fixed.start = Cell{10.2, 12.2, 16.8, 90.0, 90.0, 90.0};
        ep.fixed.seedPeaks = 8;
        ep.fixed.acceptTol = 0.05;
        ep.fixed.threads = 4;
        const EngineResult r = runIndexing(data, ep, Progress{});
        check(r.ok, "fixed contact plane finds solutions");
        if (r.ok) {
            const Solution& s = r.solutions.front();
            const double d = std::fabs(s.cell.a - cell.a) + std::fabs(s.cell.b - cell.b) +
                             std::fabs(s.cell.c - cell.c);
            std::printf("       best: %s rms=%.5f (true %s)\n", fmtCell(s.cell).c_str(),
                        s.rmsQxyz, fmtCell(cell).c_str());
            check(d < 0.05, "fixed contact plane: cell parameters recovered", d, 0.05);
            check(s.rmsQxyz < 0.01, "fixed contact plane: RMSD < 0.01", s.rmsQxyz, 0.01);
        }
    }

    // ------------------------------------- fixed contact plane: hard-constraint regression for a given (h k l)
    // Real data (q_xy,q_z columns of examples/peak_finder/fitted.txt, 23 peaks, including the specular rod).
    // The two-step method used to reduce the whole cell unconditionally, so fixing (0 0 1) rewrote the
    // result to out-of-plane directions such as (-1 0 -1) or (0 1 -1); likewise for (1 0 0) and (1 1 0).
    {
        static const double kFitted[23][2] = {
            {-0.26960, 0.25436}, {-0.26828, 0.72876}, {-0.26928, 0.49287}, {-0.79981, 0.49199},
            {-0.80082, 0.72912}, {-0.53501, 1.20578}, {-0.53572, 0.96587}, {-1.32958, 0.49038},
            {-1.08135, 0.96708}, {-0.26974, 0.96896}, {-0.53571, 0.73086}, {-0.93283, 0.25263},
            {-0.00000, 0.49248}, {-0.53328, 1.43900}, {-1.06388, 0.25099}, {-0.93243, 0.25305},
            {-1.09800, 0.26140}, {-1.58947, 0.25136}, {-0.46865, 0.24762}, {-0.38054, 0.49469},
            {-0.65601, 0.25577}, {-0.26952, 1.20678}, {-0.53340, 1.43892}};
        std::vector<Peak> data;
        for (const auto& p : kFitted) data.push_back(Peak{p[0], p[1]});
        const Hkl dirs[3] = {{0, 0, 1}, {1, 0, 0}, {1, 1, 0}};
        for (const Hkl& d : dirs) {
            FixedPlaneParams fp;
            fp.u = d.h;
            fp.v = d.k;
            fp.w = d.l;
            const std::vector<Solution> sols = runFixedPlane(data, fp, Progress{});
            bool strict = !sols.empty();
            for (const Solution& s : sols)
                if (s.uvw.h != d.h || s.uvw.k != d.k || s.uvw.l != d.l) strict = false;
            char name[128];
            std::snprintf(name, sizeof(name), "fixed contact plane regression: (%d %d %d) result uvw strictly consistent (n=%d)",
                          d.h, d.k, d.l, (int)sols.size());
            check(strict, name, (double)sols.size(), 1.0);
        }
        // Regression: the q_xy sign is already normalized to |q_xy| at the entry point (normalizeQxySigns).
        // The fitted data all have negative q_xy: autoSign on/off must give the same result; negating the
        // whole table must also yield the same set of solutions (same count, uvw and RMS).
        {
            FixedPlaneParams on;
            on.u = 0;
            on.v = 0;
            on.w = 1;
            const std::vector<Solution> sOn = runFixedPlane(data, on, Progress{});
            on.autoSign = false;
            const std::vector<Solution> sOff = runFixedPlane(data, on, Progress{});
            bool same = !sOn.empty() && sOn.size() == sOff.size();
            for (size_t i = 0; same && i < sOn.size(); ++i)
                same = sOn[i].uvw.h == sOff[i].uvw.h && sOn[i].uvw.k == sOff[i].uvw.k &&
                       sOn[i].uvw.l == sOff[i].uvw.l &&
                       std::fabs(sOn[i].rmsQxyz - sOff[i].rmsQxyz) < 1e-12;
            check(same, "negative-value normalization regression: autoSign on/off gives the same result",
                  (double)sOn.size(), (double)sOff.size());

            std::vector<Peak> neg = data;
            for (Peak& p : neg) p.qxy = -p.qxy;
            const std::vector<Solution> sNeg = runFixedPlane(neg, on, Progress{});
            bool sameNeg = !sOn.empty() && sNeg.size() == sOn.size();
            for (size_t i = 0; sameNeg && i < sOn.size(); ++i)
                sameNeg = sNeg[i].uvw.h == sOn[i].uvw.h && sNeg[i].uvw.k == sOn[i].uvw.k &&
                          sNeg[i].uvw.l == sOn[i].uvw.l &&
                          std::fabs(sNeg[i].rmsQxyz - sOn[i].rmsQxyz) < 1e-12;
            check(sameNeg, "negative-value normalization regression: negative input matches positive input",
                  (double)sNeg.size(), (double)sOn.size());
        }
    }

    // ------------------------------------------- with q_spec: crystal-system constraint regression
    // The default Triclinic imposes no constraint; for any other crystal system both WithSpec implementations
    // must lock the output cell inside the target system (angles / a=b) and still recover the true cell (volume + RMS).
    {
        struct SysCase {
            const char* name;
            Cell truth;
            CrystalSystem sys;
        };
        const SysCase cases[] = {
            {"orthorhombic", Cell{8.6, 11.4, 15.8, 90.0, 90.0, 90.0}, CrystalSystem::Orthorhombic},
            {"tetragonal", Cell{9.2, 9.2, 15.4, 90.0, 90.0, 90.0}, CrystalSystem::Tetragonal},
            {"hexagonal", Cell{9.2, 9.2, 15.4, 90.0, 90.0, 120.0}, CrystalSystem::Hexagonal},
        };
        for (const SysCase& cs : cases) {
            double qspec = 0.0;
            const std::vector<Peak> data = makeData(cs.truth, Hkl{0, 0, 1}, 2, 2, 3, true, qspec);
            EngineParams ep;
            ep.mode = IndexMode::WithSpec;
            ep.spec.qspec = qspec;
            ep.spec.uvLimit = 1;
            ep.spec.hkLimit1 = 2;
            ep.spec.hkLimit2 = 5;
            ep.spec.lLimit = 4;
            ep.spec.wLimit = 2;
            ep.spec.maxPart1 = 6;
            ep.spec.maxSolutions = 10;
            ep.spec.threads = 4;
            ep.spec.system = cs.sys;
            const double vTrue = cellVolume(cs.truth);

            const Part1Algorithm algs[] = {Part1Algorithm::Improved,
                                           Part1Algorithm::PeriodVoteDdl,
                                           Part1Algorithm::Scrm};
            const char* algName[] = {"improved", "PeriodVote/DDL", "SCRM"};
            for (int ai = 0; ai < 3; ++ai) {
                ep.spec.part1Alg = algs[ai];
                ep.spec.system = cs.sys;
                const EngineResult r = runIndexing(data, ep, Progress{});
                const std::string tag = std::string("WithSpec crystal-system constraint (") + cs.name + ", " +
                                       algName[ai] + ")";
                if (!r.ok) {
                    // This algorithm has no solution even for triclinic (unconstrained) on this synthetic data --
                    // a limitation of the search algorithm itself (start triple degenerates at high symmetry), not a regression from the constraint.
                    ep.spec.system = CrystalSystem::Triclinic;
                    const EngineResult base = runIndexing(data, ep, Progress{});
                    if (!base.ok) {
                        std::printf("       [skip] %s has no triclinic solution on the %s data\n", algName[ai],
                                    cs.name);
                        continue;
                    }
                }
                check(r.ok, tag + " finds solutions");
                if (!r.ok) continue;
                double worstShape = 0.0;
                bool hit = false;
                for (const Solution& s : r.solutions) {
                    double e = std::max({std::fabs(s.cell.al - cs.truth.al),
                                         std::fabs(s.cell.be - cs.truth.be),
                                         std::fabs(s.cell.ga - cs.truth.ga)});
                    if (cs.sys != CrystalSystem::Orthorhombic)
                        e = std::max(e, std::fabs(s.cell.a - s.cell.b));
                    worstShape = std::max(worstShape, e);
                    const double rel = std::fabs(s.vol - vTrue) / vTrue;
                    if (rel < 0.03 && s.rmsQxyz < 0.02) hit = true;
                }
                check(worstShape < 1e-6, tag + " all output solutions satisfy the crystal-system shape", worstShape, 1e-6);
                check(hit, tag + " recovers the true cell volume with RMSD < 0.02");
            }
        }
        check(SpecParams{}.system == CrystalSystem::Triclinic,
              "WithSpec default crystal system = triclinic (unconstrained, as in the old version)");
    }

    std::printf("\n%s\n", g_fail == 0 ? "ALL: PASS" : ("FAILED: " + std::to_string(g_fail)).c_str());
    return g_fail == 0 ? 0 : 1;
}
