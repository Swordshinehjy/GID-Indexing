// Temporary probe: use the 12 peak positions of dataset_DIP to check the manual solution against Part ONE / Part TWO / the fixed contact plane
#include "engine.h"
#include "fixedplane.h"
#include "refine.h"
#include "specindex.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace gid;

static const std::vector<Peak> kDip = {
    {0.000, 1.776}, {0.368, 0.091}, {0.396, 1.160}, {0.428, 1.069}, {0.689, 0.984},
    {0.734, 0.182}, {0.806, 1.232}, {1.140, 0.171}, {1.223, 0.083}, {1.529, 0.790},
    {1.625, 0.697}, {2.274, 0.953},
};

static void printSol(const Solution& s) {
    std::printf("   uvw=(%d %d %d) a=%.3f b=%.3f c=%.3f al=%.2f be=%.2f ga=%.2f V=%.1f "
                "rmsQxyz=%.4f rmsQxy=%.4f rmsQz=%.4f qspec=%.4f\n",
                s.uvw.h, s.uvw.k, s.uvw.l, s.cell.a, s.cell.b, s.cell.c, s.cell.al, s.cell.be,
                s.cell.ga, s.vol, s.rmsQxyz, s.rmsQxy, s.rmsQz, s.qspec);
}

int main() {
    double qs = 0.0;
    int rods = 0;
    detectQSpec(kDip, qs, rods);
    std::printf("qspec=%.4f rods=%d\n", qs, rods);

    // ---- Can the final manual solution be fitted under our model?
    {
        Cell mc{7.13, 8.48, 16.67, 89.4, 87.8, 89.7};
        Hkl uvw{-1, 2, 1};
        std::vector<Hkl> ah;
        const std::vector<Hkl> grid = hklGrid(6, 6, 6);
        assignHkl(mc, Orient::fromUvw(uvw), kDip, grid, ah);
        const Quality q = evaluate(mc, Orient::fromUvw(uvw), kDip, ah);
        std::printf("manual cell: rmsQxyz=%.5f rmsQxy=%.5f rmsQz=%.5f V=%.2f qspec=%.4f\n",
                    q.rmsQxyz, q.rmsQxy, q.rmsQz, cellVolume(mc), qSpecOf(reciprocalMatrix(mc), uvw));
    }

    // ---- Verify the original formula against the 3x3 solution of the tilted contact plane (u=-1,v=2)
    {
        const double u = -1.0, v = 2.0;
        struct Row { double Q, Z, h, k; };
        const Row P[3] = {{0.428, 1.069, -1, 1}, {0.689, 0.984, -1, 1}, {1.223, 0.083, 1, 1}};
        double M[3][3], b[3] = {1.0, 1.0, 1.0};
        for (int i = 0; i < 3; ++i) {
            const double H = P[i].h - u * P[i].Z / qs;
            const double K = P[i].k - v * P[i].Z / qs;
            const double iQ2 = 1.0 / (P[i].Q * P[i].Q);
            M[i][0] = H * H * iQ2 + u * u / (qs * qs);
            M[i][1] = K * K * iQ2 + v * v / (qs * qs);
            M[i][2] = -2.0 * (H * K * iQ2 + u * v / (qs * qs));
        }
        auto det3 = [](double m[3][3]) {
            return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                   m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        };
        const double D = det3(M);
        double X = 0, Y = 0, Z = 0;
        for (int c = 0; c < 3; ++c) {
            double T[3][3];
            for (int r = 0; r < 3; ++r)
                for (int cc = 0; cc < 3; ++cc) T[r][cc] = (cc == c) ? b[r] : M[r][cc];
            const double dc = det3(T);
            if (c == 0) X = dc / D; else if (c == 1) Y = dc / D; else Z = dc / D;
        }
        const double den = qs * qs - u * u * X - v * v * Y + 2.0 * u * v * Z;
        const double Lam = (X * Y - Z * Z) / den;
        const double za = std::sqrt(X + v * v * Lam), zb = std::sqrt(Y + u * u * Lam);
        const double cg = (Z + u * v * Lam) / (za * zb);
        const double sg = std::sqrt(std::max(1e-12, 1.0 - cg * cg));
        const double a = 6.283185307179586 / (za * sg), bb = 6.283185307179586 / (zb * sg);
        const double ga = std::acos(cg) * 180.0 / 3.14159265358979;
        std::printf("3x3 check uv=(-1,2): a=%.4f b=%.4f ga=%.3f\n", a, bb, ga);
    }

    // ---- Does the full Part ONE candidate set contain the manual in-plane cell
    SpecParams sp;
    sp.qspec = qs;
    sp.uvLimit = 2;
    sp.linesToPermute = 7; // include (1.223,0.083) in the seed peaks
    sp.hkLimit1 = 1;       // in-plane indices ±1 suffice (DIP only needs -1,0,1)
    sp.hkLimit2 = 4;
    sp.maxPart1 = 1000000; // no truncation, to see the true ranking
    sp.threads = -1;
    std::vector<AbgCand> c;
    {
        const auto t0 = std::chrono::steady_clock::now();
        c = improved::runPartOne(kDip, sp, Progress{});
        const double ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::printf("\nPartOne(all): %zu cands (%.0f ms)\n", c.size(), ms);
    }
    int bestHit = -1;
    {
        int hits = 0;
        for (size_t i = 0; i < c.size(); ++i) {
            const AbgCand& x = c[i];
            const double lo = std::min(x.a, x.b), hi = std::max(x.a, x.b);
            if (lo > 6.4 && lo < 8.0 && hi > 7.8 && hi < 9.6 && std::fabs(x.ga - 90.0) < 8.0) {
                if (bestHit < 0) bestHit = (int)i;
                std::printf("  HIT rank=%zu uv=(%d,%d) a=%.3f b=%.3f ga=%.2f rmsQxy=%.5f\n", i, x.u,
                            x.v, x.a, x.b, x.ga, x.rmsQxy);
                if (++hits >= 8) break;
            }
        }
        std::printf("hits(first 8)=%d bestHit=%d\n", hits, bestHit);
    }
    {
        int cnt[5][5] = {};
        for (const AbgCand& x : c)
            if (x.u >= -2 && x.u <= 2 && x.v >= -2 && x.v <= 2) ++cnt[x.u + 2][x.v + 2];
        for (int u = -2; u <= 2; ++u) {
            for (int v = -2; v <= 2; ++v) std::printf("  cnt uv=(%d,%d)=%d", u, v, cnt[u + 2][v + 2]);
            std::printf("\n");
        }
        const AbgCand* bestUV = nullptr;
        for (const AbgCand& x : c)
            if (x.u == -1 && x.v == 2 && (!bestUV || x.rmsQxy < bestUV->rmsQxy)) bestUV = &x;
        if (bestUV)
            std::printf("uv=(-1,2) best-rms: a=%.3f b=%.3f ga=%.2f rmsQxy=%.5f\n", bestUV->a,
                        bestUV->b, bestUV->ga, bestUV->rmsQxy);
        else
            std::printf("uv=(-1,2): NO candidate\n");
    }
    {
        // Candidate globally closest to the manual in-plane cell (no sorting, direct scan)
        double bestScore = 1e9;
        int bi = -1;
        for (size_t i = 0; i < c.size(); ++i) {
            const double s = std::fabs(c[i].a - 7.13) + std::fabs(c[i].b - 8.48) +
                             std::fabs(c[i].ga - 89.7) / 5.0;
            if (s < bestScore) { bestScore = s; bi = (int)i; }
        }
        if (bi >= 0) {
            std::printf("nearest manual-cell cand: uv=(%d,%d) a=%.3f b=%.3f ga=%.2f rmsQxy=%.5f\n",
                        c[bi].u, c[bi].v, c[bi].a, c[bi].b, c[bi].ga, c[bi].rmsQxy);
            if (bestScore < 0.6) bestHit = bi;
        }
    }

    // ---- Part TWO: use the hit candidate directly (no other parameters changed)
    if (bestHit >= 0) {
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<Solution> ss = improved::runPartTwo(kDip, c[bestHit], sp);
        const double ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::printf("\nPartTwo(all-on-hit): %zu sols (%.0f ms)\n", ss.size(), ms);
        for (size_t i = 0; i < ss.size() && i < 8; ++i) printSol(ss[i]);
    }

    // ---- Deterministic refinement starting from the manual cell (step 1 inside the fixed contact plane)
    {
        Cell st{7.13, 8.48, 16.67, 89.4, 87.8, 89.7};
        const Orient ori = Orient::fromUvw(Hkl{-1, 2, 1});
        const std::vector<Hkl> cands = hklGrid(4, 4, 4);
        double rms = 0.0;
        Cell out = refineAssign(kDip, ori, cands, st, CrystalSystem::Triclinic, 6, rms);
        std::vector<Hkl> ah;
        assignHkl(out, ori, kDip, cands, ah);
        const Quality q = evaluate(out, ori, kDip, ah);
        std::printf("\nrefineAssign from manual: a=%.3f b=%.3f c=%.3f al=%.2f be=%.2f ga=%.2f "
                    "rms=%.5f | eval rmsQxyz=%.5f V=%.1f qspec=%.4f\n",
                    out.a, out.b, out.c, out.al, out.be, out.ga, rms, q.rmsQxyz,
                    cellVolume(out), qSpecOf(reciprocalMatrix(out), Hkl{-1, 2, 1}));
    }

    // ---- Fixed contact plane with explicit (-1,2,1): different initial values
    {
        const Cell starts[] = {
            Cell{7.13, 8.48, 16.67, 89.4, 87.8, 89.7},
            Cell{7.5, 8.5, 16.5, 90.0, 90.0, 90.0},
            Cell{10.0, 10.0, 15.0, 90.0, 90.0, 90.0},
        };
        for (const Cell& st : starts) {
            FixedPlaneParams fp;
            fp.u = -1; fp.v = 2; fp.w = 1;
            fp.maxH = fp.maxK = fp.maxL = 6;
            fp.start = st;
            fp.seedPeaks = 8;
            fp.threads = -1;
            const auto t0 = std::chrono::steady_clock::now();
            std::vector<Solution> ss = runFixedPlane(kDip, fp, Progress{});
            const double ms = std::chrono::duration<double, std::milli>(
                                  std::chrono::steady_clock::now() - t0).count();
            std::printf("\nFixedPlane(-1,2,1) start a=%.3f b=%.3f c=%.3f: %zu sols (%.0f ms)\n",
                        st.a, st.b, st.c, ss.size(), ms);
            for (size_t i = 0; i < ss.size() && i < 10; ++i) printSol(ss[i]);
        }
    }
    return 0;
}
