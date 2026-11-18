// Temporary probe: behavior of the fixed contact plane on DIP
#include "fixedplane.h"
#include "refine.h"

#include <chrono>
#include <cstdio>
#include <vector>

using namespace gid;
using clk = std::chrono::steady_clock;

static const std::vector<Peak> kDip = {
    {0.000, 1.776}, {0.368, 0.091}, {0.396, 1.160}, {0.428, 1.069}, {0.689, 0.984},
    {0.734, 0.182}, {0.806, 1.232}, {1.140, 0.171}, {1.223, 0.083}, {1.529, 0.790},
    {1.625, 0.697}, {2.274, 0.953},
};

static void show(const char* tag, const std::vector<Solution>& s, double ms) {
    std::printf("%-34s %7.0f ms  %2zu sols\n", tag, ms, s.size());
    for (size_t i = 0; i < s.size() && i < 6; ++i) {
        const Solution& x = s[i];
        std::printf("    uvw=(%d %d %d) a=%.3f b=%.3f c=%.3f al=%.1f be=%.1f ga=%.1f V=%.0f "
                    "rms=%.5f\n",
                    x.uvw.h, x.uvw.k, x.uvw.l, x.cell.a, x.cell.b, x.cell.c, x.cell.al, x.cell.be,
                    x.cell.ga, x.vol, x.rmsQxyz);
    }
}

int main() {
    // Ground-truth solution (produced by the improved version): a=7.157 b=8.491 c=16.696 al=89.1 be=87.1 ga=89.6 rms=0.00394
    const Cell truth{7.1565, 8.4910, 16.6963, 89.142, 87.107, 89.622};
    const Hkl uvw{-1, 2, 1};

    // 0) Assignment residual of the truth cell on the h2/k2/l6 grid
    {
        std::vector<Hkl> ah;
        std::vector<double> dist;
        const double r = assignHkl(truth, Orient::fromUvw(uvw), kDip, hklGrid(2, 2, 6), ah, &dist);
        const Quality q = evaluate(truth, Orient::fromUvw(uvw), kDip, ah);
        std::printf("truth cell assignment rms=%.5f (rmsQxy=%.5f rmsQz=%.5f)\n", r, q.rmsQxy, q.rmsQz);
        for (size_t i = 0; i < kDip.size(); ++i)
            std::printf("    peak%2zu (%.3f,%.3f) -> (%d %d %d) d=%.4f\n", i, kDip[i].qxy,
                        kDip[i].qz, ah[i].h, ah[i].k, ah[i].l, i < dist.size() ? dist[i] : -1.0);
    }

    // 1) Default parameters
    {
        FixedPlaneParams fp;
        fp.u = -1; fp.v = 2; fp.w = 1;
        fp.maxH = 2; fp.maxK = 2; fp.maxL = 6;
        fp.seedPeaks = 8;
        auto a = clk::now();
        const std::vector<Solution> s = runFixedPlane(kDip, fp, Progress{});
        show("default (start 10,10,15,90,90,90)", s,
             std::chrono::duration<double, std::milli>(clk::now() - a).count());
    }

    // 2) Initial value taken near the truth
    {
        FixedPlaneParams fp;
        fp.u = -1; fp.v = 2; fp.w = 1;
        fp.maxH = 2; fp.maxK = 2; fp.maxL = 6;
        fp.start = Cell{7.5, 8.8, 17.0, 92.0, 90.0, 88.0};
        fp.seedPeaks = 8;
        auto a = clk::now();
        const std::vector<Solution> s = runFixedPlane(kDip, fp, Progress{});
        show("initial value near truth", s, std::chrono::duration<double, std::milli>(clk::now() - a).count());
    }

    // 2b) Initial value derived from |q| (three smallest d = 2π/|q|)
    {
        FixedPlaneParams fp;
        fp.u = -1; fp.v = 2; fp.w = 1;
        fp.maxH = 2; fp.maxK = 2; fp.maxL = 6;
        fp.start = Cell{16.58, 8.31, 5.45, 90.0, 90.0, 90.0};
        fp.seedPeaks = 8;
        auto a = clk::now();
        const std::vector<Solution> s = runFixedPlane(kDip, fp, Progress{});
        show("initial value = three smallest d", s, std::chrono::duration<double, std::milli>(clk::now() - a).count());
    }

    // 2c) Correct lengths, angles 90
    {
        FixedPlaneParams fp;
        fp.u = -1; fp.v = 2; fp.w = 1;
        fp.maxH = 2; fp.maxK = 2; fp.maxL = 6;
        fp.start = Cell{7.16, 8.49, 16.70, 90.0, 90.0, 90.0};
        fp.seedPeaks = 8;
        auto a = clk::now();
        const std::vector<Solution> s = runFixedPlane(kDip, fp, Progress{});
        show("initial value = true lengths/90 deg", s,
             std::chrono::duration<double, std::milli>(clk::now() - a).count());
    }

    // 3) LM directly with the truth hkl from a near-truth initial value
    {
        std::vector<Hkl> ah;
        assignHkl(truth, Orient::fromUvw(uvw), kDip, hklGrid(2, 2, 6), ah);
        Cell c{7.5, 8.8, 17.0, 92.0, 90.0, 88.0};
        const std::array<bool, 6> fm{true, true, true, true, true, true};
        const bool ok = refineCell(kDip, ah, Orient::fromUvw(uvw), c, fm, 200);
        const double r = assignHkl(c, Orient::fromUvw(uvw), kDip, hklGrid(2, 2, 6), ah);
        std::printf("truth hkl + LM: ok=%d a=%.3f b=%.3f c=%.3f al=%.1f be=%.1f ga=%.1f rms=%.5f\n",
                    (int)ok, c.a, c.b, c.c, c.al, c.be, c.ga, r);
    }

    // 4) No specular peak (drop the (0, 1.776) rod peak) -- the FP-DL deterministic search is mandatory
    {
        const std::vector<Peak> noRod(kDip.begin() + 1, kDip.end());
        const Cell starts[] = {
            Cell{10.0, 10.0, 15.0, 90.0, 90.0, 90.0},
            Cell{16.58, 8.31, 5.45, 90.0, 90.0, 90.0},
            Cell{12.5, 9.5, 20.0, 85.0, 95.0, 75.0},
        };
        for (const Cell& st : starts) {
            FixedPlaneParams fp;
            fp.u = -1; fp.v = 2; fp.w = 1;
            fp.maxH = 2; fp.maxK = 2; fp.maxL = 6;
            fp.start = st;
            auto a = clk::now();
            const std::vector<Solution> s = runFixedPlane(noRod, fp, Progress{});
            char tag[96];
            std::snprintf(tag, sizeof(tag), "no rod peak start %.1f/%.1f/%.1f", st.a, st.b, st.c);
            show(tag, s, std::chrono::duration<double, std::milli>(clk::now() - a).count());
        }
    }
    return 0;
}
