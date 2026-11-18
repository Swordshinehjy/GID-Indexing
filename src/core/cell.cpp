// ============================================================================
//  cell.cpp — see cell.h
// ============================================================================
#include "cell.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <numeric>

namespace gid {
namespace {

// ---- 3x3 row-major helpers --------------------------------------------------
double det3(const Mat3& m) {
    return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) +
           m[2] * (m[3] * m[7] - m[4] * m[6]);
}

bool inv3(const Mat3& m, Mat3& out) {
    const double d = det3(m);
    if (!(std::fabs(d) > 1e-300)) return false;
    const double id = 1.0 / d;
    out[0] = (m[4] * m[8] - m[5] * m[7]) * id;
    out[1] = (m[2] * m[7] - m[1] * m[8]) * id;
    out[2] = (m[1] * m[5] - m[2] * m[4]) * id;
    out[3] = (m[5] * m[6] - m[3] * m[8]) * id;
    out[4] = (m[0] * m[8] - m[2] * m[6]) * id;
    out[5] = (m[2] * m[3] - m[0] * m[5]) * id;
    out[6] = (m[3] * m[7] - m[4] * m[6]) * id;
    out[7] = (m[1] * m[6] - m[0] * m[7]) * id;
    out[8] = (m[0] * m[4] - m[1] * m[3]) * id;
    return true;
}

Mat3 mul3(const Mat3& A, const Mat3& B) {
    Mat3 R{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int t = 0; t < 3; ++t) s += A[3 * i + t] * B[3 * t + j];
            R[3 * i + j] = s;
        }
    return R;
}

// Quadratic form of a symmetric 3x3: v^T S v, with S stored as [G11,G22,G33,G12,G13,G23]
double quadForm(const std::array<double, 6>& g, const std::array<double, 3>& v) {
    return g[0] * v[0] * v[0] + g[1] * v[1] * v[1] + g[2] * v[2] * v[2] +
           2.0 * (g[3] * v[0] * v[1] + g[4] * v[0] * v[2] + g[5] * v[1] * v[2]);
}

double clamp1(double x) { return x > 1.0 ? 1.0 : (x < -1.0 ? -1.0 : x); }

} // namespace

// ---------------------------------------------------------------- Names
std::vector<std::string> crystalSystemNames() {
    return {"Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal",
            "Hexagonal", "Trigonal",  "Cubic"};
}

std::string crystalSystemName(CrystalSystem s) {
    switch (s) {
        case CrystalSystem::Triclinic: return "Triclinic";
        case CrystalSystem::Monoclinic: return "Monoclinic";
        case CrystalSystem::Orthorhombic: return "Orthorhombic";
        case CrystalSystem::Tetragonal: return "Tetragonal";
        case CrystalSystem::Hexagonal: return "Hexagonal";
        case CrystalSystem::Trigonal: return "Trigonal";
        case CrystalSystem::Cubic: return "Cubic";
    }
    return "Triclinic";
}

// ---------------------------------------------------------------- Basis vectors and metric
Mat3 directMatrix(const Cell& c) {
    const double ca = std::cos(c.al * RAD), cb = std::cos(c.be * RAD), cg = std::cos(c.ga * RAD);
    const double sg = std::sin(c.ga * RAD);
    const double S = std::sqrt(std::max(1e-18, 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg));
    const double s = std::max(1e-9, sg);
    Mat3 A{};
    A[0] = c.a;  A[1] = c.b * cg;  A[2] = c.c * cb;
    A[3] = 0.0;  A[4] = c.b * s;   A[5] = c.c * (ca - cb * cg) / s;
    A[6] = 0.0;  A[7] = 0.0;       A[8] = c.c * S / s;
    return A;
}

Mat3 reciprocalMatrix(const Cell& c) {
    Mat3 A = directMatrix(c);
    Mat3 Ai;
    if (!inv3(A, Ai)) return Mat3{};
    Mat3 B{};
    for (int i = 0; i < 9; ++i) B[i] = TWO_PI * Ai[i];
    return B;
}

std::array<double, 6> recipMetric(const Cell& c) {
    const Mat3 B = reciprocalMatrix(c);
    std::array<double, 6> g{};
    g[0] = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
    g[1] = B[3] * B[3] + B[4] * B[4] + B[5] * B[5];
    g[2] = B[6] * B[6] + B[7] * B[7] + B[8] * B[8];
    g[3] = B[0] * B[3] + B[1] * B[4] + B[2] * B[5];
    g[4] = B[0] * B[6] + B[1] * B[7] + B[2] * B[8];
    g[5] = B[3] * B[6] + B[4] * B[7] + B[5] * B[8];
    return g;
}

std::array<double, 6> directMetric(const Cell& c) {
    std::array<double, 6> g{};
    g[0] = c.a * c.a;
    g[1] = c.b * c.b;
    g[2] = c.c * c.c;
    g[3] = c.a * c.b * std::cos(c.ga * RAD);
    g[4] = c.a * c.c * std::cos(c.be * RAD);
    g[5] = c.b * c.c * std::cos(c.al * RAD);
    return g;
}

bool cellFromDirectMetric(const std::array<double, 6>& gd, Cell& out) {
    if (gd[0] <= 0 || gd[1] <= 0 || gd[2] <= 0) return false;
    const double a = std::sqrt(gd[0]), b = std::sqrt(gd[1]), c = std::sqrt(gd[2]);
    const double cg = clamp1(gd[3] / (a * b));
    const double cb = clamp1(gd[4] / (a * c));
    const double ca = clamp1(gd[5] / (b * c));
    out.a = a; out.b = b; out.c = c;
    out.ga = std::acos(cg) * DEG;
    out.be = std::acos(cb) * DEG;
    out.al = std::acos(ca) * DEG;
    return std::isfinite(out.a) && std::isfinite(out.al);
}

bool cellFromRecipMetric(const std::array<double, 6>& gs, Cell& out) {
    // G*^-1
    const double d = gs[0] * (gs[1] * gs[2] - gs[5] * gs[5]) -
                     gs[3] * (gs[3] * gs[2] - gs[5] * gs[4]) +
                     gs[4] * (gs[3] * gs[5] - gs[1] * gs[4]);
    if (!(std::fabs(d) > 1e-30)) return false;
    const double id = 1.0 / d;
    const double i11 = (gs[1] * gs[2] - gs[5] * gs[5]) * id;
    const double i22 = (gs[0] * gs[2] - gs[4] * gs[4]) * id;
    const double i33 = (gs[0] * gs[1] - gs[3] * gs[3]) * id;
    const double i12 = (gs[4] * gs[5] - gs[3] * gs[2]) * id;
    const double i13 = (gs[3] * gs[5] - gs[1] * gs[4]) * id;
    const double i23 = (gs[3] * gs[4] - gs[0] * gs[5]) * id;
    const double s = TWO_PI * TWO_PI;
    std::array<double, 6> gd{s * i11, s * i22, s * i33, s * i12, s * i13, s * i23};
    return cellFromDirectMetric(gd, out);
}

double cellVolume(const Cell& c) {
    const double ca = std::cos(c.al * RAD), cb = std::cos(c.be * RAD), cg = std::cos(c.ga * RAD);
    return c.a * c.b * c.c *
           std::sqrt(std::max(0.0, 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg));
}

double volumeFromRecipDet(double detGs) {
    if (!(detGs > 0.0)) return 0.0;
    const double s = TWO_PI * TWO_PI * TWO_PI;
    return s / std::sqrt(detGs);
}

// ---------------------------------------------------------------- Forward model
QPoint forwardHkl(const Mat3& recip, const Hkl& uvw, int h, int k, int l) {
    const double u = uvw.h, v = uvw.k, w = uvw.l;
    double gu[3] = {u * recip[0] + v * recip[3] + w * recip[6],
                    u * recip[1] + v * recip[4] + w * recip[7],
                    u * recip[2] + v * recip[5] + w * recip[8]};
    const double gn = std::sqrt(gu[0] * gu[0] + gu[1] * gu[1] + gu[2] * gu[2]);
    QPoint out;
    if (!(gn > 1e-12)) return out;
    const double inv = 1.0 / gn;
    for (int i = 0; i < 3; ++i) gu[i] *= inv;
    const double g[3] = {h * recip[0] + k * recip[3] + l * recip[6],
                         h * recip[1] + k * recip[4] + l * recip[7],
                         h * recip[2] + k * recip[5] + l * recip[8]};
    const double g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
    const double qz = g[0] * gu[0] + g[1] * gu[1] + g[2] * gu[2];
    out.qz = qz;
    out.qxy = std::sqrt(std::max(0.0, g2 - qz * qz));
    return out;
}

QPoint forwardHkl(const Cell& c, const Hkl& uvw, int h, int k, int l) {
    return forwardHkl(reciprocalMatrix(c), uvw, h, k, l);
}

double qSpecOf(const Mat3& recip, const Hkl& uvw) {
    const double u = uvw.h, v = uvw.k, w = uvw.l;
    const double gu[3] = {u * recip[0] + v * recip[3] + w * recip[6],
                          u * recip[1] + v * recip[4] + w * recip[7],
                          u * recip[2] + v * recip[5] + w * recip[8]};
    return std::sqrt(gu[0] * gu[0] + gu[1] * gu[1] + gu[2] * gu[2]);
}

namespace {
// Reciprocal metric G* = B B^T, element order [G11,G22,G33,G12,G13,G23]
std::array<double, 6> recipMetricOf(const Mat3& r) {
    std::array<double, 6> g{};
    g[0] = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    g[1] = r[3] * r[3] + r[4] * r[4] + r[5] * r[5];
    g[2] = r[6] * r[6] + r[7] * r[7] + r[8] * r[8];
    g[3] = r[0] * r[3] + r[1] * r[4] + r[2] * r[5];
    g[4] = r[0] * r[6] + r[1] * r[7] + r[2] * r[8];
    g[5] = r[3] * r[6] + r[4] * r[7] + r[5] * r[8];
    return g;
}
} // namespace

std::array<double, 3> nuFromUvw(const Mat3& recip, const Hkl& uvw) {
    std::array<double, 3> out{};
    const double qs = qSpecOf(recip, uvw);
    if (!(qs > 1e-12)) return out;
    const std::array<double, 3> n{static_cast<double>(uvw.h), static_cast<double>(uvw.k),
                                 static_cast<double>(uvw.l)};
    const std::array<double, 6> g = recipMetricOf(recip);
    // nu = G* n / q_spec
    out[0] = (g[0] * n[0] + g[3] * n[1] + g[4] * n[2]) / qs;
    out[1] = (g[3] * n[0] + g[1] * n[1] + g[5] * n[2]) / qs;
    out[2] = (g[4] * n[0] + g[5] * n[1] + g[2] * n[2]) / qs;
    return out;
}

bool uvwFromNu(const Mat3& recip, const std::array<double, 3>& nu, Hkl& out, int maxIndex) {
    // (u,v,w)/q_spec = (G*)^-1 nu
    const std::array<double, 6> g = recipMetricOf(recip);
    const double d = g[0] * (g[1] * g[2] - g[5] * g[5]) - g[3] * (g[3] * g[2] - g[5] * g[4]) +
                     g[4] * (g[3] * g[5] - g[1] * g[4]);
    if (!(std::fabs(d) > 1e-30)) return false;
    const double id = 1.0 / d;
    const double i11 = (g[1] * g[2] - g[5] * g[5]) * id, i22 = (g[0] * g[2] - g[4] * g[4]) * id,
                 i33 = (g[0] * g[1] - g[3] * g[3]) * id;
    const double i12 = (g[4] * g[5] - g[3] * g[2]) * id, i13 = (g[3] * g[5] - g[1] * g[4]) * id,
                 i23 = (g[3] * g[4] - g[0] * g[5]) * id;
    std::array<double, 3> dir{i11 * nu[0] + i12 * nu[1] + i13 * nu[2],
                              i12 * nu[0] + i22 * nu[1] + i23 * nu[2],
                              i13 * nu[0] + i23 * nu[1] + i33 * nu[2]};
    const double dn = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    if (!(dn > 1e-14)) return false;
    for (int i = 0; i < 3; ++i) dir[i] /= dn;

    int best[3] = {0, 0, 0};
    double bestCos = -2.0;
    for (int i = -maxIndex; i <= maxIndex; ++i)
        for (int j = -maxIndex; j <= maxIndex; ++j)
            for (int k = -maxIndex; k <= maxIndex; ++k) {
                if (i == 0 && j == 0 && k == 0) continue;
                const double n = std::sqrt((double)(i * i + j * j + k * k));
                const double cs = (i * dir[0] + j * dir[1] + k * dir[2]) / n;
                if (std::fabs(cs) > bestCos) {
                    bestCos = std::fabs(cs);
                    best[0] = i; best[1] = j; best[2] = k;
                }
            }
    if (bestCos < 1.0 - 2e-3) return false;
    if (bestCos > 0.0 && (best[0] * dir[0] + best[1] * dir[1] + best[2] * dir[2]) < 0.0) {
        best[0] = -best[0]; best[1] = -best[1]; best[2] = -best[2];
    }
    out = Hkl{best[0], best[1], best[2]};
    // Remove the common factor
    int g0 = std::abs(best[0]), g1 = std::abs(best[1]), g2 = std::abs(best[2]);
    int gc = g0 ? (g1 ? std::gcd(g0, std::gcd(g1, g2 ? g2 : g1)) : (g2 ? std::gcd(g0, g2) : g0)) : 0;
    if (g0 == 0 && g1 == 0) gc = g2;
    if (g0 == 0 && g2 == 0) gc = g1;
    if (g1 == 0 && g2 == 0) gc = g0;
    if (gc > 1) { out.h /= gc; out.k /= gc; out.l /= gc; }
    return true;
}

// ---------------------------------------------------------------- Reparameterization
void abgToZ(double a, double b, double ga, double& za, double& zb, double& C) {
    const double s = std::sin(ga * RAD);
    za = TWO_PI / (a * std::max(1e-9, s));
    zb = TWO_PI / (b * std::max(1e-9, s));
    C = za * zb * std::cos(ga * RAD);
}

void zToAbg(double za, double zb, double C, double& a, double& b, double& ga) {
    const double cg = clamp1(C / (za * zb));
    const double sg = std::sqrt(std::max(1e-18, 1.0 - cg * cg));
    ga = std::acos(cg) * DEG;
    a = TWO_PI / (za * sg);
    b = TWO_PI / (zb * sg);
}

Frame frameOf(const Cell& c) {
    Frame f;
    abgToZ(c.a, c.b, c.ga, f.za, f.zb, f.cg);
    const std::array<double, 6> g = recipMetric(c);
    f.rho = std::sqrt(std::max(1e-18, g[2]));
    f.kap = g[4] / f.rho;
    f.lam = g[5] / f.rho;
    return f;
}

bool cellFromFrame(double za, double zb, double C, double kap, double lam, double rho, Cell& out) {
    if (!(za > 0) || !(zb > 0) || !(rho > 0)) return false;
    // D^2 = kappa^2 zb^2 + 2 kappa lambda C + lambda^2 za^2 + za^2 zb^2 - C^2
    const double D2 = kap * kap * zb * zb + 2.0 * kap * lam * C + lam * lam * za * za +
                      za * za * zb * zb - C * C;
    if (!(D2 > 0.0)) return false;
    const double D = std::sqrt(D2);
    const double denom = std::sqrt(std::max(1e-18, za * za * zb * zb - C * C)); // za zb sin(gamma)
    const double cosAl = clamp1(-(kap * C / za + lam * za) / D);
    const double cosBe = clamp1(-(kap * zb + lam * C / zb) / D);
    zToAbg(za, zb, C, out.a, out.b, out.ga);
    out.c = TWO_PI * D / (rho * denom);
    out.al = std::acos(cosAl) * DEG;
    out.be = std::acos(cosBe) * DEG;
    return std::isfinite(out.c) && out.c > 0;
}

// ---------------------------------------------------------------- Crystal systems
void applySystemConstraints(CrystalSystem sys, Cell& c) {
    // Equal-length parameters take the geometric mean (instead of overwriting one branch with the
    // main one): projection is not biased toward a or b, so a correct solution on a non-default axis setting stays closer.
    auto gm2 = [](double x, double y) {
        return (x > 0.0 && y > 0.0) ? std::sqrt(x * y) : x;
    };
    switch (sys) {
        case CrystalSystem::Cubic: {
            const double m = (c.a > 0.0 && c.b > 0.0 && c.c > 0.0) ? std::cbrt(c.a * c.b * c.c) : c.a;
            c.a = c.b = c.c = m;
            c.al = c.be = c.ga = 90.0;
            break;
        }
        case CrystalSystem::Tetragonal: {
            const double m = gm2(c.a, c.b);
            c.a = c.b = m;
            c.al = c.be = c.ga = 90.0;
            break;
        }
        case CrystalSystem::Hexagonal:
        case CrystalSystem::Trigonal: {
            const double m = gm2(c.a, c.b);
            c.a = c.b = m;
            c.al = c.be = 90.0;
            c.ga = 120.0;
            break;
        }
        case CrystalSystem::Orthorhombic:
            c.al = c.be = c.ga = 90.0; break;
        case CrystalSystem::Monoclinic:
            c.al = c.ga = 90.0; break;
        case CrystalSystem::Triclinic:
        default:
            break;
    }
}

std::array<bool, 6> systemFreeMask(CrystalSystem sys) {
    switch (sys) {
        case CrystalSystem::Cubic: return {true, false, false, false, false, false};
        case CrystalSystem::Tetragonal: return {true, false, true, false, false, false};
        case CrystalSystem::Hexagonal:
        case CrystalSystem::Trigonal: return {true, false, true, false, false, false};
        case CrystalSystem::Orthorhombic: return {true, true, true, false, false, false};
        case CrystalSystem::Monoclinic: return {true, true, true, false, true, false};
        default: return {true, true, true, true, true, true};
    }
}

// ---------------------------------------------------------------- Cell reduction
Hkl applyT(const Hkl& x, const std::array<int, 9>& T) {
    // x' = x * T (T's columns = coordinates of the new basis vectors in the original basis)
    return Hkl{x.h * T[0] + x.k * T[3] + x.l * T[6], x.h * T[1] + x.k * T[4] + x.l * T[7],
               x.h * T[2] + x.k * T[5] + x.l * T[8]};
}

Reduction reduceCell(const Cell& c) {
    Reduction r;
    std::array<double, 6> g = directMetric(c);
    std::array<int, 9> T{1, 0, 0, 0, 1, 0, 0, 0, 1};

    // --- Buerger: repeatedly take the 3 shortest linearly independent (det=+-1) lattice vectors ---
    for (int iter = 0; iter < 60; ++iter) {
        struct Cand { double len; int h, k, l; };
        std::vector<Cand> cs;
        cs.reserve(125);
        const double cur = g[0] + g[1] + g[2];
        for (int h = -2; h <= 2; ++h)
            for (int k = -2; k <= 2; ++k)
                for (int l = -2; l <= 2; ++l) {
                    if (h == 0 && k == 0 && l == 0) continue;
                    const double len = quadForm(g, {(double)h, (double)k, (double)l});
                    // The summed lengths of three lattice vectors must be below cur to form a shorter basis; drop any single vector exceeding cur
                    if (!(len < cur)) continue;
                    cs.push_back({len, h, k, l});
                }
        std::sort(cs.begin(), cs.end(), [](const Cand& a, const Cand& b) { return a.len < b.len; });
        bool found = false;
        // cs is sorted by length ascending: once a triple sum exceeds cur, larger indices cannot
        // satisfy it, so all three loops can break early via the prefix-sum upper bound (the old m
        // loop used continue, scanning all 1.2e5 necessarily invalid triples, the root cause of reduceCell being slow).
        for (size_t i = 0; i + 2 < cs.size() && !found; ++i) {
            if (cs[i].len + cs[i + 1].len + cs[i + 2].len > cur - 1e-9) break;
            for (size_t j = i + 1; j + 1 < cs.size() && !found; ++j) {
                if (cs[i].len + cs[j].len + cs[j + 1].len > cur - 1e-9) break;
                for (size_t m = j + 1; m < cs.size(); ++m) {
                    if (cs[i].len + cs[j].len + cs[m].len > cur - 1e-9) break;
                    const int d = cs[i].h * (cs[j].k * cs[m].l - cs[j].l * cs[m].k) -
                                  cs[i].k * (cs[j].h * cs[m].l - cs[j].l * cs[m].h) +
                                  cs[i].l * (cs[j].h * cs[m].k - cs[j].k * cs[m].h);
                    if (std::abs(d) != 1) continue;
                    // M's columns = new basis vectors
                    std::array<int, 9> M{cs[i].h, cs[i].k, cs[i].l, cs[j].h, cs[j].k, cs[j].l,
                                         cs[m].h, cs[m].k, cs[m].l};
                    // g' = M g M^T
                    std::array<double, 6> ng{};
                    std::array<std::array<double, 3>, 3> v{
                        {{{(double)cs[i].h, (double)cs[i].k, (double)cs[i].l}},
                         {{(double)cs[j].h, (double)cs[j].k, (double)cs[j].l}},
                         {{(double)cs[m].h, (double)cs[m].k, (double)cs[m].l}}}};
                    auto bip = [&](const std::array<double, 3>& p, const std::array<double, 3>& q) {
                        return g[0] * p[0] * q[0] + g[1] * p[1] * q[1] + g[2] * p[2] * q[2] +
                               g[3] * (p[0] * q[1] + p[1] * q[0]) +
                               g[4] * (p[0] * q[2] + p[2] * q[0]) +
                               g[5] * (p[1] * q[2] + p[2] * q[1]);
                    };
                    ng[0] = bip(v[0], v[0]);
                    ng[1] = bip(v[1], v[1]);
                    ng[2] = bip(v[2], v[2]);
                    ng[3] = bip(v[0], v[1]);
                    ng[4] = bip(v[0], v[2]);
                    ng[5] = bip(v[1], v[2]);
                    // T_new = T * M
                    std::array<int, 9> NT{};
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b) {
                            int s = 0;
                            for (int t = 0; t < 3; ++t) s += T[3 * a + t] * M[3 * t + b];
                            NT[3 * a + b] = s;
                        }
                    T = NT;
                    g = ng;
                    found = true;
                    break;
                }
            }
        }
        if (!found) break;
    }

    // --- Sort by length a<=b<=c (permutation matrix det = +-1) ---
    {
        std::array<int, 3> idx{0, 1, 2};
        // Selection sort (ascending length, ties by index)
        for (int i = 0; i < 3; ++i) {
            int best = i;
            for (int j = i + 1; j < 3; ++j) {
                const double lj = g[idx[j]], lb = g[idx[best]];
                if (lj < lb - 1e-12 || (std::fabs(lj - lb) <= 1e-12 && idx[j] < idx[best])) best = j;
            }
            std::swap(idx[i], idx[best]);
        }
        // New basis order = idx: P's columns = e_{idx[i]}
        std::array<int, 9> P{};
        for (int i = 0; i < 3; ++i) P[3 * idx[i] + i] = 1;
        std::array<int, 9> NT{};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                int s = 0;
                for (int t = 0; t < 3; ++t) s += T[3 * a + t] * P[3 * t + b];
                NT[3 * a + b] = s;
            }
        T = NT;
        auto off = [&](int i, int j) -> double {
            if (i == j) return g[i];
            const int mn = std::min(i, j), mx = std::max(i, j);
            if (mn == 0 && mx == 1) return g[3];
            if (mn == 0 && mx == 2) return g[4];
            return g[5];
        };
        std::array<double, 6> ng{g[idx[0]], g[idx[1]], g[idx[2]], 0, 0, 0};
        ng[3] = off(idx[0], idx[1]);
        ng[4] = off(idx[0], idx[2]);
        ng[5] = off(idx[1], idx[2]);
        g = ng;
    }

    // --- Choose signs: prefer satisfying Niggli (|D|<=B/2, |E|<=A/2, |F|<=A/2),
    //     with angles as obtuse as possible (type I); otherwise take the least-deviating combination ---
    {
        const double A = g[0], B = g[1], C = g[2];
        double bestScore = 1e300;
        std::array<int, 3> bestS{1, 1, 1};
        for (int s0 = -1; s0 <= 1; s0 += 2)
            for (int s1 = -1; s1 <= 1; s1 += 2)
                for (int s2 = -1; s2 <= 1; s2 += 2) {
                    const double D = s1 * s2 * g[5]; // b.c cos(alpha)
                    const double E = s0 * s2 * g[4]; // a.c cos(beta)
                    const double F = s0 * s1 * g[3]; // a.b cos(gamma)
                    double pen = 0.0;
                    pen += std::max(0.0, std::fabs(D) - B / 2.0);
                    pen += std::max(0.0, std::fabs(E) - A / 2.0);
                    pen += std::max(0.0, std::fabs(F) - A / 2.0);
                    // Type consistency: all <=0 (obtuse) preferred
                    const bool allObtuse = (D <= 1e-12 && E <= 1e-12 && F <= 1e-12);
                    const bool allAcute = (D >= -1e-12 && E >= -1e-12 && F >= -1e-12);
                    if (!allObtuse && !allAcute) pen += 1.0;
                    if (allAcute && !allObtuse) pen += 0.05; // prefer type I
                    if (pen < bestScore - 1e-15) {
                        bestScore = pen;
                        bestS = {s0, s1, s2};
                    }
                }
        std::array<int, 9> S{};
        S[0] = bestS[0]; S[4] = bestS[1]; S[8] = bestS[2];
        std::array<int, 9> NT{};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                int v = 0;
                for (int t = 0; t < 3; ++t) v += T[3 * a + t] * S[3 * t + b];
                NT[3 * a + b] = v;
            }
        T = NT;
        g[3] = bestS[0] * bestS[1] * g[3];
        g[4] = bestS[0] * bestS[2] * g[4];
        g[5] = bestS[1] * bestS[2] * g[5];
    }

    Cell out;
    if (!cellFromDirectMetric(g, out)) return r;
    r.cell = out;
    r.T = T;
    r.ok = true;
    return r;
}

// ---------------------------------------------------------------- Criteria
bool niggli2D(double a, double b, double ga) {
    if (a > b + 1e-9) return false;
    return std::fabs(b * std::cos(ga * RAD)) <= a / 2.0 + 1e-9;
}

bool niggli3D(const Cell& c) {
    const double ca = std::cos(c.al * RAD), cb = std::cos(c.be * RAD), cg = std::cos(c.ga * RAD);
    const double A = c.a * c.a, B = c.b * c.b;
    // Type 1 (all <= 91°)
    const bool t1 = (c.al < 91.0 && c.be < 91.0 && c.ga < 91.0) &&
                    (c.c * ca <= c.b / 2.0 + 1e-9) && (c.c * cb <= c.a / 2.0 + 1e-9);
    // Type 2 (all >= 89°)
    const bool t2 = (c.al > 89.0 && c.be > 89.0 && c.ga > 89.0) &&
                    (c.c * std::fabs(ca) <= c.b / 2.0 + 1e-9) &&
                    (c.c * std::fabs(cb) <= c.a / 2.0 + 1e-9) &&
                    (c.b * c.c * std::fabs(ca) + c.a * c.c * std::fabs(cb) +
                         c.a * c.b * std::fabs(cg) <= (A + B) / 2.0 + 1e-9);
    // Mixed angles: accept if either branch is satisfied (the original MATLAB discarded these outright, see review B16)
    if (t1 || t2) return true;
    const bool mixed = !((c.al < 91.0 && c.be < 91.0 && c.ga < 91.0) ||
                         (c.al > 89.0 && c.be > 89.0 && c.ga > 89.0));
    if (!mixed) return false;
    const bool m1 = (c.c * ca <= c.b / 2.0 + 1e-9) && (c.c * cb <= c.a / 2.0 + 1e-9);
    const bool m2 = (c.c * std::fabs(ca) <= c.b / 2.0 + 1e-9) &&
                    (c.c * std::fabs(cb) <= c.a / 2.0 + 1e-9) &&
                    (c.b * c.c * std::fabs(ca) + c.a * c.c * std::fabs(cb) +
                         c.a * c.b * std::fabs(cg) <= (A + B) / 2.0 + 1e-9);
    return m1 || m2;
}

int latticeIndex(const std::vector<Hkl>& hkl) {
    const int n = (int)hkl.size();
    if (n < 3) return 0;
    long long g = 0;
    const int lim = std::min(n, 48);
    for (int i = 0; i < lim; ++i)
        for (int j = i + 1; j < lim; ++j)
            for (int k = j + 1; k < lim; ++k) {
                const long long a[3] = {hkl[i].h, hkl[i].k, hkl[i].l};
                const long long b[3] = {hkl[j].h, hkl[j].k, hkl[j].l};
                const long long c[3] = {hkl[k].h, hkl[k].k, hkl[k].l};
                const long long d = a[0] * (b[1] * c[2] - b[2] * c[1]) -
                                    a[1] * (b[0] * c[2] - b[2] * c[0]) +
                                    a[2] * (b[0] * c[1] - b[1] * c[0]);
                if (d == 0) continue;
                g = (g == 0) ? (d < 0 ? -d : d) : std::gcd(g, d < 0 ? -d : d);
                if (g == 1) return 1;
            }
    return (int)std::min<long long>(g, INT_MAX);
}

// ------------------------------------------------ Constrained mode: supercell / axis settings / dedup
bool acceptLatticeIndex(CrystalSystem sys, const std::vector<Hkl>& hkl) {
    const int li = latticeIndex(hkl);
    if (li <= 1) return true; // 1 = primitive cell; 0 = rank deficient (consistent with the original criterion, not additionally rejected)
    // Trigonal R centering: the observed points generate a sublattice of index 3 and all indices
    // satisfy R-centering extinction in the hexagonal setting. Then this cell is the R rhombohedral
    // cell (a 3-fold cell in hexagonal axes) and should be accepted. Other latticeIndex>1 (including P-type trigonal supercells) are still discarded as supercells, without blanket acceptance.
    if (sys != CrystalSystem::Trigonal || li != 3) return false;
    bool obverse = true, reverse = true;
    for (const Hkl& x : hkl) {
        if (x.h == 0 && x.k == 0 && x.l == 0) continue;
        if ((-x.h + x.k + x.l) % 3 != 0) obverse = false;
        if ((x.h - x.k + x.l) % 3 != 0) reverse = false;
        if (!obverse && !reverse) return false;
    }
    return obverse || reverse;
}

bool transformCell(const Cell& in, const std::array<int, 9>& T, Cell& out) {
    // A' = A*T (A's columns = a,b,c; T's columns = coordinates of the new basis vectors in the old basis)
    const Mat3 A = directMatrix(in);
    Mat3 R{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int t = 0; t < 3; ++t) s += A[3 * i + t] * (double)T[3 * t + j];
            R[3 * i + j] = s;
        }
    auto dot = [&](int i, int j) {
        return R[3 * i + 0] * R[3 * j + 0] + R[3 * i + 1] * R[3 * j + 1] +
               R[3 * i + 2] * R[3 * j + 2];
    };
    const std::array<double, 6> g{dot(0, 0), dot(1, 1), dot(2, 2), dot(0, 1), dot(0, 2), dot(1, 2)};
    return cellFromDirectMetric(g, out);
}

std::vector<std::array<int, 9>> axisSettings(CrystalSystem sys, double ga) {
    std::vector<std::array<int, 9>> out;
    auto perm = [](int i0, int i1, int i2) {
        std::array<int, 9> T{};
        T[0 * 3 + i0] = 1;
        T[1 * 3 + i1] = 1;
        T[2 * 3 + i2] = 1;
        return T;
    };
    auto mul = [](const std::array<int, 9>& A, const std::array<int, 9>& B) {
        std::array<int, 9> R{};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                int s = 0;
                for (int t = 0; t < 3; ++t) s += A[3 * i + t] * B[3 * t + j];
                R[3 * i + j] = s;
            }
        return R;
    };
    switch (sys) {
        case CrystalSystem::Cubic:
            out.push_back(perm(0, 1, 2)); // cubic has no unique axis
            break;
        case CrystalSystem::Tetragonal:
            out.push_back(perm(1, 2, 0)); // unique axis = old a
            out.push_back(perm(0, 2, 1)); // unique axis = old b
            out.push_back(perm(0, 1, 2)); // unique axis = old c
            break;
        case CrystalSystem::Orthorhombic: {
            const int P[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
            for (const auto& p : P) out.push_back(perm(p[0], p[1], p[2]));
            break;
        }
        case CrystalSystem::Monoclinic:
            // beta = (a,c) is free => the unique axis must go to b; the other two choices project the out-of-plane angle to 90°
            out.push_back(perm(1, 0, 2)); // unique axis = old a
            out.push_back(perm(0, 1, 2)); // unique axis = old b
            out.push_back(perm(0, 2, 1)); // unique axis = old c
            break;
        case CrystalSystem::Hexagonal:
        case CrystalSystem::Trigonal: {
            static const int CP[3][3] = {{1, 2, 0}, {0, 2, 1}, {0, 1, 2}};
            // gamma≈60° and gamma=120° are two settings of the same hexagonal lattice: b -> b-a.
            std::array<int, 9> M{};
            M[0] = 1; M[1] = -1; M[4] = 1; M[8] = 1; // columns = a, b-a, c
            const bool alt = (ga < 90.0);
            for (const auto& p : CP) {
                const std::array<int, 9> T = perm(p[0], p[1], p[2]);
                out.push_back(T);
                if (alt) out.push_back(mul(T, M));
            }
            break;
        }
        default:
            out.push_back(perm(0, 1, 2));
            break;
    }
    return out;
}

double axisSettingCost(CrystalSystem sys, const Cell& raw, const std::array<int, 9>& T) {
    Cell cv;
    if (!transformCell(raw, T, cv)) return 1e30;
    const Cell pre = cv; // after transform, before projection
    applySystemConstraints(sys, cv);
    auto rel = [](double x, double y) {
        const double d = std::max(std::fabs(x), std::fabs(y));
        return (d > 1e-9) ? std::fabs(x - y) / d : 0.0;
    };
    double d = rel(cv.a, pre.a) + rel(cv.b, pre.b) + rel(cv.c, pre.c);
    d += (std::fabs(cv.al - pre.al) + std::fabs(cv.be - pre.be) + std::fabs(cv.ga - pre.ga)) / 90.0;
    return d;
}

void monoclinicBranch(Cell& c, Hkl& uvw, std::vector<Hkl>& hkl) {
    if (c.be <= 90.0) return;
    const std::array<int, 9> T{1, 0, 0, 0, 1, 0, 0, 0, -1}; // new c = -old c
    Cell cv;
    if (!transformCell(c, T, cv)) return;
    c = cv;
    uvw = applyT(uvw, T);
    for (Hkl& x : hkl) x = applyT(x, T);
}

bool monoclinicSnapBeta90(Cell& c, const std::vector<Hkl>& hkl) {
    for (const Hkl& x : hkl)
        if (x.l != 0) return false; // with l!=0, beta is determined by the data; no folding
    if (std::fabs(c.be - 90.0) < 1e-6) return true;
    const double s = std::sin(c.be * RAD);
    c.a *= s; // keep a·sinβ (in-plane h reflections depend on it)
    c.c *= s; // keep c·sinβ
    c.be = 90.0;
    return true;
}

uint64_t constrainedCellKey(const Cell& c, const Hkl& uvw) {
    auto mix = [](uint64_t& h, long long x) {
        h ^= (uint64_t)(x + 0x9e3779b97f4a7c15LL);
        h *= 1099511628211ull;
    };
    auto q = [](double x, double t) -> long long { return (long long)std::llround(x / t); };
    // Canonical axis ordering: ascending edge lengths (keep original order on ties), with plane
    // indices transformed by the same permutation. This gives the same key for axis-permuted /
    // axis-setting descriptions of one lattice; lengths and angles are taken as sorted sets, independent of axis labels.
    double L[3] = {c.a, c.b, c.c};
    int idx[3] = {0, 1, 2};
    for (int i = 1; i < 3; ++i) {
        int j = i;
        while (j > 0 && L[idx[j - 1]] > L[idx[j]]) {
            std::swap(idx[j - 1], idx[j]);
            --j;
        }
    }
    uint64_t h = 1469598103934665603ull;
    for (int i = 0; i < 3; ++i) mix(h, q(L[idx[i]], 0.05));
    double A[3] = {c.al, c.be, c.ga};
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (A[j] < A[i]) std::swap(A[i], A[j]);
    for (int i = 0; i < 3; ++i) mix(h, q(A[i], 0.3));
    mix(h, q(cellVolume(c), 0.5));
    int p[3] = {(idx[0] == 0 ? uvw.h : idx[0] == 1 ? uvw.k : uvw.l),
                (idx[1] == 0 ? uvw.h : idx[1] == 1 ? uvw.k : uvw.l),
                (idx[2] == 0 ? uvw.h : idx[2] == 1 ? uvw.k : uvw.l)};
    // Plane normals n and -n are equivalent: make the first nonzero component positive
    if (p[0] == 0 && p[1] == 0) {
        if (p[2] < 0) p[2] = -p[2];
    } else if (p[0] == 0) {
        if (p[1] < 0) { p[1] = -p[1]; p[2] = -p[2]; }
    } else if (p[0] < 0) {
        p[0] = -p[0]; p[1] = -p[1]; p[2] = -p[2];
    }
    mix(h, (long long)p[0] * 10000 + (long long)p[1] * 100 + p[2]);
    return h;
}


std::string hklLabel(int h, int k, int l) {
    auto one = [](int v) {
        if (v >= 0) return std::to_string(v);
        return "-" + std::to_string(-v);
    };
    return one(h) + one(k) + one(l);
}

double dSpacing(const Cell& c, int h, int k, int l) {
    const std::array<double, 6> g = recipMetric(c);
    const double q2 = g[0] * h * h + g[1] * k * k + g[2] * l * l +
                      2.0 * (g[3] * h * k + g[4] * h * l + g[5] * k * l);
    if (!(q2 > 0.0)) return 0.0;
    return TWO_PI / std::sqrt(q2);
}

} // namespace gid
