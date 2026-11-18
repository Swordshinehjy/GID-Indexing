// ============================================================================
//  cell.h — crystallography basics: unit cell <-> reciprocal lattice <-> forward model / crystal-system constraints / cell reduction
//  All values are double; angles are in degrees externally and radians internally.
//  Conventions (consistent with GIDInd / the reference implementation):
//      Direct basis vectors as columns A = [a b c] (rows = x,y,z)
//      Reciprocal B = 2*pi * A^{-1}; row i of B is a*/b*/c*, satisfying a_j . b_i = 2*pi*delta
//      g(hkl) = (h,k,l) * B        (row vector times matrix)
//      n_hat  = G_spec/|G_spec|,  G_spec = (u,v,w) * B
//      q_z = g . n_hat,   q_xy = sqrt(|g|^2 - q_z^2)
// ============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gid {

constexpr double TWO_PI = 6.28318530717958647692;
constexpr double RAD = M_PI / 180.0;
constexpr double DEG = 180.0 / M_PI;

enum class CrystalSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Hexagonal,
    Trigonal,
    Cubic
};

std::vector<std::string> crystalSystemNames();
std::string crystalSystemName(CrystalSystem s);

struct Cell {
    double a = 10.0, b = 10.0, c = 10.0;
    double al = 90.0, be = 90.0, ga = 90.0; // alpha=(b,c) beta=(a,c) gamma=(a,b)
};

struct Hkl {
    int h = 0, k = 0, l = 0;
};

using Mat3 = std::array<double, 9>; // row-major 3x3

// ---------------------------------------------------------------- Basis vectors and metric
// Direct basis vectors (columns are a,b,c)
Mat3 directMatrix(const Cell& c);
// Reciprocal basis vectors (rows are a*,b*,c*)
Mat3 reciprocalMatrix(const Cell& c);

// Reciprocal metric G* = B B^T; element order [G11,G22,G33,G12,G13,G23]
std::array<double, 6> recipMetric(const Cell& c);
// Direct metric G = A^T A; element order [G11,G22,G33,G12,G13,G23]
std::array<double, 6> directMetric(const Cell& c);

// Recover the unit cell from the reciprocal metric: G_dir = (2pi)^2 (G*)^-1
bool cellFromRecipMetric(const std::array<double, 6>& gs, Cell& out);
// Recover the unit cell from the direct metric
bool cellFromDirectMetric(const std::array<double, 6>& gd, Cell& out);

double cellVolume(const Cell& c);
// Related to reciprocal volume / det(G*): V = (2pi)^3 / sqrt(det G*)
double volumeFromRecipDet(double detGs);

// ---------------------------------------------------------------- Forward model
struct QPoint {
    double qxy = 0.0, qz = 0.0;
};

// Given the unit cell and contact plane (u,v,w), compute (q_xy, q_z) for (h,k,l)
QPoint forwardHkl(const Mat3& recip, const Hkl& uvw, int h, int k, int l);
QPoint forwardHkl(const Cell& c, const Hkl& uvw, int h, int k, int l);

// |G_spec| = sqrt((u,v,w) G* (u,v,w)^T)
double qSpecOf(const Mat3& recip, const Hkl& uvw);

// Mutual conversion between the out-of-plane component coefficients nu = (a*.n, b*.n, c*.n) and (u,v,w)
// nu = G* (u,v,w)^T / q_spec ;  (u,v,w)/q_spec = (G*)^-1 nu
std::array<double, 3> nuFromUvw(const Mat3& recip, const Hkl& uvw);
// Find the smallest integer (u,v,w) along the nu direction (returns false if no small integer is found)
bool uvwFromNu(const Mat3& recip, const std::array<double, 3>& nu, Hkl& out, int maxIndex = 6);

// ------------------------------------------- Part ONE / TWO reparameterization (wiki §1.2)
// a* = A + kappa*e, b* = B + lambda*e, c* = rho*e (e = unit vector along c*)
// z_a = 2pi/(a sin gamma), z_b = 2pi/(b sin gamma), C = z_a z_b cos gamma
struct Frame {
    double za = 0.0, zb = 0.0, cg = 0.0; // cg = z_a z_b cos(gamma)
    double kap = 0.0, lam = 0.0, rho = 0.0;
};
// (a,b,gamma) -> (z_a, z_b, C)
void abgToZ(double a, double b, double ga, double& za, double& zb, double& C);
// (z_a, z_b, C) -> (a,b,gamma)
void zToAbg(double za, double zb, double C, double& a, double& b, double& ga);
// Full unit cell -> Frame
Frame frameOf(const Cell& c);
// (z_a,z_b,C,kappa,lambda,rho) -> full unit cell (wiki §3.4 closed form)
bool cellFromFrame(double za, double zb, double C, double kap, double lam, double rho, Cell& out);

// ---------------------------------------------------------------- Crystal-system constraints
// Snap the parameters constrained by the crystal system to their constrained values (e.g. Cubic: b=c=a, angles 90)
void applySystemConstraints(CrystalSystem sys, Cell& c);
// Free-parameter mask for the crystal system (false = not fitted)
std::array<bool, 6> systemFreeMask(CrystalSystem sys);

// ---------------------------------------------------------------- Cell reduction
// Supercell criterion in constrained mode: accept when the generated index of the observed
// reciprocal lattice is 1; trigonal (R-centered) allows index 3, but requires all indices to
// satisfy R-centering extinction under hexagonal axes (obverse: -h+k+l ≡ 0 mod 3, reverse: h-k+l ≡ 0 mod 3). Other latticeIndex > 1 are still supercells and are rejected.
bool acceptLatticeIndex(CrystalSystem sys, const std::vector<Hkl>& hkl);

// Valid axis settings for the crystal system (basis transform T: columns = integer coordinates
// of the new basis vectors in the old basis, det ±1). Tetragonal/hexagonal/trigonal: 3 choices
// of unique axis; orthorhombic: 6; monoclinic: 3; cubic: 1. ga is the original cell's gamma: for hexagonal/trigonal, gamma<90° additionally yields the b->b-a 60°<->120° equivalent setting.
std::vector<std::array<int, 9>> axisSettings(CrystalSystem sys, double ga);

// Metric cost of projecting an axis setting onto the target crystal system (relative length change
// + angle change / 90°). Used to rank axis settings before projection: the lowest-cost setting is
// usually the correct one, so a single assignment suffices instead of assigning all 6 settings.
double axisSettingCost(CrystalSystem sys, const Cell& raw, const std::array<int, 9>& T);

// Basis transform A' = A*T (T as above). hkl / uvw are transformed in sync with applyT; volume and physics unchanged.
bool transformCell(const Cell& in, const std::array<int, 9>& T, Cell& out);

// Full dedup key in constrained mode = cell fingerprint after canonical axis ordering (sorted
// lengths + sorted angles + volume) with the plane indices under that permutation; replaces the raw cell+uvw key when reduceCell is absent.
uint64_t constrainedCellKey(const Cell& c, const Hkl& uvw);
// Monoclinic beta canonical branch: for beta > 90 switch to the equivalent c -> -c setting (uvw and hkl in sync), a purely descriptive transform.
void monoclinicBranch(Cell& c, Hkl& uvw, std::vector<Hkl>& hkl);
// Monoclinic l=0 degeneracy: when all l==0 and the specular direction w==0, c and beta are
// inseparable (they enter the model only as a·sinβ, c·sinβ); snap to beta=90 keeping a·sinβ,
// c·sinβ, folding the beta family into one solution. Returns whether applied; skipped when l!=0 reflections exist.
bool monoclinicSnapBeta90(Cell& c, const std::vector<Hkl>& hkl);
// Buerger + Niggli reduction. T's columns = integer coordinates of the new basis vectors
// in the original basis. Also transforms uvw and all hkl: x' = x * T
struct Reduction {
    Cell cell;                      // reduced unit cell
    std::array<int, 9> T{1, 0, 0, 0, 1, 0, 0, 0, 1}; // columns = new basis vectors
    bool ok = false;
};
Reduction reduceCell(const Cell& c);
Hkl applyT(const Hkl& x, const std::array<int, 9>& T);

// 2D Niggli criterion (Part ONE): a<=b and b|cos gamma| <= a/2
bool niggli2D(double a, double b, double ga);
// 3D Niggli criterion (used for filtering in Part TWO)
bool niggli3D(const Cell& c);

std::string hklLabel(int h, int k, int l);
double dSpacing(const Cell& c, int h, int k, int l);

// Decide from the assigned (hkl) matrix whether the cell is a supercell:
// return the gcd of all 3x3 minors (= index in Z^3 of the sublattice generated by the observed reciprocal points).
//  1  = the observed points generate exactly this reciprocal lattice (primitive cell)
//  d>1= the observed points generate only a sublattice of index d -> this cell volume is d times the true cell (should be discarded)
//  0  = rank deficient (cannot decide)
int latticeIndex(const std::vector<Hkl>& hkl);

} // namespace gid
