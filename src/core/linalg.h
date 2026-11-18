// ============================================================================
//  linalg.h — small dense linear algebra + MATLAB-like array helpers
//  Single responsibility: numeric primitives shared by all math kernels.
// ============================================================================
#pragma once

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double DEG = 180.0 / M_PI;

// ----------------------------------------------------------------------------
//  Small dense-matrix helper (row-major)
// ----------------------------------------------------------------------------
struct Mat {
    int r = 0, c = 0;
    std::vector<double> a;

    Mat() {}
    Mat(int R, int C, double v = 0.0) : r(R), c(C), a((size_t)R * C, v) {}
    Mat(int R, int C, const std::vector<double>& d) : r(R), c(C), a(d) { a.resize((size_t)R * C, 0.0); }

    double& operator()(int i, int j) { return a[(size_t)i * c + j]; }
    double  operator()(int i, int j) const { return a[(size_t)i * c + j]; }
    double& at(int i, int j) { return a[(size_t)i * c + j]; }
    double  at(int i, int j) const { return a[(size_t)i * c + j]; }

    Mat rowSlice(int i) const {
        Mat R(1, c);
        for (int j = 0; j < c; ++j) R(0, j) = (*this)(i, j);
        return R;
    }

    static Mat eye(int n) {
        Mat I(n, n, 0.0);
        for (int i = 0; i < n; ++i) I(i, i) = 1.0;
        return I;
    }
    Mat operator*(const Mat& B) const {
        // On a dimension mismatch the inner loop reads B(k, j) striding by this->c and goes out of bounds without throwing -- this must be guarded in numeric kernels
        assert(c == B.r);
        Mat R(r, B.c, 0.0);
        for (int i = 0; i < r; ++i)
            for (int k = 0; k < c; ++k) {
                double v = (*this)(i, k);
                if (v == 0.0) continue;
                for (int j = 0; j < B.c; ++j) R(i, j) += v * B(k, j);
            }
        return R;
    }
    Mat operator*(double s) const {
        Mat R(r, c);
        for (size_t i = 0; i < a.size(); ++i) R.a[i] = a[i] * s;
        return R;
    }
    Mat T() const {
        Mat R(c, r, 0.0);
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j) R(j, i) = (*this)(i, j);
        return R;
    }
};

// 3-vector helpers
double vecnorm(const std::vector<double>& v);
std::vector<double> cross3(const std::vector<double>& u, const std::vector<double>& w);
double dot3(const std::vector<double>& u, const std::vector<double>& w);

// det of small square matrix (n<=8) by LU with partial pivoting
double matDet(const Mat& A);

// Solve A X = B (A: n x n, B: n x m) with partial pivoting. Returns false if singular.
bool matSolve(const Mat& A, const Mat& B, Mat& X);

// 3x3 inverse
bool matInv3(const Mat& A, Mat& Ai);

// Least squares solve A x = b via normal equations (A: m x n, m>=n)
bool lstsq(const Mat& A, const Mat& b, Mat& x);

// Least squares solve A x = b via Householder QR (A: m x n, m>=n).
// Numerically preferable for the overdetermined 4x3 systems of part 2.
bool lstsqQR(const Mat& A, const Mat& b, Mat& x);

// ----------------------------------------------------------------------------
//  MATLAB-like helpers
// ----------------------------------------------------------------------------

// stable lexicographic sort of rows by selected key columns (default: all)
void sortrows(Mat& M, const std::vector<int>& keys);
void sortrows(Mat& M);

// unique rows (sorted lexicographically, keep first) -> C, ia (indices into A)
void uniqueRowsSorted(const Mat& A, Mat& C, std::vector<int>& ia);

// MATLAB uniquetol(A,tol) for a vector column: default DataScale = max(abs(A)).
// Returns ia = indices of first occurrence of each unique value (values ascending).
void uniquetolVec(const std::vector<double>& A, double tol, std::vector<int>& ia);

// MATLAB round(x,n): round to n decimals, half away from zero
double roundTo(double x, int n);

// mink(v,k): indices of the k smallest values, sorted ascending by value
// (ties broken by original index).
std::vector<int> minkVec(const std::vector<double>& v, int k);

// nchoosek(1:n,3) in lexicographic order
Mat nchoosek3(int n);

// LPermutation: all n-tuples (with repetition) of hkvalues → m^n rows (m = |hkvalues|).
// Enumeration as implemented (linalg.cpp): column i (1-based) repeats each value m^(i-1)
// times, so the FIRST column cycles fastest and the LAST column stays constant the
// longest; within one block values run in DESCENDING order (hkvalues[m-1] → hkvalues[0]).
// The set of tuples is order-independent, so only row order (and hence tie-breaking in
// callers) depends on this convention.
Mat LPermutation(const std::vector<double>& hkvalues, int n);

// function_PERMVEC_NEW: all (h,k) pairs, h slowest / k fastest, -lim..lim ascending
Mat PERMVEC_NEW(int limit_fuer_hk);

// Row-filtering utilities (MATLAB logical indexing)
void removeRowsWithAnyZero(Mat& M);
void removeRowsAllZero(Mat& M);
void removeRowsAnyNonFinite(Mat& M);
