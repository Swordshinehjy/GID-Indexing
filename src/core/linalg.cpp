// ============================================================================
//  linalg.cpp — implementation of numeric primitives
// ============================================================================
#include "linalg.h"

double vecnorm(const std::vector<double>& v) {
    double s = 0; for (double x : v) s += x * x; return std::sqrt(s);
}

std::vector<double> cross3(const std::vector<double>& u, const std::vector<double>& w) {
    return { u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0] };
}

double dot3(const std::vector<double>& u, const std::vector<double>& w) {
    return u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
}

double matDet(const Mat& A) {
    assert(A.r == A.c); // traverses with A.r as the square order; a non-square matrix silently gives wrong results
    int n = A.r;
    std::vector<double> m = A.a;
    double det = 1.0;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i) if (std::fabs(m[(size_t)i * n + k]) > std::fabs(m[(size_t)p * n + k])) p = i;
        if (std::fabs(m[(size_t)p * n + k]) < 1e-300) return 0.0;
        if (p != k) {
            for (int j = 0; j < n; ++j) std::swap(m[(size_t)p * n + j], m[(size_t)k * n + j]);
            det = -det;
        }
        det *= m[(size_t)k * n + k];
        for (int i = k + 1; i < n; ++i) {
            double f = m[(size_t)i * n + k] / m[(size_t)k * n + k];
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) m[(size_t)i * n + j] -= f * m[(size_t)k * n + j];
        }
    }
    return det;
}

bool matSolve(const Mat& A, const Mat& B, Mat& X) {
    assert(A.r == A.c && B.r == A.r); // A square and B with the same row count as A, otherwise striding by n overruns
    int n = A.r, m = B.c;
    std::vector<double> a = A.a, b = B.a;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i) if (std::fabs(a[(size_t)i * n + k]) > std::fabs(a[(size_t)p * n + k])) p = i;
        if (std::fabs(a[(size_t)p * n + k]) < 1e-300) return false;
        if (p != k) {
            for (int j = 0; j < n; ++j) std::swap(a[(size_t)p * n + j], a[(size_t)k * n + j]);
            for (int j = 0; j < m; ++j) std::swap(b[(size_t)p * m + j], b[(size_t)k * m + j]);
        }
        for (int i = k + 1; i < n; ++i) {
            double f = a[(size_t)i * n + k] / a[(size_t)k * n + k];
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) a[(size_t)i * n + j] -= f * a[(size_t)k * n + j];
            for (int j = 0; j < m; ++j) b[(size_t)i * m + j] -= f * b[(size_t)k * m + j];
        }
    }
    X = Mat(n, m, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j < m; ++j) {
            double s = b[(size_t)i * m + j];
            for (int k2 = i + 1; k2 < n; ++k2) s -= a[(size_t)i * n + k2] * X(k2, j);
            X(i, j) = s / a[(size_t)i * n + i];
        }
    }
    return true;
}

bool matInv3(const Mat& A, Mat& Ai) {
    assert(A.r == 3 && A.c == 3); // addressed directly by indices 0..8 below
    double d = matDet(A);
    if (std::fabs(d) < 1e-300) return false;
    Ai = Mat(3, 3);
    const double* m = A.a.data();
    Ai(0, 0) = (m[4] * m[8] - m[5] * m[7]) / d;
    Ai(0, 1) = (m[2] * m[7] - m[1] * m[8]) / d;
    Ai(0, 2) = (m[1] * m[5] - m[2] * m[4]) / d;
    Ai(1, 0) = (m[5] * m[6] - m[3] * m[8]) / d;
    Ai(1, 1) = (m[0] * m[8] - m[2] * m[6]) / d;
    Ai(1, 2) = (m[2] * m[3] - m[0] * m[5]) / d;
    Ai(2, 0) = (m[3] * m[7] - m[4] * m[6]) / d;
    Ai(2, 1) = (m[1] * m[6] - m[0] * m[7]) / d;
    Ai(2, 2) = (m[0] * m[4] - m[1] * m[3]) / d;
    return true;
}

bool lstsq(const Mat& A, const Mat& b, Mat& x) {
    Mat AtA = A.T() * A;
    Mat Atb = A.T() * b;
    return matSolve(AtA, Atb, x);
}

bool lstsqQR(const Mat& A, const Mat& b, Mat& x) {
    int m = A.r, n = A.c;
    if (m < n || m == 0) return false;
    Mat R(A), B(b);
    std::vector<double> v(m);
    for (int k = 0; k < n; ++k) {
        double nrm2 = 0.0;
        for (int i = k; i < m; ++i) nrm2 += R(i, k) * R(i, k);
        double nrm = std::sqrt(nrm2);
        if (nrm < 1e-300) return false;
        double alpha = (R(k, k) >= 0.0) ? -nrm : nrm;
        for (int i = k; i < m; ++i) v[i] = R(i, k);
        v[k] -= alpha;
        double vtv = 0.0;
        for (int i = k; i < m; ++i) vtv += v[i] * v[i];
        if (vtv < 1e-300) return false;
        for (int j = k; j < n; ++j) {
            double dot = 0.0;
            for (int i = k; i < m; ++i) dot += v[i] * R(i, j);
            double f = 2.0 * dot / vtv;
            for (int i = k; i < m; ++i) R(i, j) -= f * v[i];
        }
        for (int j = 0; j < B.c; ++j) {
            double dot = 0.0;
            for (int i = k; i < m; ++i) dot += v[i] * B(i, j);
            double f = 2.0 * dot / vtv;
            for (int i = k; i < m; ++i) B(i, j) -= f * v[i];
        }
    }
    x = Mat(n, B.c, 0.0);
    for (int j = 0; j < B.c; ++j)
        for (int i = n - 1; i >= 0; --i) {
            double s = B(i, j);
            for (int k2 = i + 1; k2 < n; ++k2) s -= R(i, k2) * x(k2, j);
            x(i, j) = s / R(i, i);
        }
    return true;
}

void sortrows(Mat& M, const std::vector<int>& keys) {
    std::vector<int> idx(M.r);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) {
        for (int k : keys) {
            double a = M(i, k), b = M(j, k);
            if (a < b) return true;
            if (a > b) return false;
        }
        return false;
    });
    Mat R(M.r, M.c);
    for (int i = 0; i < M.r; ++i)
        for (int j = 0; j < M.c; ++j) R(i, j) = M(idx[i], j);
    M = R;
}

void sortrows(Mat& M) {
    std::vector<int> keys(M.c);
    std::iota(keys.begin(), keys.end(), 0);
    sortrows(M, keys);
}

void uniqueRowsSorted(const Mat& A, Mat& C, std::vector<int>& ia) {
    std::vector<int> idx(A.r);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) {
        for (int k = 0; k < A.c; ++k) {
            double a = A(i, k), b = A(j, k);
            if (a < b) return true;
            if (a > b) return false;
        }
        return i < j;
    });
    ia.clear();
    Mat tmp(A.r, A.c);
    int n = 0;
    for (int t = 0; t < A.r; ++t) {
        int i = idx[t];
        if (n == 0) { ia.push_back(i); for (int j = 0; j < A.c; ++j) tmp(n, j) = A(i, j); ++n; continue; }
        bool same = true;
        for (int j = 0; j < A.c; ++j) if (tmp(n - 1, j) != A(i, j)) { same = false; break; }
        if (!same) { ia.push_back(i); for (int j = 0; j < A.c; ++j) tmp(n, j) = A(i, j); ++n; }
    }
    C = Mat(n, A.c);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < A.c; ++j) C(i, j) = tmp(i, j);
}

void uniquetolVec(const std::vector<double>& A, double tol, std::vector<int>& ia) {
    double dscale = 0.0;
    for (double v : A) dscale = std::max(dscale, std::fabs(v));
    std::vector<int> idx(A.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int i, int j) {
        if (A[i] != A[j]) return A[i] < A[j];
        return i < j;
    });
    ia.clear();
    if (A.empty()) return;
    double rep = A[idx[0]];
    ia.push_back(idx[0]);
    for (size_t t = 1; t < idx.size(); ++t) {
        if (std::fabs(A[idx[t]] - rep) > tol * dscale) {
            rep = A[idx[t]];
            ia.push_back(idx[t]);
        }
    }
    // NOTE: ia is kept in ascending-VALUE order (C = A(ia) is sorted), matching
    // MATLAB uniquetol. Do NOT re-sort by index.
}

std::vector<int> minkVec(const std::vector<double>& v, int k) {
    std::vector<int> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) {
        if (v[i] != v[j]) return v[i] < v[j];
        return i < j;
    });
    if ((int)idx.size() > k) idx.resize(k);
    return idx;
}

Mat nchoosek3(int n) {
    long long total = (long long)n * (n - 1) * (n - 2) / 6;
    Mat R((int)total, 3);
    long long t = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = j + 1; k < n; ++k) {
                R((int)t, 0) = i + 1; R((int)t, 1) = j + 1; R((int)t, 2) = k + 1; ++t;
            }
    return R;
}

Mat LPermutation(const std::vector<double>& hkvalues, int n) {
    // MATLAB original (LPermutation.m):
    //   column i is filled with  m,..,m | m-1,..,m-1 | ... | 1,..,1
    //   where each value is repeated  m^(i-1)  times, and the whole pattern
    //   is repeated  m^(n-i)  times.  Reproducing this exact enumeration
    //   keeps the row order identical to the reference implementation.
    int nn = (int)hkvalues.size();
    long long rows = 1;
    for (int i = 0; i < n; ++i) rows *= nn;
    Mat L((int)rows, n);
    for (int i = 1; i <= n; ++i) {
        long long repblk = 1;                     // m^(i-1)
        for (int t = 0; t < i - 1; ++t) repblk *= nn;
        long long nblk = 1;                       // m^(n-i)
        for (int t = 0; t < n - i; ++t) nblk *= nn;
        long long r0 = 0;
        for (long long b = 0; b < nblk; ++b)
            for (int val = nn; val >= 1; --val)
                for (long long c = 0; c < repblk; ++c, ++r0)
                    L((int)r0, i - 1) = hkvalues[val - 1];
    }
    return L;
}

Mat PERMVEC_NEW(int limit_fuer_hk) {
    int nh = 2 * limit_fuer_hk + 1;
    Mat HK(nh * nh, 2);
    int t = 0;
    for (int j = 0; j < nh; ++j)
        for (int i = 0; i < nh; ++i) {
            HK(t, 0) = -limit_fuer_hk + j;
            HK(t, 1) = -limit_fuer_hk + i;
            ++t;
        }
    return HK;
}

void removeRowsWithAnyZero(Mat& M) {
    int w = 0;
    for (int i = 0; i < M.r; ++i) {
        bool hasZero = false;
        for (int j = 0; j < M.c; ++j) if (M(i, j) == 0.0) { hasZero = true; break; }
        if (!hasZero) {
            if (w != i) for (int j = 0; j < M.c; ++j) M(w, j) = M(i, j);
            ++w;
        }
    }
    M.r = w;
}

void removeRowsAllZero(Mat& M) {
    int w = 0;
    for (int i = 0; i < M.r; ++i) {
        bool allZero = true;
        for (int j = 0; j < M.c; ++j) if (M(i, j) != 0.0) { allZero = false; break; }
        if (!allZero) {
            if (w != i) for (int j = 0; j < M.c; ++j) M(w, j) = M(i, j);
            ++w;
        }
    }
    M.r = w;
}

// MATLAB round(x,n): round to n decimals, half away from zero
double roundTo(double x, int n) {
    double f = std::pow(10.0, n);
    return std::round(x * f) / f;
}

void removeRowsAnyNonFinite(Mat& M) {
    int w = 0;
    for (int i = 0; i < M.r; ++i) {
        bool bad = false;
        for (int j = 0; j < M.c; ++j) if (!std::isfinite(M(i, j))) { bad = true; break; }
        if (!bad) {
            if (w != i) for (int j = 0; j < M.c; ++j) M(w, j) = M(i, j);
            ++w;
        }
    }
    M.r = w;
}
