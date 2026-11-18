// ============================================================================
//  dataio.cpp — see dataio.h
//
//  Excel (.xlsx / .xlsm) reading is done with the xlsxio library
//  (vcpkg, x64-windows-static). xlsxio resolves shared strings / inline
//  strings and all cell types for us, so the loader is robust against the
//  OOXML quirks that a hand-rolled XML walker used to trip over (e.g. <row>
//  elements being nested under <sheetData>, not directly under <worksheet>).
//  CSV loading and CSV export are unchanged.
// ============================================================================
#include "../PathIO.hpp"
#include "dataio.h"

#include <xlsxio_read.h>

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ----------------------------------------------------------------------------
//  helpers
// ----------------------------------------------------------------------------
// Parse a double in invariant locale ("C")
static bool parseDouble(const std::string& s, double& v) {
    if (s.empty()) return false;
    const char* p = s.c_str();
    char* end = nullptr;
    double d = std::strtod(p, &end);
    if (end == p) return false;
    while (*end && std::isspace((unsigned char)*end)) ++end;
    if (*end) return false;
    v = d;
    return true;
}

static std::string lowerExt(const std::string& path) {
    size_t dot = path.find_last_of('.');
    if (dot == std::string::npos) return std::string();
    std::string e = path.substr(dot + 1);
    for (char& c : e) c = (char)std::tolower((unsigned char)c);
    return e;
}

void sortRowsByFirstColumn(Mat& M) {
    if (M.r <= 1) return;
    std::vector<int> idx(M.r);
    for (int i = 0; i < M.r; ++i) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) {
        return M(i, 0) < M(j, 0);
    });
    Mat R(M.r, M.c);
    for (int i = 0; i < M.r; ++i)
        for (int j = 0; j < M.c; ++j) R(i, j) = M(idx[i], j);
    M = R;
}

// ----------------------------------------------------------------------------
//  CSV loader
// ----------------------------------------------------------------------------
static bool loadCSV(const std::string& path, Mat& out, std::string& errorMessage) {
    auto f = pathio::input(path, std::ios::binary);
    if (!f) { errorMessage = "Cannot open file."; return false; }

    std::vector<std::vector<double>> rows;
    std::string line;
    bool first = true;
    while (std::getline(f, line)) {
        // strip UTF-8 BOM on the first line
        if (first) {
            first = false;
            if (line.size() >= 3 &&
                (unsigned char)line[0] == 0xEF && (unsigned char)line[1] == 0xBB &&
                (unsigned char)line[2] == 0xBF)
                line.erase(0, 3);
        }
        // strip \r
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();
        // skip empty lines
        bool onlySpace = true;
        for (char ch : line) if (!std::isspace((unsigned char)ch)) { onlySpace = false; break; }
        if (onlySpace) continue;

        // separator detection: ';' > ',' > tab > whitespace
        char sep = ' ';
        if (line.find(';') != std::string::npos) sep = ';';
        else if (line.find(',') != std::string::npos) sep = ',';
        else if (line.find('\t') != std::string::npos) sep = '\t';

        std::vector<double> row;
        std::string cell;
        auto flushCell = [&]() {
            double v = 0.0;
            if (!parseDouble(cell, v)) v = 0.0;   // non-numeric -> 0
            row.push_back(v);
            cell.clear();
        };
        for (char ch : line) {
            if ((sep == ' ' && std::isspace((unsigned char)ch)) || (sep != ' ' && ch == sep))
                flushCell();
            else
                cell.push_back(ch);
        }
        if (!cell.empty()) flushCell();
        if (!row.empty()) rows.push_back(std::move(row));
    }

    if (rows.empty()) { errorMessage = "No numeric rows found."; return false; }

    int cols = 0;
    for (auto& r : rows) cols = std::max(cols, (int)r.size());
    out = Mat((int)rows.size(), cols, 0.0);
    for (int i = 0; i < (int)rows.size(); ++i)
        for (int j = 0; j < (int)rows[i].size(); ++j)
            out(i, j) = rows[i][j];
    return true;
}

// ----------------------------------------------------------------------------
//  XLSX loader (xlsxio)
// ----------------------------------------------------------------------------
struct XlsxCtx {
    std::vector<std::vector<double>> rows;
    int maxRow = -1;   // 0-based
    int maxCol = -1;   // 0-based
};

static int xlsxCellCb(size_t row, size_t col, const char* value, void* data) {
    XlsxCtx* ctx = static_cast<XlsxCtx*>(data);
    int r = (row >= 1) ? (int)row - 1 : 0;
    int c = (col >= 1) ? (int)col - 1 : 0;
    if (r > ctx->maxRow) ctx->maxRow = r;
    if (c > ctx->maxCol) ctx->maxCol = c;
    while ((int)ctx->rows.size() <= r)
        ctx->rows.emplace_back();
    auto& rv = ctx->rows[r];
    if ((int)rv.size() <= c) rv.resize((size_t)c + 1, 0.0);
    double v = 0.0;
    if (value && value[0]) parseDouble(std::string(value), v);   // non-numeric -> 0
    rv[c] = v;
    return 0;   // continue
}

static int xlsxRowCb(size_t, size_t, void*) { return 0; }

static bool loadXLSX(const std::string& path, Mat& out, std::string& errorMessage) {
    xlsxioreader reader = xlsxioread_open(path.c_str());
    if (!reader) { errorMessage = "Cannot open xlsx file (not a valid .xlsx?)."; return false; }

    XlsxCtx ctx;
    int rc = xlsxioread_process(reader, nullptr, XLSXIOREAD_SKIP_NONE,
                                xlsxCellCb, xlsxRowCb, &ctx);
    xlsxioread_close(reader);
    if (rc != 0) { errorMessage = "Failed to parse xlsx worksheet."; return false; }

    int nRows = ctx.maxRow + 1;
    int nCols = ctx.maxCol + 1;
    if (nRows <= 0 || nCols <= 0) { errorMessage = "Worksheet is empty."; return false; }

    // pad rows / columns so the matrix is dense (missing cells stay 0.0)
    ctx.rows.resize((size_t)nRows);
    for (auto& rv : ctx.rows) rv.resize((size_t)nCols, 0.0);

    out = Mat(nRows, nCols, 0.0);
    for (int i = 0; i < nRows; ++i)
        for (int j = 0; j < nCols; ++j)
            out(i, j) = ctx.rows[(size_t)i][(size_t)j];
    return true;
}

bool loadDataTable(const std::string& path, Mat& out, std::string& errorMessage) {
    std::string ext = lowerExt(path);
    if (ext == "xlsx" || ext == "xlsm")
        return loadXLSX(path, out, errorMessage);
    return loadCSV(path, out, errorMessage);   // .csv and everything else
}

// ----------------------------------------------------------------------------
//  CSV export of one solution (MATLAB header layout)
// ----------------------------------------------------------------------------
static std::string fmtG(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.10g", v);
    return buf;
}
static std::string fmt4(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.4f", v);
    return buf;
}
static std::string fmt0(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.0f", v);
    return buf;
}

bool exportSolutionCSV(const std::string& path,
                       const double* par14,
                       const Mat& indicesPage,
                       const Mat& rawData,
                       const double* phipsi,
                       std::string& errorMessage) {
    auto f = pathio::output(path, std::ios::binary);
    if (!f) { errorMessage = "Cannot create file."; return false; }

    const std::string x = "xxxx";
    auto line = [&](const std::vector<std::string>& c) {
        for (size_t i = 0; i < c.size(); ++i) {
            if (i) f << ",";
            f << c[i];
        }
        f << "\r\n";
    };

    line({ "a [A]", "b [A]", "c [A]",
           phipsi ? "PHI [deg]" : x, phipsi ? "PSI [deg]" : x });
    line({ fmt4(par14[3]), fmt4(par14[4]), fmt4(par14[5]),
           phipsi ? fmtG(phipsi[0]) : x, phipsi ? fmtG(phipsi[1]) : x });
    line({ "alpha [deg]", "beta [deg]", "gamma [deg]", x, x });
    line({ fmt4(par14[6]), fmt4(par14[7]), fmt4(par14[8]), x, x });
    line({ "u", "v", "w", x, x });
    line({ fmt0(par14[0]), fmt0(par14[1]), fmt0(par14[2]), x, x });
    line({ "vol [A^3]", fmtG(par14[9]), x, x, x });
    line({ "dqxyz [1/A]", "dqz [1/A]", "dqxy [1/A]", x, x });
    line({ fmt4(par14[10]), fmt4(par14[11]), fmt4(par14[12]), x, x });
    line({ x, x, x, x, x });
    line({ "qxy [A]", "qz [A]", "h", "k", "l" });

    for (int i = 0; i < rawData.r; ++i) {
        std::vector<std::string> c = { fmtG(rawData(i, 0)), fmtG(rawData(i, 1)) };
        if (i < indicesPage.r) {
            c.push_back(fmtG(indicesPage(i, 0)));
            c.push_back(fmtG(indicesPage(i, 1)));
            c.push_back(fmtG(indicesPage(i, 2)));
        } else {
            c.push_back(""); c.push_back(""); c.push_back("");
        }
        line(c);
    }
    return true;
}
