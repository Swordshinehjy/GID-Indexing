// ============================================================================
//  dataio.h — file input/output for measurement tables
//  Load: .csv and .xlsx (first worksheet; numeric cells, missing cells = 0),
//  mirroring MATLAB readmatrix(...,'OutputType','double') closely enough for
//  the GIDInd workflows. Export: solution sheet as CSV.
//  Single responsibility: file parsing / writing. No UI code.
// ============================================================================
#pragma once

#include "linalg.h"
#include <string>

// Load a numeric table from .csv / .xlsx. Rows/columns as stored.
// Returns false and sets errorMessage on failure.
bool loadDataTable(const std::string& path, Mat& out, std::string& errorMessage);

// Stable sort of the rows by column 1 ascending (the app sorts every load).
void sortRowsByFirstColumn(Mat& M);

// Export one indexed solution as CSV, replicating the layout of the MATLAB
// ExportsolutionButtonPushed (11 header rows + [qxy qz h k l] per datapoint).
//   par14      : final_parameter_stored row (u v w a b c alpha beta gamma Vol
//                RMSD_qxyz RMSD_qxy RMSD_qz dqspec)
//   indicesPage: n x 3 (h k l), already rounded
//   rawData    : app.t (n x 2)
//   phipsi     : nullptr or [phi psi] in degrees (ppout row)
bool exportSolutionCSV(const std::string& path,
                       const double* par14,
                       const Mat& indicesPage,
                       const Mat& rawData,
                       const double* phipsi,
                       std::string& errorMessage);
