// ============================================================================
//  PeakPlot.cpp — see PeakPlot.h
// ============================================================================
#include "PeakPlot.h"

#include <wx/dcbuffer.h>
#include <wx/settings.h>

#include <algorithm>
#include <cmath>

PeakPlot::PeakPlot(wxWindow* parent)
    : wxPanel(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize,
              wxFULL_REPAINT_ON_RESIZE | wxBORDER_THEME) {
    SetBackgroundStyle(wxBG_STYLE_PAINT);
    SetBackgroundColour(*wxWHITE);
    Bind(wxEVT_PAINT, &PeakPlot::onPaint, this);
}

void PeakPlot::setPeaks(const std::vector<gid::Peak>& p) {
    m_peaks = p;
    Refresh();
}

void PeakPlot::setPredicted(const std::vector<gid::Peak>& p) {
    m_pred = p;
    Refresh();
}

void PeakPlot::onPaint(wxPaintEvent&) {
    wxAutoBufferedPaintDC dc(this);
    const wxSize sz = GetClientSize();
    dc.SetBackground(wxBrush(wxColour(253, 253, 251)));
    dc.Clear();

    const int mL = 58, mR = 16, mT = 16, mB = 38;
    if (sz.x < mL + mR + 40 || sz.y < mT + mB + 40) return;

    dc.SetTextForeground(wxColour(70, 70, 70));
    dc.SetFont(wxFont(8, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL));

    if (m_peaks.empty() && m_pred.empty()) {
        dc.DrawText(L"Import peaks to show (q_xy, q_z) distribution", mL + 12, mT + 8);
        return;
    }

    double xmax = 0.05, zmin = 0.0, zmax = 0.05;
    for (const gid::Peak& p : m_peaks) {
        xmax = std::max(xmax, p.qxy);
        zmin = std::min(zmin, p.qz);
        zmax = std::max(zmax, p.qz);
    }
    for (const gid::Peak& p : m_pred) {
        xmax = std::max(xmax, p.qxy);
        zmin = std::min(zmin, p.qz);
        zmax = std::max(zmax, p.qz);
    }
    xmax *= 1.06;
    double zspan = std::max(zmax - zmin, 1e-3) * 1.10;

    const int w = sz.x - mL - mR;
    const int h = sz.y - mT - mB;
    auto px = [&](double q) { return mL + (int)std::lround(q / xmax * (double)w); };
    auto py = [&](double q) {
        return mT + h - (int)std::lround((q - zmin) / zspan * (double)h);
    };

    // Grid
    dc.SetPen(wxPen(wxColour(232, 232, 232), 1));
    for (int i = 0; i <= 5; ++i) {
        const int x = mL + w * i / 5;
        dc.DrawLine(x, mT, x, mT + h);
        const double q = xmax * i / 5.0;
        dc.DrawText(wxString::Format("%.2f", q), x - 12, mT + h + 6);
    }
    for (int i = 0; i <= 4; ++i) {
        const int y = mT + h - h * i / 4;
        dc.DrawLine(mL, y, mL + w, y);
        const double q = zmin + zspan * i / 4.0;
        dc.DrawText(wxString::Format("%.2f", q), 8, y - 7);
    }

    // Axes
    dc.SetPen(wxPen(wxColour(110, 110, 110), 1));
    dc.DrawLine(mL, mT, mL, mT + h);
    dc.DrawLine(mL, mT + h, mL + w, mT + h);
    dc.DrawText(L"q_xy", mL + w - 26, mT + h + 20);
    dc.DrawRotatedText(L"q_z", 8, mT + 44, 90.0);

    // Predicted positions
    if (!m_pred.empty()) {
        dc.SetPen(wxPen(wxColour(200, 40, 40), 1));
        for (const gid::Peak& p : m_pred) {
            const int x = px(p.qxy), y = py(p.qz);
            dc.DrawLine(x - 4, y - 4, x + 4, y + 4);
            dc.DrawLine(x - 4, y + 4, x + 4, y - 4);
        }
    }

    // Observed peaks
    dc.SetPen(wxPen(wxColour(30, 80, 160), 1));
    dc.SetBrush(wxBrush(wxColour(60, 130, 220)));
    for (const gid::Peak& p : m_peaks) dc.DrawCircle(px(p.qxy), py(p.qz), 3);
}
