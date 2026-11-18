// ============================================================================
//  PeakPlot.h — peak scatter plot (q_xy, q_z) with predicted positions overlaid
// ============================================================================
#pragma once

#include <wx/panel.h>

#include <vector>

#include "core/types.h"

class PeakPlot : public wxPanel {
public:
    explicit PeakPlot(wxWindow* parent);

    void setPeaks(const std::vector<gid::Peak>& p);
    void setPredicted(const std::vector<gid::Peak>& p);

private:
    void onPaint(wxPaintEvent&);

    std::vector<gid::Peak> m_peaks;
    std::vector<gid::Peak> m_pred;
};
