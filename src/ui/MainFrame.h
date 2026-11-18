// ============================================================================
//  MainFrame.h — GIDIndexing main window (toolbar + menu bar)
//
//  Left: peak plot + peak table (edit q_xy / q_z directly; hkl and residuals written back)
//  Right: parameter page (switches with mode) + solution list; compute runs in background
//  via gid::Progress and can be stopped at any time
// ============================================================================
#pragma once

#include <wx/frame.h>
#include <wx/string.h>

#include "core/engine.h"
#include "core/types.h"

#include <atomic>
#include <thread>
#include <vector>

class wxBoxSizer;
class wxCheckBox;
class wxMenu;
class wxChoice;
class wxGauge;
class wxGrid;
class wxGridEvent;
class wxNotebook;
class wxPanel;
class wxScrolledWindow;
class wxSpinCtrl;
class wxSplitterWindow;
class wxStaticText;
class wxTextCtrl;
class wxToolBar;

class PeakPlot;

class MainFrame : public wxFrame {
public:
    MainFrame();
    ~MainFrame() override;

private:
    // ---------------------------------------------------------------- Build
    void buildMenu();
    void buildToolbar();
    void buildLayout();
    wxScrolledWindow* buildSpecPage(wxWindow* parent);
    wxScrolledWindow* buildNoSpecPage(wxWindow* parent);
    wxScrolledWindow* buildNoSpec2Page(wxWindow* parent);
    wxScrolledWindow* buildFixedPage(wxWindow* parent);
    wxScrolledWindow* buildSimplePage(wxWindow* parent);
    void bindEvents();

    // ---------------------------------------------------------------- Data
    void importFile(const wxString& path);
    void refreshPeakGrid();
    bool collectPeaks(); // peak table -> m_peaks
    void appendPeakRow();
    void deleteSelectedPeakRows();
    void clearPeaks();
    void autoDetectQSpec();

    // ---------------------------------------------------------------- Compute
    bool collectParams(gid::EngineParams& ep);
    void startIndexing();
    void stopIndexing();
    void onFinished(const gid::EngineResult& res);
    void setRunning(bool on);

    // ---------------------------------------------------------------- Results
    void refreshSolutionGrid();
    void showSolution(int row);
    void exportSolution();
    void syncModeUI();

    // ---------------------------------------------------------------- Events
    void onClose(wxCloseEvent& e);
    void onImport(wxCommandEvent&);
    void onExport(wxCommandEvent&);
    void onRun(wxCommandEvent&);
    void onStop(wxCommandEvent&);
    void onClear(wxCommandEvent&);
    void onModeChoice(wxCommandEvent&);
    void onModeMenu(wxCommandEvent& e);
    void onAppend(wxCommandEvent&);
    void onDelete(wxCommandEvent&);
    void onAutoQSpec(wxCommandEvent&);
    void onSolSelect(wxGridEvent& e);
    void onUsage(wxCommandEvent&);
    void onAbout(wxCommandEvent&);

private:
    // Widgets
    wxToolBar* m_toolbar = nullptr;
    wxMenu* m_menuRun = nullptr;
    wxChoice* m_modeChoice = nullptr;
    wxGauge* m_gauge = nullptr;
    wxSplitterWindow* m_splitter = nullptr;
    PeakPlot* m_plot = nullptr;
    wxGrid* m_peakGrid = nullptr;
    wxNotebook* m_note = nullptr;
    wxPanel* m_paramHost = nullptr;
    wxScrolledWindow* m_pageSpec = nullptr;
    wxScrolledWindow* m_pageNoSpec = nullptr;
    wxScrolledWindow* m_pageNoSpec2 = nullptr;
    wxScrolledWindow* m_pageFixed = nullptr;
    wxScrolledWindow* m_pageSimple = nullptr;
    wxPanel* m_solPage = nullptr;
    wxGrid* m_solGrid = nullptr;
    wxStaticText* m_solInfo = nullptr;

    // Parameters (with q_spec)
    wxTextCtrl* m_txtQSpec = nullptr;
    wxChoice* m_chP1Alg = nullptr;
    wxChoice* m_specSystem = nullptr;
    wxSpinCtrl* m_spRedundant = nullptr, *m_spUv = nullptr, *m_spHk1 = nullptr;
    wxSpinCtrl* m_spHk2 = nullptr, *m_spL = nullptr, *m_spW = nullptr;
    wxSpinCtrl* m_spLines = nullptr, *m_spMaxPart1 = nullptr;
    wxSpinCtrl* m_spRefineTop = nullptr, *m_spMaxSol = nullptr;
    wxTextCtrl* m_txtDqxy = nullptr, *m_txtDqz = nullptr;
    wxTextCtrl* m_txtDqxyz = nullptr, *m_txtDqspec = nullptr;
    wxCheckBox* m_chkNiggli = nullptr, *m_chkRestrict = nullptr;
    wxTextCtrl* m_txtALo = nullptr, *m_txtAHi = nullptr;
    wxTextCtrl* m_txtBLo = nullptr, *m_txtBHi = nullptr;
    wxTextCtrl* m_txtCLo = nullptr, *m_txtCHi = nullptr;
    wxTextCtrl* m_txtVolLo = nullptr, *m_txtVolHi = nullptr;
    wxSpinCtrl* m_spThreads = nullptr;

    // Parameters (no q_spec)
    wxChoice* m_nsMethod = nullptr;
    wxChoice* m_nsSystem = nullptr;
    wxSpinCtrl *m_nsSeedH = nullptr, *m_nsSeedK = nullptr, *m_nsSeedL = nullptr;
    wxSpinCtrl *m_nsSeedPeaks = nullptr, *m_nsMinRows = nullptr;
    wxSpinCtrl *m_nsGridH = nullptr, *m_nsGridK = nullptr, *m_nsGridL = nullptr;
    wxSpinCtrl *m_nsBeam = nullptr, *m_nsRefineTop = nullptr, *m_nsMaxSol = nullptr;
    wxTextCtrl *m_nsDqxy = nullptr, *m_nsDqz = nullptr, *m_nsDqxyz = nullptr;
    wxCheckBox* m_nsRestrict = nullptr;
    wxTextCtrl *m_nsALo = nullptr, *m_nsAHi = nullptr;
    wxTextCtrl *m_nsBLo = nullptr, *m_nsBHi = nullptr;
    wxTextCtrl *m_nsCLo = nullptr, *m_nsCHi = nullptr;
    wxTextCtrl *m_nsVolLo = nullptr, *m_nsVolHi = nullptr;
    wxSpinCtrl* m_nsThreads = nullptr;

    // Parameters (no q_spec, two-step): all search ranges, criteria and reduction switches are adjustable
    wxTextCtrl *m_ns2QsLo = nullptr, *m_ns2QsHi = nullptr;
    wxSpinCtrl *m_ns2Steps = nullptr, *m_ns2Refine = nullptr, *m_ns2MaxTau = nullptr;
    wxTextCtrl *m_ns2RodTol = nullptr, *m_ns2RodPrior = nullptr;
    wxSpinCtrl* m_ns2RodOrders = nullptr;
    wxSpinCtrl *m_ns2CoarseHk1 = nullptr, *m_ns2CoarseHk2 = nullptr;
    wxSpinCtrl *m_ns2CoarseSeedPeaks = nullptr, *m_ns2CoarseRedundant = nullptr;
    wxSpinCtrl* m_ns2CoarseMaxPart1 = nullptr;

    wxChoice* m_ns2System = nullptr;
    wxSpinCtrl *m_ns2UvLimit = nullptr, *m_ns2Hk1 = nullptr, *m_ns2Hk2 = nullptr;
    wxSpinCtrl *m_ns2Seed3Half = nullptr, *m_ns2SeedPeaks = nullptr;
    wxSpinCtrl *m_ns2Redundant = nullptr, *m_ns2MaxPart1 = nullptr;
    wxTextCtrl* m_ns2Dqxy = nullptr;

    wxSpinCtrl *m_ns2SeedL = nullptr, *m_ns2SeedLFallback = nullptr;
    wxSpinCtrl *m_ns2LLimit = nullptr, *m_ns2WLimit = nullptr;
    wxSpinCtrl *m_ns2MaxTrials = nullptr, *m_ns2RefineTop = nullptr;
    wxTextCtrl* m_ns2Dqz = nullptr;

    wxCheckBox *m_ns2Restrict = nullptr, *m_ns2Niggli = nullptr;
    wxTextCtrl *m_ns2ALo = nullptr, *m_ns2AHi = nullptr;
    wxTextCtrl *m_ns2BLo = nullptr, *m_ns2BHi = nullptr;
    wxTextCtrl *m_ns2CLo = nullptr, *m_ns2CHi = nullptr;
    wxTextCtrl *m_ns2AlLo = nullptr, *m_ns2AlHi = nullptr;
    wxTextCtrl *m_ns2BeLo = nullptr, *m_ns2BeHi = nullptr;
    wxTextCtrl *m_ns2GaLo = nullptr, *m_ns2GaHi = nullptr;
    wxTextCtrl *m_ns2VolLo = nullptr, *m_ns2VolHi = nullptr;

    wxSpinCtrl *m_ns2IpCands = nullptr, *m_ns2MaxSol = nullptr, *m_ns2Threads = nullptr;
    wxTextCtrl* m_ns2Dqxyz = nullptr;

    // Parameters (fixed contact plane)
    wxSpinCtrl* m_spU = nullptr, *m_spV = nullptr, *m_spWplane = nullptr;
    wxSpinCtrl* m_spMaxH = nullptr, *m_spMaxK = nullptr, *m_spMaxL = nullptr;
    wxChoice* m_chSystem = nullptr;
    wxTextCtrl* m_txtSa = nullptr, *m_txtSb = nullptr, *m_txtSc = nullptr;
    wxTextCtrl* m_txtSal = nullptr, *m_txtSbe = nullptr, *m_txtSga = nullptr;
    wxSpinCtrl* m_spSeed = nullptr, *m_spTrials = nullptr, *m_spTopN = nullptr;
    wxSpinCtrl* m_spRefIter = nullptr;
    wxTextCtrl* m_txtAccept = nullptr;
    wxCheckBox* m_chkAutoSign = nullptr;

    // Parameters (simple indexing, out-of-plane fixed to 001)
    wxTextCtrl *m_simpleTol = nullptr, *m_simpleTolXy = nullptr, *m_simpleTolZ = nullptr;
    wxTextCtrl* m_simpleRod = nullptr;
    wxSpinCtrl *m_simpleHk = nullptr, *m_simpleL = nullptr;
    wxSpinCtrl *m_simpleTopK = nullptr, *m_simpleGrid = nullptr;

    // State
    std::vector<gid::Peak> m_peaks;
    std::vector<gid::Solution> m_solutions;
    int m_solRow = -1;
    gid::IndexMode m_mode = gid::IndexMode::WithSpec;
    // Crystal system and Δq_xyz tolerance of the last "with q_spec" run (for the mismatch warning)
    gid::CrystalSystem m_lastSystem = gid::CrystalSystem::Triclinic;
    double m_lastDqxyz = -1.0; // <0 = no warning

    std::atomic<bool> m_cancel{false};
    std::atomic<long long> m_lastReport{0};
    bool m_running = false;
    std::thread m_worker;
};
