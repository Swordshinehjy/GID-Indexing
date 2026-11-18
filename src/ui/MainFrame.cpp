// ============================================================================
//  MainFrame.cpp — see MainFrame.h
// ============================================================================
#include "MainFrame.h"

#include "PeakPlot.h"

#include "core/dataio.h"
#include "core/linalg.h"

#include <wx/accel.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/choicdlg.h>
#include <wx/choice.h>
#include <wx/dcmemory.h>
#include <wx/display.h>
#include <wx/filedlg.h>
#include <wx/gauge.h>
#include <wx/grid.h>
#include <wx/menu.h>
#include <wx/msgdlg.h>
#include <wx/notebook.h>
#include <wx/panel.h>
#include <wx/scrolwin.h>
#include <wx/sizer.h>
#include <wx/spinctrl.h>
#include <wx/splitter.h>
#include <wx/stattext.h>
#include <wx/statusbr.h>
#include <wx/textctrl.h>
#include <wx/toolbar.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>

namespace {

enum {
    ID_IMPORT = wxID_HIGHEST + 1,
    ID_EXPORT,
    ID_RUN,
    ID_STOP,
    ID_CLEAR,
    ID_APPEND,
    ID_DELETE,
    ID_AUTO_QSPEC,
    ID_MODE_CHOICE,
    ID_MODE_SPEC,
    ID_MODE_NOSPEC,
    ID_MODE_NOSPEC2,
    ID_MODE_FIXED,
    ID_MODE_SIMPLE,
    ID_USAGE,
    ID_ABOUT
};

enum class Icon { Import, Export, Run, Stop, Clear };

// "Method" combo-box items: display name ↔ Part1Algorithm one-to-one. Both widget
// creation and reading use this table; add a new algorithm only here to keep indices aligned.
struct P1AlgItem {
    const wchar_t* name;
    gid::Part1Algorithm alg;
};
const P1AlgItem kP1AlgItems[] = {
    {L"Improved", gid::Part1Algorithm::Improved},
    {L"PeriodVote/DDL", gid::Part1Algorithm::PeriodVoteDdl},
    {L"RANSAC", gid::Part1Algorithm::Ransac},
    {L"Integer LSQ", gid::Part1Algorithm::Ills},
    {L"SCRM", gid::Part1Algorithm::Scrm},
};
constexpr int kP1AlgDefault = 1; // PeriodVote/DDL (see SpecParams::part1Alg)

// Toolbar icons: drawn in code (white background + mask), no external resources
wxBitmap makeIcon(Icon kind) {
    wxBitmap bmp(24, 24);
    wxMemoryDC dc(bmp);
    dc.SetBackground(*wxWHITE_BRUSH);
    dc.Clear();
    const wxColour ink(0x1F, 0x3B, 0x57);
    dc.SetPen(wxPen(ink, 2));
    dc.SetBrush(*wxTRANSPARENT_BRUSH);
    switch (kind) {
        case Icon::Import: // down arrow + tray
            dc.DrawLine(12, 3, 12, 13);
            dc.DrawLine(8, 9, 12, 13);
            dc.DrawLine(16, 9, 12, 13);
            dc.DrawLine(5, 15, 19, 15);
            dc.DrawLine(5, 15, 5, 20);
            dc.DrawLine(19, 15, 19, 20);
            dc.DrawLine(5, 20, 19, 20);
            break;
        case Icon::Export: // up arrow + tray
            dc.DrawLine(12, 20, 12, 10);
            dc.DrawLine(8, 14, 12, 10);
            dc.DrawLine(16, 14, 12, 10);
            dc.DrawLine(5, 4, 19, 4);
            dc.DrawLine(5, 4, 5, 9);
            dc.DrawLine(19, 4, 19, 9);
            dc.DrawLine(5, 9, 19, 9);
            break;
        case Icon::Run: { // filled triangle
            dc.SetBrush(wxBrush(ink));
            wxPoint p[3] = {wxPoint(8, 4), wxPoint(8, 20), wxPoint(20, 12)};
            dc.DrawPolygon(3, p);
            break;
        }
        case Icon::Stop: // filled square
            dc.SetBrush(wxBrush(ink));
            dc.DrawRectangle(7, 7, 10, 10);
            break;
        case Icon::Clear: // cross
            dc.DrawLine(6, 6, 18, 18);
            dc.DrawLine(18, 6, 6, 18);
            break;
    }
    dc.SelectObject(wxNullBitmap);
    bmp.SetMask(new wxMask(bmp, *wxWHITE));
    return bmp;
}

wxStaticText* addLabel(wxWindow* p, wxSizer* s, const wxString& text) {
    auto* t = new wxStaticText(p, wxID_ANY, text);
    s->Add(t, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 6);
    return t;
}

wxSpinCtrl* addSpin(wxWindow* p, wxSizer* s, const wxString& label, int val, int lo, int hi,
                    const wxString& tip = wxString()) {
    addLabel(p, s, label);
    auto* c = new wxSpinCtrl(p, wxID_ANY, wxString::Format("%d", val), wxDefaultPosition,
                             wxSize(120, -1), wxSP_ARROW_KEYS, lo, hi, val);
    if (!tip.empty()) c->SetToolTip(tip);
    s->Add(c, 0, wxEXPAND | wxBOTTOM, 3);
    return c;
}

wxTextCtrl* addText(wxWindow* p, wxSizer* s, const wxString& label, const wxString& val,
                    const wxString& tip = wxString()) {
    addLabel(p, s, label);
    auto* c = new wxTextCtrl(p, wxID_ANY, val, wxDefaultPosition, wxSize(120, -1));
    if (!tip.empty()) c->SetToolTip(tip);
    s->Add(c, 0, wxEXPAND | wxBOTTOM, 3);
    return c;
}

// Two inputs per row (lower / upper)
void addRange(wxWindow* p, wxSizer* s, const wxString& label, double lo, double hi,
              wxTextCtrl*& outLo, wxTextCtrl*& outHi) {
    addLabel(p, s, label);
    auto* hs = new wxBoxSizer(wxHORIZONTAL);
    outLo = new wxTextCtrl(p, wxID_ANY, wxString::Format("%g", lo), wxDefaultPosition, wxSize(62, -1));
    outHi = new wxTextCtrl(p, wxID_ANY, wxString::Format("%g", hi), wxDefaultPosition, wxSize(62, -1));
    hs->Add(outLo, 1, wxEXPAND | wxRIGHT, 4);
    hs->Add(new wxStaticText(p, wxID_ANY, L"~"), 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 4);
    hs->Add(outHi, 1, wxEXPAND);
    s->Add(hs, 0, wxEXPAND | wxBOTTOM, 3);
}

void addCheck(wxWindow* p, wxSizer* s, const wxString& label, bool val, wxCheckBox*& out,
              const wxString& tip = wxString()) {
    out = new wxCheckBox(p, wxID_ANY, label);
    out->SetValue(val);
    if (!tip.empty()) out->SetToolTip(tip);
    // The check box is a control in the right column (aligned with the inputs); the left cell stays empty
    s->AddSpacer(1);
    s->Add(out, 0, wxALIGN_CENTER_VERTICAL | wxBOTTOM, 3);
}

// Important: every row reserves a pair of cells (left label / right control). Section
// headers must also fill a whole row (header + one spacer), otherwise the grid shifts
// by one column and later labels land on the right while their controls land on the left.
void addSection(wxWindow* p, wxSizer* s, const wxString& text) {
    auto* t = new wxStaticText(p, wxID_ANY, text);
    t->SetFont(t->GetFont().Bold());
    s->Add(t, 0, wxALIGN_CENTER_VERTICAL | wxTOP | wxBOTTOM, 4);
    s->AddSpacer(1); // hold the right column to keep the two-cell rhythm
}

void setupGrid(wxGrid* g, const wxString* cols, int n, const int* widths) {
    g->SetRowLabelSize(0);
    g->SetColLabelSize(26);
    g->SetSelectionMode(wxGrid::wxGridSelectRows);
    g->EnableDragRowMove(false);
    g->EnableDragColMove(false);
    g->EnableDragColSize(true);
    g->EnableDragRowSize(false);
    g->SetDefaultCellFont(wxFont(9, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL));
    for (int i = 0; i < n; ++i) {
        g->SetColLabelValue(i, cols[i]);
        g->SetColSize(i, widths[i]);
    }
}

} // namespace

// ============================================================================
MainFrame::MainFrame()
    : wxFrame(nullptr, wxID_ANY, L"GIDIndexing — GIWAXS Peak Indexing", wxDefaultPosition,
              wxSize(1280, 680)) {
    buildMenu();
    buildToolbar();
    buildLayout();
    bindEvents();

    CreateStatusBar(3);
    const int widths[] = {-6, -3, -2};
    SetStatusWidths(3, widths);
    SetStatusText(L"Peaks: 0", 1);

    syncModeUI();
    SetStatusText(L"Ready: import a peak table, then click \"Run Indexing\"", 0);

    wxAcceleratorEntry acc[5];
    acc[0].Set(wxACCEL_CTRL, 'O', ID_IMPORT);
    acc[1].Set(wxACCEL_CTRL, 'S', ID_EXPORT);
    acc[2].Set(wxACCEL_NORMAL, WXK_F5, ID_RUN);
    acc[3].Set(wxACCEL_NORMAL, WXK_ESCAPE, ID_STOP);
    acc[4].Set(wxACCEL_CTRL, 'N', ID_APPEND);
    SetAcceleratorTable(wxAcceleratorTable(5, acc));

    // Small-screen / DPI scaling: fit the work area and keep the title bar on screen
    const int disp = wxDisplay::GetFromWindow(this);
    if (disp != wxNOT_FOUND) {
        const wxRect work = wxDisplay(disp).GetClientArea();
        const int margin = FromDIP(24);
        const int maxW = std::max(400, work.width - margin);
        const int maxH = std::max(300, work.height - margin);
        // Clamp the min size first, then set the size, so it cannot push the window off-screen
        SetMinSize(wxSize(std::min(FromDIP(980), maxW), std::min(FromDIP(620), maxH)));
        wxSize sz = GetSize();
        sz.x = std::min(sz.x, maxW);
        sz.y = std::min(sz.y, maxH);
        SetSize(sz);
        const wxSize actual = GetSize();
        SetPosition(wxPoint(work.x + (work.width - actual.x) / 2,
                            work.y + (work.height - actual.y) / 2));
    } else {
        SetMinSize(wxSize(980, 620));
        Centre();
    }
}

MainFrame::~MainFrame() {
    m_cancel = true;
    if (m_worker.joinable()) m_worker.join();
}

// ---------------------------------------------------------------- Menu bar
void MainFrame::buildMenu() {
    auto* mb = new wxMenuBar();

    auto* file = new wxMenu();
    file->Append(ID_IMPORT, L"Import Peaks...\tCtrl+O", L"CSV / TXT / XLSX, at least two columns (q_xy, q_z)");
    file->AppendSeparator();
    file->Append(ID_EXPORT, L"Export Solution...\tCtrl+S", L"Export cell parameters and (hkl) for all peaks");
    file->AppendSeparator();
    file->Append(wxID_EXIT, L"Exit");
    mb->Append(file, L"File");

    auto* edit = new wxMenu();
    edit->Append(ID_APPEND, L"Add Row\tCtrl+N");
    edit->Append(ID_DELETE, L"Delete Row\tCtrl+D");
    edit->AppendSeparator();
    edit->Append(ID_CLEAR, L"Clear Peaks");
    mb->Append(edit, L"Edit");

    m_menuRun = new wxMenu();
    m_menuRun->Append(ID_RUN, L"Run Indexing\tF5", L"Solve the cell with the current mode and parameters");
    m_menuRun->Append(ID_STOP, L"Stop\tEsc");
    m_menuRun->AppendSeparator();
    m_menuRun->AppendRadioItem(ID_MODE_SPEC, L"q_spec (2-step)");
    m_menuRun->AppendRadioItem(ID_MODE_NOSPEC, L"No q_spec (1-step)");
    m_menuRun->AppendRadioItem(ID_MODE_NOSPEC2, L"No q_spec (2-step)");
    m_menuRun->AppendRadioItem(ID_MODE_FIXED, L"Fixed (uvw)");
    m_menuRun->AppendRadioItem(ID_MODE_SIMPLE, L"Simple (001)");
    m_menuRun->Check(ID_MODE_SPEC, true);
    mb->Append(m_menuRun, L"Run");

    auto* help = new wxMenu();
    help->Append(ID_USAGE, L"Usage", L"Four steps to index; five indexing modes");
    help->Append(ID_ABOUT, L"About");
    mb->Append(help, L"Help");

    SetMenuBar(mb);
}

// ---------------------------------------------------------------- Toolbar
void MainFrame::buildToolbar() {
    m_toolbar = CreateToolBar(wxTB_HORIZONTAL | wxTB_FLAT | wxTB_NODIVIDER);
    m_toolbar->SetToolBitmapSize(wxSize(24, 24));

    m_toolbar->AddTool(ID_IMPORT, L"Import", makeIcon(Icon::Import), L"Import peaks (Ctrl+O)");
    m_toolbar->AddTool(ID_EXPORT, L"Export", makeIcon(Icon::Export), L"Export solution (Ctrl+S)");
    m_toolbar->AddSeparator();

    m_modeChoice = new wxChoice(m_toolbar, ID_MODE_CHOICE, wxDefaultPosition, wxSize(190, -1));
    m_modeChoice->Append(L"q_spec (2-step)");
    m_modeChoice->Append(L"No q_spec (1-step)");
    m_modeChoice->Append(L"No q_spec (2-step)");
    m_modeChoice->Append(L"Fixed (uvw)");
    m_modeChoice->Append(L"Simple (001)");
    m_modeChoice->SetSelection(0);
    m_modeChoice->SetToolTip(L"Indexing mode");
    m_toolbar->AddControl(m_modeChoice);
    m_toolbar->AddSeparator();

    m_toolbar->AddTool(ID_RUN, L"Run", makeIcon(Icon::Run), L"Run indexing (F5)");
    m_toolbar->AddTool(ID_STOP, L"Stop", makeIcon(Icon::Stop), L"Stop computation (Esc)");
    m_toolbar->AddSeparator();
    m_toolbar->AddTool(ID_CLEAR, L"Clear", makeIcon(Icon::Clear), L"Clear peak table");

    m_gauge = new wxGauge(m_toolbar, wxID_ANY, 100, wxDefaultPosition, wxSize(180, 16),
                          wxGA_HORIZONTAL | wxGA_SMOOTH);
    m_toolbar->AddStretchableSpace();
    m_toolbar->AddControl(m_gauge);

    m_toolbar->Realize();
    m_toolbar->EnableTool(ID_STOP, false);
}

// ---------------------------------------------------------------- Layout
void MainFrame::buildLayout() {
    m_splitter = new wxSplitterWindow(this, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                                      wxSP_LIVE_UPDATE | wxSP_3D);

    // ---- Left: plot + peak table
    auto* left = new wxPanel(m_splitter);
    auto* ls = new wxBoxSizer(wxVERTICAL);
    m_plot = new PeakPlot(left);
    m_peakGrid = new wxGrid(left, wxID_ANY, wxDefaultPosition, wxSize(-1, 220));
    m_peakGrid->CreateGrid(0, 8);
    const wxString pcols[] = {L"#", L"q_xy", L"q_z", L"h", L"k", L"l", L"Δq_xy", L"Δq_z"};
    const int pwidths[] = {38, 86, 86, 44, 44, 44, 76, 76};
    setupGrid(m_peakGrid, pcols, 8, pwidths);
    ls->Add(m_plot, 1, wxEXPAND | wxALL, 4);
    ls->Add(m_peakGrid, 0, wxEXPAND | wxLEFT | wxRIGHT | wxBOTTOM, 4);
    left->SetSizer(ls);

    // ---- Right: parameters / solutions
    m_note = new wxNotebook(m_splitter, wxID_ANY);

    m_paramHost = new wxPanel(m_note);
    auto* hs = new wxBoxSizer(wxVERTICAL);
    m_pageSpec = buildSpecPage(m_paramHost);
    m_pageNoSpec = buildNoSpecPage(m_paramHost);
    m_pageNoSpec2 = buildNoSpec2Page(m_paramHost);
    m_pageFixed = buildFixedPage(m_paramHost);
    m_pageSimple = buildSimplePage(m_paramHost);
    hs->Add(m_pageSpec, 1, wxEXPAND);
    hs->Add(m_pageNoSpec, 1, wxEXPAND);
    hs->Add(m_pageNoSpec2, 1, wxEXPAND);
    hs->Add(m_pageFixed, 1, wxEXPAND);
    hs->Add(m_pageSimple, 1, wxEXPAND);
    m_paramHost->SetSizer(hs);
    m_note->AddPage(m_paramHost, L"Params", true);

    m_solPage = new wxPanel(m_note);
    auto* ss = new wxBoxSizer(wxVERTICAL);
    m_solGrid = new wxGrid(m_solPage, wxID_ANY);
    m_solGrid->CreateGrid(0, 13);
    const wxString scols[] = {L"#",       L"(uvw)",  L"a",      L"b",        L"c",
                              L"α",       L"β",      L"γ",      L"V",        L"q_spec",
                              L"R(q_xyz)", L"R(q_xy)", L"R(q_z)"};
    const int swidths[] = {34, 66, 66, 66, 66, 60, 60, 60, 78, 72, 76, 72, 72};
    setupGrid(m_solGrid, scols, 13, swidths);
    m_solGrid->EnableEditing(false);
    m_solInfo = new wxStaticText(m_solPage, wxID_ANY, L"Not computed", wxDefaultPosition, wxDefaultSize,
                                 wxST_ELLIPSIZE_END);
    ss->Add(m_solGrid, 1, wxEXPAND | wxALL, 4);
    ss->Add(m_solInfo, 0, wxEXPAND | wxLEFT | wxRIGHT | wxBOTTOM, 6);
    m_solPage->SetSizer(ss);
    m_note->AddPage(m_solPage, L"Solutions");

    m_splitter->SplitVertically(left, m_note, 500);
    m_splitter->SetMinimumPaneSize(330);
    m_splitter->SetSashGravity(0.0);

    auto* top = new wxBoxSizer(wxVERTICAL);
    top->Add(m_splitter, 1, wxEXPAND);
    SetSizer(top);

    // The sash position needs the real client width, so set it after layout takes effect
    CallAfter([this]() {
        const int w = m_splitter->GetClientSize().x;
        if (w > 500) m_splitter->SetSashPosition(std::max(360, (int)(w * 0.52)));
    });
}

wxScrolledWindow* MainFrame::buildSpecPage(wxWindow* parent) {
    auto* sw = new wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                                    wxVSCROLL | wxTAB_TRAVERSAL);
    auto* s = new wxFlexGridSizer(0, 2, 4, 10);
    s->AddGrowableCol(1);

    addSection(sw, s, L"Specular Peaks");
    addLabel(sw, s, L"q_spec");
    {
        auto* hs = new wxBoxSizer(wxHORIZONTAL);
        m_txtQSpec = new wxTextCtrl(sw, wxID_ANY, L"", wxDefaultPosition, wxSize(90, -1));
        auto* b = new wxButton(sw, ID_AUTO_QSPEC, L"Auto", wxDefaultPosition, wxSize(78, -1));
        hs->Add(m_txtQSpec, 1, wxEXPAND | wxRIGHT, 6);
        hs->Add(b, 0);
        s->Add(hs, 0, wxEXPAND | wxBOTTOM, 3);
    }

    addSection(sw, s, L"Algorithm");
    addLabel(sw, s, L"Method");
    {
        m_chP1Alg = new wxChoice(sw, wxID_ANY, wxDefaultPosition, wxSize(210, -1));
        for (const P1AlgItem& it : kP1AlgItems) m_chP1Alg->Append(it.name);
        m_chP1Alg->SetSelection(kP1AlgDefault);
        m_chP1Alg->SetToolTip(
            L"Implementations for the q_spec two-step mode:\n"
            L"Improved: l triplet enumeration (slow, fallback when the whole round fails)\n"
            L"PeriodVote/DDL: difference-vector period voting (default)\n"
            L"RANSAC: 2-peak minimal sample + specular rows, deterministic with per-peak early exit\n"
            L"Integer LSQ: difference-vector lattice clustering + interval scan for exact t\n"
            L"SCRM (new): step 1 enumerates in-plane (u,v) only, without w; step 2 resolves\n"
            L"(κ,λ,ρ,w,l) from the specular constraint; standalone, shares no code with Improved");
        s->Add(m_chP1Alg, 0, wxEXPAND | wxBOTTOM, 3);
    }

    addSection(sw, s, L"System");
    addLabel(sw, s, L"System");
    m_specSystem = new wxChoice(sw, wxID_ANY, wxDefaultPosition, wxSize(120, -1));
    for (const std::string& n : gid::crystalSystemNames())
        m_specSystem->Append(wxString::FromUTF8(n.c_str()));
    m_specSystem->SetSelection(0); // triclinic = unconstrained
    m_specSystem->SetToolTip(
        L"Metric constraints for the target crystal system (higher symmetry = fewer free parameters and solutions)\n"
        L"Triclinic: unconstrained (default, same behavior as older versions)\n"
        L"Other systems: before projection, each legal axis setting is compared (3 unique axes for tetra-\n"
        L"gonal/hexagonal/trigonal, 6 permutations for orthorhombic, 3 for monoclinic); the cheapest wins.\n"
        L"LM refines only its free parameters and skips the 3D Niggli test and whole-cell reduction to keep\n"
        L"the standard setting. With 0 solutions or a clearly high best rms, the UI warns \"system may be wrong\";\n"
        L"then compare against triclinic instead of switching algorithms (switching will not fix it)");
    s->Add(m_specSystem, 0, wxEXPAND | wxBOTTOM, 3);

    addSection(sw, s, L"Search Range");
    m_spRedundant = addSpin(sw, s, L"Cand/peak R", 4, 1, 8, L"Candidates kept per peak in Part TWO");
    m_spUv = addSpin(sw, s, L"Max (u,v)", 2, 0, 3);
    m_spHk1 = addSpin(sw, s, L"Part1 (h,k)", 3, 1, 6);
    m_spHk2 = addSpin(sw, s, L"Index (h,k)", 6, 1, 8);
    m_spL = addSpin(sw, s, L"l range", 6, 1, 8);
    m_spW = addSpin(sw, s, L"Max w", 2, 0, 4);
    m_spLines = addSpin(sw, s, L"Seed triads N", 4, 3, 8);

    addSection(sw, s, L"Candidate Size");
    m_spMaxPart1 = addSpin(sw, s, L"Max Part1", 100, 1, 200);
    m_spRefineTop = addSpin(sw, s, L"Refine top", 200, 10, 2000);
    m_spMaxSol = addSpin(sw, s, L"Max sols", 50, 1, 200);

    addSection(sw, s, L"Criteria");
    m_txtDqxy = addText(sw, s, L"Δq_xy tol", L"0.05");
    m_txtDqz = addText(sw, s, L"Δq_z tol", L"0.05");
    m_txtDqxyz = addText(sw, s, L"Δq_xyz tol", L"0.05");
    m_txtDqspec = addText(sw, s, L"Δq_spec tol", L"0.01");
    addCheck(sw, s, L"Niggli", true, m_chkNiggli);
    addCheck(sw, s, L"Limit cell range", false, m_chkRestrict);
    addRange(sw, s, L"a range (Å)", 3, 30, m_txtALo, m_txtAHi);
    addRange(sw, s, L"b range (Å)", 3, 30, m_txtBLo, m_txtBHi);
    addRange(sw, s, L"c range (Å)", 3, 30, m_txtCLo, m_txtCHi);
    addRange(sw, s, L"Vol range (Å³)", 50, 100000, m_txtVolLo, m_txtVolHi);

    addSection(sw, s, L"Misc");
    m_spThreads = addSpin(sw, s, L"Threads (0=auto)", 0, 0, 32);

    sw->SetSizer(s);
    s->FitInside(sw);
    sw->SetScrollRate(0, 10);
    return sw;
}

wxScrolledWindow* MainFrame::buildNoSpecPage(wxWindow* parent) {
    auto* sw = new wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                                    wxVSCROLL | wxTAB_TRAVERSAL);
    auto* s = new wxFlexGridSizer(0, 2, 4, 10);
    s->AddGrowableCol(1);

    addSection(sw, s, L"Solver");
    addLabel(sw, s, L"Method");
    m_nsMethod = new wxChoice(sw, wxID_ANY, wxDefaultPosition, wxSize(220, -1));
    m_nsMethod->Append(L"Beam search (stable, default)");
    m_nsMethod->Append(L"LLL + sphere decoding");
    m_nsMethod->SetSelection(static_cast<int>(gid::NoSpecMethod::Beam));
    m_nsMethod->SetToolTip(L"Beam search: per-peak branch and bound, stable; LLL + sphere decoding: LLL reduction, then assignment");
    s->Add(m_nsMethod, 0, wxEXPAND | wxBOTTOM, 3);

    addSection(sw, s, L"System");
    addLabel(sw, s, L"System");
    m_nsSystem = new wxChoice(sw, wxID_ANY, wxDefaultPosition, wxSize(120, -1));
    for (const std::string& n : gid::crystalSystemNames())
        m_nsSystem->Append(wxString::FromUTF8(n.c_str()));
    m_nsSystem->SetSelection(0);
    m_nsSystem->SetToolTip(L"Metric constraints for the target system (higher symmetry = stronger constraints, fewer solutions)");
    s->Add(m_nsSystem, 0, wxEXPAND | wxBOTTOM, 3);

    addSection(sw, s, L"Index Box (h,k,l)");
    m_nsSeedH = addSpin(sw, s, L"Max h", 2, 0, 6, L"h enumeration range ± for all peaks (including seeds)");
    m_nsSeedK = addSpin(sw, s, L"Max k", 2, 0, 6);
    m_nsSeedL = addSpin(sw, s, L"Max l", 6, 1, 12,
                        L"A long c axis (e.g. alkyl-chain stacking) needs a larger l, otherwise the true solution lies outside the search space");

    addSection(sw, s, L"Beam Expansion");
    m_nsSeedPeaks = addSpin(sw, s, L"Seed peaks N", 6, 3, 10,
                            L"Take any 3 of the N peaks with the smallest |g| as seeds");
    m_nsMinRows = addSpin(sw, s, L"Min rows", 8, 3, 10,
                          L">=8 recommended for triclinic; higher symmetry needs fewer rows");
    m_nsBeam = addSpin(sw, s, L"Beam width", 8000, 500, 200000, L"States kept per layer (larger = slower)");
    m_nsRefineTop = addSpin(sw, s, L"Refine top", 300, 10, 2000);
    m_nsMaxSol = addSpin(sw, s, L"Max sols", 50, 1, 200);

    addSection(sw, s, L"Joint Reassign Grid");
    m_nsGridH = addSpin(sw, s, L"h range", 6, 1, 10);
    m_nsGridK = addSpin(sw, s, L"k range", 6, 1, 10);
    m_nsGridL = addSpin(sw, s, L"l range", 6, 1, 10);

    addSection(sw, s, L"Criteria");
    m_nsDqxy = addText(sw, s, L"Δq_xy tol", L"0.05");
    m_nsDqz = addText(sw, s, L"Δq_z tol", L"0.01");
    m_nsDqxyz = addText(sw, s, L"Δq_xyz tol", L"0.05");
    addCheck(sw, s, L"Limit cell range", false, m_nsRestrict);
    addRange(sw, s, L"a range (Å)", 3, 30, m_nsALo, m_nsAHi);
    addRange(sw, s, L"b range (Å)", 3, 30, m_nsBLo, m_nsBHi);
    addRange(sw, s, L"c range (Å)", 3, 30, m_nsCLo, m_nsCHi);
    addRange(sw, s, L"Vol range (Å³)", 50, 100000, m_nsVolLo, m_nsVolHi);

    addSection(sw, s, L"Misc");
    m_nsThreads = addSpin(sw, s, L"Threads (0=auto)", 0, 0, 32);

    sw->SetSizer(s);
    s->FitInside(sw);
    sw->SetScrollRate(0, 10);
    return sw;
}

// Fourth mode: no q_spec (two-step). It has its own nospec2::Params: every search
// range, criterion and reduction switch of both steps is visible here; no hidden defaults from SCRM.
wxScrolledWindow* MainFrame::buildNoSpec2Page(wxWindow* parent) {
    auto* sw = new wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                                    wxVSCROLL | wxTAB_TRAVERSAL);
    auto* s = new wxFlexGridSizer(0, 2, 4, 10);
    s->AddGrowableCol(1);

    addSection(sw, s, L"q_spec 1D Search");
    m_ns2QsLo = addText(sw, s, L"q_spec min (0=auto)", L"0",
                        L"Leave 0 to infer from the peak table; without specular peaks, q_s becomes the 4th unknown of step 1");
    m_ns2QsHi = addText(sw, s, L"q_spec max (0=auto)", L"0");
    m_ns2Steps = addSpin(sw, s, L"Coarse steps", 12, 4, 64,
                         L"First coarse-scan the in-plane residual R(τ), then refine each local minimum by golden section");
    m_ns2Refine = addSpin(sw, s, L"Refine iters/min", 6, 0, 24,
                          L"Golden-section iterations; 0 = use coarse grid points only");
    m_ns2MaxTau = addSpin(sw, s, L"τ into step 2", 8, 1, 24,
                          L"Several τ values may give a good in-plane fit; all must be judged by step-2 integer consistency");

    addSection(sw, s, L"τ Coarse Budget");
    m_ns2CoarseHk1 = addSpin(sw, s, L"Coarse seed (h,k)", 2, 1, 8,
                             L"The coarse scan evaluates dozens of τ, so use a light budget; it does not affect the final range");
    m_ns2CoarseHk2 = addSpin(sw, s, L"Coarse all (h,k)", 4, 1, 12);
    m_ns2CoarseSeedPeaks = addSpin(sw, s, L"Coarse seed pool", 3, 1, 12);
    m_ns2CoarseRedundant = addSpin(sw, s, L"Coarse groups/peak", 2, 1, 8);
    m_ns2CoarseMaxPart1 = addSpin(sw, s, L"Coarse max cands", 40, 1, 2000,
                                  L"Set it equal to the main Step 1 setting to give the coarse scan the full budget");

    addSection(sw, s, L"Rod Peak Prior");
    m_ns2RodTol = addText(sw, s, L"q_xy thresh (0=off)", L"0.02",
                          L"Peaks with q_xy below this value are excluded from fitting; their q_z is used only as τ candidates");
    m_ns2RodOrders = addSpin(sw, s, L"Rod order m", 3, 1, 8,
                             L"Add q_z,rod/m with m=1..this value");
    m_ns2RodPrior = addText(sw, s, L"τ score factor", L"0.5",
                            L"<1 means early admission to step 2; first-order rod peaks get the smallest factor");

    addSection(sw, s, L"Step 1: In-plane");
    addLabel(sw, s, L"System");
    m_ns2System = new wxChoice(sw, wxID_ANY, wxDefaultPosition, wxSize(120, -1));
    for (const std::string& n : gid::crystalSystemNames())
        m_ns2System->Append(wxString::FromUTF8(n.c_str()));
    m_ns2System->SetSelection(0);
    m_ns2System->SetToolTip(L"Metric constraints for the target system (higher symmetry = stronger constraints, fewer solutions)");
    s->Add(m_ns2System, 0, wxEXPAND | wxBOTTOM, 3);
    m_ns2UvLimit = addSpin(sw, s, L"Max (u,v)", 2, 0, 5,
                           L"Enumeration radius of the contact-plane direction (u,v) (half-plane representative)");
    m_ns2Hk1 = addSpin(sw, s, L"Seed (h,k)", 3, 1, 8);
    m_ns2Hk2 = addSpin(sw, s, L"All-peak (h,k)", 6, 1, 12);
    m_ns2Seed3Half = addSpin(sw, s, L"First seed (0=any)", 0, 0, 12,
                             L"If >0, shrink the (h,k) radius of the first seed to this value");
    m_ns2SeedPeaks = addSpin(sw, s, L"Seed peak pool", 4, 1, 12);
    m_ns2Redundant = addSpin(sw, s, L"Groups/peak", 4, 1, 8);
    m_ns2MaxPart1 = addSpin(sw, s, L"Step1 max cands", 200, 1, 2000,
                            L"The coarse scan and final stage use the same range parameters");
    m_ns2Dqxy = addText(sw, s, L"Δq_xy tol", L"0.05");

    addSection(sw, s, L"Step 2: Out-of-plane");
    m_ns2SeedL = addSpin(sw, s, L"Seed l radius", 1, 0, 6,
                         L"l enumeration range ± this value for the three seed peaks");
    m_ns2SeedLFallback = addSpin(sw, s, L"Fallback l (0=off)", 3, 0, 6,
                                 L"If the first τ yields no good solution, widen the seed l radius to this value");
    m_ns2LLimit = addSpin(sw, s, L"Max |l|/peak", 6, 1, 20);
    m_ns2WLimit = addSpin(sw, s, L"Max |w|", 3, 1, 12);
    m_ns2MaxTrials = addSpin(sw, s, L"3×3 trials (0=inf)", 0, 0, 200000);
    m_ns2RefineTop = addSpin(sw, s, L"Refine top", 200, 1, 2000);
    m_ns2Dqz = addText(sw, s, L"Δq_z tol", L"0.05");

    addSection(sw, s, L"Cell Criteria");
    addCheck(sw, s, L"Limit cell range", false, m_ns2Restrict);
    addRange(sw, s, L"a range (Å)", 3, 30, m_ns2ALo, m_ns2AHi);
    addRange(sw, s, L"b range (Å)", 3, 30, m_ns2BLo, m_ns2BHi);
    addRange(sw, s, L"c range (Å)", 3, 30, m_ns2CLo, m_ns2CHi);
    addRange(sw, s, L"α range (°)", 60, 120, m_ns2AlLo, m_ns2AlHi);
    addRange(sw, s, L"β range (°)", 60, 120, m_ns2BeLo, m_ns2BeHi);
    addRange(sw, s, L"γ range (°)", 60, 120, m_ns2GaLo, m_ns2GaHi);
    addRange(sw, s, L"Vol range (Å³)", 50, 100000, m_ns2VolLo, m_ns2VolHi);
    addCheck(sw, s, L"Niggli", true, m_ns2Niggli);

    addSection(sw, s, L"Scale & Output");
    m_ns2IpCands = addSpin(sw, s, L"In-plane/τ", 16, 1, 64);
    m_ns2MaxSol = addSpin(sw, s, L"Max sols", 50, 1, 200);
    m_ns2Dqxyz = addText(sw, s, L"Δq_xyz tol", L"0.05");
    m_ns2Threads = addSpin(sw, s, L"Threads (0=auto)", 0, 0, 32);

    auto* note = new wxStaticText(
        sw, wxID_ANY,
        L"These settings are the actual search ranges used by steps 1 and 2; there are no hidden defaults.\n"
        L"Rod peaks with q_xy≈0 act only as soft τ-candidate priors and are not fitted.");
    s->Add(note, 0, wxEXPAND | wxTOP | wxBOTTOM, 4);
    s->AddSpacer(1);

    sw->SetSizer(s);
    s->FitInside(sw);
    sw->SetScrollRate(0, 10);
    return sw;
}

wxScrolledWindow* MainFrame::buildFixedPage(wxWindow* parent) {
    auto* sw = new wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                                    wxVSCROLL | wxTAB_TRAVERSAL);
    auto* s = new wxFlexGridSizer(0, 2, 4, 10);
    s->AddGrowableCol(1);

    addSection(sw, s, L"Fixed Contact Plane");
    m_spU = addSpin(sw, s, L"u", 0, -8, 8);
    m_spV = addSpin(sw, s, L"v", 0, -8, 8);
    m_spWplane = addSpin(sw, s, L"w", 1, -8, 8);

    addSection(sw, s, L"System & Search");
    addLabel(sw, s, L"System");
    m_chSystem = new wxChoice(sw, wxID_ANY, wxDefaultPosition, wxSize(120, -1));
    for (const std::string& n : gid::crystalSystemNames())
        m_chSystem->Append(wxString::FromUTF8(n.c_str()));
    m_chSystem->SetSelection(0);
    s->Add(m_chSystem, 0, wxEXPAND | wxBOTTOM, 3);
    m_spMaxH = addSpin(sw, s, L"h range", 4, 1, 8);
    m_spMaxK = addSpin(sw, s, L"k range", 4, 1, 8);
    m_spMaxL = addSpin(sw, s, L"l range", 4, 1, 8);

    addSection(sw, s, L"Initial Cell");
    m_txtSa = addText(sw, s, L"a (Å)", L"10");
    m_txtSb = addText(sw, s, L"b (Å)", L"10");
    m_txtSc = addText(sw, s, L"c (Å)", L"15");
    m_txtSal = addText(sw, s, L"α (°)", L"90");
    m_txtSbe = addText(sw, s, L"β (°)", L"90");
    m_txtSga = addText(sw, s, L"γ (°)", L"90");

    addSection(sw, s, L"Deterministic Search");
    m_spSeed = addSpin(sw, s, L"Seed peaks", 8, 3, 12);
    m_spTrials = addSpin(sw, s, L"Max cands (10k)", 40, 1, 400);
    m_txtAccept = addText(sw, s, L"Accept tol (Å⁻¹)", L"0.05");
    m_spRefIter = addSpin(sw, s, L"Alt refine iters", 6, 1, 20);
    m_spTopN = addSpin(sw, s, L"Keep sols", 10, 1, 50);
    // Note: q_xy sign only reflects mirror handedness of the detection geometry; the core
    // already takes |q_xy| at data entry, so this switch no longer affects results. Kept as-is.
    addCheck(sw, s, L"Auto q_xy sign", true, m_chkAutoSign);

    sw->SetSizer(s);
    s->FitInside(sw);
    sw->SetScrollRate(0, 10);
    return sw;
}

wxScrolledWindow* MainFrame::buildSimplePage(wxWindow* parent) {
    auto* sw = new wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                                    wxVSCROLL | wxTAB_TRAVERSAL);
    auto* s = new wxFlexGridSizer(0, 2, 4, 10);
    s->AddGrowableCol(1);

    auto* note = new wxStaticText(
        sw, wxID_ANY,
        L"Simple indexing (out-of-plane fixed to (001), searched automatically).\n"
        L"Model: q_xy² = h²A + k²B − 2hkC, q_z = hP + kQ + lR.\n"
        L"In-plane hypotheses come from exact three-peak solutions + all-peak consensus; the out-of-plane\n"
        L"period R is estimated from specular rod peaks, tilt (P,Q) is found by torus Hough search, and the\n"
        L"result is refined by reassignment + linear least squares. Only 3 non-rod peaks are needed; q_spec is optional.");
    s->Add(note, 0, wxEXPAND | wxTOP | wxBOTTOM, 4);
    s->AddSpacer(1);

    addSection(sw, s, L"Match Tolerance (Å⁻¹)");
    m_simpleTol = addText(sw, s, L"Model tolerance", L"0.05",
                          L"2-D (q_xy,q_z) nearest-neighbor matching threshold; decides whether a peak is an inlier");
    m_simpleTolXy = addText(sw, s, L"In-plane Δq_xy", L"0.05", L"q_xy tolerance used for in-plane hypotheses and consensus");
    m_simpleTolZ = addText(sw, s, L"Out-plane Δq_z", L"0.05", L"q_z tolerance used by the tilt Hough search and l assignment");
    m_simpleRod = addText(sw, s, L"Rod q_xy thresh", L"0.04",
                          L"Peaks with q_xy below this value count as specular rod peaks and are used to estimate the out-of-plane period R");

    addSection(sw, s, L"Search Range");
    m_simpleHk = addSpin(sw, s, L"Max (h,k)", 2, 1, 12,
                         L"In-plane index box |h|,|k| <= this value. The staged search starts at 2 and "
                         L"stops when good enough, so raising it only costs time when needed");
    m_simpleL = addSpin(sw, s, L"Max l", 6, 1, 12, L"Out-of-plane order box |l| <= this value; also the rod order cap");

    addSection(sw, s, L"Scale & Output");
    m_simpleTopK = addSpin(sw, s, L"Max sols", 10, 1, 50);
    m_simpleGrid = addSpin(sw, s, L"Tilt Hough grid", 33, 9, 81, L"Edge length of the torus (P/R,Q/R) search grid (larger = slower)");

    sw->SetSizer(s);
    s->FitInside(sw);
    sw->SetScrollRate(0, 10);
    return sw;
}

// ---------------------------------------------------------------- Event bindings
void MainFrame::bindEvents() {
    Bind(wxEVT_CLOSE_WINDOW, &MainFrame::onClose, this);
    Bind(wxEVT_MENU, &MainFrame::onImport, this, ID_IMPORT);
    Bind(wxEVT_TOOL, &MainFrame::onImport, this, ID_IMPORT);
    Bind(wxEVT_MENU, &MainFrame::onExport, this, ID_EXPORT);
    Bind(wxEVT_TOOL, &MainFrame::onExport, this, ID_EXPORT);
    Bind(wxEVT_MENU, &MainFrame::onRun, this, ID_RUN);
    Bind(wxEVT_TOOL, &MainFrame::onRun, this, ID_RUN);
    Bind(wxEVT_MENU, &MainFrame::onStop, this, ID_STOP);
    Bind(wxEVT_TOOL, &MainFrame::onStop, this, ID_STOP);
    Bind(wxEVT_MENU, &MainFrame::onClear, this, ID_CLEAR);
    Bind(wxEVT_TOOL, &MainFrame::onClear, this, ID_CLEAR);
    Bind(wxEVT_MENU, &MainFrame::onAppend, this, ID_APPEND);
    Bind(wxEVT_MENU, &MainFrame::onDelete, this, ID_DELETE);
    Bind(wxEVT_BUTTON, &MainFrame::onAutoQSpec, this, ID_AUTO_QSPEC);
    Bind(wxEVT_MENU, &MainFrame::onModeMenu, this, ID_MODE_SPEC);
    Bind(wxEVT_MENU, &MainFrame::onModeMenu, this, ID_MODE_NOSPEC);
    Bind(wxEVT_MENU, &MainFrame::onModeMenu, this, ID_MODE_NOSPEC2);
    Bind(wxEVT_MENU, &MainFrame::onModeMenu, this, ID_MODE_FIXED);
    Bind(wxEVT_MENU, &MainFrame::onModeMenu, this, ID_MODE_SIMPLE);
    Bind(wxEVT_CHOICE, &MainFrame::onModeChoice, this, ID_MODE_CHOICE);
    Bind(wxEVT_MENU, &MainFrame::onUsage, this, ID_USAGE);
    Bind(wxEVT_MENU, &MainFrame::onAbout, this, ID_ABOUT);
    Bind(wxEVT_MENU, [this](wxCommandEvent&) { Close(true); }, wxID_EXIT);
    Bind(wxEVT_GRID_SELECT_CELL, &MainFrame::onSolSelect, this);
}

void MainFrame::onClose(wxCloseEvent& e) {
    if (m_running) {
        m_cancel = true;
        SetStatusText(L"Stopping computation…", 0);
    }
    if (m_worker.joinable()) m_worker.join();
    e.Skip();
}

// ---------------------------------------------------------------- Data
void MainFrame::importFile(const wxString& path) {
    Mat M;
    std::string err;
    if (!loadDataTable(std::string(path.utf8_str()), M, err)) {
        wxMessageBox(wxString::FromUTF8(err.c_str()), L"Import Failed", wxICON_ERROR | wxOK, this);
        return;
    }
    if (M.c < 2 || M.r < 1) {
        wxMessageBox(L"Data must have at least two columns (q_xy, q_z) and one row", L"Import Failed", wxICON_ERROR | wxOK, this);
        return;
    }
    removeRowsAnyNonFinite(M);

    int cx = 0, cz = 1;
    if (M.c > 2) {
        wxArrayString names;
        for (int j = 0; j < M.c && j < 16; ++j)
            names.Add(wxString::Format(L"Column %d   (first value %.4g)", j + 1, M(0, j)));
        wxSingleChoiceDialog d1(this, L"Select the q_xy column", L"Import Peaks", names);
        d1.SetSelection(0);
        if (d1.ShowModal() != wxID_OK) return;
        cx = d1.GetSelection();
        wxSingleChoiceDialog d2(this, L"Select the q_z column", L"Import Peaks", names);
        d2.SetSelection(M.c > 1 ? 1 : 0);
        if (d2.ShowModal() != wxID_OK) return;
        cz = d2.GetSelection();
    }

    m_peaks.clear();
    for (int i = 0; i < M.r; ++i) {
        const double qx = M(i, cx), qz = M(i, cz);
        if (!std::isfinite(qx) || !std::isfinite(qz)) continue;
        m_peaks.push_back(gid::Peak{qx, qz});
    }
    // q_xy sign only reflects mirror handedness of the detection geometry: take |q_xy| on import,
    // so computation, peak-table display, delta columns and export all use the same normalized data.
    gid::normalizeQxySigns(m_peaks);
    m_solutions.clear();
    m_solRow = -1;
    refreshPeakGrid();
    m_plot->setPeaks(m_peaks);
    m_plot->setPredicted({});
    m_solGrid->DeleteRows(0, m_solGrid->GetNumberRows());
    m_solInfo->SetLabel(L"Not computed");

    autoDetectQSpec();
    SetStatusText(wxString::Format(L"Imported %zu peaks: %s", m_peaks.size(), path), 0);
    SetStatusText(wxString::Format(L"Peaks: %zu", m_peaks.size()), 1);
}

void MainFrame::refreshPeakGrid() {
    wxGrid* g = m_peakGrid;
    const int n = (int)m_peaks.size();
    if (g->GetNumberRows() > n) g->DeleteRows(n, g->GetNumberRows() - n);
    if (g->GetNumberRows() < n) g->AppendRows(n - g->GetNumberRows());
    for (int i = 0; i < n; ++i) {
        g->SetCellValue(i, 0, wxString::Format("%d", i + 1));
        g->SetCellValue(i, 1, wxString::Format("%.5f", m_peaks[i].qxy));
        g->SetCellValue(i, 2, wxString::Format("%.5f", m_peaks[i].qz));
        for (int j = 3; j < 8; ++j) g->SetCellValue(i, j, "");
        g->SetReadOnly(i, 0, true);
        for (int j = 3; j < 8; ++j) g->SetReadOnly(i, j, true);
    }
    g->ClearSelection();
}

bool MainFrame::collectPeaks() {
    m_peaks.clear();
    wxGrid* g = m_peakGrid;
    for (int i = 0; i < g->GetNumberRows(); ++i) {
        double qx = 0.0, qz = 0.0;
        const bool okx = g->GetCellValue(i, 1).ToDouble(&qx);
        const bool okz = g->GetCellValue(i, 2).ToDouble(&qz);
        if (!okx || !okz) continue;
        if (!std::isfinite(qx) || !std::isfinite(qz)) continue;
        m_peaks.push_back(gid::Peak{qx, qz});
    }
    // Hand-entered/edited peaks are also normalized to |q_xy| before computation; the table is refreshed to match.
    gid::normalizeQxySigns(m_peaks);
    refreshPeakGrid();
    if (m_peaks.empty()) {
        wxMessageBox(L"Peak table is empty: import data or enter q_xy / q_z directly", L"Cannot Compute",
                     wxICON_WARNING | wxOK, this);
        return false;
    }
    return true;
}

void MainFrame::appendPeakRow() {
    m_peakGrid->AppendRows(1);
    const int r = m_peakGrid->GetNumberRows() - 1;
    m_peakGrid->SetCellValue(r, 0, wxString::Format("%d", r + 1));
    m_peakGrid->SetReadOnly(r, 0, true);
    for (int j = 3; j < 8; ++j) m_peakGrid->SetReadOnly(r, j, true);
    m_peakGrid->GoToCell(r, 1);
}

void MainFrame::deleteSelectedPeakRows() {
    wxGrid* g = m_peakGrid;
    wxArrayInt rows = g->GetSelectedRows();
    int r = g->GetGridCursorRow();
    if (rows.empty() && r >= 0) rows.Add(r);
    std::sort(rows.begin(), rows.end(), std::greater<int>());
    for (int i : rows)
        if (i >= 0 && i < g->GetNumberRows()) g->DeleteRows(i, 1);
    if (!collectPeaks()) return;
    m_plot->setPeaks(m_peaks);
    SetStatusText(wxString::Format(L"Peaks: %zu", m_peaks.size()), 1);
}

void MainFrame::clearPeaks() {
    m_peaks.clear();
    m_solutions.clear();
    m_solRow = -1;
    refreshPeakGrid();
    m_solGrid->DeleteRows(0, m_solGrid->GetNumberRows());
    m_solInfo->SetLabel(L"Not computed");
    m_plot->setPeaks({});
    m_plot->setPredicted({});
    SetStatusText(L"Peak table cleared", 0);
    SetStatusText(L"Peaks: 0", 1);
}

void MainFrame::autoDetectQSpec() {
    double qs = 0.0;
    int rods = 0;
    if (gid::detectQSpec(m_peaks, qs, rods)) {
        m_txtQSpec->SetValue(wxString::Format("%.5f", qs));
        SetStatusText(wxString::Format(L"Auto-detected q_spec = %.5f (%d specular rod peaks)", qs, rods), 0);
    } else {
        wxMessageBox(L"No specular rod peak with q_xy = 0 found; enter q_spec manually", L"Auto-detect",
                     wxICON_INFORMATION | wxOK, this);
    }
}

// ---------------------------------------------------------------- Parameters
bool MainFrame::collectParams(gid::EngineParams& ep) {
    ep.mode = m_mode;
    if (m_mode == gid::IndexMode::WithSpec) {
        gid::SpecParams& s = ep.spec;
        double v = 0.0;
        if (!m_txtQSpec->GetValue().ToDouble(&v) || !(v > 0.0)) {
            wxMessageBox(L"Enter q_spec first (click \"Auto\" to get it from the specular rod peaks)", L"Missing Parameter",
                         wxICON_WARNING | wxOK, this);
            return false;
        }
        s.qspec = v;
        if (m_chP1Alg) {
            const int sel = m_chP1Alg->GetSelection();
            s.part1Alg = (sel >= 0 && sel < (int)(sizeof(kP1AlgItems) / sizeof(kP1AlgItems[0])))
                             ? kP1AlgItems[sel].alg
                             : gid::Part1Algorithm::PeriodVoteDdl;
        } else {
            s.part1Alg = gid::Part1Algorithm::PeriodVoteDdl;
        }
        if (m_specSystem && m_specSystem->GetSelection() >= 0)
            s.system = static_cast<gid::CrystalSystem>(m_specSystem->GetSelection());
        s.redundant = m_spRedundant->GetValue();
        s.uvLimit = m_spUv->GetValue();
        s.hkLimit1 = m_spHk1->GetValue();
        s.hkLimit2 = m_spHk2->GetValue();
        s.lLimit = m_spL->GetValue();
        s.wLimit = m_spW->GetValue();
        s.linesToPermute = m_spLines->GetValue();
        s.maxPart1 = m_spMaxPart1->GetValue();
        s.refineTop = m_spRefineTop->GetValue();
        s.maxSolutions = m_spMaxSol->GetValue();
        s.niggli = m_chkNiggli->GetValue();
        s.restrictCell = m_chkRestrict->GetValue();
        if (m_txtDqxy->GetValue().ToDouble(&v)) s.dqxyCutoff = v;
        if (m_txtDqz->GetValue().ToDouble(&v)) s.dqzCutoff = v;
        if (m_txtDqxyz->GetValue().ToDouble(&v)) s.dqxyzCutoff = v;
        if (m_txtDqspec->GetValue().ToDouble(&v)) s.dqspec = v;
        if (s.restrictCell) {
            if (m_txtALo->GetValue().ToDouble(&v)) s.aLo = v;
            if (m_txtAHi->GetValue().ToDouble(&v)) s.aHi = v;
            if (m_txtBLo->GetValue().ToDouble(&v)) s.bLo = v;
            if (m_txtBHi->GetValue().ToDouble(&v)) s.bHi = v;
            if (m_txtCLo->GetValue().ToDouble(&v)) s.cLo = v;
            if (m_txtCHi->GetValue().ToDouble(&v)) s.cHi = v;
            if (m_txtVolLo->GetValue().ToDouble(&v)) s.volLo = v;
            if (m_txtVolHi->GetValue().ToDouble(&v)) s.volHi = v;
        }
        s.threads = m_spThreads->GetValue();
    } else if (m_mode == gid::IndexMode::NoSpec) {
        gid::NoSpecParams& n = ep.nospec;
        n.method = static_cast<gid::NoSpecMethod>(m_nsMethod->GetSelection());
        n.system = static_cast<gid::CrystalSystem>(m_nsSystem->GetSelection());
        n.seedH = m_nsSeedH->GetValue();
        n.seedK = m_nsSeedK->GetValue();
        n.seedL = m_nsSeedL->GetValue();
        n.seedPeaks = m_nsSeedPeaks->GetValue();
        n.minRows = m_nsMinRows->GetValue();
        n.beam = m_nsBeam->GetValue();
        n.refineTop = m_nsRefineTop->GetValue();
        n.maxSolutions = m_nsMaxSol->GetValue();
        n.gridH = m_nsGridH->GetValue();
        n.gridK = m_nsGridK->GetValue();
        n.gridL = m_nsGridL->GetValue();
        double v = 0.0;
        if (m_nsDqxy->GetValue().ToDouble(&v)) n.dqxyCutoff = v;
        if (m_nsDqz->GetValue().ToDouble(&v)) n.dqzCutoff = v;
        if (m_nsDqxyz->GetValue().ToDouble(&v)) n.dqxyzCutoff = v;
        n.restrictCell = m_nsRestrict->GetValue();
        if (n.restrictCell) {
            if (m_nsALo->GetValue().ToDouble(&v)) n.aLo = v;
            if (m_nsAHi->GetValue().ToDouble(&v)) n.aHi = v;
            if (m_nsBLo->GetValue().ToDouble(&v)) n.bLo = v;
            if (m_nsBHi->GetValue().ToDouble(&v)) n.bHi = v;
            if (m_nsCLo->GetValue().ToDouble(&v)) n.cLo = v;
            if (m_nsCHi->GetValue().ToDouble(&v)) n.cHi = v;
            if (m_nsVolLo->GetValue().ToDouble(&v)) n.volLo = v;
            if (m_nsVolHi->GetValue().ToDouble(&v)) n.volHi = v;
        }
        n.threads = m_nsThreads->GetValue();
        if (n.threads == 0) n.threads = -1; // 0 = auto (let parallelFor pick the hardware thread count)
    } else if (m_mode == gid::IndexMode::NoSpecTwoSteps) {
        gid::nospec2::Params& n2 = ep.nospec2;
        double v2 = 0.0;
        auto rd = [&](wxTextCtrl* c, double& out) {
            if (c->GetValue().ToDouble(&v2)) out = v2;
        };
        // q_spec 1D search
        rd(m_ns2QsLo, n2.qsLo);
        rd(m_ns2QsHi, n2.qsHi);
        n2.qsSteps = m_ns2Steps->GetValue();
        n2.qsRefine = m_ns2Refine->GetValue();
        n2.maxTau = m_ns2MaxTau->GetValue();
        n2.coarseHk1 = m_ns2CoarseHk1->GetValue();
        n2.coarseHk2 = m_ns2CoarseHk2->GetValue();
        n2.coarseSeedPeaks = m_ns2CoarseSeedPeaks->GetValue();
        n2.coarseRedundant = m_ns2CoarseRedundant->GetValue();
        n2.coarseMaxPart1 = m_ns2CoarseMaxPart1->GetValue();
        rd(m_ns2RodTol, n2.rodQxyTol);
        n2.rodOrders = m_ns2RodOrders->GetValue();
        rd(m_ns2RodPrior, n2.rodPrior);
        // Step 1: in-plane enumeration
        n2.system = static_cast<gid::CrystalSystem>(m_ns2System->GetSelection());
        n2.uvLimit = m_ns2UvLimit->GetValue();
        n2.hk1 = m_ns2Hk1->GetValue();
        n2.hk2 = m_ns2Hk2->GetValue();
        n2.seed3Half = m_ns2Seed3Half->GetValue();
        n2.seedPeaks = m_ns2SeedPeaks->GetValue();
        n2.redundant = m_ns2Redundant->GetValue();
        n2.maxPart1 = m_ns2MaxPart1->GetValue();
        rd(m_ns2Dqxy, n2.dqxyCutoff);
        // Step 2: out-of-plane enumeration
        n2.seedL = m_ns2SeedL->GetValue();
        n2.seedLFallback = m_ns2SeedLFallback->GetValue();
        n2.lLimit = m_ns2LLimit->GetValue();
        n2.wLimit = m_ns2WLimit->GetValue();
        n2.maxTrials = m_ns2MaxTrials->GetValue();
        n2.refineTop = m_ns2RefineTop->GetValue();
        rd(m_ns2Dqz, n2.dqzCutoff);
        // Cell criteria and reduction
        n2.restrictCell = m_ns2Restrict->GetValue();
        rd(m_ns2ALo, n2.aLo);   rd(m_ns2AHi, n2.aHi);
        rd(m_ns2BLo, n2.bLo);   rd(m_ns2BHi, n2.bHi);
        rd(m_ns2CLo, n2.cLo);   rd(m_ns2CHi, n2.cHi);
        rd(m_ns2AlLo, n2.alLo); rd(m_ns2AlHi, n2.alHi);
        rd(m_ns2BeLo, n2.beLo); rd(m_ns2BeHi, n2.beHi);
        rd(m_ns2GaLo, n2.gaLo); rd(m_ns2GaHi, n2.gaHi);
        rd(m_ns2VolLo, n2.volLo); rd(m_ns2VolHi, n2.volHi);
        n2.niggli = m_ns2Niggli->GetValue();
        // Scale and output
        n2.maxIpCands = m_ns2IpCands->GetValue();
        n2.maxSolutions = m_ns2MaxSol->GetValue();
        rd(m_ns2Dqxyz, n2.dqxyzCutoff);
        n2.threads = m_ns2Threads->GetValue();
        if (n2.threads == 0) n2.threads = -1;
    } else if (m_mode == gid::IndexMode::Simple) {
        gid::SimpleParams& sp = ep.simple;
        double v = 0.0;
        if (m_simpleTol->GetValue().ToDouble(&v) && v > 0.0) sp.tol = v;
        if (m_simpleTolXy->GetValue().ToDouble(&v) && v > 0.0) sp.tolXy = v;
        if (m_simpleTolZ->GetValue().ToDouble(&v) && v > 0.0) sp.tolZ = v;
        if (m_simpleRod->GetValue().ToDouble(&v) && v >= 0.0) sp.rodTol = v;
        sp.hkLimit = m_simpleHk->GetValue();
        sp.lLimit = m_simpleL->GetValue();
        sp.topK = m_simpleTopK->GetValue();
        sp.grid = m_simpleGrid->GetValue();
    } else {
        gid::FixedPlaneParams& f = ep.fixed;
        f.u = m_spU->GetValue();
        f.v = m_spV->GetValue();
        f.w = m_spWplane->GetValue();
        if (f.u == 0 && f.v == 0 && f.w == 0) {
            wxMessageBox(L"(u,v,w) cannot all be zero", L"Invalid Parameters", wxICON_WARNING | wxOK, this);
            return false;
        }
        f.system = static_cast<gid::CrystalSystem>(m_chSystem->GetSelection());
        f.maxH = m_spMaxH->GetValue();
        f.maxK = m_spMaxK->GetValue();
        f.maxL = m_spMaxL->GetValue();
        double v = 0.0;
        if (m_txtSa->GetValue().ToDouble(&v)) f.start.a = v;
        if (m_txtSb->GetValue().ToDouble(&v)) f.start.b = v;
        if (m_txtSc->GetValue().ToDouble(&v)) f.start.c = v;
        if (m_txtSal->GetValue().ToDouble(&v)) f.start.al = v;
        if (m_txtSbe->GetValue().ToDouble(&v)) f.start.be = v;
        if (m_txtSga->GetValue().ToDouble(&v)) f.start.ga = v;
        if (!(f.start.a > 0) || !(f.start.b > 0) || !(f.start.c > 0)) {
            wxMessageBox(L"Initial cell a, b, c must be positive", L"Invalid Parameters", wxICON_WARNING | wxOK, this);
            return false;
        }
        // The two former "random trials" controls now map to FP-DL: seed peak count / candidate cap
        f.seedPeaks = m_spSeed->GetValue();
        f.maxCands = m_spTrials->GetValue() * 10000;
        if (m_txtAccept->GetValue().ToDouble(&v)) f.acceptTol = v;
        f.refineIter = m_spRefIter->GetValue();
        f.topN = m_spTopN->GetValue();
        f.autoSign = m_chkAutoSign->GetValue();
        f.threads = m_spThreads->GetValue();
    }
    return true;
}

void MainFrame::syncModeUI() {
    const bool spec = (m_mode == gid::IndexMode::WithSpec);
    const bool nospec = (m_mode == gid::IndexMode::NoSpec);
    const bool nospec2 = (m_mode == gid::IndexMode::NoSpecTwoSteps);
    const bool simple = (m_mode == gid::IndexMode::Simple);
    const bool fixed = !spec && !nospec && !nospec2 && !simple;
    m_pageSpec->Show(spec);
    m_pageNoSpec->Show(nospec);
    m_pageNoSpec2->Show(nospec2);
    m_pageFixed->Show(fixed);
    m_pageSimple->Show(simple);
    // All five pages share one sizer and only one is shown: after switching, the shown page's
    // inner layout and scroll extent must be recomputed, or it keeps its previous size and misaligns.
    wxScrolledWindow* shown =
        spec ? m_pageSpec
             : (nospec ? m_pageNoSpec
                       : (nospec2 ? m_pageNoSpec2 : (simple ? m_pageSimple : m_pageFixed)));
    m_paramHost->Layout();
    if (auto* ps = shown->GetSizer()) {
        ps->Layout();
        ps->FitInside(shown);
    }
    shown->Scroll(0, 0);
    m_modeChoice->SetSelection(spec ? 0 : (nospec ? 1 : (nospec2 ? 2 : (simple ? 4 : 3))));
    if (m_menuRun) {
        m_menuRun->Check(ID_MODE_SPEC, spec);
        m_menuRun->Check(ID_MODE_NOSPEC, nospec);
        m_menuRun->Check(ID_MODE_NOSPEC2, nospec2);
        m_menuRun->Check(ID_MODE_FIXED, fixed);
        m_menuRun->Check(ID_MODE_SIMPLE, simple);
    }
    SetStatusText(wxString::FromUTF8(gid::modeName(m_mode)), 2);
}

void MainFrame::onModeChoice(wxCommandEvent&) {
    switch (m_modeChoice->GetSelection()) {
        case 0: m_mode = gid::IndexMode::WithSpec; break;
        case 1: m_mode = gid::IndexMode::NoSpec; break;
        case 2: m_mode = gid::IndexMode::NoSpecTwoSteps; break;
        case 4: m_mode = gid::IndexMode::Simple; break;
        default: m_mode = gid::IndexMode::FixedPlane; break;
    }
    syncModeUI();
}

void MainFrame::onModeMenu(wxCommandEvent& e) {
    switch (e.GetId()) {
        case ID_MODE_SPEC: m_mode = gid::IndexMode::WithSpec; break;
        case ID_MODE_NOSPEC: m_mode = gid::IndexMode::NoSpec; break;
        case ID_MODE_NOSPEC2: m_mode = gid::IndexMode::NoSpecTwoSteps; break;
        case ID_MODE_FIXED: m_mode = gid::IndexMode::FixedPlane; break;
        case ID_MODE_SIMPLE: m_mode = gid::IndexMode::Simple; break;
        default: return;
    }
    syncModeUI();
}

// ---------------------------------------------------------------- Compute
void MainFrame::startIndexing() {
    if (m_running) return;
    if (!collectPeaks()) return;
    gid::EngineParams ep;
    if (!collectParams(ep)) return;

    m_lastSystem = gid::CrystalSystem::Triclinic;
    m_lastDqxyz = -1.0;
    if (m_mode == gid::IndexMode::WithSpec) {
        m_lastSystem = ep.spec.system;
        m_lastDqxyz = ep.spec.dqxyzCutoff;
    }

    m_cancel = false;
    setRunning(true);
    const std::vector<gid::Peak> peaks = m_peaks;

    m_worker = std::thread([this, peaks, ep]() {
        gid::Progress prog;
        prog.report = [this](double v, const std::string& msg) {
            const long long now = std::chrono::duration_cast<std::chrono::milliseconds>(
                                      std::chrono::steady_clock::now().time_since_epoch())
                                      .count();
            long long prev = m_lastReport.load();
            if (now - prev < 60) return;
            if (!m_lastReport.compare_exchange_strong(prev, now)) return;
            CallAfter([this, v, msg]() {
                m_gauge->SetValue((int)(v * 100.0 + 0.5));
                SetStatusText(wxString::FromUTF8(msg.c_str()), 0);
            });
        };
        prog.cancelled = [this]() { return m_cancel.load(); };

        const gid::EngineResult res = gid::runIndexing(peaks, ep, prog);
        CallAfter([this, res]() { onFinished(res); });
    });
}

void MainFrame::stopIndexing() {
    if (!m_running) return;
    m_cancel = true;
    SetStatusText(L"Stopping…", 0);
}

void MainFrame::onFinished(const gid::EngineResult& res) {
    if (m_worker.joinable()) m_worker.join();
    setRunning(false);
    m_gauge->SetValue(res.ok ? 100 : 0);

    m_solutions = res.solutions;
    refreshSolutionGrid();

    if (!res.ok) {
        SetStatusText(wxString::FromUTF8(res.message.c_str()), 0);
        if (!res.stopped)
            wxMessageBox(wxString::FromUTF8(res.message.c_str()), L"Indexing Finished",
                         wxICON_INFORMATION | wxOK, this);
        return;
    }
    SetStatusText(wxString::FromUTF8(res.message.c_str()), 0);
    if (!m_solutions.empty()) {
        m_solGrid->SelectRow(0);
        showSolution(0);
    }
    // Warn explicitly when the chosen system does not match the data (instead of silently
    // returning an empty list): 0 solutions, or the best solution's rmsQxyz exceeds half the tolerance.
    if (m_lastDqxyz > 0.0 && m_lastSystem != gid::CrystalSystem::Triclinic && !res.stopped) {
        const wxString sysName =
            wxString::FromUTF8(gid::crystalSystemName(m_lastSystem).c_str());
        if (m_solutions.empty()) {
            wxMessageBox(
                wxString::Format(
                    L"No solution satisfies the criteria for the selected system (%s).\n\n"
                    L"This usually means the data does not match that system's symmetry, or the criteria are too strict.\n"
                    L"Try \"Triclinic\" once: if triclinic gives a clearly better solution, the system choice is wrong.",
                    sysName),
                L"System May Not Match", wxICON_WARNING | wxOK, this);
        } else if (!(m_solutions.front().rmsQxyz <= 0.5 * m_lastDqxyz)) {
            wxMessageBox(
                wxString::Format(
                    L"Best solution for the selected system (%s) has rmsQxyz = %.4f, over half the tolerance (%.4f),\n"
                    L"i.e. a forced fit. If triclinic yields a clearly lower rms, the system does not match the data.",
                    sysName, m_solutions.front().rmsQxyz, m_lastDqxyz),
                L"System May Not Match", wxICON_WARNING | wxOK, this);
        }
    }
}

void MainFrame::setRunning(bool on) {
    m_running = on;
    m_toolbar->EnableTool(ID_RUN, !on);
    m_toolbar->EnableTool(ID_STOP, on);
    m_toolbar->EnableTool(ID_IMPORT, !on);
    m_toolbar->EnableTool(ID_CLEAR, !on);
    m_modeChoice->Enable(!on);
    m_gauge->SetValue(0);
}

// ---------------------------------------------------------------- Results
void MainFrame::refreshSolutionGrid() {
    wxGrid* g = m_solGrid;
    const int n = (int)m_solutions.size();
    if (g->GetNumberRows() > n) g->DeleteRows(n, g->GetNumberRows() - n);
    if (g->GetNumberRows() < n) g->AppendRows(n - g->GetNumberRows());
    for (int i = 0; i < n; ++i) {
        const gid::Solution& s = m_solutions[i];
        g->SetCellValue(i, 0, wxString::Format("%d", i + 1));
        // P2-4: do not show a fabricated (u v w) / q_spec when there is no integer contact plane
        g->SetCellValue(i, 1, s.hasPlane
                                  ? wxString::Format("(%d %d %d)", s.uvw.h, s.uvw.k, s.uvw.l)
                                  : wxString(L"—"));
        g->SetCellValue(i, 2, wxString::Format("%.4f", s.cell.a));
        g->SetCellValue(i, 3, wxString::Format("%.4f", s.cell.b));
        g->SetCellValue(i, 4, wxString::Format("%.4f", s.cell.c));
        g->SetCellValue(i, 5, wxString::Format("%.3f", s.cell.al));
        g->SetCellValue(i, 6, wxString::Format("%.3f", s.cell.be));
        g->SetCellValue(i, 7, wxString::Format("%.3f", s.cell.ga));
        g->SetCellValue(i, 8, wxString::Format("%.2f", s.vol));
        g->SetCellValue(i, 9, s.hasPlane ? wxString::Format("%.4f", s.qspec) : wxString(L"—"));
        g->SetCellValue(i, 10, wxString::Format("%.5f", s.rmsQxyz));
        g->SetCellValue(i, 11, wxString::Format("%.5f", s.rmsQxy));
        g->SetCellValue(i, 12, wxString::Format("%.5f", s.rmsQz));
    }
    g->ClearSelection();
}

void MainFrame::showSolution(int row) {
    if (row < 0 || row >= (int)m_solutions.size()) return;
    m_solRow = row;
    const gid::Solution& s = m_solutions[row];

    std::vector<gid::Peak> pred;
    pred.reserve(m_peaks.size());
    const int n = (int)m_peaks.size();
    if (m_peakGrid->GetNumberRows() != n) refreshPeakGrid();
    for (int i = 0; i < n; ++i) {
        gid::Hkl h{0, 0, 0};
        if (i < (int)s.hkl.size()) h = s.hkl[i];
        const gid::QPoint q = gid::predictQ(s, h.h, h.k, h.l);
        m_peakGrid->SetCellValue(i, 3, wxString::Format("%d", h.h));
        m_peakGrid->SetCellValue(i, 4, wxString::Format("%d", h.k));
        m_peakGrid->SetCellValue(i, 5, wxString::Format("%d", h.l));
        // Same convention as the solver: delta_qxy is always computed from |q_xy| (sign is mirror handedness only).
        const double obsQxy = std::fabs(m_peaks[i].qxy);
        m_peakGrid->SetCellValue(i, 6, wxString::Format("%.4f", q.qxy - obsQxy));
        m_peakGrid->SetCellValue(i, 7, wxString::Format("%.4f", q.qz - m_peaks[i].qz));
        pred.push_back(gid::Peak{q.qxy, q.qz});
    }
    m_plot->setPredicted(pred);

    if (s.hasPlane) {
        m_solInfo->SetLabel(wxString::Format(
            L"(u v w)=(%d %d %d)    a=%.4f  b=%.4f  c=%.4f    α=%.3f  β=%.3f  γ=%.3f    V=%.2f Å³    "
            L"q_spec=%.4f    RMSD: q_xyz=%.5f  q_xy=%.5f  q_z=%.5f",
            s.uvw.h, s.uvw.k, s.uvw.l, s.cell.a, s.cell.b, s.cell.c, s.cell.al, s.cell.be, s.cell.ga,
            s.vol, s.qspec, s.rmsQxyz, s.rmsQxy, s.rmsQz));
    } else {
        // P2-4: nu does not correspond to any integer contact plane (common in no-q_spec modes), so (u v w)/q_spec are meaningless
        m_solInfo->SetLabel(wxString::Format(
            L"(u v w)=— (no integer contact plane)    a=%.4f  b=%.4f  c=%.4f    α=%.3f  β=%.3f  γ=%.3f    "
            L"V=%.2f Å³    q_spec=—    RMSD: q_xyz=%.5f  q_xy=%.5f  q_z=%.5f",
            s.cell.a, s.cell.b, s.cell.c, s.cell.al, s.cell.be, s.cell.ga, s.vol, s.rmsQxyz,
            s.rmsQxy, s.rmsQz));
    }
    m_solPage->Layout();
}

void MainFrame::exportSolution() {
    if (m_solRow < 0 || m_solRow >= (int)m_solutions.size()) {
        wxMessageBox(L"Select a solution on the \"Solutions\" tab first", L"Export", wxICON_INFORMATION | wxOK, this);
        return;
    }
    const gid::Solution& s = m_solutions[m_solRow];
    if (!collectPeaks()) return;

    wxFileDialog dlg(this, L"Export Indexing Solution", "", "gidindex_solution.csv",
                     L"CSV files (*.csv)|*.csv|All files (*.*)|*.*",
                     wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if (dlg.ShowModal() != wxID_OK) return;

    const double par14[14] = {(double)s.uvw.h, (double)s.uvw.k,  (double)s.uvw.l, s.cell.a,
                              s.cell.b,        s.cell.c,         s.cell.al,       s.cell.be,
                              s.cell.ga,       s.vol,            s.rmsQxyz,       s.rmsQxy,
                              s.rmsQz,         s.dqspec};
    const int n = (int)m_peaks.size();
    Mat idx(n, 3, 0.0), raw(n, 2, 0.0);
    for (int i = 0; i < n; ++i) {
        gid::Hkl h{0, 0, 0};
        if (i < (int)s.hkl.size()) h = s.hkl[i];
        idx(i, 0) = h.h;
        idx(i, 1) = h.k;
        idx(i, 2) = h.l;
        raw(i, 0) = m_peaks[i].qxy;
        raw(i, 1) = m_peaks[i].qz;
    }
    std::string err;
    if (!exportSolutionCSV(std::string(dlg.GetPath().utf8_str()), par14, idx, raw, nullptr, err))
        wxMessageBox(wxString::FromUTF8(err.c_str()), L"Export Failed", wxICON_ERROR | wxOK, this);
    else
        SetStatusText(wxString::Format(L"Exported: %s", dlg.GetPath()), 0);
}

// ---------------------------------------------------------------- Event handlers
void MainFrame::onImport(wxCommandEvent&) {
    wxFileDialog dlg(this, L"Import Peaks", "", "",
                     L"Data files (*.csv;*.txt;*.dat;*.xlsx;*.xlsm)|*.csv;*.txt;*.dat;*.xlsx;*.xlsm|"
                     L"All files (*.*)|*.*",
                     wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (dlg.ShowModal() == wxID_OK) importFile(dlg.GetPath());
}

void MainFrame::onExport(wxCommandEvent&) { exportSolution(); }
void MainFrame::onRun(wxCommandEvent&) { startIndexing(); }
void MainFrame::onStop(wxCommandEvent&) { stopIndexing(); }
void MainFrame::onClear(wxCommandEvent&) { clearPeaks(); }
void MainFrame::onAppend(wxCommandEvent&) { appendPeakRow(); }
void MainFrame::onDelete(wxCommandEvent&) { deleteSelectedPeakRows(); }
void MainFrame::onAutoQSpec(wxCommandEvent&) { autoDetectQSpec(); }

void MainFrame::onSolSelect(wxGridEvent& e) {
    if (e.GetEventObject() != m_solGrid) return;
    const int r = e.GetRow();
    if (r >= 0 && r != m_solRow && r < (int)m_solutions.size()) showSolution(r);
}

void MainFrame::onUsage(wxCommandEvent&) {
    wxMessageBox(
        L"Steps:\n\n"
        "1) Import a peak table (two columns: q_xy, q_z; CSV / TXT / XLSX), or type values directly into the table;\n"
        "2) Choose a mode:\n"
        "     · q_spec (2-step) — requires q_spec; click \"Auto\" to get it from specular rod peaks at q_xy = 0;\n"
        "     · No q_spec (1-step) — no q_spec needed; solves the reciprocal metric and orientation jointly;\n"
        "     · No q_spec (2-step) — no q_spec needed; q_spec is the 4th unknown τ searched in step 1;\n"
        "     · Fixed (uvw) — set the out-of-plane contact plane manually; needs no prior cell parameters;\n"
        "     · Simple (001) — out-of-plane defaults to (001); uses an exact three-peak solution + torus Hough\n"
        "       search, requires no q_spec, and 3 non-rod peaks suffice;\n"
        "3) Adjust search ranges and tolerances on the \"Params\" tab, then click \"Run Indexing\";\n"
        "4) Solutions are listed on the \"Solutions\" tab; click a row to see each peak's (hkl) and residual in the peak table,"
        "and save them with \"Export\".\n\n"
        "Notes on simple indexing:\n"
        "     · The model q_xy² = h²A + k²B − 2hkC, q_z = hP + kQ + lR is an ordinary\n"
        "       reciprocal metric; output solutions fix uvw = (0,0,1) and q_spec = R;\n"
        "     · In-plane (A,B,C) comes from exact three-peak solutions + all-peak consensus, with a lattice-density gate;\n"
        "     · The out-of-plane period R is estimated from specular rod peaks near q_xy≈0 (p/n scan when the order is ambiguous);\n"
        "       tilt (P,Q) is searched by Hough transform over the torus [−R/2,R/2)²;\n"
        "     · Finally, nearest-neighbor reassignment + two alternating linear least-squares refinements are applied.\n\n"
        "Notes on the no-q_spec two-step mode:\n"
        "     · The Params tab lists the τ search interval/coarse budget, the rod soft prior, and the enumeration\n"
        "       ranges and tolerances of both steps; all are editable, and q_spec need not be entered;\n"
        "     · Rod peaks near q_xy≈0 are not fitted; they only supply τ candidates (a soft prior, not a measured q_spec);\n"
        "     · The output q_spec is derived from the solution itself; at least 4 non-rod peaks are required;\n"
        "     · Without rod peaks the discrimination weakens; consider narrowing the τ range and re-running.\n\n"
        "Tip: when solutions are not unique, prefer those with small RMSD(q_xyz) and a plausible volume.",
        L"Usage", wxOK | wxICON_INFORMATION, this);
}

void MainFrame::onAbout(wxCommandEvent&) {
    wxMessageBox(L"GIDIndexing 1.0.0\n\nGIWAXS / GIXD diffraction peak indexing\n"
                 L"· q_spec (2-step)\n· No q_spec (1-step)\n· No q_spec (2-step)\n"
                 L"· Fixed (uvw) global search\n· Simple (001) automatic search",
                 L"About GIDIndexing", wxOK | wxICON_INFORMATION, this);
}
