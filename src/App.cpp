#include <wx/wx.h>

#include "ui/MainFrame.h"

class GIDIndexingApp : public wxApp {
public:
    bool OnInit() override {
        wxInitAllImageHandlers();
        auto* frame = new MainFrame();
        frame->SetIcon(wxICON(APPICON));
        frame->Show(true);
        SetTopWindow(frame);
        return true;
    }
};

wxIMPLEMENT_APP(GIDIndexingApp);
