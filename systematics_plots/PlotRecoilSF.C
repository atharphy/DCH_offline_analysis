// Plots the CMS-HTT recoil correction (MEtSys) response curves and their
// systematic uncertainties, per year.
//
// Unlike the lepton/tau SFs, this payload isn't natively a function of
// (pt, eta): it stores (1) response curves vs generator boson pT, one per
// jet-multiplicity bin, and (2) a small discrete uncertainty table indexed
// by (Response/Resolution) x (njet bin), read directly off the same ROOT
// payload MEtSysWrapper.h loads (bypassing the MEtSys class itself, since
// its response histograms and sysUnc table are private members with no
// public accessor -- this reads the identical named histograms directly).
// Only the "V" (vector boson, i.e. Z/W+jets) category is plotted: it's the
// only ProcessType this framework ever uses (see MEtSysWrapper.h's
// hardcoded MEtSys::ProcessType::BOSON).
//
// Run: root -l -b -q systematics_plots/PlotRecoilSF.C+
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>

#include "SFPlotStyle.h"

namespace {
const std::string kOutDir = "systematics_plots/output";
const std::vector<std::string> kYears = {"2016preVFP", "2016postVFP", "2017", "2018"};
const std::vector<std::string> kJetBins = {"NJet0", "NJet1", "NJetGe2"};

std::string payloadFor(const std::string& year) {
    if (year == "2016preVFP" || year == "2016postVFP") return "../HTT-utilities/RecoilCorrections/data/PFMEtSys_2016.root";
    if (year == "2017") return "../HTT-utilities/RecoilCorrections/data/PFMEtSys_2017.root";
    return "../HTT-utilities/RecoilCorrections/data/PFMEtSys_2018.root";
}
}  // namespace

void PlotRecoilSF() {
    gStyle->SetOptStat(0);
    for (const auto& year : kYears) {
        std::cout << "=== " << year << " ===" << std::endl;
        const std::string path = payloadFor(year);
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << path << std::endl; continue; }

        // --- response curves vs genVpt, one per njet bin ---
        std::vector<TH1D*> curves;
        for (auto& jb : kJetBins) {
            TH1D* h = (TH1D*)f->Get(("V_" + jb).c_str());
            if (!h) { std::cerr << "  missing V_" << jb << std::endl; continue; }
            h = (TH1D*)h->Clone(("hResp_" + jb + "_" + year).c_str());
            h->SetDirectory(nullptr);
            curves.push_back(h);
        }
        if (curves.size() == 3) {
            TCanvas c(("c_recoilResponse_" + year).c_str(), "", 1000, 750);
            c.SetLeftMargin(0.14); c.SetRightMargin(0.05); c.SetBottomMargin(0.14); c.SetTopMargin(0.12); c.SetTicks(1, 1);
            curves[0]->SetTitle("");
            curves[0]->GetXaxis()->SetTitle("Generator boson p_{T} [GeV]");
            curves[0]->GetYaxis()->SetTitle("Mean hadronic recoil response");
            curves[0]->GetXaxis()->SetTitleSize(0.045);
            curves[0]->GetYaxis()->SetTitleSize(0.045);
            curves[0]->SetLineColor(kBlack); curves[0]->SetLineWidth(2);
            curves[1]->SetLineColor(kRed + 1); curves[1]->SetLineWidth(2);
            curves[2]->SetLineColor(kAzure + 2); curves[2]->SetLineWidth(2);
            curves[0]->Draw("HIST");
            curves[1]->Draw("HIST SAME");
            curves[2]->Draw("HIST SAME");
            TLegend leg(0.55, 0.18, 0.90, 0.35);
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
            leg.AddEntry(curves[0], "N_{jet} = 0", "l");
            leg.AddEntry(curves[1], "N_{jet} = 1", "l");
            leg.AddEntry(curves[2], "N_{jet} #geq 2", "l");
            leg.Draw();
            applyCMSLabel(0.14);
            applyYearLabel(year);
            gSystem->mkdir((kOutDir + "/" + year).c_str(), kTRUE);
            c.SaveAs((kOutDir + "/" + year + "/recoilResponseCurve.png").c_str());
        }
        for (auto* h : curves) delete h;

        // --- uncertainty table: (Response/Resolution) x (njet bin), in % ---
        TH2D* hSyst = (TH2D*)f->Get("V_syst");
        if (hSyst) {
            TH2D* hUnc = new TH2D(("hUnc_recoil_" + year).c_str(), "", 2, 0, 2, 3, 0, 3);
            hUnc->SetDirectory(nullptr);
            hUnc->GetXaxis()->SetBinLabel(1, "Response");
            hUnc->GetXaxis()->SetBinLabel(2, "Resolution");
            hUnc->GetYaxis()->SetBinLabel(1, "N_{jet}=0");
            hUnc->GetYaxis()->SetBinLabel(2, "N_{jet}=1");
            hUnc->GetYaxis()->SetBinLabel(3, "N_{jet}#geq2");
            for (int xb = 0; xb < 2; ++xb)
                for (int yb = 0; yb < 3; ++yb)
                    hUnc->SetBinContent(xb + 1, yb + 1, hSyst->GetBinContent(xb + 1, yb + 1) * 100.0);

            TCanvas c(("c_recoilUnc_" + year).c_str(), "", 800, 700);
            c.SetLeftMargin(0.18); c.SetRightMargin(0.18); c.SetBottomMargin(0.14); c.SetTopMargin(0.12); c.SetTicks(1, 1);
            hUnc->SetTitle("");
            hUnc->GetZaxis()->SetTitle("Uncertainty [%]");
            hUnc->SetMinimum(0.0);
            hUnc->SetMaximum(5.0);
            hUnc->SetMarkerSize(1.6);
            gStyle->SetPaintTextFormat(".3f");
            hUnc->Draw("COLZ TEXT");
            applyCMSLabel(0.18);
            applyYearLabel(year);
            c.SaveAs((kOutDir + "/" + year + "/recoilUncertainty.png").c_str());
            delete hUnc;
        } else {
            std::cerr << "  missing V_syst" << std::endl;
        }

        f->Close();
    }
    std::cout << "Done." << std::endl;
}
