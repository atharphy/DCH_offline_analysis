// Plots the DY ZpT reweighting factor (and its derivation-side statistical
// uncertainty) vs generator Z pT, per year, directly from the already-
// derived ZPtWeights_<year>.root files (same file DCH_tight.C/DCH_tauFR.C
// read at analysis time via ZPtReweight.h's loadZPtWeights).
//
// Run: root -l -b -q systematics_plots/PlotZPtWeight.C+
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>

#include "SFPlotStyle.h"

namespace {
const std::string kOutDir = "systematics_plots/output";
const std::vector<std::string> kYears = {"2016preVFP", "2016postVFP", "2017", "2018"};
}

void PlotZPtWeight() {
    gStyle->SetOptStat(0);
    for (const auto& year : kYears) {
        std::cout << "=== " << year << " ===" << std::endl;
        const std::string fileName = "../Dependencies/zpt_weights/" + year + "/ZPtWeights_" + year + ".root";
        TFile* f = TFile::Open(fileName.c_str(), "READ");
        if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << fileName << std::endl; continue; }
        TH1D* hNom = (TH1D*)f->Get("h_zpt_weight_ee");
        if (!hNom) { std::cerr << "ERROR: missing h_zpt_weight_ee in " << fileName << std::endl; f->Close(); continue; }
        hNom = (TH1D*)hNom->Clone(("hZpt_" + year).c_str());
        hNom->SetDirectory(nullptr);
        f->Close();

        TH1D* hUp = (TH1D*)hNom->Clone(("hZptUp_" + year).c_str());
        TH1D* hDown = (TH1D*)hNom->Clone(("hZptDown_" + year).c_str());
        hUp->SetDirectory(nullptr); hDown->SetDirectory(nullptr);
        for (int b = 1; b <= hNom->GetNbinsX(); ++b) {
            const double v = hNom->GetBinContent(b);
            const double e = hNom->GetBinError(b);
            hUp->SetBinContent(b, v + e);
            hDown->SetBinContent(b, std::max(0.0, v - e));
        }

        drawSF1D(hNom, hUp, hDown, "zpt_" + year, "Generator Z p_{T} [GeV]", "ZpT reweighting factor", year,
                  kOutDir + "/" + year + "/zpt.png", 0.5, 1.5);
        delete hNom; delete hUp; delete hDown;
    }
    std::cout << "Done." << std::endl;
}
