// Plots the top-pT reweighting SF (and its Up/Down bracket) vs generator
// top/antitop pT. Unlike the other PlotXxxSF.C scripts here, this factor
// isn't derived from data and stored in a weights file -- it's the fixed
// analytic CMS Top PAG formula from DCH_modules/TopPtReweighting.h, so the
// three curves are built directly from that header instead of read back
// from a ROOT file. The SF itself doesn't depend on year, but one plot is
// still written per year (identical content) to match this directory's own
// per-year output convention.
//
// Run: root -l -b -q systematics_plots/PlotTopPtWeight.C+
#include <TH1D.h>
#include <iostream>
#include <string>
#include <vector>

#include "../include/MyBranch.C"
#include "SFPlotStyle.h"
#include "../DCH_modules/TopPtReweighting.h"

namespace {
const std::string kOutDir = "systematics_plots/output";
const std::vector<std::string> kYears = {"2016preVFP", "2016postVFP", "2017", "2018"};
}

void PlotTopPtWeight() {
    gStyle->SetOptStat(0);

    const int nBins = 60;
    const double ptMin = 0.0, ptMax = 600.0;
    TH1D* hNom = new TH1D("hTopPtNom", "", nBins, ptMin, ptMax);
    TH1D* hUp = new TH1D("hTopPtUp", "", nBins, ptMin, ptMax);
    TH1D* hDown = new TH1D("hTopPtDown", "", nBins, ptMin, ptMax);
    for (int b = 1; b <= nBins; ++b) {
        const double pt = hNom->GetXaxis()->GetBinCenter(b);
        const double sf = TopPtSF::sf(pt);
        hNom->SetBinContent(b, sf);
        hUp->SetBinContent(b, sf * sf);      // Top PAG bracket: up = SF^2
        hDown->SetBinContent(b, 1.0);        // Top PAG bracket: down = no reweighting
    }
    zeroBinErrors(hNom);

    for (const auto& year : kYears) {
        std::cout << "=== " << year << " ===" << std::endl;
        TH1D* hNomY = (TH1D*)hNom->Clone(("hTopPtNom_" + year).c_str());
        TH1D* hUpY = (TH1D*)hUp->Clone(("hTopPtUp_" + year).c_str());
        TH1D* hDownY = (TH1D*)hDown->Clone(("hTopPtDown_" + year).c_str());
        drawSF1D(hNomY, hUpY, hDownY, "topPt_" + year, "Generator top/antitop p_{T} [GeV]", "Top-p_{T} SF", year,
                  kOutDir + "/" + year + "/topPt.png", 0.6, 1.3);
        delete hNomY; delete hUpY; delete hDownY;
    }
    delete hNom; delete hUp; delete hDown;
    std::cout << "Done." << std::endl;
}
