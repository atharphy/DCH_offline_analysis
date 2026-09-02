// Plots the Rochester muon momentum correction (kScaleMC) and its
// uncertainty (kScaleMCerror) as 2D (pt, eta) maps, per year.
//
// Run: root -l -b -q systematics_plots/PlotRoccorSF.C+
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>

#include "../roccor/RoccoR.cc"
#include "SFPlotStyle.h"

namespace {
const std::string kOutDir = "systematics_plots/output";
const std::vector<std::string> kYears = {"2016preVFP", "2016postVFP", "2017", "2018"};

const std::vector<double> kPtEdges = {20, 30, 40, 50, 60, 80, 100, 150, 200};
const std::vector<double> kEtaEdges = {-2.4, -1.6, -0.8, 0.0, 0.8, 1.6, 2.4};

std::string payloadFor(const std::string& year) {
    if (year == "2016preVFP") return "RoccoR2016aUL.txt";
    if (year == "2016postVFP") return "RoccoR2016bUL.txt";
    if (year == "2017") return "RoccoR2017UL.txt";
    return "RoccoR2018UL.txt";
}
}  // namespace

void PlotRoccorSF() {
    gStyle->SetOptStat(0);
    for (const auto& year : kYears) {
        std::cout << "=== " << year << " ===" << std::endl;
        const std::string path = "../Offline_framework/offline/roccor/" + payloadFor(year);
        RoccoR rc(path);

        TH2D* hSF = new TH2D(("hSF_roccor_" + year).c_str(), "", kPtEdges.size() - 1, kPtEdges.data(), kEtaEdges.size() - 1, kEtaEdges.data());
        hSF->SetDirectory(nullptr);
        TH2D* hUnc = new TH2D(("hUnc_roccor_" + year).c_str(), "", kPtEdges.size() - 1, kPtEdges.data(), kEtaEdges.size() - 1, kEtaEdges.data());
        hUnc->SetDirectory(nullptr);

        for (int bx = 1; bx <= hSF->GetNbinsX(); ++bx) {
            const double pt = hSF->GetXaxis()->GetBinCenter(bx);
            for (int by = 1; by <= hSF->GetNbinsY(); ++by) {
                const double eta = hSF->GetYaxis()->GetBinCenter(by);
                const double sf = rc.kScaleMC(1, pt, eta, 0.0);
                const double err = rc.kScaleMCerror(1, pt, eta, 0.0);
                hSF->SetBinContent(bx, by, sf);
                hUnc->SetBinContent(bx, by, sf > 0.0 ? err / sf * 100.0 : 0.0);
            }
        }
        drawSFAndUncertainty2D(hSF, hUnc, "roccor_" + year, "p_{T} [GeV]", "#eta", year,
                                kOutDir + "/" + year + "/roccor.png", 0.98, 1.02, 0.15, true);
        delete hSF; delete hUnc;
    }
    std::cout << "Done." << std::endl;
}
