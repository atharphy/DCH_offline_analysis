#include "EtauCommon.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

std::unique_ptr<TH1D> rebinClone(TH1D* source, const std::string& newName) {
    const int nBins = static_cast<int>(ETAU_PT_EDGES.size()) - 1;
    std::unique_ptr<TH1D> result(dynamic_cast<TH1D*>(source->Rebin(nBins, newName.c_str(), ETAU_PT_EDGES.data())));
    result->SetDirectory(nullptr);
    return result;
}

}

void RebinAndPlotFakeRates() {
    gStyle->SetOptStat(0);
    gStyle->SetGridColor(kGray + 2);
    gStyle->SetGridStyle(2);
    gStyle->SetGridWidth(1);

    const std::vector<std::string> years = {"2016preVFP", "2016postVFP", "2017", "2018"};
    const std::vector<std::string> regions = {"barrel", "endcap"};

    std::map<std::string, std::unique_ptr<TH1D>> num, den, fr;

    for (const std::string& year : years) {
        const std::string inName = "cut_results/quick_check_numden_" + year + ".root";
        TFile fin(inName.c_str(), "READ");
        if (fin.IsZombie()) { std::cerr << "ERROR: cannot open " << inName << std::endl; gSystem->Exit(1); }

        for (const std::string& region : regions) {
            const std::string key = year + "_" + region;
            TH1D* numFine = (TH1D*)fin.Get(("Num_" + key).c_str());
            TH1D* denFine = (TH1D*)fin.Get(("Den_" + key).c_str());
            if (!numFine || !denFine) { std::cerr << "ERROR: missing histograms for " << key << std::endl; gSystem->Exit(2); }

            num[key] = rebinClone(numFine, "Num_" + key + "_rebin");
            den[key] = rebinClone(denFine, "Den_" + key + "_rebin");

            fr[key].reset((TH1D*)num[key]->Clone(("ETau_FR_" + key).c_str()));
            fr[key]->Divide(num[key].get(), den[key].get(), 1.0, 1.0, "B");
            fr[key]->SetDirectory(nullptr);
        }
    }

    const std::vector<std::string> allYears = {"2016preVFP", "2016postVFP", "2017", "2018", "Run2"};
    for (const std::string& region : regions) {
        const std::string run2Key = "Run2_" + region;
        std::unique_ptr<TH1D> numRun2, denRun2;
        for (const std::string& year : years) {
            const std::string key = year + "_" + region;
            if (!numRun2) {
                numRun2.reset((TH1D*)num[key]->Clone(("Num_" + run2Key).c_str()));
                denRun2.reset((TH1D*)den[key]->Clone(("Den_" + run2Key).c_str()));
            } else {
                numRun2->Add(num[key].get());
                denRun2->Add(den[key].get());
            }
        }
        fr[run2Key].reset((TH1D*)numRun2->Clone(("ETau_FR_" + run2Key).c_str()));
        fr[run2Key]->Divide(numRun2.get(), denRun2.get(), 1.0, 1.0, "B");
        num[run2Key] = std::move(numRun2);
        den[run2Key] = std::move(denRun2);
    }

    gSystem->mkdir("cut_results", kTRUE);
    const std::string outName = "cut_results/ETau_FakeRates_Final.root";
    TFile fout(outName.c_str(), "RECREATE");
    for (const std::string& year : allYears)
        for (const std::string& region : regions)
            fr[year + "_" + region]->Write();
    fout.Close();
    std::cout << "output: " << outName << std::endl;

    std::map<std::string, int> colors = {
        {"2016preVFP", kRed + 1}, {"2016postVFP", kAzure + 1}, {"2017", kGreen + 2},
        {"2018", kMagenta + 1}, {"Run2", kBlack}
    };

    TCanvas c("c_compare", "", 1800, 800);
    c.Divide(2, 1);

    for (int i = 0; i < 2; ++i) {
        const std::string& region = regions[i];
        c.cd(i + 1);
        gPad->SetGrid();
        gPad->SetTicks(1, 1);
        gPad->SetLogx();
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.10);

        TLegend* leg = new TLegend(0.55, 0.65, 0.94, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        bool first = true;
        for (const std::string& year : allYears) {
            TH1D* h = fr[year + "_" + region].get();
            h->SetLineColor(colors[year]);
            h->SetLineWidth(year == "Run2" ? 3 : 2);
            h->SetMarkerColor(colors[year]);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(year == "Run2" ? 0 : 0.8);
            h->SetTitle(";p_{T}^{#tau} [GeV];Electron-to-tau fake rate");
            h->GetYaxis()->SetRangeUser(0.0, 1.0);
            h->GetXaxis()->SetMoreLogLabels();
            h->GetXaxis()->SetNoExponent();
            h->Draw(first ? "E1" : "E1 SAME");
            leg->AddEntry(h, year.c_str(), "lep");
            first = false;
        }
        leg->Draw();

        TLatex label;
        label.SetNDC();
        label.SetTextFont(61);
        label.SetTextSize(0.05);
        label.DrawLatex(0.16, 0.93, region == "barrel" ? "Barrel (|#eta|<1.5)" : "Endcap (1.5<|#eta|<2.5)");
    }

    c.SaveAs("cut_results/ETau_FakeRates_Comparison.png");
    c.SaveAs("cut_results/ETau_FakeRates_Comparison.pdf");
    std::cout << "plot: cut_results/ETau_FakeRates_Comparison.png" << std::endl;
}
