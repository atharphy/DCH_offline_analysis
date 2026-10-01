#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void PlotDeltaRates(const std::string& year = "2018") {
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".4f");
    const std::string outdir = "DD_flip_rates/plots/5_delta_rate";
    gSystem->mkdir(outdir.c_str(), true);

    TFile* f = TFile::Open(("DD_flip_rates/" + year + "/charge_flip_delta_rate_" + year + ".root").c_str());
    if (!f || f->IsZombie()) { std::cerr << "[error] cannot open charge_flip_delta_rate_" << year << ".root" << std::endl; return; }

    std::vector<std::pair<std::string,std::string>> hists = {
        {"rate_ele_data_raw_" + year, "data (raw, no bkg subtraction)"},
        {"rate_ele_data_bkgsub_" + year, "data (bkg-subtracted)"},
        {"rate_ele_mc_" + year, "MC (DY, intrinsic)"},
        {"delta_rate_ele_" + year, "delta = max(0, data_bkgsub - MC)"},
        {"ratio_nosub_data_mc_" + year, "data(raw) / MC"},
        {"ratio_withsub_data_mc_" + year, "data(bkg-subtracted) / MC"}
    };
    for (auto& hp : hists) {
        TH2D* h = (TH2D*)f->Get(hp.first.c_str());
        if (!h) { std::cerr << "[warn] missing " << hp.first << std::endl; continue; }
        TCanvas* c = new TCanvas(("c_" + hp.first).c_str(), "", 1400, 700);
        c->SetLogx();
        c->SetRightMargin(0.14);
        c->SetLeftMargin(0.09);
        h->SetTitle((hp.second + ";p_{T} [GeV];|#eta|").c_str());
        h->GetXaxis()->SetMoreLogLabels();
        h->GetXaxis()->SetNoExponent();
        h->Draw("COLZ TEXT");
        c->SaveAs((outdir + "/" + hp.first + ".png").c_str());
        std::cout << "[saved] " << outdir << "/" << hp.first << ".png" << std::endl;
    }
    f->Close();
}
