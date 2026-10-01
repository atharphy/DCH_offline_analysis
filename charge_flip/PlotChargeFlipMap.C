#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>
void PlotChargeFlipMap() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("correction_plots", kTRUE);
    const std::vector<std::string> tags = {"2016preVFP", "2016postVFP", "2017", "2018", "Run2"};
    for (const auto& tag : tags) {
        gSystem->mkdir(("correction_plots/" + tag).c_str(), kTRUE);
        TFile* f = TFile::Open(("DD_flip_rates/" + tag + "/charge_flip_delta_rate_" + tag + ".root").c_str());
        if (!f || f->IsZombie()) { std::cout << "  [warn] missing " << tag << std::endl; continue; }
        TH2D* h = (TH2D*)f->Get(("rate_ele_data_bkgsub_" + tag).c_str());
        if (!h) { std::cout << "  [warn] no hist for " << tag << std::endl; f->Close(); continue; }
        h = (TH2D*)h->Clone(); h->SetDirectory(nullptr);
        h->SetTitle(";p_{T} [GeV];|#eta|");
        TCanvas c("c", "", 900, 700);
        c.SetLeftMargin(0.12); c.SetRightMargin(0.16); c.SetBottomMargin(0.12); c.SetTopMargin(0.10);
        h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
        h->GetZaxis()->SetTitle("charge mis-ID rate");
        h->Draw("colz");
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.04); lat.SetTextFont(62);
        lat.DrawLatex(0.12, 0.94, ("Electron charge mis-ID rate, " + tag + " data").c_str());
        c.SaveAs(("correction_plots/" + tag + "/chargeflip_map.png").c_str());
        f->Close();
    }
}
