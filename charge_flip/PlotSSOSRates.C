#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

void PlotSSOSRates() {
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".4f");

    const std::string dirRaw = "DD_flip_rates/plots/1_raw_pair_counts";
    const std::string dirRatio = "DD_flip_rates/plots/2_ratio_matrix";
    const std::string dirSolve = "DD_flip_rates/plots/3_solve_fit";
    const std::string dirRate = "DD_flip_rates/plots/4_final_rate_maps";
    gSystem->mkdir(dirRaw.c_str(), true);
    gSystem->mkdir(dirRatio.c_str(), true);
    gSystem->mkdir(dirSolve.c_str(), true);
    gSystem->mkdir(dirRate.c_str(), true);

    auto styleAxis = [](TH2D* h) {
        h->SetStats(0);
        h->GetXaxis()->SetMoreLogLabels();
        h->GetXaxis()->SetNoExponent();
        h->GetXaxis()->SetNdivisions(515);
        h->GetXaxis()->SetLabelSize(0.03);
    };
    auto drawAndSave = [&](TH2D* h, const std::string& title, const std::string& canvasName, const std::string& outPath) {
        TCanvas* c = new TCanvas(canvasName.c_str(), "", 1600, 700);
        c->SetLogx();
        c->SetRightMargin(0.14);
        c->SetLeftMargin(0.08);
        c->SetBottomMargin(0.12);
        h->SetTitle(title.c_str());
        styleAxis(h);
        h->Draw("COLZ");

        TLatex latex;
        latex.SetTextAlign(22);
        latex.SetTextAngle(90);
        latex.SetTextSize(0.028);
        latex.SetTextFont(42);
        for (int bx = 1; bx <= h->GetNbinsX(); bx++) {
            for (int by = 1; by <= h->GetNbinsY(); by++) {
                const double val = h->GetBinContent(bx, by);
                if (val == 0.0) continue;
                char buf[32];
                snprintf(buf, sizeof(buf), "%.4f", val);
                latex.DrawLatex(h->GetXaxis()->GetBinCenter(bx), h->GetYaxis()->GetBinCenter(by), buf);
            }
        }

        c->SaveAs(outPath.c_str());
    };

    auto drawMatrix = [&](TH2D* h, const std::string& title, const std::string& canvasName,
                           const std::string& outPath, bool logZ) {
        TCanvas* c = new TCanvas(canvasName.c_str(), "", 1400, 1200);
        if (logZ) c->SetLogz();
        c->SetRightMargin(0.13);
        c->SetLeftMargin(0.10);
        c->SetBottomMargin(0.10);
        h->SetStats(0);
        h->SetTitle(title.c_str());
        h->Draw("COLZ");
        c->SaveAs(outPath.c_str());
    };

    std::vector<std::string> years = {"2016preVFP", "2016postVFP", "2017", "2018"};
    for (auto& year : years) {
        TFile* fd = TFile::Open(("DD_flip_rates/" + year + "/ssos_ratio_data_" + year + ".root").c_str());
        TFile* fm = TFile::Open(("DD_flip_rates/" + year + "/ssos_ratio_mc_" + year + ".root").c_str());
        if (!fd || fd->IsZombie() || !fm || fm->IsZombie()) { std::cerr << "[error] cannot open rate files for " << year << std::endl; continue; }
        TH2D* rateData = (TH2D*)fd->Get(("rate_ele_data_" + year).c_str());
        TH2D* rateMC = (TH2D*)fm->Get(("rate_ele_mc_" + year).c_str());
        if (!rateData || !rateMC) { std::cerr << "[error] missing rate histograms for " << year << std::endl; fd->Close(); fm->Close(); continue; }

        drawAndSave(rateData, year + " electron charge mis-ID rate (data);p_{T} [GeV];|#eta|",
                    "c_data_" + year, dirRate + "/ssos_rate_data_" + year + ".png");
        drawAndSave(rateMC, year + " electron charge mis-ID rate (DY MC);p_{T} [GeV];|#eta|",
                    "c_mc_" + year, dirRate + "/ssos_rate_mc_" + year + ".png");

        TH2D* ratio = (TH2D*)rateData->Clone(("ratio_" + year).c_str());
        ratio->Divide(rateMC);
        drawAndSave(ratio, year + " charge mis-ID rate ratio (data / MC);p_{T} [GeV];|#eta|",
                    "c_ratio_" + year, dirRate + "/ssos_rate_ratio_" + year + ".png");

        for (auto& cat : {std::string("data"), std::string("mc")}) {
            TFile* f = (cat == "data") ? fd : fm;
            TH2D* nOS = (TH2D*)f->Get("n_os_ele");
            TH2D* nSS = (TH2D*)f->Get("n_ss_ele");
            if (!nOS || !nSS) { std::cerr << "[warn] missing n_os_ele/n_ss_ele for " << year << "/" << cat << std::endl; continue; }

            drawMatrix(nOS, year + " N_{OS}(cell_{i},cell_{j}) [" + cat + "];cell index i;cell index j",
                       "c_nos_" + cat + "_" + year, dirRaw + "/ssos_matrix_nos_" + cat + "_" + year + ".png", true);
            drawMatrix(nSS, year + " N_{SS}(cell_{i},cell_{j}) [" + cat + "];cell index i;cell index j",
                       "c_nss_" + cat + "_" + year, dirRaw + "/ssos_matrix_nss_" + cat + "_" + year + ".png", true);

            TH2D* rij = (TH2D*)nOS->Clone(("rij_" + cat + "_" + year).c_str());
            rij->Reset();
            for (int bx = 1; bx <= nOS->GetNbinsX(); bx++) {
                for (int by = 1; by <= nOS->GetNbinsY(); by++) {
                    const double os = nOS->GetBinContent(bx, by);
                    const double ss = nSS->GetBinContent(bx, by);
                    const double n = os + ss;
                    if (n > 0) rij->SetBinContent(bx, by, ss / n);
                }
            }
            drawMatrix(rij, year + " r_{ij} = N_{SS}/(N_{SS}+N_{OS}) [" + cat + "];cell index i;cell index j",
                       "c_rij_" + cat + "_" + year, dirRatio + "/ssos_matrix_rij_" + cat + "_" + year + ".png", false);

            TH2D* rateMap = (cat == "data") ? rateData : rateMC;
            const int nPt = rateMap->GetNbinsX();
            std::vector<double> pByCell(nOS->GetNbinsX() + 1, 0.0);
            for (int ptBin = 1; ptBin <= nPt; ptBin++)
                for (int etaBin = 1; etaBin <= 2; etaBin++)
                    pByCell[(ptBin - 1) * 2 + etaBin] = rateMap->GetBinContent(ptBin, etaBin);

            TH2D* predicted = (TH2D*)nOS->Clone(("predicted_" + cat + "_" + year).c_str());
            TH2D* residual = (TH2D*)nOS->Clone(("residual_" + cat + "_" + year).c_str());
            predicted->Reset();
            residual->Reset();
            for (int i = 1; i <= nOS->GetNbinsX(); i++) {
                for (int j = i; j <= nOS->GetNbinsY(); j++) {
                    const double measured = rij->GetBinContent(i, j);
                    if (nOS->GetBinContent(i, j) + nSS->GetBinContent(i, j) <= 0) continue;
                    const double pred = (i == j) ? 2.0 * pByCell[i] : (pByCell[i] + pByCell[j]);
                    predicted->SetBinContent(i, j, pred);
                    predicted->SetBinContent(j, i, pred);
                    residual->SetBinContent(i, j, measured - pred);
                    residual->SetBinContent(j, i, measured - pred);
                }
            }
            drawMatrix(predicted, year + " predicted r_{ij} = p_{i}+p_{j} from solved rates [" + cat + "];cell index i;cell index j",
                       "c_pred_" + cat + "_" + year, dirSolve + "/ssos_matrix_predicted_" + cat + "_" + year + ".png", false);
            drawMatrix(residual, year + " residual: measured r_{ij} - predicted [" + cat + "];cell index i;cell index j",
                       "c_resid_" + cat + "_" + year, dirSolve + "/ssos_matrix_residual_" + cat + "_" + year + ".png", false);
        }

        fd->Close();
        fm->Close();
    }
}
