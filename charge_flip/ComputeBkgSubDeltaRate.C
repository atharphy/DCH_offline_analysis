#include <TFile.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {
const std::vector<double> PT_EDGES = {
    20., 22., 24., 26., 28., 30., 32., 35., 40.,
    45., 50., 60., 70., 85., 100., 125., 150., 200.
};
const double ETA_EDGES[3] = {0.0, 1.5, 2.5};
const double MIN_PAIRS = 20.0;

TH2D* solveRate(TH2D* nOS, TH2D* nSS, const std::string& name, TH2D** errOut = nullptr) {
    const int nCells = nOS->GetNbinsX();
    std::vector<std::vector<double>> rows;
    std::vector<double> rhs, weights;
    for (int i = 1; i <= nCells; i++) {
        for (int j = i; j <= nCells; j++) {
            const double os = nOS->GetBinContent(i, j);
            const double ss = nSS->GetBinContent(i, j);
            const double nPairs = os + ss;
            if (nPairs < MIN_PAIRS) continue;
            const double r = ss / nPairs;
            if (r <= 0.0 || r >= 1.0) continue;
            const double variance = r * (1.0 - r) / nPairs;
            if (variance <= 0.0) continue;

            std::vector<double> row(nCells, 0.0);
            if (i == j) row[i - 1] = 2.0;
            else { row[i - 1] = 1.0; row[j - 1] = 1.0; }
            rows.push_back(row);
            rhs.push_back(r);
            weights.push_back(1.0 / variance);
        }
    }
    const int nEq = static_cast<int>(rows.size());
    std::cout << "  " << name << ": built " << nEq << " equations for " << nCells << " cells" << std::endl;

    TMatrixD A(nEq, nCells);
    TVectorD b(nEq);
    TVectorD w(nEq);
    for (int e = 0; e < nEq; e++) {
        for (int c = 0; c < nCells; c++) A(e, c) = rows[e][c];
        b(e) = rhs[e];
        w(e) = weights[e];
    }
    TMatrixD At(TMatrixD::kTransposed, A);
    TMatrixD AtW(At);
    for (int c = 0; c < nCells; c++)
        for (int e = 0; e < nEq; e++) AtW(c, e) *= w(e);
    TMatrixD AtWA = AtW * A;
    TVectorD AtWb = AtW * b;

    TDecompSVD svd(AtWA);
    Bool_t ok = kFALSE;
    TVectorD p = svd.Solve(AtWb, ok);
    if (!ok) std::cerr << "  [warn] SVD solve for " << name << " reported a non-exact solution" << std::endl;

    TMatrixD cov(AtWA);
    if (svd.Invert(cov)) {
        if (errOut) {
            const int nPt = static_cast<int>(PT_EDGES.size()) - 1;
            TH2D* err = new TH2D((name + "_err").c_str(), "", nPt, PT_EDGES.data(), 2, ETA_EDGES);
            err->SetDirectory(nullptr);
            for (int ptBin = 1; ptBin <= nPt; ptBin++) {
                for (int etaBin = 1; etaBin <= 2; etaBin++) {
                    const int cell = (ptBin - 1) * 2 + etaBin;
                    const double variance = cov(cell - 1, cell - 1);
                    err->SetBinContent(ptBin, etaBin, variance > 0.0 ? std::sqrt(variance) : 0.0);
                }
            }
            *errOut = err;
        }
    } else {
        std::cerr << "  [warn] could not invert (A^T W A) for " << name << " uncertainty" << std::endl;
        if (errOut) *errOut = nullptr;
    }

    const int nPt = static_cast<int>(PT_EDGES.size()) - 1;
    TH2D* rate = new TH2D(name.c_str(), "", nPt, PT_EDGES.data(), 2, ETA_EDGES);
    rate->SetDirectory(nullptr);
    for (int ptBin = 1; ptBin <= nPt; ptBin++) {
        for (int etaBin = 1; etaBin <= 2; etaBin++) {
            const int cell = (ptBin - 1) * 2 + etaBin;
            rate->SetBinContent(ptBin, etaBin, std::max(0.0, p(cell - 1)));
        }
    }
    return rate;
}
}

void ComputeBkgSubDeltaRate(const std::string& year = "2018") {
    const std::string dir = "DD_flip_rates/" + year + "/";

    TFile* fData = TFile::Open((dir + "ssos_ratio_data_" + year + ".root").c_str());
    TFile* fBkg = TFile::Open((dir + "ssos_ratio_bkg_" + year + ".root").c_str());
    TFile* fMc = TFile::Open((dir + "ssos_ratio_mc_" + year + ".root").c_str());
    if (!fData || fData->IsZombie() || !fBkg || fBkg->IsZombie() || !fMc || fMc->IsZombie()) {
        std::cerr << "[error] missing one of ssos_ratio_{data,bkg,mc}_" << year << ".root -- run MergeSSOSChunks.C for all three categories first" << std::endl;
        return;
    }

    TH2D* nOSData = (TH2D*)fData->Get("n_os_ele");
    TH2D* nSSData = (TH2D*)fData->Get("n_ss_ele");
    TH2D* nOSBkg = (TH2D*)fBkg->Get("n_os_ele");
    TH2D* nSSBkg = (TH2D*)fBkg->Get("n_ss_ele");
    TH2D* nOSMc = (TH2D*)fMc->Get("n_os_ele");
    TH2D* nSSMc = (TH2D*)fMc->Get("n_ss_ele");
    TH2D* rateDataRaw = (TH2D*)fData->Get(("rate_ele_data_" + year).c_str());
    if (!nOSData || !nSSData || !nOSBkg || !nSSBkg || !nOSMc || !nSSMc || !rateDataRaw) {
        std::cerr << "[error] missing expected histograms in one of the input files" << std::endl;
        return;
    }
    rateDataRaw = (TH2D*)rateDataRaw->Clone(("rate_ele_data_raw_" + year).c_str());
    rateDataRaw->SetDirectory(nullptr);

    std::cout << "Solving MC (intrinsic) rate..." << std::endl;
    TH2D* mcErr = nullptr;
    TH2D* rateMc = solveRate(nOSMc, nSSMc, "rate_ele_mc_" + year, &mcErr);

    TH2D* nOSSub = (TH2D*)nOSData->Clone("n_os_ele_bkgsub");
    TH2D* nSSSub = (TH2D*)nSSData->Clone("n_ss_ele_bkgsub");
    nOSSub->SetDirectory(nullptr);
    nSSSub->SetDirectory(nullptr);
    nOSSub->Add(nOSBkg, -1.0);
    nSSSub->Add(nSSBkg, -1.0);
    for (int bx = 1; bx <= nOSSub->GetNbinsX(); bx++) {
        for (int by = 1; by <= nOSSub->GetNbinsY(); by++) {
            if (nOSSub->GetBinContent(bx, by) < 0) nOSSub->SetBinContent(bx, by, 0);
            if (nSSSub->GetBinContent(bx, by) < 0) nSSSub->SetBinContent(bx, by, 0);
        }
    }

    std::cout << "Solving background-subtracted data rate..." << std::endl;
    TH2D* dataErr = nullptr;
    TH2D* rateDataBkgsub = solveRate(nOSSub, nSSSub, "rate_ele_data_bkgsub_" + year, &dataErr);

    TH2D* deltaRate = (TH2D*)rateDataBkgsub->Clone(("delta_rate_ele_" + year).c_str());
    deltaRate->SetDirectory(nullptr);
    deltaRate->Add(rateMc, -1.0);
    for (int bx = 1; bx <= deltaRate->GetNbinsX(); bx++)
        for (int by = 1; by <= deltaRate->GetNbinsY(); by++)
            if (deltaRate->GetBinContent(bx, by) < 0) deltaRate->SetBinContent(bx, by, 0);

    TH2D* deltaErr = nullptr;
    if (dataErr && mcErr) {
        deltaErr = (TH2D*)dataErr->Clone(("delta_rate_ele_" + year + "_err").c_str());
        deltaErr->SetDirectory(nullptr);
        for (int bx = 1; bx <= deltaErr->GetNbinsX(); bx++) {
            for (int by = 1; by <= deltaErr->GetNbinsY(); by++) {
                const double sData = dataErr->GetBinContent(bx, by);
                const double sMc = mcErr->GetBinContent(bx, by);
                deltaErr->SetBinContent(bx, by, std::sqrt(sData * sData + sMc * sMc));
            }
        }
    }
    TH2D* deltaRateUp = (TH2D*)deltaRate->Clone(("delta_rate_ele_" + year + "_Up").c_str());
    deltaRateUp->SetDirectory(nullptr);
    TH2D* deltaRateDown = (TH2D*)deltaRate->Clone(("delta_rate_ele_" + year + "_Down").c_str());
    deltaRateDown->SetDirectory(nullptr);
    if (deltaErr) {
        deltaRateUp->Add(deltaErr, 1.0);
        deltaRateDown->Add(deltaErr, -1.0);
        for (int bx = 1; bx <= deltaRateDown->GetNbinsX(); bx++)
            for (int by = 1; by <= deltaRateDown->GetNbinsY(); by++)
                if (deltaRateDown->GetBinContent(bx, by) < 0) deltaRateDown->SetBinContent(bx, by, 0);
    }

    TH2D* ratioNoSub = (TH2D*)rateDataRaw->Clone(("ratio_nosub_data_mc_" + year).c_str());
    ratioNoSub->SetDirectory(nullptr);
    ratioNoSub->Divide(rateMc);

    TH2D* ratioWithSub = (TH2D*)rateDataBkgsub->Clone(("ratio_withsub_data_mc_" + year).c_str());
    ratioWithSub->SetDirectory(nullptr);
    ratioWithSub->Divide(rateMc);

    gSystem->mkdir(dir.c_str(), true);
    TFile* fout = new TFile((dir + "charge_flip_delta_rate_" + year + ".root").c_str(), "RECREATE");
    rateDataRaw->Write();
    rateDataBkgsub->Write();
    TH2D* rateMcClone = (TH2D*)rateMc->Clone(("rate_ele_mc_" + year).c_str());
    rateMcClone->SetDirectory(fout);
    rateMcClone->Write();
    deltaRate->Write();
    deltaRateUp->Write();
    deltaRateDown->Write();
    ratioNoSub->Write();
    ratioWithSub->Write();
    fout->Close();

    std::cout << "=== delta_rate_ele_" << year << " = max(0, bkgsub-data - mc)  (pT, |eta|) ===" << std::endl;
    for (int bx = 1; bx <= deltaRate->GetNbinsX(); bx++) {
        for (int by = 1; by <= deltaRate->GetNbinsY(); by++) {
            std::cout << "  pt=[" << deltaRate->GetXaxis()->GetBinLowEdge(bx) << "," << deltaRate->GetXaxis()->GetBinUpEdge(bx)
                       << ") eta=[" << deltaRate->GetYaxis()->GetBinLowEdge(by) << "," << deltaRate->GetYaxis()->GetBinUpEdge(by)
                       << "): data_bkgsub=" << rateDataBkgsub->GetBinContent(bx, by)
                       << " mc=" << rateMc->GetBinContent(bx, by)
                       << " delta=" << deltaRate->GetBinContent(bx, by) << std::endl;
        }
    }
    std::cout << "Wrote " << dir << "charge_flip_delta_rate_" << year << ".root" << std::endl;
}
