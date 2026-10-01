#include <TFile.h>
#include <TH2D.h>
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <TSystem.h>

#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "../filemap/FileMap.h"

namespace {
std::string baseNameNoExt(const std::string& path) {
    size_t slash = path.find_last_of('/');
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
}

const std::vector<double> PT_EDGES = {
    20., 22., 24., 26., 28., 30., 32., 35., 40.,
    45., 50., 60., 70., 85., 100., 125., 150., 200.
};
const double ETA_EDGES[3] = {0.0, 1.5, 2.5};

const double MIN_PAIRS = 20.0;
}

void MergeSSOSChunks(const std::string& year = "2018", const std::string& category = "data") {
    std::vector<std::string> nicknames;
    if (category == "data") nicknames = {"data"};
    else if (category == "mc") nicknames = {"DY"};
    else if (category == "bkg") nicknames = {"WZ", "WW", "VVV", "ttV", "WJ", "ZZ", "ST", "TTbar", "other"};
    auto fileMap = getFileMap(year);
    std::set<std::string> expectedSamples;
    for (auto& n : nicknames) {
        if (!fileMap.count(n)) continue;
        for (auto& f : fileMap[n]) {
            if (category == "data" && f.find("EGamma") == std::string::npos && f.find("SingleElectron") == std::string::npos) continue;
            expectedSamples.insert(baseNameNoExt(f));
        }
    }

    const std::string chunkDir = "DD_flip_rates_chunks/" + year;
    void* dirp = gSystem->OpenDirectory(chunkDir.c_str());
    if (!dirp) { std::cerr << "[error] cannot open " << chunkDir << std::endl; return; }

    TH2D* nOS = nullptr;
    TH2D* nSS = nullptr;
    int nChunks = 0;
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirp))) {
        std::string name(entry);
        if (name.size() < 5 || name.substr(name.size() - 5) != ".root") continue;

        size_t chunkPos = name.rfind("_chunk");
        if (chunkPos == std::string::npos) continue;
        const std::string sampleName = name.substr(0, chunkPos);
        if (!expectedSamples.count(sampleName)) continue;

        const std::string full = chunkDir + "/" + name;
        TFile* f = TFile::Open(full.c_str());
        if (!f || f->IsZombie()) { std::cerr << "  [warn] cannot open " << full << std::endl; continue; }
        TH2D* os = (TH2D*)f->Get("n_os_ele");
        TH2D* ss = (TH2D*)f->Get("n_ss_ele");
        if (!os || !ss) { std::cerr << "  [warn] missing hists in " << full << std::endl; f->Close(); continue; }

        if (!nOS) { nOS = (TH2D*)os->Clone("n_os_ele"); nOS->SetDirectory(nullptr); }
        else nOS->Add(os);
        if (!nSS) { nSS = (TH2D*)ss->Clone("n_ss_ele"); nSS->SetDirectory(nullptr); }
        else nSS->Add(ss);

        f->Close();
        nChunks++;
    }
    gSystem->FreeDirectory(dirp);

    if (!nOS) { std::cerr << "[error] no matching chunk files found for category=" << category << std::endl; return; }
    std::cout << "Merged " << nChunks << " chunk files for " << year << "/" << category << std::endl;

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
    std::cout << "Built " << nEq << " pair-category equations for " << nCells << " cell unknowns" << std::endl;
    if (nEq < nCells) {
        std::cerr << "[warn] fewer equations (" << nEq << ") than unknowns (" << nCells
                   << ") -- some cells' rates may be underconstrained (SVD pseudo-inverse handles this,"
                   << " but treat those cells' results with caution)." << std::endl;
    }

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
    if (!ok) std::cerr << "[warn] SVD solve reported a non-exact solution (expected for a rank-deficient system)" << std::endl;

    const int nPt = static_cast<int>(PT_EDGES.size()) - 1;
    TH2D* rate = new TH2D(("rate_ele_" + category + "_" + year).c_str(), "", nPt, PT_EDGES.data(), 2, ETA_EDGES);
    rate->SetDirectory(nullptr);
    for (int ptBin = 1; ptBin <= nPt; ptBin++) {
        for (int etaBin = 1; etaBin <= 2; etaBin++) {
            const int cell = (ptBin - 1) * 2 + etaBin;
            rate->SetBinContent(ptBin, etaBin, p(cell - 1));
        }
    }

    gSystem->mkdir(("DD_flip_rates/" + year).c_str(), true);
    TFile* fout = new TFile(("DD_flip_rates/" + year + "/ssos_ratio_" + category + "_" + year + ".root").c_str(), "RECREATE");
    nOS->Write();
    nSS->Write();
    rate->Write();
    fout->Close();

    std::cout << "=== rate_ele_" << category << "_" << year << " (pT, |eta|) ===" << std::endl;
    for (int bx = 1; bx <= rate->GetNbinsX(); bx++) {
        for (int by = 1; by <= rate->GetNbinsY(); by++) {
            std::cout << "  pt=[" << rate->GetXaxis()->GetBinLowEdge(bx) << "," << rate->GetXaxis()->GetBinUpEdge(bx)
                       << ") eta=[" << rate->GetYaxis()->GetBinLowEdge(by) << "," << rate->GetYaxis()->GetBinUpEdge(by)
                       << "): p=" << rate->GetBinContent(bx, by) << std::endl;
        }
    }
    std::cout << "Wrote DD_flip_rates/" << year << "/ssos_ratio_" << category << "_" << year << ".root" << std::endl;
}
