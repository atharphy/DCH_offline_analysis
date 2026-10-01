#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TString.h>

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>

#include "../filemap/FileMap.h"
#include "../include/MyBranch.C"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#pragma GCC diagnostic ignored "-Wunused-value"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wuninitialized"
#include "../include/Kinematics.C"
#pragma GCC diagnostic pop

#include "../include/Xsections.C"
#include "../include/cat.h"
#include "../DCH_modules/ObjectAccessors.h"

namespace {
const Long64_t CHUNK_SIZE = 500000;

std::string baseNameOf(const std::string& path) {
    size_t slash = path.find_last_of('/');
    return (slash == std::string::npos) ? path : path.substr(slash + 1);
}
std::string baseNameNoExt(const std::string& path) {
    std::string base = baseNameOf(path);
    size_t dot = base.find_last_of('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
}

void restrictBranches(TTree* tree, bool isMC) {
    tree->SetBranchStatus("*", 0);
    const char* keep[] = {"cat", "pt_1", "pt_2", "eta_1", "eta_2", "phi_1", "phi_2", "q_1", "q_2"};
    for (const char* name : keep) if (tree->GetBranch(name)) tree->SetBranchStatus(name, 1);
    if (isMC && tree->GetBranch("Generator_weight")) tree->SetBranchStatus("Generator_weight", 1);
}

double getLumiFb(const std::string& year) {
    if (year == "2016preVFP") return 19520.0;
    if (year == "2016postVFP") return 16810.0;
    if (year == "2017") return 41480.0;
    if (year == "2018") return 59830.0;
    return 0.0;
}

const std::vector<double> PT_EDGES = {
    20., 22., 24., 26., 28., 30., 32., 35., 40.,
    45., 50., 60., 70., 85., 100., 125., 150., 200.
};
const double ETA_EDGES[3] = {0.0, 1.5, 2.5};
const double Z_MASS_LO = 81.0;
const double Z_MASS_HI = 101.0;

bool isElectronStream(const std::string& path) {
    return path.find("EGamma") != std::string::npos || path.find("SingleElectron") != std::string::npos;
}
}

void DeriveSSOSRatioChunkBatch(const std::string& year, int idx, int batchStart, int numInBatch) {
    const std::map<std::string, std::vector<std::string>> fileMap = getFileMap(year);
    std::vector<std::pair<std::string, std::string>> files;
    for (const auto& sample : fileMap)
        for (const std::string& fileName : sample.second) files.push_back({sample.first, fileName});
    if (idx < 0 || idx >= static_cast<int>(files.size())) {
        std::cerr << "[error] invalid idx=" << idx << ", number of files=" << files.size() << std::endl;
        return;
    }

    static const std::vector<std::string> bkgNicknames = {
        "WZ", "WW", "VVV", "ttV", "WJ", "ZZ", "ST", "TTbar", "other"
    };

    const std::string nickname = files[idx].first;
    const std::string fname = files[idx].second;
    const bool isDataFile = (nickname == "data");
    const bool isDY = (nickname == "DY");
    const bool isBkg = std::find(bkgNicknames.begin(), bkgNicknames.end(), nickname) != bkgNicknames.end();
    const bool isMC = isDY || isBkg;
    if (!((isDataFile && isElectronStream(fname)) || isMC)) {

        return;
    }
    const std::string category = isDataFile ? "data" : (isDY ? "mc" : "bkg");

    const double lumi = getLumiFb(year);
    const int nPt = static_cast<int>(PT_EDGES.size()) - 1;
    const int nCells = nPt * 2;

    TH2D* nOS = new TH2D("n_os_ele", "", nCells, 0.5, nCells + 0.5, nCells, 0.5, nCells + 0.5);
    TH2D* nSS = new TH2D("n_ss_ele", "", nCells, 0.5, nCells + 0.5, nCells, 0.5, nCells + 0.5);
    nOS->Sumw2();
    nSS->Sumw2();
    nOS->SetDirectory(nullptr);
    nSS->SetDirectory(nullptr);

    TH2D* cellFinder = new TH2D("cell_finder", "", nPt, PT_EDGES.data(), 2, ETA_EDGES);
    cellFinder->SetDirectory(nullptr);
    auto flattenCell = [&](double pt, double aeta) -> int {
        const int ptBin = cellFinder->GetXaxis()->FindFixBin(pt);
        const int etaBin = cellFinder->GetYaxis()->FindFixBin(aeta);
        return (ptBin - 1) * 2 + etaBin;
    };

    TFile* fin = TFile::Open(fname.c_str());
    if (!fin || fin->IsZombie()) { std::cerr << "[error] cannot open " << fname << std::endl; return; }
    TTree* tree = (TTree*)fin->Get("Events");
    if (!tree) { std::cerr << "[error] no Events tree in " << fname << std::endl; fin->Close(); return; }

    MyBranch(tree);
    restrictBranches(tree, isMC);
    tree->SetCacheSize(64 * 1024 * 1024);
    tree->SetCacheLearnEntries(100);
    tree->AddBranchToCache("*", true);

    double xsw = 1.0;
    if (isMC) {
        TH1D* hNWEvts = (TH1D*)fin->Get("hNWEvts");
        if (!hNWEvts) hNWEvts = (TH1D*)fin->Get("hNEvts");
        const double denom = hNWEvts ? hNWEvts->Integral() : 0.0;
        xsw = (denom > 0.0) ? lumi * XSec(baseNameOf(fname)) / denom : 0.0;
        if (xsw <= 0.0) std::cerr << "[warn] xsw<=0 for " << fname << std::endl;
    }

    const Long64_t total = tree->GetEntries();
    if (!isMC || xsw > 0.0) {
        for (int c = batchStart; c < batchStart + numInBatch; ++c) {
            const Long64_t firstEntry = static_cast<Long64_t>(c) * CHUNK_SIZE;
            if (firstEntry >= total) break;
            const Long64_t last = std::min(total, firstEntry + CHUNK_SIZE);
            for (Long64_t i = firstEntry; i < last; i++) {
                tree->GetEntry(i);
                const std::string catstr = numberToCat(cat);
                if (catstr != "ee") continue;

                const double mass = (LepV(1) + LepV(2)).M();
                if (mass < Z_MASS_LO || mass > Z_MASS_HI) continue;

                const double pt1 = ptByIndex(1), pt2 = ptByIndex(2);
                const double aeta1 = std::fabs(etaByIndex(1)), aeta2 = std::fabs(etaByIndex(2));
                if (pt1 < PT_EDGES.front() || pt2 < PT_EDGES.front()) continue;
                if (pt1 >= PT_EDGES.back() || pt2 >= PT_EDGES.back()) continue;

                const int cell1 = flattenCell(pt1, aeta1);
                const int cell2 = flattenCell(pt2, aeta2);
                const int cellLo = std::min(cell1, cell2);
                const int cellHi = std::max(cell1, cell2);

                const double weight = isMC ? (Generator_weight * xsw) : 1.0;
                const bool isOS = (q_1 * q_2 < 0);
                (isOS ? nOS : nSS)->Fill(cellLo, cellHi, weight);
            }
        }
    }
    fin->Close();

    const char* scratchEnv = gSystem->Getenv("_CONDOR_SCRATCH_DIR");
    const std::string scratch = (scratchEnv && scratchEnv[0]) ? scratchEnv : ".";
    const TString eosDir = Form("DD_flip_rates_chunks/%s", year.c_str());
    gSystem->mkdir(eosDir, kTRUE);

    const std::string sampleName = baseNameNoExt(fname);
    const std::string localOut = scratch + "/" + sampleName + "_chunk" + std::to_string(batchStart) + ".root";
    TFile* fout = new TFile(localOut.c_str(), "RECREATE");
    nOS->Write();
    nSS->Write();
    fout->Close();

    const TString eosOut = eosDir + "/" + TString(sampleName) + Form("_chunk%d.root", batchStart);
    TFile::Cp(localOut.c_str(), eosOut, kFALSE);
    std::cout << "Wrote " << eosOut << " (category=" << category << ")" << std::endl;
}
