#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <string>

#include "EtauCommon.h"

void DeriveEtauFakeRateChunk(std::string year, std::string filename, Long64_t firstEntry, Long64_t nEntries, std::string outputPath) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    ROOT::RDataFrame frame("Events", filename);
    auto ranged = frame.Range(firstEntry, firstEntry + nEntries);

    auto tag = ranged
        .Filter(etauLeg1IsElectron, {"cat"})
        .Filter("genPartFlav_1 == 1");

    auto probe = tag
        .Define("probeMatch", etauFindProbe, {"eta_1", "phi_1", "lpt", "leta", "lphi", "lflavor", "gen_match", "TauIDm", "TauIDe", "TauIDj"})
        .Define("probePt", "probeMatch[0]")
        .Define("probeEta", "probeMatch[1]")
        .Define("probeTauIDj", "probeMatch[2]")
        .Filter("probePt > 0");

    auto denom = probe.Filter("probeTauIDj > 2");
    auto num = denom.Filter("probeTauIDj > 16");

    const int nPtBins = static_cast<int>(ETAU_PT_EDGES.size()) - 1;

    auto denBarrel = denom.Filter("abs(probeEta) < 1.5").Histo1D({"denBarrel", "", nPtBins, ETAU_PT_EDGES.data()}, "probePt");
    auto denEndcap = denom.Filter("abs(probeEta) >= 1.5 && abs(probeEta) < 2.5").Histo1D({"denEndcap", "", nPtBins, ETAU_PT_EDGES.data()}, "probePt");
    auto numBarrel = num.Filter("abs(probeEta) < 1.5").Histo1D({"numBarrel", "", nPtBins, ETAU_PT_EDGES.data()}, "probePt");
    auto numEndcap = num.Filter("abs(probeEta) >= 1.5 && abs(probeEta) < 2.5").Histo1D({"numEndcap", "", nPtBins, ETAU_PT_EDGES.data()}, "probePt");

    gSystem->mkdir(gSystem->DirName(outputPath.c_str()), kTRUE);
    TFile output(outputPath.c_str(), "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "ERROR: cannot create " << outputPath << std::endl;
        gSystem->Exit(1);
        return;
    }

    denBarrel->Write("denBarrel");
    denEndcap->Write("denEndcap");
    numBarrel->Write("numBarrel");
    numEndcap->Write("numEndcap");
    output.Close();

    std::cout << "[" << year << "] " << filename << " [" << firstEntry << "," << firstEntry + nEntries << ") -> " << outputPath
              << " (den=" << denBarrel->GetEntries() + denEndcap->GetEntries()
              << ", num=" << numBarrel->GetEntries() + numEndcap->GetEntries() << ")" << std::endl;
}
