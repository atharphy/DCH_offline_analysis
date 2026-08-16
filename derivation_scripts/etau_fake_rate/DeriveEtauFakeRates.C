#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TIterator.h>
#include <TList.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "EtauCommon.h"

namespace {

void savePlot(TH1D* histogram, const std::string& outputDirectory, const std::string& name) {
    TCanvas canvas(("c_" + name).c_str(), "", 800, 700);
    canvas.SetGrid();
    histogram->SetMarkerStyle(20);
    histogram->Draw("E1");
    canvas.SaveAs((outputDirectory + "/" + name + ".png").c_str());
}

std::vector<std::string> listChunkFiles(const std::string& directory) {
    std::vector<std::string> result;
    TSystemDirectory dir("chunks", directory.c_str());
    TList* files = dir.GetListOfFiles();
    if (!files) return result;
    TIter next(files);
    while (TObject* obj = next()) {
        auto* file = dynamic_cast<TSystemFile*>(obj);
        if (!file || file->IsDirectory()) continue;
        const std::string name = file->GetName();
        if (name.size() > 5 && name.substr(name.size() - 5) == ".root") result.push_back(directory + "/" + name);
    }
    return result;
}

std::unique_ptr<TH1D> readAndSum(const std::vector<std::string>& years, const std::string& prefix,
                                  const std::string& suffix, const std::string& outputName) {
    std::unique_ptr<TH1D> sum;
    for (const std::string& inputYear : years) {
        const std::string fileName = ETAU_OUTPUT_BASE + "/" + inputYear + "/data_etau_fake_rates_" + inputYear + ".root";
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
        if (!file || file->IsZombie()) { std::cerr << "ERROR: cannot open " << fileName << std::endl; return nullptr; }
        const std::string histogramName = prefix + inputYear + suffix;
        TH1D* source = dynamic_cast<TH1D*>(file->Get(histogramName.c_str()));
        if (!source) { std::cerr << "ERROR: missing " << histogramName << " in " << fileName << std::endl; return nullptr; }
        if (!sum) sum = etauCloneHistogram(*source, outputName);
        else sum->Add(source);
    }
    return sum;
}

void mergeYear(const std::string& year) {
    const std::string chunkDir = ETAU_OUTPUT_BASE + "/" + year + "/chunks";
    const std::vector<std::string> chunkFiles = listChunkFiles(chunkDir);
    if (chunkFiles.empty()) {
        std::cerr << "ERROR: no chunk files found in " << chunkDir << std::endl;
        gSystem->Exit(2);
        return;
    }

    std::unique_ptr<TH1D> denBarrel, denEndcap, numBarrel, numEndcap;
    int nMerged = 0;
    for (const std::string& chunkFile : chunkFiles) {
        std::unique_ptr<TFile> f(TFile::Open(chunkFile.c_str(), "READ"));
        if (!f || f->IsZombie()) { std::cerr << "[warning] cannot open " << chunkFile << std::endl; continue; }
        TH1D* db = dynamic_cast<TH1D*>(f->Get("denBarrel"));
        TH1D* de = dynamic_cast<TH1D*>(f->Get("denEndcap"));
        TH1D* nb = dynamic_cast<TH1D*>(f->Get("numBarrel"));
        TH1D* ne = dynamic_cast<TH1D*>(f->Get("numEndcap"));
        if (!db || !de || !nb || !ne) { std::cerr << "[warning] missing histograms in " << chunkFile << std::endl; continue; }

        if (!denBarrel) {
            denBarrel = etauCloneHistogram(*db, "Den_et_" + year + "_barrel");
            denEndcap = etauCloneHistogram(*de, "Den_et_" + year + "_endcap");
            numBarrel = etauCloneHistogram(*nb, "Num_et_" + year + "_barrel");
            numEndcap = etauCloneHistogram(*ne, "Num_et_" + year + "_endcap");
        } else {
            denBarrel->Add(db);
            denEndcap->Add(de);
            numBarrel->Add(nb);
            numEndcap->Add(ne);
        }
        ++nMerged;
    }
    std::cout << "[" << year << "] merged " << nMerged << "/" << chunkFiles.size() << " chunk files" << std::endl;
    if (!denBarrel) { gSystem->Exit(3); return; }

    auto fakeRateBarrel = etauMakeFakeRate(*numBarrel, *denBarrel, "FR_et_" + year + "_barrel");
    auto fakeRateEndcap = etauMakeFakeRate(*numEndcap, *denEndcap, "FR_et_" + year + "_endcap");

    const int nPtBins = static_cast<int>(ETAU_PT_EDGES.size()) - 1;
    const double etaEdges[] = {0.0, 1.5, 2.5};
    TH2D payload(("Data_" + year + "_etau_fake_rate").c_str(), ";p_{T}^{#tau} [GeV];|#eta^{#tau}|;Electron-to-tau fake rate",
                 nPtBins, ETAU_PT_EDGES.data(), 2, etaEdges);
    payload.Sumw2();
    etauCopyToTwoDimensionalPayload(payload, *fakeRateBarrel, *fakeRateEndcap);

    const std::string outputDirectory = ETAU_OUTPUT_BASE + "/" + year;
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);
    const std::string outputName = outputDirectory + "/data_etau_fake_rates_" + year + ".root";

    TFile output(outputName.c_str(), "RECREATE");
    if (output.IsZombie()) { std::cerr << "ERROR: cannot create " << outputName << std::endl; gSystem->Exit(4); return; }

    payload.Write();
    fakeRateBarrel->Write();
    fakeRateEndcap->Write();
    denBarrel->Write();
    denEndcap->Write();
    numBarrel->Write();
    numEndcap->Write();
    output.Close();

    savePlot(fakeRateBarrel.get(), outputDirectory, "FR_et_" + year + "_barrel");
    savePlot(fakeRateEndcap.get(), outputDirectory, "FR_et_" + year + "_endcap");

    std::cout << "[" << year << "] output: " << outputName << std::endl;
}

void combineRun2() {
    const std::vector<std::string> years = {"2016preVFP", "2016postVFP", "2017", "2018"};

    auto denBarrel = readAndSum(years, "Den_et_", "_barrel", "Den_et_Run2_barrel");
    auto denEndcap = readAndSum(years, "Den_et_", "_endcap", "Den_et_Run2_endcap");
    auto numBarrel = readAndSum(years, "Num_et_", "_barrel", "Num_et_Run2_barrel");
    auto numEndcap = readAndSum(years, "Num_et_", "_endcap", "Num_et_Run2_endcap");

    if (!denBarrel || !denEndcap || !numBarrel || !numEndcap) { gSystem->Exit(5); return; }

    auto fakeRateBarrel = etauMakeFakeRate(*numBarrel, *denBarrel, "FR_et_Run2_barrel");
    auto fakeRateEndcap = etauMakeFakeRate(*numEndcap, *denEndcap, "FR_et_Run2_endcap");

    const int nPtBins = static_cast<int>(ETAU_PT_EDGES.size()) - 1;
    const double etaEdges[] = {0.0, 1.5, 2.5};
    TH2D payload("Data_Run2_etau_fake_rate", ";p_{T}^{#tau} [GeV];|#eta^{#tau}|;Electron-to-tau fake rate",
                 nPtBins, ETAU_PT_EDGES.data(), 2, etaEdges);
    payload.Sumw2();
    etauCopyToTwoDimensionalPayload(payload, *fakeRateBarrel, *fakeRateEndcap);

    const std::string outputDirectory = ETAU_OUTPUT_BASE + "/Run2";
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);
    const std::string outputName = ETAU_OUTPUT_BASE + "/data_etau_fake_rates_Run2.root";

    TFile output(outputName.c_str(), "RECREATE");
    if (output.IsZombie()) { std::cerr << "ERROR: cannot create " << outputName << std::endl; gSystem->Exit(6); return; }
    payload.Write();
    denBarrel->Write();
    denEndcap->Write();
    numBarrel->Write();
    numEndcap->Write();
    fakeRateBarrel->Write();
    fakeRateEndcap->Write();
    output.Close();

    savePlot(fakeRateBarrel.get(), outputDirectory, "FR_et_Run2_barrel");
    savePlot(fakeRateEndcap.get(), outputDirectory, "FR_et_Run2_endcap");

    std::cout << "[Run2] output: " << outputName << std::endl;
}

}

void DeriveEtauFakeRates(std::string year = "2018") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    if (year == "Run2") {
        combineRun2();
        return;
    }

    const std::vector<std::string> allowedYears = {"2016preVFP", "2016postVFP", "2017", "2018"};
    if (std::find(allowedYears.begin(), allowedYears.end(), year) == allowedYears.end()) {
        std::cerr << "ERROR: unsupported year " << year << std::endl;
        gSystem->Exit(1);
        return;
    }

    mergeYear(year);
}
