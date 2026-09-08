#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TauFakeRate.h"

const std::string kTauFRSystematicsBase = "systematics/results";

const std::vector<std::pair<std::string,std::string>> kTauFRSystematicSources = {
    {"dy_mc", "DyMc"},
    {"wjets_data", "WJets"},
    {"wjets_mc", "WJetsMc"},
    {"qcd_mc", "QcdMc"},
    {"tt_mc", "TtMc"}
};

inline TH2D* buildFakeRate2DFromBarrelEndcap(const TH1D& barrel, const TH1D& endcap, const std::string& name) {
    const int nPtBins = barrel.GetNbinsX();
    std::vector<double> ptEdges(nPtBins + 1);
    for (int bin = 1; bin <= nPtBins; ++bin) ptEdges[bin - 1] = barrel.GetXaxis()->GetBinLowEdge(bin);
    ptEdges[nPtBins] = barrel.GetXaxis()->GetBinUpEdge(nPtBins);

    const double etaEdges[] = {0.0, 1.5, 2.5};
    TH2D* result = new TH2D(name.c_str(), ";p_{T}^{#tau_{h}} [GeV];|#eta_{#tau_{h}}|;Fake rate", nPtBins, ptEdges.data(), 2, etaEdges);
    result->SetDirectory(nullptr);
    result->Sumw2();

    for (int xbin = 1; xbin <= nPtBins; ++xbin) {
        result->SetBinContent(xbin, 1, barrel.GetBinContent(xbin));
        result->SetBinError(xbin, 1, barrel.GetBinError(xbin));
        result->SetBinContent(xbin, 2, endcap.GetBinContent(xbin));
        result->SetBinError(xbin, 2, endcap.GetBinError(xbin));
    }
    return result;
}

inline std::shared_ptr<TauFRReader> loadSystematicTauFRReader(const std::string& sourceName, const std::string& direction, const std::string& year) {
    const std::string path = kTauFRSystematicsBase + "/" + sourceName + "/" + year + "/tau_fake_rate_syst_" + sourceName + "_" + year + ".root";
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open fake-rate systematic file: " << path << std::endl;
        return nullptr;
    }
    TH1D* barrel = dynamic_cast<TH1D*>(file->Get(("FR_DY_" + year + "_barrel_" + direction).c_str()));
    TH1D* endcap = dynamic_cast<TH1D*>(file->Get(("FR_DY_" + year + "_endcap_" + direction).c_str()));
    if (!barrel || !endcap) {
        std::cerr << "ERROR: missing barrel/endcap histograms for " << sourceName << " " << direction << " " << year << " in " << path << std::endl;
        return nullptr;
    }
    std::unique_ptr<TH2D> hist2D(buildFakeRate2DFromBarrelEndcap(*barrel, *endcap, "Data_" + year + "_tau_fake_rate_" + sourceName + "_" + direction));
    return std::make_shared<TauFRReader>(hist2D.get());
}

inline std::shared_ptr<TauFRReader> getCachedSystematicTauFRReader(const std::string& sourceName, const std::string& direction, const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<TauFRReader>> cache;
    const std::string key = sourceName + "_" + direction + "_" + year;
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    auto reader = loadSystematicTauFRReader(sourceName, direction, year);
    cache[key] = reader;
    return reader;
}
