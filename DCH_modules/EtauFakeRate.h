#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TauFakeRate.h"
#include "TauFRSystematics.h"

const std::string kEtauFakeRateFile = "Dependencies/etau_fake_rates/cut_results/ETau_FakeRates_Final.root";

inline std::shared_ptr<TauFRReader> loadEtauFRReader(const std::string& year) {
    std::unique_ptr<TFile> file(TFile::Open(kEtauFakeRateFile.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open etau fake-rate file: " << kEtauFakeRateFile << std::endl;
        return nullptr;
    }
    TH1D* barrel = dynamic_cast<TH1D*>(file->Get(("ETau_FR_" + year + "_barrel").c_str()));
    TH1D* endcap = dynamic_cast<TH1D*>(file->Get(("ETau_FR_" + year + "_endcap").c_str()));
    if (!barrel || !endcap) {
        std::cerr << "ERROR: missing ETau_FR_" << year << "_barrel/endcap in " << kEtauFakeRateFile << std::endl;
        return nullptr;
    }
    std::unique_ptr<TH2D> hist2D(buildFakeRate2DFromBarrelEndcap(*barrel, *endcap, "Data_" + year + "_etau_fake_rate"));
    return std::make_shared<TauFRReader>(hist2D.get());
}

inline std::shared_ptr<TauFRReader> getCachedEtauFRReader(const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<TauFRReader>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;
    auto reader = loadEtauFRReader(year);
    cache[year] = reader;
    return reader;
}

inline std::unique_ptr<TH1D> shiftFakeRateHistogramByError(const TH1D& source, const std::string& name, double sign) {
    std::unique_ptr<TH1D> result(dynamic_cast<TH1D*>(source.Clone(name.c_str())));
    result->SetDirectory(nullptr);
    for (int bin = 1; bin <= result->GetNbinsX(); ++bin) {
        const double shifted = std::max(0.0, std::min(1.0, source.GetBinContent(bin) + sign * source.GetBinError(bin)));
        result->SetBinContent(bin, shifted);
    }
    return result;
}

inline std::shared_ptr<TauFRReader> loadEtauFRStatReader(const std::string& year, const std::string& direction) {
    std::unique_ptr<TFile> file(TFile::Open(kEtauFakeRateFile.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open etau fake-rate file: " << kEtauFakeRateFile << std::endl;
        return nullptr;
    }
    TH1D* barrel = dynamic_cast<TH1D*>(file->Get(("ETau_FR_" + year + "_barrel").c_str()));
    TH1D* endcap = dynamic_cast<TH1D*>(file->Get(("ETau_FR_" + year + "_endcap").c_str()));
    if (!barrel || !endcap) {
        std::cerr << "ERROR: missing ETau_FR_" << year << "_barrel/endcap in " << kEtauFakeRateFile << std::endl;
        return nullptr;
    }
    const double sign = (direction == "up") ? 1.0 : -1.0;
    auto barrelShift = shiftFakeRateHistogramByError(*barrel, "ETau_FR_" + year + "_barrel_stat" + direction, sign);
    auto endcapShift = shiftFakeRateHistogramByError(*endcap, "ETau_FR_" + year + "_endcap_stat" + direction, sign);
    std::unique_ptr<TH2D> hist2D(buildFakeRate2DFromBarrelEndcap(*barrelShift, *endcapShift, "Data_" + year + "_etau_fake_rate_stat" + direction));
    return std::make_shared<TauFRReader>(hist2D.get());
}

inline std::shared_ptr<TauFRReader> getCachedEtauFRStatReader(const std::string& year, const std::string& direction) {
    static std::unordered_map<std::string, std::shared_ptr<TauFRReader>> cache;
    const std::string key = year + "_" + direction;
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    auto reader = loadEtauFRStatReader(year, direction);
    cache[key] = reader;
    return reader;
}
