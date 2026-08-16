#pragma once

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"

#include "YieldConfig.h"

struct ProcessYields {
    std::unique_ptr<TH1D> nominal;
    std::map<std::string, std::unique_ptr<TH1D>> up;
    std::map<std::string, std::unique_ptr<TH1D>> down;
};

inline std::unique_ptr<TH1D> makeYieldHistogram(const std::string& name) {
    std::unique_ptr<TH1D> histogram(new TH1D(name.c_str(), "", 7, 0.5, 7.5));
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    for (int bin = 1; bin <= 7; ++bin) histogram->GetXaxis()->SetBinLabel(bin, kBinLabels[bin - 1].c_str());
    return histogram;
}

inline void addIntegral(TH1D& destination, int outputBin, TH1D* source) {
    if (!source) return;
    double error = 0.0;
    const double value = source->IntegralAndError(0, source->GetNbinsX() + 1, error);
    destination.SetBinContent(outputBin, destination.GetBinContent(outputBin) + value);
    destination.SetBinError(outputBin, std::hypot(destination.GetBinError(outputBin), error));
}

inline bool fillProcessYields(TH1D& destination, const std::vector<TFile*>& files, const std::string& region, const std::string& suffix = "") {
    bool foundAny = false;
    for (int outputBin = 1; outputBin <= 7; ++outputBin) {
        const std::string regionName = region + "_" + kBinSuffixes[outputBin - 1];
        for (TFile* file : files) {
            if (!file) continue;
            for (const std::string& finalState : kFinalStates) {
                const std::string histogramName = "h_LT_" + finalState + "_" + regionName + suffix;
                TH1D* source = dynamic_cast<TH1D*>(file->Get(histogramName.c_str()));
                if (source) foundAny = true;
                addIntegral(destination, outputBin, source);
            }
        }
    }
    return foundAny;
}
