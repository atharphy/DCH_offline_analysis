#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "TFile.h"
#include "TH1D.h"

inline TH1D* loadZPtWeights(const std::string& year) {
    const std::string fileName = "zpt_weights/" + year + "/ZPtWeights_" + year + ".root";
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open Z pT correction file:\n  " << fileName << std::endl;
        return nullptr;
    }
    TH1D* source = dynamic_cast<TH1D*>(file->Get("h_zpt_weight_ee"));
    if (!source) {
        std::cerr << "ERROR: missing h_zpt_weight_ee in:\n  " << fileName << std::endl;
        return nullptr;
    }
    TH1D* weights = dynamic_cast<TH1D*>(source->Clone("zpt_weights_ee"));
    weights->SetDirectory(nullptr);
    return weights;
}

inline bool getZPtWeight(const TH1D* weights, double zPt, double& value, double& error) {
    value = -1.0;
    error = 0.0;
    if (!weights || !std::isfinite(zPt)) return false;
    int bin = weights->GetXaxis()->FindFixBin(zPt);
    bin = std::max(1, std::min(bin, weights->GetNbinsX()));
    const double weight = weights->GetBinContent(bin);
    if (!std::isfinite(weight) || weight < 0.0) return false;
    const double weightError = weights->GetBinError(bin);
    value = weight;
    error = std::isfinite(weightError) && weightError >= 0.0 ? weightError : 0.0;
    return true;
}
