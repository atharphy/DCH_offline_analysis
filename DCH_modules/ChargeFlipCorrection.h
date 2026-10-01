#pragma once

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH2D.h"
#include "TRandom3.h"

#include "ObjectAccessors.h"

const bool APPLY_CHARGE_FLIP_CORRECTION = true;

inline TH2D* loadChargeFlipRate(const std::string& year, const std::string& fileNameOverride = "", const std::string& histNameOverride = "") {
    const std::string fileName = !fileNameOverride.empty() ? fileNameOverride
        : "charge_flip/DD_flip_rates/" + year + "/charge_flip_delta_rate_" + year + ".root";
    const std::string histName = !histNameOverride.empty() ? histNameOverride : "delta_rate_ele_" + year;
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open charge mis-ID rate file:\n  " << fileName << std::endl;
        return nullptr;
    }
    TH2D* source = dynamic_cast<TH2D*>(file->Get(histName.c_str()));
    if (!source) {
        std::cerr << "ERROR: missing " << histName << " in:\n  " << fileName << std::endl;
        return nullptr;
    }
    TH2D* rate = dynamic_cast<TH2D*>(source->Clone(("charge_flip_rate_ele_" + histName).c_str()));
    rate->SetDirectory(nullptr);
    return rate;
}

inline double getChargeFlipRate(const TH2D& rate, double pt, double aeta) {
    const double ptClamped = std::min(pt, rate.GetXaxis()->GetXmax() - 0.1);
    const double etaClamped = std::min(aeta, rate.GetYaxis()->GetXmax() - 0.01);
    const int bx = rate.GetXaxis()->FindFixBin(ptClamped);
    const int by = rate.GetYaxis()->FindFixBin(etaClamped);
    return rate.GetBinContent(bx, by);
}

struct ChargeFlipDecision {
    int idx;
    bool flipped;
    double fNominal;
};

inline std::vector<ChargeFlipDecision> applyChargeFlipCorrection(const std::string& cat, bool isData, const TH2D& rate, TRandom3& rng) {
    std::vector<ChargeFlipDecision> decisions;
    if (isData || !APPLY_CHARGE_FLIP_CORRECTION) return decisions;
    for (int idx = 1; idx <= static_cast<int>(cat.size()) && idx <= 4; ++idx) {
        if (cat[idx - 1] != 'e') continue;
        const double pt = ptByIndex(idx);
        const double aeta = std::fabs(etaByIndex(idx));
        const double f = getChargeFlipRate(rate, pt, aeta);
        const bool flip = rng.Uniform() < f;
        decisions.push_back({idx, flip, f});
        if (!flip) continue;
        if (idx == 1) q_1 = -q_1;
        else if (idx == 2) q_2 = -q_2;
        else if (idx == 3) q_3 = -q_3;
        else if (idx == 4) q_4 = -q_4;
    }
    return decisions;
}

inline double chargeFlipWeightRatio(const std::vector<ChargeFlipDecision>& decisions, const TH2D& shiftedRate) {
    double ratio = 1.0;
    for (const auto& d : decisions) {
        const double pt = ptByIndex(d.idx);
        const double aeta = std::fabs(etaByIndex(d.idx));
        const double fShift = getChargeFlipRate(shiftedRate, pt, aeta);
        if (d.flipped) {
            if (d.fNominal > 0.0) ratio *= fShift / d.fNominal;
        } else {
            const double denom = 1.0 - d.fNominal;
            if (denom > 1e-6) ratio *= (1.0 - fShift) / denom;
        }
    }
    return ratio;
}
