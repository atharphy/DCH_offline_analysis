#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "TH1D.h"

inline double clampRate(double value) {
    if (!std::isfinite(value)) return 0.0;
    if (value < 0.0) return 0.0;
    if (value > 1.0) return 1.0;
    return value;
}

inline void makeEnvelopeVariation(const TH1D& nominal, const std::vector<const TH1D*>& alternates, const std::string& upName, const std::string& downName, std::unique_ptr<TH1D>& up, std::unique_ptr<TH1D>& down) {
    up.reset(dynamic_cast<TH1D*>(nominal.Clone(upName.c_str())));
    down.reset(dynamic_cast<TH1D*>(nominal.Clone(downName.c_str())));
    up->SetDirectory(nullptr);
    down->SetDirectory(nullptr);
    for (int bin = 1; bin <= nominal.GetNbinsX(); ++bin) {
        const double n = nominal.GetBinContent(bin);
        double maxAbsDiff = 0.0;
        for (const TH1D* alt : alternates) maxAbsDiff = std::max(maxAbsDiff, std::fabs(alt->GetBinContent(bin) - n));
        up->SetBinContent(bin, clampRate(n + maxAbsDiff));
        down->SetBinContent(bin, clampRate(n - maxAbsDiff));
    }
}

inline void makeSymmetricVariation(const TH1D& nominal, const TH1D& alt, const std::string& upName, const std::string& downName, std::unique_ptr<TH1D>& up, std::unique_ptr<TH1D>& down) {
    makeEnvelopeVariation(nominal, {&alt}, upName, downName, up, down);
}
