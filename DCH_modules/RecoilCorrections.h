#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TSystem.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "CommonConfig.h"

inline std::vector<std::string> recoilPayloadCandidates(const std::string& year) {
    const std::string directory = "HTT-utilities/RecoilCorrections/data/";
    if (USE_CUSTOM_RECOIL_CORRECTIONS) {
        if (year == "2016preVFP") return {directory + "TypeI-PFMet_2016preVFP.root"};
        if (year == "2016postVFP") return {directory + "TypeI-PFMet_2016postVFP.root"};
        if (year == "2017") return {directory + "TypeI-PFMet_2017.root"};
        if (year == "2018") return {directory + "TypeI-PFMet_2018.root"};
        return {};
    }
    if (year == "2016preVFP")
        return {directory + "TypeI-PFMet_Run2016BtoH.root", directory + "TypeI-PFMET_Run2016BtoH.root"};
    if (year == "2016postVFP")
        return {directory + "TypeI-PFMet_Run2016BtoH.root", directory + "TypeI-PFMET_Run2016BtoH.root"};
    if (year == "2017")
        return {directory + "TypeI-PFMET_2017.root", directory + "TypeI-PFMet_2017.root", directory + "TypeI-PFMet_Run2017.root", directory + "TypeI-PFMET_Run2017.root"};
    if (year == "2018")
        return {directory + "TypeI-PFMET_2018.root", directory + "TypeI-PFMet_2018.root", directory + "TypeI-PFMet_Run2018.root", directory + "TypeI-PFMET_Run2018.root"};
    return {};
}

inline std::unique_ptr<RecoilCorrector> loadRecoilCorrector(const std::string& year) {
    const std::vector<std::string> candidates = recoilPayloadCandidates(year);
    if (candidates.empty()) {
        std::cerr << "ERROR: cannot resolve recoil payload for " << year << std::endl;
        return nullptr;
    }
    for (const std::string& payload : candidates) {
        if (!gSystem->AccessPathName(payload.c_str())) {
            std::cout << "Loaded recoil corrections from:\n  " << payload << std::endl;
            return std::unique_ptr<RecoilCorrector>(new RecoilCorrector(payload.c_str()));
        }
    }
    std::cerr << "ERROR: no recoil payload found for " << year << ". Checked:" << std::endl;
    for (const std::string& payload : candidates) std::cerr << "  " << payload << std::endl;
    return nullptr;
}

inline bool applyRecoilCorrection(RecoilCorrector& corrector, double genPx, double genPy, double visPx, double visPy) {
    if (!std::isfinite(met) || !std::isfinite(metphi) || met < 0.0 || !std::isfinite(genPx) || !std::isfinite(genPy) || !std::isfinite(visPx) || !std::isfinite(visPy) || !std::isfinite(njets)) return false;
    const float metPx = static_cast<float>(met * std::cos(metphi));
    const float metPy = static_cast<float>(met * std::sin(metphi));
    float correctedMetPx = metPx;
    float correctedMetPy = metPy;
    const int recoilNJets = std::max(0, static_cast<int>(std::lround(njets)));
    corrector.CorrectByMeanResolution(metPx, metPy, static_cast<float>(genPx), static_cast<float>(genPy), static_cast<float>(visPx), static_cast<float>(visPy), recoilNJets, correctedMetPx, correctedMetPy);
    if (!std::isfinite(correctedMetPx) || !std::isfinite(correctedMetPy)) return false;
    met = std::hypot(correctedMetPx, correctedMetPy);
    metphi = std::atan2(correctedMetPy, correctedMetPx);
    return std::isfinite(met) && std::isfinite(metphi);
}
