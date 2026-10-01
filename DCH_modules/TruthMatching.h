#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "ObjectAccessors.h"

inline bool getPromptTruthDileptonXY(const std::string& catstr, double& px, double& py) {
    std::vector<int> promptLeptons;
    promptLeptons.reserve(4);
    for (int idx = 1; idx <= static_cast<int>(catstr.size()); ++idx) {
        const char flavor = catstr[idx - 1];
        if (flavor != 'e' && flavor != 'm') continue;
        if (genPartFlavByIndex(idx) != 1) continue;
        const double pt = truthPtByIndex(idx);
        const double phi = truthPhiByIndex(idx);
        if (!std::isfinite(pt) || !std::isfinite(phi) || pt <= 0.0) continue;
        promptLeptons.push_back(idx);
    }
    if (promptLeptons.size() != 2) return false;
    const int first = promptLeptons[0];
    const int second = promptLeptons[1];
    if (catstr[first - 1] != catstr[second - 1]) return false;
    px = truthPtByIndex(first) * std::cos(truthPhiByIndex(first)) + truthPtByIndex(second) * std::cos(truthPhiByIndex(second));
    py = truthPtByIndex(first) * std::sin(truthPhiByIndex(first)) + truthPtByIndex(second) * std::sin(truthPhiByIndex(second));
    return std::isfinite(px) && std::isfinite(py);
}

inline bool getPromptTruthDileptonPt(const std::string& catstr, double& zPt) {
    double px = 0.0;
    double py = 0.0;
    if (!getPromptTruthDileptonXY(catstr, px, py)) return false;
    zPt = std::hypot(px, py);
    return std::isfinite(zPt);
}

inline bool getPromptTruthLeptonXY(const std::string& catstr, double& px, double& py) {
    int promptIndex = -1;
    for (int idx = 1; idx <= static_cast<int>(catstr.size()); ++idx) {
        const char flavor = catstr[idx - 1];
        if (flavor != 'e' && flavor != 'm') continue;
        if (genPartFlavByIndex(idx) != 1) continue;
        const double pt = truthPtByIndex(idx);
        const double phi = truthPhiByIndex(idx);
        if (!std::isfinite(pt) || !std::isfinite(phi) || pt <= 0.0) continue;
        if (promptIndex != -1) return false;
        promptIndex = idx;
    }
    if (promptIndex == -1) return false;
    px = truthPtByIndex(promptIndex) * std::cos(truthPhiByIndex(promptIndex));
    py = truthPtByIndex(promptIndex) * std::sin(truthPhiByIndex(promptIndex));
    return std::isfinite(px) && std::isfinite(py);
}
