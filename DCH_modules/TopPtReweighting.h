#pragma once

#include <algorithm>
#include <cmath>

namespace TopPtSF {
constexpr double kA = 0.103;
constexpr double kB = -0.0118;
constexpr double kC = -0.000134;
constexpr double kD = 0.973;
constexpr double kPtCap = 2000.0;

inline double sf(double pt) {
    const double ptEval = std::min(pt, kPtCap);
    return kA * std::exp(kB * ptEval) + kC * ptEval + kD;
}
}

inline bool findGenTopAntitopPt(double& topPt, double& antitopPt) {
    if (!GenPart_pdgId || !GenPart_pt || !GenPart_statusFlags) return false;
    topPt = -1.0;
    antitopPt = -1.0;

    for (size_t i = 0; i < GenPart_pdgId->size(); ++i) {
        const int pdg = (*GenPart_pdgId)[i];
        if (pdg != 6 && pdg != -6) continue;

        const int flags = (*GenPart_statusFlags)[i];
        const bool fromHardProcess = (flags >> 8) & 1;
        const bool isLastCopy = (flags >> 13) & 1;
        if (!fromHardProcess || !isLastCopy) continue;

        if (pdg == 6) topPt = (*GenPart_pt)[i];
        else antitopPt = (*GenPart_pt)[i];
    }
    return topPt >= 0.0 && antitopPt >= 0.0;
}

inline double topPtWeight(bool isData, bool isTTbar) {
    if (isData || !isTTbar) return 1.0;
    double topPt, antitopPt;
    if (!findGenTopAntitopPt(topPt, antitopPt)) return 1.0;
    return std::sqrt(TopPtSF::sf(topPt) * TopPtSF::sf(antitopPt));
}

inline double topPtWeightDownVariation() { return 1.0; }

inline double topPtWeightUpVariation(bool isData, bool isTTbar) {
    const double w = topPtWeight(isData, isTTbar);
    return w * w;
}
