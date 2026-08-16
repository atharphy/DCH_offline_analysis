#pragma once

#include <string>
#include <utility>
#include <vector>

static const std::vector<std::pair<std::string, double>> kFlatXsecUncertainties = {
    {"DYJetsToLL", 0.0249},
    {"WWTo2L2Nu", 0.0577},
    {"WW_", 0.0577},
    {"WZTo2Q2L", 0.0412},
    {"WZTo3LNu", 0.0412},
    {"ZZTo2L2Nu", 0.0441},
    {"ZZTo2Q2L", 0.0441},
    {"ZZTo4L", 0.0441},
    {"ST_t-channel_antitop", 0.1882},
    {"ST_t-channel_top", 0.1429},
    {"ST_tW_top", 0.1105},
    {"ST_tW_antitop", 0.1105},
    {"ttWJets", 0.0747},
    {"ttZJets", 0.0822},
    {"ttHToTauTau", 0.0696},
    {"ttHToEE", 0.0696},
    {"ttHTo2L2Nu", 0.0696},
    {"ttHJetToNonbb", 0.0696},
    {"ZHToMuMu", 0.0410},
    {"ZHToTauTau", 0.0410},
    {"TTTo2L2Nu", 0.0319},
    {"TTToSemiLeptonic", 0.0319},
    {"TTToHadronic", 0.0319}
};

inline double getFlatXsecUncertainty(const std::string& baseName) {
    for (const auto& entry : kFlatXsecUncertainties) {
        if (baseName.rfind(entry.first, 0) == 0) return entry.second;
    }
    return 0.0;
}
