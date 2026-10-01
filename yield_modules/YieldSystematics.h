#pragma once

#include <string>
#include <vector>

struct SystSource {
    std::string name;
    std::string upSuffix;
    std::string downSuffix;
};

const std::vector<SystSource> kAllSystSources = {
    {"frStat", "_frUp", "_frDown"},
    {"frDyMc", "_frDyMcUp", "_frDyMcDown"},
    {"frWJets", "_frWJetsUp", "_frWJetsDown"},

    {"zpt", "_zptUp", "_zptDown"},
    {"pu", "_puUp", "_puDown"},
    {"l1Prefire", "_l1PrefireUp", "_l1PrefireDown"},
    {"metPhiCorr", "_metPhiCorrUp", "_metPhiCorrDown"},
    {"recoilResponse", "_recoilResponseUp", "_recoilResponseDown"},
    {"recoilResolution", "_recoilResolutionUp", "_recoilResolutionDown"},
    {"xsec", "_xsecUp", "_xsecDown"},
    {"roccor", "_roccorUp", "_roccorDown"}
};
