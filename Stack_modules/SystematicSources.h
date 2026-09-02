#pragma once

#include <string>
#include <vector>

struct SystSource {
    std::string name;
    std::string upSuffix;
    std::string downSuffix;
};

static const std::vector<SystSource> kAllKnownSystSources = {
    {"frStat", "_frUp", "_frDown"},
    {"frEtauStat", "_frEtauStatUp", "_frEtauStatDown"},
    {"frDyMc", "_frDyMcUp", "_frDyMcDown"},
    {"frWJets", "_frWJetsUp", "_frWJetsDown"},
    {"frWJetsMc", "_frWJetsMcUp", "_frWJetsMcDown"},
    {"frQcdMc", "_frQcdMcUp", "_frQcdMcDown"},
    {"frTtMc", "_frTtMcUp", "_frTtMcDown"},
    {"zpt", "_zptUp", "_zptDown"},
    {"pu", "_puUp", "_puDown"},
    {"l1Prefire", "_l1PrefireUp", "_l1PrefireDown"},
    {"recoilResponse", "_recoilResponseUp", "_recoilResponseDown"},
    {"recoilResolution", "_recoilResolutionUp", "_recoilResolutionDown"},
    {"xsec", "_xsecUp", "_xsecDown"},
    {"roccor", "_roccorUp", "_roccorDown"},
    {"eReco", "_eRecoUp", "_eRecoDown"},
    {"eIdIso", "_eIdIsoUp", "_eIdIsoDown"},
    {"eTrig", "_eTrigUp", "_eTrigDown"},
    {"muId", "_muIdUp", "_muIdDown"},
    {"muIso", "_muIsoUp", "_muIsoDown"},
    {"muTrig", "_muTrigUp", "_muTrigDown"},
    {"tauVsEle", "_tauVsEleUp", "_tauVsEleDown"},
    {"tauVsMu", "_tauVsMuUp", "_tauVsMuDown"},
    {"tauVsJet", "_tauVsJetUp", "_tauVsJetDown"},
    {"tauES", "_tauESUp", "_tauESDown"}
};

static const std::vector<SystSource> kAllSystSources = {
    {"frStat", "_frUp", "_frDown"},
    {"frEtauStat", "_frEtauStatUp", "_frEtauStatDown"},
    {"frDyMc", "_frDyMcUp", "_frDyMcDown"},
    {"frWJets", "_frWJetsUp", "_frWJetsDown"},
    // {"frWJetsMc", "_frWJetsMcUp", "_frWJetsMcDown"},
    // {"frQcdMc", "_frQcdMcUp", "_frQcdMcDown"},
    // {"frTtMc", "_frTtMcUp", "_frTtMcDown"},
    {"zpt", "_zptUp", "_zptDown"},
    {"pu", "_puUp", "_puDown"},
    {"l1Prefire", "_l1PrefireUp", "_l1PrefireDown"},
    {"recoilResponse", "_recoilResponseUp", "_recoilResponseDown"},
    {"recoilResolution", "_recoilResolutionUp", "_recoilResolutionDown"},
    {"xsec", "_xsecUp", "_xsecDown"},
    {"roccor", "_roccorUp", "_roccorDown"},
    {"eReco", "_eRecoUp", "_eRecoDown"},
    {"eIdIso", "_eIdIsoUp", "_eIdIsoDown"},
    {"eTrig", "_eTrigUp", "_eTrigDown"},
    {"muId", "_muIdUp", "_muIdDown"},
    {"muIso", "_muIsoUp", "_muIsoDown"},
    {"muTrig", "_muTrigUp", "_muTrigDown"},
    {"tauVsEle", "_tauVsEleUp", "_tauVsEleDown"},
    {"tauVsMu", "_tauVsMuUp", "_tauVsMuDown"},
    {"tauVsJet", "_tauVsJetUp", "_tauVsJetDown"},
    {"tauES", "_tauESUp", "_tauESDown"}
};

inline bool endsWith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool isSystematicVariantName(const std::string& name, const std::vector<SystSource>& sources) {
    for (const auto& src : sources) {
        if (endsWith(name, src.upSuffix) || endsWith(name, src.downSuffix)) return true;
    }
    return false;
}
