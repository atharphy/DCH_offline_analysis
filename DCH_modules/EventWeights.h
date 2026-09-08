#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "ObjectAccessors.h"
#include "SystematicPlan.h"
#include "TopPtReweighting.h"

inline double computeBaseGenWeight(double xsw) {
    double w = Generator_weight * brWeight * xsw;
    w *= L1PreFiringWeight_Nom * weightPUtruejson;
    return w;
}

inline void appendPileupWeightVariants(double evtWeight, std::vector<WeightSystematic>& weightVariants) {
    if (weightPUtruejson <= 0.0) return;
    weightVariants.push_back({"_puUp", evtWeight / weightPUtruejson * weightPUtruejson_up});
    weightVariants.push_back({"_puDown", evtWeight / weightPUtruejson * weightPUtruejson_down});
}

inline void appendL1PrefiringWeightVariants(double evtWeight, std::vector<WeightSystematic>& weightVariants) {
    if (L1PreFiringWeight_Nom <= 0.0) return;
    weightVariants.push_back({"_l1PrefireUp", evtWeight / L1PreFiringWeight_Nom * L1PreFiringWeight_Up});
    weightVariants.push_back({"_l1PrefireDown", evtWeight / L1PreFiringWeight_Nom * L1PreFiringWeight_Down});
}

inline double computeTriggerSF(const std::string& catstr) {
    const int nLep = static_cast<int>(catstr.size());
    double trigSF = 1.0;
    const bool hasElectron = catstr.find('e') != std::string::npos;
    const bool hasMuon = catstr.find('m') != std::string::npos;

    if (hasElectron && hasMuon) {
        double electronSF = 0.0;
        double muonSF = 0.0;
        for (int idx = 1; idx <= nLep; ++idx) {
            const char flavor = catstr[idx - 1];
            if (flavor == 't' || !trigPassedByIndex(idx)) continue;
            if (flavor == 'e' && electronSF == 0.0) electronSF = trigSFByIndex(idx);
            else if (flavor == 'm' && muonSF == 0.0) muonSF = trigSFByIndex(idx);
        }
        if (electronSF > 0.0 && muonSF > 0.0) trigSF = electronSF + muonSF - electronSF * muonSF;
        else if (electronSF > 0.0) trigSF = electronSF;
        else if (muonSF > 0.0) trigSF = muonSF;
    }
    else {
        if (nLep < 3) {
            if (isTrig_1 >= 1) trigSF = TrigSF_1;
            else if (isTrig_1 == -1) trigSF = TrigSF_2;
        }
        else {
            if (isTrig_1 >= 1 && isTrig_2 == 0) trigSF = TrigSF_1;
            else if (isTrig_1 == -1 && isTrig_2 == 0) trigSF = TrigSF_2;
            else if (isTrig_2 >= 1 && isTrig_1 == 0) trigSF = TrigSF_3;
            else if (isTrig_2 == -1 && isTrig_1 == 0) trigSF = TrigSF_4;
            else if (isTrig_1 == 2 && isTrig_2 == 2) trigSF = TrigSF_1;
        }
    }
    return trigSF;
}

inline double computeNominalWeight(const std::string& catstr, bool isData, bool isTTbar, double xsw) {
    if (isData) return 1.0;
    const int nLep = static_cast<int>(catstr.size());
    double objectSF = 1.0;
    for (int idx = 1; idx <= nLep; ++idx) {
        const char flavor = catstr[idx - 1];
        if (flavor == 'e' || flavor == 'm') {
            objectSF *= idSFByIndex(idx);
            objectSF *= isoSFByIndex(idx);
        }
        else if (flavor == 't') {
            objectSF *= tauEleSFByIndex(idx);
            objectSF *= tauMuSFByIndex(idx);
            objectSF *= tauJetSFByIndex(idx);
        }
    }
    return computeBaseGenWeight(xsw) * objectSF * computeTriggerSF(catstr) * topPtWeight(isData, isTTbar);
}

struct ObjectSFSystEntry { std::string suffix; double weight; };

inline std::vector<ObjectSFSystEntry> computeObjectSFSystematics(const std::string& catstr, bool isData, double baseWeight) {
    std::vector<ObjectSFSystEntry> result;
    if (isData) return result;

    double eRecoUp = 1.0, eRecoDown = 1.0, eIdIsoUp = 1.0, eIdIsoDown = 1.0;
    double muIdUp = 1.0, muIdDown = 1.0, muIsoUp = 1.0, muIsoDown = 1.0;
    double tauVsEleUp = 1.0, tauVsEleDown = 1.0, tauVsMuUp = 1.0, tauVsMuDown = 1.0, tauVsJetUp = 1.0, tauVsJetDown = 1.0;

    const int nLep = static_cast<int>(catstr.size());
    for (int idx = 1; idx <= nLep; ++idx) {
        const char flavor = catstr[idx - 1];
        const double pt = ptByIndex(idx);
        if (flavor == 'e' && pt > 20.0) {
            const double idNom = idSFByIndex(idx), idUp = idSFUpByIndex(idx), idDown = idSFDownByIndex(idx);
            if (idNom > 0.0 && idUp > 0.0 && idDown > 0.0) { eRecoUp *= idUp / idNom; eRecoDown *= idDown / idNom; }
            const double isoNom = isoSFByIndex(idx), isoUp = isoSFUpByIndex(idx), isoDown = isoSFDownByIndex(idx);
            if (isoNom > 0.0 && isoUp > 0.0 && isoDown > 0.0) { eIdIsoUp *= isoUp / isoNom; eIdIsoDown *= isoDown / isoNom; }
        }
        else if (flavor == 'm' && pt > 15.0) {
            const double idNom = idSFByIndex(idx), idUp = idSFUpByIndex(idx), idDown = idSFDownByIndex(idx);
            if (idNom > 0.0 && idUp > 0.0 && idDown > 0.0) { muIdUp *= idUp / idNom; muIdDown *= idDown / idNom; }
            const double isoNom = isoSFByIndex(idx), isoUp = isoSFUpByIndex(idx), isoDown = isoSFDownByIndex(idx);
            if (isoNom > 0.0 && isoUp > 0.0 && isoDown > 0.0) { muIsoUp *= isoUp / isoNom; muIsoDown *= isoDown / isoNom; }
        }
        else if (flavor == 't') {
            const double vE = tauEleSFByIndex(idx), vEUp = tauEleSFUpByIndex(idx), vEDown = tauEleSFDownByIndex(idx);
            if (vE > 0.0 && vEUp > 0.0 && vEDown > 0.0) { tauVsEleUp *= vEUp / vE; tauVsEleDown *= vEDown / vE; }
            const double vM = tauMuSFByIndex(idx), vMUp = tauMuSFUpByIndex(idx), vMDown = tauMuSFDownByIndex(idx);
            if (vM > 0.0 && vMUp > 0.0 && vMDown > 0.0) { tauVsMuUp *= vMUp / vM; tauVsMuDown *= vMDown / vM; }
            const double vJ = tauJetSFByIndex(idx), vJUp = tauJetSFUpByIndex(idx), vJDown = tauJetSFDownByIndex(idx);
            if (vJ > 0.0 && vJUp > 0.0 && vJDown > 0.0) { tauVsJetUp *= vJUp / vJ; tauVsJetDown *= vJDown / vJ; }
        }
    }

    result.push_back({"_eRecoUp", baseWeight * eRecoUp}); result.push_back({"_eRecoDown", baseWeight * eRecoDown});
    result.push_back({"_eIdIsoUp", baseWeight * eIdIsoUp}); result.push_back({"_eIdIsoDown", baseWeight * eIdIsoDown});
    result.push_back({"_muIdUp", baseWeight * muIdUp}); result.push_back({"_muIdDown", baseWeight * muIdDown});
    result.push_back({"_muIsoUp", baseWeight * muIsoUp}); result.push_back({"_muIsoDown", baseWeight * muIsoDown});
    result.push_back({"_tauVsEleUp", baseWeight * tauVsEleUp}); result.push_back({"_tauVsEleDown", baseWeight * tauVsEleDown});
    result.push_back({"_tauVsMuUp", baseWeight * tauVsMuUp}); result.push_back({"_tauVsMuDown", baseWeight * tauVsMuDown});
    result.push_back({"_tauVsJetUp", baseWeight * tauVsJetUp}); result.push_back({"_tauVsJetDown", baseWeight * tauVsJetDown});
    return result;
}

struct TrigSFShift { double eUp = 1.0, eDown = 1.0, muUp = 1.0, muDown = 1.0; };

inline double muonTrigPtThreshold(const std::string& year) { return year == "2017" ? 29.0 : 26.0; }

inline TrigSFShift computeTriggerSFShifts(const std::string& catstr, const std::string& year) {
    const double muTrigPtMin = muonTrigPtThreshold(year);
    const int nLep = static_cast<int>(catstr.size());
    const bool hasElectron = catstr.find('e') != std::string::npos;
    const bool hasMuon = catstr.find('m') != std::string::npos;

    TrigSFShift result;

    if (hasElectron && hasMuon) {
        double electronSF = 0.0, muonSF = 0.0;
        int eIdx = -1, mIdx = -1;
        for (int idx = 1; idx <= nLep; ++idx) {
            const char flavor = catstr[idx - 1];
            if (flavor == 't' || !trigPassedByIndex(idx)) continue;
            if (flavor == 'e' && electronSF == 0.0) { electronSF = trigSFByIndex(idx); eIdx = idx; }
            else if (flavor == 'm' && muonSF == 0.0) { muonSF = trigSFByIndex(idx); mIdx = idx; }
        }
        double nominal = 1.0;
        if (electronSF > 0.0 && muonSF > 0.0) nominal = electronSF + muonSF - electronSF * muonSF;
        else if (electronSF > 0.0) nominal = electronSF;
        else if (muonSF > 0.0) nominal = muonSF;
        result.eUp = result.eDown = result.muUp = result.muDown = nominal;

        if (eIdx > 0) {
            const double eUpSF = trigSFUpByIndex(eIdx);
            const double eDownSF = trigSFDownByIndex(eIdx);
            if (eUpSF > 0.0) result.eUp = (muonSF > 0.0) ? (eUpSF + muonSF - eUpSF * muonSF) : eUpSF;
            if (eDownSF > 0.0) result.eDown = (muonSF > 0.0) ? (eDownSF + muonSF - eDownSF * muonSF) : eDownSF;
        }
        if (mIdx > 0 && ptByIndex(mIdx) > muTrigPtMin) {
            const double mUpSF = trigSFUpByIndex(mIdx);
            const double mDownSF = trigSFDownByIndex(mIdx);
            if (mUpSF > 0.0) result.muUp = (electronSF > 0.0) ? (electronSF + mUpSF - electronSF * mUpSF) : mUpSF;
            if (mDownSF > 0.0) result.muDown = (electronSF > 0.0) ? (electronSF + mDownSF - electronSF * mDownSF) : mDownSF;
        }
        return result;
    }

    int selIdx = -1;
    double nominal = 1.0;
    if (nLep < 3) {
        if (isTrig_1 >= 1) { nominal = TrigSF_1; selIdx = 1; }
        else if (isTrig_1 == -1) { nominal = TrigSF_2; selIdx = 2; }
    }
    else {
        if (isTrig_1 >= 1 && isTrig_2 == 0) { nominal = TrigSF_1; selIdx = 1; }
        else if (isTrig_1 == -1 && isTrig_2 == 0) { nominal = TrigSF_2; selIdx = 2; }
        else if (isTrig_2 >= 1 && isTrig_1 == 0) { nominal = TrigSF_3; selIdx = 3; }
        else if (isTrig_2 == -1 && isTrig_1 == 0) { nominal = TrigSF_4; selIdx = 4; }
        else if (isTrig_1 == 2 && isTrig_2 == 2) { nominal = TrigSF_1; selIdx = 1; }
    }
    result.eUp = result.eDown = result.muUp = result.muDown = nominal;
    if (selIdx > 0 && selIdx <= nLep) {
        const char flavor = catstr[selIdx - 1];
        if (flavor == 'e') {
            const double eUpSF = trigSFUpByIndex(selIdx), eDownSF = trigSFDownByIndex(selIdx);
            if (eUpSF > 0.0) result.eUp = eUpSF;
            if (eDownSF > 0.0) result.eDown = eDownSF;
        }
        else if (flavor == 'm' && ptByIndex(selIdx) > muTrigPtMin) {
            const double mUpSF = trigSFUpByIndex(selIdx), mDownSF = trigSFDownByIndex(selIdx);
            if (mUpSF > 0.0) result.muUp = mUpSF;
            if (mDownSF > 0.0) result.muDown = mDownSF;
        }
    }
    return result;
}

inline void appendStandardWeightVariants(const std::string& catstr, bool isData, bool isTTbar, double xsecUnc,
                                          double weight, const std::string& year,
                                          std::vector<WeightSystematic>& weightVariants) {
    if (isData) return;
    appendPileupWeightVariants(weight, weightVariants);
    appendL1PrefiringWeightVariants(weight, weightVariants);
    if (xsecUnc > 0.0) {
        weightVariants.push_back({"_xsecUp", weight * (1.0 + xsecUnc)});
        weightVariants.push_back({"_xsecDown", weight * std::max(0.0, 1.0 - xsecUnc)});
    }
    for (const auto& sfVar : computeObjectSFSystematics(catstr, isData, weight)) {
        weightVariants.push_back({sfVar.suffix, sfVar.weight});
    }
    const double trigSFNominal = computeTriggerSF(catstr);
    if (trigSFNominal > 0.0) {
        const TrigSFShift trigShift = computeTriggerSFShifts(catstr, year);
        weightVariants.push_back({"_eTrigUp", weight / trigSFNominal * trigShift.eUp});
        weightVariants.push_back({"_eTrigDown", weight / trigSFNominal * trigShift.eDown});
        weightVariants.push_back({"_muTrigUp", weight / trigSFNominal * trigShift.muUp});
        weightVariants.push_back({"_muTrigDown", weight / trigSFNominal * trigShift.muDown});
    }
    if (isTTbar) {
        const double topPtNom = topPtWeight(isData, isTTbar);
        if (topPtNom > 0.0) {
            weightVariants.push_back({"_topPtUp", weight / topPtNom * topPtWeightUpVariation(isData, isTTbar)});
            weightVariants.push_back({"_topPtDown", weight / topPtNom * topPtWeightDownVariation()});
        }
    }
}
