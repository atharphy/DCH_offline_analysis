#pragma once

#include <string>
#include <vector>

#include "ObjectAccessors.h"
#include "SystematicPlan.h"

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
