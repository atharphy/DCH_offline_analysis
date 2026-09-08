#pragma once

#include <string>
#include <utility>
#include <vector>

#include "EventSelection.h"
#include "ObjectAccessors.h"

struct EventKinematics {
    double LT = 0.0;
    bool hasMZ1 = false, hasMZ2 = false, hasMH1 = false, hasMH2 = false, hasZPt = false;
    double mZ1 = 0.0, mZ2 = 0.0, mH1 = 0.0, mH2 = 0.0, zPt = 0.0;
};

inline EventKinematics buildEventKinematics(const std::string& catstr,
                                             const std::vector<std::pair<int,int>>& OS_pair,
                                             const std::vector<std::pair<int,int>>& SS_pair) {
    EventKinematics result;
    result.LT = getLT(catstr);
    if (!OS_pair.empty()) {
        const TLorentzVector zCandidate = LepV(OS_pair[0].first) + LepV(OS_pair[0].second);
        result.hasMZ1 = true;
        result.mZ1 = zCandidate.M();
        result.hasZPt = true;
        result.zPt = zCandidate.Pt();
    }
    if (OS_pair.size() > 1) {
        result.hasMZ2 = true;
        result.mZ2 = (LepV(OS_pair[1].first) + LepV(OS_pair[1].second)).M();
    }
    if (!SS_pair.empty()) {
        result.hasMH1 = true;
        result.mH1 = (LepV(SS_pair[0].first) + LepV(SS_pair[0].second)).M();
    }
    if (SS_pair.size() > 1) {
        result.hasMH2 = true;
        result.mH2 = (LepV(SS_pair[1].first) + LepV(SS_pair[1].second)).M();
    }
    return result;
}

inline std::vector<std::string> collectHistKeys(const std::string& catstr, const EventKinematics& kin,
                                                 const std::vector<std::pair<int,int>>& OS_pair) {
    std::vector<std::string> keys;
    keys.push_back(catstr);
    if (catstr.size() == 2) {
        const std::string zFlag = dilepZFlag(catstr);
        const bool isOS = q_1 * q_2 < 0;
        const bool isSS = q_1 * q_2 > 0;
        std::string dileptonKey;
        if (isOS && zFlag == "Zwin") dileptonKey = catstr + "_OS_Zwin";
        else if (isSS && zFlag == "Zwin") dileptonKey = catstr + "_SS_Zwin";
        else if (isOS && zFlag == "Zveto") dileptonKey = catstr + "_OS_Zveto";
        else if (isSS && zFlag == "Zveto") dileptonKey = catstr + "_SS_Zveto";
        if (!dileptonKey.empty()) keys.push_back(dileptonKey);
    }
    const std::string region = classifyTauRegion(catstr, kin.LT, OS_pair);
    keys.push_back(catstr + "_" + region);
    return keys;
}
