#pragma once

#include <cmath>
#include <string>

#include "ObjectAccessors.h"

inline bool passTightEvent(const std::string& catstr) {
    if (catstr.size() < 2 || catstr.size() > 4) return false;
    for (int idx = 1; idx <= static_cast<int>(catstr.size()); ++idx) {
        if (!isValidFlavor(catstr[idx - 1])) return false;
        if (ptByIndex(idx) <= 0.0) return false;
    }
    return true;
}

inline bool passChargeTopology(const std::string& catstr) {
    if (catstr.size() == 4) return q_1 + q_2 + q_3 + q_4 == 0;
    if (catstr.size() == 3) return std::abs(q_1 + q_2 + q_3) != 3;
    return true;
}

inline bool hasDuplicateTightObjects(const std::string& catstr) {
    const int nLep = static_cast<int>(catstr.size());
    for (int i = 1; i <= nLep; ++i) {
        for (int j = i + 1; j <= nLep; ++j) {
            if (getDR(etaByIndex(i), phiByIndex(i), etaByIndex(j), phiByIndex(j)) <= 0.4) return true;
        }
    }
    return false;
}

inline double getLT(const std::string& catstr) {
    double LT = 0.0;
    for (int idx = 1; idx <= static_cast<int>(catstr.size()); ++idx) LT += ptByIndex(idx);
    return LT;
}
