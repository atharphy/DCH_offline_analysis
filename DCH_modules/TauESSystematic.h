#pragma once

#include <algorithm>
#include <string>

#include "ObjectAccessors.h"

inline void computeTauESRelShift(const std::string& cat, double relUp[4], double relDown[4]) {
    const int nLep = static_cast<int>(cat.size());
    for (int idx = 1; idx <= 4; ++idx) {
        relUp[idx - 1] = relDown[idx - 1] = 0.0;
        if (idx > nLep || cat[idx - 1] != 't') continue;
        const double nom = tauESByIndex(idx);
        const double up = tauESUpByIndex(idx);
        const double down = tauESDownByIndex(idx);

        if (nom <= 0.0 || up <= 0.0 || down <= 0.0) continue;
        relUp[idx - 1] = up / nom - 1.0;
        relDown[idx - 1] = down / nom - 1.0;
    }
}

inline void shiftTauESPts(const std::string& cat, const double origPt[4], const double relArr[4], double sign) {
    for (int idx = 1; idx <= static_cast<int>(cat.size()) && idx <= 4; ++idx) {
        if (cat[idx - 1] != 't') continue;
        const double factor = std::max(0.0, 1.0 + sign * relArr[idx - 1]);
        const double shifted = origPt[idx - 1] * factor;
        if (idx == 1) pt_1 = shifted;
        else if (idx == 2) pt_2 = shifted;
        else if (idx == 3) pt_3 = shifted;
        else if (idx == 4) pt_4 = shifted;
    }
}
