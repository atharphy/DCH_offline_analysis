#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>

#include "../roccor/RoccoR.cc"

#include "ObjectAccessors.h"

const bool APPLY_ROCCOR_DATA = true;
const bool APPLY_ROCCOR_MC   = true;

inline std::unique_ptr<RoccoR> loadRoccoRCorrections(const std::string& year) {
    std::string fileName;
    if (year == "2016preVFP") fileName = "roccor/RoccoR2016aUL.txt";
    else if (year == "2016postVFP") fileName = "roccor/RoccoR2016bUL.txt";
    else if (year == "2017") fileName = "roccor/RoccoR2017UL.txt";
    else if (year == "2018") fileName = "roccor/RoccoR2018UL.txt";
    else {
        std::cerr << "ERROR: no RoccoR calibration mapped for year " << year << std::endl;
        return nullptr;
    }

    try {
        std::unique_ptr<RoccoR> rc(new RoccoR(fileName));
        std::cout << "Loaded RoccoR muon momentum corrections from:\n  " << fileName << std::endl;
        return rc;
    }
    catch (const std::exception& error) {
        std::cerr << "ERROR: cannot load RoccoR calibration from:\n  " << fileName << "\n  " << error.what() << std::endl;
        return nullptr;
    }
}

inline void applyRoccoRCorrection(const std::string& cat, bool isData, const RoccoR& rc, double relErr[4] = nullptr) {
    for (int idx = 1; idx <= static_cast<int>(cat.size()) && idx <= 4; ++idx) {
        if (relErr) relErr[idx - 1] = 0.0;
        if (cat[idx - 1] != 'm') continue;

        const int q = qByIndex(idx);
        const double pt = ptByIndex(idx);
        const double eta = etaByIndex(idx);
        const double phi = phiByIndex(idx);

        double sf = 1.0;
        double err = 0.0;
        if (isData) {
            if (!APPLY_ROCCOR_DATA) continue;
            sf = rc.kScaleDT(q, pt, eta, phi);
            err = rc.kScaleDTerror(q, pt, eta, phi);
        } else {
            if (!APPLY_ROCCOR_MC) continue;
            const double genPt = truthPtByIndex(idx);
            if (genPartFlavByIndex(idx) == 1 && genPt > 0.0) {
                sf = rc.kSpreadMC(q, pt, eta, phi, genPt);
                err = rc.kSpreadMCerror(q, pt, eta, phi, genPt);
            } else {
                sf = rc.kScaleMC(q, pt, eta, phi);

                err = 0.0;
            }
        }

        if (idx == 1) pt_1 *= sf;
        else if (idx == 2) pt_2 *= sf;
        else if (idx == 3) pt_3 *= sf;
        else if (idx == 4) pt_4 *= sf;
        if (relErr) relErr[idx - 1] = err;
    }
}

inline void shiftMuonPts(const std::string& cat, const double origPt[4], const double relErr[4], double sign) {
    for (int idx = 1; idx <= static_cast<int>(cat.size()) && idx <= 4; ++idx) {
        if (cat[idx - 1] != 'm') continue;
        const double factor = std::max(0.0, 1.0 + sign * relErr[idx - 1]);
        const double shifted = origPt[idx - 1] * factor;
        if (idx == 1) pt_1 = shifted;
        else if (idx == 2) pt_2 = shifted;
        else if (idx == 3) pt_3 = shifted;
        else if (idx == 4) pt_4 = shifted;
    }
}

inline void restorePts(const double origPt[4]) {
    pt_1 = origPt[0];
    pt_2 = origPt[1];
    pt_3 = origPt[2];
    pt_4 = origPt[3];
}
