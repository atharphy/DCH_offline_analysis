#pragma once

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "CommonConfig.h"
#include "EventKinematics.h"
#include "HistUtils.h"
#include "ObjectAccessors.h"
#include "SystematicPlan.h"

inline void deleteHists(std::unordered_map<std::string, TH1D*>& histMap) {
    for (auto& item : histMap) delete item.second;
    histMap.clear();
}

struct HistBundle {
    std::unordered_map<std::string, TH1D*> mZ1, mZ2, mH1, mH2, zPt, met, metphi, LT;
    std::unordered_map<std::string, TH1D*> pt[NlepMax], eta[NlepMax], phi[NlepMax], d0[NlepMax], dZ[NlepMax], iso[NlepMax];

    void reserveAll(size_t n) {
        mZ1.reserve(n); mZ2.reserve(n); mH1.reserve(n); mH2.reserve(n);
        zPt.reserve(n); met.reserve(n); metphi.reserve(n); LT.reserve(n);
        for (int i = 0; i < NlepMax; ++i) {
            pt[i].reserve(n); eta[i].reserve(n); phi[i].reserve(n);
            d0[i].reserve(n); dZ[i].reserve(n); iso[i].reserve(n);
        }
    }
    void writeAll() {
        writeHists(mZ1); writeHists(mZ2); writeHists(mH1); writeHists(mH2);
        writeHists(zPt); writeHists(met); writeHists(metphi); writeHists(LT);
        for (int i = 0; i < NlepMax; ++i) {
            writeHists(pt[i]); writeHists(eta[i]); writeHists(phi[i]);
            writeHists(d0[i]); writeHists(dZ[i]); writeHists(iso[i]);
        }
    }
    void deleteAll() {
        deleteHists(mZ1); deleteHists(mZ2); deleteHists(mH1); deleteHists(mH2);
        deleteHists(zPt); deleteHists(met); deleteHists(metphi); deleteHists(LT);
        for (int i = 0; i < NlepMax; ++i) {
            deleteHists(pt[i]); deleteHists(eta[i]); deleteHists(phi[i]);
            deleteHists(d0[i]); deleteHists(dZ[i]); deleteHists(iso[i]);
        }
    }
};

inline void fillMetOnly(const std::string& key, double met_, double metphi_, double weight, HistBundle& h) {
    fillOne(h.met, key, met_, weight, S_met);
    fillOne(h.metphi, key, metphi_, weight, S_metphi);
}

inline void fillAllVars(const std::string& key, const std::string& catstr, const EventKinematics& kin, double weight, HistBundle& h) {
    fillOne(h.met, key, met, weight, S_met);
    fillOne(h.metphi, key, metphi, weight, S_metphi);
    fillOne(h.LT, key, kin.LT, weight, S_LT);
    if (kin.hasMZ1) fillOne(h.mZ1, key, kin.mZ1, weight, S_mZ1);
    if (kin.hasMZ2) fillOne(h.mZ2, key, kin.mZ2, weight, S_mZ2);
    if (kin.hasMH1) fillOne(h.mH1, key, kin.mH1, weight, S_mH1);
    if (kin.hasMH2) fillOne(h.mH2, key, kin.mH2, weight, S_mH2);
    if (kin.hasZPt && std::isfinite(kin.zPt)) fillOne(h.zPt, key, kin.zPt, weight, S_zPt);
    const int nLep = std::min(static_cast<int>(catstr.size()), NlepMax);
    for (int i = 0; i < nLep; ++i) {
        const int idx = i + 1;
        fillOne(h.pt[i], key, ptByIndex(idx), weight, S_pt[i]);
        fillOne(h.eta[i], key, etaByIndex(idx), weight, S_eta[i]);
        fillOne(h.phi[i], key, phiByIndex(idx), weight, S_phi[i]);
        fillOne(h.d0[i], key, d0ByIndex(idx), weight, S_d0[i]);
        fillOne(h.dZ[i], key, dZByIndex(idx), weight, S_dZ[i]);
        fillOne(h.iso[i], key, isoByIndex(idx), weight, S_iso[i]);
    }
}

inline void fillEventHistograms(const std::vector<std::string>& histKeys, const std::string& catstr, const EventKinematics& kin,
                                 double weight, const std::vector<WeightSystematic>& weightVariants,
                                 const std::vector<MetSystematic>& metVariants, HistBundle& h) {
    for (const auto& key : histKeys) fillAllVars(key, catstr, kin, weight, h);
    for (const auto& variant : weightVariants)
        for (const auto& key : histKeys) fillAllVars(key + variant.suffix, catstr, kin, variant.weight, h);
    for (const auto& metVar : metVariants) {
        double shiftedMet, shiftedMetPhi;
        cartesianToPolar(metVar.px, metVar.py, shiftedMet, shiftedMetPhi);
        for (const auto& key : histKeys) fillMetOnly(key + metVar.suffix, shiftedMet, shiftedMetPhi, weight, h);
    }
}
