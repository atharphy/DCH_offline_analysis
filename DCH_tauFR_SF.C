// DCH_tauFR_SF.C
//
// Jet->tau_h fake background via the SF-correction method (CMS-AN-24-174
// Sec. 8.2 / CMS-AN-19-111 Appendix J style), instead of the tight-to-loose
// promotion used by DCH_tauFR.C. Every event keeps exactly the tau legs it
// was already reconstructed with (no loose-collection candidates, no 2^N
// promotion): a reco-tight tau that fails MC-truth match to a genuine
// hadronic tau (genPartFlav != 5) has its event weight multiplied by
// SF(pT,eta) = f_data/f_DY_MC (built by
// fake_rates/tau_fake_rates/BuildTauFakeRateSF.C); a genuine tau is left
// untouched. For events with 2 or 3 fake tau legs, each leg's SF is
// multiplied in independently (same factorization assumption the fake-
// factor method already makes everywhere else).
//
// Scope note: e->tau fakes (genPartFlav==1) are deliberately left
// untouched (SF=1) in this driver -- there is no independent data-driven
// measurement to build a data/MC SF from for that channel (the existing
// e->tau rate is itself MC-truth only), so reusing it here would conflate
// two different quantities. This is an open item, not an oversight.
//
// Compile:
//   root -l -b -q 'DCH_tauFR_SF.C++("2018",0,-1,1)'
//
// Run one file:
//   root -l -b -q 'DCH_tauFR_SF.C+("2018",0,1,1)'
//
// Run many files in one ROOT job with process parallelism:
//   root -l -b -q 'DCH_tauFR_SF.C+("2018",0,-1,8)'

#if !defined(__CLING__)
#pragma GCC optimize("O3,unroll-loops")
#endif

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Compression.h"
#include "ROOT/TProcessExecutor.hxx"
#include "correction.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#if defined(__CLING__)
R__LOAD_LIBRARY(libcorrectionlib.so)
R__LOAD_LIBRARY(libHTT-utilitiesRecoilCorrections.so)
#endif

#include "filemap/FileMap.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "include/MyBranch.C"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#pragma GCC diagnostic ignored "-Wunused-value"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wuninitialized"
#include "include/Kinematics.C"
#pragma GCC diagnostic pop

#include "include/Xsections.C"

#include "DCH_modules/CommonConfig.h"
#include "DCH_modules/HistUtils.h"
#include "DCH_modules/ObjectAccessors.h"
#include "DCH_modules/METCorrections.h"
#include "DCH_modules/RecoilCorrections.h"
#include "DCH_modules/RoccoRCorrections.h"
#include "DCH_modules/ZPtReweight.h"
#include "DCH_modules/GenZMatching.h"
#include "DCH_modules/GenBosonMomentum.h"
#include "DCH_modules/FileJob.h"
#include "DCH_modules/TauFakeSFReader.h"
#include "DCH_modules/EventWeights.h"
#include "DCH_modules/PairBuilder.h"
#include "DCH_modules/MEtSysWrapper.h"
#include "DCH_modules/SystematicPlan.h"
#include "DCH_modules/FlatXsecSystematic.h"
#include "DCH_modules/TauESSystematic.h"

using std::string;
using std::vector;
using std::unordered_map;

const bool SKIP_DATA = false;
const bool APPLY_TAU_FAKE_RATE = true;

const string FR_SF_FILE = "Dependencies/fake_rates/tau_fake_rates/DY_tau_fake_rate_SF_2D.root";

const vector<string> tauRegions = {
    "DYCR_0tau", "DYCR_1tau", "DYveto_0tau", "DYveto_1tau",
    "CR_0tau", "CR_1tau", "CR_2tau", "CR_3tau", "CR_3lep0tau", "CR_3lep1tau", "CR_3lep2tau",
    "VR_0tau", "VR_1tau", "VR_2tau", "VR_3tau", "VR_3lep0tau", "VR_3lep1tau", "VR_3lep2tau",
    "SR_0tau", "SR_1tau", "SR_2tau", "SR_3tau", "SR_3lep0tau", "SR_3lep1tau", "SR_3lep2tau"
};

struct Obj {
    char flav = 'x';
    double pt = -99;
    double eta = -99;
    double phi = -99;
    double mass = 0;
    int q = 0;
    double d0 = -99;
    double dZ = -99;
    double iso = -99;
    bool fromTight = false;
    int tightIdx = -1;
    int objectOrder = 999;
};

bool isDilep(const string& ch) { return ch.size() == 2; }

bool getRecoZPt(const string& catstr, double& zPt) {
    vector<std::pair<int,int>> Z_pair, Zv_pair, Ztt_pair, SS_pair, OSDF_pair;
    processPairs(catstr.c_str(), Z_pair, Zv_pair, Ztt_pair, SS_pair, OSDF_pair);
    auto pairPtGreater = [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
        return (LepV(a.first) + LepV(a.second)).Pt() > (LepV(b.first) + LepV(b.second)).Pt();
    };
    std::sort(Z_pair.begin(), Z_pair.end(), pairPtGreater);
    std::sort(Zv_pair.begin(), Zv_pair.end(), pairPtGreater);
    std::sort(Ztt_pair.begin(), Ztt_pair.end(), pairPtGreater);
    Zv_pair = removeOverlap(Zv_pair, static_cast<int>(catstr.size()));
    vector<std::pair<int,int>> OSSF_pair;
    OSSF_pair.insert(OSSF_pair.end(), Z_pair.begin(), Z_pair.end());
    OSSF_pair.insert(OSSF_pair.end(), Zv_pair.begin(), Zv_pair.end());
    OSSF_pair = removeOverlap(OSSF_pair, static_cast<int>(catstr.size()));
    vector<std::pair<int,int>> OS_pair;
    OS_pair.insert(OS_pair.end(), OSSF_pair.begin(), OSSF_pair.end());
    OS_pair.insert(OS_pair.end(), Ztt_pair.begin(), Ztt_pair.end());
    OS_pair.insert(OS_pair.end(), OSDF_pair.begin(), OSDF_pair.end());
    OS_pair = removeOverlap(OS_pair, static_cast<int>(catstr.size()));
    if (OS_pair.empty()) return false;
    const int first = OS_pair[0].first;
    const int second = OS_pair[0].second;
    if (first < 1 || second < 1 || first > static_cast<int>(catstr.size()) || second > static_cast<int>(catstr.size())) return false;
    const char firstFlavor = catstr[first - 1];
    const char secondFlavor = catstr[second - 1];
    if (firstFlavor != secondFlavor) return false;
    if (firstFlavor != 'e' && firstFlavor != 'm') return false;
    if (qByIndex(first) * qByIndex(second) >= 0) return false;
    const TLorentzVector candidate = LepV(first) + LepV(second);
    const double mass = candidate.M();
    if (!std::isfinite(mass) || std::fabs(mass - 91.2) > 10.0) return false;
    zPt = candidate.Pt();
    return std::isfinite(zPt);
}

void resetSlots() { pt_1=pt_2=pt_3=pt_4=-99; eta_1=eta_2=eta_3=eta_4=-99; phi_1=phi_2=phi_3=phi_4=-99; m_1=m_2=m_3=m_4=0; q_1=q_2=q_3=q_4=0; d0_1=d0_2=d0_3=d0_4=-99; dZ_1=dZ_2=dZ_3=dZ_4=-99; iso_1=iso_2=iso_3=iso_4=-99; }

void setSlot(int idx, const Obj& o) {
    if (idx == 1) { pt_1=o.pt; eta_1=o.eta; phi_1=o.phi; m_1=o.mass; q_1=o.q; d0_1=o.d0; dZ_1=o.dZ; iso_1=o.iso; }
    else if (idx == 2) { pt_2=o.pt; eta_2=o.eta; phi_2=o.phi; m_2=o.mass; q_2=o.q; d0_2=o.d0; dZ_2=o.dZ; iso_2=o.iso; }
    else if (idx == 3) { pt_3=o.pt; eta_3=o.eta; phi_3=o.phi; m_3=o.mass; q_3=o.q; d0_3=o.d0; dZ_3=o.dZ; iso_3=o.iso; }
    else if (idx == 4) { pt_4=o.pt; eta_4=o.eta; phi_4=o.phi; m_4=o.mass; q_4=o.q; d0_4=o.d0; dZ_4=o.dZ; iso_4=o.iso; }
}

string catFromObjects(const vector<Obj>& objs) { string s; for (const auto& o : objs) s += o.flav; return s; }

void loadObjectsIntoGlobals(const vector<Obj>& objs) { resetSlots(); for (int i = 0; i < (int)objs.size() && i < 4; ++i) setSlot(i + 1, objs[i]); }

bool passChargeTopologyConfig(const vector<Obj>& objs) {
    if (objs.size() == 4) {
        int qsum = 0;
        for (const auto& o : objs) qsum += o.q;
        return qsum == 0;
    }
    if (objs.size() == 3) {
        int qsum = 0;
        for (const auto& o : objs) qsum += o.q;
        return std::abs(qsum) != 3;
    }
    return true;
}

bool hasDuplicateObjects(const vector<Obj>& objs) {
    for (size_t first = 0; first < objs.size(); ++first) {
        for (size_t second = first + 1; second < objs.size(); ++second) {
            if (getDR(objs[first].eta, objs[first].phi, objs[second].eta, objs[second].phi) <= 0.4) return true;
        }
    }
    return false;
}

double getLTFromObjects(const vector<Obj>& objs) { double lt = 0.0; for (const auto& o : objs) lt += o.pt; return lt; }

double computeEventWeightWithoutObjectSF(const string& originalCat, bool isData, double xsw) {
    if (isData) return 1.0;
    return computeBaseGenWeight(xsw) * computeTriggerSF(originalCat);
}

double computeConfigurationObjectSF(const vector<Obj>& objs) {
    double objectSF = 1.0;
    for (const auto& obj : objs) {
        if (!obj.fromTight || obj.tightIdx < 1 || obj.tightIdx > 4) continue;
        if (obj.flav == 'e' || obj.flav == 'm') {
            objectSF *= idSFByIndex(obj.tightIdx);
            objectSF *= isoSFByIndex(obj.tightIdx);
        }
        else if (obj.flav == 't' && genPartFlavByIndex(obj.tightIdx) == 5) {
            objectSF *= tauEleSFByIndex(obj.tightIdx);
            objectSF *= tauMuSFByIndex(obj.tightIdx);
            objectSF *= tauJetSFByIndex(obj.tightIdx);
        }
    }
    return objectSF;
}

struct ObjectSFSystEntry { string suffix; double objectSF; };

// Reads the skim's own precomputed Up/Down SF branches directly (see
// DCH_modules/ObjectAccessors.h) instead of recomputing a relative
// uncertainty via correctionlib at analysis time. `year` is unused now but
// kept in the signature to avoid touching the call site.
vector<ObjectSFSystEntry> computeConfigurationObjectSFSystematics(const vector<Obj>& objs, const string& year) {
    (void)year;
    double eRecoUp = 1.0, eRecoDown = 1.0, eIdIsoUp = 1.0, eIdIsoDown = 1.0;
    double muIdUp = 1.0, muIdDown = 1.0, muIsoUp = 1.0, muIsoDown = 1.0;
    double tauVsEleUp = 1.0, tauVsEleDown = 1.0, tauVsMuUp = 1.0, tauVsMuDown = 1.0, tauVsJetUp = 1.0, tauVsJetDown = 1.0;

    for (const auto& obj : objs) {
        if (!obj.fromTight || obj.tightIdx < 1 || obj.tightIdx > 4) continue;
        const int idx = obj.tightIdx;
        if (obj.flav == 'e' && obj.pt > 20.0) {
            const double idNom = idSFByIndex(idx), idUp = idSFUpByIndex(idx), idDown = idSFDownByIndex(idx);
            if (idNom > 0.0 && idUp > 0.0 && idDown > 0.0) { eRecoUp *= idUp / idNom; eRecoDown *= idDown / idNom; }
            const double isoNom = isoSFByIndex(idx), isoUp = isoSFUpByIndex(idx), isoDown = isoSFDownByIndex(idx);
            if (isoNom > 0.0 && isoUp > 0.0 && isoDown > 0.0) { eIdIsoUp *= isoUp / isoNom; eIdIsoDown *= isoDown / isoNom; }
        }
        else if (obj.flav == 'm' && obj.pt > 15.0) {
            const double idNom = idSFByIndex(idx), idUp = idSFUpByIndex(idx), idDown = idSFDownByIndex(idx);
            if (idNom > 0.0 && idUp > 0.0 && idDown > 0.0) { muIdUp *= idUp / idNom; muIdDown *= idDown / idNom; }
            const double isoNom = isoSFByIndex(idx), isoUp = isoSFUpByIndex(idx), isoDown = isoSFDownByIndex(idx);
            if (isoNom > 0.0 && isoUp > 0.0 && isoDown > 0.0) { muIsoUp *= isoUp / isoNom; muIsoDown *= isoDown / isoNom; }
        }
        else if (obj.flav == 't' && genPartFlavByIndex(idx) == 5) {
            const double vE = tauEleSFByIndex(idx), vEUp = tauEleSFUpByIndex(idx), vEDown = tauEleSFDownByIndex(idx);
            if (vE > 0.0 && vEUp > 0.0 && vEDown > 0.0) { tauVsEleUp *= vEUp / vE; tauVsEleDown *= vEDown / vE; }
            const double vM = tauMuSFByIndex(idx), vMUp = tauMuSFUpByIndex(idx), vMDown = tauMuSFDownByIndex(idx);
            if (vM > 0.0 && vMUp > 0.0 && vMDown > 0.0) { tauVsMuUp *= vMUp / vM; tauVsMuDown *= vMDown / vM; }
            const double vJ = tauJetSFByIndex(idx), vJUp = tauJetSFUpByIndex(idx), vJDown = tauJetSFDownByIndex(idx);
            if (vJ > 0.0 && vJUp > 0.0 && vJDown > 0.0) { tauVsJetUp *= vJUp / vJ; tauVsJetDown *= vJDown / vJ; }
        }
    }

    const double nominal = computeConfigurationObjectSF(objs);
    return {
        {"_eRecoUp", nominal * eRecoUp}, {"_eRecoDown", nominal * eRecoDown},
        {"_eIdIsoUp", nominal * eIdIsoUp}, {"_eIdIsoDown", nominal * eIdIsoDown},
        {"_muIdUp", nominal * muIdUp}, {"_muIdDown", nominal * muIdDown},
        {"_muIsoUp", nominal * muIsoUp}, {"_muIsoDown", nominal * muIsoDown},
        {"_tauVsEleUp", nominal * tauVsEleUp}, {"_tauVsEleDown", nominal * tauVsEleDown},
        {"_tauVsMuUp", nominal * tauVsMuUp}, {"_tauVsMuDown", nominal * tauVsMuDown},
        {"_tauVsJetUp", nominal * tauVsJetUp}, {"_tauVsJetDown", nominal * tauVsJetDown},
    };
}

struct TrigSFShift { double eUp = 1.0, eDown = 1.0, muUp = 1.0, muDown = 1.0; };

double muonTrigPtThreshold(const string& year) {
    if (year == "2017") return 29.0;
    return 26.0;
}

// Reads the skim's precomputed TrigSF_Up/Down branches directly (absolute
// shifted values) instead of recomputing a relative uncertainty via
// correctionlib/TGraphAsymmErrors. `year` unused now, kept for call-site
// compatibility.
TrigSFShift computeTriggerSFShifts(const string& catstr, const string& year) {
    const double muTrigPtMin = muonTrigPtThreshold(year);

    const int nLep = static_cast<int>(catstr.size());
    const bool hasElectron = catstr.find('e') != string::npos;
    const bool hasMuon = catstr.find('m') != string::npos;

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

void fillAllVars(const string& key, const vector<Obj>& objs, double weight, unordered_map<string,TH1D*>& h_mZ1, unordered_map<string,TH1D*>& h_mZ2, unordered_map<string,TH1D*>& h_mH1, unordered_map<string,TH1D*>& h_mH2, unordered_map<string,TH1D*>& h_zPt, unordered_map<string,TH1D*>& h_met, unordered_map<string,TH1D*>& h_metphi, unordered_map<string,TH1D*>& h_LT, unordered_map<string,TH1D*> h_pt[NlepMax], unordered_map<string,TH1D*> h_eta[NlepMax], unordered_map<string,TH1D*> h_phi[NlepMax], unordered_map<string,TH1D*> h_d0[NlepMax], unordered_map<string,TH1D*> h_dZ[NlepMax], unordered_map<string,TH1D*> h_iso[NlepMax])
{
    loadObjectsIntoGlobals(objs);
    string catstr = catFromObjects(objs);
    vector<std::pair<int,int>> OS_pair, SS_pair;
    buildPairs(catstr, OS_pair, SS_pair);
    fillOne(h_met, key, met, weight, S_met);
    fillOne(h_metphi, key, metphi, weight, S_metphi);
    fillOne(h_LT, key, getLTFromObjects(objs), weight, S_LT);
    if (!OS_pair.empty()) {
        const TLorentzVector zCandidate = LepV(OS_pair[0].first) + LepV(OS_pair[0].second);
        fillOne(h_mZ1, key, (LepV(OS_pair[0].first) + LepV(OS_pair[0].second)).M(), weight, S_mZ1);
        fillOne(h_zPt, key, zCandidate.Pt(), weight, S_zPt);
    }
    if (OS_pair.size() > 1) fillOne(h_mZ2, key, (LepV(OS_pair[1].first) + LepV(OS_pair[1].second)).M(), weight, S_mZ2);
    if (!SS_pair.empty()) fillOne(h_mH1, key, (LepV(SS_pair[0].first) + LepV(SS_pair[0].second)).M(), weight, S_mH1);
    if (SS_pair.size() > 1) fillOne(h_mH2, key, (LepV(SS_pair[1].first) + LepV(SS_pair[1].second)).M(), weight, S_mH2);
    for (int i = 0; i < (int)objs.size() && i < NlepMax; ++i) {
        fillOne(h_pt[i], key, objs[i].pt, weight, S_pt[i]);
        fillOne(h_eta[i], key, objs[i].eta, weight, S_eta[i]);
        fillOne(h_phi[i], key, objs[i].phi, weight, S_phi[i]);
        fillOne(h_d0[i], key, objs[i].d0, weight, S_d0[i]);
        fillOne(h_dZ[i], key, objs[i].dZ, weight, S_dZ[i]);
        fillOne(h_iso[i], key, objs[i].iso, weight, S_iso[i]);
    }
}

vector<string> resolveHistKeys(const vector<Obj>& objs, const std::unordered_set<string>& finalStateSet) {
    vector<string> keys;
    if (objs.size() < 2 || objs.size() > 4) return keys;
    if (!passChargeTopologyConfig(objs)) return keys;
    if (hasDuplicateObjects(objs)) return keys;
    string catstr = catFromObjects(objs);
    if (!finalStateSet.count(catstr)) return keys;
    loadObjectsIntoGlobals(objs);
    vector<std::pair<int,int>> Z_pair, Zv_pair, Ztt_pair, SS_pair, OSDF_pair;
    processPairs(catstr.c_str(), Z_pair, Zv_pair, Ztt_pair, SS_pair, OSDF_pair);
    Zv_pair = removeOverlap(Zv_pair, catstr.size());
    vector<std::pair<int,int>> OSSF_pair;
    OSSF_pair.insert(OSSF_pair.end(), Z_pair.begin(), Z_pair.end());
    OSSF_pair.insert(OSSF_pair.end(), Zv_pair.begin(), Zv_pair.end());
    OSSF_pair = removeOverlap(OSSF_pair, catstr.size());
    vector<std::pair<int,int>> OS_pair;
    OS_pair.insert(OS_pair.end(), OSSF_pair.begin(), OSSF_pair.end());
    OS_pair.insert(OS_pair.end(), Ztt_pair.begin(), Ztt_pair.end());
    OS_pair.insert(OS_pair.end(), OSDF_pair.begin(), OSDF_pair.end());
    OS_pair = removeOverlap(OS_pair, catstr.size());
    keys.push_back(catstr);
    if (catstr.size() == 2) {
        string zflag = dilepZFlag(catstr);
        bool isOS = q_1 * q_2 < 0;
        bool isSS = q_1 * q_2 > 0;
        if      (isOS && zflag == "Zwin")   keys.push_back(catstr + "_OS_Zwin");
        else if (isSS && zflag == "Zwin")   keys.push_back(catstr + "_SS_Zwin");
        else if (isOS && zflag == "Zveto")  keys.push_back(catstr + "_OS_Zveto");
        else if (isSS && zflag == "Zveto")  keys.push_back(catstr + "_SS_Zveto");
    }
    string region = classifyTauRegion(catstr, getLTFromObjects(objs), OS_pair);
    keys.push_back(catstr + "_" + region);
    return keys;
}

void fillConfiguration(const vector<Obj>& objs, double weight, const std::unordered_set<string>& finalStateSet, unordered_map<string,TH1D*>& h_mZ1, unordered_map<string,TH1D*>& h_mZ2, unordered_map<string,TH1D*>& h_mH1, unordered_map<string,TH1D*>& h_mH2, unordered_map<string,TH1D*>& h_zPt, unordered_map<string,TH1D*>& h_met, unordered_map<string,TH1D*>& h_metphi, unordered_map<string,TH1D*>& h_LT, unordered_map<string,TH1D*> h_pt[NlepMax], unordered_map<string,TH1D*> h_eta[NlepMax], unordered_map<string,TH1D*> h_phi[NlepMax], unordered_map<string,TH1D*> h_d0[NlepMax], unordered_map<string,TH1D*> h_dZ[NlepMax], unordered_map<string,TH1D*> h_iso[NlepMax], const string& keySuffix = "")
{
    for (const auto& key : resolveHistKeys(objs, finalStateSet))
        fillAllVars(key + keySuffix, objs, weight, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
}

void fillConfigurationMetOnly(const vector<Obj>& objs, double met_, double metphi_, double weight, const std::unordered_set<string>& finalStateSet, unordered_map<string,TH1D*>& h_met, unordered_map<string,TH1D*>& h_metphi, const string& keySuffix)
{
    for (const auto& key : resolveHistKeys(objs, finalStateSet)) {
        fillOne(h_met, key + keySuffix, met_, weight, S_met);
        fillOne(h_metphi, key + keySuffix, metphi_, weight, S_metphi);
    }
}

void enableBranchesTauFRSF(TTree* t, bool isData) {
    t->SetBranchStatus("*", 0);
    const char* brs[] = {
        "cat", "run", "met", "metphi", "nPV", "nPVGood", "njets",
        "pt_1", "pt_2", "pt_3", "pt_4",
        "eta_1", "eta_2", "eta_3", "eta_4",
        "phi_1", "phi_2", "phi_3", "phi_4",
        "m_1", "m_2", "m_3", "m_4",
        "q_1", "q_2", "q_3", "q_4",
        "d0_1", "d0_2", "d0_3", "d0_4",
        "dZ_1", "dZ_2", "dZ_3", "dZ_4",
        "iso_1", "iso_2", "iso_3", "iso_4",
        "genPartFlav_1", "genPartFlav_2", "genPartFlav_3", "genPartFlav_4",
        "pt_1_tr", "pt_2_tr", "pt_3_tr", "pt_4_tr",
        "phi_1_tr", "phi_2_tr", "phi_3_tr", "phi_4_tr",
        "IDSF_1", "IDSF_2", "IDSF_3", "IDSF_4",
        "ISOSF_1", "ISOSF_2", "ISOSF_3", "ISOSF_4",
        "TrigSF_1", "TrigSF_2", "TrigSF_3", "TrigSF_4",
        "isTrig_1", "isTrig_2",
        "TauVsEleIDSF_1", "TauVsEleIDSF_2", "TauVsEleIDSF_3", "TauVsEleIDSF_4",
        "TauVsMuIDSF_1", "TauVsMuIDSF_2", "TauVsMuIDSF_3", "TauVsMuIDSF_4",
        "TauVsJetIDSF_1", "TauVsJetIDSF_2", "TauVsJetIDSF_3", "TauVsJetIDSF_4",
        "TauES_1", "TauES_2", "TauES_3", "TauES_4",
        "IDSF_Up_1", "IDSF_Up_2", "IDSF_Up_3", "IDSF_Up_4",
        "IDSF_Down_1", "IDSF_Down_2", "IDSF_Down_3", "IDSF_Down_4",
        "ISOSF_Up_1", "ISOSF_Up_2", "ISOSF_Up_3", "ISOSF_Up_4",
        "ISOSF_Down_1", "ISOSF_Down_2", "ISOSF_Down_3", "ISOSF_Down_4",
        "TrigSF_Up_1", "TrigSF_Up_2", "TrigSF_Up_3", "TrigSF_Up_4",
        "TrigSF_Down_1", "TrigSF_Down_2", "TrigSF_Down_3", "TrigSF_Down_4",
        "TauVsEleIDSF_Up_1", "TauVsEleIDSF_Up_2", "TauVsEleIDSF_Up_3", "TauVsEleIDSF_Up_4",
        "TauVsEleIDSF_Down_1", "TauVsEleIDSF_Down_2", "TauVsEleIDSF_Down_3", "TauVsEleIDSF_Down_4",
        "TauVsMuIDSF_Up_1", "TauVsMuIDSF_Up_2", "TauVsMuIDSF_Up_3", "TauVsMuIDSF_Up_4",
        "TauVsMuIDSF_Down_1", "TauVsMuIDSF_Down_2", "TauVsMuIDSF_Down_3", "TauVsMuIDSF_Down_4",
        "TauVsJetIDSF_Up_1", "TauVsJetIDSF_Up_2", "TauVsJetIDSF_Up_3", "TauVsJetIDSF_Up_4",
        "TauVsJetIDSF_Down_1", "TauVsJetIDSF_Down_2", "TauVsJetIDSF_Down_3", "TauVsJetIDSF_Down_4",
        "TauES_Up_1", "TauES_Up_2", "TauES_Up_3", "TauES_Up_4",
        "TauES_Down_1", "TauES_Down_2", "TauES_Down_3", "TauES_Down_4",
        "Generator_weight", "brWeight", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Down",
        "weightPUtruejson", "weightPUtruejson_up", "weightPUtruejson_down",
        "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",
        "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_status", "GenPart_statusFlags",
        "GenVisTau_genPartIdxMother", "GenVisTau_eta", "GenVisTau_phi",
    };
    for (auto br : brs) if (t->GetBranch(br)) t->SetBranchStatus(br, 1);
}

Obj makeTightObj(char flav, int idx) {
    Obj o;
    o.flav = flav;
    o.pt = ptByIndex(idx);
    o.eta = etaByIndex(idx);
    o.phi = phiByIndex(idx);
    o.mass = massByIndex(idx);
    o.q = qByIndex(idx);
    o.d0 = d0ByIndex(idx);
    o.dZ = dZByIndex(idx);
    o.iso = isoByIndex(idx);
    o.fromTight = true;
    o.tightIdx = idx;
    o.objectOrder = idx;
    return o;
}

struct TauFakeSFResult { double sf = 1.0; double relErr = 0.0; bool hasFakeTau = false; };

// Applies SF(pT,eta) = f_data/f_DY_MC to every reco-tight tau leg that fails
// MC-truth match to a genuine hadronic tau (genPartFlav != 5), factorizing
// across however many such legs the event has (1, 2 or 3). genPartFlav==1
// (e->tau) legs are deliberately skipped -- see file header.
TauFakeSFResult computeTauFakeSF(const vector<Obj>& objs, const TauFakeSFReader& sfReader) {
    TauFakeSFResult result;
    double relVarSum = 0.0;
    for (const auto& o : objs) {
        if (o.flav != 't' || !o.fromTight || o.tightIdx < 1 || o.tightIdx > 4) continue;
        const int gpf = genPartFlavByIndex(o.tightIdx);
        if (gpf == 5) continue;
        if (gpf == 1) continue;
        double sf = 1.0, err = 0.0;
        sfReader.get(o.pt, o.eta, sf, err);
        result.sf *= sf;
        result.hasFakeTau = true;
        if (sf > 1e-6) { const double rel = err / sf; relVarSum += rel * rel; }
    }
    result.relErr = std::sqrt(relVarSum);
    return result;
}

int ProcessTauFRSFFile(FileJob job) {
    std::cout << "\n[file " << job.index + 1 << "/" << job.total << "] " << job.fileName << std::endl;
    TFile* fin = TFile::Open(job.fileName.c_str(), "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "  [skip] cannot open file" << std::endl; if (fin) fin->Close(); return 1; }
    TTree* tree = (TTree*)fin->Get("Events");
    if (!tree) { fin->Close(); return 2; }
    string baseName = gSystem->BaseName(job.fileName.c_str());
    bool isData = (XSec(baseName) == 1);
    const bool isDY = job.process == "DY" || job.process == "DY10_50";
    const bool isWJ = job.process == "WJ";
    const double xsecUnc = getFlatXsecUncertainty(baseName);
    if (SKIP_DATA && isData) { std::cout << "  [skip data]" << std::endl; fin->Close(); return 0; }
    std::unique_ptr<TH1D> zPtWeights;
    if (APPLY_ZPT_REWEIGHTING && isDY) {
        zPtWeights.reset(loadZPtWeights(job.year));
        if (!zPtWeights) { std::cerr << "  [skip] Z pT correction could not be loaded for " << job.year << std::endl; fin->Close(); return 5; }
    }
    std::unique_ptr<RecoilCorrector> recoilCorrector;
    const bool useRecoilCorrection = !isData && ((APPLY_DY_RECOIL_CORRECTION && isDY) || (APPLY_WJ_RECOIL_CORRECTION && isWJ));
    if (useRecoilCorrection) {
        recoilCorrector = loadRecoilCorrector(job.year);
        if (!recoilCorrector) { fin->Close(); return 6; }
    }
    std::unique_ptr<OfficialMETCorrections> metCorrections;
    if (APPLY_OFFICIAL_MET_CORRECTION) {
        metCorrections = loadOfficialMETCorrections(job.year, isData);
        if (!metCorrections) { fin->Close(); return 4; }
    }
    std::unique_ptr<RoccoR> roccor = loadRoccoRCorrections(job.year);
    if (!roccor) { fin->Close(); return 7; }
    enableBranchesTauFRSF(tree, isData);
    MyBranch(tree);
    if (!isData && !tree->GetBranch("IDSF_Up_1")) {
        std::cerr << "  [skip] MC file lacks precomputed SF Up/Down branches: " << job.fileName << std::endl;
        fin->Close();
        return 8;
    }
    tree->SetCacheSize(200 * 1024 * 1024);
    tree->AddBranchToCache("*", kTRUE);
    std::shared_ptr<TauFakeSFReader> sfReader;
    if (!isData && APPLY_TAU_FAKE_RATE) sfReader = getCachedTauFakeSFReader(FR_SF_FILE, job.year);
    std::shared_ptr<MEtSys> metSys;
    if (useRecoilCorrection) {
        metSys = getCachedMEtSys(job.year);
        if (!metSys) std::cerr << "  [warn] MEtSys payload not found for " << job.year << "; recoil systematic variations will be skipped" << std::endl;
    }
    TH1D* hNWEvts = (TH1D*)fin->Get("hNWEvts");
    if (!hNWEvts) hNWEvts = (TH1D*)fin->Get("hNEvts");
    double denom = hNWEvts ? hNWEvts->Integral() : 0.0;
    double xsw = (!isData && denom > 0.0) ? job.lumi * XSec(baseName) / denom : 1.0;
    Long64_t nEnt = tree->GetEntriesFast();
    TString outName = Form("%s/hist_%s", job.outDir.c_str(), baseName.c_str());
    TFile* fout = new TFile(outName, "RECREATE");
    fout->SetCompressionSettings(ROOT::CompressionSettings(ROOT::kZSTD, 1));
    if (hNWEvts) hNWEvts->Write();
    unordered_map<string,TH1D*> h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT;
    unordered_map<string,TH1D*> h_pt[NlepMax], h_eta[NlepMax], h_phi[NlepMax], h_d0[NlepMax], h_dZ[NlepMax], h_iso[NlepMax];
    std::unordered_set<string> finalStateSet(finalStates.begin(), finalStates.end());
    std::cout << "  entries in tree: " << nEnt << std::endl;
    Long64_t nPassInput = 0;
    Long64_t nZPtReweighted = 0;
    Long64_t nZPtNoCandidate = 0;
    Long64_t nZPtInvalid = 0;
    Long64_t nRecoilCorrected = 0;
    Long64_t nRecoilNoCandidate = 0;
    Long64_t nRecoilInvalid = 0;
    for (Long64_t i = 0; i < nEnt; ++i) {
        tree->GetEntry(i);
        if (i > 0 && i % 10000000 == 0) std::cout << "    processed " << i << " / " << nEnt << std::endl;
        string originalCat = numberToCat(cat);
        if (originalCat.size() > 4) continue;
        if (isData && originalCat.empty()) continue;
        applyTauES(originalCat);
        double roccorRelErr[4] = {0.0, 0.0, 0.0, 0.0};
        applyRoccoRCorrection(originalCat, isData, *roccor, roccorRelErr);
        // No dedicated systematic is derived from the official MET-xy
        // correction: the JME POG's own correctionlib payload for it
        // ships only nominal pt/phi outputs, no up/down/syst nodes --
        // unlike JES/JER, which always ship a Total/Up/Down source when a
        // systematic is intended. See DCH_tight.C for the full reasoning.
        if (APPLY_OFFICIAL_MET_CORRECTION && !applyOfficialMETCorrection(*metCorrections)) continue;
        RecoilSystShift recoilShift;
        bool hasRecoilVariation = false;
        if (useRecoilCorrection) {
            double genPx = 0.0;
            double genPy = 0.0;
            double visPx = 0.0;
            double visPy = 0.0;
            bool hasCandidate = false;
            if (isDY || isWJ) {
                hasCandidate = getGenBosonMomentum(genPx, genPy, visPx, visPy);
            }
            if (!hasCandidate) ++nRecoilNoCandidate;
            else if (!applyRecoilCorrection(*recoilCorrector, genPx, genPy, visPx, visPy)) ++nRecoilInvalid;
            else {
                ++nRecoilCorrected;
                if (metSys) {
                    const int recoilNJets = std::max(0, static_cast<int>(std::lround(njets)));
                    const float correctedPx = static_cast<float>(met * std::cos(metphi));
                    const float correctedPy = static_cast<float>(met * std::sin(metphi));
                    recoilShift = computeRecoilSystShift(*metSys, correctedPx, correctedPy, static_cast<float>(genPx), static_cast<float>(genPy), static_cast<float>(visPx), static_cast<float>(visPy), recoilNJets);
                    hasRecoilVariation = recoilShift.valid;
                }
            }
        }
        if (job.year == "2018") {
            if (isData && run >= 319077 && applyHEMveto(originalCat) == "yes") continue;
        }
        double evtWeight = computeEventWeightWithoutObjectSF(originalCat, isData, xsw);
        if (job.year == "2018" && !isData && applyHEMveto(originalCat) == "yes") evtWeight *= 0.35;
        vector<WeightSystematic> weightVariants;
        if (APPLY_ZPT_REWEIGHTING && isDY) {
            GenZMatch genZ;
            vector<std::pair<int,int>> zPtOSPair, zPtSSPair;
            buildPairs(originalCat, zPtOSPair, zPtSSPair);
            if (!zPtOSPair.empty()) {
                const TLorentzVector lep1 = LepV(zPtOSPair[0].first);
                const TLorentzVector lep2 = LepV(zPtOSPair[0].second);
                genZ = matchGenZ(lep1.Pt(), lep1.Eta(), lep1.Phi(), lep2.Pt(), lep2.Eta(), lep2.Phi());
            }
            if (genZ.matched) {
                double zPtWeight = 0.0, zPtWeightErr = 0.0;
                if (!getZPtWeight(zPtWeights.get(), genZ.pt, zPtWeight, zPtWeightErr)) { ++nZPtInvalid; continue; }
                evtWeight *= zPtWeight;
                ++nZPtReweighted;
                if (zPtWeight > 0.0) {
                    weightVariants.push_back({"_zptUp", evtWeight / zPtWeight * (zPtWeight + zPtWeightErr)});
                    weightVariants.push_back({"_zptDown", evtWeight / zPtWeight * std::max(0.0, zPtWeight - zPtWeightErr)});
                }
            }
            else ++nZPtNoCandidate;
        }
        if (!isData) {
            appendPileupWeightVariants(evtWeight, weightVariants);
            appendL1PrefiringWeightVariants(evtWeight, weightVariants);
            if (xsecUnc > 0.0) {
                weightVariants.push_back({"_xsecUp", evtWeight * (1.0 + xsecUnc)});
                weightVariants.push_back({"_xsecDown", evtWeight * std::max(0.0, 1.0 - xsecUnc)});
            }
            const double trigSFNominal = computeTriggerSF(originalCat);
            if (trigSFNominal > 0.0) {
                const TrigSFShift trigShift = computeTriggerSFShifts(originalCat, job.year);
                weightVariants.push_back({"_eTrigUp", evtWeight / trigSFNominal * trigShift.eUp});
                weightVariants.push_back({"_eTrigDown", evtWeight / trigSFNominal * trigShift.eDown});
                weightVariants.push_back({"_muTrigUp", evtWeight / trigSFNominal * trigShift.muUp});
                weightVariants.push_back({"_muTrigDown", evtWeight / trigSFNominal * trigShift.muDown});
            }
        }
        vector<MetSystematic> metVariants;
        if (hasRecoilVariation) {
            metVariants.push_back({"_recoilResponseUp", recoilShift.responseUpPx, recoilShift.responseUpPy});
            metVariants.push_back({"_recoilResponseDown", recoilShift.responseDownPx, recoilShift.responseDownPy});
            metVariants.push_back({"_recoilResolutionUp", recoilShift.resolutionUpPx, recoilShift.resolutionUpPy});
            metVariants.push_back({"_recoilResolutionDown", recoilShift.resolutionDownPx, recoilShift.resolutionDownPy});
        }
        vector<Obj> objs;
        for (int idx = 1; idx <= (int)originalCat.size(); ++idx) {
            char flav = originalCat[idx - 1];
            if (!isValidFlavor(flav)) continue;
            if (ptByIndex(idx) <= 0) continue;
            objs.push_back(makeTightObj(flav, idx));
        }
        if (objs.empty()) continue;
        ++nPassInput;

        TauFakeSFResult fakeSF;
        if (!isData && APPLY_TAU_FAKE_RATE && sfReader) fakeSF = computeTauFakeSF(objs, *sfReader);

        vector<Obj> objsRoccorUp = objs;
        vector<Obj> objsRoccorDown = objs;
        // Gated on originalCat containing a muon (structural), NOT on
        // whether roccorRelErr happens to be nonzero for one -- see
        // DCH_tight.C/DCH_tauFR.C for why the per-object-availability gate
        // silently drops events from the Up/Down histograms.
        const bool hasRoccorVariation = originalCat.find('m') != string::npos;
        for (auto& o : objsRoccorUp) {
            if (o.flav != 'm' || o.tightIdx < 1 || o.tightIdx > 4) continue;
            const double delta = roccorRelErr[o.tightIdx - 1];
            if (delta <= 0.0) continue;
            o.pt *= (1.0 + delta);
        }
        for (auto& o : objsRoccorDown) {
            if (o.flav != 'm' || o.tightIdx < 1 || o.tightIdx > 4) continue;
            const double delta = roccorRelErr[o.tightIdx - 1];
            if (delta <= 0.0) continue;
            o.pt *= std::max(0.0, 1.0 - delta);
        }

        double tauESRelUp[4], tauESRelDown[4];
        computeTauESRelShift(originalCat, tauESRelUp, tauESRelDown);
        // Gated on originalCat containing a tau (structural), NOT on
        // whether a shift was actually available -- see DCH_tight.C for
        // the full explanation of why the availability gate silently
        // collapsed tau-channel _tauESUp/_tauESDown histograms.
        const bool hasTauESVariation = originalCat.find('t') != string::npos;
        vector<Obj> objsTauESUp = objs;
        vector<Obj> objsTauESDown = objs;
        for (auto& o : objsTauESUp) {
            if (o.flav != 't' || o.tightIdx < 1 || o.tightIdx > 4) continue;
            o.pt *= (1.0 + tauESRelUp[o.tightIdx - 1]);
        }
        for (auto& o : objsTauESDown) {
            if (o.flav != 't' || o.tightIdx < 1 || o.tightIdx > 4) continue;
            o.pt *= (1.0 + tauESRelDown[o.tightIdx - 1]);
        }

        const double objectSF = isData ? 1.0 : computeConfigurationObjectSF(objs);
        const double w = evtWeight * fakeSF.sf * objectSF;

        fillConfiguration(objs, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);

        for (const auto& metVar : metVariants) {
            double shiftedMet, shiftedMetPhi;
            cartesianToPolar(metVar.px, metVar.py, shiftedMet, shiftedMetPhi);
            fillConfigurationMetOnly(objs, shiftedMet, shiftedMetPhi, w, finalStateSet, h_met, h_metphi, metVar.suffix);
        }
        if (hasRoccorVariation) {
            fillConfiguration(objsRoccorUp, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_roccorUp");
            fillConfiguration(objsRoccorDown, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_roccorDown");
        }
        if (hasTauESVariation) {
            fillConfiguration(objsTauESUp, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_tauESUp");
            fillConfiguration(objsTauESDown, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_tauESDown");
        }
        if (fakeSF.hasFakeTau) {
            const double wSFUp = evtWeight * fakeSF.sf * (1.0 + fakeSF.relErr) * objectSF;
            const double wSFDown = evtWeight * fakeSF.sf * std::max(0.0, 1.0 - fakeSF.relErr) * objectSF;
            fillConfiguration(objs, wSFUp, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_frSFUp");
            fillConfiguration(objs, wSFDown, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_frSFDown");
        }
        if (!isData) {
            for (const auto& sfVar : computeConfigurationObjectSFSystematics(objs, job.year)) {
                const double wSFVar = evtWeight * fakeSF.sf * sfVar.objectSF;
                fillConfiguration(objs, wSFVar, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, sfVar.suffix);
            }
        }
        for (const auto& variant : weightVariants) {
            const double wVariant = variant.weight * fakeSF.sf * objectSF;
            fillConfiguration(objs, wVariant, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, variant.suffix);
        }
    }
    fout->cd();
    writeHists(h_mZ1);
    writeHists(h_mZ2);
    writeHists(h_mH1);
    writeHists(h_mH2);
    writeHists(h_zPt);
    writeHists(h_met);
    writeHists(h_metphi);
    writeHists(h_LT);
    for (int k = 0; k < NlepMax; ++k) {
        writeHists(h_pt[k]);
        writeHists(h_eta[k]);
        writeHists(h_phi[k]);
        writeHists(h_d0[k]);
        writeHists(h_dZ[k]);
        writeHists(h_iso[k]);
    }
    fout->Close();
    fin->Close();
    std::cout << "  [done] input_passed=" << nPassInput;
    if (APPLY_ZPT_REWEIGHTING && isDY) {
        std::cout << " zpt_reweighted=" << nZPtReweighted << " zpt_no_candidate=" << nZPtNoCandidate << " zpt_invalid=" << nZPtInvalid;
    }
    if (useRecoilCorrection) {
        std::cout << " recoil_corrected=" << nRecoilCorrected << " recoil_no_candidate=" << nRecoilNoCandidate << " recoil_invalid=" << nRecoilInvalid;
    }
    std::cout << " written -> " << outName << std::endl;
    return 0;
}

void DCH_tauFR_SF(string inYear="2018", int firstFile=0, int nFilesToRun=-1, int nProc=1, string processFilter="")
{
    TString outDirBase = "/eos/user/a/atahmad/DCH_offline_analysis/hists/run2_hists_tauFR_SF";
    TString outDir = Form("%s/%s", outDirBase.Data(), inYear.c_str());
    gSystem->mkdir(outDir, kTRUE);
    std::map<string, vector<string>> fileMap = getFileMap(inYear);
    vector<std::pair<string,string>> files;
    for (const auto& item : fileMap) {
        if (!processFilter.empty() && item.first != processFilter) continue;
        for (const auto& fName : item.second) files.push_back({item.first, fName});
    }
    if (files.empty()) { std::cerr << "ERROR: no files found for year " << inYear << std::endl; return; }
    int lastFile = files.size();
    if (nFilesToRun >= 0) lastFile = std::min((int)files.size(), firstFile + nFilesToRun);
    if (firstFile < 0 || firstFile >= (int)files.size()) { std::cout << "No file for firstFile=" << firstFile << std::endl; return; }
    double lumi = (inYear=="2016preVFP") ? 19520.0 : (inYear=="2016postVFP") ? 16810.0 : (inYear=="2017") ? 41480.0 : (inYear=="2018") ? 59830.0 : 19520.0 + 16810.0 + 41480.0 + 59830.0;
    vector<FileJob> jobs;
    for (int jf = firstFile; jf < lastFile; ++jf) {
        FileJob job;
        job.year = inYear;
        job.process = files[jf].first;
        job.fileName = files[jf].second;
        job.outDir = outDir.Data();
        job.lumi = lumi;
        job.index = jf;
        job.total = files.size();
        jobs.push_back(job);
    }
    if (nProc <= 1 || jobs.size() <= 1) {
        for (auto& job : jobs) ProcessTauFRSFFile(job);
    }
    else {
        ROOT::TProcessExecutor pool(nProc);
        auto results = pool.Map(ProcessTauFRSFFile, jobs);
        (void)results;
    }
}
