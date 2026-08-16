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
#include "DCH_modules/TruthMatching.h"
#include "DCH_modules/ZPtReweight.h"
#include "DCH_modules/FileJob.h"
#include "DCH_modules/TauFakeRate.h"
#include "DCH_modules/EtauFakeRate.h"
#include "DCH_modules/EventWeights.h"
#include "DCH_modules/PairBuilder.h"
#include "DCH_modules/TauFRSystematics.h"
#include "DCH_modules/MEtSysWrapper.h"
#include "DCH_modules/SystematicPlan.h"
#include "DCH_modules/FlatXsecSystematic.h"

using std::string;
using std::vector;
using std::unordered_map;

const bool SKIP_DATA = false;
const bool APPLY_TAU_FAKE_RATE = true;
const bool APPLY_ETAU_FAKE_RATE = true;

const string FR_FILE = "Dependencies/fake_rates/tau_fake_rates/DY_data_tau_fake_rate_2D.root";

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

struct TauCand : public Obj {
    double fr = 0.0;
    double frErr = 0.0;
    bool isEtauSource = false;
};

struct FRSystSource {
    string suffix;
    const TauFRReader* reader;
    bool isEtau;
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

vector<string> allHistKeys() {
    vector<string> keys;
    for (const auto& ch : finalStates) {
        keys.push_back(ch);
        if (isDilep(ch)) {
            keys.push_back(ch + "_OS_Zwin");
            keys.push_back(ch + "_SS_Zwin");
            keys.push_back(ch + "_OS_Zveto");
            keys.push_back(ch + "_SS_Zveto");
        }
        for (const auto& reg : tauRegions) keys.push_back(ch + "_" + reg);
    }
    return keys;
}

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

void enableBranchesTauFR(TTree* t, bool isData) {
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
        "Generator_weight", "brWeight", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Down",
        "weightPUtruejson", "weightPUtruejson_up", "weightPUtruejson_down",
    };
    if (!isData) {
        const char* looseBrs[] = {"lpt", "leta", "lphi", "lmass", "lq", "lflavor", "liso", "ld0", "ldZ", "gen_match"};
        for (auto br : looseBrs) if (t->GetBranch(br)) t->SetBranchStatus(br, 1);
    }
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

TauCand makeLooseTauCand(size_t i, const TauFRReader& frReader, bool isEtau) {
    TauCand c;
    c.flav = 't';
    c.pt = lpt->at(i);
    c.eta = leta->at(i);
    c.phi = lphi->at(i);
    c.mass = lmass ? lmass->at(i) : 0.0;
    c.q = lq->at(i);
    c.d0 = ld0 ? ld0->at(i) : -99.0;
    c.dZ = ldZ ? ldZ->at(i) : -99.0;
    c.iso = liso ? liso->at(i) : -99.0;
    c.fromTight = false;
    c.tightIdx = -1;
    c.objectOrder = 100 + static_cast<int>(i);
    c.isEtauSource = isEtau;
    frReader.get(c.pt, c.eta, c.fr, c.frErr);
    return c;
}

TauCand makeTightFakeTauCand(const Obj& o, const TauFRReader& frReader, bool isEtau) {
    TauCand c;
    c.flav = 't';
    c.pt = o.pt;
    c.eta = o.eta;
    c.phi = o.phi;
    c.mass = o.mass;
    c.q = o.q;
    c.d0 = o.d0;
    c.dZ = o.dZ;
    c.iso = o.iso;
    c.fromTight = true;
    c.tightIdx = o.tightIdx;
    c.objectOrder = o.objectOrder;
    c.isEtauSource = isEtau;
    frReader.get(c.pt, c.eta, c.fr, c.frErr);
    return c;
}

bool overlapsAny(const Obj& x, const vector<Obj>& objs) {
    for (const auto& o : objs) if (getDR(x.eta, x.phi, o.eta, o.phi) <= 0.4) return true;
    return false;
}

bool overlapsAnyTauCand(const Obj& x, const vector<TauCand>& objs) {
    for (const auto& o : objs) if (getDR(x.eta, x.phi, o.eta, o.phi) <= 0.4) return true;
    return false;
}

void generateAndFillConfigs(const vector<Obj>& baseObjects, const vector<TauCand>& tauCands, double evtWeight, const std::unordered_set<string>& finalStateSet, const vector<FRSystSource>& frSourceVariants, const vector<WeightSystematic>& weightVariants, const vector<MetSystematic>& metVariants, unordered_map<string,TH1D*>& h_mZ1, unordered_map<string,TH1D*>& h_mZ2, unordered_map<string,TH1D*>& h_mH1, unordered_map<string,TH1D*>& h_mH2, unordered_map<string,TH1D*>& h_zPt, unordered_map<string,TH1D*>& h_met, unordered_map<string,TH1D*>& h_metphi, unordered_map<string,TH1D*>& h_LT, unordered_map<string,TH1D*> h_pt[NlepMax], unordered_map<string,TH1D*> h_eta[NlepMax], unordered_map<string,TH1D*> h_phi[NlepMax], unordered_map<string,TH1D*> h_d0[NlepMax], unordered_map<string,TH1D*> h_dZ[NlepMax], unordered_map<string,TH1D*> h_iso[NlepMax], const vector<Obj>& baseObjectsRoccorUp, const vector<Obj>& baseObjectsRoccorDown, bool hasRoccorVariation)
{
    const int nCand = tauCands.size();
    const int nBase = baseObjects.size();
    const int maxPick = std::min(4 - nBase, nCand);
    if (nBase > 4) return;
    if (maxPick < 0) return;
    const int nStates = 1 << nCand;
    for (int mask = 0; mask < nStates; ++mask) {
        int nPick = __builtin_popcount((unsigned)mask);
        if (nPick > maxPick) continue;
        int nTotal = nBase + nPick;
        if (nTotal < 2 || nTotal > 4) continue;
        double wFR = 1.0;
        double relVarSum = 0.0;
        vector<Obj> cfg = baseObjects;
        vector<Obj> cfgRoccorUp = baseObjectsRoccorUp;
        vector<Obj> cfgRoccorDown = baseObjectsRoccorDown;
        for (int i = 0; i < nCand; ++i) {
            double f = tauCands[i].fr;
            double errF = tauCands[i].frErr;
            double factor;
            if (mask & (1 << i)) {
                factor = f;
                wFR *= f;
                cfg.push_back(tauCands[i]);
                cfgRoccorUp.push_back(tauCands[i]);
                cfgRoccorDown.push_back(tauCands[i]);
            }
            else {
                factor = 1.0 - f;
                wFR *= factor;
            }
            if (factor > 1e-6) {
                double relErr = errF / factor;
                relVarSum += relErr * relErr;
            }
        }
        std::stable_sort(cfg.begin(), cfg.end(), [](const Obj& a, const Obj& b) { return a.objectOrder < b.objectOrder; });
        if (hasRoccorVariation) {
            std::stable_sort(cfgRoccorUp.begin(), cfgRoccorUp.end(), [](const Obj& a, const Obj& b) { return a.objectOrder < b.objectOrder; });
            std::stable_sort(cfgRoccorDown.begin(), cfgRoccorDown.end(), [](const Obj& a, const Obj& b) { return a.objectOrder < b.objectOrder; });
        }
        if (wFR <= 0.0) continue;
        const double relErrTotal = std::sqrt(relVarSum);
        const double wFRUp = wFR * (1.0 + relErrTotal);
        const double wFRDown = wFR * std::max(0.0, 1.0 - relErrTotal);
        const double objectSF = computeConfigurationObjectSF(cfg);
        const double w = evtWeight * wFR * objectSF;
        const double wUp = evtWeight * wFRUp * objectSF;
        const double wDown = evtWeight * wFRDown * objectSF;
        fillConfiguration(cfg, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "");
        fillConfiguration(cfg, wUp, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_frUp");
        fillConfiguration(cfg, wDown, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_frDown");
        if (hasRoccorVariation) {
            fillConfiguration(cfgRoccorUp, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_roccorUp");
            fillConfiguration(cfgRoccorDown, w, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_roccorDown");
        }
        for (const auto& source : frSourceVariants) {
            double wFRSource = 1.0;
            for (int i = 0; i < nCand; ++i) {
                double f;
                if (tauCands[i].isEtauSource == source.isEtau) {
                    double altF = 0.0, altErr = 0.0;
                    source.reader->get(tauCands[i].pt, tauCands[i].eta, altF, altErr);
                    f = altF;
                } else {
                    f = tauCands[i].fr;
                }
                wFRSource *= (mask & (1 << i)) ? f : (1.0 - f);
            }
            if (wFRSource <= 0.0) continue;
            const double wSource = evtWeight * wFRSource * objectSF;
            fillConfiguration(cfg, wSource, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, source.suffix);
        }
        for (const auto& variant : weightVariants) {
            const double wVariant = variant.weight * wFR * objectSF;
            fillConfiguration(cfg, wVariant, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, variant.suffix);
        }
        for (const auto& metVar : metVariants) {
            double shiftedMet, shiftedMetPhi;
            cartesianToPolar(metVar.px, metVar.py, shiftedMet, shiftedMetPhi);
            fillConfigurationMetOnly(cfg, shiftedMet, shiftedMetPhi, w, finalStateSet, h_met, h_metphi, metVar.suffix);
        }
    }
}

int ProcessTauFRFile(FileJob job) {
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
    enableBranchesTauFR(tree, isData);
    MyBranch(tree);
    tree->SetCacheSize(200 * 1024 * 1024);
    tree->AddBranchToCache("*", kTRUE);
    std::shared_ptr<TauFRReader> frReader;
    if (!isData && APPLY_TAU_FAKE_RATE) frReader = getCachedTauFRReader(FR_FILE, job.year);
    std::shared_ptr<TauFRReader> etauFrReader;
    if (!isData && APPLY_TAU_FAKE_RATE && APPLY_ETAU_FAKE_RATE) etauFrReader = getCachedEtauFRReader(job.year);
    vector<std::shared_ptr<TauFRReader>> frSystReadersOwned;
    vector<FRSystSource> frSystReaders;
    if (!isData && APPLY_TAU_FAKE_RATE) {
        for (const auto& source : kTauFRSystematicSources) {
            auto up = getCachedSystematicTauFRReader(source.first, "up", job.year);
            auto down = getCachedSystematicTauFRReader(source.first, "down", job.year);
            if (up) frSystReadersOwned.push_back(up);
            if (down) frSystReadersOwned.push_back(down);
            if (up) frSystReaders.push_back({"_fr" + source.second + "Up", up.get(), false});
            if (down) frSystReaders.push_back({"_fr" + source.second + "Down", down.get(), false});
        }
        if (APPLY_ETAU_FAKE_RATE && etauFrReader) {
            auto etauStatUp = getCachedEtauFRStatReader(job.year, "up");
            auto etauStatDown = getCachedEtauFRStatReader(job.year, "down");
            if (etauStatUp) frSystReadersOwned.push_back(etauStatUp);
            if (etauStatDown) frSystReadersOwned.push_back(etauStatDown);
            if (etauStatUp) frSystReaders.push_back({"_frEtauStatUp", etauStatUp.get(), true});
            if (etauStatDown) frSystReaders.push_back({"_frEtauStatDown", etauStatDown.get(), true});
        }
    }
    std::shared_ptr<MEtSys> metSys;
    if (useRecoilCorrection) {
        metSys = getCachedMEtSys(job.year);
        if (!metSys) std::cerr << "  [warn] MEtSys payload not found for " << job.year << "; recoil systematic variations will be skipped" << std::endl;
    }
    TH1D* hNWEvts = (TH1D*)fin->Get("hNWEvts");
    if (!hNWEvts) hNWEvts = (TH1D*)fin->Get("hNEvts");
    double denom = hNWEvts ? hNWEvts->Integral() : 0.0;
    double xsw = (!isData && denom > 0.0) ? job.lumi * XSec(baseName) / denom : 1.0;
    TString outName = Form("%s/hist_%s", job.outDir.c_str(), baseName.c_str());
    TFile* fout = new TFile(outName, "RECREATE");
    fout->SetCompressionSettings(ROOT::CompressionSettings(ROOT::kZSTD, 1));
    if (hNWEvts) hNWEvts->Write();
    unordered_map<string,TH1D*> h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT;
    unordered_map<string,TH1D*> h_pt[NlepMax], h_eta[NlepMax], h_phi[NlepMax], h_d0[NlepMax], h_dZ[NlepMax], h_iso[NlepMax];
    std::unordered_set<string> finalStateSet(finalStates.begin(), finalStates.end());
    Long64_t nEnt = tree->GetEntriesFast();
    std::cout << "  entries in tree: " << nEnt << std::endl;
    Long64_t nPassInput = 0;
    Long64_t nFilledSeed = 0;
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
        const double rawMetPx = met * std::cos(metphi);
        const double rawMetPy = met * std::sin(metphi);
        if (APPLY_OFFICIAL_MET_CORRECTION && !applyOfficialMETCorrection(*metCorrections)) continue;
        bool hasMetPhiVariation = false;
        double metPhiDeltaPx = 0.0, metPhiDeltaPy = 0.0;
        if (APPLY_OFFICIAL_MET_CORRECTION) {
            metPhiDeltaPx = met * std::cos(metphi) - rawMetPx;
            metPhiDeltaPy = met * std::sin(metphi) - rawMetPy;
            hasMetPhiVariation = true;
        }
        RecoilSystShift recoilShift;
        bool hasRecoilVariation = false;
        if (useRecoilCorrection) {
            double genPx = 0.0;
            double genPy = 0.0;
            double visPx = 0.0;
            double visPy = 0.0;
            bool hasCandidate = false;
            if (isDY) {
                hasCandidate = getPromptTruthDileptonXY(originalCat, genPx, genPy);
                visPx = genPx;
                visPy = genPy;
            }
            else if (isWJ) {
                hasCandidate = getPromptTruthLeptonXY(originalCat, visPx, visPy);
                if (hasCandidate) {
                    genPx = visPx + met * std::cos(metphi);
                    genPy = visPy + met * std::sin(metphi);
                }
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
            double promptTruthZPt = 0.0;
            if (getPromptTruthDileptonPt(originalCat, promptTruthZPt)) {
                double zPtWeight = 0.0, zPtWeightErr = 0.0;
                if (!getZPtWeight(zPtWeights.get(), promptTruthZPt, zPtWeight, zPtWeightErr)) { ++nZPtInvalid; continue; }
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
        }
        vector<MetSystematic> metVariants;
        if (hasMetPhiVariation) {
            const double finalPx = met * std::cos(metphi);
            const double finalPy = met * std::sin(metphi);
            metVariants.push_back({"_metPhiCorrUp", finalPx + metPhiDeltaPx, finalPy + metPhiDeltaPy});
            metVariants.push_back({"_metPhiCorrDown", finalPx - metPhiDeltaPx, finalPy - metPhiDeltaPy});
        }
        if (hasRecoilVariation) {
            metVariants.push_back({"_recoilResponseUp", recoilShift.responseUpPx, recoilShift.responseUpPy});
            metVariants.push_back({"_recoilResponseDown", recoilShift.responseDownPx, recoilShift.responseDownPy});
            metVariants.push_back({"_recoilResolutionUp", recoilShift.resolutionUpPx, recoilShift.resolutionUpPy});
            metVariants.push_back({"_recoilResolutionDown", recoilShift.resolutionDownPx, recoilShift.resolutionDownPy});
        }
        vector<Obj> baseObjects;
        vector<TauCand> tauCands;
        for (int idx = 1; idx <= (int)originalCat.size(); ++idx) {
            char flav = originalCat[idx - 1];
            if (!isValidFlavor(flav)) continue;
            if (ptByIndex(idx) <= 0) continue;
            Obj o = makeTightObj(flav, idx);
            if (flav != 't') { baseObjects.push_back(o); continue; }
            if (isData || !APPLY_TAU_FAKE_RATE) baseObjects.push_back(o);
            else {
                int gpf = genPartFlavByIndex(idx);
                if (gpf == 5) baseObjects.push_back(o);
                else {
                    bool useEtau = (gpf == 1 && etauFrReader != nullptr);
                    const TauFRReader& reader = useEtau ? *etauFrReader : *frReader;
                    tauCands.push_back(makeTightFakeTauCand(o, reader, useEtau));
                }
            }
        }
        if (!isData && APPLY_TAU_FAKE_RATE && lpt && lflavor && gen_match) {
            for (size_t il = 0; il < lpt->size(); ++il) {
                if (std::abs(lflavor->at(il)) != 15) continue;
                if (lpt->at(il) <= 0) continue;
                if (il >= gen_match->size()) continue;
                int gm = gen_match->at(il);
                if (gm == 5) continue;
                bool useEtau = (gm == 1 && etauFrReader != nullptr);
                const TauFRReader& reader = useEtau ? *etauFrReader : *frReader;
                TauCand c = makeLooseTauCand(il, reader, useEtau);
                if (overlapsAny(c, baseObjects)) continue;
                if (overlapsAnyTauCand(c, tauCands)) continue;
                tauCands.push_back(c);
            }
        }
        if (baseObjects.empty() && tauCands.empty()) continue;
        ++nPassInput;
        vector<Obj> baseObjectsRoccorUp = baseObjects;
        vector<Obj> baseObjectsRoccorDown = baseObjects;
        bool hasRoccorVariation = false;
        for (auto& o : baseObjectsRoccorUp) {
            if (o.flav != 'm' || o.tightIdx < 1 || o.tightIdx > 4) continue;
            const double delta = roccorRelErr[o.tightIdx - 1];
            if (delta <= 0.0) continue;
            hasRoccorVariation = true;
            o.pt *= (1.0 + delta);
        }
        for (auto& o : baseObjectsRoccorDown) {
            if (o.flav != 'm' || o.tightIdx < 1 || o.tightIdx > 4) continue;
            const double delta = roccorRelErr[o.tightIdx - 1];
            if (delta <= 0.0) continue;
            o.pt *= std::max(0.0, 1.0 - delta);
        }
        if (isData || !APPLY_TAU_FAKE_RATE) {
            fillConfiguration(baseObjects, evtWeight, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
            for (const auto& metVar : metVariants) {
                double shiftedMet, shiftedMetPhi;
                cartesianToPolar(metVar.px, metVar.py, shiftedMet, shiftedMetPhi);
                fillConfigurationMetOnly(baseObjects, shiftedMet, shiftedMetPhi, evtWeight, finalStateSet, h_met, h_metphi, metVar.suffix);
            }
            if (hasRoccorVariation) {
                fillConfiguration(baseObjectsRoccorUp, evtWeight, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_roccorUp");
                fillConfiguration(baseObjectsRoccorDown, evtWeight, finalStateSet, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, "_roccorDown");
            }
        }
        else generateAndFillConfigs(baseObjects, tauCands, evtWeight, finalStateSet, frSystReaders, weightVariants, metVariants, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso, baseObjectsRoccorUp, baseObjectsRoccorDown, hasRoccorVariation);
        ++nFilledSeed;
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
    std::cout << "  [done] input_passed=" << nPassInput << " seeds=" << nFilledSeed;
    if (APPLY_ZPT_REWEIGHTING && isDY) {
        std::cout << " zpt_reweighted=" << nZPtReweighted << " zpt_no_candidate=" << nZPtNoCandidate << " zpt_invalid=" << nZPtInvalid;
    }
    if (useRecoilCorrection) {
        std::cout << " recoil_corrected=" << nRecoilCorrected << " recoil_no_candidate=" << nRecoilNoCandidate << " recoil_invalid=" << nRecoilInvalid;
    }
    std::cout << " written -> " << outName << std::endl;
    return 0;
}

void DCH_tauFR(string inYear="2018", int firstFile=0, int nFilesToRun=-1, int nProc=1, string processFilter="")
{
    TString outDirBase = "hists/run2_hists_tauFR_etau_roccor";
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
        for (auto& job : jobs) ProcessTauFRFile(job);
    }
    else {
        ROOT::TProcessExecutor pool(nProc);
        auto results = pool.Map(ProcessTauFRFile, jobs);
        (void)results;
    }
}
