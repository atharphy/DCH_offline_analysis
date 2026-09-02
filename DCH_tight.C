// DCH_tight.C
//
// Tight-object nominal histograms for data and MC.
// Histograms are split by final state and CR/VR/SR region.
// Dilepton channels additionally get OS/SS x Zwin/Zveto categories.
//
// Compile and run one file:
//   root -l -b -q 'DCH_tight.C+("2018",0,1,1)'
//
// Run all files sequentially:
//   root -l -b -q 'DCH_tight.C+("2018",0,-1,1)'
//
// Run all files with process parallelism:
//   root -l -b -q 'DCH_tight.C+("2018",0,-1,12)'
//
// Run only one process group (e.g. just the signal MC), all its files:
//   root -l -b -q 'DCH_tight.C+("2018",0,-1,1,"signal")'
//
// Arguments:
//   inYear        = year
//   firstFile     = first file index
//   nFilesToRun   = number of files, -1 means all remaining files
//   nProc         = number of parallel worker processes
//   processFilter = if non-empty, only run files whose FileMap.h process
//                   key exactly matches this (e.g. "signal"); default ""
//                   runs every process as before

#if !defined(__CLING__)
#pragma GCC optimize("O3,unroll-loops")
#endif

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1D.h"
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
#include <utility>
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
#include "DCH_modules/EventWeights.h"
#include "DCH_modules/PairBuilder.h"
#include "DCH_modules/MEtSysWrapper.h"
#include "DCH_modules/SystematicPlan.h"
#include "DCH_modules/FlatXsecSystematic.h"
#include "DCH_modules/TauESSystematic.h"

using std::string;
using std::unordered_map;
using std::vector;

const int OUTPUT_COMPRESSION = ROOT::CompressionSettings(ROOT::kZSTD, 1);
const Long64_t TREE_CACHE_SIZE = 200LL * 1024LL * 1024LL;
const Long64_t CHUNK_SIZE = 500000LL;

struct EventKinematics {
    double LT = 0.0;
    bool hasMZ1 = false, hasMZ2 = false, hasMH1 = false, hasMH2 = false, hasZPt = false;
    double mZ1 = 0.0, mZ2 = 0.0, mH1 = 0.0, mH2 = 0.0, zPt = 0.0;
};

bool passTightEvent(const string& catstr) {
    if (catstr.size() < 2 || catstr.size() > 4) return false;
    for (int idx=1; idx<=static_cast<int>(catstr.size()); ++idx) {
        if (!isValidFlavor(catstr[idx-1])) return false;
        if (ptByIndex(idx) <= 0.0) return false;
    }
    return true;
}

bool passChargeTopology(const string& catstr) {
    if (catstr.size()==4) return q_1 + q_2 + q_3 + q_4 == 0;
    if (catstr.size()==3) return std::abs(q_1 + q_2 + q_3) != 3;
    return true;
}

bool hasDuplicateTightObjects(const string& catstr) {
    const int nLep = static_cast<int>(catstr.size());
    for (int i=1; i<=nLep; ++i) {
        for (int j=i+1; j<=nLep; ++j) {
            if (getDR(etaByIndex(i), phiByIndex(i), etaByIndex(j), phiByIndex(j)) <= 0.4) return true;
        }
    }
    return false;
}

double getLT(const string& catstr) { double LT = 0.0; for (int idx=1; idx<=static_cast<int>(catstr.size()); ++idx) LT += ptByIndex(idx); return LT; }

void deleteHists(unordered_map<string,TH1D*>& histMap) {
    for (auto& item : histMap) delete item.second;
    histMap.clear();
}

double computeNominalWeight(const string& catstr, bool isData, double xsw) {
    if (isData) return 1.0;
    const int nLep = static_cast<int>(catstr.size());
    double objectSF = 1.0;
    for (int idx = 1; idx <= nLep; ++idx) {
        const char flavor = catstr[idx - 1];
        if (flavor == 'e') {
            objectSF *= idSFByIndex(idx);
            objectSF *= isoSFByIndex(idx);
        }
        else if (flavor == 'm') {
            objectSF *= idSFByIndex(idx);
            objectSF *= isoSFByIndex(idx);
        }
        else if (flavor == 't') {
            objectSF *= tauEleSFByIndex(idx);
            objectSF *= tauMuSFByIndex(idx);
            objectSF *= tauJetSFByIndex(idx);
        }
    }
    return computeBaseGenWeight(xsw) * objectSF * computeTriggerSF(catstr);
}

struct ObjectSFSystEntry { string suffix; double weight; };

vector<ObjectSFSystEntry> computeObjectSFSystematics(const string& catstr, bool isData, double baseWeight, const string& year) {
    (void)year;
    vector<ObjectSFSystEntry> result;
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

double muonTrigPtThreshold(const string& year) {
    if (year == "2017") return 29.0;
    return 26.0;
}

// Reads the skim's precomputed TrigSF_Up/Down branches directly (absolute
// shifted values) instead of recomputing a relative uncertainty via
// correctionlib/TGraphAsymmErrors.
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

EventKinematics buildEventKinematics(const string& catstr, const vector<std::pair<int,int>>& OS_pair, const vector<std::pair<int,int>>& SS_pair) {
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

vector<string> collectHistKeys(const string& catstr, const EventKinematics& kin, const vector<std::pair<int,int>>& OS_pair) {
    vector<string> keys;
    keys.push_back(catstr);
    if (catstr.size() == 2) {
        const string zFlag = dilepZFlag(catstr);
        const bool isOS = q_1 * q_2 < 0;
        const bool isSS = q_1 * q_2 > 0;
        string dileptonKey;
        if (isOS && zFlag=="Zwin") dileptonKey = catstr + "_OS_Zwin";
        else if (isSS && zFlag=="Zwin") dileptonKey = catstr + "_SS_Zwin";
        else if (isOS && zFlag=="Zveto") dileptonKey = catstr + "_OS_Zveto";
        else if (isSS && zFlag=="Zveto") dileptonKey = catstr + "_SS_Zveto";
        if (!dileptonKey.empty()) keys.push_back(dileptonKey);
    }
    const string region = classifyTauRegion(catstr, kin.LT, OS_pair);
    keys.push_back(catstr + "_" + region);
    return keys;
}

void fillMetOnly(const string& key, double met_, double metphi_, double weight, unordered_map<string,TH1D*>& h_met, unordered_map<string,TH1D*>& h_metphi) {
    fillOne(h_met, key, met_, weight, S_met);
    fillOne(h_metphi, key, metphi_, weight, S_metphi);
}

void fillAllVars(const string& key, const string& catstr, const EventKinematics& kin, double weight, unordered_map<string,TH1D*>& h_mZ1, unordered_map<string,TH1D*>& h_mZ2, unordered_map<string,TH1D*>& h_mH1, unordered_map<string,TH1D*>& h_mH2, unordered_map<string,TH1D*>& h_zPt, unordered_map<string,TH1D*>& h_met, unordered_map<string,TH1D*>& h_metphi, unordered_map<string,TH1D*>& h_LT, unordered_map<string,TH1D*> h_pt[NlepMax], unordered_map<string,TH1D*> h_eta[NlepMax], unordered_map<string,TH1D*> h_phi[NlepMax], unordered_map<string,TH1D*> h_d0[NlepMax], unordered_map<string,TH1D*> h_dZ[NlepMax], unordered_map<string,TH1D*> h_iso[NlepMax])
{
    fillOne(h_met, key, met, weight, S_met);
    fillOne(h_metphi, key, metphi, weight, S_metphi);
    fillOne(h_LT, key, kin.LT, weight, S_LT);
    if (kin.hasMZ1) fillOne(h_mZ1, key, kin.mZ1, weight, S_mZ1);
    if (kin.hasMZ2) fillOne(h_mZ2, key, kin.mZ2, weight, S_mZ2);
    if (kin.hasMH1) fillOne(h_mH1, key, kin.mH1, weight, S_mH1);
    if (kin.hasMH2) fillOne(h_mH2, key, kin.mH2, weight, S_mH2);
    if (kin.hasZPt && std::isfinite(kin.zPt)) fillOne(h_zPt, key, kin.zPt, weight, S_zPt);
    const int nLep = std::min(static_cast<int>(catstr.size()), NlepMax);
    for (int i = 0; i < nLep; ++i) {
        const int idx = i + 1;
        fillOne(h_pt[i], key, ptByIndex(idx), weight, S_pt[i]);
        fillOne(h_eta[i], key, etaByIndex(idx), weight, S_eta[i]);
        fillOne(h_phi[i], key, phiByIndex(idx), weight, S_phi[i]);
        fillOne(h_d0[i], key, d0ByIndex(idx), weight, S_d0[i]);
        fillOne(h_dZ[i], key, dZByIndex(idx), weight, S_dZ[i]);
        fillOne(h_iso[i], key, isoByIndex(idx), weight, S_iso[i]);
    }
}

void enableBranches(TTree* tree, bool isData) {
    tree->SetBranchStatus("*", 0);
    const char* commonBranches[] = {
        "cat","run","met","metphi","nPV","nPVGood",
        "pt_1","pt_2","pt_3","pt_4",
        "eta_1","eta_2","eta_3","eta_4",
        "phi_1","phi_2","phi_3","phi_4",
        "m_1","m_2","m_3","m_4",
        "q_1","q_2","q_3","q_4",
        "d0_1","d0_2","d0_3","d0_4",
        "dZ_1","dZ_2","dZ_3","dZ_4",
        "iso_1","iso_2","iso_3","iso_4",
        "TauES_1","TauES_2","TauES_3","TauES_4"
    };
    for (const char* branchName : commonBranches) if (tree->GetBranch(branchName)) tree->SetBranchStatus(branchName, 1);
    if (!isData) {
        const char* mcBranches[] = {
            "Generator_weight","brWeight","L1PreFiringWeight_Nom","L1PreFiringWeight_Up","L1PreFiringWeight_Down",
            "weightPUtruejson","weightPUtruejson_up","weightPUtruejson_down",
            "genPartFlav_1","genPartFlav_2","genPartFlav_3","genPartFlav_4",
            "pt_1_tr","pt_2_tr","pt_3_tr","pt_4_tr",
            "eta_1_tr","eta_2_tr","eta_3_tr","eta_4_tr",
            "phi_1_tr","phi_2_tr","phi_3_tr","phi_4_tr",
            "IDSF_1","IDSF_2","IDSF_3","IDSF_4",
            "ISOSF_1","ISOSF_2","ISOSF_3","ISOSF_4",
            "TrigSF_1","TrigSF_2","TrigSF_3","TrigSF_4",
            "isTrig_1","isTrig_2",
            "TauVsEleIDSF_1","TauVsEleIDSF_2","TauVsEleIDSF_3","TauVsEleIDSF_4",
            "TauVsMuIDSF_1","TauVsMuIDSF_2","TauVsMuIDSF_3","TauVsMuIDSF_4",
            "TauVsJetIDSF_1","TauVsJetIDSF_2","TauVsJetIDSF_3","TauVsJetIDSF_4",
            "IDSF_Up_1","IDSF_Up_2","IDSF_Up_3","IDSF_Up_4",
            "IDSF_Down_1","IDSF_Down_2","IDSF_Down_3","IDSF_Down_4",
            "ISOSF_Up_1","ISOSF_Up_2","ISOSF_Up_3","ISOSF_Up_4",
            "ISOSF_Down_1","ISOSF_Down_2","ISOSF_Down_3","ISOSF_Down_4",
            "TrigSF_Up_1","TrigSF_Up_2","TrigSF_Up_3","TrigSF_Up_4",
            "TrigSF_Down_1","TrigSF_Down_2","TrigSF_Down_3","TrigSF_Down_4",
            "TauVsEleIDSF_Up_1","TauVsEleIDSF_Up_2","TauVsEleIDSF_Up_3","TauVsEleIDSF_Up_4",
            "TauVsEleIDSF_Down_1","TauVsEleIDSF_Down_2","TauVsEleIDSF_Down_3","TauVsEleIDSF_Down_4",
            "TauVsMuIDSF_Up_1","TauVsMuIDSF_Up_2","TauVsMuIDSF_Up_3","TauVsMuIDSF_Up_4",
            "TauVsMuIDSF_Down_1","TauVsMuIDSF_Down_2","TauVsMuIDSF_Down_3","TauVsMuIDSF_Down_4",
            "TauVsJetIDSF_Up_1","TauVsJetIDSF_Up_2","TauVsJetIDSF_Up_3","TauVsJetIDSF_Up_4",
            "TauVsJetIDSF_Down_1","TauVsJetIDSF_Down_2","TauVsJetIDSF_Down_3","TauVsJetIDSF_Down_4",
            "TauES_Up_1","TauES_Up_2","TauES_Up_3","TauES_Up_4",
            "TauES_Down_1","TauES_Down_2","TauES_Down_3","TauES_Down_4",
            "GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass",
            "GenPart_pdgId","GenPart_genPartIdxMother","GenPart_status","GenPart_statusFlags",
            "GenVisTau_genPartIdxMother","GenVisTau_eta","GenVisTau_phi"
        };
        for (const char* branchName : mcBranches) if (tree->GetBranch(branchName)) tree->SetBranchStatus(branchName, 1);
    }
}

int ProcessTightFile(FileJob job) {
    std::cout << "\n[file " << job.index + 1 << "/" << job.total << "] " << job.fileName << std::endl;
    TFile* fin = TFile::Open(job.fileName.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "  [skip] cannot open file" << std::endl;
        if (fin) fin->Close();
        delete fin;
        return 1;
    }
    TTree* tree = dynamic_cast<TTree*>(fin->Get("Events"));
    if (!tree) {
        std::cerr << "  [skip] no Events tree" << std::endl;
        fin->Close();
        delete fin;
        return 2;
    }
    const string baseName = gSystem->BaseName(job.fileName.c_str());
    const bool isData = XSec(baseName) == 1;
    const bool isDY = job.process == "DY" || job.process == "DY10_50";
    const bool isWJ = job.process == "WJ";
    const double xsecUnc = getFlatXsecUncertainty(baseName);
    std::unique_ptr<TH1D> zPtWeights;
    if (APPLY_ZPT_REWEIGHTING && isDY) {
        zPtWeights.reset(loadZPtWeights(job.year));
        if (!zPtWeights) { fin->Close(); delete fin; return 5; }
    }
    std::unique_ptr<RecoilCorrector> recoilCorrector;
    const bool useRecoilCorrection = !isData && ((APPLY_DY_RECOIL_CORRECTION && isDY) || (APPLY_WJ_RECOIL_CORRECTION && isWJ));
    if (useRecoilCorrection) {
        recoilCorrector = loadRecoilCorrector(job.year);
        if (!recoilCorrector) { fin->Close(); delete fin; return 6; }
    }
    std::shared_ptr<MEtSys> metSys;
    if (useRecoilCorrection) {
        metSys = getCachedMEtSys(job.year);
        if (!metSys) std::cerr << "  [warn] MEtSys payload not found for " << job.year << "; recoil systematic variations will be skipped" << std::endl;
    }
    std::unique_ptr<OfficialMETCorrections> metCorrections;
    if (APPLY_OFFICIAL_MET_CORRECTION) {
        metCorrections = loadOfficialMETCorrections(job.year, isData);
        if (!metCorrections) {
            std::cerr << "  [skip] official MET corrections could not be loaded for " << job.year << std::endl;
            fin->Close();
            delete fin;
            return 4;
        }
    }
    std::unique_ptr<RoccoR> roccor = loadRoccoRCorrections(job.year);
    if (!roccor) {
        std::cerr << "  [skip] RoccoR muon momentum corrections could not be loaded for " << job.year << std::endl;
        fin->Close();
        delete fin;
        return 7;
    }
    enableBranches(tree, isData);
    MyBranch(tree);
    if (!isData && !tree->GetBranch("IDSF_Up_1")) {
        std::cerr << "  [skip] MC file lacks precomputed SF Up/Down branches: " << job.fileName << std::endl;
        fin->Close();
        delete fin;
        return 8;
    }
    tree->SetCacheSize(TREE_CACHE_SIZE);
    tree->SetCacheLearnEntries(10);
    tree->AddBranchToCache("*", kTRUE);
    TH1D* hNWEvts = dynamic_cast<TH1D*>(fin->Get("hNWEvts"));
    if (!hNWEvts) hNWEvts = dynamic_cast<TH1D*>(fin->Get("hNEvts"));
    const double denominator = hNWEvts ? hNWEvts->Integral() : 0.0;
    const double xsw = (!isData && denominator > 0.0) ? job.lumi * XSec(baseName) / denominator : 1.0;
    const Long64_t nEntriesFull = tree->GetEntriesFast();
    const bool isChunked = job.startEntry >= 0 && job.endEntry >= 0;
    const Long64_t loopStart = isChunked ? job.startEntry : 0;
    const Long64_t loopEnd = isChunked ? std::min(job.endEntry, nEntriesFull) : nEntriesFull;
    std::string stem = baseName;
    const size_t dotRoot = stem.rfind(".root");
    if (dotRoot != std::string::npos) stem = stem.substr(0, dotRoot);
    const TString outName = isChunked
        ? Form("%s/hist_%s_chunk%d.root", job.outDir.c_str(), stem.c_str(), job.chunkIndex)
        : Form("%s/hist_%s", job.outDir.c_str(), baseName.c_str());
    TFile* fout = TFile::Open(outName, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "  [skip] cannot create output file: " << outName << std::endl;
        if (fout) fout->Close();
        delete fout;
        fin->Close();
        delete fin;
        return 3;
    }
    fout->SetCompressionSettings(OUTPUT_COMPRESSION);
    fout->cd();
    if (hNWEvts) hNWEvts->Write();
    unordered_map<string,TH1D*> h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT;
    unordered_map<string,TH1D*> h_pt[NlepMax], h_eta[NlepMax], h_phi[NlepMax], h_d0[NlepMax], h_dZ[NlepMax], h_iso[NlepMax];
    h_mZ1.reserve(64);
    h_mZ2.reserve(64);
    h_mH1.reserve(64);
    h_mH2.reserve(64);
    h_zPt.reserve(64);
    h_met.reserve(64);
    h_metphi.reserve(64);
    h_LT.reserve(64);
    for (int i=0; i<NlepMax; ++i) {
        h_pt[i].reserve(64);
        h_eta[i].reserve(64);
        h_phi[i].reserve(64);
        h_d0[i].reserve(64);
        h_dZ[i].reserve(64);
        h_iso[i].reserve(64);
    }
    const std::unordered_set<string> finalStateSet(finalStates.begin(), finalStates.end());
    const Long64_t nEntries = nEntriesFull;
    std::cout << "  entries in tree: " << nEntries << std::endl;
    Long64_t nPassed = 0;
    Long64_t nZPtReweighted = 0;
    Long64_t nZPtNoCandidate = 0;
    Long64_t nZPtInvalid = 0;
    Long64_t nRecoilCorrected = 0;
    Long64_t nRecoilNoCandidate = 0;
    Long64_t nRecoilInvalid = 0;
    for (Long64_t entry=loopStart; entry<loopEnd; ++entry) {
        tree->GetEntry(entry);
        if (entry > loopStart && (entry - loopStart) % 10000000LL == 0) std::cout << "    processed " << (entry - loopStart) << " / " << (loopEnd - loopStart) << std::endl;
        const string catstr = numberToCat(cat);
        if (!finalStateSet.count(catstr)) continue;
        if (!passTightEvent(catstr)) continue;
        if (!passChargeTopology(catstr)) continue;
        if (job.year=="2018" && isData && run>=319077 && applyHEMveto(catstr)=="yes") continue;
        applyTauES(catstr);
        double roccorRelErr[4] = {0.0, 0.0, 0.0, 0.0};
        applyRoccoRCorrection(catstr, isData, *roccor, roccorRelErr);
        // No dedicated systematic is derived from the official MET-xy
        // correction: the JME POG's own correctionlib payload for it
        // (met.json.gz) ships only nominal pt/phi outputs, no up/down/syst
        // nodes -- unlike JES/JER, which always ship a Total/Up/Down source
        // when a systematic is intended. Inventing one (e.g. mirroring the
        // full raw-vs-corrected delta as a +-1sigma envelope) is not the
        // standard convention and produced a wildly oversized band.
        if (APPLY_OFFICIAL_MET_CORRECTION && !applyOfficialMETCorrection(*metCorrections)) continue;
        if (hasDuplicateTightObjects(catstr)) continue;
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
        double weight = computeNominalWeight(catstr, isData, xsw);
        if (job.year=="2018" && !isData && applyHEMveto(catstr)=="yes") weight *= 0.35;
        vector<std::pair<int,int>> OS_pair, SS_pair;
        buildPairs(catstr, OS_pair, SS_pair);
        const EventKinematics kin = buildEventKinematics(catstr, OS_pair, SS_pair);
        vector<WeightSystematic> weightVariants;
        if (APPLY_ZPT_REWEIGHTING && isDY) {
            GenZMatch genZ;
            if (!OS_pair.empty()) {
                const TLorentzVector lep1 = LepV(OS_pair[0].first);
                const TLorentzVector lep2 = LepV(OS_pair[0].second);
                genZ = matchGenZ(lep1.Pt(), lep1.Eta(), lep1.Phi(), lep2.Pt(), lep2.Eta(), lep2.Phi());
            }
            if (genZ.matched) {
                double zPtWeight = 0.0, zPtWeightErr = 0.0;
                if (!getZPtWeight(zPtWeights.get(), genZ.pt, zPtWeight, zPtWeightErr)) { ++nZPtInvalid; continue; }
                weight *= zPtWeight;
                ++nZPtReweighted;
                if (zPtWeight > 0.0) {
                    weightVariants.push_back({"_zptUp", weight / zPtWeight * (zPtWeight + zPtWeightErr)});
                    weightVariants.push_back({"_zptDown", weight / zPtWeight * std::max(0.0, zPtWeight - zPtWeightErr)});
                }
            }
            else ++nZPtNoCandidate;
        }
        if (!isData) {
            appendPileupWeightVariants(weight, weightVariants);
            appendL1PrefiringWeightVariants(weight, weightVariants);
            if (xsecUnc > 0.0) {
                weightVariants.push_back({"_xsecUp", weight * (1.0 + xsecUnc)});
                weightVariants.push_back({"_xsecDown", weight * std::max(0.0, 1.0 - xsecUnc)});
            }
            for (const auto& sfVar : computeObjectSFSystematics(catstr, isData, weight, job.year)) {
                weightVariants.push_back({sfVar.suffix, sfVar.weight});
            }
            const double trigSFNominal = computeTriggerSF(catstr);
            if (trigSFNominal > 0.0) {
                const TrigSFShift trigShift = computeTriggerSFShifts(catstr, job.year);
                weightVariants.push_back({"_eTrigUp", weight / trigSFNominal * trigShift.eUp});
                weightVariants.push_back({"_eTrigDown", weight / trigSFNominal * trigShift.eDown});
                weightVariants.push_back({"_muTrigUp", weight / trigSFNominal * trigShift.muUp});
                weightVariants.push_back({"_muTrigDown", weight / trigSFNominal * trigShift.muDown});
            }
        }
        vector<MetSystematic> metVariants;
        if (hasRecoilVariation) {
            metVariants.push_back({"_recoilResponseUp", recoilShift.responseUpPx, recoilShift.responseUpPy});
            metVariants.push_back({"_recoilResponseDown", recoilShift.responseDownPx, recoilShift.responseDownPy});
            metVariants.push_back({"_recoilResolutionUp", recoilShift.resolutionUpPx, recoilShift.resolutionUpPy});
            metVariants.push_back({"_recoilResolutionDown", recoilShift.resolutionDownPx, recoilShift.resolutionDownPy});
        }
        const vector<string> histKeys = collectHistKeys(catstr, kin, OS_pair);
        for (const auto& key : histKeys) fillAllVars(key, catstr, kin, weight, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
        for (const auto& variant : weightVariants) {
            for (const auto& key : histKeys) fillAllVars(key + variant.suffix, catstr, kin, variant.weight, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
        }
        for (const auto& metVar : metVariants) {
            double shiftedMet, shiftedMetPhi;
            cartesianToPolar(metVar.px, metVar.py, shiftedMet, shiftedMetPhi);
            for (const auto& key : histKeys) fillMetOnly(key + metVar.suffix, shiftedMet, shiftedMetPhi, weight, h_met, h_metphi);
        }
        // Gated on catstr containing a muon (structural, always correct),
        // NOT on whether roccorRelErr happens to be nonzero for one: the
        // shift itself already no-ops safely per-slot when unavailable
        // (factor=1), so gating on per-event availability would silently
        // drop events with an unavailable/zero roccor error from the
        // Up/Down histograms entirely instead of including them unshifted
        // -- exactly the bug found and fixed for TauES below.
        const bool hasRoccorVariation = catstr.find('m') != string::npos;
        if (hasRoccorVariation) {
            const double origPt[4] = {pt_1, pt_2, pt_3, pt_4};
            shiftMuonPts(catstr, origPt, roccorRelErr, 1.0);
            const EventKinematics kinRoccorUp = buildEventKinematics(catstr, OS_pair, SS_pair);
            for (const auto& key : histKeys) fillAllVars(key + "_roccorUp", catstr, kinRoccorUp, weight, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
            shiftMuonPts(catstr, origPt, roccorRelErr, -1.0);
            const EventKinematics kinRoccorDown = buildEventKinematics(catstr, OS_pair, SS_pair);
            for (const auto& key : histKeys) fillAllVars(key + "_roccorDown", catstr, kinRoccorDown, weight, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
            restorePts(origPt);
        }
        double tauESRelUp[4], tauESRelDown[4];
        computeTauESRelShift(catstr, tauESRelUp, tauESRelDown);
        // Gated on catstr containing a tau (structural), NOT on whether a
        // shift was actually available for that tau: shiftTauESPts already
        // no-ops per-slot when relUp/relDown are 0 (e.g. the TauES_Up/Down
        // sentinel case), so this correctly includes every tau-channel
        // event in the Up/Down histograms -- shifted where available,
        // unchanged where not -- instead of the previous data-availability
        // gate, which silently dropped an event from _tauESUp/_tauESDown
        // entirely whenever its tau's decay-mode-aware shift was
        // unavailable. That was collapsing tau-channel _tauESUp/_tauESDown
        // histograms to near-zero relative to nominal.
        const bool hasTauESVariation = catstr.find('t') != string::npos;
        if (hasTauESVariation) {
            const double origPt[4] = {pt_1, pt_2, pt_3, pt_4};
            shiftTauESPts(catstr, origPt, tauESRelUp, 1.0);
            const EventKinematics kinTauESUp = buildEventKinematics(catstr, OS_pair, SS_pair);
            for (const auto& key : histKeys) fillAllVars(key + "_tauESUp", catstr, kinTauESUp, weight, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
            shiftTauESPts(catstr, origPt, tauESRelDown, 1.0);  // relDown already signed correctly
            const EventKinematics kinTauESDown = buildEventKinematics(catstr, OS_pair, SS_pair);
            for (const auto& key : histKeys) fillAllVars(key + "_tauESDown", catstr, kinTauESDown, weight, h_mZ1, h_mZ2, h_mH1, h_mH2, h_zPt, h_met, h_metphi, h_LT, h_pt, h_eta, h_phi, h_d0, h_dZ, h_iso);
            restorePts(origPt);
        }
        ++nPassed;
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
    for (int i=0; i<NlepMax; ++i) {
        writeHists(h_pt[i]);
        writeHists(h_eta[i]);
        writeHists(h_phi[i]);
        writeHists(h_d0[i]);
        writeHists(h_dZ[i]);
        writeHists(h_iso[i]);
    }
    fout->Write();
    fout->Close();
    fin->Close();
    deleteHists(h_mZ1);
    deleteHists(h_mZ2);
    deleteHists(h_mH1);
    deleteHists(h_mH2);
    deleteHists(h_zPt);
    deleteHists(h_met);
    deleteHists(h_metphi);
    deleteHists(h_LT);
    for (int i=0; i<NlepMax; ++i) {
        deleteHists(h_pt[i]);
        deleteHists(h_eta[i]);
        deleteHists(h_phi[i]);
        deleteHists(h_d0[i]);
        deleteHists(h_dZ[i]);
        deleteHists(h_iso[i]);
    }
    delete fout;
    delete fin;
    std::cout << "  [done] passed=" << nPassed;
    if (APPLY_ZPT_REWEIGHTING && isDY) {
        std::cout << " zpt_reweighted=" << nZPtReweighted << " zpt_no_candidate=" << nZPtNoCandidate << " zpt_invalid=" << nZPtInvalid;
    }
    if (useRecoilCorrection) {
        std::cout << " recoil_corrected=" << nRecoilCorrected << " recoil_no_candidate=" << nRecoilNoCandidate << " recoil_invalid=" << nRecoilInvalid;
    }
    std::cout << " written -> " << outName << std::endl;
    return 0;
}

void DCH_tight(string inYear="2018", int firstFile=0, int nFilesToRun=-1, int nProc=1, string processFilter="", bool skipData=false) {
    const TString outDirBase = "/eos/user/a/atahmad/DCH_offline_analysis/new_hists/run2_noFR_metphi_zpt_recoil_roccor_v2";
    const TString outDir = Form("%s/%s", outDirBase.Data(), inYear.c_str());
    gSystem->mkdir(outDir, kTRUE);
    const std::map<string,vector<string>> fileMap = getFileMap(inYear);
    vector<std::pair<string,string>> files;
    for (const auto& sample : fileMap) {
        if (!processFilter.empty() && sample.first != processFilter) continue;
        if (skipData && sample.first == "data") continue;
        for (const string& fileName : sample.second) files.push_back({sample.first, fileName});
    }
    if (files.empty()) { std::cerr << "ERROR: no files found for year " << inYear << std::endl; return; }
    if (firstFile < 0 || firstFile >= static_cast<int>(files.size())) {
        std::cerr << "ERROR: invalid firstFile=" << firstFile << ", number of files=" << files.size() << std::endl;
        return;
    }
    int lastFile = static_cast<int>(files.size());
    if (nFilesToRun >= 0) lastFile = std::min(static_cast<int>(files.size()), firstFile + nFilesToRun);
    const double lumi = inYear=="2016preVFP" ? 19520.0 : inYear=="2016postVFP" ? 16810.0 : inYear=="2017" ? 41480.0 : inYear=="2018" ? 59830.0 : 19520.0 + 16810.0 + 41480.0 + 59830.0;
    vector<FileJob> jobs;
    jobs.reserve(lastFile - firstFile);
    for (int fileIndex=firstFile; fileIndex<lastFile; ++fileIndex) {
        FileJob job;
        job.year = inYear;
        job.process = files[fileIndex].first;
        job.fileName = files[fileIndex].second;
        job.outDir = outDir.Data();
        job.lumi = lumi;
        job.index = fileIndex;
        job.total = static_cast<int>(files.size());
        jobs.push_back(job);
    }
    if (jobs.empty()) { std::cout << "No files selected." << std::endl; return; }
    if (nProc < 1) nProc = 1;
    nProc = std::min(nProc, static_cast<int>(jobs.size()));
    std::cout << "Selected " << jobs.size() << " files, using " << nProc << " process" << (nProc==1 ? "" : "es") << "." << std::endl;
    if (nProc == 1) {
        for (const FileJob& job : jobs) ProcessTightFile(job);
    }
    else {
        ROOT::TProcessExecutor pool(nProc);
        const auto results = pool.Map(ProcessTightFile, jobs);
        int nFailed = 0;
        for (const int result : results) if (result != 0) ++nFailed;
        if (nFailed > 0) std::cerr << "WARNING: " << nFailed << " file jobs returned a nonzero status." << std::endl;
    }
}

// Processes exactly one 1M-event chunk of one (year, file-index) file,
// mirroring the one-job-per-file condor convention already used by
// run_dch_perfile.sh, just with an added chunk index. Chunk outputs are
// staged in a chunks/ subdirectory of the normal per-year output dir and
// are meant to be combined back into hist_<sample>.root by
// MergeChunkHists.C before Stackhist.C ever sees them.
void DCH_tight_Chunk(string inYear="2018", int idx=0, int chunkIdx=0) {
    const TString outDir = Form("/eos/user/a/atahmad/DCH_offline_analysis/new_hists/chunks/run2_noFR_metphi_zpt_recoil_roccor_v2/%s", inYear.c_str());
    gSystem->mkdir(outDir, kTRUE);
    const std::map<string,vector<string>> fileMap = getFileMap(inYear);
    vector<std::pair<string,string>> files;
    for (const auto& sample : fileMap)
        for (const string& fileName : sample.second) files.push_back({sample.first, fileName});
    if (idx < 0 || idx >= static_cast<int>(files.size())) {
        std::cerr << "ERROR: invalid idx=" << idx << ", number of files=" << files.size() << std::endl;
        return;
    }
    const double lumi = inYear=="2016preVFP" ? 19520.0 : inYear=="2016postVFP" ? 16810.0 : inYear=="2017" ? 41480.0 : inYear=="2018" ? 59830.0 : 19520.0 + 16810.0 + 41480.0 + 59830.0;
    FileJob job;
    job.year = inYear;
    job.process = files[idx].first;
    job.fileName = files[idx].second;
    job.outDir = outDir.Data();
    job.lumi = lumi;
    job.index = idx;
    job.total = static_cast<int>(files.size());
    job.startEntry = static_cast<Long64_t>(chunkIdx) * CHUNK_SIZE;
    job.endEntry = job.startEntry + CHUNK_SIZE;
    job.chunkIndex = chunkIdx;
    ProcessTightFile(job);
}
