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
#include "DCH_modules/TruthMatching.h"
#include "DCH_modules/ZPtReweight.h"
#include "DCH_modules/FileJob.h"
#include "DCH_modules/EventWeights.h"
#include "DCH_modules/PairBuilder.h"
#include "DCH_modules/MEtSysWrapper.h"
#include "DCH_modules/SystematicPlan.h"
#include "DCH_modules/FlatXsecSystematic.h"

using std::string;
using std::unordered_map;
using std::vector;

const int OUTPUT_COMPRESSION = ROOT::CompressionSettings(ROOT::kZSTD, 1);
const Long64_t TREE_CACHE_SIZE = 200LL * 1024LL * 1024LL;

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
            "TauVsJetIDSF_1","TauVsJetIDSF_2","TauVsJetIDSF_3","TauVsJetIDSF_4"
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
    tree->SetCacheSize(TREE_CACHE_SIZE);
    tree->SetCacheLearnEntries(10);
    tree->AddBranchToCache("*", kTRUE);
    TH1D* hNWEvts = dynamic_cast<TH1D*>(fin->Get("hNWEvts"));
    if (!hNWEvts) hNWEvts = dynamic_cast<TH1D*>(fin->Get("hNEvts"));
    const double denominator = hNWEvts ? hNWEvts->Integral() : 0.0;
    const double xsw = (!isData && denominator > 0.0) ? job.lumi * XSec(baseName) / denominator : 1.0;
    const TString outName = Form("%s/hist_%s", job.outDir.c_str(), baseName.c_str());
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
    const Long64_t nEntries = tree->GetEntriesFast();
    std::cout << "  entries in tree: " << nEntries << std::endl;
    Long64_t nPassed = 0;
    Long64_t nZPtReweighted = 0;
    Long64_t nZPtNoCandidate = 0;
    Long64_t nZPtInvalid = 0;
    Long64_t nRecoilCorrected = 0;
    Long64_t nRecoilNoCandidate = 0;
    Long64_t nRecoilInvalid = 0;
    for (Long64_t entry=0; entry<nEntries; ++entry) {
        tree->GetEntry(entry);
        if (entry > 0 && entry % 10000000LL == 0) std::cout << "    processed " << entry << " / " << nEntries << std::endl;
        const string catstr = numberToCat(cat);
        if (!finalStateSet.count(catstr)) continue;
        if (!passTightEvent(catstr)) continue;
        if (!passChargeTopology(catstr)) continue;
        if (job.year=="2018" && isData && run>=319077 && applyHEMveto(catstr)=="yes") continue;
        applyTauES(catstr);
        double roccorRelErr[4] = {0.0, 0.0, 0.0, 0.0};
        applyRoccoRCorrection(catstr, isData, *roccor, roccorRelErr);
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
        if (hasDuplicateTightObjects(catstr)) continue;
        RecoilSystShift recoilShift;
        bool hasRecoilVariation = false;
        if (useRecoilCorrection) {
            double genPx = 0.0;
            double genPy = 0.0;
            double visPx = 0.0;
            double visPy = 0.0;
            bool hasCandidate = false;
            if (isDY) {
                hasCandidate = getPromptTruthDileptonXY(catstr, genPx, genPy);
                visPx = genPx;
                visPy = genPy;
            }
            else if (isWJ) {
                hasCandidate = getPromptTruthLeptonXY(catstr, visPx, visPy);
                if (hasCandidate) {
                    const double metPx = met * std::cos(metphi);
                    const double metPy = met * std::sin(metphi);
                    genPx = visPx + metPx;
                    genPy = visPy + metPy;
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
        double weight = computeNominalWeight(catstr, isData, xsw);
        if (job.year=="2018" && !isData && applyHEMveto(catstr)=="yes") weight *= 0.35;
        vector<std::pair<int,int>> OS_pair, SS_pair;
        buildPairs(catstr, OS_pair, SS_pair);
        const EventKinematics kin = buildEventKinematics(catstr, OS_pair, SS_pair);
        vector<WeightSystematic> weightVariants;
        if (APPLY_ZPT_REWEIGHTING && isDY) {
            double promptTruthZPt = 0.0;
            if (getPromptTruthDileptonPt(catstr, promptTruthZPt)) {
                double zPtWeight = 0.0, zPtWeightErr = 0.0;
                if (!getZPtWeight(zPtWeights.get(), promptTruthZPt, zPtWeight, zPtWeightErr)) { ++nZPtInvalid; continue; }
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
        bool hasRoccorVariation = false;
        for (int idx = 1; idx <= (int)catstr.size() && idx <= 4; ++idx) {
            if (catstr[idx - 1] == 'm' && roccorRelErr[idx - 1] > 0.0) hasRoccorVariation = true;
        }
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

void DCH_tight(string inYear="2018", int firstFile=0, int nFilesToRun=-1, int nProc=1, string processFilter="") {
    const TString outDirBase = "hists/run2_hists_noFR_reccor";
    const TString outDir = Form("%s/%s", outDirBase.Data(), inYear.c_str());
    gSystem->mkdir(outDir, kTRUE);
    const std::map<string,vector<string>> fileMap = getFileMap(inYear);
    vector<std::pair<string,string>> files;
    for (const auto& sample : fileMap) {
        if (!processFilter.empty() && sample.first != processFilter) continue;
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
