#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Compression.h"
#include "TFile.h"
#include "TTree.h"

#include "CommonConfig.h"
#include "EventSelection.h"
#include "EventWeights.h"
#include "FileJob.h"
#include "FileSetup.h"
#include "GenBosonMomentum.h"
#include "GenZMatching.h"
#include "MEtSysWrapper.h"
#include "METCorrections.h"
#include "PairBuilder.h"
#include "RecoilCorrections.h"
#include "RoccoRCorrections.h"
#include "../Mass_reco_modules/SignalShapeSelection.h"
#include "SystematicPlan.h"
#include "TauESSystematic.h"
#include "ZPtReweight.h"

const int TREEONLY_OUTPUT_COMPRESSION = ROOT::CompressionSettings(ROOT::kZSTD, 1);
const Long64_t TREEONLY_TREE_CACHE_SIZE = 200LL * 1024LL * 1024LL;
const Long64_t TREEONLY_CHUNK_SIZE = 500000LL;

inline int ProcessTightFileTreeOnly(FileJob job) {
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

    const std::string baseName = gSystem->BaseName(job.fileName.c_str());
    FileContext ctx;
    int errorCode = 0;
    if (!setupFileContext(job, fin, tree, baseName, ctx, errorCode)) {
        std::cerr << "  [skip] setup failed (code " << errorCode << ") for " << job.fileName << std::endl;
        fin->Close();
        delete fin;
        return errorCode;
    }

    tree->SetCacheSize(TREEONLY_TREE_CACHE_SIZE);
    tree->SetCacheLearnEntries(10);
    tree->AddBranchToCache("*", kTRUE);

    const Long64_t nEntriesFull = tree->GetEntriesFast();
    const bool isChunked = job.startEntry >= 0 && job.endEntry >= 0;
    const Long64_t loopStart = isChunked ? job.startEntry : 0;
    const Long64_t loopEnd = isChunked ? std::min(job.endEntry, nEntriesFull) : nEntriesFull;

    std::string stem = baseName;
    const size_t dotRoot = stem.rfind(".root");
    if (dotRoot != std::string::npos) stem = stem.substr(0, dotRoot);
    const TString outNameMass = isChunked
        ? Form("%s/masstree_%s_chunk%d.root", job.treeOutDir.c_str(), stem.c_str(), job.chunkIndex)
        : Form("%s/masstree_%s", job.treeOutDir.c_str(), baseName.c_str());
    TFile* foutMass = TFile::Open(outNameMass, "RECREATE");
    if (!foutMass || foutMass->IsZombie()) {
        std::cerr << "  [skip] cannot create mass-tree output file: " << outNameMass << std::endl;
        if (foutMass) foutMass->Close();
        delete foutMass;
        fin->Close();
        delete fin;
        return 3;
    }
    foutMass->SetCompressionSettings(TREEONLY_OUTPUT_COMPRESSION);

    SignalShapeHistBundle sigHists;
    sigHists.reserveAll(4);

    const std::unordered_set<std::string> finalStateSet(finalStates.begin(), finalStates.end());
    std::cout << "  entries in tree: " << nEntriesFull << std::endl;
    Long64_t nPassed = 0;

    for (Long64_t entry = loopStart; entry < loopEnd; ++entry) {
        tree->GetEntry(entry);
        if (entry > loopStart && (entry - loopStart) % 10000000LL == 0)
            std::cout << "    processed " << (entry - loopStart) << " / " << (loopEnd - loopStart) << std::endl;

        const std::string catstr = numberToCat(cat);
        if (!finalStateSet.count(catstr)) continue;
        if (!passTightEvent(catstr)) continue;
        if (!passChargeTopology(catstr)) continue;
        if (job.year == "2018" && ctx.isData && run >= 319077 && applyHEMveto(catstr) == "yes") continue;
        if (catstr.size() != 3 && catstr.size() != 4) continue;

        applyTauES(catstr);
        double roccorRelErr[4] = {0.0, 0.0, 0.0, 0.0};
        applyRoccoRCorrection(catstr, ctx.isData, *ctx.roccor, roccorRelErr);
        std::vector<ChargeFlipDecision> chargeFlipDecisions;
        if (!ctx.isData) chargeFlipDecisions = applyChargeFlipCorrection(catstr, ctx.isData, *ctx.chargeFlipRate, ctx.chargeFlipRng);
        const double chargeFlipRatioUp = ctx.isData ? 1.0 : chargeFlipWeightRatio(chargeFlipDecisions, *ctx.chargeFlipRateUp);
        const double chargeFlipRatioDown = ctx.isData ? 1.0 : chargeFlipWeightRatio(chargeFlipDecisions, *ctx.chargeFlipRateDown);

        if (APPLY_OFFICIAL_MET_CORRECTION && !applyOfficialMETCorrection(*ctx.metCorrections)) continue;
        if (hasDuplicateTightObjects(catstr)) continue;

        RecoilSystShift recoilShift;
        bool hasRecoilVariation = false;
        if (ctx.useRecoilCorrection) {
            double genPx = 0.0, genPy = 0.0, visPx = 0.0, visPy = 0.0;
            const bool hasCandidate = (ctx.isDY || ctx.isWJ) && getGenBosonMomentum(genPx, genPy, visPx, visPy);
            if (hasCandidate && applyRecoilCorrection(*ctx.recoilCorrector, genPx, genPy, visPx, visPy) && ctx.metSys) {
                const int recoilNJets = std::max(0, static_cast<int>(std::lround(njets)));
                const float correctedPx = static_cast<float>(met * std::cos(metphi));
                const float correctedPy = static_cast<float>(met * std::sin(metphi));
                recoilShift = computeRecoilSystShift(*ctx.metSys, correctedPx, correctedPy,
                                                      static_cast<float>(genPx), static_cast<float>(genPy),
                                                      static_cast<float>(visPx), static_cast<float>(visPy), recoilNJets);
                hasRecoilVariation = recoilShift.valid;
            }
        }

        double weight = computeNominalWeight(catstr, ctx.isData, ctx.isTTbar, ctx.xsw);
        if (job.year == "2018" && !ctx.isData && applyHEMveto(catstr) == "yes") weight *= 0.35;

        std::vector<std::pair<int,int>> OS_pair, SS_pair;
        buildPairs(catstr, OS_pair, SS_pair);

        std::vector<WeightSystematic> weightVariants;
        if (APPLY_ZPT_REWEIGHTING && ctx.isDY) {
            GenZMatch genZ;
            if (!OS_pair.empty()) {
                const TLorentzVector lep1 = LepV(OS_pair[0].first);
                const TLorentzVector lep2 = LepV(OS_pair[0].second);
                genZ = matchGenZ(lep1.Pt(), lep1.Eta(), lep1.Phi(), lep2.Pt(), lep2.Eta(), lep2.Phi());
            }
            if (genZ.matched) {
                double zPtWeight = 0.0, zPtWeightErr = 0.0;
                if (!getZPtWeight(ctx.zPtWeights.get(), genZ.pt, zPtWeight, zPtWeightErr)) continue;
                weight *= zPtWeight;
                if (zPtWeight > 0.0) {
                    weightVariants.push_back({"_zptUp", weight / zPtWeight * (zPtWeight + zPtWeightErr)});
                    weightVariants.push_back({"_zptDown", weight / zPtWeight * std::max(0.0, zPtWeight - zPtWeightErr)});
                }
            }
        }
        appendStandardWeightVariants(catstr, ctx.isData, ctx.isTTbar, ctx.xsecUnc, weight, job.year, weightVariants);
        if (!ctx.isData && !chargeFlipDecisions.empty()) {
            weightVariants.push_back({"_chargeFlipUp", weight * chargeFlipRatioUp});
            weightVariants.push_back({"_chargeFlipDown", weight * chargeFlipRatioDown});
        }

        const SignalShapeVariables ssNominal = computeSignalShapeVariables(catstr, cat, SS_pair);
        fillSignalShapeVariant(ssNominal, weight, "", sigHists);
        for (const auto& variant : weightVariants) fillSignalShapeVariant(ssNominal, variant.weight, variant.suffix, sigHists);
        if (hasRecoilVariation) {
            const std::pair<const char*, std::pair<double,double>> recoilVars[4] = {
                {"_recoilResponseUp",    {recoilShift.responseUpPx,   recoilShift.responseUpPy}},
                {"_recoilResponseDown",  {recoilShift.responseDownPx, recoilShift.responseDownPy}},
                {"_recoilResolutionUp",  {recoilShift.resolutionUpPx, recoilShift.resolutionUpPy}},
                {"_recoilResolutionDown",{recoilShift.resolutionDownPx, recoilShift.resolutionDownPy}}};
            for (const auto& rv : recoilVars)
                fillSignalShapeVariant(ssNominal, weight, rv.first, sigHists,
                                       std::hypot(rv.second.first, rv.second.second), std::atan2(rv.second.second, rv.second.first));
        }

        if (catstr.find('m') != std::string::npos) {
            const double origPt[4] = {pt_1, pt_2, pt_3, pt_4};
            shiftMuonPts(catstr, origPt, roccorRelErr, 1.0);
            fillSignalShapeVariant(computeSignalShapeVariables(catstr, cat, SS_pair), weight, "_roccorUp", sigHists);
            shiftMuonPts(catstr, origPt, roccorRelErr, -1.0);
            fillSignalShapeVariant(computeSignalShapeVariables(catstr, cat, SS_pair), weight, "_roccorDown", sigHists);
            restorePts(origPt);
        }

        double tauESRelUp[4], tauESRelDown[4];
        computeTauESRelShift(catstr, tauESRelUp, tauESRelDown);
        if (catstr.find('t') != std::string::npos) {
            const double origPt[4] = {pt_1, pt_2, pt_3, pt_4};
            shiftTauESPts(catstr, origPt, tauESRelUp, 1.0);
            fillSignalShapeVariant(computeSignalShapeVariables(catstr, cat, SS_pair), weight, "_tauESUp", sigHists);
            shiftTauESPts(catstr, origPt, tauESRelDown, 1.0);
            fillSignalShapeVariant(computeSignalShapeVariables(catstr, cat, SS_pair), weight, "_tauESDown", sigHists);
            restorePts(origPt);
        }

        ++nPassed;
    }

    foutMass->cd();
    sigHists.writeAll();
    foutMass->Write();
    foutMass->Close();
    delete foutMass;
    sigHists.deleteAll();

    fin->Close();
    delete fin;

    std::cout << "  [done] passed=" << nPassed << " written -> " << outNameMass << std::endl;
    return 0;
}
