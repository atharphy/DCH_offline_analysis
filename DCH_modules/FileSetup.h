#pragma once

#include <memory>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "BranchSetup.h"
#include "CommonConfig.h"
#include "FileJob.h"
#include "FlatXsecSystematic.h"
#include "MEtSysWrapper.h"
#include "METCorrections.h"
#include "RecoilCorrections.h"
#include "RoccoRCorrections.h"
#include "ZPtReweight.h"

struct FileContext {
    bool isData = false, isDY = false, isWJ = false, isTTbar = false;
    double xsw = 1.0, xsecUnc = 0.0;
    bool useRecoilCorrection = false;

    std::unique_ptr<TH1D> zPtWeights;
    std::unique_ptr<RecoilCorrector> recoilCorrector;
    std::shared_ptr<MEtSys> metSys;
    std::unique_ptr<OfficialMETCorrections> metCorrections;
    std::unique_ptr<RoccoR> roccor;
};

inline bool setupFileContext(const FileJob& job, TFile* fin, TTree* tree, const std::string& baseName, FileContext& ctx, int& errorCode) {
    ctx.isData = XSec(baseName) == 1;
    ctx.isDY = job.process == "DY" || job.process == "DY10_50";
    ctx.isWJ = job.process == "WJ";
    ctx.isTTbar = job.process == "TTbar";
    ctx.xsecUnc = getFlatXsecUncertainty(baseName);

    if (APPLY_ZPT_REWEIGHTING && ctx.isDY) {
        ctx.zPtWeights.reset(loadZPtWeights(job.year));
        if (!ctx.zPtWeights) { errorCode = 5; return false; }
    }

    ctx.useRecoilCorrection = !ctx.isData && ((APPLY_DY_RECOIL_CORRECTION && ctx.isDY) || (APPLY_WJ_RECOIL_CORRECTION && ctx.isWJ));
    if (ctx.useRecoilCorrection) {
        ctx.recoilCorrector = loadRecoilCorrector(job.year);
        if (!ctx.recoilCorrector) { errorCode = 6; return false; }
        ctx.metSys = getCachedMEtSys(job.year);
        if (!ctx.metSys) std::cerr << "  [warn] MEtSys payload not found for " << job.year << "; recoil systematic variations will be skipped" << std::endl;
    }

    if (APPLY_OFFICIAL_MET_CORRECTION) {
        ctx.metCorrections = loadOfficialMETCorrections(job.year, ctx.isData);
        if (!ctx.metCorrections) { errorCode = 4; return false; }
    }

    ctx.roccor = loadRoccoRCorrections(job.year);
    if (!ctx.roccor) { errorCode = 7; return false; }

    enableBranches(tree, ctx.isData);
    MyBranch(tree);
    if (!ctx.isData && !tree->GetBranch("IDSF_Up_1")) { errorCode = 8; return false; }

    TH1D* hNWEvts = dynamic_cast<TH1D*>(fin->Get("hNWEvts"));
    if (!hNWEvts) hNWEvts = dynamic_cast<TH1D*>(fin->Get("hNEvts"));
    const double denominator = hNWEvts ? hNWEvts->Integral() : 0.0;
    ctx.xsw = (!ctx.isData && denominator > 0.0) ? job.lumi * XSec(baseName) / denominator : 1.0;
    return true;
}
