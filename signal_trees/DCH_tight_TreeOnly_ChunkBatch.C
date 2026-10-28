#if !defined(__CLING__)
#pragma GCC optimize("O3,unroll-loops")
#endif

#include "correction.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#if defined(__CLING__)
R__LOAD_LIBRARY(libcorrectionlib.so)
R__LOAD_LIBRARY(libHTT-utilitiesRecoilCorrections.so)
#endif

#include "filemap/FileMapHTQCD.h"
#include "include/MyBranch.C"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#pragma GCC diagnostic ignored "-Wunused-value"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wuninitialized"
#include "include/Kinematics.C"
#pragma GCC diagnostic pop

#include "include/Xsections.C"

#include "DCH_modules/TightTreeOnlyFileProcessor.h"
#include "MergeChunkHists.C"

using std::string;
using std::vector;

void DCH_tight_TreeOnly_ChunkBatch(string inYear = "2018", int idx = 0, int batchStart = 0, int numInBatch = 15) {
    const char* scratchEnv = gSystem->Getenv("_CONDOR_SCRATCH_DIR");
    const std::string scratch = (scratchEnv && scratchEnv[0]) ? scratchEnv : ".";
    const std::string localTreeDir = scratch + "/tree_chunks";
    gSystem->mkdir(localTreeDir.c_str(), kTRUE);

    const char* fmEnv = gSystem->Getenv("DCH_FILEMAP");
    const bool qcdHT = fmEnv && std::string(fmEnv) == "qcdht";
    const TString treeChunkDir = Form("output/trees/chunks/%s/%s",
        qcdHT ? "run2_noFR_metphi_zpt_recoil_roccor_toppt_chargeflip_LT_met_qcdht" : "run2_noFR_metphi_zpt_recoil_roccor_toppt_chargeflip_LT", inYear.c_str());
    gSystem->mkdir(treeChunkDir, kTRUE);

    const std::map<string, vector<string>> fileMap = qcdHT ? getFileMapHTQCD(inYear) : getFileMap(inYear);
    vector<std::pair<string,string>> files;
    for (const auto& sample : fileMap)
        for (const string& fileName : sample.second) files.push_back({sample.first, fileName});
    if (idx < 0 || idx >= static_cast<int>(files.size())) {
        std::cerr << "ERROR: invalid idx=" << idx << ", number of files=" << files.size() << std::endl;
        return;
    }

    const double lumi = inYear=="2016preVFP" ? 19520.0 : inYear=="2016postVFP" ? 16810.0 : inYear=="2017" ? 41480.0 : inYear=="2018" ? 59830.0 : 19520.0 + 16810.0 + 41480.0 + 59830.0;

    for (int c = batchStart; c < batchStart + numInBatch; ++c) {
        FileJob job;
        job.year = inYear;
        job.process = files[idx].first;
        job.fileName = files[idx].second;
        job.treeOutDir = localTreeDir;
        job.lumi = lumi;
        job.index = idx;
        job.total = static_cast<int>(files.size());
        job.startEntry = static_cast<Long64_t>(c) * TREEONLY_CHUNK_SIZE;
        job.endEntry = job.startEntry + TREEONLY_CHUNK_SIZE;
        job.chunkIndex = c;
        const int rc = ProcessTightFileTreeOnly(job);
        if (rc != 0) { std::cerr << "  [batch] chunk " << c << " failed (code " << rc << "), stopping batch" << std::endl; break; }
    }

    const std::string localMergedTreeDir = scratch + "/merged_tree";
    gSystem->mkdir(localMergedTreeDir.c_str(), kTRUE);
    MergeChunkHists(localTreeDir, localMergedTreeDir);

    TSystemDirectory lister("d", localMergedTreeDir.c_str());
    TList* mergedFiles = lister.GetListOfFiles();
    if (!mergedFiles) return;
    TIter next(mergedFiles);
    TSystemFile* sf;
    while ((sf = (TSystemFile*)next())) {
        TString n = sf->GetName();
        if (!n.EndsWith(".root")) continue;
        TString tagged = n; tagged.ReplaceAll(".root", Form("_chunk%d.root", batchStart));
        TFile::Cp(Form("%s/%s", localMergedTreeDir.c_str(), n.Data()), treeChunkDir + "/" + tagged, kFALSE);
    }
}
