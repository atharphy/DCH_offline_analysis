// DCH_tight_ChunkBatch.C
//
// Processes a small BATCH of consecutive 500k-event chunks of one (year,
// file-index) file, writing each chunk's intermediate output to the job's
// own local condor scratch directory (never touching EOS), merging that
// batch locally, then uploading only the one small merged-batch file to
// EOS's normal chunks/ staging directory -- using the exact same
// "hist_<sample>_chunk<N>.root" naming MergeChunkHists.C already expects,
// tagged by the batch's own starting chunk index. No change to
// MergeChunkHists.C or the final per-sample merge step is needed: they
// just see fewer, larger "chunk" files.
//
// This exists because staging all ~6550 raw chunks directly on EOS (for
// both the main histogram output and the mass-reco tree output) peaked at
// well over 100GB before the final merge shrank it back down to <20GB --
// a real risk against a ~265GB EOS headroom. Batching locally on
// per-job scratch cuts that peak EOS footprint by roughly the batch size.
//
// Run: root -l -b -q 'DCH_tight_ChunkBatch.C+("2018",0,0,15)'

#if !defined(__CLING__)
#pragma GCC optimize("O3,unroll-loops")
#endif

#include "correction.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#if defined(__CLING__)
R__LOAD_LIBRARY(libcorrectionlib.so)
R__LOAD_LIBRARY(libHTT-utilitiesRecoilCorrections.so)
#endif

#include "filemap/FileMap.h"
#include "include/MyBranch.C"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#pragma GCC diagnostic ignored "-Wunused-value"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wuninitialized"
#include "include/Kinematics.C"
#pragma GCC diagnostic pop

#include "include/Xsections.C"

#include "DCH_modules/TightFileProcessor.h"
#include "MergeChunkHists.C"

using std::string;
using std::vector;

void DCH_tight_ChunkBatch(string inYear = "2018", int idx = 0, int batchStart = 0, int numInBatch = 15) {
    const char* scratchEnv = gSystem->Getenv("_CONDOR_SCRATCH_DIR");
    const std::string scratch = (scratchEnv && scratchEnv[0]) ? scratchEnv : ".";
    const std::string localHistDir = scratch + "/hist_chunks";
    const std::string localTreeDir = scratch + "/tree_chunks";
    gSystem->mkdir(localHistDir.c_str(), kTRUE);
    gSystem->mkdir(localTreeDir.c_str(), kTRUE);

    const TString eosHistChunkDir = Form("/eos/user/a/atahmad/DCH_offline_analysis/new_hists/chunks/run2_noFR_metphi_zpt_recoil_roccor_toppt_syst_v2/%s", inYear.c_str());
    const TString eosTreeChunkDir = Form("/eos/user/a/atahmad/DCH_offline_analysis/new_trees/chunks/run2_noFR_metphi_zpt_recoil_roccor_toppt_syst_v2/%s", inYear.c_str());
    gSystem->mkdir(eosHistChunkDir, kTRUE);
    gSystem->mkdir(eosTreeChunkDir, kTRUE);

    const std::map<string, vector<string>> fileMap = getFileMap(inYear);
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
        job.outDir = localHistDir;
        job.treeOutDir = localTreeDir;
        job.lumi = lumi;
        job.index = idx;
        job.total = static_cast<int>(files.size());
        job.startEntry = static_cast<Long64_t>(c) * CHUNK_SIZE;
        job.endEntry = job.startEntry + CHUNK_SIZE;
        job.chunkIndex = c;
        const int rc = ProcessTightFile(job);
        // A chunk past the end of the file's real entry count still opens
        // fine and just processes zero events (TTree::GetEntry on an
        // out-of-range index is a no-op) -- only a genuine open/setup
        // failure returns non-zero, and once one chunk in the batch hits
        // that, every later one in the same batch would too (same file).
        if (rc != 0) { std::cerr << "  [batch] chunk " << c << " failed (code " << rc << "), stopping batch" << std::endl; break; }
    }

    // Merge the batch's chunks LOCALLY first (target is still scratch, not
    // EOS) -- with many batch jobs for the same sample running at once,
    // merging straight into a shared EOS "hist_<sample>.root" name (then
    // renaming) would race: two jobs could momentarily write the same
    // filename at once. Merging locally into a batch-private directory and
    // uploading the single already-uniquely-named result sidesteps that
    // entirely, since the EOS destination filename is unique per batch
    // from the first byte written.
    const std::string localMergedHistDir = scratch + "/merged_hist";
    const std::string localMergedTreeDir = scratch + "/merged_tree";
    gSystem->mkdir(localMergedHistDir.c_str(), kTRUE);
    gSystem->mkdir(localMergedTreeDir.c_str(), kTRUE);
    MergeChunkHists(localHistDir, localMergedHistDir);
    MergeChunkHists(localTreeDir, localMergedTreeDir);

    auto uploadTagged = [&](const std::string& localMergedDir, const TString& eosChunkDir) {
        TSystemDirectory lister("d", localMergedDir.c_str());
        TList* mergedFiles = lister.GetListOfFiles();
        if (!mergedFiles) return;
        TIter next(mergedFiles);
        TSystemFile* sf;
        while ((sf = (TSystemFile*)next())) {
            TString n = sf->GetName();
            if (!n.EndsWith(".root")) continue;
            TString tagged = n; tagged.ReplaceAll(".root", Form("_chunk%d.root", batchStart));
            TFile::Cp(Form("%s/%s", localMergedDir.c_str(), n.Data()), eosChunkDir + "/" + tagged, kFALSE);
        }
    };
    uploadTagged(localMergedHistDir, eosHistChunkDir);
    uploadTagged(localMergedTreeDir, eosTreeChunkDir);
}
