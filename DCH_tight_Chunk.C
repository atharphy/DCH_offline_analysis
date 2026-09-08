// DCH_tight_Chunk.C
//
// Processes exactly one 1M-event chunk of one (year, file-index) file from
// DCH_tight.C's own file list, mirroring run_dch_perfile.sh's one-job-per-file
// condor convention with an added chunk index. Chunk outputs are staged in a
// chunks/ subdirectory of the normal per-year output dir and are meant to be
// combined back into hist_<sample>.root by MergeChunkHists.C before
// Stackhist.C ever sees them.
//
// Run: root -l -b -q 'DCH_tight_Chunk.C+("2018",0,0)'

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

using std::string;
using std::vector;

void DCH_tight_Chunk(string inYear="2018", int idx=0, int chunkIdx=0) {
    const TString outDir = Form("/eos/user/a/atahmad/DCH_offline_analysis/new_hists/chunks/run2_noFR_metphi_zpt_recoil_roccor_toppt_syst_v2/%s", inYear.c_str());
    gSystem->mkdir(outDir, kTRUE);
    const TString treeOutDir = Form("/eos/user/a/atahmad/DCH_offline_analysis/new_trees/chunks/run2_noFR_metphi_zpt_recoil_roccor_toppt_syst_v2/%s", inYear.c_str());
    gSystem->mkdir(treeOutDir, kTRUE);

    const std::map<string, vector<string>> fileMap = getFileMap(inYear);
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
    job.treeOutDir = treeOutDir.Data();
    job.lumi = lumi;
    job.index = idx;
    job.total = static_cast<int>(files.size());
    job.startEntry = static_cast<Long64_t>(chunkIdx) * CHUNK_SIZE;
    job.endEntry = job.startEntry + CHUNK_SIZE;
    job.chunkIndex = chunkIdx;
    ProcessTightFile(job);
}
