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
// Chunked per-year batch processing lives in DCH_tight_Chunk.C; merging
// chunk output back into per-sample files is MergeChunkHists.C.

#if !defined(__CLING__)
#pragma GCC optimize("O3,unroll-loops")
#endif

#include "correction.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "ROOT/TProcessExecutor.hxx"

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

void DCH_tight(string inYear="2018", int firstFile=0, int nFilesToRun=-1, int nProc=1, string processFilter="", bool skipData=false) {
    const TString outDirBase = "/eos/user/a/atahmad/DCH_offline_analysis/new_hists/run2_noFR_metphi_zpt_recoil_roccor_toppt_syst_v2";
    const TString outDir = Form("%s/%s", outDirBase.Data(), inYear.c_str());
    gSystem->mkdir(outDir, kTRUE);

    const std::map<string, vector<string>> fileMap = getFileMap(inYear);
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

    nProc = std::max(1, std::min(nProc, static_cast<int>(jobs.size())));
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
