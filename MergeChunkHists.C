// Merges per-chunk histogram files (hist_<sample>_chunk<N>.root) living in
// chunkDir back into hist_<sample>.root written into targetDir, matching
// exactly the file each driver (DCH_tight.C / DCH_tauFR.C) would have
// produced unchunked -- Stackhist.C needs zero changes downstream.
//
// hNWEvts/hNEvts is copied once (not summed): every chunk of the same
// source file carries an identical copy of that file's full normalization
// histogram (it's read straight from the source file, not restricted to
// the chunk's entry range), so summing it across chunks would multiply the
// normalization by the chunk count.
//
// After a successful merge, the chunk files for that sample are deleted.
//
// Run from Updated_offline_framework/:
//   root -l -b -q 'MergeChunkHists.C("/eos/user/a/atahmad/DCH_offline_analysis/new_hists/chunks/run2_noFR_metphi_zpt_recoil_roccor/2018", "/eos/user/a/atahmad/DCH_offline_analysis/new_hists/run2_noFR_metphi_zpt_recoil_roccor/2018")'

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

void MergeChunkHists(std::string chunkDir, std::string targetDir) {
    gSystem->mkdir(targetDir.c_str(), kTRUE);
    TSystemDirectory dir("chunks", chunkDir.c_str());
    TList* entries = dir.GetListOfFiles();
    if (!entries) { std::cerr << "ERROR: cannot list " << chunkDir << std::endl; return; }

    std::map<std::string, std::vector<std::string>> groups;  // base sample name -> chunk file paths
    std::regex chunkPattern("^(.*)_chunk[0-9]+\\.root$");
    TIter next(entries);
    TObject* obj;
    while ((obj = next())) {
        std::string name = obj->GetName();
        if (name.size() < 5 || name.substr(name.size() - 5) != ".root") continue;
        std::smatch m;
        if (!std::regex_match(name, m, chunkPattern)) continue;
        std::string base = m[1].str() + ".root";
        groups[base].push_back(chunkDir + "/" + name);
    }
    delete entries;

    std::cout << "Found " << groups.size() << " sample(s) to merge in " << chunkDir << std::endl;

    for (auto& kv : groups) {
        const std::string& base = kv.first;
        std::vector<std::string>& chunkFiles = kv.second;
        std::sort(chunkFiles.begin(), chunkFiles.end());
        std::cout << "  merging " << chunkFiles.size() << " chunk(s) -> " << base << std::endl;

        std::map<std::string, TH1*> merged;
        TH1* normHist = nullptr;

        for (const auto& path : chunkFiles) {
            TFile* f = TFile::Open(path.c_str(), "READ");
            if (!f || f->IsZombie()) { std::cerr << "    [warn] cannot open " << path << std::endl; continue; }
            TIter keyIter(f->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)keyIter())) {
                TObject* keyObj = key->ReadObj();
                TH1* h = dynamic_cast<TH1*>(keyObj);
                if (!h) { delete keyObj; continue; }
                std::string hname = h->GetName();
                if (hname == "hNWEvts" || hname == "hNEvts") {
                    if (!normHist) { normHist = (TH1*)h->Clone(hname.c_str()); normHist->SetDirectory(nullptr); }
                    delete h;
                    continue;
                }
                auto it = merged.find(hname);
                if (it == merged.end()) {
                    TH1* clone = (TH1*)h->Clone(hname.c_str());
                    clone->SetDirectory(nullptr);
                    merged[hname] = clone;
                    delete h;
                } else {
                    it->second->Add(h);
                    delete h;
                }
            }
            f->Close();
            delete f;
        }

        std::string outPath = targetDir + "/" + base;
        TFile* fout = TFile::Open(outPath.c_str(), "RECREATE");
        fout->cd();
        if (normHist) normHist->Write();
        for (auto& hkv : merged) hkv.second->Write();
        fout->Close();
        delete fout;
        std::cout << "    wrote " << outPath << " (" << merged.size() << " histograms)" << std::endl;

        // These were heap-cloned above and are never attached to a TFile
        // directory (SetDirectory(nullptr)), so nothing but us owns them --
        // without this, one merge job processing an entire year's worth of
        // samples in a single process leaks every sample's histograms and
        // accumulates multiple GB by the time it reaches the later samples.
        delete normHist;
        for (auto& hkv : merged) delete hkv.second;

        for (auto& path : chunkFiles) gSystem->Unlink(path.c_str());
        std::cout << "    deleted " << chunkFiles.size() << " chunk file(s)" << std::endl;
    }
}
