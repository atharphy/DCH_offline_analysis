#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TChain.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <vector>

void MergeChunkHists(std::string chunkDir, std::string targetDir, std::vector<std::string> allowedBases = {}) {
    std::set<std::string> allowedSet(allowedBases.begin(), allowedBases.end());
    gSystem->mkdir(targetDir.c_str(), kTRUE);
    TSystemDirectory dir("chunks", chunkDir.c_str());
    TList* entries = dir.GetListOfFiles();
    if (!entries) { std::cerr << "ERROR: cannot list " << chunkDir << std::endl; return; }

    std::map<std::string, std::vector<std::string>> groups;
    std::regex chunkPattern("^(.*)_chunk[0-9]+\\.root$");
    TIter next(entries);
    TObject* obj;
    while ((obj = next())) {
        std::string name = obj->GetName();
        if (name.size() < 5 || name.substr(name.size() - 5) != ".root") continue;
        std::smatch m;
        if (!std::regex_match(name, m, chunkPattern)) continue;
        std::string base = m[1].str() + ".root";
        if (!allowedSet.empty() && allowedSet.find(base) == allowedSet.end()) continue;
        groups[base].push_back(chunkDir + "/" + name);
    }
    delete entries;

    std::cout << "Found " << groups.size() << " sample(s) to merge in " << chunkDir
               << (allowedSet.empty() ? "" : " (filtered slice)") << std::endl;

    for (auto& kv : groups) {
        const std::string& base = kv.first;
        std::vector<std::string>& chunkFiles = kv.second;
        std::sort(chunkFiles.begin(), chunkFiles.end());
        std::cout << "  merging " << chunkFiles.size() << " chunk(s) -> " << base << std::endl;

        std::map<std::string, TH1*> merged;
        std::set<std::string> treeNames;
        TH1* normHist = nullptr;

        for (const auto& path : chunkFiles) {
            TFile* f = TFile::Open(path.c_str(), "READ");
            if (!f || f->IsZombie()) { std::cerr << "    [warn] cannot open " << path << std::endl; continue; }
            TIter keyIter(f->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)keyIter())) {
                TObject* keyObj = key->ReadObj();
                TH1* h = dynamic_cast<TH1*>(keyObj);
                if (!h) {
                    TTree* t = dynamic_cast<TTree*>(keyObj);
                    if (t) treeNames.insert(t->GetName());
                    delete keyObj;
                    continue;
                }
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
        TFile* fout = TFile::Open(outPath.c_str(), "RECREATE", "", ROOT::CompressionSettings(ROOT::kLZMA, 9));
        if (!fout || fout->IsZombie()) {
            std::cerr << "    [skip] cannot create output file: " << outPath << std::endl;
            if (fout) fout->Close();
            delete fout;
            for (auto& hkv : merged) delete hkv.second;
            delete normHist;
            continue;
        }
        fout->cd();
        if (normHist) normHist->Write();
        for (auto& hkv : merged) hkv.second->Write();
        for (const auto& tname : treeNames) {
            TChain chain(tname.c_str());
            for (const auto& path : chunkFiles) chain.Add(path.c_str());
            fout->cd();
            TTree* mergedTree = chain.CloneTree(-1);
            mergedTree->Write();
        }
        fout->Close();
        delete fout;
        std::cout << "    wrote " << outPath << " (" << merged.size() << " histograms, " << treeNames.size() << " tree(s))" << std::endl;

        delete normHist;
        for (auto& hkv : merged) delete hkv.second;

        for (auto& path : chunkFiles) gSystem->Unlink(path.c_str());
        std::cout << "    deleted " << chunkFiles.size() << " chunk file(s)" << std::endl;
    }
}
