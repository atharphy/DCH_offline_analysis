#include "../../filemap/FileMap.h"
#include <fstream>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

void gen_chunklist() {
    const Long64_t CHUNK_SIZE = 500000LL;
    std::ofstream out("batch/dch_tight/chunklist.txt");
    for (std::string year : {"2016preVFP", "2016postVFP", "2017", "2018"}) {
        auto fileMap = getFileMap(year);
        int idx = 0;
        Long64_t totalChunks = 0;
        for (auto& entry : fileMap) {
            for (auto& fileName : entry.second) {
                if (fileName.empty() || fileName.back() == '/') continue;
                TFile* f = TFile::Open(fileName.c_str(), "READ");
                Long64_t nEnt = 0;
                if (f && !f->IsZombie()) {
                    TTree* t = dynamic_cast<TTree*>(f->Get("Events"));
                    if (t) nEnt = t->GetEntries();
                    f->Close();
                }
                delete f;
                const int nChunks = std::max(1, static_cast<int>((nEnt + CHUNK_SIZE - 1) / CHUNK_SIZE));
                for (int c = 0; c < nChunks; ++c) out << year << "," << idx << "," << c << "\n";
                std::cout << year << " idx=" << idx << " " << entry.first << " entries=" << nEnt << " chunks=" << nChunks << std::endl;
                totalChunks += nChunks;
                ++idx;
            }
        }
        std::cout << year << ": " << idx << " files, " << totalChunks << " total chunks" << std::endl;
    }
    out.close();
    std::cout << "Wrote batch/dch_tight/chunklist.txt" << std::endl;
}
