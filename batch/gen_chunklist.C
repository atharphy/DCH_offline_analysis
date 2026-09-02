// Regenerates chunklist.txt from FileMap.h + each file's real entry count:
// one line per (year, file index, chunk index), chunks of CHUNK_SIZE (1M)
// events each. Needs a valid grid proxy (reads every file's Events tree
// over xrootd just for GetEntries(), no event loop).
//
// Run from Updated_offline_framework/:
//   root -l -b -q batch/gen_chunklist.C

#include "../filemap/FileMap.h"
#include <fstream>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

void gen_chunklist() {
    const Long64_t CHUNK_SIZE = 500000LL;
    std::ofstream out("batch/chunklist.txt");
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
    std::cout << "Wrote batch/chunklist.txt" << std::endl;
}
