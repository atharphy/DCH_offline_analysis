#include "../filemap/FileMapHTQCD.h"
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
void gen_treeonly_qcdht_list() {
    const Long64_t CHUNK_SIZE = 500000LL; const int BATCH = 15;
    std::ofstream out("batch/batched_chunklist_treeonly_qcdht.txt"), ml("batch/mergelist_treeonly_qcdht.txt");
    for (std::string year : {"2016preVFP", "2016postVFP", "2017", "2018"}) {
        auto fm = getFileMapHTQCD(year); int idx = 0;
        for (auto& e : fm) for (auto& fn : e.second) {
            if (fn.empty() || fn.back()=='/') continue;
            TFile* f = TFile::Open(fn.c_str(), "READ"); Long64_t n = 0;
            if (f && !f->IsZombie()) { TTree* t = dynamic_cast<TTree*>(f->Get("Events")); if (t) n = t->GetEntries(); f->Close(); }
            delete f;
            int nc = std::max(1, (int)((n + CHUNK_SIZE - 1)/CHUNK_SIZE));
            for (int s = 0; s < nc; s += BATCH) out << year << "," << idx << "," << s << "," << std::min(BATCH, nc - s) << "\n";
            std::string b = gSystem->BaseName(fn.c_str()); b = b.substr(0, b.rfind(".root"));
            ml << year << "," << b << "\n";
            std::cout << year << " " << idx << " " << e.first << " " << b << " n=" << n << " chunks=" << nc << std::endl;
            ++idx;
        }
    }
}
