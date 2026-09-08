// Regenerates mergelist_persample.txt from FileMap.h: one "year,sample"
// line per (year, file) pair -- "sample" here is the file's own basename
// stem, matching the hist_<stem>/masstree_<stem> naming ProcessTightFile
// actually writes (not the FileMap process key, which can map to several
// files, e.g. "TTbar" -> TTTo2L2Nu/TTToSemiLeptonic/TTToHadronic).
//
// Run from Updated_offline_framework/:
//   root -l -b -q batch/gen_mergelist_persample.C

#include "../filemap/FileMap.h"
#include <fstream>
#include <iostream>

void gen_mergelist_persample() {
    std::ofstream out("batch/mergelist_persample.txt");
    for (const std::string year : {"2016preVFP", "2016postVFP", "2017", "2018"}) {
        auto fileMap = getFileMap(year);
        for (auto& entry : fileMap) {
            for (auto& fileName : entry.second) {
                if (fileName.empty() || fileName.back() == '/') continue;
                std::string base = gSystem->BaseName(fileName.c_str());
                const size_t dotRoot = base.rfind(".root");
                if (dotRoot != std::string::npos) base = base.substr(0, dotRoot);
                out << year << "," << base << "\n";
            }
        }
    }
    std::cout << "Wrote batch/mergelist_persample.txt" << std::endl;
}
