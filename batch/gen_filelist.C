// Regenerates filelist.txt from FileMap.h directly, so the condor per-file
// submissions ($(year),$(idx) pairs) always match the sample lists in
// FileMap.h exactly -- no manually-maintained/stale snapshot.
//
// Run from Updated_offline_framework/:
//   root -l -b -q batch/gen_filelist.C

#include "../filemap/FileMap.h"
#include <fstream>
#include <iostream>

void gen_filelist() {
    std::ofstream out("batch/filelist.txt");
    for (std::string year : {"2016preVFP", "2016postVFP", "2017", "2018"}) {
        auto fileMap = getFileMap(year);
        int idx = 0;
        for (auto& entry : fileMap) {
            for (auto& fileName : entry.second) {
                if (fileName.empty() || fileName.back() == '/') continue;
                out << year << "," << idx << "\n";
                ++idx;
            }
        }
        std::cout << year << ": " << idx << " files" << std::endl;
    }
    out.close();
    std::cout << "Wrote batch/filelist.txt" << std::endl;
}
