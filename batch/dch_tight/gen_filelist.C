#include "../../filemap/FileMap.h"
#include <fstream>
#include <iostream>

void gen_filelist() {
    std::ofstream out("batch/dch_tight/filelist.txt");
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
    std::cout << "Wrote batch/dch_tight/filelist.txt" << std::endl;
}
