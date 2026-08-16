#include "../filemap/FileMap.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void MakeFileList() {
    const std::vector<std::string> years = {"2016preVFP", "2016postVFP", "2017", "2018"};
    std::ofstream out("filelist.txt");

    int total = 0;
    for (const std::string& year : years) {
        auto fileMap = getFileMap(year);
        int n = 0;
        for (const auto& item : fileMap) n += item.second.size();
        for (int idx = 0; idx < n; ++idx) out << year << "," << idx << "\n";
        std::cout << year << ": " << n << " files" << std::endl;
        total += n;
    }
    out.close();
    std::cout << "wrote filelist.txt with " << total << " entries" << std::endl;
}
