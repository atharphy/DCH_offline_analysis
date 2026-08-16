#pragma once

#include <string>

struct FileJob {
    std::string year, process, fileName, outDir;
    double lumi = 1.0;
    int index = 0, total = 0;
};
