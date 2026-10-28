#pragma once

#include <string>

#include "Rtypes.h"

struct FileJob {
    std::string year, process, fileName, outDir;

    std::string treeOutDir;
    double lumi = 1.0;
    int index = 0, total = 0;

    Long64_t startEntry = -1, endEntry = -1;

    int chunkIndex = -1;
};
