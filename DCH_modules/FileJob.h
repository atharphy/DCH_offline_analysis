#pragma once

#include <string>

#include "Rtypes.h"

struct FileJob {
    std::string year, process, fileName, outDir;
    double lumi = 1.0;
    int index = 0, total = 0;
    // -1 means "process the whole file" (default, existing behavior for
    // every script using FileJob). Set both to a non-negative range to
    // process only entries [startEntry, endEntry) of the file instead.
    Long64_t startEntry = -1, endEntry = -1;
    // Set alongside startEntry/endEntry when chunking: used only to name
    // the per-chunk output file distinctly (e.g. hist_<sample>_chunk3.root).
    int chunkIndex = -1;
};
