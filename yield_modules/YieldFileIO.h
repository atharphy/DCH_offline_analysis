#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TSystem.h"

#include "YieldConfig.h"

struct YieldOpenResult {
    int nOpened = 0;
    int nMissing = 0;
};

inline std::string ensureTrailingSlash(const std::string& dir) {
    if (dir.empty() || dir.back() == '/') return dir;
    return dir + "/";
}

inline YieldOpenResult openYieldInputFiles(const std::string& year, const std::string& histogramFolder, std::map<std::string, std::vector<TFile*>>& filesByProcess) {
    YieldOpenResult result;
    const std::string base = ensureTrailingSlash(histogramFolder);

    for (const std::string& inputYear : yearsFor(year)) {
        const auto fileMap = getFileMap(inputYear);
        const std::string inputDirectory = base + inputYear;

        for (const auto& process : fileMap) {
            for (const std::string& sourceName : process.second) {
                const std::string sampleName = std::string(gSystem->BaseName(sourceName.c_str()));
                const std::string fileName = inputDirectory + "/hist_" + sampleName;

                if (gSystem->AccessPathName(fileName.c_str())) {
                    ++result.nMissing;
                    continue;
                }

                TFile* file = TFile::Open(fileName.c_str(), "READ");
                if (!file || file->IsZombie()) {
                    if (file) delete file;
                    ++result.nMissing;
                    continue;
                }

                filesByProcess[process.first].push_back(file);
                ++result.nOpened;
            }
        }
    }

    return result;
}
