#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TSystem.h"

struct OpenFilesResult {
    int nOpened = 0;
    int nMissing = 0;
    int nBad = 0;
};

inline std::string ensureTrailingSlash(const std::string& dir) {
    if (dir.empty() || dir.back() == '/') return dir;
    return dir + "/";
}

inline std::vector<std::string> resolveYearsToRead(const std::string& inYear) {
    if (inYear == "Run2") return {"2016preVFP", "2016postVFP", "2017", "2018"};
    if (inYear == "2016") return {"2016preVFP", "2016postVFP"};
    if (inYear == "2016preVFP" || inYear == "2016postVFP" || inYear == "2017" || inYear == "2018") return {inYear};
    std::cerr << "ERROR: unsupported year argument: " << inYear << std::endl;
    return {};
}

inline OpenFilesResult openInputFiles(const std::string& inYear, const std::string& histBase, const std::string& filePrefix, std::map<std::string, std::vector<TFile*>>& handles) {
    OpenFilesResult result;
    const std::vector<std::string> yearsToRead = resolveYearsToRead(inYear);

    for (const std::string& inputYear : yearsToRead) {
        const std::map<std::string, std::vector<std::string>> files = getFileMap(inputYear);
        if (files.empty()) {
            std::cerr << "[warning] FileMap is empty for " << inputYear << std::endl;
            continue;
        }

        const std::string inputDir = histBase + inputYear + "/";
        std::cout << "\nReading year : " << inputYear << std::endl;
        std::cout << "Input path   : " << inputDir << std::endl;

        for (const auto& processEntry : files) {
            const std::string& process = processEntry.first;

            for (const std::string& originalFile : processEntry.second) {
                const std::string baseName = gSystem->BaseName(originalFile.c_str());

                if (process == "signal" && !SIGNAL_MASS_FILTER.empty() && baseName.find(SIGNAL_MASS_FILTER) == std::string::npos) continue;

                const std::string fullName = inputDir + filePrefix + baseName;

                if (gSystem->AccessPathName(fullName.c_str())) {
                    std::cout << "[missing file] " << fullName << std::endl;
                    ++result.nMissing;
                    continue;
                }

                TFile* input = TFile::Open(fullName.c_str(), "READ");
                if (!input || input->IsZombie()) {
                    std::cerr << "[bad file] " << fullName << std::endl;
                    if (input) { input->Close(); delete input; }
                    ++result.nBad;
                    continue;
                }

                handles[process].push_back(input);
                ++result.nOpened;
            }
        }
    }

    return result;
}
