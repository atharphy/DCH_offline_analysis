#pragma once

#include <cstdlib>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "TSystem.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"

inline std::vector<std::string> metSysPayloadCandidates(const std::string& year) {
    const std::string directory = "HTT-utilities/RecoilCorrections/data/";
    if (year == "2016preVFP" || year == "2016postVFP") return {directory + "PFMEtSys_2016.root"};
    if (year == "2017") return {directory + "PFMEtSys_2017.root", directory + "MEtSys_2017.root"};
    if (year == "2018") return {directory + "PFMEtSys_2018.root"};
    return {};
}

inline std::shared_ptr<MEtSys> loadMEtSys(const std::string& year) {
    const std::vector<std::string> candidates = metSysPayloadCandidates(year);
    for (const std::string& payload : candidates) {
        if (!gSystem->AccessPathName(payload.c_str())) return std::make_shared<MEtSys>(payload.c_str());
    }
    return nullptr;
}

inline std::shared_ptr<MEtSys> getCachedMEtSys(const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<MEtSys>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;
    auto sys = loadMEtSys(year);
    cache[year] = sys;
    return sys;
}

struct RecoilSystShift {
    bool valid = false;
    float responseUpPx = 0, responseUpPy = 0;
    float responseDownPx = 0, responseDownPy = 0;
    float resolutionUpPx = 0, resolutionUpPy = 0;
    float resolutionDownPx = 0, resolutionDownPy = 0;
};

inline RecoilSystShift computeRecoilSystShift(MEtSys& metSys, float correctedMetPx, float correctedMetPy, float genPx, float genPy, float visPx, float visPy, int njets) {
    RecoilSystShift result;
    metSys.ApplyMEtSys(correctedMetPx, correctedMetPy, genPx, genPy, visPx, visPy, njets, MEtSys::ProcessType::BOSON, MEtSys::SysType::Response, MEtSys::SysShift::Up, result.responseUpPx, result.responseUpPy);
    metSys.ApplyMEtSys(correctedMetPx, correctedMetPy, genPx, genPy, visPx, visPy, njets, MEtSys::ProcessType::BOSON, MEtSys::SysType::Response, MEtSys::SysShift::Down, result.responseDownPx, result.responseDownPy);
    metSys.ApplyMEtSys(correctedMetPx, correctedMetPy, genPx, genPy, visPx, visPy, njets, MEtSys::ProcessType::BOSON, MEtSys::SysType::Resolution, MEtSys::SysShift::Up, result.resolutionUpPx, result.resolutionUpPy);
    metSys.ApplyMEtSys(correctedMetPx, correctedMetPy, genPx, genPy, visPx, visPy, njets, MEtSys::ProcessType::BOSON, MEtSys::SysType::Resolution, MEtSys::SysShift::Down, result.resolutionDownPx, result.resolutionDownPy);
    result.valid = true;
    return result;
}
