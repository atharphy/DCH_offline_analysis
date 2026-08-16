#pragma once

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

inline void buildPairs(const std::string& catstr, std::vector<std::pair<int,int>>& OS_pair, std::vector<std::pair<int,int>>& SS_pair) {
    std::vector<std::pair<int,int>> Z_pair, Zv_pair, Ztt_pair, OSDF_pair;
    processPairs(catstr.c_str(), Z_pair, Zv_pair, Ztt_pair, SS_pair, OSDF_pair);

    auto pairPtGreater = [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
        return (LepV(a.first)+LepV(a.second)).Pt() > (LepV(b.first)+LepV(b.second)).Pt();
    };
    std::sort(Z_pair.begin(), Z_pair.end(), pairPtGreater);
    std::sort(Zv_pair.begin(), Zv_pair.end(), pairPtGreater);
    std::sort(Ztt_pair.begin(), Ztt_pair.end(), pairPtGreater);
    std::sort(SS_pair.begin(), SS_pair.end(), pairPtGreater);

    Zv_pair = removeOverlap(Zv_pair, static_cast<int>(catstr.size()));

    std::vector<std::pair<int,int>> OSSF_pair;
    OSSF_pair.reserve(Z_pair.size() + Zv_pair.size());
    OSSF_pair.insert(OSSF_pair.end(), Z_pair.begin(), Z_pair.end());
    OSSF_pair.insert(OSSF_pair.end(), Zv_pair.begin(), Zv_pair.end());
    OSSF_pair = removeOverlap(OSSF_pair, static_cast<int>(catstr.size()));

    OS_pair.clear();
    OS_pair.reserve(OSSF_pair.size() + Ztt_pair.size() + OSDF_pair.size());
    OS_pair.insert(OS_pair.end(), OSSF_pair.begin(), OSSF_pair.end());
    OS_pair.insert(OS_pair.end(), Ztt_pair.begin(), Ztt_pair.end());
    OS_pair.insert(OS_pair.end(), OSDF_pair.begin(), OSDF_pair.end());
    OS_pair = removeOverlap(OS_pair, static_cast<int>(catstr.size()));
}
