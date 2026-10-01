#pragma once

#include <map>
#include <string>
#include <vector>

const std::vector<std::string> kFinalStates = {
    "ee", "em", "et", "mm", "mt", "tt",
    "eee", "eem", "eet", "eme", "emm", "emt",
    "ete", "etm", "ett", "mme", "mmm", "mmt",
    "mte", "mtm", "mtt", "tte", "ttm", "ttt",
    "eeee", "eeem", "eeet", "eemm", "eemt", "eett",
    "emem", "emet", "emmm", "emmt", "emtt", "etet",
    "etmm", "etmt", "ettt", "mmmm", "mmmt", "mmtt",
    "mtmt", "mttt", "tttt"
};

const std::vector<std::string> kProcesses = {
    "DY", "DY10_50", "WZ", "WW", "ZZ", "TTbar", "ttV",
    "WJ", "ST", "VVV", "QCD", "other", "signal"
};

const std::map<std::string, int> kColors = {
    {"DY", 7}, {"DY10_50", 20}, {"WZ", 8}, {"WW", 40},
    {"ZZ", 5}, {"TTbar", 46}, {"ttV", 4},
    {"WJ", 9}, {"ST", 38}, {"VVV", 6},
    {"QCD", kOrange - 3}, {"other", 3}, {"signal", 2}
};

const std::map<std::string, std::string> kLabels = {
    {"DY", "DY"}, {"DY10_50", "DY10_50"}, {"WZ", "WZ"}, {"WW", "WW"},
    {"ZZ", "ZZ"}, {"TTbar", "TTbar"}, {"ttV", "ttV"},
    {"WJ", "WJ"}, {"ST", "ST"}, {"VVV", "VVV"},
    {"QCD", "QCD"}, {"other", "other"},
    {"signal", "signal"}, {"data", "Data"}
};

const std::vector<std::string> kBinSuffixes = {
    "0tau", "1tau", "2tau", "3tau",
    "3lep0tau", "3lep1tau", "3lep2tau"
};

const std::vector<std::string> kBinLabels = {
    "4l,0#tau", "4l,1#tau", "4l,2#tau", "4l,3#tau",
    "3l,0#tau", "3l,1#tau", "3l,2#tau"
};

inline std::vector<std::string> yearsFor(const std::string& year) {
    if (year == "Run2") return {"2016preVFP", "2016postVFP", "2017", "2018"};
    if (year == "2016") return {"2016preVFP", "2016postVFP"};
    return {year};
}

inline std::string periodLabel(const std::string& year) {
    if (year == "2016preVFP") return "2016 preVFP (19.5 fb^{-1})";
    if (year == "2016postVFP") return "2016 postVFP (16.8 fb^{-1})";
    if (year == "2016") return "2016 (36.3 fb^{-1})";
    if (year == "2017") return "2017 (41.5 fb^{-1})";
    if (year == "2018") return "2018 (59.7 fb^{-1})";
    if (year == "Run2") return "Run 2 (137.6 fb^{-1})";
    return year;
}
