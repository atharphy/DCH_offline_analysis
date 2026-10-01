#pragma once

#include <map>
#include <string>
#include <vector>

static std::string year = "Run2";

static std::string SIGNAL_MASS_FILTER = "HppM500";

struct BinSpec {
    int nbins;
    double xmin;
    double xmax;
    std::vector<double> edges;
};

static std::map<std::string, BinSpec> binning = {
    {"mZ1", {40, 70, 110, {}}},
    {"mZ1:VR", {50, 0, 500, {}}},
    {"mZ1:SR", {50, 0, 1000, {}}},
    {"mZ2", {50, 0, 1000, {}}},
    {"mZ2:SR", {50, 0, 1000, {}}},
    {"mH1", {20, 0, 1000, {}}},
    {"mH2", {20, 0, 1000, {}}},
    {"zPt", {40, 0, 120, {}}},
    {"zPt:SR", {50, 0, 500, {}}},
    {"met", {50, 0, 200, {}}},
    {"met:SR", {50, 0, 1000, {}}},
    {"metphi", {16, -3.2, 3.2, {}}},
    {"LT",  {50, 0, 1000, {}}},
    {"LT:SR",  {100, 0, 2000, {}}},

    {"pt1", {40, 0, 200, {}}},
    {"pt2", {40, 0, 200, {}}},
    {"pt3", {40, 0, 200, {}}},
    {"pt4", {40, 0, 200, {}}},
    {"pt1:SR", {50, 0, 1000, {}}},
    {"pt2:SR", {50, 0, 1000, {}}},
    {"pt3:SR", {50, 0, 1000, {}}},
    {"pt4:SR", {50, 0, 1000, {}}},

    {"eta1", {10, -3, 3, {}}},
    {"eta2", {10, -3, 3, {}}},
    {"eta3", {10, -3, 3, {}}},
    {"eta4", {10, -3, 3, {}}},

    {"phi1", {10, -3.5, 3.5, {}}},
    {"phi2", {10, -3.5, 3.5, {}}},
    {"phi3", {10, -3.5, 3.5, {}}},
    {"phi4", {10, -3.5, 3.5, {}}},

    {"d01", {10, -0.045, 0.045, {}}},
    {"d02", {10, -0.045, 0.045, {}}},
    {"d03", {10, -0.045, 0.045, {}}},
    {"d04", {10, -0.045, 0.045, {}}},

    {"dZ1", {10, -0.1, 0.1, {}}},
    {"dZ2", {10, -0.1, 0.1, {}}},
    {"dZ3", {10, -0.1, 0.1, {}}},
    {"dZ4", {10, -0.1, 0.1, {}}},

    {"iso1", {10, 0, 0.6, {}}},
    {"iso2", {10, 0, 0.6, {}}},
    {"iso3", {10, 0, 0.6, {}}},
    {"iso4", {10, 0, 0.6, {}}}
};

static std::vector<std::string> allVariables = {
    "mZ1","mZ2","mH1","mH2","zPt","met","metphi","LT",
    "pt1","pt2","pt3","pt4",
    "eta1","eta2","eta3","eta4",
    "phi1","phi2","phi3","phi4",
    "d01","d02","d03","d04",
    "dZ1","dZ2","dZ3","dZ4",
    "iso1","iso2","iso3","iso4"
};

static std::vector<std::string> finalStates = {
    "ee","em","et","mm","mt","tt",

    "eee","eem","eet",
    "eme","emm","emt",
    "ete","etm","ett",
    "mme","mmm","mmt",
    "mte","mtm","mtt",
    "tte","ttm","ttt",

    "eeee","eeem","eeet","eemm","eemt","eett",
    "emem","emet","emmm","emmt","emtt",
    "etet","etmm","etmt","ettt",
    "mmmm","mmmt","mmtt",
    "mtmt","mttt",
    "tttt"
};

static std::vector<std::string> regions = {
    "DYCR_0tau", "DYCR_1tau", "DYveto_0tau", "DYveto_1tau",

    "CR_0tau", "CR_1tau", "CR_2tau", "CR_3tau",
    "CR_3lep0tau", "CR_3lep1tau", "CR_3lep2tau",

    "VR_0tau", "VR_1tau", "VR_2tau", "VR_3tau",
    "VR_3lep0tau", "VR_3lep1tau", "VR_3lep2tau",

    "SR_0tau", "SR_1tau", "SR_2tau", "SR_3tau",
    "SR_3lep0tau", "SR_3lep1tau", "SR_3lep2tau"
};
