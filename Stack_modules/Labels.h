#pragma once

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "TLatex.h"

#include "StackConfig.h"

inline std::string getVariableFromHname(const std::string& hname) {
    std::vector<std::string> vars = allVariables;
    std::sort(vars.begin(), vars.end(), [](const std::string& a, const std::string& b) { return a.size() > b.size(); });
    for (const auto& var : vars) {
        std::string pattern = "h_" + var + "_";
        if (hname.find(pattern) == 0) return var;
    }
    return hname;
}

inline std::string getRegionGroupFromHname(const std::string& hname) {
    std::vector<std::string> regs = regions;
    std::sort(regs.begin(), regs.end(), [](const std::string& a, const std::string& b) { return a.size() > b.size(); });
    for (const auto& reg : regs) {
        const std::string suffix = "_" + reg;
        if (hname.size() >= suffix.size() && hname.compare(hname.size() - suffix.size(), suffix.size(), suffix) == 0) {
            const size_t us = reg.find('_');
            return us == std::string::npos ? reg : reg.substr(0, us);
        }
    }
    return "";
}

inline std::string getXaxisLabel(const std::string& hname) {
    std::string var = getVariableFromHname(hname);

    static std::map<std::string, std::string> labels = {
        {"mZ1", "m_{Z,1} [GeV]"},
        {"mZ2", "m_{Z,2} [GeV]"},
        {"mH1", "m_{H,1} [GeV]"},
        {"mH2", "m_{H,2} [GeV]"},
        {"zPt", "Z p_{T} [GeV]"},
        {"met", "E_{T}^{miss} [GeV]"},
        {"metphi", "#phi(E_{T}^{miss})"},
        {"LT",  "L_{T} [GeV]"},

        {"pt1", "p_{T}^{1} [GeV]"},
        {"pt2", "p_{T}^{2} [GeV]"},
        {"pt3", "p_{T}^{3} [GeV]"},
        {"pt4", "p_{T}^{4} [GeV]"},

        {"eta1", "#eta_{1}"},
        {"eta2", "#eta_{2}"},
        {"eta3", "#eta_{3}"},
        {"eta4", "#eta_{4}"},

        {"phi1", "#phi_{1}"},
        {"phi2", "#phi_{2}"},
        {"phi3", "#phi_{3}"},
        {"phi4", "#phi_{4}"},

        {"d01", "d_{0}^{1} [cm]"},
        {"d02", "d_{0}^{2} [cm]"},
        {"d03", "d_{0}^{3} [cm]"},
        {"d04", "d_{0}^{4} [cm]"},

        {"dZ1", "d_{z}^{1} [cm]"},
        {"dZ2", "d_{z}^{2} [cm]"},
        {"dZ3", "d_{z}^{3} [cm]"},
        {"dZ4", "d_{z}^{4} [cm]"},

        {"iso1", "I_{rel}^{1}"},
        {"iso2", "I_{rel}^{2}"},
        {"iso3", "I_{rel}^{3}"},
        {"iso4", "I_{rel}^{4}"}
    };

    if (labels.count(var)) return labels[var];
    return var;
}

inline std::string getLumiLabel(const std::string& y) {
    if (y == "2016preVFP")  return "2016 preVFP (19.5 fb^{-1})";
    if (y == "2016postVFP") return "2016 postVFP (16.8 fb^{-1})";
    if (y == "2016")        return "2016 (36.3 fb^{-1})";
    if (y == "2017")        return "2017 (41.5 fb^{-1})";
    if (y == "2018")        return "2018 (59.7 fb^{-1})";
    if (y == "Run2")        return "Run 2 (137.6 fb^{-1})";
    return y;
}

inline void drawCMSAndLumi() {
    TLatex cms;
    cms.SetNDC();
    cms.SetTextFont(61);
    cms.SetTextSize(0.068);
    cms.DrawLatex(0.17, 0.905, "CMS");

    TLatex prelim;
    prelim.SetNDC();
    prelim.SetTextFont(52);
    prelim.SetTextSize(0.064);
    prelim.DrawLatex(0.285, 0.905, "Preliminary");

    TLatex lumi;
    lumi.SetNDC();
    lumi.SetTextFont(42);
    lumi.SetTextSize(0.060);
    lumi.SetTextAlign(31);
    lumi.DrawLatex(0.96, 0.905, getLumiLabel(year).c_str());
}
