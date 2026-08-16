#pragma once

#include <string>
#include <unordered_map>

#include "TH1D.h"

#include "CommonConfig.h"

struct HistSpec {
    std::string name, title;
    int bins = 0;
    double xmin = 0.0, xmax = 0.0;
};

inline TH1D* getHist(std::unordered_map<std::string,TH1D*>& histMap, const std::string& key, const HistSpec& spec) {
    auto found = histMap.find(key);
    if (found != histMap.end()) return found->second;
    const std::string histName = "h_" + spec.name + "_" + key;
    const std::string histTitle = spec.title + " " + key;
    TH1D* hist = new TH1D(histName.c_str(), histTitle.c_str(), spec.bins, spec.xmin, spec.xmax);
    hist->Sumw2();
    hist->SetDirectory(nullptr);
    histMap[key] = hist;
    return hist;
}

inline void fillOne(std::unordered_map<std::string,TH1D*>& histMap, const std::string& key, double value, double weight, const HistSpec& spec) {
    getHist(histMap, key, spec)->Fill(value, weight);
}

inline void writeHists(std::unordered_map<std::string,TH1D*>& histMap) {
    for (auto& item : histMap) {
        TH1D* hist = item.second;
        if (!hist) continue;
        if (hist->GetEntries()==0 && hist->Integral()==0.0) continue;
        hist->Write();
    }
}

const HistSpec S_mZ1 {"mZ1","mZ1",1000,0.0,1000.0};
const HistSpec S_mZ2 {"mZ2","mZ2",1000,0.0,1000.0};
const HistSpec S_mH1 {"mH1","mH1",1000,0.0,1000.0};
const HistSpec S_mH2 {"mH2","mH2",1000,0.0,1000.0};
const HistSpec S_zPt {"zPt","Z p_{T}",500,0.0,500.0};
const HistSpec S_met {"met","MET",1000,0.0,1000.0};
const HistSpec S_metphi {"metphi", "MET #phi", 64, -3.2, 3.2};
const HistSpec S_LT {"LT","LT",1000,0.0,2000.0};

const HistSpec S_pt[NlepMax] = {{"pt1","pt1",1000,0.0,1000.0},{"pt2","pt2",1000,0.0,1000.0},{"pt3","pt3",1000,0.0,1000.0},{"pt4","pt4",1000,0.0,1000.0}};
const HistSpec S_eta[NlepMax] = {{"eta1","eta1",60,-3.0,3.0},{"eta2","eta2",60,-3.0,3.0},{"eta3","eta3",60,-3.0,3.0},{"eta4","eta4",60,-3.0,3.0}};
const HistSpec S_phi[NlepMax] = {{"phi1","phi1",64,-3.2,3.2},{"phi2","phi2",64,-3.2,3.2},{"phi3","phi3",64,-3.2,3.2},{"phi4","phi4",64,-3.2,3.2}};
const HistSpec S_d0[NlepMax] = {{"d01","d01",60,-0.06,0.06},{"d02","d02",60,-0.06,0.06},{"d03","d03",60,-0.06,0.06},{"d04","d04",60,-0.06,0.06}};
const HistSpec S_dZ[NlepMax] = {{"dZ1","dZ1",60,-0.15,0.15},{"dZ2","dZ2",60,-0.15,0.15},{"dZ3","dZ3",60,-0.15,0.15},{"dZ4","dZ4",60,-0.15,0.15}};
const HistSpec S_iso[NlepMax] = {{"iso1","iso1",60,0.0,0.6},{"iso2","iso2",60,0.0,0.6},{"iso3","iso3",60,0.0,0.6},{"iso4","iso4",60,0.0,0.6}};
