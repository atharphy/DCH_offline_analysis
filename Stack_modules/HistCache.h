#pragma once

#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"

#include "StackConfig.h"
#include "Labels.h"

inline std::string makeOldRegionName(const std::string& hname) {
    for (const auto& reg : regions) {
        std::string tag = "_" + reg;
        size_t pos = hname.find(tag);
        if (pos != std::string::npos) return hname.substr(0, pos) + "__" + reg;
    }
    return hname;
}

inline std::string normalizeRegionName(const std::string& hname) {
    for (const auto& reg : regions) {
        const std::string oldTag = "__" + reg;
        const size_t pos = hname.find(oldTag);
        if (pos != std::string::npos) return hname.substr(0, pos) + "_" + reg + hname.substr(pos + oldTag.size());
    }
    return hname;
}

static std::unordered_map<TFile*, std::unordered_set<std::string>> fileKeys;

inline const std::unordered_set<std::string>& keysInFile(TFile* f) {
    auto found = fileKeys.find(f);
    if (found != fileKeys.end()) return found->second;

    std::unordered_set<std::string> keys;
    if (f && !f->IsZombie()) {
        TIter next(f->GetListOfKeys());
        while (TKey* key = dynamic_cast<TKey*>(next())) keys.insert(key->GetName());
    }
    return fileKeys.emplace(f, std::move(keys)).first->second;
}

inline TH1D* fetchAndClone(TFile* f, const std::string& hname) {
    if (!f) return nullptr;

    const auto& keys = keysInFile(f);
    std::string storedName = hname;

    if (!keys.count(storedName)) {
        std::string oldName = makeOldRegionName(hname);
        if (oldName == hname || !keys.count(oldName)) return nullptr;
        storedName = oldName;
    }

    TH1D* h = dynamic_cast<TH1D*>(f->Get(storedName.c_str()));
    if (!h) return nullptr;

    TH1D* c = (TH1D*)h->Clone();
    c->SetDirectory(nullptr);
    c->Sumw2();

    std::string var = getVariableFromHname(hname);
    const std::string regionKey = var + ":" + getRegionGroupFromHname(hname);
    std::string key;
    if (binning.count(regionKey)) key = regionKey;
    else if (binning.count(var)) key = var;
    else return c;

    BinSpec b = binning[key];
    TH1D* reb = nullptr;

    if (!b.edges.empty()) {
        reb = new TH1D((std::string(c->GetName()) + "_reb").c_str(), c->GetTitle(), b.edges.size() - 1, b.edges.data());
    } else {
        reb = new TH1D((std::string(c->GetName()) + "_reb").c_str(), c->GetTitle(), b.nbins, b.xmin, b.xmax);
    }
    reb->Sumw2();

    for (int ib = 1; ib <= c->GetNbinsX() + 1; ++ib) {
        const double x = c->GetXaxis()->GetBinCenter(ib);
        const double w = c->GetBinContent(ib);
        const double e = c->GetBinError(ib);

        int newbin = reb->FindBin(x);
        if (newbin < 1) continue;
        if (newbin > reb->GetNbinsX()) newbin = reb->GetNbinsX();

        reb->SetBinContent(newbin, reb->GetBinContent(newbin) + w);
        reb->SetBinError(newbin, std::hypot(reb->GetBinError(newbin), e));
    }

    delete c;
    return reb;
}
