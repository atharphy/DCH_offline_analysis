#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "TTree.h"

struct SignalShapeTreeVars {
    double mll1 = 0, mll2 = 0, mH1 = 0, mH2 = 0, mDCH1 = 0, mDCH2 = 0, evtwt = 0, LT = 0, met = 0, metphi = 0;
    int cat = 0, gencat = 0;
    std::string region;
};

struct SignalShapeHistBundle {
    std::unordered_map<std::string, TTree*> tree;
    std::unordered_map<std::string, std::unique_ptr<SignalShapeTreeVars>> treeVars;

    void reserveAll(size_t n) {
        tree.reserve(n);
        treeVars.reserve(n);
    }

    TTree* getTree(const std::string& region) {
        auto it = tree.find(region);
        if (it != tree.end()) return it->second;
        TTree* t = new TTree(("Events_" + region).c_str(), "");
        t->SetDirectory(nullptr);
        auto vars = std::make_unique<SignalShapeTreeVars>();
        t->Branch("mll1", &vars->mll1);
        t->Branch("mll2", &vars->mll2);
        t->Branch("mH1", &vars->mH1);
        t->Branch("mH2", &vars->mH2);
        t->Branch("mDCH1", &vars->mDCH1);
        t->Branch("mDCH2", &vars->mDCH2);
        t->Branch("evtwt", &vars->evtwt);
        t->Branch("LT", &vars->LT);
        t->Branch("met", &vars->met);
        t->Branch("metphi", &vars->metphi);
        t->Branch("cat", &vars->cat);
        t->Branch("gencat", &vars->gencat);
        t->Branch("region", &vars->region);
        tree[region] = t;
        treeVars[region] = std::move(vars);
        return t;
    }

    void writeAll() {
        for (auto& kv : tree) if (kv.second->GetEntries() > 0) kv.second->Write();
    }

    void deleteAll() {
        for (auto& kv : tree) delete kv.second;
        tree.clear();
        treeVars.clear();
    }
};
