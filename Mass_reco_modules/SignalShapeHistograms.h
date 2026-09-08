#pragma once

// Tree-only bundle for the 4-lepton/3-lepton signal-shape mass-reconstruction
// output (0tau/1tau/2tau/3tau/3lep0tau/3lep1tau/3lep2tau regions, plus a
// systematic-suffixed tree per region for each weight/shape variation --
// e.g. "Events_0tau_roccorUp"). No histograms: this module used to also
// fill a full ~30-variable x per-systematic histogram set (mirroring the
// main HistBundle), but that duplicated what the main hist_<sample>.root
// output already has and was most of why these files ballooned to
// thousands of keys. Combine only ever needs the discriminant itself,
// which lives in the tree branches -- bin whichever branch is needed into
// a histogram per systematic at datacard-build time instead.

#include <memory>
#include <string>
#include <unordered_map>

#include "TTree.h"

// mH1/mH2 = get_legs_comb's masses (Mdch01/Mdch02); mDCH1/mDCH2 =
// run_reco_event's "fit"-policy masses (MdchFit1/MdchFit2) -- kept as
// separate branches so either reconstruction method's mass is available
// downstream without re-deriving it.
struct SignalShapeTreeVars {
    double mll1 = 0, mll2 = 0, mH1 = 0, mH2 = 0, mDCH1 = 0, mDCH2 = 0, evtwt = 0;
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
