#pragma once

// Relies on include/cat.h, include/Kinematics.C (numberToCat, processPairs,
// removeOverlap, LepV, deltaR, deltaPhi, calculateMT, calculateMTtot) and
// include/MyBranch.C's global branch variables already being in scope --
// same assumption PairBuilder.h/HistFilling.h already make.

#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "TLorentzVector.h"

inline int catToNumber(const std::string& catStr) {
    for (int i = 0; i < 49; ++i) {
        if (catStr == numberToCat(i)) return i;
    }
    return -1;
}

// Builds the OSSF pair list matching guru_offline/DCH_presel_tree.C's own
// convention -- its local processPairs mapped "Zwindow", "Zv", AND
// "ZttPair" all into the same OSSF_pair bucket (it called its own
// processPairs with the same output vector passed as both the Z_pair and
// Zv_pair args, and its "ZttPair" case pushed into Zv_pair too). That's
// different from PairBuilder.h::buildPairs' own internal OSSF_pair (which
// deliberately keeps Ztt_pair out, folding it into OS_pair separately) --
// using that narrower convention here silently emptied OSSF_pair for any
// category whose only same-flavor OS pair is tau-tau (e.g. ettt/mttt have
// no e-e or m-m pair to fall back on), dropping those events entirely.
inline void buildSignalShapeOSSFPairs(const std::string& catstr, std::vector<std::pair<int,int>>& OSSF_pair) {
    std::vector<std::pair<int,int>> Z_pair, Zv_pair, Ztt_pair, SS_pair_unused, OSDF_pair_unused;
    processPairs(catstr.c_str(), Z_pair, Zv_pair, Ztt_pair, SS_pair_unused, OSDF_pair_unused);

    auto pairPtGreater = [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
        return (LepV(a.first) + LepV(a.second)).Pt() > (LepV(b.first) + LepV(b.second)).Pt();
    };
    std::sort(Z_pair.begin(), Z_pair.end(), pairPtGreater);
    std::sort(Zv_pair.begin(), Zv_pair.end(), pairPtGreater);
    std::sort(Ztt_pair.begin(), Ztt_pair.end(), pairPtGreater);
    Zv_pair = removeOverlap(Zv_pair, static_cast<int>(catstr.size()));

    OSSF_pair.clear();
    OSSF_pair.reserve(Z_pair.size() + Zv_pair.size() + Ztt_pair.size());
    OSSF_pair.insert(OSSF_pair.end(), Z_pair.begin(), Z_pair.end());
    OSSF_pair.insert(OSSF_pair.end(), Zv_pair.begin(), Zv_pair.end());
    OSSF_pair.insert(OSSF_pair.end(), Ztt_pair.begin(), Ztt_pair.end());
    OSSF_pair = removeOverlap(OSSF_pair, static_cast<int>(catstr.size()));
}

// Same logic as Kinematics.C's remaining_idx, but returns -1 instead of
// falling off the end (undefined behavior) when every lepton already
// appears in `pairs` -- can happen once restricted to 4-lepton events,
// where the source script's own remaining_idx has no such guard.
inline int safeRemainingIdx(const std::vector<std::pair<int,int>>& pairs, int nLep) {
    std::set<int> used;
    for (const auto& p : pairs) { used.insert(p.first); used.insert(p.second); }
    for (int i = 1; i <= nLep; ++i) if (!used.count(i)) return i;
    return -1;
}

constexpr double SIGNAL_SHAPE_SENTINEL = -1.0;

// Pair accessor with the same out-of-bounds guard applied uniformly to every
// OSSF_pair[k] access the source script makes unconditionally (it only ever
// checks OSSF_pair.size() >= 1, then reads indices 0..3) -- returns the
// sentinel instead of undefined behavior when fewer than idx+1 OSSF pairs
// exist, which is routine for non-quad-flavor categories (e.g. eeem-derived
// realcats never reach here since we're 4-lepton-only, but eemm/emem-style
// splits can still yield fewer than 4 OSSF pairs).
inline std::pair<int,int> ossfAt(const std::vector<std::pair<int,int>>& OSSF_pair, size_t idx) {
    if (idx < OSSF_pair.size()) return OSSF_pair[idx];
    return {-1, -1};
}

#include <map>
#include "SignalShapeReco.h"
#include "SignalShapeHistograms.h"

// Holds every derived quantity for one event under one kinematic hypothesis
// (nominal, or lepton momenta already shifted for a roccor/tauES variation)
// -- computing this once and filling multiple weight-only variants from it
// avoids re-running the expensive TMinuit-based reconstruction for
// systematics that don't touch lepton momenta (zpt, pu, l1Prefire, xsec,
// topPt, lepton/tau ID-SF). Shape systematics that do shift momenta
// (roccor, tauES) need their own computeSignalShapeVariables call after the
// shift, since LepV() picks up whatever pt_1..4 currently hold.
struct SignalShapeVariables {
    bool eligible = false;
    std::string region;
    double mll1 = 0, mll2 = 0;
    int catNum = 0;
    SignalShapeMasses masses;
};

// Reconstructs the ~30 kinematic variables for one event under whatever
// lepton momenta (pt_1..4/eta.../phi...) currently hold, mirroring
// guru_offline/DCH_presel_tree.C's "Pre" selection block. Handles both
// 4-lepton (cat 1-21) and 3-lepton (cat 22-39) final states.
//
// The source script's own 3-lepton branch was unreachable dead code -- an
// earlier, unconditional "if (cat > 21) continue;" skipped every 3-lepton
// event before the mass reconstruction or hist_variable_map were ever
// built, so its "3lep0tau/3lep1tau/3lep2tau" trees/histograms were always
// empty. There is no working original behavior to replicate for the
// 3-lepton case; this reuses the same 4-slot reconstruction machinery the
// 4-lepton path uses, with the (non-existent) second lepton of the second
// "leg" replaced by a zero four-vector via the existing LepV(-1) convention
// -- L[2] is the one leftover lepton (not in the same-sign pair), L[3] is
// the phantom zero-momentum slot. get_legs_comb/get_legs_tau read past the
// end of a 3-character realcat as the string's null terminator, which
// safely reads as "not a tau" and excludes the phantom slot from neutrino
// consideration.
//
// SS_pair (and its pairing membership) is passed in already computed at
// nominal kinematics, matching the main HistBundle's own convention of
// freezing pairing membership and only letting LepV() pick up shifted
// momenta for systematic variations -- not recomputing which leptons pair
// together under each shift.
inline SignalShapeVariables computeSignalShapeVariables(const std::string& catstr, int catNumber,
                                                          const std::vector<std::pair<int,int>>& SS_pair) {
    SignalShapeVariables out;
    const bool is4Lep = catstr.size() == 4;
    const bool is3Lep = catstr.size() == 3;
    if (!is4Lep && !is3Lep) return out;
    if (is4Lep && (catNumber < 1 || catNumber > 21)) return out;
    if (is3Lep && (catNumber < 22 || catNumber > 39)) return out;
    if (SS_pair.empty()) return out;
    if (is4Lep && SS_pair.size() < 2) return out;

    std::vector<std::pair<int,int>> OSSF_pair;
    buildSignalShapeOSSFPairs(catstr, OSSF_pair);
    if (OSSF_pair.empty()) return out;

    int H1_idx1 = SS_pair[0].first, H1_idx2 = SS_pair[0].second;
    int H2_idx1, H2_idx2;
    if (is4Lep) {
        H2_idx1 = SS_pair[1].first;
        H2_idx2 = SS_pair[1].second;
        if ((LepV(H1_idx1) + LepV(H1_idx2)).M() <= (LepV(H2_idx1) + LepV(H2_idx2)).M()) {
            std::swap(H1_idx1, H2_idx1);
            std::swap(H1_idx2, H2_idx2);
        }
    } else {
        H2_idx1 = safeRemainingIdx({SS_pair[0]}, static_cast<int>(catstr.size()));
        H2_idx2 = -1;  // no second lepton for this leg -- LepV(-1) gives a zero four-vector
    }

    const double mll1 = (LepV(H1_idx1) + LepV(H1_idx2)).M();
    const double mll2 = (LepV(H2_idx1) + LepV(H2_idx2)).M();

    // Only mZ1-4 are still needed here (for the mZ-window veto below) --
    // dR/mT/mTtot/W-leg variables the source script also derived at this
    // point only ever fed the now-removed histogram set, not any cut or
    // tree branch, so they're gone rather than computed and discarded.
    auto ossfMass = [&](size_t idx) {
        auto p = ossfAt(OSSF_pair, idx);
        return p.first < 0 ? SIGNAL_SHAPE_SENTINEL : (LepV(p.first) + LepV(p.second)).M();
    };
    const double mZ1 = ossfMass(0), mZ2 = ossfMass(1), mZ3 = ossfMass(2), mZ4 = ossfMass(3);

    TLorentzVector MET; MET.SetPtEtaPhiM(met, 0, metphi, 0);

    const double st = is4Lep ? (pt_1 + pt_2 + pt_3 + pt_4) : (pt_1 + pt_2 + pt_3);

    std::string realcat;
    realcat.push_back(catstr[H1_idx1 - 1]);
    realcat.push_back(catstr[H1_idx2 - 1]);
    realcat.push_back(catstr[H2_idx1 - 1]);
    if (is4Lep) realcat.push_back(catstr[H2_idx2 - 1]);

    std::array<TLorentzVector,4> L = {LepV(H1_idx1), LepV(H1_idx2), LepV(H2_idx1), LepV(H2_idx2)};
    mat2 metcov;
    metcov[0][0] = metcov00; metcov[0][1] = metcov01;
    metcov[1][0] = metcov10; metcov[1][1] = metcov11;
    const SignalShapeMasses masses = computeSignalShapeMasses(realcat, L, MET, metcov);

    const int Nlep = cat_lepCount(catstr, 'e', 'm');
    const int Ntau = static_cast<int>(catstr.size()) - Nlep;
    const int maxNtau = is4Lep ? 3 : 2;  // Ntau==maxNtau+1 (tttt/ttt) has no destination region, same as the source script
    if (Ntau < 0 || Ntau > maxNtau) return out;

    const bool applyMzVeto = Ntau < maxNtau;
    if (applyMzVeto && (std::abs(mZ1 - 91.2) < 10 || std::abs(mZ2 - 91.2) < 10 ||
                         (is4Lep && (std::abs(mZ3 - 91.2) < 10 || std::abs(mZ4 - 91.2) < 10)))) return out;

    std::string region;
    double stThreshold;
    if (is4Lep) {
        region = (Ntau == 0) ? "0tau" : (Ntau == 1) ? "1tau" : (Ntau == 2) ? "2tau" : "3tau";
        stThreshold = (Ntau == 3) ? 100.0 : 400.0;
    } else {
        region = (Ntau == 0) ? "3lep0tau" : (Ntau == 1) ? "3lep1tau" : "3lep2tau";
        stThreshold = (Ntau == 0) ? 300.0 : (Ntau == 1) ? 200.0 : 100.0;
    }
    // 3lep0tau's M_ll1+M_ll2 check is a real cut (>=250) in the source script;
    // every other case (all 4-lepton regions, 3lep1tau, 3lep2tau) leaves it
    // hardcoded true.
    const bool mll1Plus2Cut = (!is4Lep && Ntau == 0) ? (mll1 + mll2 >= 250.0) : true;

    const std::map<std::string, bool> cutbools = {
        {"ST", st >= stThreshold}, {"dRll", true}, {"M_ll1", true}, {"M_ll1+M_ll2 >= 500", mll1Plus2Cut}, {"dRlplm", true}
    };
    bool pass = true;
    for (const auto& kv : cutbools) pass = pass && kv.second;
    if (!pass) return out;

    out.eligible = true;
    out.region = region;
    out.mll1 = mll1; out.mll2 = mll2;
    out.catNum = catToNumber(realcat);
    out.masses = masses;
    return out;
}

// Fills one entry into the region-keyed signal-shape tree from an already
// -computed SignalShapeVariables, under the given weight and key suffix
// (e.g. "_roccorUp", or "" for nominal). Cheap: no reconstruction happens
// here, so this is safe to call once per weight-only systematic variant
// sharing the same (nominal) SignalShapeVariables. No histograms -- the
// discriminant lives in the tree branches; Combine only ever needs those,
// binned into a histogram per systematic at datacard-build time, and the
// main hist_<sample>.root output separately already has everything else
// this used to duplicate.
inline void fillSignalShapeVariant(const SignalShapeVariables& c, double weight, const std::string& keySuffix,
                                    SignalShapeHistBundle& sh) {
    if (!c.eligible) return;
    const std::string region = c.region + keySuffix;

    TTree* t = sh.getTree(region);
    SignalShapeTreeVars& v = *sh.treeVars[region];
    v.mll1 = c.mll1; v.mll2 = c.mll2;
    v.mH1 = c.masses.Mdch01; v.mH2 = c.masses.Mdch02;
    v.mDCH1 = c.masses.MdchFit1; v.mDCH2 = c.masses.MdchFit2;
    v.evtwt = weight; v.cat = c.catNum; v.gencat = gen_cat; v.region = c.region;
    t->Fill();
}
