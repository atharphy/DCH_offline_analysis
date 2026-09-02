#pragma once
//
// Tau energy scale (TauES) up/down kinematic systematic.
//
// The skim (Online_framework/funcs/outTuple_test.py) writes TauES_Up_1..4 /
// TauES_Down_1..4 -- absolute shifted per-lepton-slot values on the same
// scale as the nominal TauES_1..4 branch, computed with the real per-lepton
// decay mode available at production time (unlike an on-the-fly correction-
// lib recompute at analysis time, which has no per-slot decay-mode branch
// to work with for cat-tight slots -- see sf_uncertainty/LeptonTauSFUncertainty.h's
// own comments on this exact limitation).
//
// computeTauESRelShift() turns those into a per-slot RELATIVE shift; the
// callers in DCH_tauFR.C / DCH_tight.C apply it to `pt` only (mass is
// deliberately left unshifted -- see each call site for why: mass has zero
// downstream effect in either driver's histograms today).
//
// IMPORTANT: relDown is typically NEGATIVE (Down lowers the SF). Any
// "does this slot have a variation" check must be `!= 0.0`, never a
// sign-based `> 0.0` guard (that would silently discard every Down
// variation) -- see ObjectAccessors.h's tauES*ByIndex for the raw values.

#include <algorithm>
#include <string>

#include "ObjectAccessors.h"

inline void computeTauESRelShift(const std::string& cat, double relUp[4], double relDown[4]) {
    const int nLep = static_cast<int>(cat.size());
    for (int idx = 1; idx <= 4; ++idx) {
        relUp[idx - 1] = relDown[idx - 1] = 0.0;
        if (idx > nLep || cat[idx - 1] != 't') continue;
        const double nom = tauESByIndex(idx);
        const double up = tauESUpByIndex(idx);
        const double down = tauESDownByIndex(idx);
        // TauES_Up/Down_N fall back to the -99 "not available" sentinel
        // (same convention as ptByIndex etc.) when the skim couldn't compute
        // a per-slot decay-mode-aware shift for this lepton/sample. Treating
        // that as a genuine SF ratio would produce a ~-10000% "shift",
        // zeroing the tau's pt outright -- so require both variants to be
        // positive (physically sane) before trusting them.
        if (nom <= 0.0 || up <= 0.0 || down <= 0.0) continue;
        relUp[idx - 1] = up / nom - 1.0;
        relDown[idx - 1] = down / nom - 1.0;  // typically NEGATIVE
    }
}

// For DCH_tight.C's flat-global-shift style (mirrors RoccoRCorrections.h's
// shiftMuonPts exactly, but for 't'-flavor slots and a caller-supplied
// per-slot relative array instead of a fixed roccor error array). Shifts
// pt_N in place; caller is responsible for restorePts(origPt) afterward
// (mass is deliberately left unshifted -- see this header's top comment).
inline void shiftTauESPts(const std::string& cat, const double origPt[4], const double relArr[4], double sign) {
    for (int idx = 1; idx <= static_cast<int>(cat.size()) && idx <= 4; ++idx) {
        if (cat[idx - 1] != 't') continue;
        const double factor = std::max(0.0, 1.0 + sign * relArr[idx - 1]);
        const double shifted = origPt[idx - 1] * factor;
        if (idx == 1) pt_1 = shifted;
        else if (idx == 2) pt_2 = shifted;
        else if (idx == 3) pt_3 = shifted;
        else if (idx == 4) pt_4 = shifted;
    }
}
