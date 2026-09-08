#pragma once

#include <array>
#include <string>

#include "TLorentzVector.h"
#include "TVector2.h"

#include "massreco/rjm_test.C"   // get_legs_comb, mat2
#include "massreco/rjm_guru.C"   // get_legs_tau
#include "massreco/rjm_final.C"  // run_reco_event, Result

struct SignalShapeMasses {
    double Mdch01 = 0.0, Mdch02 = 0.0;      // get_legs_comb (feeds h_mT1/h_mT2 and the tree's mH1/mH2)
    double Mdch1 = 0.0, Mdch2 = 0.0;        // get_legs_tau (feeds h_mTtot1/h_mTtot2)
    double MdchFit1 = 0.0, MdchFit2 = 0.0;  // run_reco_event "fit" policy (feeds h_mDCH1/h_mDCH2)
};

inline SignalShapeMasses computeSignalShapeMasses(const std::string& realcat, const std::array<TLorentzVector,4>& L,
                                                   TLorentzVector& MET, const mat2& metcov) {
    SignalShapeMasses out;

    auto [nlegs0, Mdch01, Mdch02, bestdPhi, bestChi2, isOpp] = get_legs_comb(realcat, L, MET, metcov, false);
    out.Mdch01 = Mdch01;
    out.Mdch02 = Mdch02;

    auto [nlegs, Mdch1, Mdch2] = get_legs_tau(realcat, L, MET);
    out.Mdch1 = Mdch1;
    out.Mdch2 = Mdch2;

    TVector2 met_xy(MET.X(), MET.Y());
    Result res_fit = run_reco_event(L, met_xy, "fit", true);
    out.MdchFit1 = res_fit.MH1;
    out.MdchFit2 = res_fit.MH2;

    return out;
}
