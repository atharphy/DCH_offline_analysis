#pragma once

#include <array>
#include <string>

#include "TLorentzVector.h"
#include "TVector2.h"

#include "massreco/rjm_test.C"
#include "massreco/rjm_guru.C"
#include "massreco/rjm_final.C"

struct SignalShapeMasses {
    double Mdch01 = 0.0, Mdch02 = 0.0;
    double Mdch1 = 0.0, Mdch2 = 0.0;
    double MdchFit1 = 0.0, MdchFit2 = 0.0;
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
