// met_reco_with_dotcheck_fastpath.cc
// Compile:
// g++ -O2 -std=c++17 `root-config --cflags --libs` met_reco_with_dotcheck_fastpath.cc -o met_reco_with_dotcheck_fastpath

#include <iostream>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMinuit.h"
#include "TMath.h"

// ---------- Globals for Minuit (alphas fit) ----------
static std::array<TLorentzVector,4> gLep_ref, gLep_lab;
static TVector2 gMET_ref;
static TVector2 g_uT[4];
static double g_reg_alpha = 1e-3; // tiny regularizer weight

// alpha_i = x_i^2 to enforce non-negativity; Minuit params are x_i.
void fcn_alphas(int &npar, double *, double &f, double *x, int) {
    double a[4];
    for (int i=0;i<4;++i) a[i] = x[i]*x[i];

    // MET closure term (GeV^2)
    TVector2 met_fit(0,0);
    for (int i=0;i<4;++i) met_fit += a[i] * g_uT[i];
    TVector2 r = gMET_ref - met_fit;
    double mismatch = r.Mod2();
	
	if (mismatch > 1e-9) {
        f = 1e12 * mismatch;  // forbid non-closure
        return;
    }
    // Regularizer to pick small alphas among continuum
    double term_reg = 0;
    for (int i=0;i<4;++i) term_reg += a[i]*a[i];

    f = mismatch + g_reg_alpha *term_reg;
}

// ---------- Helpers ----------
TLorentzVector make_collinear_nu_ref(const TLorentzVector &lep_ref, double a) {
    // a is transverse magnitude for neutrino in REF frame (GeV)
    double ptL = lep_ref.Pt();
    if (ptL <= 1e-12 || a <= 0) return TLorentzVector(0,0,0,0);
    double k = a / ptL;
    TVector3 pnu = k * lep_ref.Vect(); // collinear in 3D
    return TLorentzVector(pnu.X(), pnu.Y(), pnu.Z(), pnu.Mag()); // massless
}

// rapidity of a TLorentzVector (works for massless or massive)
double rapidity(const TLorentzVector &v) {
    double E = v.E(), pz = v.Pz();
    if (E <= fabs(pz) + 1e-12) return 0.0;
    return 0.5 * log((E + pz) / (E - pz));
}

// Signature:
// returns pair<MH1, MH2>. If zero_opposite==true, any neutrino whose 3-vector
// dot lepton_3vec <= eps will be treated as zero for the mass computation.
// Optionally fills zeroedIndices (size 4) and nu_lab_out (size 4).
std::pair<double,double> masses_from_alphas_with_pzpolicy(
    const double a_ref[4],
    const std::string &pz_policy,
    const std::vector<double> &s_fit, // for "fit" policy; size==4 expected
    bool zero_opposite,
    std::array<bool,4>* zeroedIndices = nullptr,
    std::array<TLorentzVector,4>* nu_lab_out = nullptr
) {
    // build neutrinos in REF
    TLorentzVector nu_ref[4];
    for (int i=0;i<4;++i) {
        double pt_n = a_ref[i];
        if (pt_n <= 0) { nu_ref[i].SetPxPyPzE(0,0,0,0); continue; }
        double ux = g_uT[i].X(), uy = g_uT[i].Y();
        double px = pt_n * ux, py = pt_n * uy;
        double pz = 0.0;
        if (pz_policy == "collinear") {
            double ptL = gLep_ref[i].Pt();
            double pzL = gLep_ref[i].Pz();
            double k = (ptL>1e-12) ? (pt_n / ptL) : 0.0;
            pz = k * pzL;
        } else if (pz_policy == "equal_rapidity") {
            double yL = rapidity(gLep_ref[i]);
            pz = pt_n * std::sinh(yL);
        } else if (pz_policy == "fit") {
            if ((int)s_fit.size() == 4) {
                double pzL = gLep_ref[i].Pz();
                pz = s_fit[i] * pzL;
            } else pz = 0.0;
        } else {
            // fallback to collinear
            double ptL = gLep_ref[i].Pt();
            double pzL = gLep_ref[i].Pz();
            double k = (ptL>1e-12) ? (pt_n / ptL) : 0.0;
            pz = k * pzL;
        }
        double E = std::sqrt(px*px + py*py + pz*pz);
        nu_ref[i].SetPxPyPzE(px,py,pz,E);
    }

    // If requested, zero neutrinos that are anti-aligned with their lepton:
    std::array<bool,4> zeroed = {false,false,false,false};
    if (zero_opposite) {
        const double eps = 1e-8;
        for (int i=0;i<4;++i) {
            if (nu_ref[i].P() <= 0) { zeroed[i] = true; nu_ref[i].SetPxPyPzE(0,0,0,0); continue; }
            double dot = nu_ref[i].Vect().Dot(gLep_ref[i].Vect());
            if (dot <= eps) { // not pointing in same direction -> discard neutrino
                zeroed[i] = true;
                nu_ref[i].SetPxPyPzE(0,0,0,0);
            }
        }
        if (zeroedIndices) *zeroedIndices = zeroed;
    }

    // compute betaZ to boost back to lab
    TLorentzVector Pvis_lab(0,0,0,0);
    for (int i=0;i<4;++i) Pvis_lab += gLep_lab[i];
    double betaZ = (std::abs(Pvis_lab.E())>1e-12) ? (-Pvis_lab.Pz()/Pvis_lab.E()) : 0.0;

    // boost nu_ref -> nu_lab and optionally fill nu_lab_out
    TLorentzVector nu_lab[4];
    for (int i=0;i<4;++i) { nu_lab[i] = nu_ref[i]; nu_lab[i].Boost(0,0,-betaZ); }
    if (nu_lab_out) for (int i=0;i<4;++i) (*nu_lab_out)[i] = nu_lab[i];

    // build H1=(0,1), H2=(2,3)
    TLorentzVector H1 = gLep_lab[0] + gLep_lab[1] + nu_lab[0] + nu_lab[1];
    TLorentzVector H2 = gLep_lab[2] + gLep_lab[3] + nu_lab[2] + nu_lab[3];

    return { H1.M(), H2.M() };
}

// ---------- Fit for pz scale factors (fit policy) ----------
static std::array<double,4> g_alphas_for_pzfit; // store alphas for the fit FCN
static double g_pzfit_reg = 1e-2;
void fcn_pzfit(int &npar, double*, double &f, double *par, int) {
    std::vector<double> s(4,0.0);
    for (int i=0;i<4;++i) s[i] = par[i];
    auto [m1,m2] = masses_from_alphas_with_pzpolicy(g_alphas_for_pzfit.data(), "fit", s, /*zero_opposite=*/false, nullptr, nullptr);
    double term_massdiff = (m1 - m2)*(m1 - m2);
    double term_reg = 0.0;
    for (int i=0;i<4;++i) term_reg += s[i]*s[i];
    f = term_massdiff + g_pzfit_reg * term_reg;
}

// ---------- NEW: analytic solve for up to 2 candidates ----------
// Solve single candidate: alpha = u·MET (u is unit)
static double solve_single_alpha(const TVector2 &u, const TVector2 &MET) {
    return u.X()*MET.X() + u.Y()*MET.Y(); // since u unit, dot = alpha
}

// Solve two-vector system: a1*u1 + a2*u2 = MET (u1,u2 not necessarily orth)
// Return (ok, a1, a2, residual)
// If system nearly singular, returns ok=false.
static std::tuple<bool,double,double,double>
solve_two_alphas(const TVector2 &u1, const TVector2 &u2, const TVector2 &MET) {
    const double u1u1 = u1.Mod2(); // should be 1 if unit, but keep general
    const double u2u2 = u2.Mod2();
    const double u1u2 = u1.X()*u2.X() + u1.Y()*u2.Y();
    const double b1 = u1.X()*MET.X() + u1.Y()*MET.Y();
    const double b2 = u2.X()*MET.X() + u2.Y()*MET.Y();

    const double det = u1u1*u2u2 - u1u2*u1u2;
    if (std::abs(det) < 1e-8) return {false, 0.0, 0.0, 1e9};

    double a1 = ( u2u2*b1 - u1u2*b2) / det;
    double a2 = (-u1u2*b1 + u1u1*b2) / det;

    TVector2 fit = a1*u1 + a2*u2;
    double residual = (MET - fit).Mod();
    return {true, a1, a2, residual};
}
/*static std::tuple<bool,double,double,double>
solve_two_alphas(const TVector2 &u1, const TVector2 &u2, const TVector2 &MET) {
 double denom = u1.X()*u2.Y()-u2.X()*u1.Y();
    if (fabs(denom) < 1e-6) {
        std::cerr << "Lepton directions nearly collinear! No unique solution.\n";
        return {false ,0.0, 0.0, 1e9};
    }

    double alpha1 = (MET.X()*u2.Y()-u2.X()*MET.Y()) / denom;
    double alpha2 = (u1.X()*MET.Y()-MET.X()*u1.Y()) / denom;
	
	TVector2 fit = alpha1*u1 + alpha2*u2;
    double residual = (MET - fit).Mod();
    
    return {true, alpha1, alpha2, residual};
}*/

// ---------- Main routine: runs full flow for one event ----------
struct Result {
    bool ok = false;
    double alpha[4]={0,0,0,0}; // in REF frame (final, possibly zeroed)
    double residual = 0;
    double MH1 = 0, MH2 = 0;
    TLorentzVector nu_lab[4];
};

Result run_reco_event(const std::array<TLorentzVector,4> &leps_lab, const TVector2 &MET_lab,
                      const std::string &pz_policy = "collinear",
                      bool run_pz_fit = false)
{
    Result out;
    gLep_lab = leps_lab;

    // 1) boost to REF where Pz_vis = 0
    TLorentzVector Pvis_lab(0,0,0,0);
    for (int i=0;i<4;++i) Pvis_lab += gLep_lab[i];
    double betaZ = (std::abs(Pvis_lab.E())>1e-12) ? (-Pvis_lab.Pz()/Pvis_lab.E()) : 0.0;

    gLep_ref = gLep_lab;
    for (int i=0;i<4;++i) gLep_ref[i].Boost(0,0,betaZ);

    // MET_ref unchanged under pure z-boost
    gMET_ref = MET_lab;

    // unit transverse directions
    for (int i=0;i<4;++i) {
        double pt = gLep_ref[i].Pt();
        if (pt > 1e-12) g_uT[i] = TVector2(gLep_ref[i].Px()/pt, gLep_ref[i].Py()/pt);
        else g_uT[i] = TVector2(0,0);
    }

    // 2) quick check: projections of MET along uT
    std::vector<int> candidates;
    std::vector<double> proj(4,0.0);
    for (int i=0;i<4;++i) {
        if (g_uT[i].Mod2() <= 0) { proj[i]=0.0; continue; }
        proj[i] = g_uT[i].X()*gMET_ref.X() + g_uT[i].Y()*gMET_ref.Y(); // u·MET
        if (proj[i] > 1e-3) candidates.push_back(i);
    }

    // If up to 2 candidates, solve analytically
    double alphas[4] = {0,0,0,0};
    if ((int)candidates.size() <= 2) {
        if (candidates.empty()) {
            // all zeros -- trivial
            for (int i=0;i<4;++i) alphas[i]=0.0;
        } else if ((int)candidates.size() == 1) {
            int i = candidates[0];
            double a = solve_single_alpha(g_uT[i], gMET_ref);
            alphas[i] = std::max(0.0, a);
        } else {
            int i = candidates[0], j = candidates[1];
            auto [ok, a_i, a_j, res] = solve_two_alphas(g_uT[i], g_uT[j], gMET_ref);
            if (!ok) {
                // nearly collinear or numerical issue: fall back to choosing bigger projection
                double bi = proj[i], bj = proj[j];
                if (bi >= bj) { alphas[i] = std::max(0.0, bi); alphas[j]=0.0; }
                else { alphas[j] = std::max(0.0, bj); alphas[i]=0.0; }
            } else {
                // enforce non-negativity: if any negative, project to single-vector solution
                if (a_i >= 0 && a_j >= 0) {
                    alphas[i] = a_i; alphas[j] = a_j;
                } else if (a_i < 0 && a_j >= 0) {
                    alphas[i] = 0.0; alphas[j] = std::max(0.0, solve_single_alpha(g_uT[j], gMET_ref));
                } else if (a_j < 0 && a_i >= 0) {
                    alphas[j] = 0.0; alphas[i] = std::max(0.0, solve_single_alpha(g_uT[i], gMET_ref));
                } else {
                    // both negative (unlikely) -> pick largest positive projection if any
                    double bi = proj[i], bj = proj[j];
                    if (bi >= bj) { alphas[i] = std::max(0.0, bi); alphas[j]=0.0; }
                    else           { alphas[j] = std::max(0.0, bj); alphas[i]=0.0; }
                }
            }
        }
    } else {
        // 3 or 4 candidates -> use Minuit fit (global alpha fit) as before
        TMinuit minuit_a(4);
        minuit_a.SetFCN(fcn_alphas);
        minuit_a.SetPrintLevel(-1);

        double start[4];
        double step[4] = {0.5,0.5,0.5,0.5};
        double metmag = gMET_ref.Mod();
        for (int i=0;i<4;++i) start[i] = std::sqrt(std::max(0.0, metmag/4.0));

        for (int i=0;i<4;++i) minuit_a.DefineParameter(i, Form("x%d",i), start[i], step[i], 0, 0);

        minuit_a.Migrad();

        double xi, ex;
        for (int i=0;i<4;++i) {
            minuit_a.GetParameter(i, xi, ex);
            alphas[i] = xi*xi;
        }
    }

    // compute residual MET in REF (before dot-check)
    TVector2 met_fit_before(0,0);
    for (int i=0;i<4;++i) met_fit_before += alphas[i] * g_uT[i];
    double resid_before = (gMET_ref - met_fit_before).Mod();

    // 3) If using fit-policy for pz, run small Minuit to find s_i
    std::vector<double> sfit_vec = {};
    if ((int)candidates.size() <= 2) {
		if (pz_policy == "fit" || run_pz_fit) {
		    for (int i=0;i<4;++i) g_alphas_for_pzfit[i] = alphas[i];
		    TMinuit minuit_pz(4);
		    minuit_pz.SetFCN(fcn_pzfit);
		    minuit_pz.SetPrintLevel(-1);
		    for (int i=0;i<4;++i) {
		        double ptL = gLep_ref[i].Pt();
		        double sstart = (ptL>1e-12) ? (alphas[i]/ptL) : 0.0;
		        double sstep = std::max(0.1, fabs(sstart)*0.2 + 0.1);
		        minuit_pz.DefineParameter(i, Form("s%d",i), sstart, sstep, -10.0, 10.0);
		    }
		    minuit_pz.Migrad();
		    double sfit[4], err;
		    for (int i=0;i<4;++i) { minuit_pz.GetParameter(i, sfit[i], err); }
		    sfit_vec.assign(sfit, sfit+4);
		}
	}

    // 4) Final mass computation with dot-check (zero neutrinos anti-aligned with lepton)
    std::array<bool,4> zeroed = {false,false,false,false};
    std::array<TLorentzVector,4> nu_lab_out;
    auto masses = masses_from_alphas_with_pzpolicy(alphas, pz_policy, sfit_vec, /*zero_opposite=*/true, &zeroed, &nu_lab_out);
    out.MH1 = masses.first; out.MH2 = masses.second;

    // Now set final alphas and nu_lab: zero alphas where zeroed==true
    for (int i=0;i<4;++i) {
        out.alpha[i] = zeroed[i] ? 0.0 : alphas[i];
        out.nu_lab[i] = nu_lab_out[i];
    }

    // Recompute residual MET with zeroed alphas
    TVector2 met_fit_after(0,0);
    for (int i=0;i<4;++i) met_fit_after += out.alpha[i] * g_uT[i];
    out.residual = (gMET_ref - met_fit_after).Mod();

    out.ok = true;
    return out;
}
/*// --------------------- Demo main ---------------------
int main(){
    // Toy event: lepton four-vectors in LAB (pt,eta,phi,m)
    std::array<TLorentzVector,4> L;
    L[0].SetPtEtaPhiM(40, 0.2, 0.1, 0.000511);
    L[1].SetPtEtaPhiM(30,-0.3, 1.0, 0.1057);
    L[2].SetPtEtaPhiM(45, 0.4,-2.0, 0.000511);
    L[3].SetPtEtaPhiM(25,-0.5,-1.2, 0.1057);

    TVector2 MET(35, -20);

    // Try collinear with dot-check and fast-path analytic solve
    auto res_col = run_reco_event(L, MET, "collinear", false);
    std::cout << "--- collinear + dot-check + fast-path ---\n";
    for (int i=0;i<4;++i) std::cout << "alpha["<<i<<"]="<<res_col.alpha[i] << "  ";
    std::cout << "\nresidual MET ref = " << res_col.residual << " GeV\n";
    std::cout << "MH1="<<res_col.MH1<<"  MH2="<<res_col.MH2<<"\n";

    // Try fit-based pz adjustment with dot-check (still uses analytic path if <=2 candidates)
    auto res_fit = run_reco_event(L, MET, "fit", true);
    std::cout << "\n--- fit + dot-check + fast-path ---\n";
    for (int i=0;i<4;++i) std::cout << "alpha["<<i<<"]="<<res_fit.alpha[i] << "  ";
    std::cout << "\nresidual MET ref = " << res_fit.residual << " GeV\n";
    std::cout << "MH1="<<res_fit.MH1<<"  MH2="<<res_fit.MH2<<"\n";

    return 0;
}*/
