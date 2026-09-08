#include <array>
#include <cmath>
#include <algorithm>

using vec2 = std::array<double,2>;
using vec3 = std::array<double,3>;
using mat2 = std::array<std::array<double,2>,2>;

inline vec2 v2(double x, double y){ return {x,y}; }
inline vec2 operator+(const vec2&a,const vec2&b){ return {a[0]+b[0], a[1]+b[1]}; }
inline vec2 operator-(const vec2&a,const vec2&b){ return {a[0]-b[0], a[1]-b[1]}; }
inline vec2 operator*(double s, const vec2&v){ return {s*v[0], s*v[1]}; }
//inline double dot(const vec2&a,const vec2&b){ return a[0]*b[0] + a[1]*b[1]; }
inline double norm2(const vec2&v){ return dot(v,v); }
inline double norm(const vec2&v){ return std::sqrt(norm2(v)); }

inline bool inv2(const mat2& S, mat2& Sinv){
    double det = S[0][0]*S[1][1] - S[0][1]*S[1][0];
    if (std::fabs(det) < 1e-14) return false;
    double id = 1.0/det;
    Sinv = {{
        {  S[1][1]*id, -S[0][1]*id },
        { -S[1][0]*id,  S[0][0]*id }
    }};
    return true;
}

/*inline vec2 mat2_vec(const mat2& M, const vec2& v){
    return { M[0][0]*v[0] + M[0][1]*v[1],
             M[1][0]*v[0] + M[1][1]*v[1] };
}*/

// stable softplus
inline double softplus(double y){
    if (y > 30.0)  return y;            // log(1+exp(y)) ~ y
    if (y < -30.0) return std::exp(y);  // log(1+exp(y)) ~ exp(y)
    return std::log1p(std::exp(y));
}

// x(y) and derivative dx/dy = sigmoid(y)
inline double sigmoid(double y){
    if (y >= 0) {
        double e = std::exp(-y);
        return 1.0 / (1.0 + e);
    } else {
        double e = std::exp(y);
        return e / (1.0 + e);
    }
}

inline vec3 y_to_x(const vec3& y){
    return { softplus(y[0]), softplus(y[1]), softplus(y[2]) };
}

// map global lepton index -> x slot (0..2), or -1 if no neutrino
inline int nu_slot(int lepIdx, int idx1, int idx2, int idx3){
    if (lepIdx==idx1) return 0;
    if (lepIdx==idx2) return 1;
    if (lepIdx==idx3) return 2;
    return -1;
}

inline vec2 met_pred_from_x(const vec3& x, const std::array<vec2,3>& uhat){
    return x[0]*uhat[0] + x[1]*uhat[1] + x[2]*uhat[2];
}

inline vec2 sumNuT_inPair(const vec3& x, const std::array<vec2,3>& uhat,
                          int a, int b, int idx1, int idx2, int idx3){
    vec2 s{0,0};
    int sa = nu_slot(a,idx1,idx2,idx3);
    int sb = nu_slot(b,idx1,idx2,idx3);
    if (sa>=0) s = s + x[sa]*uhat[sa];
    if (sb>=0) s = s + x[sb]*uhat[sb];
    return s;
}

struct FitOut {
    bool ok;
    vec3 x;     // neutrino pT magnitudes along lepton directions (GeV)
    double chi2;
    double mA;  // proxy/full mass (GeV)
    double mB;  // proxy/full mass (GeV)
};


FitOut Fit3Nu_Option1_MTproxy(
    const std::array<TLorentzVector,4>& lep,
    const TVector3& met,
    const mat2& metcov,
    int idx1,int idx2,int idx3,          // 3 leptons that carry neutrinos
    int hA1,int hA2,int hB1,      // pairing into Higgs A (hA1,hA2) and Higgs B (hB1,hB2)
    double sigma_bal = 20.0,             // GeV, tune from |pT(HH)| proxy
    double sigma_m   = 15.0,             // GeV, tune from expected mass closeness
    double w_bal = 1.0,
    double w_m   = 1.0,
    int maxIter = 200
){
    FitOut out; out.ok=false; out.x = {0,0,0}; out.chi2=1e99; out.mA=0; out.mB=0;

    mat2 W; if (!inv2(metcov, W)) return out;

    // uhat for neutrino leptons
    std::array<int,3> idx = {idx1,idx2,idx3};
    std::array<vec2,3> uhat;
    for(int k=0;k<3;k++){
        vec2 lT{lep[idx[k]].Px(), lep[idx[k]].Py()};
        double pt = norm(lT);
        if (pt < 1e-9) return out;
        uhat[k] = (1.0/pt)*lT;
    }
	int hB2 = 0+1+2+3 - (idx1 + idx2 + idx3);//leftover idx
    vec2 metObs{met.X(), met.Y()};

    // visible transverse sums
    vec2 VA{ lep[hA1].Px()+lep[hA2].Px(), lep[hA1].Py()+lep[hA2].Py() };
    vec2 VB{ lep[hB1].Px()+lep[hB2].Px(), lep[hB1].Py()+lep[hB2].Py() };

    // visible invariant masses for the two dileptons
    double mvisA2 = (lep[hA1]+lep[hA2]).M2();
    double mvisB2 = (lep[hB1]+lep[hB2]).M2();
    if (mvisA2 < 0) mvisA2 = 0;
    if (mvisB2 < 0) mvisB2 = 0;

    auto compute_MT = [&](const vec2& V, double mvis2, const vec2& N, const vec2& H){
        double ETvis = std::sqrt(mvis2 + norm2(V));
        double ETnu  = norm(N);
        double MT2 = (ETvis + ETnu)*(ETvis + ETnu) - norm2(H);
        return std::sqrt(std::max(0.0, MT2));
    };

    auto chi2_total = [&](const vec3& x, double* mA=nullptr, double* mB=nullptr)->double{
        // MET term
        vec2 pred = met_pred_from_x(x, uhat);
        vec2 r = metObs - pred;
        double chi_met = dot(r, mat2_vec(W,r));

        // balance term
        vec2 NA = sumNuT_inPair(x, uhat, hA1,hA2, idx1,idx2,idx3);
        vec2 NB = sumNuT_inPair(x, uhat, hB1,hB2, idx1,idx2,idx3);
        vec2 HA = VA + NA;
        vec2 HB = VB + NB;
        vec2 sumH = HA + HB;
        double chi_bal = norm2(sumH) / (sigma_bal*sigma_bal);

        // mass proxy term
        double mTa = compute_MT(VA, mvisA2, NA, HA);
        double mTb = compute_MT(VB, mvisB2, NB, HB);
        double chi_m = (mTa - mTb)*(mTa - mTb) / (sigma_m*sigma_m);

        if (mA) *mA = mTa;
        if (mB) *mB = mTb;

        return chi_met + w_bal*chi_bal + w_m*chi_m;
    };

    // --- Initialize: solve least squares for MET with x>=0 by simple projection
    // Start with x=0 then do a few projected gradient steps on MET-only
    vec3 x{0,0,0};
    {
        double step = 0.5;
        for(int it=0; it<40; ++it){
            // grad of chi_met is: -2 Uhat^T W r
            vec2 pred = met_pred_from_x(x, uhat);
            vec2 r = metObs - pred;
            vec2 Wr = mat2_vec(W,r);
            vec3 g{
                -2.0*dot(uhat[0], Wr),
                -2.0*dot(uhat[1], Wr),
                -2.0*dot(uhat[2], Wr)
            };
            for(int k=0;k<3;k++){
                x[k] -= step*g[k];
                if (x[k] < 0) x[k] = 0;
            }
        }
    }

    // --- Projected gradient on full chi2 with backtracking
    double f = chi2_total(x, &out.mA, &out.mB);
    out.chi2 = f; out.x = x;

    const double eps = 1e-4;
    double step0 = 0.5;

    for(int it=0; it<maxIter; ++it){
        // numerical gradient of total objective
        vec3 g{0,0,0};
        for(int k=0;k<3;k++){
            vec3 xp = x;
            xp[k] += eps;
            double fp = chi2_total(xp);
            g[k] = (fp - f)/eps;
        }

        // stop if small
        double gn = std::sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
        if (gn < 1e-4) break;

        // backtracking line search with projection
        double step = step0;
        bool accepted = false;
        for(int ls=0; ls<20; ++ls){
            vec3 xn = x;
            for(int k=0;k<3;k++){
                xn[k] = x[k] - step*g[k];
                if (xn[k] < 0) xn[k] = 0;
            }
            double fn = chi2_total(xn);
            if (fn < f){
                x = xn; f = fn;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if (!accepted) break;

        if (f < out.chi2){
            out.chi2 = f; out.x = x;
            chi2_total(x, &out.mA, &out.mB);
        }
    }

    out.ok = true;
    return out;
}

FitOut Fit3Nu_Option2_FullMass(
    const std::array<TLorentzVector,4>& lep,
    const TVector3& met,
    const mat2& metcov,
    int idx1,int idx2,int idx3,
    int hA1,int hA2,int hB1,
    double sigma_bal = 20.0,
    double sigma_m   = 15.0,
    double w_bal = 1.0,
    double w_m   = 1.0,
    int maxIter = 200
){
    FitOut out; out.ok=false; out.x={0,0,0}; out.chi2=1e99; out.mA=0; out.mB=0;

    mat2 W; if (!inv2(metcov, W)) return out;

    std::array<int,3> idx = {idx1,idx2,idx3};
    std::array<vec2,3> uhat;
    std::array<double,3> ptL;
    
	int hB2 = 0+1+2+3 - (idx1 + idx2 + idx3);//leftover idx
	
    for(int k=0;k<3;k++){
        vec2 lT{lep[idx[k]].Px(), lep[idx[k]].Py()};
        double pt = norm(lT);
        if (pt < 1e-9) return out;
        ptL[k] = pt;
        uhat[k] = (1.0/pt)*lT;
    }

    vec2 metObs{met.X(), met.Y()};
    vec2 VA{ lep[hA1].Px()+lep[hA2].Px(), lep[hA1].Py()+lep[hA2].Py() };
    vec2 VB{ lep[hB1].Px()+lep[hB2].Px(), lep[hB1].Py()+lep[hB2].Py() };

    auto buildNu4 = [&](int lepIndex, const vec3& x)->TLorentzVector{
        int s = nu_slot(lepIndex, idx1,idx2,idx3);
        if (s < 0) return TLorentzVector(0,0,0,0);
        double alpha = x[s] / ptL[s];
        // nu = alpha * lepton 4-vector (collinear)
        TLorentzVector nu = lep[lepIndex];
        nu *= alpha;
        return nu;
    };

    auto chi2_total = [&](const vec3& x, double* mA=nullptr, double* mB=nullptr)->double{
        // MET term
        vec2 pred = met_pred_from_x(x, uhat);
        vec2 r = metObs - pred;
        double chi_met = dot(r, mat2_vec(W,r));

        // balance term
        vec2 NA = sumNuT_inPair(x, uhat, hA1,hA2, idx1,idx2,idx3);
        vec2 NB = sumNuT_inPair(x, uhat, hB1,hB2, idx1,idx2,idx3);
        vec2 HA = VA + NA;
        vec2 HB = VB + NB;
        vec2 sumH = HA + HB;
        double chi_bal = norm2(sumH) / (sigma_bal*sigma_bal);

        // full masses
        TLorentzVector HA4 = lep[hA1] + lep[hA2] + buildNu4(hA1,x) + buildNu4(hA2,x);
        TLorentzVector HB4 = lep[hB1] + lep[hB2] + buildNu4(hB1,x) + buildNu4(hB2,x);

        double mHa = HA4.M();
        double mHb = HB4.M();
        double chi_m = (mHa - mHb)*(mHa - mHb) / (sigma_m*sigma_m);

        if (mA) *mA = mHa;
        if (mB) *mB = mHb;

        return chi_met + w_bal*chi_bal + w_m*chi_m;
    };

    // Init with MET-only projected gradient (same as Option 1)
    vec3 x{0,0,0};
    {
        double step = 0.5;
        for(int it=0; it<40; ++it){
            vec2 pred = met_pred_from_x(x, uhat);
            vec2 r = metObs - pred;
            vec2 Wr = mat2_vec(W,r);
            vec3 g{
                -2.0*dot(uhat[0], Wr),
                -2.0*dot(uhat[1], Wr),
                -2.0*dot(uhat[2], Wr)
            };
            for(int k=0;k<3;k++){
                x[k] -= step*g[k];
                if (x[k] < 0) x[k] = 0;
            }
        }
    }

    double f = chi2_total(x, &out.mA, &out.mB);
    out.chi2 = f; out.x = x;

    const double eps = 1e-4;
    double step0 = 0.5;

    for(int it=0; it<maxIter; ++it){
        vec3 g{0,0,0};
        for(int k=0;k<3;k++){
            vec3 xp = x;
            xp[k] += eps;
            double fp = chi2_total(xp);
            g[k] = (fp - f)/eps;
        }
        double gn = std::sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
        if (gn < 1e-4) break;

        double step = step0;
        bool accepted = false;
        for(int ls=0; ls<20; ++ls){
            vec3 xn = x;
            for(int k=0;k<3;k++){
                xn[k] = x[k] - step*g[k];
                if (xn[k] < 0) xn[k] = 0;
            }
            double fn = chi2_total(xn);
            if (fn < f){
                x = xn; f = fn;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if (!accepted) break;

        if (f < out.chi2){
            out.chi2 = f; out.x = x;
            chi2_total(x, &out.mA, &out.mB);
        }
    }

    out.ok = true;
    return out;
}

