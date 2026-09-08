template <typename Chi2Func>
static bool minimize_in_y(
    vec3& y,
    Chi2Func chi2_of_y,
    double& bestChi2,
    int maxIter = 200,
    double epsGrad = 1e-4
){
    const double eps = 1e-4;  // finite-diff step in y
    double f = chi2_of_y(y);
    bestChi2 = f;

    for (int it=0; it<maxIter; ++it){
        // numerical gradient in y
        vec3 g{0,0,0};
        for(int k=0;k<3;k++){
            vec3 yp = y;
            yp[k] += eps;
            double fp = chi2_of_y(yp);
            g[k] = (fp - f)/eps;
        }

        double gn = std::sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
        if (gn < epsGrad) return true;

        // backtracking line search
        double step = 0.5;
        bool accepted = false;
        for(int ls=0; ls<20; ++ls){
            vec3 yn = { y[0]-step*g[0], y[1]-step*g[1], y[2]-step*g[2] };
            double fn = chi2_of_y(yn);
            if (fn < f){
                y = yn;
                f = fn;
                if (f < bestChi2) bestChi2 = f;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if (!accepted) return false;
    }
    return true;
}

FitOut Fit3Nu_Option1_MTproxy_YBarrier(
    const std::array<TLorentzVector,4>& lep,
    const TVector3& met,
    const mat2& metcov,
    int idx1,int idx2,int idx3,
    int hA1,int hA2,int hB1,
    double sigma_bal = 20.0,
    double sigma_m   = 15.0,
    double w_bal = 1.0,
    double w_m   = 1.0,
    double mu_barrier = 0.1,    // strength (small!)
    double eps_barrier = 1e-3,   // GeV
    int maxIter = 200
){
    FitOut out; out.ok=false; out.x={0,0,0}; out.chi2=1e99; out.mA=0; out.mB=0;

    mat2 W; if (!inv2(metcov, W)) return out;

	int hB2 = 0+1+2+3 - (idx1 + idx2 + idx3);//leftover idx
	
    // neutrino lepton directions
    std::array<int,3> idx = {idx1,idx2,idx3};
    std::array<vec2,3> uhat;
    for(int k=0;k<3;k++){
        vec2 lT{lep[idx[k]].Px(), lep[idx[k]].Py()};
        double pt = norm(lT);
        if (pt < 1e-9) return out;
        uhat[k] = (1.0/pt)*lT;
    }

    vec2 metObs{met.X(), met.Y()};

    // visible transverse sums
    vec2 VA{ lep[hA1].Px()+lep[hA2].Px(), lep[hA1].Py()+lep[hA2].Py() };
    vec2 VB{ lep[hB1].Px()+lep[hB2].Px(), lep[hB1].Py()+lep[hB2].Py() };

    double mvisA2 = (lep[hA1]+lep[hA2]).M2(); if (mvisA2 < 0) mvisA2 = 0;
    double mvisB2 = (lep[hB1]+lep[hB2]).M2(); if (mvisB2 < 0) mvisB2 = 0;

    auto compute_MT = [&](const vec2& V, double mvis2, const vec2& N, const vec2& H){
        double ETvis = std::sqrt(mvis2 + norm2(V));
        double ETnu  = norm(N);
        double MT2 = (ETvis + ETnu)*(ETvis + ETnu) - norm2(H);
        return std::sqrt(std::max(0.0, MT2));
    };

    auto chi2_of_y = [&](const vec3& y, double* mA=nullptr, double* mB=nullptr)->double{
        vec3 x = y_to_x(y);

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

        // mass proxy
        double mTa = compute_MT(VA, mvisA2, NA, HA);
        double mTb = compute_MT(VB, mvisB2, NB, HB);
        double chi_m = (mTa - mTb)*(mTa - mTb) / (sigma_m*sigma_m);

        // barrier: discourage x -> 0
        double chi_bar = -mu_barrier * ( std::log(x[0] + eps_barrier)
                                       + std::log(x[1] + eps_barrier)
                                       + std::log(x[2] + eps_barrier) );

        if (mA) *mA = mTa;
        if (mB) *mB = mTb;

        return chi_met + w_bal*chi_bal + w_m*chi_m + chi_bar;
    };

    // initialize y from small positive x (~1 GeV)
    vec3 y{0,0,0}; // softplus(0)=0.693 -> x ~ 0.7 GeV seed
    double bestChi2;
    bool ok = minimize_in_y(y,
        [&](const vec3& yy){ return chi2_of_y(yy); },
        bestChi2, maxIter, 1e-4);

    if (!ok) return out;

    vec3 x = y_to_x(y);
    double mA, mB;
    double chi2 = chi2_of_y(y, &mA, &mB);

    out.ok = true;
    out.x = x;
    out.chi2 = chi2;
    out.mA = mA;
    out.mB = mB;
    return out;
}

FitOut Fit3Nu_Option2_FullMass_YBarrier(
    const std::array<TLorentzVector,4>& lep,
    const TVector3& met,
    const mat2& metcov,
    int idx1,int idx2,int idx3,
    int hA1,int hA2,int hB1,
    double sigma_bal = 15,
    double sigma_m   = 5.0,
    double w_bal = 1.0,
    double w_m   = 5.0,
    double mu_barrier = 0.1,
    double eps_barrier = 1e-2,
    int maxIter = 200
){
    FitOut out; out.ok=false; out.x={0,0,0}; out.chi2=1e99; out.mA=0; out.mB=0;

    mat2 W; if (!inv2(metcov, W)) return out;
    
	int hB2 = 0+1+2+3 - (idx1 + idx2 + idx3);//leftover idx
	
    std::array<int,3> idx = {idx1,idx2,idx3};
    std::array<vec2,3> uhat;
    std::array<double,3> ptL;

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
        TLorentzVector nu = lep[lepIndex];
        nu *= alpha;
        return nu;
    };

    auto chi2_of_y = [&](const vec3& y, double* mA=nullptr, double* mB=nullptr)->double{
        vec3 x = y_to_x(y);

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

        // barrier
        double chi_bar = -mu_barrier * ( std::log(x[0] + eps_barrier)
                                       + std::log(x[1] + eps_barrier)
                                       + std::log(x[2] + eps_barrier) );

        if (mA) *mA = mHa;
        if (mB) *mB = mHb;

        return chi_met + w_bal*chi_bal + w_m*chi_m + chi_bar;
    };

    vec3 y{0,0,0}; // seed
    double bestChi2;
    bool ok = minimize_in_y(y,
        [&](const vec3& yy){ return chi2_of_y(yy); },
        bestChi2, maxIter, 1e-4);

    if (!ok) return out;

    vec3 x = y_to_x(y);
    double mA, mB;
    double chi2 = chi2_of_y(y, &mA, &mB);

    out.ok = true;
    out.x = x;
    out.chi2 = chi2;
    out.mA = mA;
    out.mB = mB;
    return out;
}
