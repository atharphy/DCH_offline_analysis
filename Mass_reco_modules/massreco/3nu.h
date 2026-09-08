#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"

struct AlphaSolution {
    bool ok;
    double alpha1, alpha2, alpha3;
    double chi2;
    double mHA, mHB;
};

AlphaSolution SolveThreeNuMETCov(
    const TLorentzVector& lep1,
    const TLorentzVector& lep2,
    const TLorentzVector& lep3,
    const TVector3& MET,
    const TMatrixD& SigmaMET,   // 2x2
    double m2_ll1,
    double m2_ll2,
    /*double alp1, double alp2,*/
    int maxIter = 20,
    double tol = 1e-6
) {
    AlphaSolution sol;
    sol.ok = false;
    TLorentzVector l1=lep1, l2=lep2, l3=lep3;
	//double boost = (LepV(1)+LepV(2)+LepV(3)+LepV(4)).Pz()/(LepV(1)+LepV(2)+LepV(3)+LepV(4)).E();
	//double boost = (l1+l2+l3).Pz()/(l1+l2+l3).E();
	//l1.Boost(0,0,-boost);l2.Boost(0,0,-boost);l3.Boost(0,0,-boost);
    // ----------------------------------------------------
    // Section 3: pz constraint → eliminate alpha3
    // ----------------------------------------------------
    if (std::abs(l3.Pz()) < 1e-6) return sol;

    double c1 = -l1.Pz() / l3.Pz();
    double c2 = -l2.Pz() / l3.Pz();

    // ----------------------------------------------------
    // Section 4: reduced MET system
    // q1 = k1T + c1 k3T, q2 = k2T + c2 k3T
    // ----------------------------------------------------
    TVector2 q1(l1.Px() + c1*l3.Px(),
                l1.Py() + c1*l3.Py());
    TVector2 q2(l2.Px() + c2*l3.Px(),
                l2.Py() + c2*l3.Py());

    TMatrixD Q(2,2);
    Q(0,0) = q1.X(); Q(0,1) = q2.X();
    Q(1,0) = q1.Y(); Q(1,1) = q2.Y();

    double detQ = Q.Determinant();
    if (std::abs(detQ) < 1e-8) return sol;

    TMatrixD Qinv = Q;
    Qinv.Invert();

    // MET-only solution α̂
    TVectorD met(2);
    met[0] = MET.X();
    met[1] = MET.Y();

    TVectorD alphaHat = Qinv * met;

    // Σ_α = Q⁻¹ Σ Q⁻ᵀ
    TMatrixD SigmaAlpha = Qinv * SigmaMET * TMatrixD(TMatrixD::kTransposed, Qinv);
    TMatrixD W = SigmaAlpha;
    W.Invert();

    // ----------------------------------------------------
    // Section 6–7: χ² minimization
    // ----------------------------------------------------
    TVectorD alpha(2);
    //alpha = (alp1,alp2);
	alpha = alphaHat;// initialize at MET-only
	double lambda = 0.0000000001;
    
    for (int iter = 0; iter < maxIter; ++iter) {

        double a1 = alpha[0];
        double a2 = alpha[1];
        double a3 = c1*a1 + c2*a2;

        // Δ(α) = A(1+a1)(1+a2) − B(1 + c1 a1 + c2 a2)
        double Delta =
            m2_ll1 * (1.0 + a1)*(1.0 + a2)
          - m2_ll2 * (1.0 + a3);

        // ∇Δ
        TVectorD gradDelta(2);
        gradDelta[0] = m2_ll1*(1.0 + a2) - m2_ll2*c1;
        gradDelta[1] = m2_ll1*(1.0 + a1) - m2_ll2*c2;
		
        // Gradient g = 2W(α−α̂) + 2λΔ∇Δ
        TVectorD g = 2.0 * (W * (alpha - alphaHat))
                   + 2.0 * lambda * Delta * gradDelta;

        // Hessian
        TMatrixD H = 2.0 * W;
		
        // ∇Δ ∇Δᵀ term
        TMatrixD outer(2,2);
        for (int i=0;i<2;i++)
            for (int j=0;j<2;j++)
                outer(i,j) = gradDelta[i]*gradDelta[j];

        H += 2.0 * lambda * outer;
        //cout<<"HEHEHEHEH "<<H(1,0)<<"\t"<<H(0,1)<<"\t"<<H(0,0)<<"\t"<<H(1,1)<<endl;
        // + Δ ∇²Δ
        TMatrixD HessDelta(2,2);
        HessDelta(0,0)=0; HessDelta(0,1)=m2_ll1;
        HessDelta(1,0)=m2_ll1; HessDelta(1,1)=0;

        H += 2.0 * lambda * Delta * HessDelta;
		//cout<<"GEGEGEGEG "<<H(1,0)<<"\t"<<H(0,1)<<"\t"<<H(0,0)<<"\t"<<H(1,1)<<endl;
        // Newton step
        TMatrixD Hinv = H.Invert();

        if (Hinv.Determinant() < 1e-12) continue;
		//cout<<"Hdet: "<<H.Determinant()<<endl;
        //std::cout << "iter " << iter<< " alpha = " << alpha[0] << ", " << alpha[1]<< std::endl;
          
        TVectorD step = Hinv * g;
        alpha -= step;
        //if (step.Norm2Sqr() < tol*tol) break;
    }

    // ----------------------------------------------------
    // Final solution
    // ----------------------------------------------------
    double a1 = alpha[0];
    double a2 = alpha[1];
    double a3 = c1*a1 + c2*a2;
	cout<<" ALPHAS: "<<a1<<"\t"<<a2<<"\t"<<a3<<endl;
    if (a1 < 0 || a2 < 0 || a3 < 0) return sol;

    sol.alpha1 = a1;
    sol.alpha2 = a2;
    sol.alpha3 = a3;

    double m2HA = m2_ll1*(1+a1)*(1+a2);
    double m2HB = m2_ll2*(1+a3);

    sol.mHA = std::sqrt(std::max(0.0,m2HA));
    sol.mHB = std::sqrt(std::max(0.0,m2HB));

    // χ²
    TVectorD diff = alpha - alphaHat;
    sol.chi2 = diff*(W*diff)
             + lambda*(m2HA - m2HB)*(m2HA - m2HB);
	cout<<"SOL chi2: "<<sol.chi2<<" MASSES: "<<sol.mHA<<"\t"<<sol.mHB<<endl;
    sol.ok = true;
    return sol;
}
