#include <array>
#include <cmath>
#include <iostream>

using vec3 = std::array<double,3>;
using mat3 = std::array<std::array<double,3>,3>;

void BoostLeptonsZ(
    std::array<TLorentzVector, 4>& leptons,
    double betaZ
) {
    TVector3 boost(0,0,-betaZ);
    for (auto& l : leptons)
        l.Boost(boost);
}

std::tuple<double,double, double> FitAlphasNuPzZero(
    const std::array<TLorentzVector, 4>& leptons_in, // size = 4
    const TVector3& met,
    const mat2& metcov,
    double m2_ll_A, double m2_ll_B,
    vec3& alpha, // initial guess IN, solution OUT
    double& chi2min,
    int idx1,int idx2, int idx3,
    double lambda = 1.0,
    int maxIter = 5
) {
	TMatrixD SigmaMET(2,2);
	SigmaMET(0,0)=metcov[0][0];
	SigmaMET(0,1)=metcov[0][1];
	SigmaMET(1,0)=metcov[1][0];
	SigmaMET(1,1)=metcov[1][1];
    constexpr double eps = 1e-5;
    //=====================INIT alphas======
	int idx4 = 0+1+2+3 - (idx1 + idx2 + idx3);//leftover idx
	TVector3 Lxy[4];
	for (int i=0;i<4;++i) {
		TVector3 v(leptons_in[i].Vect().X(), leptons_in[i].Vect().Y(), 0);
		Lxy[i] = v;
	}
				
	double A1 = Lxy[idx1].Dot(Lxy[idx3]);
	double A2 = Lxy[idx2].Dot(Lxy[idx3]);
	double K = Lxy[idx3].Dot(Lxy[idx3]);
	double C1 = (Lxy[idx1].Cross(Lxy[idx3])).Mag();
	double C2 = (Lxy[idx2].Cross(Lxy[idx3])).Mag();
	double N3 = (met.Cross(Lxy[idx3])).Mag();
	double S3 = met.Dot(Lxy[idx3]);
	
	double r = m2_ll_A/m2_ll_B;
	double A = C1*K*r;
	double B = A2*C1-C2*A1-K*N3*r;
	double C = -A2*N3+C2*(K*(1-r)+S3);
			
	alpha[0] = (-B+sqrt(B*B-4*A*C))*0.5/A;
	alpha[1] = (N3-alpha[0]*C1)/C2;
	alpha[2] = r*(1+alpha[0])*(1+alpha[1])-1;
	//alpha = {1,1,1};
	//======================================
	double betaZ = 0;
    for (int iter = 0; iter < maxIter; ++iter) {
        // ---- 1) compute neutrino boost    
        double Pznu = alpha[0]*leptons_in[idx1].Pz() + alpha[1]*leptons_in[idx2].Pz() + alpha[2]*leptons_in[idx3].Pz();
        double Enu = alpha[0]*leptons_in[idx1].E() + alpha[1]*leptons_in[idx2].E() + alpha[2]*leptons_in[idx3].E();
		
		if(Enu > 1e-6) {
        	betaZ = Pznu/Enu;
		}
		//else betaZ = 0;
		else betaZ = (leptons_in[idx1].Pz()+leptons_in[idx3].Pz()+leptons_in[idx3].Pz()+leptons_in[idx4].Pz())/(leptons_in[idx1].E() + leptons_in[idx2].E() + leptons_in[idx3].E() + leptons_in[idx4].E());
        // ---- 2) boost leptons
        std::array<TLorentzVector, 4> leptons = leptons_in;
        BoostLeptonsZ(leptons, betaZ);

		AlphaSolution sol = SolveThreeNuMETCov(
					leptons[idx1], leptons[idx2], leptons[idx3],
					met,
					SigmaMET,
					m2_ll_A,
					m2_ll_B,
					alpha[0],
					alpha[1]
				);
		//if (sol.alpha1 < 0 || sol.alpha2 < 0 || sol.alpha3 < 0) continue;
		alpha = {sol.alpha1, sol.alpha2, sol.alpha3};
		if(!sol.ok) continue;
        if (sol.chi2 >= chi2min) continue;
		chi2min = sol.chi2;
		alpha = {sol.alpha1, sol.alpha2, sol.alpha3};
		cout<<"betaZ: "<<betaZ<<" CHI2 "<<chi2min<<" ALP: "<<alpha[0]<<"\t"<<alpha[1]<<"\t"<<alpha[2]<<endl;
	}
	double MDCH1 = sqrt(m2_ll_A*(1+alpha[0])*(1+alpha[1]));
	double MDCH2 = sqrt(m2_ll_B*(1+alpha[2]));
	cout<<"mH1: "<<MDCH1<<"\t mH2: "<<MDCH2<<endl;
    return {MDCH1, MDCH2, betaZ};
}
