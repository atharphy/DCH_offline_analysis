#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMinuit.h"
#include "TMath.h"
#include "fir_alphas_metcov.h"
#include "fir_alphas_metcovGN.h"
#include "3nu_new.h"
#include "3nu_good.h"
#include "3nu_good1.h"

std::tuple<int, double, double, double, double, bool> get_legs_comb(std::string cat, const std::array<TLorentzVector,4> &Lep_lab, TLorentzVector &MET_lab, const mat2 &metcov, bool doPrint=false){
	std::cout << std::fixed;
    std::cout << std::setprecision(6);

	double betaZ = 0, sumE = 0;
	std::array<TLorentzVector,4> Lep_ref;
	TVector3 MET_ref; TVector3 Lxy[4];
	double alpha[4] = {0,0,0,0};
	for (int i=0;i<4;++i) {
		sumE += Lep_lab[i].E();
		betaZ += Lep_lab[i].Pz();
	}
	betaZ = betaZ/sumE;
	for (int i=0;i<4;++i) {
		Lep_ref[i] = Lep_lab[i];

		MET_ref = MET_lab.Vect();
		TVector3 v(Lep_ref[i].Vect().X(), Lep_ref[i].Vect().Y(), 0);

		Lxy[i] = v;
	}
	std::vector<int> ileg;
	for (int i=0;i<4;++i) {
        if (cat[i]=='t') ileg.push_back(i);
    }
    double mll_1 = (Lep_ref[0]+Lep_ref[1]).M();
    double mll_2 = (Lep_ref[2]+Lep_ref[3]).M();

    double bestChi2 = 1e15; double best_Chi2 = 1e15; double bestchi2 = 1e15; double chee2 = 999999999;
	double chi2min; double chi2_mass;
	double M2star;
	double bestdPhi = 999;

	float dmH = 99999;
	bool isOpp = true;
	AlphaSolution sol; FitOut out;
    for (int i=0;i<4;i++){
    	for (int j=i+1;j<4;j++){

    		bool i_in = std::find(ileg.begin(), ileg.end(), i) != ileg.end();
            bool j_in = std::find(ileg.begin(), ileg.end(), j) != ileg.end();

			double alpha1 = MET_ref.Cross(Lxy[j].Unit())*Lxy[i].Unit().Cross(Lxy[j].Unit())/Lxy[i].Unit().Cross(Lxy[j].Unit()).Mag2();
			double alpha2 = MET_ref.Cross(Lxy[i].Unit())*Lxy[j].Unit().Cross(Lxy[i].Unit())/Lxy[j].Unit().Cross(Lxy[i].Unit()).Mag2();

			alpha1 = alpha1/Lxy[i].Mag(); alpha2 = alpha2/Lxy[j].Mag();

			double mHa =  -99, mHb = -999;

			vec2 alpha_star = {{0,0}};
			mat2 cov_alpha_star = {{{0, 0},{0, 0}}};
			vec2 alpha_init = {{alpha1,alpha2}};
			if ((i==0 and j==2) or (i==0 and j==3) or (i==1 and j==2) or (i==1 and j==3)){

				bool ok1 = FitAlphasWithMETCov(
				Lxy[i], Lxy[j],
				MET_ref,
				metcov,
				(Lep_ref[0]+Lep_ref[1]).M2(), (Lep_ref[2]+Lep_ref[3]).M2(),
				alpha1, alpha2,
				alpha_star,
				cov_alpha_star,
				chi2min,
				M2star
				);if(!ok1) continue;

				if (alpha_star[0] <= 0 or alpha_star[1] <= 0) continue;

				if (chi2min > bestChi2) continue;
				bestChi2 = chi2min; bestdPhi = deltaPhi(Lep_ref[i], Lep_ref[j]);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = alpha_star[0];
				alpha[j] = alpha_star[1];

				mHa = sqrt(1+alpha_star[0])*mll_1;
				mHb = sqrt(1+alpha_star[1])*mll_2;

				isOpp = true;
			}
			else if (i==0 and j==1){
				bool ok2 = FitAlphasMultiplicativeGN(
				Lxy[i], Lxy[j],
				MET_ref,
				metcov,
				mll_1, mll_2,
				alpha_star,
				cov_alpha_star,
				chi2min,
				alpha_init,
				20,
				1e-6
				);if(!ok2) continue;

				if (alpha_star[0] <= 0 or alpha_star[1] <= 0) continue;

				if (chi2min > bestChi2) continue;
				bestChi2 = chi2min;bestdPhi = deltaPhi(Lep_ref[i], Lep_ref[j]);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = alpha_star[0];
				alpha[j] = alpha_star[1];

				mHa = sqrt((1+alpha_star[0])*(1+alpha_star[1]))*mll_1;
			 	mHb = mll_2;

			 	isOpp = false;
			}
			else if (i==2 and j==3){
				bool ok2 = FitAlphasMultiplicativeGN(
				Lxy[i], Lxy[j],
				MET_ref,
				metcov,
				mll_2, mll_1,
				alpha_star,
				cov_alpha_star,
				chi2min,
				alpha_init,
				20,
				1e-6
				);if(!ok2) continue;

				if (alpha_star[0] <= 0 or alpha_star[1] <= 0) continue;

				if (chi2min > bestChi2) continue;
				bestChi2 = chi2min;bestdPhi = deltaPhi(Lep_ref[i], Lep_ref[j]);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = alpha_star[0];
				alpha[j] = alpha_star[1];

				mHa = mll_1;
				mHb = sqrt((1+alpha_star[0])*(1+alpha_star[1]))*mll_2;

				isOpp = false;
			}

    	}
    }

    TLorentzVector nu[4]; int nlegs =0;
    for (int i=0;i<4;++i) {

    	nu[i] = alpha[i]*Lep_ref[i];
    	if (nu[i].E()>0) nlegs +=1;
    }

    double MDCH1 = (nu[0]+nu[1]+Lep_ref[0]+Lep_ref[1]).M();
	double MDCH2 = (nu[2]+nu[3]+Lep_ref[2]+Lep_ref[3]).M();

	if (doPrint){
		for (int i=0;i<4;++i) {cout<<"nu"<<i<<": "; print4Vec(nu[i]);}
		for (int i=0;i<4;++i) {cout<<"Lep"<<i<<": "; print4Vec(Lep_ref[i]);}
		cout<<"MET "<<MET_ref.Mag()<<"\t"<<"Sum (nu pT) "<<(nu[0]+nu[1]+nu[2]+nu[3]).Pt()<<endl;
		cout<<"Mll "<< mll_1 <<"\t"<< mll_2 <<endl;
		cout<<"DCH masses "<< MDCH1 <<"\t"<< MDCH2<<endl;
		cout<<"SOL masses "<< sol.mHA <<"\t"<< sol.mHB<<endl;
	}

    return {nlegs, MDCH1, MDCH2, bestdPhi, bestChi2, isOpp};
}
