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
//#include "3nu_best.h"

std::tuple<int, double, double, double, double, bool> get_legs_comb(std::string cat, const std::array<TLorentzVector,4> &Lep_lab, TLorentzVector &MET_lab, const mat2 &metcov, bool doPrint=false){//this gives valid lepton legs with neutrinos
	std::cout << std::fixed;
    std::cout << std::setprecision(6);
    
	//Boosting to have Pz_tot(vis) = 0
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
		//Lep_ref[i].Boost(0,0,-betaZ);
		//print4Vec(Lep_ref[i]);
		MET_ref = MET_lab.Vect(); //no need to boost met
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
	//for 2 alphas
	float dmH = 99999;//abs(mll_1 - mll_2) ;
	bool isOpp = true;
	AlphaSolution sol; FitOut out;
    for (int i=0;i<4;i++){
    	for (int j=i+1;j<4;j++){
    		//if (int(ileg.size())==1){if (ileg[0] != i and ileg[0] != j) continue;}
    		//if (int(ileg.size())==2){if ((ileg[0] != i or ileg[1] != j) and ( ileg[1] != i or ileg[0] != j))continue;}
    		
    		bool i_in = std::find(ileg.begin(), ileg.end(), i) != ileg.end();
            bool j_in = std::find(ileg.begin(), ileg.end(), j) != ileg.end();
            //if (ileg.size()==1 and !(i_in or j_in)) continue;
			//else if (ileg.size()>1 and !(i_in and j_in)) continue;
    		
    		
			double alpha1 = MET_ref.Cross(Lxy[j].Unit())*Lxy[i].Unit().Cross(Lxy[j].Unit())/Lxy[i].Unit().Cross(Lxy[j].Unit()).Mag2();
			double alpha2 = MET_ref.Cross(Lxy[i].Unit())*Lxy[j].Unit().Cross(Lxy[i].Unit())/Lxy[j].Unit().Cross(Lxy[i].Unit()).Mag2();
			//auto [alpha1, alpha2] = solve_two_alphas(Lxy[i].Unit(),Lxy[j].Unit(),MET_ref);
			alpha1 = alpha1/Lxy[i].Mag(); alpha2 = alpha2/Lxy[j].Mag();
			
			//if (doPrint)cout<<"alpha1: "<<sqrt(1+alpha1)*mll_1<<"\t alpha2: "<<sqrt(1+alpha2)*mll_2<<endl;
			//if (alpha1 <= 0 or alpha2 <= 0) continue;
			
			double mHa =  -99, mHb = -999; 
			/*if ((i==0 and j==2) or (i==0 and j==3) or (i==1 and j==2) or (i==1 and j==3)){ 231ALPHAS  -1876470.6993	161999.3888	457010.5042

				mHa = sqrt(1+alpha1)*mll_1;
				mHb = sqrt(1+alpha2)*mll_2;
			}
			else if (i==0 and j==1){
				mHa = sqrt((1+alpha1)*(1+alpha2))*mll_1;
			 	mHb = mll_2;
			}
			else if (i==2 and j==3){
				mHa = mll_1;
				mHb = sqrt((1+alpha1)*(1+alpha2))*mll_2;
			}
			
			cout<<"mH1: "<<mHa<<"\t mH2: "<<mHb<<endl;
			if (abs(mHa - mHb) >= dmH) continue;
			dmH = abs(mHa - mHb);
			for(int l=0; l<4; l++) alpha[l]=0;
			alpha[i] = alpha1;
			alpha[j] = alpha2;*/
			
			vec2 alpha_star = {{0,0}};
			mat2 cov_alpha_star = {{{0, 0},{0, 0}}};
			vec2 alpha_init = {{alpha1,alpha2}};	
			if ((i==0 and j==2) or (i==0 and j==3) or (i==1 and j==2) or (i==1 and j==3)){ 
				//continue;
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
				//if (doPrint) cout<<i<<j<<"\t chi2: "<<chi2min<<" alpha1: "<<alpha_star[0]<<"\t alpha2: "<<alpha_star[1]<<"\t mH(cov): "<<sqrt(M2star)<<" dPhi: "<<deltaPhi(Lep_ref[i], Lep_ref[j])<<endl;
				if (alpha_star[0] <= 0 or alpha_star[1] <= 0) continue;
				
				if (chi2min > bestChi2) continue;
				bestChi2 = chi2min; bestdPhi = deltaPhi(Lep_ref[i], Lep_ref[j]);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = alpha_star[0];
				alpha[j] = alpha_star[1];
				
				mHa = sqrt(1+alpha_star[0])*mll_1;
				mHb = sqrt(1+alpha_star[1])*mll_2;
				//cout<<mHa<<"\t"<<mHb<<endl;
				isOpp = true;
			}
			else if (i==0 and j==1){ //continue;
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
				//if (doPrint) cout<<i<<j<<"\t chi2: "<<chi2min<<" alpha1: "<<alpha_star[0]<<"\t alpha2: "<<alpha_star[1]<<"\t dPhi: "<<deltaPhi(Lep_ref[i], Lep_ref[j])<<endl;
				if (alpha_star[0] <= 0 or alpha_star[1] <= 0) continue;
				
				if (chi2min > bestChi2) continue;
				bestChi2 = chi2min;bestdPhi = deltaPhi(Lep_ref[i], Lep_ref[j]);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = alpha_star[0];
				alpha[j] = alpha_star[1];
				//cout<<"cov mH: "<<sqrt(M2star)<<endl;
				
				mHa = sqrt((1+alpha_star[0])*(1+alpha_star[1]))*mll_1;
			 	mHb = mll_2;
			 	//cout<<mHa<<"\t"<<mHb<<endl;
			 	isOpp = false;
			}
			else if (i==2 and j==3){//continue;
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
				//if (doPrint) cout<<i<<j<<"\t chi2: "<<chi2min<<" alpha1: "<<alpha_star[0]<<"\t alpha2: "<<alpha_star[1]<<"\t dPhi: "<<deltaPhi(Lep_ref[i], Lep_ref[j])<<endl;
				if (alpha_star[0] <= 0 or alpha_star[1] <= 0) continue;
				
				if (chi2min > bestChi2) continue;
				bestChi2 = chi2min;bestdPhi = deltaPhi(Lep_ref[i], Lep_ref[j]);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = alpha_star[0];
				alpha[j] = alpha_star[1];
				//cout<<"cov mH: "<<sqrt(M2star)<<endl;
				
				mHa = mll_1;
				mHb = sqrt((1+alpha_star[0])*(1+alpha_star[1]))*mll_2;
				//cout<<mHa<<"\t"<<mHb<<endl;
				isOpp = false;
			}

			//TESTING 3nu.h
			/*TMatrixD SigmaMET(2,2);
			SigmaMET(0,0)=metcov[0][0];
			SigmaMET(0,1)=metcov[0][1];
			SigmaMET(1,0)=metcov[1][0];
			SigmaMET(1,1)=metcov[1][1];
			for(int k=j+1;k<4;k++){cout<<"INDEX: "<<i<<j<<k<<endl;
				if (i==0 and j==1) sol = SolveThreeNuMETCov(
					Lep_ref[i],Lep_ref[j],Lep_ref[k],
					MET_ref,
					SigmaMET,
					mll_1*mll_1,
					mll_2*mll_2
				); 
				else if (j==2 and k==3) sol = SolveThreeNuMETCov(
					Lep_ref[j],Lep_ref[k],Lep_ref[i],
					MET_ref,
					SigmaMET,
					mll_2*mll_2,
					mll_1*mll_1
				);
				
				if(!sol.ok) continue;
				if (sol.alpha1 <= 0 or sol.alpha2 <= 0 or sol.alpha3 <=0) continue;
				if (sol.chi2 > bestchi2) continue;
				bestchi2 = sol.chi2;
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = sol.alpha1;
				alpha[j] = sol.alpha2;
				alpha[k] = sol.alpha3;
			}*/
			// Testing 3nu_new.h
			/*for(int k=j+1;k<4;k++){cout<<"INDEX: "<<i<<j<<k<<endl;
				bool k_in = std::find(ileg.begin(), ileg.end(), k) != ileg.end();
				//if(!k_in) continue;
				vec3 alp;
				bool okay;
				if (i==0 and j==1)  std::tie(mHa, mHb, betaZ) = FitAlphasNuPzZero(
					Lep_ref, // size = 4
					MET_ref,
					metcov,
					mll_1*mll_1, mll_2*mll_2,
					alp,
					chee2,
					i,j,k
				);
				else if (j==2 and k==3) std::tie(mHb, mHa, betaZ) = FitAlphasNuPzZero(
					Lep_ref, // size = 4
					MET_ref,
					metcov,
					mll_2*mll_2, mll_1*mll_1,
					alp,
					chee2,
					j,k,i
				);
				
				cout<<"cheee2 "<<chee2<<" betaZ: "<<betaZ<<endl;

				cout<<"alphas: "<<alp[0]<<"\t"<<alp[1]<<"\t"<<alp[2]<<endl;
				//if (alp[0] < 0 or alp[1] < 0 or alp[2] < 0) continue;
				if (chee2 >= bestchi2) continue;
				//if (abs(mHb-mHa) > dmH) continue;
				bestchi2 = chee2;
				dmH = abs(mHb-mHa);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i] = alp[0];
				alpha[j] = alp[1];
				alpha[k] = alp[2];
				cout<<"BEST CHI2: "<<bestchi2<<" "<<chee2<<endl;
				//cout<<"alphas: "<<alpha[0]<<"\t"<<alpha[1]<<"\t"<<alpha[2]<<"\t"<<alpha[3]<<endl;
				cout<<"mH1: "<<mHa<<"\t mH2: "<<mHb<<endl;
			}*/
			//testing 3nu_good.h
			/*for(int k=j+1;k<4;k++){cout<<"INDEX: "<<i<<j<<k<<endl;
				bool k_in = std::find(ileg.begin(), ileg.end(), k) != ileg.end();
				//if(!k_in) continue;
				vec3 alp;
				bool okay;
				
				out = Fit3Nu_Option2_FullMass_YBarrier(
					Lep_ref, // size = 4
					MET_ref,
					metcov,
					i,j,k,
					i,j,k
				);
				
				cout<<"cheee2 "<<out.chi2 <<"\t"<<out.x[0]<<"\t"<<out.x[1]<<"\t"<<out.x[2]<<endl;
				if (out.x[0]*out.x[1]*out.x[2] ==0)continue;
				if (out.chi2 >= bestchi2) continue;
				bestchi2 = out.chi2;
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i]=out.x[0]/Lep_ref[i].Pt();
				alpha[j]=out.x[1]/Lep_ref[j].Pt();
				alpha[k]=out.x[2]/Lep_ref[k].Pt();

				cout<<"mH1: "<<out.mA<<"\t mH2: "<<out.mB<<endl;
			}*/
			//testing 3nu_best.h
			/*for(int k=j+1;k<4;k++){cout<<"INDEX: "<<i<<j<<k<<endl;
				bool k_in = std::find(ileg.begin(), ileg.end(), k) != ileg.end();
				//if(!k_in) continue;
				
				std::array<int,3> idxNu = {0,1,2}; // leptons 0,1,2 have neutrinos
				std::array<int,2> pairA = {0,1};  // Higgs A from leptons 0+1
				std::array<int,2> pairB = {2,3};  // Higgs B from leptons 2+3

				if (i==0 and j==1)  out = Fit3Neutrinos_GlobalLocal(
					Lep_ref, // size = 4
					MET_ref,
					metcov,
					i,j,k,
					i,j,k
				);
				else if (j==2 and k==3) 
				
				out = Fit3Neutrinos_GlobalLocal(
					Lep_ref, // size = 4
					MET_ref,
					metcov,
					i,j,k,
					i,j,k
				);
				
				cout<<"cheee2 "<<out.chi2 <<"\t"<<out.x[0]<<"\t"<<out.x[1]<<"\t"<<out.x[2]<<endl;
				if (out.x[0]*out.x[1]*out.x[2] ==0)continue;
				if (out.chi2 >= bestchi2) continue;
				bestchi2 = out.chi2;
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[i]=out.x[0]/Lep_ref[i].Pt();
				alpha[j]=out.x[1]/Lep_ref[j].Pt();
				alpha[k]=out.x[2]/Lep_ref[k].Pt();

				cout<<"mH1: "<<out.mA<<"\t mH2: "<<out.mB<<endl;
			}*/
			
			//cout<<"mH1: "<<mHa<<"\t mH2: "<<mHb<<"\t M2star: "<<M2star<<endl;
    	}
    }

	/*//3 neutrinos
	float dmH = 99999;//abs(mll_1 - mll_2) ;
    for (int i=0;i<4;i++){
    	for (int j=i+1;j<4;j++){
    		for (int k=j+1;k<4;k++){
				bool i_in = std::find(ileg.begin(), ileg.end(), i) != ileg.end();
		        bool j_in = std::find(ileg.begin(), ileg.end(), j) != ileg.end();
		        bool k_in = std::find(ileg.begin(), ileg.end(), k) != ileg.end();
		        if (ileg.size()==1 and !(i_in or j_in or k_in)) continue;
				if (ileg.size()==2 and !((i_in and j_in) or (j_in and k_in) or (k_in and i_in))) continue;
				if (ileg.size()>2 and !((i_in and j_in and k_in))) continue;
				
				int ii = -9, jj = -9, kk = -9;
				if(i==0 and j==1){ii = i; jj = j; kk = k;}//k = 2 or 3
				else if(j==2 and k==3){ii = j; jj = k; kk = i;}//i = 0 or 1
				int kk_1 = 0+1+2+3 - (ii + jj + kk);//leftover idx

				double A1 = Lxy[ii].Dot(Lxy[kk]);
				double A2 = Lxy[jj].Dot(Lxy[kk]);
				double K = Lxy[kk].Dot(Lxy[kk]);
				double C1 = (Lxy[ii].Cross(Lxy[kk])).Mag();
				double C2 = (Lxy[jj].Cross(Lxy[kk])).Mag();
				double N3 = (MET_ref.Cross(Lxy[kk])).Mag();
				double S3 = MET_ref.Dot(Lxy[kk]);
				
				double r = pow((Lep_ref[ii]+Lep_ref[jj]).M()/(Lep_ref[kk]+Lep_ref[kk_1]).M(),2);
				double A = C1*K*r;
				double B = A2*C1-C2*A1-K*N3*r;
				double C = -A2*N3+C2*(K*(1-r)+S3);
				
				double alpha1 = (-B+sqrt(B*B-4*A*C))*0.5/A;
				double alpha2 = (N3-alpha1*C1)/C2;
				double alpha3 = r*(1+alpha1)*(1+alpha2)-1;
				
				if (C1< 0.001 or C2 < 0.001) continue;
				if (alpha1 < 0 or alpha2 < 0 or alpha3 < 0) continue;
				cout<<"ALPHAS "<<alpha1<<"\t"<<alpha2<<"\t"<<alpha3<<endl;
				double mHa =  -99, mHb = -999;
				if (i==0 and j==1){
					mHa = sqrt((1+alpha1)*(1+alpha2))*(Lep_ref[ii]+Lep_ref[jj]).M();
				 	mHb = sqrt(1+alpha3)*(Lep_ref[kk]+Lep_ref[kk_1]).M();
				}
				else if (j==2 and k==3){
					mHa = sqrt(1+alpha3)*(Lep_ref[kk]+Lep_ref[kk_1]).M();
					mHb = sqrt((1+alpha1)*(1+alpha2))*(Lep_ref[ii]+Lep_ref[jj]).M();
				}
				
				cout<<"mH1: "<<mHa<<"\t mH2: "<<mHb<<endl;
				if (abs(mHa - mHb) >= dmH) continue;
				dmH = abs(mHa - mHb);
				for(int l=0; l<4; l++) alpha[l]=0;
				alpha[ii] = alpha1;
				alpha[jj] = alpha2;
				alpha[kk] = alpha3;
			}
    	}
    }*/
    
    //BoostLeptonsZ(Lep_ref, betaZ);
    TLorentzVector nu[4]; int nlegs =0;
    for (int i=0;i<4;++i) {//making neutrinos
    	//nu[i].SetPxPyPzE(alpha[i]*Lep_ref[i].Px(), alpha[i]*Lep_ref[i].Py(), alpha[i]*Lep_ref[i].Pz(), alpha[i]*Lep_ref[i].E());
    	nu[i] = alpha[i]*Lep_ref[i]; 
    	if (nu[i].E()>0) nlegs +=1;
    }
    
    double MDCH1 = (nu[0]+nu[1]+Lep_ref[0]+Lep_ref[1]).M();
	double MDCH2 = (nu[2]+nu[3]+Lep_ref[2]+Lep_ref[3]).M();

	//if(nlegs!=3) return {nlegs, MDCH1, MDCH2, bestdPhi, bestChi2, isOpp};
	if (doPrint){
		for (int i=0;i<4;++i) {cout<<"nu"<<i<<": "; print4Vec(nu[i]);}
		for (int i=0;i<4;++i) {cout<<"Lep"<<i<<": "; print4Vec(Lep_ref[i]);}
		cout<<"MET "<<MET_ref.Mag()<<"\t"<<"Sum (nu pT) "<<(nu[0]+nu[1]+nu[2]+nu[3]).Pt()<<endl;
		cout<<"Mll "<< mll_1 <<"\t"<< mll_2 <<endl;
		cout<<"DCH masses "<< MDCH1 <<"\t"<< MDCH2<<endl;
		cout<<"SOL masses "<< sol.mHA <<"\t"<< sol.mHB<<endl;
	}
    //cout<<nlegs<<endl;
    return {nlegs, MDCH1, MDCH2, bestdPhi, bestChi2, isOpp};
}

	
