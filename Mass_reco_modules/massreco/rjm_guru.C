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

//-----------if one and two alphas, solve for empherical solutions---------------------
double solve_one_alpha(const TVector3 &uLxy, const TVector3 &MET) {
	double alpha = MET.Mag();//uLxy.Dot(MET);
    //std::cout<<"ONE "<<alpha <<std::endl; //since u is unit, u.MET = alpha
    return alpha;
}

std::tuple<double, double> solve_two_alphas(const TVector3 &uL1xy, const TVector3 &uL2xy, const TVector3 &MET) {
    double alpha1 = MET.Cross(uL2xy)*uL1xy.Cross(uL2xy)/uL1xy.Cross(uL2xy).Mag2();
    double alpha2 = MET.Cross(uL1xy)*uL2xy.Cross(uL1xy)/uL2xy.Cross(uL1xy).Mag2();

	//std::cout << "TWO "<<"\t"<<alpha1<<"\t"<<alpha2<<std::endl; 
	//eliminate negative alphas and fall back to one alpha.
	if (alpha1 < 0 and alpha2 < 0) return {0,0};
	else if (alpha1 < 0){
    		alpha1 = 0;
    		alpha2 = solve_one_alpha(uL2xy ,MET);
    }
    else if (alpha2 < 0){
    		alpha1 = solve_one_alpha(uL1xy ,MET);
    		alpha2 = 0;
    }
    return {alpha1, alpha2};
}

std::tuple<double, double, double> solve_three_alphas(const TVector3 &L1, const TVector3 &L2, const TVector3 &L3, const TVector3 &MET) {
	double c1 = -L1.Z()/L3.Z();
	double c2 = -L2.Z()/L3.Z();
	 
	// Build transverse vectors (TVector3)
    TVector3 p1T(L1.X(), L1.Y(), 0);
    TVector3 p2T(L2.X(), L2.Y(), 0);
    TVector3 p3T(L3.X(), L3.Y(), 0);
    
    TVector3 q1 = p1T.Unit() + c1 * p3T.Unit();
    TVector3 q2 = p2T.Unit() + c2 * p3T.Unit();
	 
	auto cross2 = [](const TVector3 &a, const TVector3 &b){
        return a.X()*b.Y() - a.Y()*b.X();
    };
    
    // denominator
    double denom = cross2(q1, q2);
    if (std::abs(denom) < 1e-9) {
        // q1 and q2 are (nearly) collinear -> no unique solution
        return {0,0,0};
    }
    // Now solve for alpha1 and alpha2
    // alpha1 = (met x q2) / (q1 x q2)
    // alpha2 = (met x q1) / (q2 x q1) = - (met x q1)/(q1 x q2)
    TVector3 metT(MET.X(), MET.Y(), 0);
    double met_x_q2 = cross2(metT, q2);
    double met_x_q1 = cross2(metT, q1);

    double alpha1 = met_x_q2 / denom;
    double alpha2 = - met_x_q1 / denom; // because q2 x q1 = - (q1 x q2)
    double alpha3 = c1 * alpha1 + c2 * alpha2;
    //std::cout<< "THREE "<<"\t"<<alpha1<<"\t"<<alpha2<<"\t"<<alpha3<<std::endl; 
    
	if (alpha1 < 0 and alpha2 < 0){
 		alpha1 = 0;
 		alpha2 = 0;
   		alpha3 = solve_one_alpha(p3T.Unit(), MET);
    }
    else if (alpha1 < 0 and alpha3 < 0){
 		alpha1 = 0;
 		alpha3 = 0;
   		alpha2 = solve_one_alpha(p2T.Unit(), MET);
    }
    else if (alpha2 < 0 and alpha3 < 0){
 		alpha2 = 0;
 		alpha3 = 0;
   		alpha1 = solve_one_alpha(p1T.Unit(), MET);
    }
    else if (alpha1 < 0){
 		alpha1 = 0;
   		alpha2 = std::get<0>(solve_two_alphas(p2T.Unit(), p3T.Unit(),MET));
   		alpha3 = std::get<1>(solve_two_alphas(p2T.Unit(), p3T.Unit(),MET));
    }
    else if (alpha2 < 0){
 		alpha2 = 0;
   		alpha1 = std::get<0>(solve_two_alphas(p1T.Unit(), p3T.Unit(),MET));
   		alpha3 = std::get<1>(solve_two_alphas(p1T.Unit(), p3T.Unit(),MET));
    }
    else if (alpha3 < 0){
 		alpha3 = 0;
   		alpha1 = std::get<0>(solve_two_alphas(p1T.Unit(), p2T.Unit(),MET));
   		alpha2 = std::get<1>(solve_two_alphas(p1T.Unit(), p2T.Unit(),MET));
    }
    return {alpha1, alpha2, alpha3};
}
//------------------------------------------------------------------------------
struct FitContext {
    std::array<TVector3,4> Lxy;
    std::array<TLorentzVector,4> Lep_ref;
    TVector2 MET;
    std::vector<int> ileg;
    double w_MET = 1.0;
    double w_mass = 1.0;
};

FitContext g_ctx;

// --------------------- Chi^2 function for Minuit ---------------------
void chi2Fcn(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t iflag) {
    const auto &Lxy = g_ctx.Lxy;
    const auto &Lep_ref = g_ctx.Lep_ref;
    const auto &ileg = g_ctx.ileg;

    double sumx = 0, sumy = 0;
    TLorentzVector nu[4];
    for (int i=0;i<4;++i) nu[i].SetPxPyPzE(0,0,0,0);

    for (size_t i=0; i<ileg.size(); ++i) {
        int idx = ileg[i];
        TVector3 v = Lxy[idx].Unit();
        sumx += par[i] * v.X();
        sumy += par[i] * v.Y();

        TLorentzVector tmp;
        tmp.SetPxPyPzE(par[i]*v.X(), par[i]*v.Y(), 0, par[i]); // massless neutrino
        nu[idx] = tmp;
    }

    // MET residual term
    double dx = sumx - g_ctx.MET.X();
    double dy = sumy - g_ctx.MET.Y();
    double chi2_MET = dx*dx + dy*dy;

    // DCH masses
    double MDCH1 = (Lep_ref[0]+Lep_ref[1]+nu[0]+nu[1]).M();
    double MDCH2 = (Lep_ref[2]+Lep_ref[3]+nu[2]+nu[3]).M();
    double chi2_mass = (MDCH1 - MDCH2)*(MDCH1 - MDCH2);

    fval = g_ctx.w_MET * chi2_MET + g_ctx.w_mass * chi2_mass;
}

std::tuple<int, double, double> get_legs(const std::array<TLorentzVector,4> &Lep_lab, TLorentzVector &MET_lab){//this gives valid lepton legs with neutrinos
	//Boosting to have Pz_tot(vis) = 0
	double betaZ = 0; 
	std::array<TLorentzVector,4> Lep_ref; 
	TVector3 MET_ref; TVector3 Lxy[4];
	double alpha[4] = {0,0,0,0}; 
	for (int i=0;i<4;++i) betaZ += Lep_lab[i].Pz();  
	for (int i=0;i<4;++i) {
		Lep_ref[i] = Lep_lab[i]; 
		//Lep_ref[i].Boost(0,0,-betaZ);
		MET_ref = MET_lab.Vect(); //no need to boost met
		TVector3 v(Lep_ref[i].Vect().X(), Lep_ref[i].Vect().Y(), 0);
		Lxy[i] = v;
	}
	std::vector<int> ileg;
	for (int i=0;i<4;++i) {
        if (Lxy[i].Mag() <= 0){ continue;}
        double proj = Lxy[i].Dot(MET_ref);
        if (proj > 0.001) ileg.push_back(i);
    }
    if ((int)ileg.size() ==1) {
    	alpha[ileg[0]] = solve_one_alpha(Lxy[ileg[0]].Unit(),MET_ref);
    }
    else if ((int)ileg.size() ==2){
    	auto [alpha1, alpha2] = solve_two_alphas(Lxy[ileg[0]].Unit(),Lxy[ileg[1]].Unit(),MET_ref);
    	alpha[ileg[0]] = alpha1;
    	alpha[ileg[1]] = alpha2;
    }
    else if ((int)ileg.size() ==3){
    	auto [alpha1, alpha2, alpha3] = solve_three_alphas(Lep_ref[ileg[0]].Vect(), Lep_ref[ileg[1]].Vect(), Lep_ref[ileg[2]].Vect(), MET_ref);
    	alpha[ileg[0]] = alpha1;
    	alpha[ileg[1]] = alpha2;
    	alpha[ileg[2]] = alpha3;	
    }
    //else if (ileg.size() >= 3) {
    else if (ileg.size() > 3) {
        // setup fit context
        for (int i=0;i<4;++i) g_ctx.Lxy[i] = Lxy[i];
        g_ctx.Lep_ref = Lep_ref;
        g_ctx.MET = TVector2(MET_ref.X(), MET_ref.Y());
        g_ctx.ileg = ileg;
        g_ctx.w_MET = 1.0;
        g_ctx.w_mass = 0;  

        TMinuit minuit(ileg.size());
        minuit.SetPrintLevel(-1);
        minuit.SetFCN(chi2Fcn);

        for (size_t i=0; i<ileg.size(); ++i) {
            double start = 1.0, step = 0.1;
            minuit.DefineParameter(i, Form("alpha%d", i), start, step, 0.0, 1e3);
        }

        minuit.Migrad();

        for (size_t i=0; i<ileg.size(); ++i) {
            double val, err;
            minuit.GetParameter(i, val, err);
            alpha[ileg[i]] = std::max(0.0, val);
        }
    }
    TLorentzVector nu[4]; double nlegs =0;
    for (int i=0;i<ileg.size();++i) {//making neutrinos
    	int idx = ileg[i];
    	TLorentzVector tmp;
    	//tmp.SetVectM(Lep_ref[idx].Vect(),0);
    	tmp = Lep_ref[idx];
    	//tmp.SetVectM(Lxy[idx],0.0);
    	nu[idx] = (alpha[idx]/Lxy[idx].Mag())*tmp;
    	if (nu[idx].E()>0) nlegs +=1;
    }
    for (int i=0;i<4;++i) {cout<<"nu"<<i<<": "; print4Vec(nu[i]);}
    double mll_1 = (Lep_ref[0]+Lep_ref[1]).M();
    double mll_2 = (Lep_ref[2]+Lep_ref[3]).M();
    
    double MDCH1 = (nu[0]+nu[1]+Lep_ref[0]+Lep_ref[1]).M();
    double MDCH2 = (nu[2]+nu[3]+Lep_ref[2]+Lep_ref[3]).M();
    cout<<"MET "<<MET_ref.Mag()<<"\t"<<"Sum (nu pT) "<<(nu[0]+nu[1]+nu[2]+nu[3]).Pt()<<endl;
    cout<<"Mll "<< mll_1 <<"\t"<< mll_2 <<endl;
    cout<<"DCH masses "<< (nu[0]+nu[1]+Lep_ref[0]+Lep_ref[1]).M() <<"\t"<< (nu[2]+nu[3]+Lep_ref[2]+Lep_ref[3]).M()<<endl;
    cout<<nlegs<<endl;
    return {nlegs, MDCH1, MDCH2};
}

std::tuple<int, double, double> get_legs_tau(std::string cat, const std::array<TLorentzVector,4> &Lep_lab, TLorentzVector &MET_lab){//this gives valid lepton legs with neutrinos
	//Boosting to have Pz_tot(vis) = 0
	double betaZ = 0; 
	std::array<TLorentzVector,4> Lep_ref; 
	TVector3 MET_ref; TVector3 Lxy[4];
	double alpha[4] = {0,0,0,0}; 
	//for (int i=0;i<4;++i) betaZ += Lep_lab[i].Pz();  
	for (int i=0;i<4;++i) {
		Lep_ref[i] = Lep_lab[i]; 
		//Lep_ref[i].Boost(0,0,-betaZ);
		MET_ref = MET_lab.Vect(); //no need to boost met
		TVector3 v(Lep_ref[i].Vect().X(), Lep_ref[i].Vect().Y(), 0);
		Lxy[i] = v;
	}
	std::vector<int> ileg;
	for (int i=0;i<4;++i) {
        if (cat[i]=='t') ileg.push_back(i);
    }
    if (int(ileg.size())==1) {
   		alpha[ileg[0]] = solve_one_alpha(Lxy[ileg[0]].Unit(),MET_ref);
    }
    else if (int(ileg.size())==2) {
    	auto [alpha1, alpha2] = solve_two_alphas(Lxy[ileg[0]].Unit(),Lxy[ileg[1]].Unit(),MET_ref);
    	alpha[ileg[0]] = alpha1;
    	alpha[ileg[1]] = alpha2;
    }
    else if ((int)ileg.size() ==3){
    	auto [alpha1, alpha2, alpha3] = solve_three_alphas(Lep_ref[ileg[0]].Vect(), Lep_ref[ileg[1]].Vect(), Lep_ref[ileg[2]].Vect(), MET_ref);
    	alpha[ileg[0]] = alpha1;
    	alpha[ileg[1]] = alpha2;
    	alpha[ileg[2]] = alpha3;	
    }
    
    TLorentzVector nu[4]; int nlegs =0;
    for (int i=0;i<4;++i) {//making neutrinos
    	TLorentzVector tmp;
    	//tmp.SetVectM(Lep_ref[i].Vect(),0);
    	tmp = Lep_ref[i];
    	//tmp.SetVectM(Lxy[i],0.0);
    	nu[i] = (alpha[i]/Lxy[i].Mag())*tmp;
    	if (nu[i].E()>0) nlegs +=1;
    }
    //for (int i=0;i<4;++i) {cout<<"nu"<<i<<": "; print4Vec(nu[i]);}
    double mll_1 = (Lep_ref[0]+Lep_ref[1]).M();
    double mll_2 = (Lep_ref[2]+Lep_ref[3]).M();
    double MDCH1 = (nu[0]+nu[1]+Lep_ref[0]+Lep_ref[1]).M();
    double MDCH2 = (nu[2]+nu[3]+Lep_ref[2]+Lep_ref[3]).M();
    /*cout<<"MET "<<MET_ref.Mag()<<"\t"<<"Sum (nu pT) "<<(nu[0]+nu[1]+nu[2]+nu[3]).Pt()<<endl;
    cout<<"Mll "<< mll_1 <<"\t"<< mll_2 <<endl;
    cout<<"DCH masses "<< (nu[0]+nu[1]+Lep_ref[0]+Lep_ref[1]).M() <<"\t"<< (nu[2]+nu[3]+Lep_ref[2]+Lep_ref[3]).M()<<endl;
    cout<<nlegs<<endl;*/
    return {nlegs, MDCH1, MDCH2};
}

	
