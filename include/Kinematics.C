#include "TLorentzVector.h"
#include "TVector3.h"
#include <TMinuit.h>

#include "Kinematics.h"
#include "cat.h"

void print4Vec(TLorentzVector Vec){
	cout<< Vec.Pt()<<"\t"<<Vec.Eta()<<"\t"<<Vec.Phi()<<"\t"<<Vec.M()<<endl;
}

std::vector<double> SortPt(){
	double pt[] = {pt_1,pt_2,pt_3,pt_4};
	std::vector<double> sorted_pt;
	for (int j = 0; j < 4; j++){
		sorted_pt.push_back(pt[j]);
	}
	std::sort(sorted_pt.begin(), sorted_pt.end(), greater());
	return sorted_pt;
}

string applyHEMveto(string cat){
    if (cat.size() > 0 && cat[0]=='e' && eta_1 > -3.0 && eta_1 < -1.3 && phi_1 > -1.57 && phi_1 < -0.87 && pt_1 >15) return "yes";
    if (cat.size() > 1 && cat[1]=='e' && eta_2 > -3.0 && eta_2 < -1.3 && phi_2 > -1.57 && phi_2 < -0.87 && pt_2 >15) return "yes";
    if (cat.size() > 2 && cat[2]=='e' && eta_3 > -3.0 && eta_3 < -1.3 && phi_3 > -1.57 && phi_3 < -0.87 && pt_3 >15) return "yes";
    if (cat.size() > 3 && cat[3]=='e' && eta_4 > -3.0 && eta_4 < -1.3 && phi_4 > -1.57 && phi_4 < -0.87 && pt_4 >15) return "yes";
    return "no";
}

void applyTauES(string cat){
    if (cat.size() > 0 && cat[0] == 't'){ pt_1 *= TauES_1; m_1 *= TauES_1; }
    if (cat.size() > 1 && cat[1] == 't'){ pt_2 *= TauES_2; m_2 *= TauES_2; }
    if (cat.size() > 2 && cat[2] == 't'){ pt_3 *= TauES_3; m_3 *= TauES_3; }
    if (cat.size() > 3 && cat[3] == 't'){ pt_4 *= TauES_4; m_4 *= TauES_4; }
}

double ST(string cat){
	double st = 0;
	if (cat.find("e")+1==1 or cat.find("m")+1==1)
		st += pt_1;
	if (cat.find("e")+1==2 or cat.find("m")+1==2)
		st += pt_2;
	if (cat.find("e")+1==3 or cat.find("m")+1==3)
		st += pt_3;
	if (cat.find("e")+1==4 or cat.find("m")+1==4)
		st += pt_4;
	return st;
}

int remaining_idx(vector<pair<int, int>>& pair, std::string cat_name){
	int catSize = cat_name.length();
    std::set<int> used_indices;

    for (const auto& p : pair) {
        used_indices.insert(p.first);
        used_indices.insert(p.second);
    }
    for (int i = 1; i <= catSize; ++i) {
        if (used_indices.find(i) == used_indices.end()) return i;
    }
}

double dPhi(double phi1, double phi2){
	double pi = TMath::Pi();
    double DPhi = fmin(fabs(phi2 - phi1), 2.0 * pi - fabs(phi2 - phi1));
    return DPhi;
}
double deltaPhi(const TLorentzVector& v1, const TLorentzVector& v2){
	return dPhi(v1.Phi(), v2.Phi());
}

double getDR(double eta1, double phi1, double eta2, double phi2) {
    double DR = sqrt(pow(dPhi(phi1, phi2), 2) + pow(eta2 - eta1, 2));
    return DR;
}

double deltaR(const TLorentzVector& v1, const TLorentzVector& v2) {
    double eta1 = v1.Eta();
    double phi1 = v1.Phi();
    double eta2 = v2.Eta();
    double phi2 = v2.Phi();
    return getDR(eta1, phi1, eta2, phi2);
}

double calculateMT(const TLorentzVector& V, const TLorentzVector& MET) {

    double pt_total = (V+MET).Pt();
    return std::sqrt(std::pow(V.Et() + MET.Et(), 2) - pt_total * pt_total);
}

double calculateMTtot(const TLorentzVector& l1, const TLorentzVector& l2) {
	TLorentzVector LLpair = l1+l2;

    double px_ll = LLpair.Px();
    double py_ll = LLpair.Py();

    double px_met = met * std::cos(metphi);
    double py_met = met * std::sin(metphi);

    double px_total = px_ll + px_met;
    double py_total = py_ll + py_met;

    double Et_ll = LLpair.Et();
    double pt_total = std::sqrt(px_total * px_total + py_total * py_total);

    double mt_tot = std::sqrt(std::pow(Et_ll + met, 2) - pt_total * pt_total);

    return mt_tot;
}

std::pair<double, double> ComputeDCHMasses(
    const TLorentzVector& lep1,
    const TLorentzVector& lep2,
    const TLorentzVector& lep3,
    const TLorentzVector& lep4,
    const TLorentzVector& MET)
{

    TLorentzVector vis = lep1 + lep2 + lep3 + lep4;

    double beta_z = -vis.Pz() / vis.E();
    TVector3 boostVec(0, 0, beta_z);

    TLorentzVector l1 = lep1;
    TLorentzVector l2 = lep2;
    TLorentzVector l3 = lep3;
    TLorentzVector l4 = lep4;

    l1.Boost(boostVec);
    l2.Boost(boostVec);
    l3.Boost(boostVec);
    l4.Boost(boostVec);

    TLorentzVector met4 = MET;
	met4.Boost(boostVec);

    TLorentzVector leg1_vis = l1 + l2;
    TLorentzVector leg2_vis = l3 + l4;

    TVector3 pT_leg1_miss = 0.5 * (met4.Vect() + (leg2_vis.Vect() - leg1_vis.Vect()));
    pT_leg1_miss.SetZ(0);

    TVector3 pT_leg2_miss = met4.Vect() - pT_leg1_miss;
    pT_leg2_miss.SetZ(0);

    double a_guess = 0;
    auto mass_diff = [&](double a){
        TLorentzVector miss1(pT_leg1_miss.X(), pT_leg1_miss.Y(), a, sqrt(pT_leg1_miss.Mag2() + a*a));
        TLorentzVector miss2(pT_leg2_miss.X(), pT_leg2_miss.Y(), -a, sqrt(pT_leg2_miss.Mag2() + a*a));
        double m1 = (leg1_vis + miss1).M();
        double m2 = (leg2_vis + miss2).M();
        return m1 - m2;
    };

    double lo = -2000, hi = 2000;
    for (int i=0; i<100; i++) {
        double mid = 0.5*(lo+hi);
        if (mass_diff(lo)*mass_diff(mid) <= 0) hi = mid;
        else lo = mid;
    }
    double a_sol = 0.5*(lo+hi);

    TLorentzVector miss1(pT_leg1_miss.X(), pT_leg1_miss.Y(), a_sol, sqrt(pT_leg1_miss.Mag2() + a_sol*a_sol));
    TLorentzVector miss2(pT_leg2_miss.X(), pT_leg2_miss.Y(), -a_sol, sqrt(pT_leg2_miss.Mag2() + a_sol*a_sol));

    TLorentzVector dch1 = leg1_vis + miss1;
    TLorentzVector dch2 = leg2_vis + miss2;

    return std::make_pair(dch1.M(), dch2.M());
}

struct InputData {
    TLorentzVector l1, l2, l3, l4;
    TVector2 MET;
};

static TLorentzVector VA, VB;
static TVector2 METvec;
static double minMass = 1e9;
static double best_pzA = 0, best_pzB = 0;

void fcn(int &npar, double *gin, double &f, double *par, int iflag) {

    double pzA = par[0];
    double pzB = par[1];
    double pxA = par[2];
    double pyA = par[3];

    double pxB = METvec.X() - pxA;
    double pyB = METvec.Y() - pyA;

    double pA_E = sqrt(pxA*pxA + pyA*pyA + pzA*pzA);
    double pB_E = sqrt(pxB*pxB + pyB*pyB + pzB*pzB);

    TLorentzVector IA(pxA, pyA, pzA, pA_E);
    TLorentzVector IB(pxB, pyB, pzB, pB_E);

    double MA = (VA + IA).M();
    double MB = (VB + IB).M();

	f = std::max(MA, MB);

    if (f < minMass) {
        minMass = f;
        best_pzA = pzA;
        best_pzB = pzB;
    }
}

std::pair<double, double>ReconstructMass(const InputData &in) {
    VA = in.l1 + in.l2;
    VB = in.l3 + in.l4;
    METvec = in.MET;

    TMinuit minuit(4);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(fcn);

    double step[4] = {1., 1., 1., 1.};
    double start[4] = {0., 0., 0., METvec.X()/2.};
    double bound[4] = {0.};

    for (int i = 0; i < 4; i++) {
        minuit.DefineParameter(i, Form("par%d", i), start[i], step[i], 0, 0);
    }

    minuit.Migrad();

    double pxA = start[3];
    double pyA = start[3];
    double pzA = best_pzA;

    double pxB = METvec.X() - pxA;
    double pyB = METvec.Y() - pyA;
    double pzB = best_pzB;

    TLorentzVector IA(pxA, pyA, pzA, sqrt(pxA*pxA + pyA*pyA + pzA*pzA));
    TLorentzVector IB(pxB, pyB, pzB, sqrt(pxB*pxB + pyB*pyB + pzB*pzB));

    return std::make_pair((VA + IA).M(), (VB + IB).M());
}

bool isDuplicate(const Lepton &lepton1, const Lepton &lepton2){
	 if ( getDR(lepton1.eta, lepton1.phi, lepton2.eta, lepton2.phi) <= 0.4)
	 	return true;
	 else return false;
}

TLorentzVector LepV(int n);

string pairFunc(int m, int n, string cat, double Zwindow){
	TLorentzVector lep1, lep2;
	int c1, c2;
	if (m==1){
		lep1 = LepV(1);
		c1 = q_1;
	}
	else if (m==2){
		lep1 = LepV(2);
		c1 = q_2;
	}
	else if (m==3){
		lep1 = LepV(3);
		c1 = q_3;
	}
	else if (m==4){
		lep1 = LepV(4);
		c1 = q_4;
	}
	if (n==1){
		lep2 = LepV(1);
		c2 = q_1;
	}
	else if (n==2){
		lep2 = LepV(2);
		c2 = q_2;
	}
	else if (n==3){
		lep2 = LepV(3);
		c2 = q_3;
	}
	else if (n==4){
		lep2 = LepV(4);
		c2 = q_4;
	}
	if (c1 == c2 ) return "SSpair";
	else if (c1 == -c2){
		if (cat[m-1] == cat[n-1]){
			if (cat[m-1] == 'e' or cat[m-1] == 'm' ){
				return (abs((lep1+lep2).M()-91.2) < Zwindow) ? "Zwindow" : "Zv";
			}
			else if (cat[m-1] == 't') return "ZttPair" ;
		}
		else return "OSDFpair";
	}
	else "messed up";
}

TLorentzVector LepV(int n){
	TLorentzVector lepV;
	if (n==1)
		lepV.SetPtEtaPhiM(pt_1, eta_1, phi_1, 0);
	else if (n==2)
		lepV.SetPtEtaPhiM(pt_2, eta_2, phi_2, 0);
	else if (n==3)
		lepV.SetPtEtaPhiM(pt_3, eta_3, phi_3, 0);
	else if (n==4)
		lepV.SetPtEtaPhiM(pt_4, eta_4, phi_4, 0);
	else lepV.SetPtEtaPhiM(0,0,0,0);
	return lepV;
}

void processPairs( const char* cat_name, vector<pair<int, int>>& Z_pair, vector<pair<int, int>>& Zv_pair, vector<pair<int, int>>& Ztt_pair, vector<pair<int, int>>& SS_pair, vector<pair<int,int>>& OSDF_pair) {
	int len = strlen(cat_name);

    vector<bool> used(len+1,false);

    struct Zcand {
        int i,j;
        double diff;
    };
    vector<Zcand> Zcands;

    for (int m=1; m<=len; ++m){
        for (int n=m+1; n<=len; ++n){

            string pair_name = pairFunc(m,n,cat_name,10);

            if(pair_name=="Zwindow"){
                double mass = (LepV(m)+LepV(n)).M();
                double diff = fabs(mass - 91.2);
                Zcands.push_back({m,n,diff});
            }
            else if(pair_name=="ZttPair") Ztt_pair.push_back({m,n});
            else if(pair_name=="OSDFpair") OSDF_pair.push_back({m,n});
            else if(pair_name=="SSpair") SS_pair.push_back({m,n});
        }
    }

    sort(Zcands.begin(), Zcands.end(),
         [](const Zcand& a,const Zcand& b){ return a.diff < b.diff; });

    for(auto &zcand : Zcands){
        if(!used[zcand.i] && !used[zcand.j]){
            Z_pair.push_back({zcand.i,zcand.j});
            used[zcand.i] = true;
            used[zcand.j] = true;
        }
    }

    for (int m=1; m<=len; ++m){
    	if(used[m]) continue;
        for (int n=m+1; n<=len; ++n){
            if(used[n]) continue;
            string pair_name = pairFunc(m,n,cat_name,10);
            if(pair_name=="Zv"){
                Zv_pair.push_back({m,n});
                used[m] = true;
                used[n] = true;
            }
        }
    }

}

vector<pair<int,int>> removeOverlap(const vector<pair<int,int>>& pairs, int nLep)
{
    vector<pair<int,int>> clean;
    vector<bool> used(nLep + 1,false);

    for(const auto& p : pairs){
        if(!used[p.first] && !used[p.second]){
            clean.push_back(p);
            used[p.first]  = true;
            used[p.second] = true;
        }
    }
    return clean;
}

std::string classifyTauRegion(std::string cat_name, double LT, vector<pair<int, int>> OS_pair){
	int Ntau = cat_name.size()-cat_lepCount(cat_name,'e','m');
	double Zwindow = 10, Zmass = 91.2;
	double mPair = 0;
	int numZ =0;
	for (int i = 0; i<OS_pair.size(); i++){
		int ia = OS_pair[i].first, ib = OS_pair[i].second;
		if (cat_name[ia-1]=='t' or cat_name[ia-1] != cat_name[ib-1])continue;
		mPair = (LepV(ia)+LepV(ib)).M();
		if (abs(mPair- Zmass) <= Zwindow) numZ += 1;

	}
	if (cat_name.size()==2){
		if (Ntau==0) return (numZ==1) ? "DYCR_0tau" : "DYveto_0tau";
		else if (Ntau==1) return (numZ==1) ? "DYCR_1tau" : "DYveto_1tau";
	}
	if (numZ==0 and cat_name.size()==4){
		if (Ntau==0) return (LT < 400) ? "VR_0tau" : "SR_0tau";
		else if (Ntau==1) return (LT < 400) ? "VR_1tau" : "SR_1tau";
		else if (Ntau==2) return (LT < 400) ? "VR_2tau" : "SR_2tau";
		else if (Ntau==3) return (LT < 100) ? "VR_3tau" : "SR_3tau";
	}
	else if (numZ==0 and cat_name.size()==3){
		if (Ntau==0) return (LT < 300) ? "VR_3lep0tau" : "SR_3lep0tau";
		else if (Ntau==1) return (LT < 200) ? "VR_3lep1tau" : "SR_3lep1tau";
		else if (Ntau==2) return (LT < 100) ? "VR_3lep2tau" : "SR_3lep2tau";
	}
	else if (numZ > 0 and cat_name.size()==3){
		if (Ntau==0) return "CR_3lep0tau";
		else if (Ntau==1) return "CR_3lep1tau";
		else if (Ntau==2) return "CR_3lep2tau";
	}
	else if (numZ > 0 and cat_name.size()==4){
		if (Ntau==0) return "CR_0tau";
		else if (Ntau==1) return "CR_1tau";
		else if (Ntau==2) return "CR_2tau";
		else if (Ntau==3) return "CR_3tau";
	}
	else {cout<<"CRAAAPPP "<<cat_name<<"\t"<<numZ<<endl; return "crap";}
}

std::string classifyLepRegion(std::string cat_name, vector<pair<int, int>> pairVec){
	std::string lepRegion = "";
	int cat_size = cat_name.size();

	int nPairs = pairVec.size();
	if(nPairs == 0) return cat_name;
	if(cat_size == 2) return cat_name;

	vector<bool> used(cat_name.size(), false);

    for (auto &p : pairVec) {

        lepRegion += cat_name[p.first-1];
		lepRegion += cat_name[p.second-1];
        used[p.first-1] = true;
        used[p.second-1] = true;
    }

    string extra = "";
    for (int i = 0; i < cat_name.size(); i++) {
        if (used[i]) continue;
        extra += cat_name[i];
    }

    if (extra == "me") extra = "em";
    else if (extra == "te") extra = "et";
    else if (extra == "tm") extra = "mt";
    lepRegion += extra;

    return lepRegion;
}
