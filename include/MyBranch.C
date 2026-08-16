#ifndef MYBRANCH_C
#define MYBRANCH_C

#include <TTree.h>
#include <vector>

using std::vector;

Int_t nElectron, nMuon, nTau;
Int_t evt;
ULong64_t run;

Int_t cat, gen_cat;
Int_t nPV;
Int_t nPVGood;
Int_t genPartFlav_1, genPartFlav_2, genPartFlav_3, genPartFlav_4;

Double_t brWeight;
Double_t weight;
Double_t weightPU;
Double_t weightPUtrue;
Double_t weightPUtruejson;
Double_t weightPUtruejson_up;
Double_t weightPUtruejson_down;
Double_t Generator_weight;

Double_t L1PreFiringWeight_Nom;
Double_t L1PreFiringWeight_Up;
Double_t L1PreFiringWeight_Down;

Double_t met, metphi;
Double_t metNoCor, metphiNoCor;
Double_t metNoTauES, metphiNoTauES;
Double_t njets;

Double_t pt_1, pt_2, pt_3, pt_4;
Double_t eta_1, eta_2, eta_3, eta_4;
Double_t phi_1, phi_2, phi_3, phi_4;
Double_t pt_1_tr, pt_2_tr, pt_3_tr, pt_4_tr;
Double_t eta_1_tr, eta_2_tr, eta_3_tr, eta_4_tr;
Double_t phi_1_tr, phi_2_tr, phi_3_tr, phi_4_tr;
Double_t m_1, m_2, m_3, m_4;
Double_t q_1, q_2, q_3, q_4;
Double_t iso_1, iso_2, iso_3, iso_4;
Double_t d0_1, d0_2, d0_3, d0_4;
Double_t dZ_1, dZ_2, dZ_3, dZ_4;

Double_t IDSF_1, IDSF_2, IDSF_3, IDSF_4;
Double_t ISOSF_1, ISOSF_2, ISOSF_3, ISOSF_4;
Double_t TrigSF_1, TrigSF_2, TrigSF_3, TrigSF_4;
Double_t isTrig_1, isTrig_2, isDoubleTrig;

Double_t TauVsEleIDSF_1, TauVsEleIDSF_2, TauVsEleIDSF_3, TauVsEleIDSF_4;
Double_t TauVsMuIDSF_1,  TauVsMuIDSF_2,  TauVsMuIDSF_3,  TauVsMuIDSF_4;
Double_t TauVsJetIDSF_1, TauVsJetIDSF_2, TauVsJetIDSF_3, TauVsJetIDSF_4;

Double_t TauES_1, TauES_2, TauES_3, TauES_4;

vector<double> *LHEScaleWeights = nullptr;

vector<double> *lpt     = nullptr;
vector<double> *leta    = nullptr;
vector<double> *lphi    = nullptr;
vector<double> *lmass   = nullptr;
vector<int>    *lq      = nullptr;
vector<int>    *lflavor = nullptr;
vector<double> *liso    = nullptr;
vector<double> *ld0     = nullptr;
vector<double> *ldZ     = nullptr;

vector<double> *EleID      = nullptr;
vector<double> *EleID_WP90 = nullptr;
vector<double> *EleID_WPL  = nullptr;

vector<int> *MuID      = nullptr;
vector<int> *TauIDe    = nullptr;
vector<int> *TauIDm    = nullptr;
vector<int> *TauIDj    = nullptr;
vector<int> *decayMode = nullptr;
vector<int> *gen_match = nullptr;

void MyBranch(TTree* tree)
{
    if (!tree) return;

    tree->SetBranchAddress("nElectron", &nElectron);
    tree->SetBranchAddress("nMuon",     &nMuon);
    tree->SetBranchAddress("nTau",      &nTau);
    tree->SetBranchAddress("evt",       &evt);
    tree->SetBranchAddress("run",       &run);

    tree->SetBranchAddress("cat",     &cat);
    tree->SetBranchAddress("gen_cat", &gen_cat);

    tree->SetBranchAddress("genPartFlav_1", &genPartFlav_1);
    tree->SetBranchAddress("genPartFlav_2", &genPartFlav_2);
    tree->SetBranchAddress("genPartFlav_3", &genPartFlav_3);
    tree->SetBranchAddress("genPartFlav_4", &genPartFlav_4);

    tree->SetBranchAddress("brWeight",                  &brWeight);
    tree->SetBranchAddress("weight",                    &weight);
    tree->SetBranchAddress("weightPU",                  &weightPU);
    tree->SetBranchAddress("weightPUtrue",              &weightPUtrue);
    tree->SetBranchAddress("weightPUtruejson",          &weightPUtruejson);
    tree->SetBranchAddress("weightPUtruejson_up",       &weightPUtruejson_up);
    tree->SetBranchAddress("weightPUtruejson_down",     &weightPUtruejson_down);
    tree->SetBranchAddress("Generator_weight",          &Generator_weight);

    tree->SetBranchAddress("L1PreFiringWeight_Nom",  &L1PreFiringWeight_Nom);
    tree->SetBranchAddress("L1PreFiringWeight_Up",   &L1PreFiringWeight_Up);
    tree->SetBranchAddress("L1PreFiringWeight_Down", &L1PreFiringWeight_Down);

    tree->SetBranchAddress("met",          &met);
    tree->SetBranchAddress("nPV",       &nPV);
    tree->SetBranchAddress("nPVGood",   &nPVGood);
    tree->SetBranchAddress("metphi",       &metphi);
    tree->SetBranchAddress("metNoCor",     &metNoCor);
    tree->SetBranchAddress("metphiNoCor",  &metphiNoCor);
    tree->SetBranchAddress("metNoTauES",   &metNoTauES);
    tree->SetBranchAddress("metphiNoTauES",&metphiNoTauES);
    tree->SetBranchAddress("njets",        &njets);

    tree->SetBranchAddress("pt_1",  &pt_1);
    tree->SetBranchAddress("pt_2",  &pt_2);
    tree->SetBranchAddress("pt_3",  &pt_3);
    tree->SetBranchAddress("pt_4",  &pt_4);

    tree->SetBranchAddress("eta_1", &eta_1);
    tree->SetBranchAddress("eta_2", &eta_2);
    tree->SetBranchAddress("eta_3", &eta_3);
    tree->SetBranchAddress("eta_4", &eta_4);

    tree->SetBranchAddress("phi_1", &phi_1);
    tree->SetBranchAddress("phi_2", &phi_2);
    tree->SetBranchAddress("phi_3", &phi_3);
    tree->SetBranchAddress("phi_4", &phi_4);

    tree->SetBranchAddress("pt_1_tr",  &pt_1_tr);
    tree->SetBranchAddress("pt_2_tr",  &pt_2_tr);
    tree->SetBranchAddress("pt_3_tr",  &pt_3_tr);
    tree->SetBranchAddress("pt_4_tr",  &pt_4_tr);

    tree->SetBranchAddress("eta_1_tr", &eta_1_tr);
    tree->SetBranchAddress("eta_2_tr", &eta_2_tr);
    tree->SetBranchAddress("eta_3_tr", &eta_3_tr);
    tree->SetBranchAddress("eta_4_tr", &eta_4_tr);

    tree->SetBranchAddress("phi_1_tr", &phi_1_tr);
    tree->SetBranchAddress("phi_2_tr", &phi_2_tr);
    tree->SetBranchAddress("phi_3_tr", &phi_3_tr);
    tree->SetBranchAddress("phi_4_tr", &phi_4_tr);

    tree->SetBranchAddress("m_1",   &m_1);
    tree->SetBranchAddress("m_2",   &m_2);
    tree->SetBranchAddress("m_3",   &m_3);
    tree->SetBranchAddress("m_4",   &m_4);

    tree->SetBranchAddress("q_1",   &q_1);
    tree->SetBranchAddress("q_2",   &q_2);
    tree->SetBranchAddress("q_3",   &q_3);
    tree->SetBranchAddress("q_4",   &q_4);

    tree->SetBranchAddress("iso_1", &iso_1);
    tree->SetBranchAddress("iso_2", &iso_2);
    tree->SetBranchAddress("iso_3", &iso_3);
    tree->SetBranchAddress("iso_4", &iso_4);

    tree->SetBranchAddress("d0_1",  &d0_1);
    tree->SetBranchAddress("d0_2",  &d0_2);
    tree->SetBranchAddress("d0_3",  &d0_3);
    tree->SetBranchAddress("d0_4",  &d0_4);

    tree->SetBranchAddress("dZ_1",  &dZ_1);
    tree->SetBranchAddress("dZ_2",  &dZ_2);
    tree->SetBranchAddress("dZ_3",  &dZ_3);
    tree->SetBranchAddress("dZ_4",  &dZ_4);

    tree->SetBranchAddress("IDSF_1", &IDSF_1);
    tree->SetBranchAddress("IDSF_2", &IDSF_2);
    tree->SetBranchAddress("IDSF_3", &IDSF_3);
    tree->SetBranchAddress("IDSF_4", &IDSF_4);

    tree->SetBranchAddress("ISOSF_1", &ISOSF_1);
    tree->SetBranchAddress("ISOSF_2", &ISOSF_2);
    tree->SetBranchAddress("ISOSF_3", &ISOSF_3);
    tree->SetBranchAddress("ISOSF_4", &ISOSF_4);

    tree->SetBranchAddress("TrigSF_1", &TrigSF_1);
    tree->SetBranchAddress("TrigSF_2", &TrigSF_2);
    tree->SetBranchAddress("TrigSF_3", &TrigSF_3);
    tree->SetBranchAddress("TrigSF_4", &TrigSF_4);

    tree->SetBranchAddress("isTrig_1",      &isTrig_1);
    tree->SetBranchAddress("isTrig_2",      &isTrig_2);
    tree->SetBranchAddress("isDoubleTrig",  &isDoubleTrig);

    tree->SetBranchAddress("TauVsEleIDSF_1", &TauVsEleIDSF_1);
    tree->SetBranchAddress("TauVsEleIDSF_2", &TauVsEleIDSF_2);
    tree->SetBranchAddress("TauVsEleIDSF_3", &TauVsEleIDSF_3);
    tree->SetBranchAddress("TauVsEleIDSF_4", &TauVsEleIDSF_4);

    tree->SetBranchAddress("TauVsMuIDSF_1", &TauVsMuIDSF_1);
    tree->SetBranchAddress("TauVsMuIDSF_2", &TauVsMuIDSF_2);
    tree->SetBranchAddress("TauVsMuIDSF_3", &TauVsMuIDSF_3);
    tree->SetBranchAddress("TauVsMuIDSF_4", &TauVsMuIDSF_4);

    tree->SetBranchAddress("TauVsJetIDSF_1", &TauVsJetIDSF_1);
    tree->SetBranchAddress("TauVsJetIDSF_2", &TauVsJetIDSF_2);
    tree->SetBranchAddress("TauVsJetIDSF_3", &TauVsJetIDSF_3);
    tree->SetBranchAddress("TauVsJetIDSF_4", &TauVsJetIDSF_4);

    tree->SetBranchAddress("TauES_1", &TauES_1);
    tree->SetBranchAddress("TauES_2", &TauES_2);
    tree->SetBranchAddress("TauES_3", &TauES_3);
    tree->SetBranchAddress("TauES_4", &TauES_4);

    tree->SetBranchAddress("LHEScaleWeights", &LHEScaleWeights);

    tree->SetBranchAddress("lpt",     &lpt);
    tree->SetBranchAddress("leta",    &leta);
    tree->SetBranchAddress("lphi",    &lphi);
    tree->SetBranchAddress("lmass",   &lmass);
    tree->SetBranchAddress("lq",      &lq);
    tree->SetBranchAddress("lflavor", &lflavor);
    tree->SetBranchAddress("liso",    &liso);
    tree->SetBranchAddress("ld0",     &ld0);
    tree->SetBranchAddress("ldZ",     &ldZ);

    tree->SetBranchAddress("EleID",      &EleID);
    tree->SetBranchAddress("EleID_WP90", &EleID_WP90);
    tree->SetBranchAddress("EleID_WPL",  &EleID_WPL);

    tree->SetBranchAddress("MuID",      &MuID);
    tree->SetBranchAddress("TauIDe",    &TauIDe);
    tree->SetBranchAddress("TauIDm",    &TauIDm);
    tree->SetBranchAddress("TauIDj",    &TauIDj);
    tree->SetBranchAddress("decayMode", &decayMode);
    tree->SetBranchAddress("gen_match", &gen_match);
}

#endif
