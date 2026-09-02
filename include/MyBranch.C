#ifndef MYBRANCH_C
#define MYBRANCH_C

#include <TTree.h>
#include <vector>

using std::vector;

Int_t nElectron, nMuon, nTau;
Int_t nLooseElectron, nLooseMuon, nLooseTau;
Int_t evt;
Int_t lumiBlock;
ULong64_t run;

Int_t cat, gen_cat;
Int_t nPV;
Int_t nPVGood;
Int_t nPU, nPUEOOT, nPULOOT;
Double_t nPUtrue;
Int_t nJet;
Int_t genPartFlav_1, genPartFlav_2, genPartFlav_3, genPartFlav_4;

Int_t electronTriggerWord, muonTriggerWord;
Int_t whichTriggerWord, whichTriggerWordSubL;

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
Double_t MET_pt_UnclUp, MET_phi_UnclUp;
Double_t MET_pt_UnclDown, MET_phi_UnclDown;
Double_t met_UnclX, met_UnclY;
Double_t MET_T1Smear_pt, MET_T1Smear_phi;
Double_t metcov00, metcov01, metcov10, metcov11;
Double_t njets;
Double_t nbtagL, nbtagM, nbtagT;

Double_t pt_1, pt_2, pt_3, pt_4;
Double_t eta_1, eta_2, eta_3, eta_4;
Double_t phi_1, phi_2, phi_3, phi_4;
Double_t pt_1_tr, pt_2_tr, pt_3_tr, pt_4_tr;
Double_t eta_1_tr, eta_2_tr, eta_3_tr, eta_4_tr;
Double_t phi_1_tr, phi_2_tr, phi_3_tr, phi_4_tr;
Double_t m_1_tr, m_2_tr;
Double_t m_1, m_2, m_3, m_4;
Double_t q_1, q_2, q_3, q_4;
Double_t iso_1, iso_2, iso_3, iso_4;
Double_t d0_1, d0_2, d0_3, d0_4;
Double_t dZ_1, dZ_2, dZ_3, dZ_4;
Int_t decayMode_1, decayMode_2, decayMode_3, decayMode_4;
Int_t nTrackerLayers_1, nTrackerLayers_2, nTrackerLayers_3, nTrackerLayers_4;

Double_t mDCH1_sv, mtDCH1_sv, mDCH2_sv, mtDCH2_sv;
Double_t mDCH1_sv4, mtDCH1_sv4, mDCH2_sv4, mtDCH2_sv4;

Double_t IDSF_1, IDSF_2, IDSF_3, IDSF_4;
Double_t ISOSF_1, ISOSF_2, ISOSF_3, ISOSF_4;
Double_t TrigSF_1, TrigSF_2, TrigSF_3, TrigSF_4;
Double_t isTrig_1, isTrig_2, isDoubleTrig;

Double_t TauVsEleIDSF_1, TauVsEleIDSF_2, TauVsEleIDSF_3, TauVsEleIDSF_4;
Double_t TauVsMuIDSF_1,  TauVsMuIDSF_2,  TauVsMuIDSF_3,  TauVsMuIDSF_4;
Double_t TauVsJetIDSF_1, TauVsJetIDSF_2, TauVsJetIDSF_3, TauVsJetIDSF_4;

Double_t TauES_1, TauES_2, TauES_3, TauES_4;

Double_t IDSF_Up_1, IDSF_Up_2, IDSF_Up_3, IDSF_Up_4;
Double_t IDSF_Down_1, IDSF_Down_2, IDSF_Down_3, IDSF_Down_4;
Double_t ISOSF_Up_1, ISOSF_Up_2, ISOSF_Up_3, ISOSF_Up_4;
Double_t ISOSF_Down_1, ISOSF_Down_2, ISOSF_Down_3, ISOSF_Down_4;
Double_t TrigSF_Up_1, TrigSF_Up_2, TrigSF_Up_3, TrigSF_Up_4;
Double_t TrigSF_Down_1, TrigSF_Down_2, TrigSF_Down_3, TrigSF_Down_4;

Double_t TauVsEleIDSF_Up_1, TauVsEleIDSF_Up_2, TauVsEleIDSF_Up_3, TauVsEleIDSF_Up_4;
Double_t TauVsEleIDSF_Down_1, TauVsEleIDSF_Down_2, TauVsEleIDSF_Down_3, TauVsEleIDSF_Down_4;
Double_t TauVsMuIDSF_Up_1, TauVsMuIDSF_Up_2, TauVsMuIDSF_Up_3, TauVsMuIDSF_Up_4;
Double_t TauVsMuIDSF_Down_1, TauVsMuIDSF_Down_2, TauVsMuIDSF_Down_3, TauVsMuIDSF_Down_4;
Double_t TauVsJetIDSF_Up_1, TauVsJetIDSF_Up_2, TauVsJetIDSF_Up_3, TauVsJetIDSF_Up_4;
Double_t TauVsJetIDSF_Down_1, TauVsJetIDSF_Down_2, TauVsJetIDSF_Down_3, TauVsJetIDSF_Down_4;

Double_t TauES_Up_1, TauES_Up_2, TauES_Up_3, TauES_Up_4;
Double_t TauES_Down_1, TauES_Down_2, TauES_Down_3, TauES_Down_4;

vector<double> *LHEScaleWeights = nullptr;

vector<double> *GenPart_pt               = nullptr;
vector<double> *GenPart_eta              = nullptr;
vector<double> *GenPart_phi              = nullptr;
vector<double> *GenPart_mass             = nullptr;
vector<int>    *GenPart_pdgId            = nullptr;
vector<int>    *GenPart_genPartIdxMother = nullptr;
vector<int>    *GenPart_status           = nullptr;
vector<int>    *GenPart_statusFlags      = nullptr;
vector<int>    *lGenPart_status          = nullptr;
vector<int>    *lGenPart_statusFlags     = nullptr;

vector<double> *GenVisTau_pt               = nullptr;
vector<double> *GenVisTau_eta              = nullptr;
vector<double> *GenVisTau_phi              = nullptr;
vector<double> *GenVisTau_mass             = nullptr;
vector<int>    *GenVisTau_charge           = nullptr;
vector<int>    *GenVisTau_genPartIdxMother = nullptr;
vector<int>    *GenVisTau_status           = nullptr;

Double_t GenMET_pt, GenMET_phi;

vector<double> *Jet_pt            = nullptr;
vector<double> *Jet_eta           = nullptr;
vector<double> *Jet_phi           = nullptr;
vector<double> *Jet_mass          = nullptr;
vector<int>    *Jet_jetId         = nullptr;
vector<int>    *Jet_partonFlavour = nullptr;
vector<int>    *Jet_hadronFlavour = nullptr;
vector<double> *Jet_btagDeepB     = nullptr;

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
    tree->SetBranchAddress("nLooseElectron", &nLooseElectron);
    tree->SetBranchAddress("nLooseMuon",     &nLooseMuon);
    tree->SetBranchAddress("nLooseTau",      &nLooseTau);
    tree->SetBranchAddress("evt",       &evt);
    tree->SetBranchAddress("lumi",      &lumiBlock);
    tree->SetBranchAddress("run",       &run);

    tree->SetBranchAddress("cat",     &cat);
    tree->SetBranchAddress("gen_cat", &gen_cat);

    tree->SetBranchAddress("nPU",     &nPU);
    tree->SetBranchAddress("nPUEOOT", &nPUEOOT);
    tree->SetBranchAddress("nPULOOT", &nPULOOT);
    tree->SetBranchAddress("nPUtrue", &nPUtrue);
    tree->SetBranchAddress("nJet",    &nJet);

    tree->SetBranchAddress("electronTriggerWord",   &electronTriggerWord);
    tree->SetBranchAddress("muonTriggerWord",       &muonTriggerWord);
    tree->SetBranchAddress("whichTriggerWord",      &whichTriggerWord);
    tree->SetBranchAddress("whichTriggerWordSubL",  &whichTriggerWordSubL);

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
    tree->SetBranchAddress("MET_pt_UnclUp",    &MET_pt_UnclUp);
    tree->SetBranchAddress("MET_phi_UnclUp",   &MET_phi_UnclUp);
    tree->SetBranchAddress("MET_pt_UnclDown",  &MET_pt_UnclDown);
    tree->SetBranchAddress("MET_phi_UnclDown", &MET_phi_UnclDown);
    tree->SetBranchAddress("met_UnclX",        &met_UnclX);
    tree->SetBranchAddress("met_UnclY",        &met_UnclY);
    tree->SetBranchAddress("MET_T1Smear_pt",   &MET_T1Smear_pt);
    tree->SetBranchAddress("MET_T1Smear_phi",  &MET_T1Smear_phi);
    tree->SetBranchAddress("metcov00", &metcov00);
    tree->SetBranchAddress("metcov01", &metcov01);
    tree->SetBranchAddress("metcov10", &metcov10);
    tree->SetBranchAddress("metcov11", &metcov11);
    tree->SetBranchAddress("njets",        &njets);
    tree->SetBranchAddress("nbtagL", &nbtagL);
    tree->SetBranchAddress("nbtagM", &nbtagM);
    tree->SetBranchAddress("nbtagT", &nbtagT);

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

    tree->SetBranchAddress("m_1_tr", &m_1_tr);
    tree->SetBranchAddress("m_2_tr", &m_2_tr);

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

    tree->SetBranchAddress("decayMode_1", &decayMode_1);
    tree->SetBranchAddress("decayMode_2", &decayMode_2);
    tree->SetBranchAddress("decayMode_3", &decayMode_3);
    tree->SetBranchAddress("decayMode_4", &decayMode_4);

    tree->SetBranchAddress("nTrackerLayers_1", &nTrackerLayers_1);
    tree->SetBranchAddress("nTrackerLayers_2", &nTrackerLayers_2);
    tree->SetBranchAddress("nTrackerLayers_3", &nTrackerLayers_3);
    tree->SetBranchAddress("nTrackerLayers_4", &nTrackerLayers_4);

    tree->SetBranchAddress("mDCH1_sv",   &mDCH1_sv);
    tree->SetBranchAddress("mtDCH1_sv",  &mtDCH1_sv);
    tree->SetBranchAddress("mDCH2_sv",   &mDCH2_sv);
    tree->SetBranchAddress("mtDCH2_sv",  &mtDCH2_sv);
    tree->SetBranchAddress("mDCH1_sv4",  &mDCH1_sv4);
    tree->SetBranchAddress("mtDCH1_sv4", &mtDCH1_sv4);
    tree->SetBranchAddress("mDCH2_sv4",  &mDCH2_sv4);
    tree->SetBranchAddress("mtDCH2_sv4", &mtDCH2_sv4);

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

    tree->SetBranchAddress("IDSF_Up_1", &IDSF_Up_1);
    tree->SetBranchAddress("IDSF_Up_2", &IDSF_Up_2);
    tree->SetBranchAddress("IDSF_Up_3", &IDSF_Up_3);
    tree->SetBranchAddress("IDSF_Up_4", &IDSF_Up_4);
    tree->SetBranchAddress("IDSF_Down_1", &IDSF_Down_1);
    tree->SetBranchAddress("IDSF_Down_2", &IDSF_Down_2);
    tree->SetBranchAddress("IDSF_Down_3", &IDSF_Down_3);
    tree->SetBranchAddress("IDSF_Down_4", &IDSF_Down_4);

    tree->SetBranchAddress("ISOSF_Up_1", &ISOSF_Up_1);
    tree->SetBranchAddress("ISOSF_Up_2", &ISOSF_Up_2);
    tree->SetBranchAddress("ISOSF_Up_3", &ISOSF_Up_3);
    tree->SetBranchAddress("ISOSF_Up_4", &ISOSF_Up_4);
    tree->SetBranchAddress("ISOSF_Down_1", &ISOSF_Down_1);
    tree->SetBranchAddress("ISOSF_Down_2", &ISOSF_Down_2);
    tree->SetBranchAddress("ISOSF_Down_3", &ISOSF_Down_3);
    tree->SetBranchAddress("ISOSF_Down_4", &ISOSF_Down_4);

    tree->SetBranchAddress("TrigSF_Up_1", &TrigSF_Up_1);
    tree->SetBranchAddress("TrigSF_Up_2", &TrigSF_Up_2);
    tree->SetBranchAddress("TrigSF_Up_3", &TrigSF_Up_3);
    tree->SetBranchAddress("TrigSF_Up_4", &TrigSF_Up_4);
    tree->SetBranchAddress("TrigSF_Down_1", &TrigSF_Down_1);
    tree->SetBranchAddress("TrigSF_Down_2", &TrigSF_Down_2);
    tree->SetBranchAddress("TrigSF_Down_3", &TrigSF_Down_3);
    tree->SetBranchAddress("TrigSF_Down_4", &TrigSF_Down_4);

    tree->SetBranchAddress("TauVsEleIDSF_Up_1", &TauVsEleIDSF_Up_1);
    tree->SetBranchAddress("TauVsEleIDSF_Up_2", &TauVsEleIDSF_Up_2);
    tree->SetBranchAddress("TauVsEleIDSF_Up_3", &TauVsEleIDSF_Up_3);
    tree->SetBranchAddress("TauVsEleIDSF_Up_4", &TauVsEleIDSF_Up_4);
    tree->SetBranchAddress("TauVsEleIDSF_Down_1", &TauVsEleIDSF_Down_1);
    tree->SetBranchAddress("TauVsEleIDSF_Down_2", &TauVsEleIDSF_Down_2);
    tree->SetBranchAddress("TauVsEleIDSF_Down_3", &TauVsEleIDSF_Down_3);
    tree->SetBranchAddress("TauVsEleIDSF_Down_4", &TauVsEleIDSF_Down_4);

    tree->SetBranchAddress("TauVsMuIDSF_Up_1", &TauVsMuIDSF_Up_1);
    tree->SetBranchAddress("TauVsMuIDSF_Up_2", &TauVsMuIDSF_Up_2);
    tree->SetBranchAddress("TauVsMuIDSF_Up_3", &TauVsMuIDSF_Up_3);
    tree->SetBranchAddress("TauVsMuIDSF_Up_4", &TauVsMuIDSF_Up_4);
    tree->SetBranchAddress("TauVsMuIDSF_Down_1", &TauVsMuIDSF_Down_1);
    tree->SetBranchAddress("TauVsMuIDSF_Down_2", &TauVsMuIDSF_Down_2);
    tree->SetBranchAddress("TauVsMuIDSF_Down_3", &TauVsMuIDSF_Down_3);
    tree->SetBranchAddress("TauVsMuIDSF_Down_4", &TauVsMuIDSF_Down_4);

    tree->SetBranchAddress("TauVsJetIDSF_Up_1", &TauVsJetIDSF_Up_1);
    tree->SetBranchAddress("TauVsJetIDSF_Up_2", &TauVsJetIDSF_Up_2);
    tree->SetBranchAddress("TauVsJetIDSF_Up_3", &TauVsJetIDSF_Up_3);
    tree->SetBranchAddress("TauVsJetIDSF_Up_4", &TauVsJetIDSF_Up_4);
    tree->SetBranchAddress("TauVsJetIDSF_Down_1", &TauVsJetIDSF_Down_1);
    tree->SetBranchAddress("TauVsJetIDSF_Down_2", &TauVsJetIDSF_Down_2);
    tree->SetBranchAddress("TauVsJetIDSF_Down_3", &TauVsJetIDSF_Down_3);
    tree->SetBranchAddress("TauVsJetIDSF_Down_4", &TauVsJetIDSF_Down_4);

    tree->SetBranchAddress("TauES_Up_1", &TauES_Up_1);
    tree->SetBranchAddress("TauES_Up_2", &TauES_Up_2);
    tree->SetBranchAddress("TauES_Up_3", &TauES_Up_3);
    tree->SetBranchAddress("TauES_Up_4", &TauES_Up_4);
    tree->SetBranchAddress("TauES_Down_1", &TauES_Down_1);
    tree->SetBranchAddress("TauES_Down_2", &TauES_Down_2);
    tree->SetBranchAddress("TauES_Down_3", &TauES_Down_3);
    tree->SetBranchAddress("TauES_Down_4", &TauES_Down_4);

    tree->SetBranchAddress("LHEScaleWeights", &LHEScaleWeights);

    tree->SetBranchAddress("GenPart_pt",               &GenPart_pt);
    tree->SetBranchAddress("GenPart_eta",              &GenPart_eta);
    tree->SetBranchAddress("GenPart_phi",              &GenPart_phi);
    tree->SetBranchAddress("GenPart_mass",             &GenPart_mass);
    tree->SetBranchAddress("GenPart_pdgId",            &GenPart_pdgId);
    tree->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    tree->SetBranchAddress("GenPart_status",           &GenPart_status);
    tree->SetBranchAddress("GenPart_statusFlags",      &GenPart_statusFlags);
    tree->SetBranchAddress("lGenPart_status",          &lGenPart_status);
    tree->SetBranchAddress("lGenPart_statusFlags",     &lGenPart_statusFlags);

    tree->SetBranchAddress("GenVisTau_pt",               &GenVisTau_pt);
    tree->SetBranchAddress("GenVisTau_eta",              &GenVisTau_eta);
    tree->SetBranchAddress("GenVisTau_phi",              &GenVisTau_phi);
    tree->SetBranchAddress("GenVisTau_mass",             &GenVisTau_mass);
    tree->SetBranchAddress("GenVisTau_charge",           &GenVisTau_charge);
    tree->SetBranchAddress("GenVisTau_genPartIdxMother", &GenVisTau_genPartIdxMother);
    tree->SetBranchAddress("GenVisTau_status",           &GenVisTau_status);

    tree->SetBranchAddress("GenMET_pt",  &GenMET_pt);
    tree->SetBranchAddress("GenMET_phi", &GenMET_phi);

    tree->SetBranchAddress("Jet_pt",            &Jet_pt);
    tree->SetBranchAddress("Jet_eta",           &Jet_eta);
    tree->SetBranchAddress("Jet_phi",           &Jet_phi);
    tree->SetBranchAddress("Jet_mass",          &Jet_mass);
    tree->SetBranchAddress("Jet_jetId",         &Jet_jetId);
    tree->SetBranchAddress("Jet_partonFlavour", &Jet_partonFlavour);
    tree->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
    tree->SetBranchAddress("Jet_btagDeepB",     &Jet_btagDeepB);

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
