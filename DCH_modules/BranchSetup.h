#pragma once

#include "TTree.h"

inline void enableBranches(TTree* tree, bool isData) {
    tree->SetBranchStatus("*", 0);
    const char* commonBranches[] = {
        "cat","run","met","metphi","nPV","nPVGood",
        "pt_1","pt_2","pt_3","pt_4",
        "eta_1","eta_2","eta_3","eta_4",
        "phi_1","phi_2","phi_3","phi_4",
        "m_1","m_2","m_3","m_4",
        "q_1","q_2","q_3","q_4",
        "d0_1","d0_2","d0_3","d0_4",
        "dZ_1","dZ_2","dZ_3","dZ_4",
        "iso_1","iso_2","iso_3","iso_4",
        "TauES_1","TauES_2","TauES_3","TauES_4",
        "metcov00","metcov01","metcov10","metcov11"
    };
    for (const char* branchName : commonBranches) if (tree->GetBranch(branchName)) tree->SetBranchStatus(branchName, 1);
    if (isData) return;

    const char* mcBranches[] = {
        "Generator_weight","brWeight","L1PreFiringWeight_Nom","L1PreFiringWeight_Up","L1PreFiringWeight_Down",
        "weightPUtruejson","weightPUtruejson_up","weightPUtruejson_down",
        "genPartFlav_1","genPartFlav_2","genPartFlav_3","genPartFlav_4",
        "pt_1_tr","pt_2_tr","pt_3_tr","pt_4_tr",
        "eta_1_tr","eta_2_tr","eta_3_tr","eta_4_tr",
        "phi_1_tr","phi_2_tr","phi_3_tr","phi_4_tr",
        "IDSF_1","IDSF_2","IDSF_3","IDSF_4",
        "ISOSF_1","ISOSF_2","ISOSF_3","ISOSF_4",
        "TrigSF_1","TrigSF_2","TrigSF_3","TrigSF_4",
        "isTrig_1","isTrig_2",
        "TauVsEleIDSF_1","TauVsEleIDSF_2","TauVsEleIDSF_3","TauVsEleIDSF_4",
        "TauVsMuIDSF_1","TauVsMuIDSF_2","TauVsMuIDSF_3","TauVsMuIDSF_4",
        "TauVsJetIDSF_1","TauVsJetIDSF_2","TauVsJetIDSF_3","TauVsJetIDSF_4",
        "IDSF_Up_1","IDSF_Up_2","IDSF_Up_3","IDSF_Up_4",
        "IDSF_Down_1","IDSF_Down_2","IDSF_Down_3","IDSF_Down_4",
        "ISOSF_Up_1","ISOSF_Up_2","ISOSF_Up_3","ISOSF_Up_4",
        "ISOSF_Down_1","ISOSF_Down_2","ISOSF_Down_3","ISOSF_Down_4",
        "TrigSF_Up_1","TrigSF_Up_2","TrigSF_Up_3","TrigSF_Up_4",
        "TrigSF_Down_1","TrigSF_Down_2","TrigSF_Down_3","TrigSF_Down_4",
        "TauVsEleIDSF_Up_1","TauVsEleIDSF_Up_2","TauVsEleIDSF_Up_3","TauVsEleIDSF_Up_4",
        "TauVsEleIDSF_Down_1","TauVsEleIDSF_Down_2","TauVsEleIDSF_Down_3","TauVsEleIDSF_Down_4",
        "TauVsMuIDSF_Up_1","TauVsMuIDSF_Up_2","TauVsMuIDSF_Up_3","TauVsMuIDSF_Up_4",
        "TauVsMuIDSF_Down_1","TauVsMuIDSF_Down_2","TauVsMuIDSF_Down_3","TauVsMuIDSF_Down_4",
        "TauVsJetIDSF_Up_1","TauVsJetIDSF_Up_2","TauVsJetIDSF_Up_3","TauVsJetIDSF_Up_4",
        "TauVsJetIDSF_Down_1","TauVsJetIDSF_Down_2","TauVsJetIDSF_Down_3","TauVsJetIDSF_Down_4",
        "TauES_Up_1","TauES_Up_2","TauES_Up_3","TauES_Up_4",
        "TauES_Down_1","TauES_Down_2","TauES_Down_3","TauES_Down_4",
        "GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass",
        "GenPart_pdgId","GenPart_genPartIdxMother","GenPart_status","GenPart_statusFlags",
        "GenVisTau_genPartIdxMother","GenVisTau_eta","GenVisTau_phi",
        "gen_cat"
    };
    for (const char* branchName : mcBranches) if (tree->GetBranch(branchName)) tree->SetBranchStatus(branchName, 1);
}
