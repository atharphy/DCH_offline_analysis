#pragma once

#include <cmath>
#include <string>

inline bool isValidFlavor(char flavor) { return flavor == 'e' || flavor == 'm' || flavor == 't'; }

inline double ptByIndex(int idx) { if (idx==1) return pt_1; if (idx==2) return pt_2; if (idx==3) return pt_3; if (idx==4) return pt_4; return -99.0; }
inline double etaByIndex(int idx) { if (idx==1) return eta_1; if (idx==2) return eta_2; if (idx==3) return eta_3; if (idx==4) return eta_4; return -99.0; }
inline double phiByIndex(int idx) { if (idx==1) return phi_1; if (idx==2) return phi_2; if (idx==3) return phi_3; if (idx==4) return phi_4; return -99.0; }
inline double massByIndex(int idx) { if (idx==1) return m_1; if (idx==2) return m_2; if (idx==3) return m_3; if (idx==4) return m_4; return 0.0; }
inline double d0ByIndex(int idx) { if (idx==1) return d0_1; if (idx==2) return d0_2; if (idx==3) return d0_3; if (idx==4) return d0_4; return -99.0; }
inline double dZByIndex(int idx) { if (idx==1) return dZ_1; if (idx==2) return dZ_2; if (idx==3) return dZ_3; if (idx==4) return dZ_4; return -99.0; }
inline double isoByIndex(int idx) { if (idx==1) return iso_1; if (idx==2) return iso_2; if (idx==3) return iso_3; if (idx==4) return iso_4; return -99.0; }
inline int qByIndex(int idx) { if (idx==1) return q_1; if (idx==2) return q_2; if (idx==3) return q_3; if (idx==4) return q_4; return 0; }

inline int genPartFlavByIndex(int idx) { if (idx==1) return genPartFlav_1; if (idx==2) return genPartFlav_2; if (idx==3) return genPartFlav_3; if (idx==4) return genPartFlav_4; return -99; }
inline double truthPtByIndex(int idx) { if (idx==1) return pt_1_tr; if (idx==2) return pt_2_tr; if (idx==3) return pt_3_tr; if (idx==4) return pt_4_tr; return -99.0; }
inline double truthPhiByIndex(int idx) { if (idx==1) return phi_1_tr; if (idx==2) return phi_2_tr; if (idx==3) return phi_3_tr; if (idx==4) return phi_4_tr; return -99.0; }

inline double idSFByIndex(int idx) { if (idx==1) return IDSF_1; if (idx==2) return IDSF_2; if (idx==3) return IDSF_3; if (idx==4) return IDSF_4; return 1.0; }
inline double isoSFByIndex(int idx) { if (idx==1) return ISOSF_1; if (idx==2) return ISOSF_2; if (idx==3) return ISOSF_3; if (idx==4) return ISOSF_4; return 1.0; }
inline double tauEleSFByIndex(int idx) { if (idx==1) return TauVsEleIDSF_1; if (idx==2) return TauVsEleIDSF_2; if (idx==3) return TauVsEleIDSF_3; if (idx==4) return TauVsEleIDSF_4; return 1.0; }
inline double tauMuSFByIndex(int idx) { if (idx==1) return TauVsMuIDSF_1; if (idx==2) return TauVsMuIDSF_2; if (idx==3) return TauVsMuIDSF_3; if (idx==4) return TauVsMuIDSF_4; return 1.0; }
inline double tauJetSFByIndex(int idx) { if (idx==1) return TauVsJetIDSF_1; if (idx==2) return TauVsJetIDSF_2; if (idx==3) return TauVsJetIDSF_3; if (idx==4) return TauVsJetIDSF_4; return 1.0; }

inline bool trigPassedByIndex(int idx) { if (idx==1) return isTrig_1==1 || isTrig_1==2; if (idx==2) return isTrig_1==-1 || isTrig_1==2; if (idx==3) return isTrig_2==1 || isTrig_2==2; if (idx==4) return isTrig_2==-1 || isTrig_2==2; return false; }
inline double trigSFByIndex(int idx) { if (idx==1) return TrigSF_1; if (idx==2) return TrigSF_2; if (idx==3) return TrigSF_3; if (idx==4) return TrigSF_4; return 1.0; }

inline std::string dilepZFlag(const std::string& catstr) {
    if (catstr.size() != 2) return "";
    if (catstr[0] != catstr[1] || catstr[0] == 't') return "Zveto";
    return std::fabs((LepV(1) + LepV(2)).M() - 91.2) <= 10.0 ? "Zwin" : "Zveto";
}
