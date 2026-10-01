#pragma once

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "TFile.h"
#include "correction.h"

struct OfficialMETCorrections {
    std::unique_ptr<correction::CorrectionSet> set;
    correction::Correction::Ref pt;
    correction::Correction::Ref phi;
};

inline std::unique_ptr<OfficialMETCorrections> loadOfficialMETCorrections(const std::string& year, bool isData) {
    const std::string fileName = "json_files/met_" + year + ".json.gz";
    const std::string suffix = isData ? "data" : "mc";
    try {
        std::unique_ptr<OfficialMETCorrections> result(new OfficialMETCorrections());
        result->set = correction::CorrectionSet::from_file(fileName);
        result->pt = result->set->at("pt_metphicorr_pfmet_" + suffix);
        result->phi = result->set->at("phi_metphicorr_pfmet_" + suffix);
        std::cout << "Loaded official PF MET corrections from:\n  " << fileName << std::endl;
        return result;
    }
    catch (const std::exception& error) {
        std::cerr << "ERROR: cannot load official PF MET corrections from:\n  " << fileName << "\n  " << error.what() << std::endl;
        return nullptr;
    }
}

inline bool applyOfficialMETCorrection(const OfficialMETCorrections& corrections) {
    static int nDiag = 0;
    if (!corrections.pt || !corrections.phi || !std::isfinite(met) || !std::isfinite(metphi) || met < 0.0 || nPV < 0) {
        if (nDiag < 10) {
            ++nDiag;
            std::cerr << "  [metcorr-diag] pre-check failed: pt=" << (bool)corrections.pt << " phi=" << (bool)corrections.phi
                      << " met=" << met << " metphi=" << metphi << " nPV=" << nPV << " nPVGood=" << nPVGood << " run=" << run << std::endl;
        }
        return false;
    }
    const double rawMet = met;
    const double rawMetPhi = metphi;
    const double npvs = static_cast<double>(nPV);
    const double runNumber = static_cast<double>(run);
    try {
        const double correctedMet = corrections.pt->evaluate({rawMet, rawMetPhi, npvs, runNumber});
        const double correctedMetPhi = corrections.phi->evaluate({rawMet, rawMetPhi, npvs, runNumber});
        if (!std::isfinite(correctedMet) || !std::isfinite(correctedMetPhi) || correctedMet < 0.0) {
            if (nDiag < 10) {
                ++nDiag;
                std::cerr << "  [metcorr-diag] post-check failed: correctedMet=" << correctedMet << " correctedMetPhi=" << correctedMetPhi
                          << " input met=" << rawMet << " metphi=" << rawMetPhi << " nPV=" << npvs << " run=" << runNumber << std::endl;
            }
            return false;
        }
        met = correctedMet;
        metphi = correctedMetPhi;
        return true;
    }
    catch (const std::exception& error) {
        if (nDiag < 10) {
            ++nDiag;
            std::cerr << "  [metcorr-diag] evaluate() threw: " << error.what()
                      << " (met=" << rawMet << " metphi=" << rawMetPhi << " nPV=" << npvs << " run=" << runNumber << ")" << std::endl;
        }
        return false;
    }
}
