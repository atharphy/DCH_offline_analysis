// Standalone smoke test for LeptonTauSFUncertainty.h -- not part of the
// analysis pipeline, just a compile+sanity check.
//
// Run (from the offline/ top-level directory):
//   root -l -b -q 'sf_uncertainty/TestSFUncertainty.C+("2018")'

#include "correction.h"
#if defined(__CLING__)
R__LOAD_LIBRARY(libcorrectionlib.so)
#endif

#include "LeptonTauSFUncertainty.h"

#include <iostream>

void TestSFUncertainty(std::string year = "2018") {
    const std::string jsonDir = "sf_uncertainty/json";
    const std::string rootDir = "sf_uncertainty/root";

    auto eReader = getCachedElectronSFUncReader(jsonDir, year);
    auto eTrigReader = getCachedElectronTrigSFUncReader(rootDir, year);
    auto muReader = getCachedMuonSFUncReader(jsonDir, year);
    auto tauReader = getCachedTauSFUncReader(jsonDir, year);

    std::cout << "\n=== Electron (pt=40, eta=1.2) ===" << std::endl;
    auto eReco = eReader->getRecoUnc(40.0, 1.2);
    auto eIdIso = eReader->getIdIsoUnc(40.0, 1.2, "wp90noiso");
    auto eTrig = eTrigReader->getTrigUnc(40.0, 1.2);
    std::cout << "  reco:   +" << eReco.up << " / -" << eReco.down << std::endl;
    std::cout << "  idIso:  +" << eIdIso.up << " / -" << eIdIso.down << std::endl;
    std::cout << "  trig:   +" << eTrig.up << " / -" << eTrig.down << std::endl;

    std::cout << "\n=== Muon (pt=35, eta=-0.8) ===" << std::endl;
    auto muId = muReader->getIdUnc(35.0, -0.8);
    auto muIso = muReader->getIsoUnc(35.0, -0.8);
    auto muTrig = muReader->getTrigUnc(35.0, -0.8);
    std::cout << "  id:     +" << muId.up << " / -" << muId.down << std::endl;
    std::cout << "  iso:    +" << muIso.up << " / -" << muIso.down << std::endl;
    std::cout << "  trig:   +" << muTrig.up << " / -" << muTrig.down << std::endl;

    std::cout << "\n=== Tau (pt=45, eta=1.0, dm=1, genmatch=5 [genuine]) ===" << std::endl;
    auto tVsE = tauReader->getVsEleUnc(1.0, 5);
    auto tVsMu = tauReader->getVsMuUnc(1.0, 5);
    auto tVsJet = tauReader->getVsJetUnc(45.0, 1, 5);
    auto tES = tauReader->getEnergyScaleUnc(45.0, 1.0, 1, 5);
    std::cout << "  vsEle:  +" << tVsE.up << " / -" << tVsE.down << std::endl;
    std::cout << "  vsMu:   +" << tVsMu.up << " / -" << tVsMu.down << std::endl;
    std::cout << "  vsJet:  +" << tVsJet.up << " / -" << tVsJet.down << std::endl;
    std::cout << "  ES:     +" << tES.up << " / -" << tES.down << std::endl;

    std::cout << "\n=== Tau (pt=45, eta=1.0, dm=1, genmatch=1 [e->tau fake]) ===" << std::endl;
    auto tVsEFake = tauReader->getVsEleUnc(1.0, 1);
    std::cout << "  vsEle:  +" << tVsEFake.up << " / -" << tVsEFake.down << std::endl;

    std::cout << "\ndone" << std::endl;
}
