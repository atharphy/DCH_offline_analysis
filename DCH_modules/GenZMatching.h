#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include "../include/Kinematics.h"

struct GenZMatch {
    bool matched = false;
    double pt = -1.0, eta = 0.0, phi = 0.0, mass = -1.0;
};

inline GenZMatch matchGenZ(double recoPt1, double recoEta1, double recoPhi1,
                            double recoPt2, double recoEta2, double recoPhi2,
                            double maxDeltaR = 0.3)
{
    GenZMatch result;
    if (!GenPart_pdgId || !GenPart_genPartIdxMother || !GenPart_statusFlags ||
        !GenPart_pt || !GenPart_eta || !GenPart_phi || !GenPart_mass) {
        return result;
    }

    for (size_t i = 0; i < GenPart_pdgId->size(); ++i) {
        if ((*GenPart_pdgId)[i] != 23) continue;
        if ((((*GenPart_statusFlags)[i] >> 13) & 1) == 0) continue;

        std::vector<std::pair<double,double>> daughters;
        for (size_t j = 0; j < GenPart_genPartIdxMother->size(); ++j) {
            if ((*GenPart_genPartIdxMother)[j] != (int)i) continue;
            int dpdg = std::abs((*GenPart_pdgId)[j]);
            if (dpdg == 11 || dpdg == 13) {
                daughters.push_back({(*GenPart_eta)[j], (*GenPart_phi)[j]});
            } else if (dpdg == 15 && GenVisTau_genPartIdxMother && GenVisTau_eta && GenVisTau_phi) {
                for (size_t k = 0; k < GenVisTau_genPartIdxMother->size(); ++k) {
                    if ((*GenVisTau_genPartIdxMother)[k] == (int)j) {
                        daughters.push_back({(*GenVisTau_eta)[k], (*GenVisTau_phi)[k]});
                        break;
                    }
                }
            }
        }

        if (daughters.size() != 2) continue;

        const bool pairingA = getDR(recoEta1, recoPhi1, daughters[0].first, daughters[0].second) < maxDeltaR &&
                               getDR(recoEta2, recoPhi2, daughters[1].first, daughters[1].second) < maxDeltaR;
        const bool pairingB = getDR(recoEta1, recoPhi1, daughters[1].first, daughters[1].second) < maxDeltaR &&
                               getDR(recoEta2, recoPhi2, daughters[0].first, daughters[0].second) < maxDeltaR;

        if (pairingA || pairingB) {
            result.matched = true;
            result.pt   = (*GenPart_pt)[i];
            result.eta  = (*GenPart_eta)[i];
            result.phi  = (*GenPart_phi)[i];
            result.mass = (*GenPart_mass)[i];
            return result;
        }
    }

    return result;
}
