#pragma once

#include <cmath>

inline bool getGenBosonMomentum(double& fullPx, double& fullPy, double& visPx, double& visPy) {
    if (!GenPart_pdgId || !GenPart_status || !GenPart_statusFlags || !GenPart_pt || !GenPart_phi) return false;

    fullPx = fullPy = visPx = visPy = 0.0;
    bool found = false;

    for (size_t i = 0; i < GenPart_pdgId->size(); ++i) {
        const int pdg = std::abs((*GenPart_pdgId)[i]);
        const bool isLepton = (pdg == 11 || pdg == 13);
        const bool isNeutrino = (pdg == 12 || pdg == 14 || pdg == 16);
        if (!isLepton && !isNeutrino) continue;

        const bool fromHardProcessFinalState =
            (*GenPart_status)[i] == 1 && (((*GenPart_statusFlags)[i] >> 8) & 1);
        const bool isDirectHardProcessTauDecayProduct =
            (((*GenPart_statusFlags)[i] >> 10) & 1);

        if (!((fromHardProcessFinalState && (isLepton || isNeutrino)) || isDirectHardProcessTauDecayProduct)) continue;

        const double px = (*GenPart_pt)[i] * std::cos((*GenPart_phi)[i]);
        const double py = (*GenPart_pt)[i] * std::sin((*GenPart_phi)[i]);
        fullPx += px;
        fullPy += py;
        found = true;

        if ((fromHardProcessFinalState && isLepton) || (isDirectHardProcessTauDecayProduct && !isNeutrino)) {
            visPx += px;
            visPy += py;
        }
    }
    return found;
}
