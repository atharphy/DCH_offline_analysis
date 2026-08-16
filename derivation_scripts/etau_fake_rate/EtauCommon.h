#pragma once

#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TSystem.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "../../filemap/FileMap.h"
#include "../../include/cat.h"

using ROOT::VecOps::RVec;

const std::string ETAU_OUTPUT_BASE = "results";

const std::vector<double> ETAU_PT_EDGES = {
    20., 22., 24., 26., 28., 30., 32., 35., 40.,
    45., 50., 60., 70., 85., 100., 125., 150., 200.
};

const double ETAU_DR_SEPARATION = 0.5;

inline double etauDPhiWrap(double phi1, double phi2) {
    const double pi = TMath::Pi();
    return std::min(std::fabs(phi2 - phi1), 2.0 * pi - std::fabs(phi2 - phi1));
}

inline double etauDR(double eta1, double phi1, double eta2, double phi2) {
    return std::sqrt(std::pow(etauDPhiWrap(phi1, phi2), 2) + std::pow(eta2 - eta1, 2));
}

inline bool etauLeg1IsElectron(int cat) {
    const std::string name = numberToCat(cat);
    return !name.empty() && name[0] == 'e';
}

inline std::array<double, 3> etauFindProbe(double tagEta, double tagPhi,
                                            const RVec<double>& lpt, const RVec<double>& leta, const RVec<double>& lphi,
                                            const RVec<int>& lflavor, const RVec<int>& genMatch,
                                            const RVec<int>& TauIDm, const RVec<int>& TauIDe, const RVec<int>& TauIDj) {
    const std::size_t size = std::min({lpt.size(), leta.size(), lphi.size(), lflavor.size(), genMatch.size(), TauIDm.size(), TauIDe.size(), TauIDj.size()});
    for (std::size_t i = 0; i < size; ++i) {
        if (std::abs(lflavor[i]) != 15) continue;
        if (genMatch[i] != 1) continue;
        if (TauIDm[i] <= 4) continue;
        if (TauIDe[i] <= 8) continue;
        if (etauDR(leta[i], lphi[i], tagEta, tagPhi) <= ETAU_DR_SEPARATION) continue;
        return {lpt[i], leta[i], static_cast<double>(TauIDj[i])};
    }
    return {-1.0, -99.0, -1.0};
}

inline bool etauIsDYFile(const std::string& fileName) {
    const std::string baseName = gSystem->BaseName(fileName.c_str());
    return baseName.find("DYJetsToLL") != std::string::npos;
}

inline std::vector<std::string> etauDYFiles(const std::string& year) {
    const auto fileMap = getFileMap(year);
    std::vector<std::string> result;
    for (const std::string& group : {"DY", "DY10_50"}) {
        const auto found = fileMap.find(group);
        if (found == fileMap.end()) continue;
        for (const std::string& fileName : found->second) {
            if (!fileName.empty()) result.push_back(fileName);
        }
    }
    return result;
}

inline std::unique_ptr<TH1D> etauCloneHistogram(const TH1D& source, const std::string& name) {
    std::unique_ptr<TH1D> result(dynamic_cast<TH1D*>(source.Clone(name.c_str())));
    result->SetDirectory(nullptr);
    return result;
}

inline std::unique_ptr<TH1D> etauMakeFakeRate(const TH1D& tight, const TH1D& denominator, const std::string& name) {
    std::unique_ptr<TH1D> result = etauCloneHistogram(tight, name);
    result->Divide(&tight, &denominator, 1.0, 1.0, "B");
    result->SetTitle("");
    result->GetXaxis()->SetTitle("p_{T}^{#tau} [GeV]");
    result->GetYaxis()->SetTitle("Electron-to-tau fake rate");
    result->SetMinimum(0.0);
    result->SetMaximum(1.0);
    return result;
}

inline void etauCopyToTwoDimensionalPayload(TH2D& destination, const TH1D& barrel, const TH1D& endcap) {
    for (int xbin = 1; xbin <= destination.GetNbinsX(); ++xbin) {
        destination.SetBinContent(xbin, 1, barrel.GetBinContent(xbin));
        destination.SetBinError(xbin, 1, barrel.GetBinError(xbin));
        destination.SetBinContent(xbin, 2, endcap.GetBinContent(xbin));
        destination.SetBinError(xbin, 2, endcap.GetBinError(xbin));
    }
}
