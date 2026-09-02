// Plots tau ID/energy-scale scale factors and their uncertainties:
//   - vsEle, vsMu ID SF: 1D vs eta (their only kinematic dependence at the
//     "pt"-flag branch this framework uses -- see LeptonTauSFUncertainty.h's
//     own comment on this).
//   - vsJet ID SF and energy scale: 2D vs (pt, decay mode), since both
//     genuinely depend on decay mode.
// All evaluated at genmatch=5 (genuine hadronic tau), matching the only
// case DCH_tight.C/DCH_tauFR.C ever apply these SFs in (see
// computeConfigurationObjectSFSystematics's genPartFlavByIndex(idx)==5
// gate).
//
// Run: root -l -b -q systematics_plots/PlotTauSF.C+
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>

#include "correction.h"
#include "SFPlotStyle.h"

namespace {
const std::string kJsonDir = "sf_uncertainty/json";
const std::string kOutDir = "systematics_plots/output";
const std::vector<std::string> kYears = {"2016preVFP", "2016postVFP", "2017", "2018"};
const int kGenMatch = 5;  // genuine hadronic tau -- the only case these SFs are ever applied in this framework

std::string eraFileTag(const std::string& year) { return year == "2016postVFP" ? "2016" : year; }

const std::vector<double> kEtaEdges = {-2.3, -1.6, -0.8, 0.0, 0.8, 1.6, 2.3};
const std::vector<double> kPtEdgesFine = {20, 25, 30, 35, 40, 50, 60, 80, 100, 150, 200};  // for the 1D vsJet curve
const std::vector<double> kPtEdges = {20, 30, 40, 60, 100, 150, 200};  // for the 2D tauES x decay-mode grid
// CMS tau decay-mode codes: 0=1-prong, 1=1-prong+pi0, 10=3-prong, 11=3-prong+pi0
const std::vector<int> kDecayModes = {0, 1, 10, 11};
}  // namespace

void PlotTauSF() {
    gStyle->SetOptStat(0);
    for (const auto& year : kYears) {
        std::cout << "=== " << year << " ===" << std::endl;
        const std::string fileName = kJsonDir + "/tau_" + eraFileTag(year) + ".json.gz";
        std::unique_ptr<correction::CorrectionSet> set;
        try { set = correction::CorrectionSet::from_file(fileName); }
        catch (const std::exception& e) { std::cerr << "ERROR loading " << fileName << ": " << e.what() << std::endl; continue; }
        auto vsE = set->at("DeepTau2017v2p1VSe");
        auto vsMu = set->at("DeepTau2017v2p1VSmu");
        auto vsJet = set->at("DeepTau2017v2p1VSjet");
        auto es = set->at("tau_energy_scale");

        // --- vsEle, vsMu: 1D vs eta ---
        for (auto& kv : std::vector<std::pair<std::string, correction::Correction::Ref>>{{"tauVsEle", vsE}, {"tauVsMu", vsMu}}) {
            const std::string& label = kv.first;
            const std::string wp = (label == "tauVsEle") ? "VVLoose" : "Medium";
            TH1D* hNom = new TH1D(("hNom_" + label + "_" + year).c_str(), "", kEtaEdges.size() - 1, kEtaEdges.data());
            TH1D* hUp = new TH1D(("hUp_" + label + "_" + year).c_str(), "", kEtaEdges.size() - 1, kEtaEdges.data());
            TH1D* hDown = new TH1D(("hDown_" + label + "_" + year).c_str(), "", kEtaEdges.size() - 1, kEtaEdges.data());
            hNom->SetDirectory(nullptr); hUp->SetDirectory(nullptr); hDown->SetDirectory(nullptr);
            auto& corr = kv.second;
            for (int b = 1; b <= hNom->GetNbinsX(); ++b) {
                const double eta = hNom->GetXaxis()->GetBinCenter(b);
                try {
                    const double nom = corr->evaluate({eta, kGenMatch, wp, "nom"});
                    const double up = corr->evaluate({eta, kGenMatch, wp, "up"});
                    const double down = corr->evaluate({eta, kGenMatch, wp, "down"});
                    hNom->SetBinContent(b, nom); hUp->SetBinContent(b, up); hDown->SetBinContent(b, down);
                } catch (const std::exception&) {}
            }
            zeroBinErrors(hNom);  // deterministic JSON lookup, not a statistical sample -- see zeroBinErrors' own comment
            drawSF1D(hNom, hUp, hDown, label + "_" + year, "#eta", "Scale factor", year,
                      kOutDir + "/" + year + "/" + label + ".png", 0.5, 1.5);
            delete hNom; delete hUp; delete hDown;
        }

        // --- vsJet ID SF: 1D vs pt only. The "pt"-flag branch this
        // framework uses genuinely does not consult decay mode at all
        // (verified directly against the JSON, see
        // LeptonTauSFUncertainty.h's own comment on this) -- plotting it
        // as a 2D grid vs decay mode would just repeat the same row 4
        // times, which is misleading rather than informative.
        {
            TH1D* hNom = new TH1D(("hNom_tauVsJet_" + year).c_str(), "", kPtEdgesFine.size() - 1, kPtEdgesFine.data());
            TH1D* hUp = new TH1D(("hUp_tauVsJet_" + year).c_str(), "", kPtEdgesFine.size() - 1, kPtEdgesFine.data());
            TH1D* hDown = new TH1D(("hDown_tauVsJet_" + year).c_str(), "", kPtEdgesFine.size() - 1, kPtEdgesFine.data());
            hNom->SetDirectory(nullptr); hUp->SetDirectory(nullptr); hDown->SetDirectory(nullptr);
            for (int b = 1; b <= hNom->GetNbinsX(); ++b) {
                const double pt = hNom->GetXaxis()->GetBinCenter(b);
                try {
                    const double nom = vsJet->evaluate({pt, 0, kGenMatch, "Medium", "VVLoose", "nom", "pt"});
                    const double up = vsJet->evaluate({pt, 0, kGenMatch, "Medium", "VVLoose", "up", "pt"});
                    const double down = vsJet->evaluate({pt, 0, kGenMatch, "Medium", "VVLoose", "down", "pt"});
                    hNom->SetBinContent(b, nom); hUp->SetBinContent(b, up); hDown->SetBinContent(b, down);
                } catch (const std::exception&) {}
            }
            zeroBinErrors(hNom);  // deterministic JSON lookup, not a statistical sample -- see zeroBinErrors' own comment
            drawSF1D(hNom, hUp, hDown, "tauVsJet_" + year, "p_{T} [GeV]", "Scale factor", year,
                      kOutDir + "/" + year + "/tauVsJet.png", 0.5, 1.5);
            delete hNom; delete hUp; delete hDown;
        }

        // --- Tau energy scale: 2D vs (pt, decay mode) ---
        {
            TH2D* hSF = new TH2D(("hSF_tauES_" + year).c_str(), "", kPtEdges.size() - 1, kPtEdges.data(), kDecayModes.size(), 0, kDecayModes.size());
            TH2D* hUnc = new TH2D(("hUnc_tauES_" + year).c_str(), "", kPtEdges.size() - 1, kPtEdges.data(), kDecayModes.size(), 0, kDecayModes.size());
            hSF->SetDirectory(nullptr); hUnc->SetDirectory(nullptr);
            for (size_t iy = 0; iy < kDecayModes.size(); ++iy) {
                hSF->GetYaxis()->SetBinLabel(iy + 1, Form("DM%d", kDecayModes[iy]));
                hUnc->GetYaxis()->SetBinLabel(iy + 1, Form("DM%d", kDecayModes[iy]));
            }
            for (int bx = 1; bx <= hSF->GetNbinsX(); ++bx) {
                const double pt = hSF->GetXaxis()->GetBinCenter(bx);
                for (size_t iy = 0; iy < kDecayModes.size(); ++iy) {
                    try {
                        const double nom = es->evaluate({pt, 0.0, kDecayModes[iy], kGenMatch, "DeepTau2017v2p1", "nom"});
                        const double up = es->evaluate({pt, 0.0, kDecayModes[iy], kGenMatch, "DeepTau2017v2p1", "up"});
                        const double down = es->evaluate({pt, 0.0, kDecayModes[iy], kGenMatch, "DeepTau2017v2p1", "down"});
                        hSF->SetBinContent(bx, iy + 1, nom);
                        const double relUnc = nom > 0.0 ? 0.5 * (std::fabs(up - nom) + std::fabs(nom - down)) / nom * 100.0 : 0.0;
                        hUnc->SetBinContent(bx, iy + 1, relUnc);
                    } catch (const std::exception&) {}
                }
            }
            drawSFAndUncertainty2D(hSF, hUnc, "tauES_" + year, "p_{T} [GeV]", "Decay mode", year,
                                    kOutDir + "/" + year + "/tauES.png", 0.9, 1.1, 5.0, true);
            delete hSF; delete hUnc;
        }
    }
    std::cout << "Done." << std::endl;
}
