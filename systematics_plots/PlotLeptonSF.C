// Plots electron (Reco, ID+Iso) and muon (ID, Iso, Trigger) scale factors
// and their relative uncertainties as 2D (pt, eta) maps, one canvas per
// source per year (SF on the left, uncertainty in % on the right).
//
// Reads the same correctionlib JSON files already vendored for
// sf_uncertainty/ (self-contained, no CMSSW_BASE dependence beyond
// correctionlib itself). Evaluates correctionlib directly (same call
// signatures as sf_uncertainty/LeptonTauSFUncertainty.h) rather than
// reusing those reader classes, since this needs the nominal SF value too
// (those classes only expose the derived relative uncertainty).
//
// Run: root -l -b -q systematics_plots/PlotLeptonSF.C+
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

std::string eraFileTag(const std::string& year) { return year == "2016postVFP" ? "2016" : year; }
std::string muonYearKey(const std::string& year) { return year + "_UL"; }

// pt/eta grid: deliberately sparse so bin-value text stays legible, and
// eta bin CENTERS are chosen to avoid the ECAL barrel-endcap crack
// (|eta| in [1.444,1.566], where the electron SF is not meaningfully
// defined -- putting a bin center inside it previously produced a
// nonsensical "100% uncertainty" artifact from dividing by a near-zero SF).
const std::vector<double> kPtEdges = {15, 25, 40, 60, 100, 150, 200};
const std::vector<double> kEtaEdgesE = {-2.5, -2.0, -1.6, -0.8, 0.8, 1.6, 2.0, 2.5};
const std::vector<double> kEtaEdgesMu = {0.0, 0.9, 1.2, 2.1, 2.4};

TH2D* makeGrid(const std::string& name, const std::vector<double>& ptEdges, const std::vector<double>& etaEdges) {
    TH2D* h = new TH2D(name.c_str(), "", ptEdges.size() - 1, ptEdges.data(), etaEdges.size() - 1, etaEdges.data());
    h->SetDirectory(nullptr);
    return h;
}

void plotElectron(const std::string& year) {
    const std::string fileName = kJsonDir + "/electron_" + eraFileTag(year) + ".json.gz";
    std::unique_ptr<correction::CorrectionSet> set;
    try { set = correction::CorrectionSet::from_file(fileName); }
    catch (const std::exception& e) { std::cerr << "ERROR loading " << fileName << ": " << e.what() << std::endl; return; }
    auto idSf = set->at("UL-Electron-ID-SF");

    for (const auto& wpEntry : std::vector<std::pair<std::string,std::string>>{{"RecoAbove20", "eReco"}, {"wp90noiso", "eIdIso"}}) {
        const std::string& wp = wpEntry.first;
        const std::string& label = wpEntry.second;
        TH2D* hSF = makeGrid("hSF_" + label + "_" + year, kPtEdges, kEtaEdgesE);
        TH2D* hUnc = makeGrid("hUnc_" + label + "_" + year, kPtEdges, kEtaEdgesE);
        for (int bx = 1; bx <= hSF->GetNbinsX(); ++bx) {
            double pt = hSF->GetXaxis()->GetBinCenter(bx);
            if (wp == "RecoAbove20" && pt <= 20.0) pt = 20.5;  // RecoAbove20 map starts at pt=20
            for (int by = 1; by <= hSF->GetNbinsY(); ++by) {
                double eta = hSF->GetYaxis()->GetBinCenter(by);
                try {
                    const double nom = idSf->evaluate({year, "sf", wp, eta, pt});
                    const double up = idSf->evaluate({year, "sfup", wp, eta, pt});
                    const double down = idSf->evaluate({year, "sfdown", wp, eta, pt});
                    hSF->SetBinContent(bx, by, nom);
                    const double relUnc = nom > 0.0 ? 0.5 * (std::fabs(up - nom) + std::fabs(nom - down)) / nom * 100.0 : 0.0;
                    hUnc->SetBinContent(bx, by, relUnc);
                } catch (const std::exception&) { hSF->SetBinContent(bx, by, 0.0); hUnc->SetBinContent(bx, by, 0.0); }
            }
        }
        drawSFAndUncertainty2D(hSF, hUnc, label + "_" + year, "p_{T} [GeV]", "#eta", year,
                                kOutDir + "/" + year + "/" + label + ".png", 0.85, 1.10, 8.0, true);
        delete hSF; delete hUnc;
    }
}

void plotMuon(const std::string& year) {
    const std::string fileName = kJsonDir + "/muon_Z_" + eraFileTag(year) + ".json.gz";
    std::unique_ptr<correction::CorrectionSet> set;
    try { set = correction::CorrectionSet::from_file(fileName); }
    catch (const std::exception& e) { std::cerr << "ERROR loading " << fileName << ": " << e.what() << std::endl; return; }
    auto id = set->at("NUM_TightID_DEN_TrackerMuons");
    auto iso = set->at("NUM_TightRelIso_DEN_TightIDandIPCut");
    auto track = set->at("NUM_TrackerMuons_DEN_genTracks");
    const std::string trigName = (year == "2016preVFP" || year == "2016postVFP")
        ? "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"
        : (year == "2017") ? "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"
        : "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";
    auto trig = set->at(trigName);
    const std::string yk = muonYearKey(year);

    struct Source { std::string label; correction::Correction::Ref corr1; correction::Correction::Ref corr2; };
    std::vector<Source> sources = {{"muId", id, nullptr}, {"muIso", iso, track}, {"muTrig", trig, nullptr}};

    for (auto& src : sources) {
        TH2D* hSF = makeGrid("hSF_" + src.label + "_" + year, kPtEdges, kEtaEdgesMu);
        TH2D* hUnc = makeGrid("hUnc_" + src.label + "_" + year, kPtEdges, kEtaEdgesMu);
        for (int bx = 1; bx <= hSF->GetNbinsX(); ++bx) {
            const double pt = hSF->GetXaxis()->GetBinCenter(bx);
            for (int by = 1; by <= hSF->GetNbinsY(); ++by) {
                const double ae = std::fabs(hSF->GetYaxis()->GetBinCenter(by));
                try {
                    double nom, up, down;
                    if (src.corr2) {
                        nom = src.corr1->evaluate({yk, ae, pt, "sf"}) * src.corr2->evaluate({yk, ae, pt, "sf"});
                        up = src.corr1->evaluate({yk, ae, pt, "systup"}) * src.corr2->evaluate({yk, ae, pt, "systup"});
                        down = src.corr1->evaluate({yk, ae, pt, "systdown"}) * src.corr2->evaluate({yk, ae, pt, "systdown"});
                    } else {
                        nom = src.corr1->evaluate({yk, ae, pt, "sf"});
                        up = src.corr1->evaluate({yk, ae, pt, "systup"});
                        down = src.corr1->evaluate({yk, ae, pt, "systdown"});
                    }
                    hSF->SetBinContent(bx, by, nom);
                    const double relUnc = nom > 0.0 ? 0.5 * (std::fabs(up - nom) + std::fabs(nom - down)) / nom * 100.0 : 0.0;
                    hUnc->SetBinContent(bx, by, relUnc);
                } catch (const std::exception&) { hSF->SetBinContent(bx, by, 0.0); hUnc->SetBinContent(bx, by, 0.0); }
            }
        }
        drawSFAndUncertainty2D(hSF, hUnc, src.label + "_" + year, "p_{T} [GeV]", "|#eta|", year,
                                kOutDir + "/" + year + "/" + src.label + ".png", 0.85, 1.10, 8.0, true);
        delete hSF; delete hUnc;
    }
}

}  // namespace

void PlotLeptonSF() {
    gStyle->SetOptStat(0);
    for (const auto& year : kYears) {
        std::cout << "=== " << year << " ===" << std::endl;
        plotElectron(year);
        plotMuon(year);
    }
    std::cout << "Done. Plots written under " << kOutDir << "/<year>/" << std::endl;
}
