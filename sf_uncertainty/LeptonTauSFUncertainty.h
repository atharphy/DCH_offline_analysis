#pragma once
//
// Electron/muon/tau scale-factor UNCERTAINTY lookup, on top of the already-
// baked-in nominal SF branches (IDSF_*, ISOSF_*, TrigSF_*, TauVsEleIDSF_*,
// TauVsMuIDSF_*, TauVsJetIDSF_*) that the skims (Online_framework) write.
// Those branches carry no matching uncertainty, so this reads the SAME
// correctionlib JSON files Online_framework/funcs/outTuple.py uses to
// derive the nominal values, and evaluates the "up"/"down" (or "sfup"/
// "sfdown", "systup"/"systdown" -- naming differs per POG) systematic
// variants directly, returning a RELATIVE uncertainty:
//   relUp   = sf(up)/sf(nom)   - 1
//   relDown = 1 - sf(down)/sf(nom)
// which the caller applies to the ntuple's stored nominal value, e.g.
//   weightUp   = weightNominal * (1.0 + relUp)
//   weightDown = weightNominal * (1.0 - relDown)
// exactly like every other systematic source already applied in
// DCH_tauFR.C/DCH_tauFR_AltF.C (fake-rate, Z pT, pileup, ...). This keeps
// the result correct even though our own recomputed "nominal" (needed as
// the denominator for the ratio) may differ at the last few digits from
// the ntuple's stored nominal -- only the ratio matters.
//
// Source of truth for the exact correction names / input signatures /
// syst-string spelling below: Online_framework/funcs/outTuple.py
// (getIDISOTrigSF, getTauIDSFs) and the JSON files it loads from
// Online_framework/tools/{electron,muon_Z,tau}_<year>.json.gz. This
// header is self-contained -- offline/sf_uncertainty/json/ already holds
// its own copies of those same 12 files (4 years x electron/muon_Z/tau),
// and offline/sf_uncertainty/root/ holds the 4 electron-trigger ROOT
// files, so nothing here depends on Online_framework/ being present at
// runtime.
//
// Electron trigger SF is the one exception to "same source as
// Online_framework": Online_framework computes it from TGraphAsymmErrors
// efficiency curves (funcs/ScaleFactor.py) and never derives an
// uncertainty for it at all. ElectronTrigSFUncReader below derives one
// directly from those same graphs' own per-point asymmetric errors (see
// its class comment) -- this is method not present anywhere in
// Online_framework, not a translation of existing code.
//
// Known gap: the DeepTau2017v2p1VSjet "up"/"down" total-uncertainty keys
// (flag="pt" branch, what this header uses) do NOT depend on decay mode --
// verified directly against the JSON: the "dm" input is only consulted on
// the flag="dm" branch, which this header never selects. So a placeholder
// dm=0 is safe there. tau_energy_scale genuinely DOES depend on decay
// mode, though, so its uncertainty is only meaningful when a real decay
// mode is available (currently only for loose-vector-sourced tau
// candidates, via the `decayMode` branch in include/MyBranch.C -- cat-tight
// candidate slots have no per-slot decay-mode branch in the ntuple today).

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "TAxis.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"

#include "correction.h"

struct SFUnc {
    double up = 0.0;    // relative: sf(up)/sf(nom) - 1
    double down = 0.0;  // relative: 1 - sf(down)/sf(nom)
};

namespace SFUncDetail {

inline SFUnc fromVariants(double nom, double up, double down) {
    SFUnc result;
    if (!(nom > 0.0) || !std::isfinite(nom) || !std::isfinite(up) || !std::isfinite(down)) return result;
    result.up = up / nom - 1.0;
    result.down = 1.0 - down / nom;
    return result;
}

// electron_<year>.json.gz / muon_Z_<year>.json.gz filenames use this
// mapping (verified against Online_framework/tools/): the postVFP file is
// literally named "..._2016.json.gz", not "..._2016postVFP...".
inline std::string eraFileTag(const std::string& year) {
    if (year == "2016postVFP") return "2016";
    return year; // "2016preVFP", "2017", "2018" match directly
}

// String actually passed as the JSON's "year" category key -- this is
// DIFFERENT from the filename tag above for postVFP (see
// Online_framework/funcs/outTuple.py's `yearin` derivation).
inline std::string electronYearKey(const std::string& year) { return year; } // "2016preVFP"/"2016postVFP"/"2017"/"2018"
inline std::string muonYearKey(const std::string& year) { return year + "_UL"; }

}

// ---------------------------------------------------------------------
// Electron: single correction "UL-Electron-ID-SF" in electron_<year>.json.gz
// evaluate(year, ValType, WorkingPoint, eta, pt); ValType: sf/sfup/sfdown.
// ---------------------------------------------------------------------
class ElectronSFUncReader {
public:
    ElectronSFUncReader(const std::string& jsonDir, const std::string& year) : year_(year) {
        const std::string fileName = jsonDir + "/electron_" + SFUncDetail::eraFileTag(year) + ".json.gz";
        try {
            set_ = correction::CorrectionSet::from_file(fileName);
            idSf_ = set_->at("UL-Electron-ID-SF");
            std::cout << "[eSF] loaded UL-Electron-ID-SF from " << fileName << std::endl;
        }
        catch (const std::exception& error) {
            std::cerr << "ERROR: cannot load electron SF from " << fileName << ": " << error.what() << std::endl;
        }
    }

    // Reco SF uncertainty (only meaningful/applied for pt > 20 upstream, matching Online_framework).
    SFUnc getRecoUnc(double pt, double eta) const { return getWpUnc(pt, eta, "RecoAbove20"); }
    // ID+iso SF uncertainty. wp: "wp90noiso" (single-lepton legs) or "wp90iso" (dilepton legs), matching outTuple.py/outTuple2Lep.py.
    SFUnc getIdIsoUnc(double pt, double eta, const std::string& wp = "wp90noiso") const { return getWpUnc(pt, eta, wp); }

private:
    SFUnc getWpUnc(double pt, double eta, const std::string& wp) const {
        if (!idSf_) return SFUnc();
        const std::string yk = SFUncDetail::electronYearKey(year_);
        try {
            const double nom = idSf_->evaluate({yk, "sf", wp, eta, pt});
            const double up = idSf_->evaluate({yk, "sfup", wp, eta, pt});
            const double down = idSf_->evaluate({yk, "sfdown", wp, eta, pt});
            return SFUncDetail::fromVariants(nom, up, down);
        }
        catch (const std::exception& error) {
            std::cerr << "  [eSF-diag] evaluate() failed for wp=" << wp << " pt=" << pt << " eta=" << eta << ": " << error.what() << std::endl;
            return SFUnc();
        }
    }

    std::string year_;
    std::unique_ptr<correction::CorrectionSet> set_;
    correction::Correction::Ref idSf_;
};

// ---------------------------------------------------------------------
// Electron TRIGGER: not correctionlib-based in Online_framework at all --
// it's a data/MC efficiency ratio read off TGraphAsymmErrors curves in
// Electron_RunUL<year>_*.root (see funcs/ScaleFactor.py's SFs class).
// Online_framework never propagates an uncertainty for this SF (it only
// ever calls get_EfficiencyData/get_EfficiencyMC and divides), but the
// TGraphAsymmErrors themselves DO carry real per-point asymmetric y errors
// -- that's the whole reason that class is used instead of plain TGraph.
// This reader extracts those errors directly and propagates them through
// the ratio SF = effData/effMC via standard error propagation:
//   relErrHigh = sqrt((errDataHigh/effData)^2 + (errMCLow/effMC)^2)
//   relErrLow  = sqrt((errDataLow /effData)^2 + (errMCHigh/effMC)^2)
// Bin lookup (eta label, then pT bin within that eta's graph) exactly
// replicates ScaleFactor.py's SetAxisBins/FindPtBin/FindEtaLabel logic.
// ---------------------------------------------------------------------
class ElectronTrigSFUncReader {
public:
    ElectronTrigSFUncReader(const std::string& rootDir, const std::string& year) {
        std::string fileName;
        if (year == "2016preVFP") fileName = "Electron_RunUL2016preVFP_Ele25_EtaLt2p1.root";
        else if (year == "2016postVFP") fileName = "Electron_RunUL2016postVFP_Ele25_EtaLt2p1.root";
        else if (year == "2017") fileName = "Electron_RunUL2017_Ele35.root";
        else fileName = "Electron_RunUL2018_Ele35.root"; // 2018 (and default, matching outTuple.py)
        const std::string path = rootDir + "/" + fileName;

        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: cannot open electron trigger SF file: " << path << std::endl;
            if (f) f->Close();
            return;
        }
        TH1* h0 = (TH1*)f->Get("etaBinsH");
        if (!h0) {
            std::cerr << "ERROR: missing etaBinsH in " << path << std::endl;
            f->Close();
            return;
        }
        etaBinsH_.reset((TH1*)h0->Clone("etaBinsH_clone"));
        etaBinsH_->SetDirectory(nullptr);

        for (int ib = 1; ib <= etaBinsH_->GetNbinsX(); ++ib) {
            std::string rawLabel = etaBinsH_->GetXaxis()->GetBinLabel(ib);
            std::string label = (rawLabel.rfind("Eta", 0) == 0) ? rawLabel.substr(3) : rawLabel; // strip leading "Eta"
            TGraphAsymmErrors* gData = (TGraphAsymmErrors*)f->Get(("ZMassEta" + label + "_Data").c_str());
            TGraphAsymmErrors* gMC = (TGraphAsymmErrors*)f->Get(("ZMassEta" + label + "_MC").c_str());
            if (!gData || !gMC) {
                std::cerr << "  [eTrigSF-diag] missing graph for eta label '" << label << "' in " << path << std::endl;
                continue;
            }
            auto dClone = std::shared_ptr<TGraphAsymmErrors>((TGraphAsymmErrors*)gData->Clone());
            auto mClone = std::shared_ptr<TGraphAsymmErrors>((TGraphAsymmErrors*)gMC->Clone());
            setAxisBins(dClone.get());
            setAxisBins(mClone.get());
            dataGraphs_[label] = dClone;
            mcGraphs_[label] = mClone;
        }
        f->Close();
        std::cout << "[eTrigSF] loaded electron trigger efficiency graphs from " << path << std::endl;
    }

    SFUnc getTrigUnc(double pt, double eta) const {
        if (!etaBinsH_) return SFUnc();
        const int binNumber = etaBinsH_->GetXaxis()->FindFixBin(std::fabs(eta));
        std::string rawLabel = etaBinsH_->GetXaxis()->GetBinLabel(binNumber);
        const std::string label = (rawLabel.rfind("Eta", 0) == 0) ? rawLabel.substr(3) : rawLabel;

        auto itD = dataGraphs_.find(label);
        auto itM = mcGraphs_.find(label);
        if (itD == dataGraphs_.end() || itM == mcGraphs_.end()) return SFUnc();

        double effData, errDataLow, errDataHigh;
        double effMC, errMCLow, errMCHigh;
        if (!pointAt(itD->second.get(), pt, effData, errDataLow, errDataHigh)) return SFUnc();
        if (!pointAt(itM->second.get(), pt, effMC, errMCLow, errMCHigh)) return SFUnc();
        if (!(effData > 0.0) || !(effMC > 0.0)) return SFUnc();

        SFUnc result;
        result.up = std::sqrt(std::pow(errDataHigh / effData, 2) + std::pow(errMCLow / effMC, 2));
        result.down = std::sqrt(std::pow(errDataLow / effData, 2) + std::pow(errMCHigh / effMC, 2));
        return result;
    }

private:
    static void setAxisBins(TGraphAsymmErrors* g) {
        const int n = g->GetN();
        if (n <= 0) return;
        std::vector<double> edges(n + 1);
        for (int i = 0; i < n; ++i) edges[i] = g->GetX()[i] - g->GetErrorXlow(i);
        edges[n] = g->GetX()[n - 1] + g->GetErrorXhigh(n - 1);
        g->GetXaxis()->Set(n, edges.data());
    }

    // Replicates FindPtBin + the eff>1/eff<0 sentinel clamps from
    // ScaleFactor.py, but returns false (caller treats as "no info, no
    // uncertainty") instead of silently substituting eff=1 or eff=-1.
    static bool pointAt(const TGraphAsymmErrors* g, double pt, double& eff, double& errLow, double& errHigh) {
        const int n = g->GetN();
        if (n <= 0) return false;
        const double ptMax = g->GetX()[n - 1] + g->GetErrorXhigh(n - 1);
        const double ptMin = g->GetX()[0] - g->GetErrorXlow(0);
        int bin;
        if (pt >= ptMax) bin = n;
        else if (pt < ptMin) return false;
        else bin = g->GetXaxis()->FindFixBin(pt);
        if (bin < 1 || bin > n) return false;
        eff = g->GetY()[bin - 1];
        errLow = g->GetErrorYlow(bin - 1);
        errHigh = g->GetErrorYhigh(bin - 1);
        if (!std::isfinite(eff) || eff <= 0.0 || eff > 1.0) return false;
        if (!std::isfinite(errLow) || errLow < 0.0) errLow = 0.0;
        if (!std::isfinite(errHigh) || errHigh < 0.0) errHigh = 0.0;
        return true;
    }

    std::unique_ptr<TH1> etaBinsH_;
    std::unordered_map<std::string, std::shared_ptr<TGraphAsymmErrors>> dataGraphs_;
    std::unordered_map<std::string, std::shared_ptr<TGraphAsymmErrors>> mcGraphs_;
};

// ---------------------------------------------------------------------
// Muon: several corrections in muon_Z_<year>.json.gz.
// evaluate(year_UL, abseta, pt, ValType); ValType: sf/systup/systdown.
// Trigger correction name is year-dependent (matches outTuple.py exactly).
// ---------------------------------------------------------------------
class MuonSFUncReader {
public:
    MuonSFUncReader(const std::string& jsonDir, const std::string& year) : year_(year) {
        const std::string fileName = jsonDir + "/muon_Z_" + SFUncDetail::eraFileTag(year) + ".json.gz";
        try {
            set_ = correction::CorrectionSet::from_file(fileName);
            id_ = set_->at("NUM_TightID_DEN_TrackerMuons");
            iso_ = set_->at("NUM_TightRelIso_DEN_TightIDandIPCut");
            track_ = set_->at("NUM_TrackerMuons_DEN_genTracks");
            const std::string trigName = (year == "2016preVFP" || year == "2016postVFP")
                ? "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"
                : (year == "2017") ? "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"
                : "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"; // 2018
            trig_ = set_->at(trigName);
            std::cout << "[muSF] loaded muon SF corrections from " << fileName << std::endl;
        }
        catch (const std::exception& error) {
            std::cerr << "ERROR: cannot load muon SF from " << fileName << ": " << error.what() << std::endl;
        }
    }

    SFUnc getIdUnc(double pt, double eta) const { return evalOne(id_, pt, eta); }

    // Combined iso*tracking uncertainty. Online_framework computes the
    // nominal ISOSF as the PRODUCT of these two corrections; the up/down
    // variants here are formed the same way (product of each correction's
    // own up/down), i.e. treating the two sources as fully correlated --
    // a simplification, not the only defensible choice, but it mirrors how
    // the nominal is actually built.
    SFUnc getIsoUnc(double pt, double eta) const {
        if (!iso_ || !track_) return SFUnc();
        const std::string yk = SFUncDetail::muonYearKey(year_);
        const double ae = std::fabs(eta);
        try {
            const double nom = iso_->evaluate({yk, ae, pt, "sf"}) * track_->evaluate({yk, ae, pt, "sf"});
            const double up = iso_->evaluate({yk, ae, pt, "systup"}) * track_->evaluate({yk, ae, pt, "systup"});
            const double down = iso_->evaluate({yk, ae, pt, "systdown"}) * track_->evaluate({yk, ae, pt, "systdown"});
            return SFUncDetail::fromVariants(nom, up, down);
        }
        catch (const std::exception& error) {
            std::cerr << "  [muSF-diag] iso evaluate() failed pt=" << pt << " eta=" << eta << ": " << error.what() << std::endl;
            return SFUnc();
        }
    }

    SFUnc getTrigUnc(double pt, double eta) const { return evalOne(trig_, pt, eta); }

private:
    SFUnc evalOne(const correction::Correction::Ref& corr, double pt, double eta) const {
        if (!corr) return SFUnc();
        const std::string yk = SFUncDetail::muonYearKey(year_);
        const double ae = std::fabs(eta);
        try {
            const double nom = corr->evaluate({yk, ae, pt, "sf"});
            const double up = corr->evaluate({yk, ae, pt, "systup"});
            const double down = corr->evaluate({yk, ae, pt, "systdown"});
            return SFUncDetail::fromVariants(nom, up, down);
        }
        catch (const std::exception& error) {
            std::cerr << "  [muSF-diag] evaluate() failed pt=" << pt << " eta=" << eta << ": " << error.what() << std::endl;
            return SFUnc();
        }
    }

    std::string year_;
    std::unique_ptr<correction::CorrectionSet> set_;
    correction::Correction::Ref id_, iso_, track_, trig_;
};

// ---------------------------------------------------------------------
// Tau: 4 corrections in tau_<year>.json.gz, all DeepTau2017v2p1-based.
// WPs are hardcoded to match Online_framework's getTauIDSFs exactly
// (VVLoose vs-e, Medium vs-mu, Medium vs-jet, "pt"-flag vs-jet branch).
// ---------------------------------------------------------------------
class TauSFUncReader {
public:
    TauSFUncReader(const std::string& jsonDir, const std::string& year) {
        const std::string fileName = jsonDir + "/tau_" + SFUncDetail::eraFileTag(year) + ".json.gz";
        try {
            set_ = correction::CorrectionSet::from_file(fileName);
            vsE_ = set_->at("DeepTau2017v2p1VSe");
            vsMu_ = set_->at("DeepTau2017v2p1VSmu");
            vsJet_ = set_->at("DeepTau2017v2p1VSjet");
            es_ = set_->at("tau_energy_scale");
            std::cout << "[tauSF] loaded tau ID/ES SF corrections from " << fileName << std::endl;
        }
        catch (const std::exception& error) {
            std::cerr << "ERROR: cannot load tau SF from " << fileName << ": " << error.what() << std::endl;
        }
    }

    SFUnc getVsEleUnc(double eta, int genmatch) const {
        if (!vsE_) return SFUnc();
        try {
            const double nom = vsE_->evaluate({eta, genmatch, "VVLoose", "nom"});
            const double up = vsE_->evaluate({eta, genmatch, "VVLoose", "up"});
            const double down = vsE_->evaluate({eta, genmatch, "VVLoose", "down"});
            return SFUncDetail::fromVariants(nom, up, down);
        }
        catch (const std::exception& error) {
            std::cerr << "  [tauSF-diag] vsE evaluate() failed eta=" << eta << " gm=" << genmatch << ": " << error.what() << std::endl;
            return SFUnc();
        }
    }

    SFUnc getVsMuUnc(double eta, int genmatch) const {
        if (!vsMu_) return SFUnc();
        try {
            const double nom = vsMu_->evaluate({eta, genmatch, "Medium", "nom"});
            const double up = vsMu_->evaluate({eta, genmatch, "Medium", "up"});
            const double down = vsMu_->evaluate({eta, genmatch, "Medium", "down"});
            return SFUncDetail::fromVariants(nom, up, down);
        }
        catch (const std::exception& error) {
            std::cerr << "  [tauSF-diag] vsMu evaluate() failed eta=" << eta << " gm=" << genmatch << ": " << error.what() << std::endl;
            return SFUnc();
        }
    }

    // dm: pass a real decay mode when available (loose-vector taus, via the
    // `decayMode` branch); a placeholder (e.g. 0) is safe here specifically
    // because the "pt"-flag branch this calls never consults dm (verified
    // against the JSON's node tree -- dm only matters on the "dm"-flag
    // branch, which is not used by Online_framework's nominal computation
    // either, so this matches it exactly).
    SFUnc getVsJetUnc(double pt, int dm, int genmatch) const {
        if (!vsJet_) return SFUnc();
        try {
            const double nom = vsJet_->evaluate({pt, dm, genmatch, "Medium", "VVLoose", "nom", "pt"});
            const double up = vsJet_->evaluate({pt, dm, genmatch, "Medium", "VVLoose", "up", "pt"});
            const double down = vsJet_->evaluate({pt, dm, genmatch, "Medium", "VVLoose", "down", "pt"});
            return SFUncDetail::fromVariants(nom, up, down);
        }
        catch (const std::exception& error) {
            std::cerr << "  [tauSF-diag] vsJet evaluate() failed pt=" << pt << " dm=" << dm << " gm=" << genmatch << ": " << error.what() << std::endl;
            return SFUnc();
        }
    }

    // dm here genuinely matters (tau energy scale depends on decay mode) --
    // only call this with a real decay mode, not a placeholder.
    SFUnc getEnergyScaleUnc(double pt, double eta, int dm, int genmatch) const {
        if (!es_) return SFUnc();
        try {
            const double nom = es_->evaluate({pt, eta, dm, genmatch, "DeepTau2017v2p1", "nom"});
            const double up = es_->evaluate({pt, eta, dm, genmatch, "DeepTau2017v2p1", "up"});
            const double down = es_->evaluate({pt, eta, dm, genmatch, "DeepTau2017v2p1", "down"});
            return SFUncDetail::fromVariants(nom, up, down);
        }
        catch (const std::exception& error) {
            std::cerr << "  [tauSF-diag] ES evaluate() failed pt=" << pt << " eta=" << eta << " dm=" << dm << " gm=" << genmatch << ": " << error.what() << std::endl;
            return SFUnc();
        }
    }

private:
    std::unique_ptr<correction::CorrectionSet> set_;
    correction::Correction::Ref vsE_, vsMu_, vsJet_, es_;
};

// ---------------------------------------------------------------------
// Per-year caching, matching the getCachedTauFRReader(...) convention
// used elsewhere in DCH_modules.
// ---------------------------------------------------------------------
inline std::shared_ptr<ElectronSFUncReader> getCachedElectronSFUncReader(const std::string& jsonDir, const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<ElectronSFUncReader>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;
    auto reader = std::make_shared<ElectronSFUncReader>(jsonDir, year);
    cache[year] = reader;
    return reader;
}

inline std::shared_ptr<ElectronTrigSFUncReader> getCachedElectronTrigSFUncReader(const std::string& rootDir, const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<ElectronTrigSFUncReader>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;
    auto reader = std::make_shared<ElectronTrigSFUncReader>(rootDir, year);
    cache[year] = reader;
    return reader;
}

inline std::shared_ptr<MuonSFUncReader> getCachedMuonSFUncReader(const std::string& jsonDir, const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<MuonSFUncReader>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;
    auto reader = std::make_shared<MuonSFUncReader>(jsonDir, year);
    cache[year] = reader;
    return reader;
}

inline std::shared_ptr<TauSFUncReader> getCachedTauSFUncReader(const std::string& jsonDir, const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<TauSFUncReader>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;
    auto reader = std::make_shared<TauSFUncReader>(jsonDir, year);
    cache[year] = reader;
    return reader;
}
