#pragma once

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
    double up = 0.0;
    double down = 0.0;
};

namespace SFUncDetail {

inline SFUnc fromVariants(double nom, double up, double down) {
    SFUnc result;
    if (!(nom > 0.0) || !std::isfinite(nom) || !std::isfinite(up) || !std::isfinite(down)) return result;
    result.up = up / nom - 1.0;
    result.down = 1.0 - down / nom;
    return result;
}

inline std::string eraFileTag(const std::string& year) {
    if (year == "2016postVFP") return "2016";
    return year;
}

inline std::string electronYearKey(const std::string& year) { return year; }
inline std::string muonYearKey(const std::string& year) { return year + "_UL"; }

}

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

    SFUnc getRecoUnc(double pt, double eta) const { return getWpUnc(pt, eta, "RecoAbove20"); }

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

class ElectronTrigSFUncReader {
public:
    ElectronTrigSFUncReader(const std::string& rootDir, const std::string& year) {
        std::string fileName;
        if (year == "2016preVFP") fileName = "Electron_RunUL2016preVFP_Ele25_EtaLt2p1.root";
        else if (year == "2016postVFP") fileName = "Electron_RunUL2016postVFP_Ele25_EtaLt2p1.root";
        else if (year == "2017") fileName = "Electron_RunUL2017_Ele35.root";
        else fileName = "Electron_RunUL2018_Ele35.root";
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
            std::string label = (rawLabel.rfind("Eta", 0) == 0) ? rawLabel.substr(3) : rawLabel;
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
                : "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";
            trig_ = set_->at(trigName);
            std::cout << "[muSF] loaded muon SF corrections from " << fileName << std::endl;
        }
        catch (const std::exception& error) {
            std::cerr << "ERROR: cannot load muon SF from " << fileName << ": " << error.what() << std::endl;
        }
    }

    SFUnc getIdUnc(double pt, double eta) const { return evalOne(id_, pt, eta); }

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
