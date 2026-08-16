#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

#include "TFile.h"
#include "TH2D.h"

class TauFRReader {
public:
    TauFRReader(const std::string& frFile, const std::string& year) {
        TFile* f = TFile::Open(frFile.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: cannot open tau FR file: " << frFile << std::endl;
            std::exit(1);
        }
        std::string hname = "Data_" + year + "_tau_fake_rate";
        TH2D* h0 = (TH2D*)f->Get(hname.c_str());
        if (!h0) {
            std::cerr << "ERROR: missing FR histogram: " << hname << std::endl;
            std::exit(1);
        }
        hFR = (TH2D*)h0->Clone((hname + "_clone").c_str());
        hFR->SetDirectory(0);
        f->Close();
        std::cout << "[FR] loaded " << hname << std::endl;
    }

    explicit TauFRReader(TH2D* preBuilt) {
        if (!preBuilt) return;
        hFR = (TH2D*)preBuilt->Clone((std::string(preBuilt->GetName()) + "_clone").c_str());
        hFR->SetDirectory(0);
    }

    void get(double pt, double eta, double& value, double& error) const {
        value = 0.0;
        error = 0.0;
        if (!hFR) return;

        TAxis* xAxis = hFR->GetXaxis();
        TAxis* yAxis = hFR->GetYaxis();

        const double xMin = xAxis->GetXmin();
        const double xMax = xAxis->GetXmax();
        const double yMin = yAxis->GetXmin();
        const double yMax = yAxis->GetXmax();

        const double xEps = 1e-3 * (xMax - xMin);
        const double yEps = 1e-3 * (yMax - yMin);

        double ptc = pt;
        if (ptc < xMin) ptc = xMin + xEps;
        if (ptc > xMax) ptc = xMax - xEps;

        double aeta = std::fabs(eta);
        if (aeta < yMin) aeta = yMin + yEps;
        if (aeta > yMax) aeta = yMax - yEps;

        int bx = xAxis->FindBin(ptc);
        int by = yAxis->FindBin(aeta);

        bx = std::max(1, std::min(bx, hFR->GetNbinsX()));
        by = std::max(1, std::min(by, hFR->GetNbinsY()));

        double f = hFR->GetBinContent(bx, by);
        double err = hFR->GetBinError(bx, by);

        if (!std::isfinite(f) || f < 0.0) f = 0.0;
        if (f > 0.999) f = 0.999;
        if (!std::isfinite(err) || err < 0.0) err = 0.0;

        value = f;
        error = err;
    }

private:
    TH2D* hFR = nullptr;
};

inline std::shared_ptr<TauFRReader> getCachedTauFRReader(const std::string& frFile, const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<TauFRReader>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;

    auto reader = std::make_shared<TauFRReader>(frFile, year);
    cache[year] = reader;
    return reader;
}
