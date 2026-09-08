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

class TauFakeSFReader {
public:
    TauFakeSFReader(const std::string& sfFile, const std::string& year) {
        TFile* f = TFile::Open(sfFile.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: cannot open tau fake-rate SF file: " << sfFile << std::endl;
            std::exit(1);
        }
        std::string hname = "Data_" + year + "_tau_fake_rate";
        TH2D* h0 = (TH2D*)f->Get(hname.c_str());
        if (!h0) {
            std::cerr << "ERROR: missing SF histogram: " << hname << std::endl;
            std::exit(1);
        }
        hSF = (TH2D*)h0->Clone((hname + "_sf_clone").c_str());
        hSF->SetDirectory(0);
        f->Close();
        std::cout << "[FakeSF] loaded " << hname << std::endl;
    }

    void get(double pt, double eta, double& value, double& error) const {
        value = 1.0;
        error = 0.0;
        if (!hSF) return;

        TAxis* xAxis = hSF->GetXaxis();
        TAxis* yAxis = hSF->GetYaxis();

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

        bx = std::max(1, std::min(bx, hSF->GetNbinsX()));
        by = std::max(1, std::min(by, hSF->GetNbinsY()));

        double sf = hSF->GetBinContent(bx, by);
        double err = hSF->GetBinError(bx, by);

        if (!std::isfinite(sf) || sf < 0.0) sf = 1.0;
        if (!std::isfinite(err) || err < 0.0) err = 0.0;

        value = sf;
        error = err;
    }

private:
    TH2D* hSF = nullptr;
};

inline std::shared_ptr<TauFakeSFReader> getCachedTauFakeSFReader(const std::string& sfFile, const std::string& year) {
    static std::unordered_map<std::string, std::shared_ptr<TauFakeSFReader>> cache;
    auto it = cache.find(year);
    if (it != cache.end()) return it->second;

    auto reader = std::make_shared<TauFakeSFReader>(sfFile, year);
    cache[year] = reader;
    return reader;
}
