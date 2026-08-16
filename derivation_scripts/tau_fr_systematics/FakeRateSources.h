#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"

namespace {

const std::string kFRRatesBase = "../../Dependencies/fake_rates/tau_fake_rates/rates";

const std::vector<std::string> kYears = {"2016preVFP", "2016postVFP", "2017", "2018"};

std::unique_ptr<TH1D> loadConsolidatedRate(const std::string& fileName, const std::string& period, const std::string& region) {
    const std::string path = kFRRatesBase + "/" + fileName;
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie()) { std::cerr << "ERROR: cannot open " << path << std::endl; return nullptr; }
    TH1D* source = dynamic_cast<TH1D*>(file->Get(("FR_" + period + "_" + region).c_str()));
    if (!source) { std::cerr << "ERROR: missing FR_" << period << "_" << region << " in " << path << std::endl; return nullptr; }
    std::unique_ptr<TH1D> clone(dynamic_cast<TH1D*>(source->Clone((fileName + "_" + period + "_" + region).c_str())));
    clone->SetDirectory(nullptr);
    return clone;
}

std::unique_ptr<TH1D> loadDataRate(const std::string& period, const std::string& region) { return loadConsolidatedRate("DY_data_tau_fake_rate.root", period, region); }
std::unique_ptr<TH1D> loadDYMCRate(const std::string& period, const std::string& region) { return loadConsolidatedRate("DY_MC_tau_fake_rate.root", period, region); }
std::unique_ptr<TH1D> loadWJetsMCRate(const std::string& period, const std::string& region) { return loadConsolidatedRate("WJ_MC_tau_fake_rate.root", period, region); }
std::unique_ptr<TH1D> loadQCDMCRate(const std::string& period, const std::string& region) { return loadConsolidatedRate("QCD_MC_tau_fake_rate.root", period, region); }
std::unique_ptr<TH1D> loadTTMCRate(const std::string& period, const std::string& region) { return loadConsolidatedRate("TT_MC_tau_fake_rate.root", period, region); }
std::unique_ptr<TH1D> loadWJetsDataRate(const std::string& period, const std::string& region) { return loadConsolidatedRate("WJ_data_tau_fake_rate.root", period, region); }

}
