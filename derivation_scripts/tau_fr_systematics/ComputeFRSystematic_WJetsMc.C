#include "FakeRateSources.h"
#include "SystematicUtils.h"

#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

const std::string kOutputBase = "../../Dependencies/systematics/results/wjets_mc";

bool processOne(const std::string& period, const std::string& region, TFile& out) {
    auto data = loadDataRate(period, region);
    auto wjMc = loadWJetsMCRate(period, region);
    if (!data || !wjMc) { std::cerr << "Skipping " << period << " " << region << " (missing DY data or W+jets MC input)" << std::endl; return false; }
    const std::string base = "FR_DY_" + period + "_" + region;
    std::unique_ptr<TH1D> up, down;
    makeSymmetricVariation(*data, *wjMc, base + "_up", base + "_down", up, down);
    out.cd();
    std::unique_ptr<TH1D> nominal(dynamic_cast<TH1D*>(data->Clone(base.c_str())));
    nominal->SetDirectory(nullptr);
    nominal->Write();
    up->Write();
    down->Write();
    return true;
}

}

void ComputeFRSystematic_WJetsMc() {
    gROOT->SetBatch(kTRUE);
    std::vector<std::string> periods = kYears;
    periods.push_back("Run2");
    int written = 0;
    for (const std::string& period : periods) {
        const std::string outDir = kOutputBase + "/" + period;
        gSystem->mkdir(outDir.c_str(), kTRUE);
        const std::string outPath = outDir + "/tau_fake_rate_syst_wjets_mc_" + period + ".root";
        std::unique_ptr<TFile> out(TFile::Open(outPath.c_str(), "RECREATE"));
        if (!out || out->IsZombie()) { std::cerr << "ERROR: cannot create " << outPath << std::endl; continue; }
        for (const std::string& region : {"barrel", "endcap"}) if (processOne(period, region, *out)) ++written;
        out->Close();
    }
    std::cout << "Wrote " << written << " of " << (periods.size() * 2) << " (period, region) W+jets-MC systematic variations under " << kOutputBase << std::endl;
}
