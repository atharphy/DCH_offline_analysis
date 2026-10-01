#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "filemap/FileMap.h"

static bool DRAW_BANDS = false;
static bool USE_LOG_Y = true;
static bool blind_SR = true;
static std::string INPUT_DIR = "hists/run2_hists_tauFR_roccor/";
static std::string OUTPUT_DIR = "yield_plot/run2_plots_tauFR_roccor_updated_SF/";

#include "yield_modules/YieldConfig.h"
#include "yield_modules/YieldSystematics.h"
#include "yield_modules/YieldHistUtils.h"
#include "yield_modules/YieldFileIO.h"
#include "yield_modules/YieldDraw.h"

void YieldPlots(std::string year = "Run2") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    const std::string outputDirectory = ensureTrailingSlash(OUTPUT_DIR) + year;
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);

    std::map<std::string, std::vector<TFile*>> filesByProcess;
    const YieldOpenResult opened = openYieldInputFiles(year, INPUT_DIR, filesByProcess);

    std::cout << "Histogram folder : " << INPUT_DIR << std::endl;
    std::cout << "Output folder    : " << outputDirectory << std::endl;
    std::cout << "Opened files     : " << opened.nOpened << std::endl;
    std::cout << "Missing files    : " << opened.nMissing << std::endl;

    if (opened.nOpened == 0) {
        std::cerr << "ERROR: no input ROOT files were opened" << std::endl;
        return;
    }

    auto readSFFromCSV = [](const std::string& csvPath) {
        static std::map<std::string, std::map<std::string, double>> cache;
        auto& yearMap = cache[csvPath];
        if (yearMap.empty()) {
            std::ifstream fin(csvPath);
            std::string line;
            if (fin) {
                std::getline(fin, line);
                while (std::getline(fin, line)) {
                    size_t p1 = line.find(',');
                    size_t p2 = line.find(',', p1 == std::string::npos ? p1 : p1 + 1);
                    if (p1 == std::string::npos || p2 == std::string::npos) continue;
                    yearMap[line.substr(0, p1)] = std::stod(line.substr(p1 + 1, p2 - p1 - 1));
                }
            }
        }
        return yearMap;
    };
    auto scaleYieldGroup = [&](ProcessYields& py, const std::string& csvPath) {
        const auto yearMap = readSFFromCSV(csvPath);
        auto it = yearMap.find(year);
        if (it == yearMap.end()) {
            std::cerr << "[warning] no scale factor for year " << year << " in " << csvPath << ", leaving unscaled" << std::endl;
            return;
        }
        if (py.nominal) py.nominal->Scale(it->second);
        for (auto& u : py.up) u.second->Scale(it->second);
        for (auto& d : py.down) d.second->Scale(it->second);
    };

    for (const std::string& region : {"CR", "VR", "SR"}) {
        std::map<std::string, ProcessYields> yields;

        for (const std::string& process : kProcesses) {
            ProcessYields py;
            py.nominal = makeYieldHistogram("yield_" + process + "_" + region);
            auto found = filesByProcess.find(process);
            if (found != filesByProcess.end()) {
                fillProcessYields(*py.nominal, found->second, region);
                if (DRAW_BANDS) {
                    for (const auto& src : kAllSystSources) {
                        auto up = makeYieldHistogram("yield_" + process + "_" + region + "_" + src.name + "Up");
                        auto down = makeYieldHistogram("yield_" + process + "_" + region + "_" + src.name + "Down");
                        const bool foundUp = fillProcessYields(*up, found->second, region, src.upSuffix);
                        const bool foundDown = fillProcessYields(*down, found->second, region, src.downSuffix);
                        if (foundUp) py.up[src.name] = std::move(up);
                        if (foundDown) py.down[src.name] = std::move(down);
                    }
                }
            }
            yields[process] = std::move(py);
        }

        if (yields.count("WZ")) scaleYieldGroup(yields["WZ"], "wz_scale_factors.csv");
        if (yields.count("ZZ")) scaleYieldGroup(yields["ZZ"], "zz_scale_factors.csv");

        ProcessYields dataYields;
        dataYields.nominal = makeYieldHistogram("yield_data_" + region);
        auto dataFiles = filesByProcess.find("data");
        if (dataFiles != filesByProcess.end()) fillProcessYields(*dataYields.nominal, dataFiles->second, region);
        yields["data"] = std::move(dataYields);

        const bool blind = (region == "SR") && blind_SR;
        drawRegion(region, year, USE_LOG_Y, DRAW_BANDS, yields, outputDirectory, blind);
    }

    for (auto& process : filesByProcess) {
        for (TFile* file : process.second) {
            file->Close();
            delete file;
        }
    }

    std::cout << "Finished CR, VR, and SR yield plots" << std::endl;
}
