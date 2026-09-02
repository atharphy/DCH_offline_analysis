// Stackhist.C
//
// Stacked-histogram plotter. Reads whatever "hist_*.root" files exist under
// INPUT_DIR (per-year subfolders, same layout the DCH scripts write to) and
// writes plots under OUTPUT_DIR. Produces one plot per (variable, channel,
// region), discovered directly from the histogram keys actually present.
//
// Run:
//   root -l -b -q 'Stackhist.C+("2018")'
//
// Edit INPUT_DIR / OUTPUT_DIR / DRAW_BANDS / USE_LOG_Y / EVENTS_PER_BIN_WIDTH
// below and recompile to change what is read, where plots go, or how they
// are drawn.
//
// Arguments:
//   inYear          = year, or "2016"/"Run2" for combined periods
//   firstVariable   = first variable index
//   nVariablesToRun = number of variables, -1 means all remaining variables

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TError.h>
#include <TGaxis.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "filemap/FileMap.h"

#include "Stack_modules/StackConfig.h"

static bool DRAW_BANDS = true;
static bool USE_LOG_Y = false;
static bool EVENTS_PER_BIN_WIDTH = false;
static bool blind_SR = true;
static std::string INPUT_DIR = "/eos/user/a/atahmad/DCH_offline_analysis/new_hists/run2_noFR_metphi_zpt_recoil_roccor_v2/";
static std::string OUTPUT_DIR = "/eos/user/a/atahmad/DCH_offline_analysis/new_plots/run2_noFR_metphi_zpt_recoil_roccor_v2/";

#include "Stack_modules/Labels.h"
#include "Stack_modules/HistCache.h"
#include "Stack_modules/SystematicSources.h"
#include "Stack_modules/HistMerge.h"
#include "Stack_modules/PlotUtils.h"
#include "Stack_modules/FileDiscovery.h"
#include "Stack_modules/StackDraw.h"

void Stackhist(std::string inYear = "2018", int firstVariable = 0, int nVariablesToRun = -1) {
    year = inYear;

    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);
    gStyle->SetGridColor(kGray + 2);
    gStyle->SetGridStyle(2);
    gStyle->SetGridWidth(1);
    TGaxis::SetExponentOffset(-0.07, 0.02, "y");

    const std::string histBase = ensureTrailingSlash(INPUT_DIR);
    const std::string filePrefix = "hist_";
    const std::string outdir = ensureTrailingSlash(OUTPUT_DIR) + year + "/";

    gSystem->mkdir(outdir.c_str(), kTRUE);

    std::map<std::string, int> fill_colors = {
        {"DY", 7}, {"DY10_50", 20}, {"WZ", 8}, {"WW", 40}, {"ZZ", 5}, {"TTbar", 46}, {"ttV", 4},
        {"WJ", 9}, {"ST", 38}, {"VVV", 6}, {"QCD", kOrange - 3}, {"other", 3},
        {"signal", 2}, {"data", kBlack}
    };

    std::map<std::string, std::vector<TFile*>> handles;
    const OpenFilesResult opened = openInputFiles(inYear, histBase, filePrefix, handles);

    std::cout << "\nInput path        : " << histBase << std::endl;
    std::cout << "Requested period  : " << year << std::endl;
    std::cout << "Output path       : " << outdir << std::endl;
    std::cout << "Opened ROOT files : " << opened.nOpened << std::endl;
    std::cout << "Missing files     : " << opened.nMissing << std::endl;
    std::cout << "Bad ROOT files    : " << opened.nBad << std::endl;

    if (opened.nOpened == 0) {
        std::cerr << "ERROR: no histogram ROOT files were opened." << std::endl;
        return;
    }

    const bool useLog = USE_LOG_Y;

    const int first = std::max(0, firstVariable);
    const int last = nVariablesToRun < 0 ? static_cast<int>(allVariables.size()) : std::min(static_cast<int>(allVariables.size()), first + nVariablesToRun);
    if (first >= last) {
        std::cerr << "ERROR: empty variable range [" << first << ", " << last << ")" << std::endl;
        return;
    }

    std::unordered_set<std::string> selectedVariables;
    for (int index = first; index < last; ++index) selectedVariables.insert(allVariables[index]);

    std::cout << "Pre-summing input histograms by process..." << std::endl;
    mergeInputFilesByProcess(handles, selectedVariables, kAllKnownSystSources, "merged_");
    std::cout << "Input pre-summing complete." << std::endl;

    std::set<std::string> plotNames;
    for (auto& processEntry : handles) {
        for (TFile* input : processEntry.second) {
            for (const std::string& storedName : keysInFile(input)) {
                if (storedName.rfind("h_", 0) != 0) continue;
                if (isSystematicVariantName(storedName, kAllKnownSystSources)) continue;
                const std::string hname = normalizeRegionName(storedName);
                if (!selectedVariables.count(getVariableFromHname(hname))) continue;
                plotNames.insert(hname);
            }
        }
    }

    std::cout << "Variables [" << first << ", " << last << "), existing nominal histograms to plot: " << plotNames.size() << std::endl;

    std::string previousVariable;
    for (const std::string& hname : plotNames) {
        const std::string var = getVariableFromHname(hname);
        if (var != previousVariable) {
            std::cout << "\n=== Starting variable: " << var << " ===" << std::endl;
            previousVariable = var;
        }
        const bool isSR = hname.find("_SR_") != std::string::npos;
        drawAndSave(hname, {hname}, handles, fill_colors, outdir, useLog, kAllSystSources, true, /*blind=*/isSR && blind_SR);
    }

    for (auto& processEntry : handles) {
        for (TFile* input : processEntry.second) {
            if (!input) continue;
            input->Close();
            delete input;
        }
        processEntry.second.clear();
    }

    handles.clear();
    fileKeys.clear();
    CleanUpROOTMemory();
}
