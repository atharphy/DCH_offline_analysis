#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

const std::string kInputBase = "output";

const std::vector<std::string> kYears = {
    "2016preVFP", "2016postVFP", "2017", "2018"
};

struct Quantity {
    std::string component;
    std::string statistic;
    std::string yTitle;
    std::string label;
};

const std::vector<Quantity> kQuantities = {
    {"upar",  "mean",  "#mu(U_{1}) [GeV]",    "Parallel response"},
    {"upar",  "sigma", "#sigma(U_{1}) [GeV]", "Parallel resolution"},
    {"uperp", "mean",  "#mu(U_{2}) [GeV]",    "Perpendicular mean"},
    {"uperp", "sigma", "#sigma(U_{2}) [GeV]", "Perpendicular resolution"}
};

std::string periodLabel(const std::string& year) {
    if (year == "2016preVFP")  return "2016 preVFP (19.5 fb^{-1})";
    if (year == "2016postVFP") return "2016 postVFP (16.8 fb^{-1})";
    if (year == "2017")        return "2017 (41.5 fb^{-1})";
    if (year == "2018")        return "2018 (59.7 fb^{-1})";
    return year;
}

std::string jetLabel(int jet) {
    if (jet == 0) return "N_{jets} = 0";
    if (jet == 1) return "N_{jets} = 1";
    return "N_{jets} #geq 2";
}

std::unique_ptr<TH1D> makeLogAxisCopy(
    const TH1D* source,
    const std::string& name
) {
    const double edges[] = {1., 10., 20., 30., 50., 1000.};
    auto output = std::make_unique<TH1D>(name.c_str(), "", 5, edges);
    output->SetDirectory(nullptr);
    for (int bin = 1; bin <= 5; ++bin) {
        output->SetBinContent(bin, source->GetBinContent(bin));
        output->SetBinError(bin, source->GetBinError(bin));
    }
    return output;
}

void setRange(TH1D* data, TH1D* simulation, bool isSigma) {
    double minimum = 0.0;
    double maximum = 0.0;
    bool initialized = false;
    for (TH1D* histogram : {data, simulation}) {
        for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
            const double value = histogram->GetBinContent(bin);
            const double error = histogram->GetBinError(bin);
            if (!std::isfinite(value) || !std::isfinite(error)) continue;
            const double low = value - error;
            const double high = value + error;
            if (!initialized) {
                minimum = low;
                maximum = high;
                initialized = true;
            } else {
                minimum = std::min(minimum, low);
                maximum = std::max(maximum, high);
            }
        }
    }
    if (!initialized) {
        minimum = 0.0;
        maximum = 1.0;
    }
    const double span = std::max(maximum - minimum, 1.0);
    if (isSigma) {
        data->SetMinimum(0.0);
        data->SetMaximum(std::max(1.25 * maximum, 1.0));
    } else {
        data->SetMinimum(minimum - 0.25 * span);
        data->SetMaximum(maximum + 0.35 * span);
    }
}

bool drawPlot(
    TFile& input,
    const std::string& year,
    const Quantity& quantity,
    int jet,
    const std::string& outputDirectory
) {
    const std::string stem = quantity.component + "_" + quantity.statistic;
    const std::string dataName =
        "h_" + stem + "_data_njet" + std::to_string(jet);
    const std::string mcName =
        "h_" + stem + "_mc_njet" + std::to_string(jet);

    TH1D* sourceData = dynamic_cast<TH1D*>(input.Get(dataName.c_str()));
    TH1D* sourceMC = dynamic_cast<TH1D*>(input.Get(mcName.c_str()));
    if (!sourceData || !sourceMC) {
        std::cerr << "Missing " << dataName << " or " << mcName << std::endl;
        return false;
    }

    auto data = makeLogAxisCopy(sourceData, dataName + "_plot");
    auto simulation = makeLogAxisCopy(sourceMC, mcName + "_plot");

    data->SetLineColor(kBlack);
    data->SetMarkerColor(kBlack);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(1.15);
    data->SetLineWidth(2);
    simulation->SetLineColor(kRed + 1);
    simulation->SetMarkerColor(kRed + 1);
    simulation->SetMarkerStyle(24);
    simulation->SetMarkerSize(1.15);
    simulation->SetLineWidth(2);

    data->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
    data->GetYaxis()->SetTitle(quantity.yTitle.c_str());
    data->GetXaxis()->SetMoreLogLabels();
    data->GetXaxis()->SetNoExponent();
    data->GetXaxis()->SetTitleSize(0.050);
    data->GetYaxis()->SetTitleSize(0.050);
    data->GetXaxis()->SetLabelSize(0.043);
    data->GetYaxis()->SetLabelSize(0.043);
    data->GetYaxis()->SetTitleOffset(1.15);
    setRange(data.get(), simulation.get(), quantity.statistic == "sigma");

    TCanvas canvas("c_recoil", "", 1200, 800);
    canvas.SetLeftMargin(0.13);
    canvas.SetRightMargin(0.04);
    canvas.SetTopMargin(0.10);
    canvas.SetBottomMargin(0.13);
    canvas.SetLogx();
    canvas.SetGridx();
    canvas.SetGridy();
    canvas.SetTicks(1, 1);

    data->Draw("E1 P");
    simulation->Draw("E1 P SAME");

    TLegend legend(0.68, 0.72, 0.93, 0.87);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.040);
    legend.AddEntry(data.get(), "Data - non-DY", "lep");
    legend.AddEntry(simulation.get(), "DY simulation", "lep");
    legend.Draw();

    TLatex text;
    text.SetNDC();
    text.SetTextFont(62);
    text.SetTextSize(0.050);
    text.DrawLatex(0.13, 0.925, "CMS");
    text.SetTextFont(52);
    text.SetTextSize(0.040);
    text.DrawLatex(0.225, 0.925, "Preliminary");
    text.SetTextFont(42);
    text.SetTextAlign(31);
    text.DrawLatex(0.96, 0.925, periodLabel(year).c_str());
    text.SetTextAlign(11);
    text.SetTextSize(0.041);
    text.DrawLatex(0.17, 0.84, quantity.label.c_str());
    text.DrawLatex(0.17, 0.78, jetLabel(jet).c_str());

    const std::string outputStem = outputDirectory + "/" + stem +
        "_njet" + std::to_string(jet);
    canvas.SaveAs((outputStem + ".png").c_str());
    canvas.SaveAs((outputStem + ".pdf").c_str());
    return true;
}

void plotYear(const std::string& year) {
    const std::string inputName = kInputBase + "/" + year +
        "/RecoilCorrections_" + year + ".root";
    std::unique_ptr<TFile> input(TFile::Open(inputName.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "Cannot open " << inputName << std::endl;
        return;
    }

    const std::string outputDirectory =
        kInputBase + "/" + year + "/parameter_plots_logx";
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);
    int created = 0;
    for (const Quantity& quantity : kQuantities) {
        for (int jet = 0; jet < 3; ++jet) {
            if (drawPlot(*input, year, quantity, jet, outputDirectory)) ++created;
        }
    }
    std::cout << "[" << year << "] created " << created
              << " plots in " << outputDirectory << std::endl;
}

}

void PlotRecoilParameters(std::string period = "all") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(5);
    if (period == "all" || period == "Run2") {
        for (const std::string& year : kYears) plotYear(year);
        return;
    }
    if (std::find(kYears.begin(), kYears.end(), period) == kYears.end()) {
        std::cerr << "Unsupported period: " << period << std::endl;
        return;
    }
    plotYear(period);
}
