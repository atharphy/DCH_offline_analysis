#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
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

const std::string kInputBase =
    "../../Dependencies/zpt_weights";

const std::vector<std::string> kPeriods = {
    "2016preVFP", "2016postVFP", "2017", "2018", "Run2"
};

std::string periodLabel(const std::string& period) {
    if (period == "2016preVFP")  return "2016 preVFP (19.5 fb^{-1})";
    if (period == "2016postVFP") return "2016 postVFP (16.8 fb^{-1})";
    if (period == "2017")        return "2017 (41.5 fb^{-1})";
    if (period == "2018")        return "2018 (59.7 fb^{-1})";
    if (period == "Run2")        return "Run 2 (137.6 fb^{-1})";
    return period;
}

std::unique_ptr<TH1D> makeLogCopy(
    const TH1D* source,
    const std::string& name
) {
    std::vector<double> edges(source->GetNbinsX() + 1);
    for (int edge = 0; edge <= source->GetNbinsX(); ++edge)
        edges[edge] = source->GetXaxis()->GetBinLowEdge(edge + 1);
    edges.front() = 0.1;

    auto output = std::make_unique<TH1D>(
        name.c_str(), "", source->GetNbinsX(), edges.data()
    );
    output->SetDirectory(nullptr);
    for (int bin = 1; bin <= source->GetNbinsX(); ++bin) {
        output->SetBinContent(bin, source->GetBinContent(bin));
        output->SetBinError(bin, source->GetBinError(bin));
    }
    return output;
}

void styleAxes(TH1D* histogram, bool ratio) {
    histogram->GetXaxis()->SetMoreLogLabels();
    histogram->GetXaxis()->SetNoExponent();
    histogram->GetXaxis()->SetTitle("reconstructed p_{T}^{Z} [GeV]");
    if (!ratio) {
        histogram->GetYaxis()->SetTitle("Normalized events");
        histogram->GetXaxis()->SetLabelSize(0.0);
        histogram->GetXaxis()->SetTitleSize(0.0);
        histogram->GetYaxis()->SetTitleSize(0.060);
        histogram->GetYaxis()->SetLabelSize(0.050);
        histogram->GetYaxis()->SetTitleOffset(1.05);
    } else {
        histogram->GetYaxis()->SetTitle("Data / DY");
        histogram->GetXaxis()->SetTitleSize(0.125);
        histogram->GetXaxis()->SetLabelSize(0.105);
        histogram->GetXaxis()->SetTitleOffset(1.05);
        histogram->GetYaxis()->SetTitleSize(0.105);
        histogram->GetYaxis()->SetLabelSize(0.090);
        histogram->GetYaxis()->SetTitleOffset(0.55);
        histogram->GetYaxis()->SetNdivisions(505);
    }
}

bool drawChannel(
    TFile& input,
    const std::string& period,
    const std::string& channel,
    const std::string& outputDirectory
) {
    TH1D* sourceData = dynamic_cast<TH1D*>(
        input.Get(("h_data_shape_" + channel).c_str())
    );
    TH1D* sourceDY = dynamic_cast<TH1D*>(
        input.Get(("h_dy_shape_" + channel).c_str())
    );
    TH1D* sourceRatio = dynamic_cast<TH1D*>(
        input.Get(("h_zpt_weight_" + channel).c_str())
    );
    if (!sourceData || !sourceDY || !sourceRatio) {
        std::cerr << "Missing Z pT histograms for " << period << " "
                  << channel << std::endl;
        return false;
    }

    auto data = makeLogCopy(sourceData, "data_plot_" + period + channel);
    auto dy = makeLogCopy(sourceDY, "dy_plot_" + period + channel);
    auto ratio = makeLogCopy(sourceRatio, "ratio_plot_" + period + channel);

    data->SetMarkerStyle(20);
    data->SetMarkerSize(1.05);
    data->SetMarkerColor(kBlack);
    data->SetLineColor(kBlack);
    data->SetLineWidth(2);
    dy->SetMarkerStyle(24);
    dy->SetMarkerSize(1.05);
    dy->SetMarkerColor(kRed + 1);
    dy->SetLineColor(kRed + 1);
    dy->SetLineWidth(2);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(0.95);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetLineWidth(2);

    styleAxes(data.get(), false);
    styleAxes(ratio.get(), true);
    const double shapeMaximum = std::max(data->GetMaximum(), dy->GetMaximum());
    data->SetMinimum(0.0);
    data->SetMaximum(shapeMaximum > 0.0 ? 1.38 * shapeMaximum : 1.0);
    ratio->SetMinimum(0.0);
    ratio->SetMaximum(std::max(2.0, 1.25 * ratio->GetMaximum()));

    TCanvas canvas(("c_zpt_" + period + channel).c_str(), "", 1200, 900);
    TPad top(("top_" + period + channel).c_str(), "", 0.0, 0.30, 1.0, 1.0);
    TPad bottom(("bottom_" + period + channel).c_str(), "", 0.0, 0.0, 1.0, 0.30);
    top.SetLeftMargin(0.13);
    top.SetRightMargin(0.04);
    top.SetTopMargin(0.12);
    top.SetBottomMargin(0.025);
    bottom.SetLeftMargin(0.13);
    bottom.SetRightMargin(0.04);
    bottom.SetTopMargin(0.025);
    bottom.SetBottomMargin(0.34);
    top.SetLogx();
    bottom.SetLogx();
    top.SetGridx();
    top.SetGridy();
    bottom.SetGridx();
    bottom.SetGridy();
    top.SetTicks(1, 1);
    bottom.SetTicks(1, 1);
    top.Draw();
    bottom.Draw();

    top.cd();
    data->Draw("E1 P");
    dy->Draw("E1 P SAME");
    data->Draw("E1 P SAME");

    TLegend legend(0.66, 0.70, 0.93, 0.86);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.043);
    legend.AddEntry(data.get(), "Data - non-DY", "lep");
    legend.AddEntry(dy.get(), "DY simulation", "lep");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(62);
    label.SetTextSize(0.055);
    label.DrawLatex(0.13, 0.92, "CMS");
    label.SetTextFont(52);
    label.SetTextSize(0.043);
    label.DrawLatex(0.225, 0.92, "Preliminary");
    label.SetTextFont(42);
    label.SetTextAlign(31);
    label.DrawLatex(0.96, 0.92, periodLabel(period).c_str());
    label.SetTextAlign(11);
    label.SetTextSize(0.044);
    label.DrawLatex(
        0.17, 0.81,
        channel == "ee" ? "Z #rightarrow ee" : "Z #rightarrow #mu#mu"
    );

    bottom.cd();
    ratio->Draw("E1 P");
    TLine unity(0.1, 1.0, 500.0, 1.0);
    unity.SetLineStyle(2);
    unity.SetLineWidth(2);
    unity.Draw("SAME");

    const std::string outputStem = outputDirectory +
        "/ZPtWeight_" + channel + "_logx_ratio";
    canvas.SaveAs((outputStem + ".png").c_str());
    canvas.SaveAs((outputStem + ".pdf").c_str());
    return true;
}

void plotPeriod(const std::string& period) {
    const std::string inputName = kInputBase + "/" + period +
        "/ZPtWeights_" + period + ".root";
    std::unique_ptr<TFile> input(TFile::Open(inputName.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "Cannot open " << inputName << std::endl;
        return;
    }
    const std::string outputDirectory =
        kInputBase + "/" + period + "/formatted_plots_logx";
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);
    int created = 0;
    for (const std::string& channel : {std::string("ee"), std::string("mm")})
        if (drawChannel(*input, period, channel, outputDirectory)) ++created;
    std::cout << "[" << period << "] created " << created
              << " plots in " << outputDirectory << std::endl;
}

void drawYearComparison(const std::string& channel) {
    const std::vector<std::string> years = {
        "2016preVFP", "2016postVFP", "2017", "2018"
    };
    const int colors[] = {kAzure + 1, kGreen + 2, kOrange + 7, kRed + 1};
    const int markers[] = {20, 21, 22, 23};
    std::vector<std::unique_ptr<TFile>> files;
    std::vector<std::unique_ptr<TH1D>> weights;

    for (size_t index = 0; index < years.size(); ++index) {
        const std::string inputName = kInputBase + "/" + years[index] +
            "/ZPtWeights_" + years[index] + ".root";
        std::unique_ptr<TFile> file(TFile::Open(inputName.c_str(), "READ"));
        if (!file || file->IsZombie()) {
            std::cerr << "Cannot open " << inputName << std::endl;
            return;
        }
        TH1D* source = dynamic_cast<TH1D*>(
            file->Get(("h_zpt_weight_" + channel).c_str())
        );
        if (!source) {
            std::cerr << "Missing h_zpt_weight_" << channel << " in "
                      << inputName << std::endl;
            return;
        }
        auto weight = makeLogCopy(
            source, "weight_comparison_" + years[index] + channel
        );
        weight->SetLineColor(colors[index]);
        weight->SetMarkerColor(colors[index]);
        weight->SetMarkerStyle(markers[index]);
        weight->SetMarkerSize(1.05);
        weight->SetLineWidth(2);
        files.push_back(std::move(file));
        weights.push_back(std::move(weight));
    }

    double maximum = 0.0;
    for (const auto& weight : weights) {
        for (int bin = 1; bin <= weight->GetNbinsX(); ++bin)
            maximum = std::max(
                maximum,
                weight->GetBinContent(bin) + weight->GetBinError(bin)
            );
    }

    TH1D* axis = weights.front().get();
    axis->GetXaxis()->SetTitle("reconstructed p_{T}^{Z} [GeV]");
    axis->GetYaxis()->SetTitle("Z p_{T} correction factor");
    axis->GetXaxis()->SetMoreLogLabels();
    axis->GetXaxis()->SetNoExponent();
    axis->GetXaxis()->SetTitleSize(0.050);
    axis->GetYaxis()->SetTitleSize(0.050);
    axis->GetXaxis()->SetLabelSize(0.043);
    axis->GetYaxis()->SetLabelSize(0.043);
    axis->GetYaxis()->SetTitleOffset(1.15);
    axis->SetMinimum(0.0);
    axis->SetMaximum(std::max(2.0, 1.25 * maximum));

    TCanvas canvas(("c_year_comparison_" + channel).c_str(), "", 1200, 800);
    canvas.SetLeftMargin(0.13);
    canvas.SetRightMargin(0.04);
    canvas.SetTopMargin(0.10);
    canvas.SetBottomMargin(0.13);
    canvas.SetLogx();
    canvas.SetGridx();
    canvas.SetGridy();
    canvas.SetTicks(1, 1);

    axis->Draw("E1 P");
    for (size_t index = 1; index < weights.size(); ++index)
        weights[index]->Draw("E1 P SAME");

    TLine unity(0.1, 1.0, 500.0, 1.0);
    unity.SetLineStyle(2);
    unity.SetLineWidth(2);
    unity.Draw("SAME");

    TLegend legend(0.69, 0.66, 0.93, 0.87);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.039);
    legend.AddEntry(weights[0].get(), "2016 preVFP", "lep");
    legend.AddEntry(weights[1].get(), "2016 postVFP", "lep");
    legend.AddEntry(weights[2].get(), "2017", "lep");
    legend.AddEntry(weights[3].get(), "2018", "lep");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(62);
    label.SetTextSize(0.050);
    label.DrawLatex(0.13, 0.925, "CMS");
    label.SetTextFont(52);
    label.SetTextSize(0.040);
    label.DrawLatex(0.225, 0.925, "Preliminary");
    label.SetTextFont(42);
    label.SetTextAlign(31);
    label.DrawLatex(0.96, 0.925, "Run 2 (13 TeV)");
    label.SetTextAlign(11);
    label.SetTextSize(0.041);
    label.DrawLatex(
        0.17, 0.83,
        channel == "ee" ? "Z #rightarrow ee" : "Z #rightarrow #mu#mu"
    );

    const std::string outputDirectory =
        kInputBase + "/year_comparison_logx";
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);
    const std::string outputStem = outputDirectory +
        "/ZPtCorrectionFactor_year_comparison_" + channel;
    canvas.SaveAs((outputStem + ".png").c_str());
    canvas.SaveAs((outputStem + ".pdf").c_str());
    std::cout << "Created " << outputStem << ".{png,pdf}" << std::endl;
}

}

void PlotZPtWeightsFormatted(std::string period = "all") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(5);
    if (period == "all") {
        for (const std::string& item : kPeriods) plotPeriod(item);
        drawYearComparison("ee");
        drawYearComparison("mm");
        return;
    }
    if (period == "comparison") {
        drawYearComparison("ee");
        drawYearComparison("mm");
        return;
    }
    if (std::find(kPeriods.begin(), kPeriods.end(), period) == kPeriods.end()) {
        std::cerr << "Unsupported period: " << period << std::endl;
        return;
    }
    plotPeriod(period);
}
