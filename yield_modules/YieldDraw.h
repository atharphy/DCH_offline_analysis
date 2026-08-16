#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TSystem.h"

#include "YieldConfig.h"
#include "YieldHistUtils.h"
#include "YieldSystematics.h"

inline void drawRegion(const std::string& region, const std::string& year, bool useLog, bool drawBands, const std::map<std::string, ProcessYields>& yields, const std::string& outputDirectory, bool blind = false) {
    TCanvas canvas(("c_" + region).c_str(), "", 1000, 800);

    TPad mainPad("mainPad", "", 0.00, 0.32, 0.72, 1.00);
    TPad ratioPad("ratioPad", "", 0.00, 0.00, 0.72, 0.32);
    TPad legendPad("legendPad", "", 0.72, 0.00, 1.00, 1.00);

    mainPad.SetLeftMargin(0.16);
    mainPad.SetRightMargin(0.04);
    mainPad.SetTopMargin(0.12);
    mainPad.SetBottomMargin(0.03);
    mainPad.SetLogy(useLog);
    mainPad.SetTicks(1, 1);
    mainPad.SetGridx(1);
    mainPad.SetGridy(1);

    ratioPad.SetLeftMargin(0.16);
    ratioPad.SetRightMargin(0.04);
    ratioPad.SetTopMargin(0.05);
    ratioPad.SetBottomMargin(0.35);
    ratioPad.SetTicks(1, 1);
    ratioPad.SetGrid(1, 1);

    legendPad.SetLeftMargin(0.0);
    legendPad.SetRightMargin(0.05);
    legendPad.SetTopMargin(0.05);
    legendPad.SetBottomMargin(0.05);

    mainPad.Draw();
    ratioPad.Draw();
    legendPad.Draw();

    mainPad.cd();
    THStack stack(("stack_" + region).c_str(), "");
    std::unique_ptr<TH1D> totalMC = makeYieldHistogram("total_mc_" + region);
    std::vector<std::string> activeProcesses;

    for (const std::string& process : kProcesses) {
        auto found = yields.find(process);
        if (found == yields.end() || found->second.nominal->Integral() <= 0.0) continue;
        TH1D* histogram = found->second.nominal.get();
        histogram->SetFillColor(kColors.at(process));
        histogram->SetLineColor(kBlack);
        histogram->SetLineWidth(1);
        stack.Add(histogram);
        totalMC->Add(histogram);
        activeProcesses.push_back(process);
    }

    const TH1D* data = nullptr;
    auto dataFound = yields.find("data");
    if (dataFound != yields.end() && !blind) data = dataFound->second.nominal.get();

    double maximum = totalMC->GetMaximum();
    if (data) maximum = std::max(maximum, data->GetMaximum());
    if (!(maximum > 0.0)) maximum = 1.0;

    stack.Draw("HIST");
    stack.GetYaxis()->SetTitle("Events");
    stack.GetYaxis()->SetTitleSize(0.052);
    stack.GetYaxis()->SetTitleOffset(1.5);
    stack.GetYaxis()->SetLabelSize(0.045);
    stack.GetXaxis()->SetLabelSize(0.0);
    stack.GetXaxis()->SetTickLength(0.0);
    double logMinimum = 1.0;
    if (maximum < 10.0) logMinimum = 0.01;
    if (maximum < 1.0) logMinimum = 0.001;
    stack.SetMinimum(useLog ? logMinimum : 0.0);
    stack.SetMaximum(useLog ? std::max(10.0 * logMinimum, 50.0 * maximum) : maximum * 1.5);

    std::map<std::string, std::unique_ptr<TH1D>> totalUp, totalDown;
    if (drawBands) {
        for (const auto& src : kAllSystSources) {
            auto up = makeYieldHistogram("total_up_" + src.name + "_" + region);
            auto down = makeYieldHistogram("total_down_" + src.name + "_" + region);
            for (const std::string& process : activeProcesses) {
                const ProcessYields& py = yields.at(process);
                TH1D* upSrc = py.up.count(src.name) ? py.up.at(src.name).get() : py.nominal.get();
                TH1D* downSrc = py.down.count(src.name) ? py.down.at(src.name).get() : py.nominal.get();
                up->Add(upSrc);
                down->Add(downSrc);
            }
            totalUp[src.name] = std::move(up);
            totalDown[src.name] = std::move(down);
        }
    }

    auto bandErrorAt = [&](int bin, double nom) {
        double stat = totalMC->GetBinError(bin);
        double sumSqHigh = stat * stat;
        double sumSqLow = stat * stat;
        for (const auto& src : kAllSystSources) {
            if (!totalUp.count(src.name)) continue;
            double up = totalUp.at(src.name)->GetBinContent(bin);
            double down = totalDown.at(src.name)->GetBinContent(bin);
            if (up < nom) up = nom;
            if (down > nom) down = nom;
            if (down < 0.0) down = 0.0;
            const double devHigh = up - nom;
            const double devLow = nom - down;
            sumSqHigh += devHigh * devHigh;
            sumSqLow += devLow * devLow;
        }
        return std::make_pair(std::sqrt(sumSqLow), std::sqrt(sumSqHigh));
    };

    std::unique_ptr<TGraphAsymmErrors> bandMain;
    if (drawBands) {
        bandMain.reset(new TGraphAsymmErrors(7));
        for (int bin = 1; bin <= 7; ++bin) {
            const double nom = totalMC->GetBinContent(bin);
            const auto err = bandErrorAt(bin, nom);
            bandMain->SetPoint(bin - 1, bin, nom);
            bandMain->SetPointError(bin - 1, 0.5, 0.5, err.first, err.second);
        }
        bandMain->SetFillColor(kGray + 2);
        bandMain->SetFillStyle(3004);
        bandMain->SetLineColor(kGray + 2);
        bandMain->SetMarkerSize(0);
        bandMain->Draw("E2 SAME");
    }

    if (data) {
        const_cast<TH1D*>(data)->SetMarkerStyle(20);
        const_cast<TH1D*>(data)->SetMarkerSize(1.1);
        const_cast<TH1D*>(data)->SetLineColor(kBlack);
        const_cast<TH1D*>(data)->Draw("E1 SAME");
    }

    TLatex text;
    text.SetNDC();
    text.SetTextFont(61);
    text.SetTextSize(0.068);
    text.DrawLatex(0.17, 0.905, "CMS");
    text.SetTextFont(52);
    text.SetTextSize(0.064);
    text.DrawLatex(0.285, 0.905, "Preliminary");
    text.SetTextFont(42);
    text.SetTextAlign(31);
    text.SetTextSize(0.060);
    text.DrawLatex(0.96, 0.905, periodLabel(year).c_str());

    ratioPad.cd();
    std::unique_ptr<TH1D> ratio = makeYieldHistogram("ratio_" + region);
    std::unique_ptr<TH1D> ratioBand = makeYieldHistogram("ratio_band_" + region);

    for (int bin = 1; bin <= 7; ++bin) {
        const double prediction = totalMC->GetBinContent(bin);
        const double predictionError = totalMC->GetBinError(bin);
        const double observation = data ? data->GetBinContent(bin) : 0.0;
        const double observationError = data ? data->GetBinError(bin) : 0.0;

        if (prediction > 0.0) {
            ratio->SetBinContent(bin, observation / prediction);
            ratio->SetBinError(bin, observationError / prediction);
            ratioBand->SetBinContent(bin, 1.0);
            ratioBand->SetBinError(bin, predictionError / prediction);
        }
    }

    ratioBand->SetTitle("");
    ratioBand->SetMinimum(0.0);
    ratioBand->SetMaximum(2.0);
    ratioBand->GetYaxis()->SetTitle(blind ? "Ratio" : "Data/Simulation");
    ratioBand->GetYaxis()->SetNdivisions(505);
    ratioBand->GetYaxis()->SetTitleSize(0.1);
    ratioBand->GetYaxis()->SetTitleOffset(0.5);
    ratioBand->GetYaxis()->SetLabelSize(0.08);
    ratioBand->GetXaxis()->SetLabelSize(0.12);
    ratioBand->GetXaxis()->SetLabelOffset(0.02);

    ratioBand->Draw("");

    std::unique_ptr<TGraphAsymmErrors> bandRatio;
    if (drawBands) {
        bandRatio.reset(new TGraphAsymmErrors(7));
        for (int bin = 1; bin <= 7; ++bin) {
            const double nom = totalMC->GetBinContent(bin);
            double eyHigh = 0.0, eyLow = 0.0;
            if (nom > 0.0) {
                const auto err = bandErrorAt(bin, nom);
                eyLow = err.first / nom;
                eyHigh = err.second / nom;
            }
            bandRatio->SetPoint(bin - 1, bin, 1.0);
            bandRatio->SetPointError(bin - 1, 0.5, 0.5, eyLow, eyHigh);
        }
        bandRatio->SetFillColor(kGray + 2);
        bandRatio->SetFillStyle(3004);
        bandRatio->SetLineColor(kGray + 2);
        bandRatio->SetMarkerSize(0);
        bandRatio->Draw("E2 SAME");
    }

    if (!blind) {
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(1.0);
        ratio->SetLineColor(kBlack);
        ratio->Draw("E1 SAME");
    }

    TLine unity(0.5, 1.0, 7.5, 1.0);
    unity.SetLineStyle(2);
    unity.Draw();
    ratioPad.RedrawAxis();
    ratioPad.Modified();
    ratioPad.Update();

    legendPad.cd();
    TLegend legend(0.05, 0.05, 0.95, 0.95);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.08);

    for (const std::string& process : activeProcesses) {
        legend.AddEntry(yields.at(process).nominal.get(), kLabels.at(process).c_str(), "f");
    }
    if (data) legend.AddEntry(data, "Data", "lep");
    if (drawBands && bandMain) legend.AddEntry(bandMain.get(), "Stat. #oplus syst. unc.", "f");
    legend.Draw();

    canvas.cd();
    const std::string suffix = useLog ? "_log" : "_linear";
    const std::string baseName = outputDirectory + "/yield_" + region + "_" + year + suffix;
    canvas.SaveAs((baseName + ".png").c_str());
    canvas.SaveAs((baseName + ".pdf").c_str());
}
