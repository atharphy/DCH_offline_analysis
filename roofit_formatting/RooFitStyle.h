#pragma once
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TObject.h"
#include "TPad.h"
#include "TStyle.h"
#include "RooPlot.h"

struct RooFitLegendEntry {
    TObject* obj;
    std::string label;
    std::string drawOpt; // "f", "lep", "l", ...
};

inline std::string roofitLumiLabel(const std::string& year) {
    if (year == "2016preVFP")  return "2016 preVFP (19.5 fb^{-1})";
    if (year == "2016postVFP") return "2016 postVFP (16.8 fb^{-1})";
    if (year == "2016")        return "2016 (36.3 fb^{-1})";
    if (year == "2017")        return "2017 (41.5 fb^{-1})";
    if (year == "2018")        return "2018 (59.7 fb^{-1})";
    if (year == "Run2")        return "Run 2 (137.6 fb^{-1})";
    return year;
}

inline void applyRooFitGlobalStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetGridColor(kGray + 2);
    gStyle->SetGridStyle(2);
    gStyle->SetGridWidth(1);
}

inline void drawRooFitCMSAndLumi(const std::string& year) {
    TLatex cms;
    cms.SetNDC();
    cms.SetTextFont(61);
    cms.SetTextSize(0.046);
    cms.DrawLatex(0.17, 0.905, "CMS");

    TLatex prelim;
    prelim.SetNDC();
    prelim.SetTextFont(52);
    prelim.SetTextSize(0.044);
    prelim.DrawLatex(0.285, 0.905, "Preliminary");

    TLatex lumi;
    lumi.SetNDC();
    lumi.SetTextFont(42);
    lumi.SetTextSize(0.041);
    lumi.SetTextAlign(31);
    lumi.DrawLatex(0.96, 0.905, roofitLumiLabel(year).c_str());
}

inline TPad* formatRooFitCanvas(TCanvas& c, RooPlot* frame, const std::string& xTitle, const std::string& year, const std::vector<RooFitLegendEntry>& legendEntries, const std::string& fitResultText = "") {
    applyRooFitGlobalStyle();
    c.cd();

    TPad* p1 = new TPad((std::string(c.GetName()) + "_p1").c_str(), "", 0.00, 0.00, 0.72, 1.00);
    TPad* p3 = new TPad((std::string(c.GetName()) + "_p3").c_str(), "", 0.72, 0.00, 1.00, 1.00);

    p1->SetTopMargin(0.12); p1->SetBottomMargin(0.14); p1->SetLeftMargin(0.16); p1->SetRightMargin(0.04);
    p1->SetTicks(1, 1); p1->SetGridx(1); p1->SetGridy(1);
    p3->SetLeftMargin(0.0); p3->SetRightMargin(0.05); p3->SetTopMargin(0.05); p3->SetBottomMargin(0.05);

    p1->Draw();
    p3->Draw();

    p1->cd();
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Events");
    frame->GetYaxis()->SetTitleSize(0.052);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitle(xTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->Draw();
    p1->RedrawAxis();
    drawRooFitCMSAndLumi(year);

    p3->cd();
    TLegend* leg = new TLegend(0.0, 0.0, 1.0, 1.0);
    leg->SetTextSize(0.08);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (const auto& e : legendEntries) leg->AddEntry(e.obj, e.label.c_str(), e.drawOpt.c_str());
    leg->Draw();

    p1->cd();
    if (!fitResultText.empty()) {
        TLatex fitText;
        fitText.SetNDC();
        fitText.SetTextFont(42);
        fitText.SetTextSize(0.04);
        fitText.DrawLatex(0.40, 0.72, fitResultText.c_str());
    }

    c.cd();
    return p1;
}
