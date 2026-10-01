#pragma once

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <string>

inline void applyCMSLabel(double x, double textSize = 0.048) {
    TLatex cms;
    cms.SetNDC();
    cms.SetTextFont(61);
    cms.SetTextSize(textSize);
    cms.DrawLatex(x, 0.93, "CMS");
    TLatex prelim;
    prelim.SetNDC();
    prelim.SetTextFont(52);
    prelim.SetTextSize(textSize * 0.9);
    prelim.DrawLatex(x + 0.09, 0.93, "Preliminary");
}

inline void applyYearLabel(const std::string& year, double textSize = 0.042) {
    TLatex yr;
    yr.SetNDC();
    yr.SetTextFont(42);
    yr.SetTextAlign(31);
    yr.SetTextSize(textSize);
    yr.DrawLatex(0.90, 0.93, year.c_str());
}

inline void drawSFAndUncertainty2D(TH2D* hSF, TH2D* hUnc, const std::string& title,
                                    const std::string& xTitle, const std::string& yTitle,
                                    const std::string& year, const std::string& outPath,
                                    double sfMin = 0.85, double sfMax = 1.10,
                                    double uncMax = 10.0, bool logX = false) {
    TCanvas c(("c_" + title).c_str(), "", 1700, 750);
    c.Divide(2, 1);

    c.cd(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.18);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.12);
    gPad->SetTicks(1, 1);
    if (logX) gPad->SetLogx();
    hSF->SetTitle("");
    hSF->GetXaxis()->SetTitle(xTitle.c_str());
    hSF->GetYaxis()->SetTitle(yTitle.c_str());
    hSF->GetZaxis()->SetTitle("Scale factor");
    hSF->GetXaxis()->SetTitleSize(0.045);
    hSF->GetYaxis()->SetTitleSize(0.045);
    hSF->GetZaxis()->SetTitleSize(0.040);
    hSF->SetMinimum(sfMin);
    hSF->SetMaximum(sfMax);
    hSF->SetMarkerSize(1.1);
    gStyle->SetPaintTextFormat(".3f");
    hSF->Draw("COLZ TEXT");
    applyCMSLabel(0.14);
    applyYearLabel(year);

    c.cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.18);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.12);
    gPad->SetTicks(1, 1);
    if (logX) gPad->SetLogx();
    hUnc->SetTitle("");
    hUnc->GetXaxis()->SetTitle(xTitle.c_str());
    hUnc->GetYaxis()->SetTitle(yTitle.c_str());
    hUnc->GetZaxis()->SetTitle("Uncertainty [%]");
    hUnc->GetXaxis()->SetTitleSize(0.045);
    hUnc->GetYaxis()->SetTitleSize(0.045);
    hUnc->GetZaxis()->SetTitleSize(0.040);
    hUnc->SetMinimum(0.0);
    hUnc->SetMaximum(uncMax);
    hUnc->SetMarkerSize(1.1);
    gStyle->SetPaintTextFormat(".3f");
    hUnc->Draw("COLZ TEXT");
    applyCMSLabel(0.14);
    applyYearLabel(year);

    gSystem->mkdir(gSystem->DirName(outPath.c_str()), kTRUE);
    c.SaveAs(outPath.c_str());
}

inline void zeroBinErrors(TH1D* h) {
    if (!h) return;
    for (int b = 1; b <= h->GetNbinsX(); ++b) h->SetBinError(b, 0.0);
}

inline void drawSF1D(TH1D* hNom, TH1D* hUp, TH1D* hDown, const std::string& title,
                      const std::string& xTitle, const std::string& yTitle,
                      const std::string& year, const std::string& outPath,
                      double yMin, double yMax) {
    TCanvas c(("c_" + title).c_str(), "", 1000, 750);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.05);
    c.SetBottomMargin(0.14);
    c.SetTopMargin(0.12);
    c.SetTicks(1, 1);

    hNom->SetTitle("");
    hNom->GetXaxis()->SetTitle(xTitle.c_str());
    hNom->GetYaxis()->SetTitle(yTitle.c_str());
    hNom->GetXaxis()->SetTitleSize(0.045);
    hNom->GetYaxis()->SetTitleSize(0.045);
    hNom->SetMinimum(yMin);
    hNom->SetMaximum(yMax);
    hNom->SetLineColor(kBlack);
    hNom->SetLineWidth(2);
    hNom->SetMarkerStyle(20);
    hNom->SetMarkerSize(0.9);
    hNom->Draw("E1");

    if (hUp) { hUp->SetLineColor(kRed + 1); hUp->SetLineStyle(2); hUp->SetLineWidth(2); hUp->Draw("HIST SAME"); }
    if (hDown) { hDown->SetLineColor(kAzure + 2); hDown->SetLineStyle(2); hDown->SetLineWidth(2); hDown->Draw("HIST SAME"); }

    TLegend* leg = new TLegend(0.62, 0.68, 0.93, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);

    leg->AddEntry(hNom, "Nominal", "lp");
    if (hUp) leg->AddEntry(hUp, "Up", "l");
    if (hDown) leg->AddEntry(hDown, "Down", "l");
    leg->Draw();

    applyCMSLabel(0.14);
    applyYearLabel(year);

    gSystem->mkdir(gSystem->DirName(outPath.c_str()), kTRUE);
    c.SaveAs(outPath.c_str());
    delete leg;
}
