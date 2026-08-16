#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TSystem.h"

#include "HistCache.h"
#include "PlotUtils.h"
#include "StackConfig.h"
#include "SystematicSources.h"
#include "Labels.h"

struct GroupHists {
    TH1D* nominal = nullptr;
    std::map<std::string, TH1D*> up;
    std::map<std::string, TH1D*> down;
};

inline void drawAndSave(const std::string& hname, const std::vector<std::string>& sourceNames, std::map<std::string, std::vector<TFile*>>& handles, std::map<std::string, int>& fill_colors, const std::string& outdir, bool useLog, const std::vector<SystSource>& activeSources, bool groupByChannel, bool blind = false) {
    std::map<std::string, GroupHists> G;

    for (auto& kv : handles) {
        const std::string& grp = kv.first;
        GroupHists gh;

        for (TFile* f : kv.second) {
            if (!f || f->IsZombie()) continue;

            for (const std::string& sourceName : sourceNames) {
                TH1D* h = fetchAndClone(f, sourceName);
                if (h) {
                    if (!gh.nominal) {
                        gh.nominal = h;
                        gh.nominal->SetName((grp + "_" + hname).c_str());
                        h = nullptr;
                    } else {
                        gh.nominal->Add(h);
                    }
                    delete h;
                }

                if (!DRAW_BANDS) continue;

                for (const auto& sys : activeSources) {
                    TH1D* hup = fetchAndClone(f, sourceName + sys.upSuffix);
                    if (hup) {
                        TH1D*& acc = gh.up[sys.name];
                        if (!acc) {
                            acc = hup;
                            acc->SetName((grp + "_" + hname + "_" + sys.name + "Up").c_str());
                            hup = nullptr;
                        } else {
                            acc->Add(hup);
                        }
                        delete hup;
                    }

                    TH1D* hdn = fetchAndClone(f, sourceName + sys.downSuffix);
                    if (hdn) {
                        TH1D*& acc = gh.down[sys.name];
                        if (!acc) {
                            acc = hdn;
                            acc->SetName((grp + "_" + hname + "_" + sys.name + "Down").c_str());
                            hdn = nullptr;
                        } else {
                            acc->Add(hdn);
                        }
                        delete hdn;
                    }
                }
            }
        }

        if (gh.nominal) {
            G[grp] = gh;
        } else {
            for (auto& u : gh.up) delete u.second;
            for (auto& d : gh.down) delete d.second;
        }
    }

    if (G.empty()) {
        std::cout << "[skip] " << hname << std::endl;
        return;
    }

    auto scaleGroup = [&](const std::string& grp, double factor) {
        if (!G.count(grp)) return;
        G[grp].nominal->Scale(factor);
        for (auto& u : G[grp].up) u.second->Scale(factor);
        for (auto& d : G[grp].down) d.second->Scale(factor);
    };

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
    auto scaleGroupFromCSV = [&](const std::string& grp, const std::string& csvPath) {
        if (!G.count(grp)) return;
        const auto yearMap = readSFFromCSV(csvPath);
        auto it = yearMap.find(year);
        if (it == yearMap.end()) {
            std::cerr << "[warning] no " << grp << " scale factor for year " << year << " in " << csvPath << ", leaving unscaled" << std::endl;
            return;
        }
        scaleGroup(grp, it->second);
    };
    scaleGroupFromCSV("WZ", "wz_scale_factors.csv");
    scaleGroupFromCSV("ZZ", "zz_scale_factors.csv");

    auto clampNegativeBins = [](TH1D* h) {
        if (!h) return;
        for (int ib = 0; ib <= h->GetNbinsX() + 1; ++ib) {
            if (h->GetBinContent(ib) < 0.0) {
                h->SetBinContent(ib, 0.0);
                h->SetBinError(ib, 0.0);
            }
        }
    };
    for (auto& kv : G) {
        if (kv.first == "data") continue;
        clampNegativeBins(kv.second.nominal);
        for (auto& u : kv.second.up) clampNegativeBins(u.second);
        for (auto& d : kv.second.down) clampNegativeBins(d.second);
    }

    std::map<std::string, double> rawIntegral;
    for (auto& kv : G) rawIntegral[kv.first] = kv.second.nominal->Integral(1, kv.second.nominal->GetNbinsX());

    if (EVENTS_PER_BIN_WIDTH) {
        for (auto& kv : G) {
            divideByBinWidth(kv.second.nominal);
            for (auto& u : kv.second.up) divideByBinWidth(u.second);
            for (auto& d : kv.second.down) divideByBinWidth(d.second);
        }
    }

    const std::vector<std::string> bkgOrder = {"DY", "DY10_50", "WZ", "WW", "ZZ", "TTbar", "ttV", "WJ", "ST", "VVV", "QCD", "other", "signal"};

    THStack* st = new THStack(("st_" + hname).c_str(), "");
    TH1D* h_bkg = nullptr;

    for (const std::string& grp : bkgOrder) {
        if (!G.count(grp)) continue;

        G[grp].nominal->SetFillColor(fill_colors.count(grp) ? fill_colors[grp] : kGray);
        G[grp].nominal->SetLineColor(kBlack);
        G[grp].nominal->SetTitle("");

        st->Add(G[grp].nominal, "hist");

        if (!h_bkg) {
            h_bkg = (TH1D*)G[grp].nominal->Clone(("bkg_" + hname).c_str());
            h_bkg->SetDirectory(nullptr);
        } else {
            h_bkg->Add(G[grp].nominal);
        }
    }

    if (!h_bkg) {
        std::cout << "[skip no bkg] " << hname << std::endl;
        delete st;
        for (auto& kv : G) {
            delete kv.second.nominal;
            for (auto& u : kv.second.up) delete u.second;
            for (auto& d : kv.second.down) delete d.second;
        }
        return;
    }

    std::map<std::string, TH1D*> h_bkg_up_by_source, h_bkg_down_by_source;
    if (DRAW_BANDS) {
        for (const auto& sys : activeSources) {
            TH1D* bup = nullptr;
            TH1D* bdn = nullptr;
            for (const std::string& grp : bkgOrder) {
                if (!G.count(grp)) continue;
                TH1D* upSrc = G[grp].up.count(sys.name) ? G[grp].up[sys.name] : G[grp].nominal;
                TH1D* dnSrc = G[grp].down.count(sys.name) ? G[grp].down[sys.name] : G[grp].nominal;

                if (!bup) { bup = (TH1D*)upSrc->Clone(("bkg_up_" + sys.name + "_" + hname).c_str()); bup->SetDirectory(nullptr); }
                else bup->Add(upSrc);

                if (!bdn) { bdn = (TH1D*)dnSrc->Clone(("bkg_down_" + sys.name + "_" + hname).c_str()); bdn->SetDirectory(nullptr); }
                else bdn->Add(dnSrc);
            }
            if (bup) h_bkg_up_by_source[sys.name] = bup;
            if (bdn) h_bkg_down_by_source[sys.name] = bdn;
        }
    }

    double ymax = h_bkg->GetMaximum();
    if (G.count("data")) ymax = std::max(ymax, G["data"].nominal->GetMaximum());

    if (ymax <= 0.0) {
        std::cout << "[skip empty] " << hname << std::endl;
        delete st;
        delete h_bkg;
        for (auto& kv : h_bkg_up_by_source) delete kv.second;
        for (auto& kv : h_bkg_down_by_source) delete kv.second;
        for (auto& kv : G) {
            delete kv.second.nominal;
            for (auto& u : kv.second.up) delete u.second;
            for (auto& d : kv.second.down) delete d.second;
        }
        return;
    }

    auto bandErrorAt = [&](int ib, double nom) {
        double stat = h_bkg->GetBinError(ib);
        double sumSqHigh = stat * stat;
        double sumSqLow = stat * stat;
        for (const auto& sys : activeSources) {
            double up = h_bkg_up_by_source.count(sys.name) ? h_bkg_up_by_source[sys.name]->GetBinContent(ib) : nom;
            double down = h_bkg_down_by_source.count(sys.name) ? h_bkg_down_by_source[sys.name]->GetBinContent(ib) : nom;
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

    TGraphAsymmErrors* g_unc_band = nullptr;
    if (DRAW_BANDS) {
        g_unc_band = new TGraphAsymmErrors(h_bkg->GetNbinsX());
        for (int ib = 1; ib <= h_bkg->GetNbinsX(); ++ib) {
            const double x = h_bkg->GetXaxis()->GetBinCenter(ib);
            const double ex = 0.5 * h_bkg->GetXaxis()->GetBinWidth(ib);
            const double nom = h_bkg->GetBinContent(ib);
            const auto err = bandErrorAt(ib, nom);
            g_unc_band->SetPoint(ib - 1, x, nom);
            g_unc_band->SetPointError(ib - 1, ex, ex, err.first, err.second);
        }
        g_unc_band->SetFillColor(kGray + 2);
        g_unc_band->SetFillStyle(3004);
        g_unc_band->SetLineColor(kGray + 2);
        g_unc_band->SetMarkerSize(0);
    }

    bool haveData = !blind && G.count("data") && h_bkg && h_bkg->Integral() > 0.0;

    TCanvas* c = new TCanvas(("c_" + hname + (blind ? "_blinded" : "")).c_str(), "", 1000, 800);
    TPad* p1 = nullptr;
    TPad* p2 = nullptr;
    TPad* p3 = nullptr;

    p1 = new TPad(("p1_" + hname).c_str(), "", 0.00, 0.32, 0.72, 1.00);
    p2 = new TPad(("p2_" + hname).c_str(), "", 0.00, 0.00, 0.72, 0.32);
    p3 = new TPad(("p3_" + hname).c_str(), "", 0.72, 0.00, 1.00, 1.00);

    p1->SetTopMargin(0.12); p1->SetBottomMargin(0.03); p1->SetLeftMargin(0.16); p1->SetRightMargin(0.04);
    p1->SetTicks(1, 1); p1->SetGridx(1); p1->SetGridy(1);

    p2->SetTopMargin(0.05); p2->SetBottomMargin(0.35); p2->SetLeftMargin(0.16); p2->SetRightMargin(0.04);
    p2->SetTicks(1, 1); p2->SetGridx(1); p2->SetGridy(1);

    p3->SetLeftMargin(0.0); p3->SetRightMargin(0.05); p3->SetTopMargin(0.05); p3->SetBottomMargin(0.05);

    p1->Draw(); p2->Draw(); p3->Draw();
    p1->cd();

    if (useLog) {
        double logMin = 1.0;
        if (ymax < 10.0) logMin = 0.01;
        if (ymax < 1.0) logMin = 0.001;
        st->SetMinimum(logMin);
        st->SetMaximum(std::max(10.0 * logMin, 50.0 * ymax));
        p1->SetLogy(1);
    } else {
        st->SetMinimum(0.0);
        st->SetMaximum(1.5 * ymax);
    }

    st->Draw("hist");
    st->SetTitle("");
    st->GetYaxis()->SetTitle(EVENTS_PER_BIN_WIDTH ? "Events / bin" : "Events");
    st->GetYaxis()->SetTitleSize(0.052);
    st->GetYaxis()->SetTitleOffset(1.5);
    st->GetYaxis()->SetLabelSize(0.045);

    st->GetXaxis()->SetTitle("");
    st->GetXaxis()->SetLabelSize(0.0);
    st->GetXaxis()->SetTitleSize(0.0);
    st->GetXaxis()->SetTickLength(0.0);

    if (DRAW_BANDS && g_unc_band) g_unc_band->Draw("E2 SAME");

    if (!blind && G.count("data")) {
        G["data"].nominal->SetMarkerStyle(20);
        G["data"].nominal->SetMarkerColor(kBlack);
        G["data"].nominal->SetLineColor(kBlack);
        G["data"].nominal->SetTitle("");
        G["data"].nominal->Draw("E SAME");
    }

    TLegend* leg = new TLegend(0.05, 0.05, 0.95, 0.95);
    leg->SetTextSize(0.08);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    auto labelWithIntegral = [&](const std::string& grp) {
        std::ostringstream ss;
        ss.setf(std::ios::fixed);
        ss.precision(1);
        ss << grp << " (" << (rawIntegral.count(grp) ? rawIntegral[grp] : 0.0) << ")";
        return ss.str();
    };

    for (const std::string& grp : bkgOrder) {
        if (G.count(grp)) leg->AddEntry(G[grp].nominal, labelWithIntegral(grp).c_str(), "f");
    }

    if (!blind && G.count("data")) {
        std::ostringstream ds;
        ds.setf(std::ios::fixed);
        ds.precision(1);
        ds << "Data (" << rawIntegral["data"] << ")";
        leg->AddEntry(G["data"].nominal, ds.str().c_str(), "lep");
    }

    if (DRAW_BANDS && g_unc_band) leg->AddEntry(g_unc_band, "Stat. #oplus syst. unc.", "f");

    TLegend* legClone = nullptr;
    if (p3) {
        p3->cd();
        legClone = (TLegend*)leg->Clone();
        legClone->SetX1NDC(0.0);
        legClone->SetX2NDC(1.0);
        legClone->SetY1NDC(0.0);
        legClone->SetY2NDC(1.0);
        legClone->Draw();
    }

    p1->cd(); p1->RedrawAxis(); drawCMSAndLumi();

    TH1D* ratio = nullptr;
    TH1D* ratioFrame = nullptr;
    TGraphAsymmErrors* g_ratio_band = nullptr;

    {
        p2->cd();

        ratioFrame = (TH1D*)h_bkg->Clone(("ratioframe_" + hname).c_str());
        ratioFrame->SetDirectory(nullptr);
        ratioFrame->Reset();
        ratioFrame->SetTitle("");
        ratioFrame->GetYaxis()->SetTitle(haveData ? "Data/Simulation" : "Ratio");
        ratioFrame->GetYaxis()->SetNdivisions(505);
        ratioFrame->GetYaxis()->SetTitleSize(0.1);
        ratioFrame->GetYaxis()->SetTitleOffset(0.5);
        ratioFrame->GetYaxis()->SetLabelSize(0.08);
        ratioFrame->GetXaxis()->SetTitle(getXaxisLabel(hname).c_str());
        ratioFrame->GetXaxis()->SetTitleSize(0.12);
        ratioFrame->GetXaxis()->SetTitleOffset(1.05);
        ratioFrame->GetXaxis()->SetLabelSize(0.09);
        ratioFrame->SetMinimum(0.0);
        ratioFrame->SetMaximum(2.0);

        ratioFrame->Draw("");

        if (DRAW_BANDS) {
            g_ratio_band = new TGraphAsymmErrors(h_bkg->GetNbinsX());
            for (int ib = 1; ib <= h_bkg->GetNbinsX(); ++ib) {
                const double nom = h_bkg->GetBinContent(ib);
                const double x = h_bkg->GetXaxis()->GetBinCenter(ib);
                const double ex = 0.5 * h_bkg->GetXaxis()->GetBinWidth(ib);

                double eyHigh = 0.0, eyLow = 0.0;
                if (nom > 0.0) {
                    const auto err = bandErrorAt(ib, nom);
                    eyLow = err.first / nom;
                    eyHigh = err.second / nom;
                }

                g_ratio_band->SetPoint(ib - 1, x, 1.0);
                g_ratio_band->SetPointError(ib - 1, ex, ex, eyLow, eyHigh);
            }
            g_ratio_band->SetFillColor(kGray + 2);
            g_ratio_band->SetFillStyle(3004);
            g_ratio_band->SetLineColor(kGray + 2);
            g_ratio_band->SetMarkerSize(0);
            g_ratio_band->Draw("E2 SAME");
        }

        if (haveData) {
            ratio = (TH1D*)G["data"].nominal->Clone(("ratio_" + hname).c_str());
            ratio->SetDirectory(nullptr);
            ratio->Divide(h_bkg);
            ratio->SetTitle("");
            ratio->SetMarkerStyle(20);
            ratio->SetMarkerColor(kBlack);
            ratio->SetLineColor(kBlack);
            ratio->Draw("E SAME");
        }

        TLine L;
        L.SetLineStyle(2);
        L.DrawLine(ratioFrame->GetXaxis()->GetXmin(), 1.0, ratioFrame->GetXaxis()->GetXmax(), 1.0);
        p2->RedrawAxis();
    }

    std::string outsub = outdir;
    if (groupByChannel) {
        outsub = outdir + "/" + channelFromHname(hname);
        gSystem->mkdir(outsub.c_str(), kTRUE);
    }

    std::string png = outsub + "/" + hname + (blind ? "_blinded" : "") + ".png";
    c->SaveAs(png.c_str());

    if (!gSystem->AccessPathName(png.c_str())) std::cout << "[saved] " << png << std::endl;
    else std::cerr << "[save failed] " << png << std::endl;

    delete ratio;
    delete ratioFrame;
    delete g_ratio_band;
    delete g_unc_band;
    delete h_bkg;
    for (auto& kv : h_bkg_up_by_source) delete kv.second;
    for (auto& kv : h_bkg_down_by_source) delete kv.second;

    delete legClone;
    delete leg;
    delete st;

    for (auto& kv : G) {
        delete kv.second.nominal;
        for (auto& u : kv.second.up) delete u.second;
        for (auto& d : kv.second.down) delete d.second;
    }

    delete p1;
    delete p2;
    delete p3;

    c->Close();
    delete c;

    gROOT->cd();
}
