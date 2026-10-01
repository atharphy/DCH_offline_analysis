#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../../filemap/FileMap.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ZPtRW {

const std::vector<double> edges = {
    0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,
    10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.,
    30., 40., 50., 70., 100., 150., 200., 300., 500.
};

const std::vector<std::string> backgrounds = {
    "DY10_50", "VV", "ZZ", "TTbar", "ttV", "WJ", "ST",
    "VVV", "QCD", "other", "ttH", "signal"
};

double groupScale(const std::string &group) {

  if (group == "VV") return 1.12631;
  if (group == "ZZ") return 1.31366;
  return 1.0;
}

bool isDataFileForChannel(const std::string &mappedName,
                          const std::string &channel) {
  const std::string base = gSystem->BaseName(mappedName.c_str());
  if (channel == "ee")
    return base.find("SingleElectron") != std::string::npos ||
           base.find("EGamma") != std::string::npos;
  if (channel == "mm")
    return base.find("SingleMuon") != std::string::npos;
  return false;
}

std::vector<std::string> yearsToRead(const std::string &year) {
  if (year == "Run2")
    return {"2016preVFP", "2016postVFP", "2017", "2018"};
  if (year == "2016")
    return {"2016preVFP", "2016postVFP"};
  return {year};
}

TH1D *makeEmpty(const std::string &name, const std::string &title = "") {
  TH1D *h = new TH1D(name.c_str(), title.c_str(),
                     static_cast<int>(edges.size()) - 1, edges.data());
  h->Sumw2();
  h->SetDirectory(nullptr);
  return h;
}

TH1D *readRebinned(TFile *file, const std::string &histName,
                  const std::string &newName) {
  if (!file) return nullptr;
  TH1 *source = dynamic_cast<TH1 *>(file->Get(histName.c_str()));
  if (!source) return nullptr;

  TH1D *result = makeEmpty(newName, source->GetTitle());
  const int last = result->GetNbinsX();

  for (int ib = 1; ib <= source->GetNbinsX() + 1; ++ib) {
    const double content = source->GetBinContent(ib);
    const double error = source->GetBinError(ib);
    if (content == 0.0 && error == 0.0) continue;

    int target = result->FindBin(source->GetXaxis()->GetBinCenter(ib));
    if (target < 1) continue;
    if (target > last) target = last;

    result->SetBinContent(target, result->GetBinContent(target) + content);
    result->SetBinError(target,
                        std::hypot(result->GetBinError(target), error));
  }
  return result;
}

void addHistogram(TH1D *&sum, TH1D *input, const std::string &name) {
  if (!input) return;
  if (!sum) {
    sum = dynamic_cast<TH1D *>(input->Clone(name.c_str()));
    sum->SetDirectory(nullptr);
  } else {
    sum->Add(input);
  }
}

struct Inputs {
  std::unique_ptr<TH1D> data;
  std::unique_ptr<TH1D> dy;
  std::unique_ptr<TH1D> background;
  int filesOpened = 0;
  int histogramsRead = 0;
};

Inputs loadInputs(const std::string &year, const std::string &channel,
                  const std::string &histBase) {
  Inputs out;
  TH1D *data = nullptr;
  TH1D *dy = nullptr;
  TH1D *background = nullptr;
  const std::string histName = "h_zPt_" + channel + "_DYCR_0tau";

  for (const std::string &inputYear : yearsToRead(year)) {
    const auto fileMap = getFileMap(inputYear);
    const std::string directory = histBase + "/" + inputYear + "/";

    for (const auto &entry : fileMap) {
      const std::string &group = entry.first;
      const bool isData = group == "data";
      const bool isDY = group == "DY";
      const bool isBackground =
          std::find(backgrounds.begin(), backgrounds.end(), group) !=
          backgrounds.end();
      if (!isData && !isDY && !isBackground) continue;

      for (const std::string &mappedName : entry.second) {
        if (mappedName.empty()) continue;

        if (isData && !isDataFileForChannel(mappedName, channel)) continue;
        const std::string inputName =
            directory + "hist_" + gSystem->BaseName(mappedName.c_str());
        std::unique_ptr<TFile> file(TFile::Open(inputName.c_str(), "READ"));
        if (!file || file->IsZombie()) {
          std::cerr << "[missing/bad] " << inputName << std::endl;
          continue;
        }
        ++out.filesOpened;

        std::unique_ptr<TH1D> h(readRebinned(
            file.get(), histName,
            Form("tmp_%s_%s_%d", channel.c_str(), group.c_str(),
                 out.histogramsRead)));
        if (!h) {
          std::cerr << "[missing histogram] " << histName << " in "
                    << inputName << std::endl;
          continue;
        }
        ++out.histogramsRead;
        h->Scale(groupScale(group));

        if (isData)
          addHistogram(data, h.get(), "h_data_" + channel);
        else if (isDY)
          addHistogram(dy, h.get(), "h_dy_" + channel);
        else
          addHistogram(background, h.get(), "h_nonDY_" + channel);
      }
    }
  }

  out.data.reset(data);
  out.dy.reset(dy);
  out.background.reset(background ? background
                                  : makeEmpty("h_nonDY_" + channel));
  return out;
}

struct Products {
  std::unique_ptr<TH1D> dataSubtracted;
  std::unique_ptr<TH1D> dataShape;
  std::unique_ptr<TH1D> dyShape;
  std::unique_ptr<TH1D> weight;
  std::unique_ptr<TH1D> weightUp;
  std::unique_ptr<TH1D> weightDown;
};

Products derive(const Inputs &in, const std::string &channel) {
  Products p;
  p.dataSubtracted.reset(dynamic_cast<TH1D *>(
      in.data->Clone(("h_data_minus_nonDY_" + channel).c_str())));
  p.dataSubtracted->SetDirectory(nullptr);
  p.dataSubtracted->Add(in.background.get(), -1.0);

  for (int ib = 1; ib <= p.dataSubtracted->GetNbinsX(); ++ib) {
    if (p.dataSubtracted->GetBinContent(ib) < 0.0) {
      std::cerr << "[warning] negative background-subtracted yield in "
                << channel << " bin " << ib << "; setting central value to 0"
                << std::endl;
      p.dataSubtracted->SetBinContent(ib, 0.0);

    }
  }

  p.dataShape.reset(dynamic_cast<TH1D *>(
      p.dataSubtracted->Clone(("h_data_shape_" + channel).c_str())));
  p.dyShape.reset(dynamic_cast<TH1D *>(
      in.dy->Clone(("h_dy_shape_" + channel).c_str())));
  p.dataShape->SetDirectory(nullptr);
  p.dyShape->SetDirectory(nullptr);

  const double dataIntegral = p.dataShape->Integral();
  const double dyIntegral = p.dyShape->Integral();
  if (!(dataIntegral > 0.0) || !(dyIntegral > 0.0)) return p;
  p.dataShape->Scale(1.0 / dataIntegral);
  p.dyShape->Scale(1.0 / dyIntegral);

  p.weight.reset(makeEmpty("h_zpt_weight_" + channel,
                           ";reconstructed p_{T}^{Z} [GeV];Z p_{T} weight"));
  p.weightUp.reset(makeEmpty("h_zpt_weight_" + channel + "_up"));
  p.weightDown.reset(makeEmpty("h_zpt_weight_" + channel + "_down"));

  for (int ib = 1; ib <= p.weight->GetNbinsX(); ++ib) {
    const double d = p.dataShape->GetBinContent(ib);
    const double de = p.dataShape->GetBinError(ib);
    const double m = p.dyShape->GetBinContent(ib);
    const double me = p.dyShape->GetBinError(ib);

    double value = 1.0;
    double error = 0.0;
    if (m > 0.0) {
      value = d / m;
      error = std::sqrt((de / m) * (de / m) +
                        (d * me / (m * m)) * (d * me / (m * m)));
    } else {
      std::cerr << "[warning] zero DY denominator in " << channel << " bin "
                << ib << "; using weight 1" << std::endl;
    }

    p.weight->SetBinContent(ib, value);
    p.weight->SetBinError(ib, error);
    p.weightUp->SetBinContent(ib, value + error);
    p.weightDown->SetBinContent(ib, std::max(0.0, value - error));
  }
  return p;
}

void drawDiagnostic(const std::string &channel, const Products &p,
                    const std::string &outputName) {
  if (!p.dataShape || !p.dyShape || !p.weight) return;
  TCanvas canvas(("c_zpt_" + channel).c_str(), "", 850, 850);
  TPad top(("top_" + channel).c_str(), "", 0., 0.31, 1., 1.);
  TPad bottom(("bottom_" + channel).c_str(), "", 0., 0., 1., 0.31);
  top.SetBottomMargin(0.03);
  bottom.SetTopMargin(0.04);
  bottom.SetBottomMargin(0.32);
  top.Draw();
  bottom.Draw();

  top.cd();
  p.dataShape->SetTitle("");
  p.dataShape->GetYaxis()->SetTitle("Normalized events");
  p.dataShape->GetXaxis()->SetLabelSize(0.0);
  p.dataShape->SetMarkerStyle(20);
  p.dataShape->SetLineColor(kBlack);
  p.dyShape->SetLineColor(kRed + 1);
  p.dyShape->SetLineWidth(2);
  p.dataShape->SetMaximum(1.35 * std::max(p.dataShape->GetMaximum(),
                                          p.dyShape->GetMaximum()));
  p.dataShape->Draw("E1");
  p.dyShape->Draw("HIST SAME");
  p.dataShape->Draw("E1 SAME");
  TLegend legend(0.55, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(p.dataShape.get(), "Data - non-DY", "lep");
  legend.AddEntry(p.dyShape.get(), "DY simulation", "l");
  legend.Draw();

  bottom.cd();
  p.weight->SetMarkerStyle(20);
  p.weight->SetLineColor(kBlack);
  p.weight->GetXaxis()->SetTitle("reconstructed p_{T}^{Z} [GeV]");
  p.weight->GetYaxis()->SetTitle("Shape weight");
  p.weight->GetYaxis()->SetNdivisions(505);
  p.weight->GetXaxis()->SetTitleSize(0.12);
  p.weight->GetXaxis()->SetLabelSize(0.10);
  p.weight->GetYaxis()->SetTitleSize(0.10);
  p.weight->GetYaxis()->SetLabelSize(0.09);
  p.weight->GetYaxis()->SetTitleOffset(0.55);
  p.weight->SetMinimum(0.0);
  p.weight->SetMaximum(std::max(2.0, 1.25 * p.weight->GetMaximum()));
  p.weight->Draw("E1");
  TLine unity;
  unity.SetLineStyle(2);
  unity.DrawLine(edges.front(), 1.0, edges.back(), 1.0);
  canvas.SaveAs(outputName.c_str());
}

void writeTable(const TH1D *weight, const std::string &fileName) {
  std::ofstream out(fileName);
  out << "# reco_zpt_low reco_zpt_high nominal up down\n";
  out << std::fixed << std::setprecision(8);
  for (int ib = 1; ib <= weight->GetNbinsX(); ++ib) {
    const double value = weight->GetBinContent(ib);
    const double error = weight->GetBinError(ib);
    out << weight->GetXaxis()->GetBinLowEdge(ib) << " "
        << weight->GetXaxis()->GetBinUpEdge(ib) << " " << value << " "
        << value + error << " " << std::max(0.0, value - error) << "\n";
  }
}

}

void DeriveZPtWeights(
    std::string year = "2018",
    std::string histBase = "hists/run2_hists_noFR_all_channels_zpt",
    std::string outputBase = "../../Dependencies/zpt_weights") {
  gROOT->SetBatch(kTRUE);
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  const std::vector<std::string> allowed = {
      "2016preVFP", "2016postVFP", "2016", "2017", "2018", "Run2"};
  if (std::find(allowed.begin(), allowed.end(), year) == allowed.end()) {
    std::cerr << "ERROR: unsupported year " << year << std::endl;
    return;
  }

  const std::string outputDir = outputBase + "/" + year;
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  const std::string rootName = outputDir + "/ZPtWeights_" + year + ".root";
  TFile output(rootName.c_str(), "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "ERROR: cannot create " << rootName << std::endl;
    return;
  }

  int channelsWritten = 0;

  for (const std::string &channel : {std::string("ee"), std::string("mm")}) {
    ZPtRW::Inputs inputs = ZPtRW::loadInputs(year, channel, histBase);
    if (!inputs.data || !inputs.dy) {
      std::cerr << "ERROR: missing summed Data or DY histogram for " << channel
                << std::endl;
      continue;
    }
    std::cout << "[" << channel << "] opened files=" << inputs.filesOpened
              << ", histograms=" << inputs.histogramsRead
              << ", Data=" << inputs.data->Integral()
              << ", non-DY=" << inputs.background->Integral()
              << ", DY=" << inputs.dy->Integral() << std::endl;

    ZPtRW::Products products = ZPtRW::derive(inputs, channel);
    if (!products.weight) {
      std::cerr << "ERROR: non-positive Data-subtracted or DY integral for "
                << channel << std::endl;
      continue;
    }

    output.cd();
    inputs.data->Write();
    inputs.background->Write();
    inputs.dy->Write();
    products.dataSubtracted->Write();
    products.dataShape->Write();
    products.dyShape->Write();
    products.weight->Write();
    products.weightUp->Write();
    products.weightDown->Write();

    ZPtRW::drawDiagnostic(channel, products,
                          outputDir + "/ZPtWeight_" + channel + ".png");
    ZPtRW::writeTable(products.weight.get(),
                      outputDir + "/ZPtWeight_" + channel + ".txt");
    ++channelsWritten;
  }

  output.Close();

  if (channelsWritten != 2) {
    std::cerr << "ERROR: produced " << channelsWritten
              << " of 2 requested Z pT correction channels for " << year
              << std::endl;
    gSystem->Exit(2);
    return;
  }

  std::cout << "Z pT correction outputs written to " << outputDir << std::endl;
}
