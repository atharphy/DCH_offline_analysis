// One-off companion to BuildQCDFakeRates.C: measures the QCD-enriched-region
// electron and muon fake rates in 10 SIGNED |eta|<2.5 bins (-2.5 to 2.5),
// pT-marginalized, combined across all 4 Run2 eras for statistics. The main
// (adopted) BuildQCDFakeRates.C stays untouched -- this is a diagnostic
// eta-shape check, reusing that script's exact region cuts and electron
// conversion-correction logic verbatim (see its own header comment for the
// physics rationale), just with different binning and a different input
// path (reads directly from the per-sample-group merged files under
// QCD_FakeRate_merged/<year>/{Data,MC/<category>}/*.root via RDataFrame's
// native multi-file support, instead of requiring merge_qcd_fakerate_skims.py's
// single-file-per-category output, which doesn't exist at this path yet).
//
// Run: root -l -b -q BuildQCDFakeRatesEtaScan.C+
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

const std::string SKIM_BASE = "/eos/user/a/atahmad/DCH_offline_analysis/QCD_FakeRate_merged_crossflavor";
const std::vector<std::string> YEARS = {"2016preVFP", "2016postVFP", "2017", "2018"};
const std::string ELECTRON_CORRECTION_DONOR = "QCD_EMEnriched";

// Same AN-19-111/X53-AN2016 baseline region definition as BuildQCDFakeRates.C.
const double MET_MAX = 25.0;
const double MT_MAX = 25.0;
const double ZVETO_LOW = 81.1;
const double ZVETO_HIGH = 101.1;
const std::map<std::string, double> DEEPJET_MEDIUM_WP = {
    {"2016preVFP", 0.2598}, {"2016postVFP", 0.2489}, {"2017", 0.3040}, {"2018", 0.2783},
};

// pT edges kept the same as the main script (not the focus of this scan);
// eta is now 10 SIGNED bins spanning the full -2.5..2.5 range for both
// flavors, per request -- note muon's real acceptance is |eta|<2.4, so its
// outermost bins will be sparse/empty by physics, not a bug.
const std::vector<double> PT_EDGES = {20, 30, 45, 70, 100, 200};
const std::vector<double> ETA_EDGES = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

struct FlavorDef { std::string label; int catCode; };
const std::vector<FlavorDef> FLAVORS = {{"electron", 1}, {"muon", 2}};

std::vector<std::string> globFiles(const std::string& pattern) {
    std::vector<std::string> out;
    TString cmd = Form("ls %s 2>/dev/null", pattern.c_str());
    FILE* p = gSystem->OpenPipe(cmd.Data(), "r");
    if (!p) return out;
    char buf[4096];
    while (fgets(buf, sizeof(buf), p)) {
        std::string s(buf);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        if (!s.empty()) out.push_back(s);
    }
    gSystem->ClosePipe(p);
    return out;
}

std::unique_ptr<TH2D> makeHist2D(const std::string& name) {
    auto h = std::make_unique<TH2D>(name.c_str(), "", (int)PT_EDGES.size() - 1, PT_EDGES.data(),
                                     (int)ETA_EDGES.size() - 1, ETA_EDGES.data());
    h->SetDirectory(nullptr);
    h->Sumw2();
    return h;
}

void addInPlace(TH2D& target, const TH2D& add) { target.Add(&add); }

void subtractInPlace(TH2D& target, const TH2D& toSubtract) {
    for (int bx = 1; bx <= target.GetNbinsX(); ++bx)
        for (int by = 1; by <= target.GetNbinsY(); ++by) {
            const double v = std::max(0.0, target.GetBinContent(bx, by) - toSubtract.GetBinContent(bx, by));
            target.SetBinContent(bx, by, v);
        }
}

struct FlavorHistos {
    std::unique_ptr<TH2D> loose, tight, looseReal, tightReal;
};

// Accumulates one (sample, flavor) across all 4 years into a single combined
// FlavorHistos, by summing each year's freshly-filled histos in place.
FlavorHistos fillFlavorHistosOneYear(const std::vector<std::string>& files, const FlavorDef& fd,
                                      const std::string& sampleLabel, const std::string& year, bool isMC) {
    FlavorHistos out;
    if (files.empty()) return out;
    ROOT::RDataFrame frame("Events", files);

    auto regionCut = frame.Filter([catCode = fd.catCode](int cat) { return cat == catCode; }, {"cat"})
                          .Filter(
                              [](float metPt, double mt, double closestMlj) {
                                  if (metPt >= MET_MAX) return false;
                                  if (mt >= MT_MAX) return false;
                                  if (closestMlj > ZVETO_LOW && closestMlj < ZVETO_HIGH) return false;
                                  return true;
                              },
                              {"MET_pt", "mt", "closest_mlj"});

    auto base = regionCut.Filter(
        [btagWP = DEEPJET_MEDIUM_WP.at(year)](const ROOT::VecOps::RVec<float>& btag) {
            for (float v : btag) if (v > btagWP) return false;
            return true;
        },
        {"Jet_btagDeepFlavB"});

    auto withWeight = isMC ? base.Define("frWeight", "genWeight >= 0 ? 1.0 : -1.0")
                            : base.Define("frWeight", "1.0");

    const std::string tag = sampleLabel + "_" + fd.label + "_" + year;
    out.loose = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
        withWeight.Histo2D({("loose_" + tag).c_str(), "", (int)PT_EDGES.size() - 1, PT_EDGES.data(),
                             (int)ETA_EDGES.size() - 1, ETA_EDGES.data()}, "pt", "eta", "frWeight")
            .GetPtr()->Clone()));
    out.loose->SetDirectory(nullptr);
    out.tight = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
        withWeight.Filter("passTight")
            .Histo2D({("tight_" + tag).c_str(), "", (int)PT_EDGES.size() - 1, PT_EDGES.data(),
                      (int)ETA_EDGES.size() - 1, ETA_EDGES.data()}, "pt", "eta", "frWeight")
            .GetPtr()->Clone()));
    out.tight->SetDirectory(nullptr);

    if (isMC && fd.label == "electron") {
        auto realBase = withWeight.Filter("gen_match != 0");
        out.looseReal = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
            realBase.Histo2D({("looseReal_" + tag).c_str(), "", (int)PT_EDGES.size() - 1, PT_EDGES.data(),
                               (int)ETA_EDGES.size() - 1, ETA_EDGES.data()}, "pt", "eta", "frWeight")
                .GetPtr()->Clone()));
        out.looseReal->SetDirectory(nullptr);
        out.tightReal = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
            realBase.Filter("passTight")
                .Histo2D({("tightReal_" + tag).c_str(), "", (int)PT_EDGES.size() - 1, PT_EDGES.data(),
                          (int)ETA_EDGES.size() - 1, ETA_EDGES.data()}, "pt", "eta", "frWeight")
                .GetPtr()->Clone()));
        out.tightReal->SetDirectory(nullptr);
    }
    return out;
}

}  // namespace

void BuildQCDFakeRatesEtaScan() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    const std::string outDir = "results";
    gSystem->mkdir(outDir.c_str(), kTRUE);
    TFile out((outDir + "/QCD_fake_rates_EtaScan_Run2.root").c_str(), "RECREATE");

    const std::vector<std::pair<std::string, bool>> samples = {
        {"Data", false}, {"QCD_HT", true}, {"QCD_MuEnriched", true}, {"QCD_EMEnriched", true}
    };

    // [sample][flavor] combined across all 4 years.
    std::map<std::string, std::map<std::string, FlavorHistos>> combined;
    for (const auto& [sampleLabel, isMC] : samples)
        for (const auto& fd : FLAVORS)
            combined[sampleLabel][fd.label] = FlavorHistos{makeHist2D("comb_loose_" + sampleLabel + "_" + fd.label),
                                                             makeHist2D("comb_tight_" + sampleLabel + "_" + fd.label),
                                                             nullptr, nullptr};
    combined["QCD_EMEnriched"]["electron"].looseReal = makeHist2D("comb_looseReal_EM_electron");
    combined["QCD_EMEnriched"]["electron"].tightReal = makeHist2D("comb_tightReal_EM_electron");

    for (const auto& year : YEARS) {
        std::cout << "\n=== " << year << " ===" << std::endl;
        const std::string base = SKIM_BASE + "/" + year;
        for (const auto& [sampleLabel, isMC] : samples) {
            std::string pattern = (sampleLabel == "Data") ? (base + "/Data/*.root")
                                                            : (base + "/MC/" + sampleLabel + "/*.root");
            auto files = globFiles(pattern);
            std::cout << "  " << sampleLabel << ": " << files.size() << " file(s)" << std::endl;
            if (files.empty()) continue;
            for (const auto& fd : FLAVORS) {
                auto fh = fillFlavorHistosOneYear(files, fd, sampleLabel, year, isMC);
                if (!fh.loose) continue;
                addInPlace(*combined[sampleLabel][fd.label].loose, *fh.loose);
                addInPlace(*combined[sampleLabel][fd.label].tight, *fh.tight);
                if (isMC && fd.label == "electron" && fh.looseReal) {
                    addInPlace(*combined[sampleLabel][fd.label].looseReal, *fh.looseReal);
                    addInPlace(*combined[sampleLabel][fd.label].tightReal, *fh.tightReal);
                }
            }
        }
    }

    // Electron conversion correction, derived from the QCD_EMEnriched donor,
    // exactly as in BuildQCDFakeRates.C (see that file's header for why).
    auto& donor = combined[ELECTRON_CORRECTION_DONOR]["electron"];
    std::unique_ptr<TH2D> eleRealFracLoose, eleRealFracTight;
    if (donor.looseReal && donor.tightReal && donor.looseReal->Integral() > 0) {
        eleRealFracLoose = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(donor.looseReal->Clone("eleRealFracLoose")));
        eleRealFracLoose->Divide(donor.loose.get());
        eleRealFracTight = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(donor.tightReal->Clone("eleRealFracTight")));
        eleRealFracTight->Divide(donor.tight.get());
    } else {
        std::cerr << "[WARNING] no electron correction available, electron rates will be uncorrected." << std::endl;
    }

    std::map<std::string, std::unique_ptr<TH1D>> frVsEta;  // key: sampleLabel+"_"+flavor

    for (const auto& [sampleLabel, isMC] : samples) {
        for (const auto& fd : FLAVORS) {
            auto& fh = combined[sampleLabel][fd.label];
            if (fh.loose->Integral() <= 0) continue;

            auto looseClean = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(fh.loose->Clone()));
            auto tightClean = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(fh.tight->Clone()));

            if (fd.label == "electron" && eleRealFracLoose) {
                if (isMC && fh.looseReal) {
                    subtractInPlace(*looseClean, *fh.looseReal);
                    subtractInPlace(*tightClean, *fh.tightReal);
                } else if (!isMC) {
                    for (int bx = 1; bx <= looseClean->GetNbinsX(); ++bx)
                        for (int by = 1; by <= looseClean->GetNbinsY(); ++by) {
                            const double keepL = 1.0 - eleRealFracLoose->GetBinContent(bx, by);
                            const double keepT = 1.0 - eleRealFracTight->GetBinContent(bx, by);
                            looseClean->SetBinContent(bx, by, looseClean->GetBinContent(bx, by) * std::max(0.0, keepL));
                            tightClean->SetBinContent(bx, by, tightClean->GetBinContent(bx, by) * std::max(0.0, keepT));
                        }
                }
            }

            out.cd();
            looseClean->Write(("loose_clean_" + sampleLabel + "_" + fd.label).c_str());
            tightClean->Write(("tight_clean_" + sampleLabel + "_" + fd.label).c_str());

            // Marginalize over pT (ProjectionY sums along the x/pT axis),
            // then binomial-divide for the eta-only fake rate + proper errors.
            std::unique_ptr<TH1D> loose1D(tightClean->ProjectionY(("loose1D_" + sampleLabel + "_" + fd.label).c_str()));
            loose1D->Reset();
            loose1D->Add(looseClean->ProjectionY("_tmpL"));
            std::unique_ptr<TH1D> tight1D(tightClean->ProjectionY(("tight1D_" + sampleLabel + "_" + fd.label).c_str()));

            std::string key = sampleLabel + "_" + fd.label;
            frVsEta[key] = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(tight1D->Clone(("fr_eta_" + key).c_str())));
            frVsEta[key]->SetDirectory(nullptr);
            frVsEta[key]->Divide(tight1D.get(), loose1D.get(), 1.0, 1.0, "B");
            frVsEta[key]->SetTitle(";#eta;Fake rate");
            frVsEta[key]->Write();

            std::cout << "[" << key << "] N_loose=" << looseClean->Integral() << " N_tight=" << tightClean->Integral()
                       << " inclusive_rate=" << (looseClean->Integral() > 0 ? tightClean->Integral() / looseClean->Integral() : 0.0) << std::endl;
            for (int b = 1; b <= frVsEta[key]->GetNbinsX(); ++b) {
                std::cout << "    eta[" << frVsEta[key]->GetXaxis()->GetBinLowEdge(b) << ","
                          << frVsEta[key]->GetXaxis()->GetBinUpEdge(b) << "]: rate="
                          << frVsEta[key]->GetBinContent(b) << " +/- " << frVsEta[key]->GetBinError(b)
                          << "  (N_loose=" << loose1D->GetBinContent(b) << ")" << std::endl;
            }
        }
    }

    // Plot: Data-based electron and muon fake rate vs eta (the actual
    // QCD-derived rate used in the analysis), overlaid on one canvas.
    TCanvas c("c", "", 900, 700);
    c.SetLeftMargin(0.13);
    c.SetGridy();
    auto* hE = frVsEta.count("Data_electron") ? frVsEta["Data_electron"].get() : nullptr;
    auto* hM = frVsEta.count("Data_muon") ? frVsEta["Data_muon"].get() : nullptr;
    if (hE) {
        hE->SetLineColor(kRed + 1);
        hE->SetMarkerColor(kRed + 1);
        hE->SetMarkerStyle(20);
        hE->GetYaxis()->SetTitle("Fake rate");
        hE->GetXaxis()->SetTitle("#eta");
        double ymax = 0;
        for (int b = 1; b <= hE->GetNbinsX(); ++b) ymax = std::max(ymax, hE->GetBinContent(b) + hE->GetBinError(b));
        if (hM) for (int b = 1; b <= hM->GetNbinsX(); ++b) ymax = std::max(ymax, hM->GetBinContent(b) + hM->GetBinError(b));
        hE->GetYaxis()->SetRangeUser(0, ymax * 1.4);
        hE->Draw("E1");
    }
    if (hM) {
        hM->SetLineColor(kAzure + 2);
        hM->SetMarkerColor(kAzure + 2);
        hM->SetMarkerStyle(21);
        hM->Draw(hE ? "E1 SAME" : "E1");
    }
    TLegend leg(0.6, 0.75, 0.88, 0.88);
    leg.SetBorderSize(0);
    if (hE) leg.AddEntry(hE, "Electron (Data, QCD-enriched CR)", "lep");
    if (hM) leg.AddEntry(hM, "Muon (Data, QCD-enriched CR)", "lep");
    leg.Draw();
    gSystem->mkdir((outDir + "/plots").c_str(), kTRUE);
    c.Print((outDir + "/plots/QCD_fake_rate_vs_eta_Run2.png").c_str());

    out.Close();
    std::cout << "\nWrote " << outDir << "/QCD_fake_rates_EtaScan_Run2.root" << std::endl;
    std::cout << "Wrote " << outDir << "/plots/QCD_fake_rate_vs_eta_Run2.png" << std::endl;
}
