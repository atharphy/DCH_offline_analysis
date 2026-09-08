// Builds jet->lepton fake rates f(pT,|eta|) = N_tight/N_loose for electrons,
// muons, and taus (cat 1/2/3) from the QCD-enriched-region skims produced by
// Online_framework/DCH/qcd_fakerate_skim.py.
//
// IMPORTANT: the skim itself only gates on "one loose object" + "an away
// jet exists" -- MET/mT/Z-veto/trigger/b-tag are all stored as branches,
// not applied there (see that file's docstring for why: so region cuts can
// be retuned without re-skimming ~10k condor jobs each time). This script
// is where the QCD-enriched region actually gets defined: MET_MAX/MT_MAX/
// ZVETO_LOW/ZVETO_HIGH below (AN-19-111/X53-AN2016 baseline values) are
// always applied; the b-veto (DEEPJET_MEDIUM_WP) is applied by default but
// toggleable via the applyBVeto argument, since it's our own addition on
// top of the AN's own region, not part of the literature baseline.
//
// Run on the merged files from merge_qcd_fakerate_skims.py, one call per
// year, producing Data_<flavor>_fake_rate and
// <QCD_HT|QCD_MuEnriched|QCD_EMEnriched>_<flavor>_fake_rate 2D histograms
// (pT x |eta|) in a single output file per year.
//
// tight is a strict subset of loose by construction, so N_tight/N_loose gets
// binomial ("B") error propagation via TH1::Divide, same convention as
// DeriveDYTauFakeRates.C's makeFakeRate(). That subset property survives the
// electron correction below too (see there).
//
// QCD MC categories (QCD_HT/QCD_MuEnriched/QCD_EMEnriched) are built and
// stored SEPARATELY, never summed -- different filter efficiencies. Pick
// whichever category is the appropriate f_MC comparison for a given flavor
// at the point of use.
//
// --- electron conversion/prompt correction ---
// Diagnosed directly against QCD_HT 2018: 29.5% of loose electrons in pure
// QCD MC are gen-matched to something real (mostly genPartFlav==22, i.e.
// photon-conversion electrons -- pi0->gg->e+e- inside jets), and that real
// population passes tight ID at ~59% vs ~15% for genuine fakes, since
// conversions are real electron tracks and pass MVA ID/isolation far more
// readily than a jet fragment does. Checked directly: convVeto and lostHits
// are already maxed out for 100% of these (no reco-level cut headroom left
// -- see the conversation this was diagnosed in), so this cannot be fixed by
// tightening object selection. Muons show no analogous gap (0.216 fake vs
// 0.224 real pass rate -- statistically indistinguishable, no correction
// applied there).
//
// Matches AN-19-111's own handling exactly: measure the raw rate, then
// apply a *derived, multiplicative* correction rather than a truth-level
// cut, so the correction is usable on data (which has no gen truth) too.
// For MC, the correction is self-derived per sample from that sample's own
// gen_match branch (genPartFlav != 0 => "real", subtracted from both tight
// and loose bin counts). For Data, gen_match doesn't exist, so the
// correction fraction is instead borrowed from QCD_EMEnriched (the
// electron-flavor-filtered MC, the most statistically appropriate donor)
// and applied as a per-bin multiplicative scale-down of Data's raw counts.
// Subtracting the same (gen-matched-real) population from both tight and
// loose preserves tight-clean subset-of-loose-clean, so "B" division stays
// valid on the corrected histograms too.
//
// Run: root -l -b -q BuildQCDFakeRates.C+'("2018")'
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

const std::string MERGED_BASE = "/eos/user/a/atahmad/DCH_offline_analysis/QCD_FakeRate_merged";
const std::string ELECTRON_CORRECTION_DONOR = "QCD_EMEnriched";

// The skim itself only gates on "one loose object" + "an away jet exists"
// (see qcd_fakerate_skim.py's docstring) -- everything else needed to
// actually define the QCD-enriched region is stored as branches and
// applied HERE instead, so it can be retuned without re-skimming. These
// four are the AN-19-111/X53-AN2016 baseline values.
const double MET_MAX = 25.0;
const double MT_MAX = 25.0;
const double ZVETO_LOW = 81.1;
const double ZVETO_HIGH = 101.1;

// BTV POG-recommended DeepJet (DeepFlavB) Medium WPs, per UL year. Not part
// of the AN's own region definition -- added after diagnosing that
// semileptonic heavy-flavor decay (b/c->l) is the dominant source of "real"
// (non-fake) contamination in the loose sample (see the electron
// conversion-correction note below for the analogous electron-side issue).
// Toggle via BuildQCDFakeRates's applyBVeto argument.
const std::map<std::string, double> DEEPJET_MEDIUM_WP = {
    {"2016preVFP", 0.2598}, {"2016postVFP", 0.2489}, {"2017", 0.3040}, {"2018", 0.2783},
};

struct FlavorBinning {
    std::string label;   // "electron", "muon", "tau"
    int catCode;          // 1, 2, 3
    std::vector<double> ptEdges;
    std::vector<double> etaEdges;  // |eta|
};

// pT binning kept moderately coarse by default given this is a QCD-enriched
// control region (limited statistics, especially for data electrons/muons
// at high pT) -- retune once actual per-bin yields from a full year are in
// hand. |eta| splits follow each object's own standard barrel/endcap
// convention: ECAL crack sits at the electron split (1.479, the usual
// EGamma POG boundary -- edges are fine there since it's the crack sitting
// *between* bins, not a bin center inside it, which is the failure mode).
const std::vector<FlavorBinning> FLAVORS = {
    {"electron", 1, {20, 30, 45, 70, 100, 200}, {0.0, 1.479, 2.5}},
    {"muon",     2, {20, 30, 45, 70, 100, 200}, {0.0, 1.2, 2.4}},
    {"tau",      3, {20, 30, 45, 70, 100, 200}, {0.0, 1.5, 2.3}},
};

std::unique_ptr<TH2D> makeHist2D(const std::string& name, const FlavorBinning& fb) {
    auto h = std::make_unique<TH2D>(name.c_str(), "",
        static_cast<int>(fb.ptEdges.size()) - 1, fb.ptEdges.data(),
        static_cast<int>(fb.etaEdges.size()) - 1, fb.etaEdges.data());
    h->SetDirectory(nullptr);
    h->Sumw2();
    return h;
}

// Divides out(bin) by "correction stays >=0" safety: never lets a bin go
// negative from the subtraction (would happen only from a statistical
// fluctuation in a very low-stat bin, not a real physics case).
void subtractInPlace(TH2D& target, const TH2D& toSubtract) {
    for (int bx = 1; bx <= target.GetNbinsX(); ++bx) {
        for (int by = 1; by <= target.GetNbinsY(); ++by) {
            const double v = std::max(0.0, target.GetBinContent(bx, by) - toSubtract.GetBinContent(bx, by));
            target.SetBinContent(bx, by, v);
            // conservative error: keep the original (uncorrected) uncertainty
            // rather than trying to propagate a subtraction of two
            // correlated-in-origin MC-derived quantities precisely.
        }
    }
}

std::unique_ptr<TH2D> divideBinomial(const TH2D& tight, const TH2D& loose, const std::string& name) {
    std::unique_ptr<TH2D> fr(dynamic_cast<TH2D*>(tight.Clone(name.c_str())));
    fr->SetDirectory(nullptr);
    fr->Divide(&tight, &loose, 1.0, 1.0, "B");
    fr->SetTitle(";p_{T} [GeV];|#eta|;Fake rate");
    fr->GetZaxis()->SetTitle("Fake rate");
    return fr;
}

struct FlavorHistos {
    std::unique_ptr<TH2D> loose;
    std::unique_ptr<TH2D> tight;
    std::unique_ptr<TH2D> looseReal;  // MC only: gen_match != 0 subset
    std::unique_ptr<TH2D> tightReal;  // MC only
};

FlavorHistos fillFlavorHistos(
    ROOT::RDataFrame& frame,
    const FlavorBinning& fb,
    const std::string& sampleLabel,
    const std::string& year,
    bool isMC,
    bool applyBVeto
) {
    auto regionCut = frame.Filter(
        [catCode = fb.catCode](int cat) { return cat == catCode; }, {"cat"}
    ).Filter(
        [](double metPt, double mt, double closestMlj) {
            if (metPt >= MET_MAX) return false;
            if (mt >= MT_MAX) return false;
            if (closestMlj > ZVETO_LOW && closestMlj < ZVETO_HIGH) return false;
            return true;
        },
        {"MET_pt", "mt", "closest_mlj"}
    );

    // Always inserts one Filter node (RDataFrame's node type changes with
    // each .Filter()/.Define() call, so a ternary between "filtered" and
    // "unfiltered" branches won't compile) -- made a no-op when applyBVeto
    // is false by just short-circuiting inside the lambda instead.
    auto base = regionCut.Filter(
        [applyBVeto, btagWP = DEEPJET_MEDIUM_WP.at(year)](const ROOT::VecOps::RVec<float>& btag) {
            if (!applyBVeto) return true;
            for (float v : btag) if (v > btagWP) return false;
            return true;
        },
        {"Jet_btagDeepFlavB"}
    ).Define("absEta", "std::fabs(eta)");

    auto withWeight = isMC
        ? base.Define("frWeight", "genWeight >= 0 ? 1.0 : -1.0")
        : base.Define("frWeight", "1.0");

    const std::string tag = sampleLabel + "_" + fb.label + "_" + year;
    FlavorHistos out;
    out.loose = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
        withWeight.Histo2D({("loose_" + tag).c_str(), "", static_cast<int>(fb.ptEdges.size()) - 1,
                             fb.ptEdges.data(), static_cast<int>(fb.etaEdges.size()) - 1, fb.etaEdges.data()},
                            "pt", "absEta", "frWeight")
            .GetPtr()
            ->Clone()));
    out.loose->SetDirectory(nullptr);
    out.tight = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
        withWeight.Filter("passTight")
            .Histo2D({("tight_" + tag).c_str(), "", static_cast<int>(fb.ptEdges.size()) - 1, fb.ptEdges.data(),
                      static_cast<int>(fb.etaEdges.size()) - 1, fb.etaEdges.data()},
                     "pt", "absEta", "frWeight")
            .GetPtr()
            ->Clone()));
    out.tight->SetDirectory(nullptr);

    if (isMC && fb.label == "electron") {
        auto realBase = withWeight.Filter("gen_match != 0");
        out.looseReal = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
            realBase
                .Histo2D({("looseReal_" + tag).c_str(), "", static_cast<int>(fb.ptEdges.size()) - 1,
                          fb.ptEdges.data(), static_cast<int>(fb.etaEdges.size()) - 1, fb.etaEdges.data()},
                         "pt", "absEta", "frWeight")
                .GetPtr()
                ->Clone()));
        out.looseReal->SetDirectory(nullptr);
        out.tightReal = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(
            realBase.Filter("passTight")
                .Histo2D({("tightReal_" + tag).c_str(), "", static_cast<int>(fb.ptEdges.size()) - 1,
                          fb.ptEdges.data(), static_cast<int>(fb.etaEdges.size()) - 1, fb.etaEdges.data()},
                         "pt", "absEta", "frWeight")
                .GetPtr()
                ->Clone()));
        out.tightReal->SetDirectory(nullptr);
    }

    return out;
}

}  // namespace

void BuildQCDFakeRates(std::string year = "2018", bool applyBVeto = true) {
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);

    const std::string outDir = "results";
    gSystem->mkdir(outDir.c_str(), kTRUE);
    const std::string outName = outDir + "/QCD_fake_rates_" + year + ".root";
    TFile out(outName.c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "ERROR: cannot create " << outName << std::endl;
        gSystem->Exit(1);
        return;
    }

    const std::string base = MERGED_BASE + "/" + year;
    const std::vector<std::pair<std::string, bool>> samples = {
        {"Data", false}, {"QCD_HT", true}, {"QCD_MuEnriched", true}, {"QCD_EMEnriched", true}
    };

    // Pass 1: fill every flavor's histograms for every sample, and along the
    // way capture the electron correction donor's real-fraction (as a
    // fraction of that donor's own loose/tight, so it can be reapplied
    // multiplicatively to Data's raw counts, which have no gen_match of
    // their own).
    std::map<std::string, std::map<std::string, FlavorHistos>> histos;  // [sample][flavor]
    for (const auto& [sampleLabel, isMC] : samples) {
        const std::string filePath = base + "/" + sampleLabel + "_" + year + ".root";
        if (gSystem->AccessPathName(filePath.c_str())) {
            std::cerr << "[WARNING] missing merged file, skipping: " << filePath << std::endl;
            continue;
        }
        ROOT::RDataFrame frame("Events", filePath);
        for (const auto& fb : FLAVORS) {
            histos[sampleLabel][fb.label] = fillFlavorHistos(frame, fb, sampleLabel, year, isMC, applyBVeto);
        }
    }

    const FlavorBinning* electronBinning = nullptr;
    for (const auto& fb : FLAVORS)
        if (fb.label == "electron") electronBinning = &fb;

    // Electron real-fraction from the donor MC, as a per-bin fraction (not
    // an absolute count), so it can be rescaled onto Data's own statistics.
    std::unique_ptr<TH2D> eleRealFracLoose, eleRealFracTight;
    if (electronBinning && histos.count(ELECTRON_CORRECTION_DONOR) &&
        histos[ELECTRON_CORRECTION_DONOR].count("electron")) {
        auto& donor = histos[ELECTRON_CORRECTION_DONOR]["electron"];
        if (donor.looseReal && donor.tightReal) {
            eleRealFracLoose = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(donor.looseReal->Clone("eleRealFracLoose")));
            eleRealFracLoose->SetDirectory(nullptr);
            eleRealFracLoose->Divide(donor.loose.get());
            eleRealFracTight = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(donor.tightReal->Clone("eleRealFracTight")));
            eleRealFracTight->SetDirectory(nullptr);
            eleRealFracTight->Divide(donor.tight.get());
            std::cout << "[electron correction] derived from " << ELECTRON_CORRECTION_DONOR
                      << " " << year << ":" << std::endl;
            for (int bx = 1; bx <= eleRealFracLoose->GetNbinsX(); ++bx) {
                for (int by = 1; by <= eleRealFracLoose->GetNbinsY(); ++by) {
                    std::cout << "  pt[" << eleRealFracLoose->GetXaxis()->GetBinLowEdge(bx) << ","
                              << eleRealFracLoose->GetXaxis()->GetBinUpEdge(bx) << "] eta bin " << by
                              << ": loose real-frac=" << eleRealFracLoose->GetBinContent(bx, by)
                              << " tight real-frac=" << eleRealFracTight->GetBinContent(bx, by) << std::endl;
                }
            }
        }
    }
    if (!eleRealFracLoose) {
        std::cerr << "[WARNING] could not derive electron conversion correction (missing "
                  << ELECTRON_CORRECTION_DONOR << "), electron rates will be uncorrected." << std::endl;
    }

    // Pass 2: apply correction (electron only) and write everything out.
    for (const auto& [sampleLabel, isMC] : samples) {
        if (!histos.count(sampleLabel)) continue;
        for (const auto& fb : FLAVORS) {
            if (!histos[sampleLabel].count(fb.label)) continue;
            auto& fh = histos[sampleLabel][fb.label];

            auto looseClean = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(fh.loose->Clone()));
            looseClean->SetDirectory(nullptr);
            auto tightClean = std::unique_ptr<TH2D>(dynamic_cast<TH2D*>(fh.tight->Clone()));
            tightClean->SetDirectory(nullptr);

            if (fb.label == "electron" && eleRealFracLoose) {
                if (isMC && fh.looseReal && fh.tightReal) {
                    // self-derived: subtract this sample's OWN gen-matched-real counts
                    subtractInPlace(*looseClean, *fh.looseReal);
                    subtractInPlace(*tightClean, *fh.tightReal);
                } else if (!isMC) {
                    // Data: no gen_match of its own -- scale down by the donor's
                    // real-fraction instead of subtracting an absolute count.
                    for (int bx = 1; bx <= looseClean->GetNbinsX(); ++bx) {
                        for (int by = 1; by <= looseClean->GetNbinsY(); ++by) {
                            const double keepLoose = 1.0 - eleRealFracLoose->GetBinContent(bx, by);
                            const double keepTight = 1.0 - eleRealFracTight->GetBinContent(bx, by);
                            looseClean->SetBinContent(bx, by, looseClean->GetBinContent(bx, by) * std::max(0.0, keepLoose));
                            tightClean->SetBinContent(bx, by, tightClean->GetBinContent(bx, by) * std::max(0.0, keepTight));
                        }
                    }
                }
            }

            const double nLooseRaw = fh.loose->Integral();
            const double nTightRaw = fh.tight->Integral();
            const double nLooseClean = looseClean->Integral();
            const double nTightClean = tightClean->Integral();
            std::cout << "[" << sampleLabel << " " << year << " " << fb.label << "] "
                      << "raw: N_loose=" << nLooseRaw << " N_tight=" << nTightRaw
                      << " rate=" << (nLooseRaw > 0 ? nTightRaw / nLooseRaw : 0.0)
                      << "  |  corrected: N_loose=" << nLooseClean << " N_tight=" << nTightClean
                      << " rate=" << (nLooseClean > 0 ? nTightClean / nLooseClean : 0.0) << std::endl;

            out.cd();
            fh.loose->Write(("raw_loose_" + sampleLabel + "_" + fb.label + "_" + year).c_str());
            fh.tight->Write(("raw_tight_" + sampleLabel + "_" + fb.label + "_" + year).c_str());
            auto frRaw = divideBinomial(*fh.tight, *fh.loose, sampleLabel + "_" + fb.label + "_fake_rate_raw");
            frRaw->Write();
            auto frClean = divideBinomial(*tightClean, *looseClean, sampleLabel + "_" + fb.label + "_fake_rate");
            frClean->Write();
        }
    }

    out.Close();
    std::cout << "Wrote " << outName << std::endl;
}
