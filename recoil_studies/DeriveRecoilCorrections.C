#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <correction.h>

#if defined(__CLING__)
R__LOAD_LIBRARY(libcorrectionlib.so)
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../filemap/FileMap.h"
#include "../include/Xsections.C"

namespace {

const bool APPLY_ZPT_REWEIGHTING = true;
const std::string OUTPUT_BASE =
    "/eos/user/a/atahmad/DCH_offline_analysis/recoil_studies";

const std::vector<double> ZPT_EDGES = {
    0., 10., 20., 30., 50., 1000.
};

enum SampleGroup { DATA = 0, BACKGROUND = 1, DY = 2, N_GROUPS = 3 };
enum Component { PARALLEL = 0, PERPENDICULAR = 1, N_COMPONENTS = 2 };

struct Moment {
    double sumW = 0.0;
    double sumWX = 0.0;
    double sumWX2 = 0.0;
    double sumW2 = 0.0;

    void fill(double x, double weight) {
        sumW += weight;
        sumWX += weight * x;
        sumWX2 += weight * x * x;
        sumW2 += weight * weight;
    }
};

struct Estimate {
    bool valid = false;
    double yield = 0.0;
    double mean = 0.0;
    double sigma = 0.0;
    double meanError = 0.0;
    double sigmaError = 0.0;
    double effectiveEntries = 0.0;
};

struct OfficialMETCorrections {
    std::unique_ptr<correction::CorrectionSet> set;
    correction::Correction::Ref pt;
    correction::Correction::Ref phi;
};

struct Branches {
    Int_t cat = 0;
    Int_t nPV = 0;
    Double_t njets = 0.0;
    ULong64_t run = 0;
    Double_t met = 0.0;
    Double_t metphi = 0.0;
    Double_t pt[2] = {};
    Double_t eta[2] = {};
    Double_t phi[2] = {};
    Double_t mass[2] = {};
    Double_t charge[2] = {};
    std::vector<double>* genPartPt = nullptr;
    std::vector<double>* genPartPhi = nullptr;
    std::vector<int>* genPartPdgId = nullptr;
    std::vector<int>* genPartStatus = nullptr;
    std::vector<int>* genPartStatusFlags = nullptr;
    Double_t idSF[2] = {1.0, 1.0};
    Double_t isoSF[2] = {1.0, 1.0};
    Double_t trigSF[2] = {1.0, 1.0};
    Double_t isTrig1 = 0.0;
    Double_t generatorWeight = 1.0;
    Double_t brWeight = 1.0;
    Double_t prefiringWeight = 1.0;
    Double_t pileupWeight = 1.0;
};

struct Accumulator {
    std::vector<std::vector<std::vector<Moment>>> moments;
    std::vector<std::vector<std::vector<TH1D*>>> distributions;

    Accumulator() {
        const int nPtBins = static_cast<int>(ZPT_EDGES.size()) - 1;
        moments.resize(N_GROUPS);
        distributions.resize(N_GROUPS);

        for (int group = 0; group < N_GROUPS; ++group) {
            moments[group].resize(N_COMPONENTS);
            distributions[group].resize(N_COMPONENTS);
            for (int component = 0; component < N_COMPONENTS; ++component) {
                moments[group][component].resize(3 * nPtBins);
                distributions[group][component].resize(3 * nPtBins, nullptr);
            }
        }
    }

    ~Accumulator() {
        for (auto& byGroup : distributions) {
            for (auto& byComponent : byGroup) {
                for (TH1D* histogram : byComponent) delete histogram;
            }
        }
    }
};

double luminosity(const std::string& year) {
    if (year == "2016preVFP") return 19520.0;
    if (year == "2016postVFP") return 16810.0;
    if (year == "2017") return 41480.0;
    if (year == "2018") return 59830.0;
    return 0.0;
}

std::string categoryString(Int_t category) {
    if (category == 40) return "ee";
    if (category == 43) return "mm";
    return "";
}

// Follows the CMS-HTT RecoilCorrections README exactly (same algorithm as
// DCH_modules/GenBosonMomentum.h's getGenBosonMomentum, self-contained here
// since this script binds its own branches rather than sharing MyBranch.C's
// globals):
//   full:    (fromHardProcessFinalState && (isMuon||isElectron||isNeutrino)) || isDirectHardProcessTauDecayProduct
//   visible: (fromHardProcessFinalState && (isMuon||isElectron)) || (isDirectHardProcessTauDecayProduct && !isNeutrino)
// fromHardProcessFinalState = status==1 && statusFlags bit 8 (fromHardProcess).
// isDirectHardProcessTauDecayProduct = statusFlags bit 10.
bool findGenBosonAxis(
    const Branches& b,
    double& fullPx,
    double& fullPy,
    double& visPx,
    double& visPy
) {
    if (!b.genPartPdgId || !b.genPartStatus || !b.genPartStatusFlags ||
        !b.genPartPt || !b.genPartPhi) return false;

    fullPx = fullPy = visPx = visPy = 0.0;
    bool found = false;

    for (size_t i = 0; i < b.genPartPdgId->size(); ++i) {
        const int pdg = std::abs((*b.genPartPdgId)[i]);
        const bool isLepton = (pdg == 11 || pdg == 13);
        const bool isNeutrino = (pdg == 12 || pdg == 14 || pdg == 16);
        if (!isLepton && !isNeutrino) continue;

        const bool fromHardProcessFinalState =
            (*b.genPartStatus)[i] == 1 && (((*b.genPartStatusFlags)[i] >> 8) & 1);
        const bool isDirectHardProcessTauDecayProduct =
            (((*b.genPartStatusFlags)[i] >> 10) & 1);

        if (!((fromHardProcessFinalState && (isLepton || isNeutrino)) || isDirectHardProcessTauDecayProduct)) continue;

        const double px = (*b.genPartPt)[i] * std::cos((*b.genPartPhi)[i]);
        const double py = (*b.genPartPt)[i] * std::sin((*b.genPartPhi)[i]);
        fullPx += px;
        fullPy += py;
        found = true;

        if ((fromHardProcessFinalState && isLepton) || (isDirectHardProcessTauDecayProduct && !isNeutrino)) {
            visPx += px;
            visPy += py;
        }
    }
    return found;
}

int ptBin(double pt) {
    if (!std::isfinite(pt) || pt < ZPT_EDGES.front() || pt >= ZPT_EDGES.back()) {
        return -1;
    }
    const auto upper = std::upper_bound(ZPT_EDGES.begin(), ZPT_EDGES.end(), pt);
    return static_cast<int>(upper - ZPT_EDGES.begin()) - 1;
}

int jetBin(double njets) {
    if (!std::isfinite(njets) || njets < 0.0) return -1;
    if (njets < 0.5) return 0;
    if (njets < 1.5) return 1;
    return 2;
}

int flatBin(int jet, int pt) {
    return jet * (static_cast<int>(ZPT_EDGES.size()) - 1) + pt;
}

TH1D* getDistribution(
    Accumulator& accumulator,
    int group,
    int component,
    int jet,
    int pt
) {
    TH1D*& histogram =
        accumulator.distributions[group][component][flatBin(jet, pt)];
    if (histogram) return histogram;

    const char* groupName[] = {"data", "background", "dy"};
    const char* componentName[] = {"upar", "uperp"};
    const double xmin = -400.0;
    const double xmax = 400.0;

    histogram = new TH1D(
        Form("h_%s_%s_njet%d_pt%d",
             groupName[group], componentName[component], jet, pt),
        "",
        320,
        xmin,
        xmax
    );
    histogram->Sumw2();
    histogram->SetDirectory(nullptr);
    return histogram;
}

void fillAccumulator(
    Accumulator& accumulator,
    int group,
    int jet,
    int pt,
    double uParallel,
    double uPerpendicular,
    double weight
) {
    const int index = flatBin(jet, pt);
    accumulator.moments[group][PARALLEL][index].fill(uParallel, weight);
    accumulator.moments[group][PERPENDICULAR][index].fill(uPerpendicular, weight);
    getDistribution(accumulator, group, PARALLEL, jet, pt)
        ->Fill(uParallel, weight);
    getDistribution(accumulator, group, PERPENDICULAR, jet, pt)
        ->Fill(uPerpendicular, weight);
}

Estimate estimate(const Moment& moment) {
    Estimate result;
    result.yield = moment.sumW;

    if (!(moment.sumW > 0.0) || !(moment.sumW2 > 0.0)) return result;

    result.mean = moment.sumWX / moment.sumW;
    const double variance = moment.sumWX2 / moment.sumW - result.mean * result.mean;
    if (!(variance > 0.0) || !std::isfinite(variance)) return result;

    result.sigma = std::sqrt(variance);
    result.effectiveEntries = moment.sumW * moment.sumW / moment.sumW2;
    if (!(result.effectiveEntries > 2.0)) return result;

    result.meanError = result.sigma / std::sqrt(result.effectiveEntries);
    result.sigmaError = result.sigma /
        std::sqrt(2.0 * (result.effectiveEntries - 1.0));
    result.valid = std::isfinite(result.mean) && std::isfinite(result.sigma);
    return result;
}

Moment subtractBackground(const Moment& data, const Moment& background) {
    Moment result;
    result.sumW = data.sumW - background.sumW;
    result.sumWX = data.sumWX - background.sumWX;
    result.sumWX2 = data.sumWX2 - background.sumWX2;
    result.sumW2 = data.sumW2 + background.sumW2;
    return result;
}

std::unique_ptr<OfficialMETCorrections> loadMETCorrections(
    const std::string& year,
    bool isData
) {
    const std::string fileName = "json_files/met_" + year + ".json.gz";
    const std::string suffix = isData ? "data" : "mc";

    try {
        std::unique_ptr<OfficialMETCorrections> result(new OfficialMETCorrections());
        result->set = correction::CorrectionSet::from_file(fileName);
        result->pt = result->set->at("pt_metphicorr_pfmet_" + suffix);
        result->phi = result->set->at("phi_metphicorr_pfmet_" + suffix);
        return result;
    }
    catch (const std::exception& error) {
        std::cerr << "ERROR loading " << fileName << ": " << error.what() << std::endl;
        return nullptr;
    }
}

bool correctedMET(
    const OfficialMETCorrections& correction,
    const Branches& b,
    double& correctedPt,
    double& correctedPhi
) {
    if (!std::isfinite(b.met) || !std::isfinite(b.metphi) ||
        b.met < 0.0 || b.nPV < 0) return false;

    try {
        correctedPt = correction.pt->evaluate({
            b.met,
            b.metphi,
            static_cast<double>(b.nPV),
            static_cast<double>(b.run)
        });
        correctedPhi = correction.phi->evaluate({
            b.met,
            b.metphi,
            static_cast<double>(b.nPV),
            static_cast<double>(b.run)
        });
        return std::isfinite(correctedPt) && std::isfinite(correctedPhi) &&
               correctedPt >= 0.0;
    }
    catch (const std::exception&) {
        return false;
    }
}

std::unique_ptr<TH1D> loadZPtWeights(const std::string& year) {
    const std::string fileName =
        "/eos/user/a/atahmad/DCH_offline_analysis/zpt_studies/" +
        year + "/ZPtWeights_" + year + ".root";
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    if (!file || file->IsZombie()) return nullptr;

    TH1D* source = dynamic_cast<TH1D*>(file->Get("h_zpt_weight_ee"));
    if (!source) return nullptr;

    std::unique_ptr<TH1D> result(
        dynamic_cast<TH1D*>(source->Clone("h_zpt_weights_for_recoil"))
    );
    result->SetDirectory(nullptr);
    return result;
}

double zPtWeight(const TH1D* histogram, double pt) {
    if (!histogram || !std::isfinite(pt)) return 1.0;
    int bin = histogram->GetXaxis()->FindFixBin(pt);
    bin = std::max(1, std::min(bin, histogram->GetNbinsX()));
    const double weight = histogram->GetBinContent(bin);
    return std::isfinite(weight) && weight >= 0.0 ? weight : 1.0;
}

bool bindBranches(TTree* tree, Branches& b, bool isData) {
    tree->SetBranchStatus("*", 0);

    auto bind = [tree](const char* name, void* address) {
        if (!tree->GetBranch(name)) {
            std::cerr << "ERROR: missing branch " << name << std::endl;
            return false;
        }
        tree->SetBranchStatus(name, 1);
        tree->SetBranchAddress(name, address);
        return true;
    };

    bool valid = true;
    valid &= bind("cat", &b.cat);
    valid &= bind("run", &b.run);
    valid &= bind("nPV", &b.nPV);
    valid &= bind("njets", &b.njets);
    valid &= bind("met", &b.met);
    valid &= bind("metphi", &b.metphi);

    for (int index = 0; index < 2; ++index) {
        const std::string suffix = "_" + std::to_string(index + 1);
        valid &= bind(("pt" + suffix).c_str(), &b.pt[index]);
        valid &= bind(("eta" + suffix).c_str(), &b.eta[index]);
        valid &= bind(("phi" + suffix).c_str(), &b.phi[index]);
        valid &= bind(("m" + suffix).c_str(), &b.mass[index]);
        valid &= bind(("q" + suffix).c_str(), &b.charge[index]);
    }

    if (!isData) {
        valid &= bind("Generator_weight", &b.generatorWeight);
        valid &= bind("brWeight", &b.brWeight);
        valid &= bind("L1PreFiringWeight_Nom", &b.prefiringWeight);
        valid &= bind("weightPUtruejson", &b.pileupWeight);
        valid &= bind("isTrig_1", &b.isTrig1);
        valid &= bind("GenPart_pt", &b.genPartPt);
        valid &= bind("GenPart_phi", &b.genPartPhi);
        valid &= bind("GenPart_pdgId", &b.genPartPdgId);
        valid &= bind("GenPart_status", &b.genPartStatus);
        valid &= bind("GenPart_statusFlags", &b.genPartStatusFlags);

        for (int index = 0; index < 2; ++index) {
            const std::string suffix = "_" + std::to_string(index + 1);
            valid &= bind(("IDSF" + suffix).c_str(), &b.idSF[index]);
            valid &= bind(("ISOSF" + suffix).c_str(), &b.isoSF[index]);
            valid &= bind(("TrigSF" + suffix).c_str(), &b.trigSF[index]);
        }
    }

    tree->SetCacheSize(128LL * 1024LL * 1024LL);
    tree->SetCacheLearnEntries(20);
    tree->AddBranchToCache("*", kTRUE);
    return valid;
}

double nominalWeight(const Branches& b, double crossSectionWeight) {
    double weight = b.generatorWeight * b.brWeight * crossSectionWeight;
    weight *= b.prefiringWeight * b.pileupWeight;
    weight *= b.idSF[0] * b.isoSF[0] * b.idSF[1] * b.isoSF[1];

    double triggerWeight = 1.0;
    if (b.isTrig1 >= 1.0) triggerWeight = b.trigSF[0];
    else if (b.isTrig1 == -1.0) triggerWeight = b.trigSF[1];
    return weight * triggerWeight;
}

bool acceptsDataChannel(const std::string& baseName, const std::string& category) {
    if (category == "ee") {
        return baseName.find("SingleElectron") != std::string::npos ||
               baseName.find("EGamma") != std::string::npos;
    }
    if (category == "mm") {
        return baseName.find("SingleMuon") != std::string::npos;
    }
    return false;
}

int processFile(
    const std::string& year,
    const std::string& process,
    const std::string& fileName,
    double lumi,
    Accumulator& accumulator,
    const OfficialMETCorrections& dataMET,
    const OfficialMETCorrections& mcMET,
    const TH1D* zPtWeights
) {
    std::unique_ptr<TFile> inputFile(TFile::Open(fileName.c_str(), "READ"));
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "[skip] cannot open " << fileName << std::endl;
        return 1;
    }

    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if (!tree) {
        std::cerr << "[skip] no Events tree in " << fileName << std::endl;
        return 2;
    }

    const bool isData = process == "data";
    const bool isDY = process == "DY";
    const int group = isData ? DATA : (isDY ? DY : BACKGROUND);
    const std::string baseName = gSystem->BaseName(fileName.c_str());

    Branches b;
    if (!bindBranches(tree, b, isData)) return 3;

    double crossSectionWeight = 1.0;
    if (!isData) {
        TH1D* eventCounter = dynamic_cast<TH1D*>(inputFile->Get("hNWEvts"));
        if (!eventCounter) eventCounter = dynamic_cast<TH1D*>(inputFile->Get("hNEvts"));
        const double denominator = eventCounter ? eventCounter->Integral() : 0.0;
        if (!(denominator > 0.0)) return 4;
        crossSectionWeight = lumi * XSec(baseName) / denominator;
    }

    Long64_t selected = 0;
    const Long64_t entries = tree->GetEntriesFast();

    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);

        const std::string category = categoryString(b.cat);
        if (category.empty()) continue;
        if (isData && !acceptsDataChannel(baseName, category)) continue;
        if (!(b.pt[0] > 0.0) || !(b.pt[1] > 0.0)) continue;
        if (!(b.charge[0] * b.charge[1] < 0.0)) continue;

        TLorentzVector firstLepton;
        TLorentzVector secondLepton;
        firstLepton.SetPtEtaPhiM(b.pt[0], b.eta[0], b.phi[0], b.mass[0]);
        secondLepton.SetPtEtaPhiM(b.pt[1], b.eta[1], b.phi[1], b.mass[1]);
        const TLorentzVector boson = firstLepton + secondLepton;

        if (!std::isfinite(boson.M()) || std::fabs(boson.M() - 91.2) > 10.0) continue;

        double axisPx = boson.Px();
        double axisPy = boson.Py();
        double axisPt = boson.Pt();
        double genPx = axisPx;
        double genPy = axisPy;
        double visPx = axisPx;
        double visPy = axisPy;

        if (isDY) {
            if (!findGenBosonAxis(b, genPx, genPy, visPx, visPy)) continue;
            axisPx = genPx;
            axisPy = genPy;
            axisPt = std::hypot(axisPx, axisPy);
        }

        const int zPtBin = ptBin(axisPt);
        const int nJetBin = jetBin(b.njets);
        if (zPtBin < 0 || nJetBin < 0 || !(axisPt > 0.0)) continue;

        double metPt = 0.0;
        double metPhi = 0.0;
        const OfficialMETCorrections& metCorrection = isData ? dataMET : mcMET;
        if (!correctedMET(metCorrection, b, metPt, metPhi)) continue;

        const double metX = metPt * std::cos(metPhi);
        const double metY = metPt * std::sin(metPhi);
        const double recoilX = metX + visPx - genPx;
        const double recoilY = metY + visPy - genPy;
        const double cosPhi = axisPx / axisPt;
        const double sinPhi = axisPy / axisPt;
        const double uParallel =
            recoilX * cosPhi + recoilY * sinPhi;
        const double uPerpendicular =
            -recoilX * sinPhi + recoilY * cosPhi;

        double weight = isData ? 1.0 : nominalWeight(b, crossSectionWeight);
        if (APPLY_ZPT_REWEIGHTING && isDY) {
            weight *= zPtWeight(zPtWeights, axisPt);
        }

        if (!std::isfinite(weight)) continue;
        fillAccumulator(
            accumulator,
            group,
            nJetBin,
            zPtBin,
            uParallel,
            uPerpendicular,
            weight
        );
        ++selected;
    }

    std::cout << "[done] " << baseName << " selected=" << selected << std::endl;
    return 0;
}

std::unique_ptr<TH1D> makeBinnedHistogram(
    const std::string& name,
    const std::string& yTitle
) {
    std::unique_ptr<TH1D> histogram(new TH1D(
        name.c_str(),
        (";p_{T}^{Z} [GeV];" + yTitle).c_str(),
        static_cast<int>(ZPT_EDGES.size()) - 1,
        ZPT_EDGES.data()
    ));
    histogram->SetDirectory(nullptr);
    return histogram;
}

void drawSummary(
    TH1D* data,
    TH1D* dy,
    const std::string& outputName
) {
    TCanvas canvas("summary", "summary", 900, 750);
    canvas.SetLeftMargin(0.14);
    canvas.SetBottomMargin(0.12);
    canvas.SetGrid();

    data->SetLineColor(kBlack);
    data->SetMarkerColor(kBlack);
    data->SetMarkerStyle(20);
    dy->SetLineColor(kRed + 1);
    dy->SetMarkerColor(kRed + 1);
    dy->SetMarkerStyle(24);

    data->Draw("E1");
    dy->Draw("E1 SAME");

    TLegend legend(0.62, 0.76, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.AddEntry(data, "Data - non-DY", "lep");
    legend.AddEntry(dy, "DY simulation", "lep");
    legend.Draw();
    canvas.SaveAs(outputName.c_str());
}

Estimate usableEstimate(const Estimate& preferred, const Estimate& fallback) {
    if (preferred.valid && preferred.sigma > 0.0) return preferred;
    if (fallback.valid && fallback.sigma > 0.0) return fallback;

    Estimate result;
    result.valid = true;
    result.mean = 0.0;
    result.sigma = 30.0;
    return result;
}

void writeGaussian(
    TFile& output,
    const std::string& name,
    double mean,
    double sigma
) {
    sigma = std::max(0.5, sigma);
    const double xmin = mean - 8.0 * sigma;
    const double xmax = mean + 8.0 * sigma;
    TF1 function(name.c_str(), "gaus", xmin, xmax);
    function.SetParameters(
        1.0 / (std::sqrt(2.0 * TMath::Pi()) * sigma),
        mean,
        sigma
    );
    output.cd();
    function.Write(name.c_str(), TObject::kOverwrite);
}

void writeCompatiblePayload(
    const std::string& year,
    const std::string& outputDirectory,
    Accumulator& accumulator
) {
    const std::string fileName =
        "TypeI-PFMet_" + year + ".root";
    const std::string outputName = outputDirectory + "/" + fileName;
    TFile output(outputName.c_str(), "RECREATE");

    TH1D projection("projH", "", 2, -1.0, 1.0);
    projection.GetXaxis()->SetBinLabel(1, "recoilZPerp");
    projection.GetXaxis()->SetBinLabel(2, "recoilZParal");
    projection.Write();

    TH1D zPtBins(
        "ZPtBinsH",
        "",
        static_cast<int>(ZPT_EDGES.size()) - 1,
        ZPT_EDGES.data()
    );
    const char* ptLabels[] = {
        "Pt0to10", "Pt10to20", "Pt20to30", "Pt30to50", "PtGt50"
    };
    for (int pt = 0; pt < 5; ++pt)
        zPtBins.GetXaxis()->SetBinLabel(pt + 1, ptLabels[pt]);
    zPtBins.Write();

    TH1D jetBins("nJetBinsH", "", 3, -0.5, 2.5);
    const char* jetLabels[] = {"NJet0", "NJet1", "NJetGe2"};
    for (int jet = 0; jet < 3; ++jet)
        jetBins.GetXaxis()->SetBinLabel(jet + 1, jetLabels[jet]);
    jetBins.Write();

    const char* componentLabels[] = {"recoilZParal", "recoilZPerp"};
    const int nPtBins = static_cast<int>(ZPT_EDGES.size()) - 1;

    for (int jet = 0; jet < 3; ++jet) {
        for (int pt = 0; pt < nPtBins; ++pt) {
            const int index = flatBin(jet, pt);

            for (int component = 0; component < N_COMPONENTS; ++component) {
                const Moment targetMoment = subtractBackground(
                    accumulator.moments[DATA][component][index],
                    accumulator.moments[BACKGROUND][component][index]
                );
                const Estimate rawTarget = estimate(targetMoment);
                const Estimate rawMC = estimate(
                    accumulator.moments[DY][component][index]
                );
                const Estimate target = usableEstimate(rawTarget, rawMC);
                const Estimate simulation = usableEstimate(rawMC, rawTarget);

                const double targetMean =
                    component == PERPENDICULAR ? 0.0 : target.mean;
                const double simulationMean =
                    component == PERPENDICULAR ? 0.0 : simulation.mean;
                const std::string baseName =
                    std::string(componentLabels[component]) + "_" +
                    jetLabels[jet] + ptLabels[pt];

                writeGaussian(
                    output,
                    baseName + "_data",
                    targetMean,
                    target.sigma
                );
                writeGaussian(
                    output,
                    baseName + "_mc",
                    simulationMean,
                    simulation.sigma
                );

                TH1D* dataDistribution =
                    accumulator.distributions[DATA][component][index];
                TH1D* backgroundDistribution =
                    accumulator.distributions[BACKGROUND][component][index];
                TH1D* mcDistribution =
                    accumulator.distributions[DY][component][index];

                if (dataDistribution) {
                    std::unique_ptr<TH1D> targetHistogram(
                        dynamic_cast<TH1D*>(dataDistribution->Clone(
                            (baseName + "_hist_data").c_str()
                        ))
                    );
                    targetHistogram->SetDirectory(nullptr);
                    if (backgroundDistribution)
                        targetHistogram->Add(backgroundDistribution, -1.0);
                    output.cd();
                    targetHistogram->Write();
                }

                if (mcDistribution) {
                    std::unique_ptr<TH1D> mcHistogram(
                        dynamic_cast<TH1D*>(mcDistribution->Clone(
                            (baseName + "_hist_mc").c_str()
                        ))
                    );
                    mcHistogram->SetDirectory(nullptr);
                    output.cd();
                    mcHistogram->Write();
                }
            }
        }
    }

    output.Close();
    std::cout << "Wrote RecoilCorrector-compatible payload to "
              << outputName << std::endl;

    const char* cmsswBase = std::getenv("CMSSW_BASE");
    if (cmsswBase) {
        const std::string installedName =
            std::string(cmsswBase) +
            "/src/HTT-utilities/RecoilCorrections/data/" + fileName;
        if (gSystem->CopyFile(
            outputName.c_str(),
            installedName.c_str(),
            true
        ) == 0) {
            std::cout << "Installed recoil payload at "
                      << installedName << std::endl;
        }
        else {
            std::cerr << "WARNING: could not copy recoil payload to "
                      << installedName << std::endl;
        }
    }
}

void writePayload(
    const std::string& year,
    const std::string& outputDirectory,
    Accumulator& accumulator
) {
    const std::string outputName =
        outputDirectory + "/RecoilCorrections_" + year + ".root";
    TFile outputFile(outputName.c_str(), "RECREATE");

    const int nPtBins = static_cast<int>(ZPT_EDGES.size()) - 1;
    const char* componentName[] = {"upar", "uperp"};

    for (int jet = 0; jet < 3; ++jet) {
        for (int component = 0; component < N_COMPONENTS; ++component) {
            auto meanData = makeBinnedHistogram(
                Form("h_%s_mean_data_njet%d", componentName[component], jet),
                Form("#mu(%s) [GeV]", componentName[component])
            );
            auto meanDY = makeBinnedHistogram(
                Form("h_%s_mean_mc_njet%d", componentName[component], jet),
                Form("#mu(%s) [GeV]", componentName[component])
            );
            auto sigmaData = makeBinnedHistogram(
                Form("h_%s_sigma_data_njet%d", componentName[component], jet),
                Form("#sigma(%s) [GeV]", componentName[component])
            );
            auto sigmaDY = makeBinnedHistogram(
                Form("h_%s_sigma_mc_njet%d", componentName[component], jet),
                Form("#sigma(%s) [GeV]", componentName[component])
            );
            auto shift = makeBinnedHistogram(
                Form("h_%s_shift_njet%d", componentName[component], jet),
                "Data-MC mean shift [GeV]"
            );
            auto scale = makeBinnedHistogram(
                Form("h_%s_scale_njet%d", componentName[component], jet),
                "Data/MC resolution scale"
            );
            auto offset = makeBinnedHistogram(
                Form("h_%s_offset_njet%d", componentName[component], jet),
                "Affine offset [GeV]"
            );
            auto valid = makeBinnedHistogram(
                Form("h_%s_valid_njet%d", componentName[component], jet),
                "Valid correction bin"
            );

            for (int pt = 0; pt < nPtBins; ++pt) {
                const int index = flatBin(jet, pt);
                const Moment targetMoment = subtractBackground(
                    accumulator.moments[DATA][component][index],
                    accumulator.moments[BACKGROUND][component][index]
                );
                const Estimate target = estimate(targetMoment);
                const Estimate simulation = estimate(
                    accumulator.moments[DY][component][index]
                );
                const int bin = pt + 1;

                if (target.valid) {
                    meanData->SetBinContent(bin, target.mean);
                    meanData->SetBinError(bin, target.meanError);
                    sigmaData->SetBinContent(bin, target.sigma);
                    sigmaData->SetBinError(bin, target.sigmaError);
                }
                if (simulation.valid) {
                    meanDY->SetBinContent(bin, simulation.mean);
                    meanDY->SetBinError(bin, simulation.meanError);
                    sigmaDY->SetBinContent(bin, simulation.sigma);
                    sigmaDY->SetBinError(bin, simulation.sigmaError);
                }

                const bool correctionValid =
                    target.valid && simulation.valid && simulation.sigma > 0.0;
                if (correctionValid) {
                    const double shiftValue = target.mean - simulation.mean;
                    const double scaleValue = target.sigma / simulation.sigma;
                    const double shiftError = std::hypot(
                        target.meanError,
                        simulation.meanError
                    );
                    const double relativeTarget = target.sigmaError / target.sigma;
                    const double relativeSimulation =
                        simulation.sigmaError / simulation.sigma;

                    shift->SetBinContent(bin, shiftValue);
                    shift->SetBinError(bin, shiftError);
                    scale->SetBinContent(bin, scaleValue);
                    scale->SetBinError(
                        bin,
                        scaleValue * std::hypot(relativeTarget, relativeSimulation)
                    );
                    offset->SetBinContent(
                        bin,
                        target.mean - scaleValue * simulation.mean
                    );
                    valid->SetBinContent(bin, 1.0);
                }
                else {
                    shift->SetBinContent(bin, 0.0);
                    scale->SetBinContent(bin, 1.0);
                    offset->SetBinContent(bin, 0.0);
                    valid->SetBinContent(bin, 0.0);
                }

                TH1D* dataDistribution =
                    accumulator.distributions[DATA][component][index];
                TH1D* backgroundDistribution =
                    accumulator.distributions[BACKGROUND][component][index];
                TH1D* dyDistribution =
                    accumulator.distributions[DY][component][index];

                if (dataDistribution) dataDistribution->Write();
                if (backgroundDistribution) backgroundDistribution->Write();
                if (dyDistribution) dyDistribution->Write();

                if (dataDistribution) {
                    std::unique_ptr<TH1D> targetDistribution(
                        dynamic_cast<TH1D*>(dataDistribution->Clone(
                            Form("h_target_%s_njet%d_pt%d",
                                 componentName[component], jet, pt)
                        ))
                    );
                    targetDistribution->SetDirectory(nullptr);
                    if (backgroundDistribution) targetDistribution->Add(backgroundDistribution, -1.0);
                    targetDistribution->Write();
                }
            }

            meanData->Write();
            meanDY->Write();
            sigmaData->Write();
            sigmaDY->Write();
            shift->Write();
            scale->Write();
            offset->Write();
            valid->Write();

            drawSummary(
                meanData.get(),
                meanDY.get(),
                outputDirectory + "/" + componentName[component] +
                    "_mean_njet" + std::to_string(jet) + ".png"
            );
            drawSummary(
                sigmaData.get(),
                sigmaDY.get(),
                outputDirectory + "/" + componentName[component] +
                    "_sigma_njet" + std::to_string(jet) + ".png"
            );
        }
    }

    outputFile.Close();
    std::cout << "Wrote recoil payload to " << outputName << std::endl;
}

std::vector<std::pair<std::string,std::string>> buildFileList(const std::string& year) {
    std::vector<std::pair<std::string,std::string>> files;
    const std::map<std::string, std::vector<std::string>> fileMap = getFileMap(year);
    for (const auto& processEntry : fileMap) {
        for (const std::string& fileName : processEntry.second) {
            if (fileName.empty() || fileName.back() == '/') continue;
            files.push_back({processEntry.first, fileName});
        }
    }
    return files;
}

void writeAccumulator(Accumulator& accumulator, TFile& output) {
    output.cd();
    const int nPtBins = static_cast<int>(ZPT_EDGES.size()) - 1;

    int group = 0;
    int component = 0;
    int index = 0;
    double sumW = 0.0;
    double sumWX = 0.0;
    double sumWX2 = 0.0;
    double sumW2 = 0.0;

    TTree momentsTree("moments", "moments");
    momentsTree.Branch("group", &group);
    momentsTree.Branch("component", &component);
    momentsTree.Branch("index", &index);
    momentsTree.Branch("sumW", &sumW);
    momentsTree.Branch("sumWX", &sumWX);
    momentsTree.Branch("sumWX2", &sumWX2);
    momentsTree.Branch("sumW2", &sumW2);

    for (group = 0; group < N_GROUPS; ++group) {
        for (component = 0; component < N_COMPONENTS; ++component) {
            for (index = 0; index < 3 * nPtBins; ++index) {
                const Moment& moment = accumulator.moments[group][component][index];
                sumW = moment.sumW;
                sumWX = moment.sumWX;
                sumWX2 = moment.sumWX2;
                sumW2 = moment.sumW2;
                momentsTree.Fill();

                TH1D* histogram = accumulator.distributions[group][component][index];
                if (histogram) histogram->Write();
            }
        }
    }
    momentsTree.Write();
}

void readAndAddAccumulator(Accumulator& accumulator, const std::string& fileName) {
    std::unique_ptr<TFile> input(TFile::Open(fileName.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "[skip] cannot open perfile accumulator: " << fileName << std::endl;
        return;
    }

    TTree* momentsTree = dynamic_cast<TTree*>(input->Get("moments"));
    if (!momentsTree) {
        std::cerr << "[skip] no moments tree in: " << fileName << std::endl;
        return;
    }

    int group = 0;
    int component = 0;
    int index = 0;
    double sumW = 0.0;
    double sumWX = 0.0;
    double sumWX2 = 0.0;
    double sumW2 = 0.0;
    momentsTree->SetBranchAddress("group", &group);
    momentsTree->SetBranchAddress("component", &component);
    momentsTree->SetBranchAddress("index", &index);
    momentsTree->SetBranchAddress("sumW", &sumW);
    momentsTree->SetBranchAddress("sumWX", &sumWX);
    momentsTree->SetBranchAddress("sumWX2", &sumWX2);
    momentsTree->SetBranchAddress("sumW2", &sumW2);

    const Long64_t nEntries = momentsTree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        momentsTree->GetEntry(entry);
        Moment& target = accumulator.moments[group][component][index];
        target.sumW += sumW;
        target.sumWX += sumWX;
        target.sumWX2 += sumWX2;
        target.sumW2 += sumW2;
    }

    const char* groupName[] = {"data", "background", "dy"};
    const char* componentName[] = {"upar", "uperp"};
    const int nPtBins = static_cast<int>(ZPT_EDGES.size()) - 1;

    for (int g = 0; g < N_GROUPS; ++g) {
        for (int c = 0; c < N_COMPONENTS; ++c) {
            for (int jet = 0; jet < 3; ++jet) {
                for (int pt = 0; pt < nPtBins; ++pt) {
                    const int idx = flatBin(jet, pt);
                    const std::string name = Form("h_%s_%s_njet%d_pt%d", groupName[g], componentName[c], jet, pt);
                    TH1D* source = dynamic_cast<TH1D*>(input->Get(name.c_str()));
                    if (!source) continue;
                    TH1D*& target = accumulator.distributions[g][c][idx];
                    if (!target) {
                        target = dynamic_cast<TH1D*>(source->Clone(name.c_str()));
                        target->SetDirectory(nullptr);
                    }
                    else {
                        target->Add(source);
                    }
                }
            }
        }
    }
}

}

void DeriveRecoilCorrections(std::string year = "2018") {
    if (year == "Run2") {
        const std::vector<std::string> years = {
            "2016preVFP", "2016postVFP", "2017", "2018"
        };
        for (const std::string& runYear : years)
            DeriveRecoilCorrections(runYear);
        return;
    }

    const double lumi = luminosity(year);
    if (!(lumi > 0.0)) {
        std::cerr << "ERROR: unsupported year " << year << std::endl;
        return;
    }

    const std::map<std::string, std::vector<std::string>> fileMap =
        getFileMap(year);
    if (fileMap.empty()) {
        std::cerr << "ERROR: empty FileMap for " << year << std::endl;
        return;
    }

    std::unique_ptr<OfficialMETCorrections> dataMET =
        loadMETCorrections(year, true);
    std::unique_ptr<OfficialMETCorrections> mcMET =
        loadMETCorrections(year, false);
    if (!dataMET || !mcMET) return;

    std::unique_ptr<TH1D> zPtWeights;
    if (APPLY_ZPT_REWEIGHTING) {
        zPtWeights = loadZPtWeights(year);
        if (!zPtWeights) {
            std::cerr << "ERROR: cannot load Z pT weights for " << year << std::endl;
            return;
        }
    }

    const std::string outputDirectory = OUTPUT_BASE + "/" + year;
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);

    Accumulator accumulator;
    int failed = 0;
    int succeeded = 0;

    for (const auto& processEntry : fileMap) {
        for (const std::string& fileName : processEntry.second) {
            if (fileName.empty() || fileName.back() == '/') continue;
            const int status = processFile(
                year,
                processEntry.first,
                fileName,
                lumi,
                accumulator,
                *dataMET,
                *mcMET,
                zPtWeights.get()
            );
            if (status == 0) ++succeeded;
            else ++failed;
        }
    }

    if (succeeded == 0) {
        std::cerr << "ERROR: no input files were processed successfully for "
                  << year << "; no recoil payload was written" << std::endl;
        gSystem->Exit(2);
        return;
    }

    writePayload(year, outputDirectory, accumulator);
    writeCompatiblePayload(year, outputDirectory, accumulator);
    std::cout << "Finished " << year << " with " << succeeded
              << " successful and " << failed << " failed input files"
              << std::endl;
}

void DeriveRecoilCorrectionsPerFile(std::string year, int fileIndex) {
    const double lumi = luminosity(year);
    if (!(lumi > 0.0)) {
        std::cerr << "ERROR: unsupported year " << year << std::endl;
        return;
    }

    const std::vector<std::pair<std::string,std::string>> files = buildFileList(year);
    if (fileIndex < 0 || fileIndex >= static_cast<int>(files.size())) {
        std::cerr << "ERROR: invalid fileIndex=" << fileIndex << ", total files=" << files.size() << std::endl;
        return;
    }

    std::unique_ptr<OfficialMETCorrections> dataMET = loadMETCorrections(year, true);
    std::unique_ptr<OfficialMETCorrections> mcMET = loadMETCorrections(year, false);
    if (!dataMET || !mcMET) return;

    std::unique_ptr<TH1D> zPtWeights;
    if (APPLY_ZPT_REWEIGHTING) {
        zPtWeights = loadZPtWeights(year);
        if (!zPtWeights) {
            std::cerr << "ERROR: cannot load Z pT weights for " << year << std::endl;
            return;
        }
    }

    const std::string outputDirectory = OUTPUT_BASE + "/" + year + "/perfile";
    gSystem->mkdir(outputDirectory.c_str(), kTRUE);

    Accumulator accumulator;
    const int status = processFile(
        year,
        files[fileIndex].first,
        files[fileIndex].second,
        lumi,
        accumulator,
        *dataMET,
        *mcMET,
        zPtWeights.get()
    );

    const std::string baseName = gSystem->BaseName(files[fileIndex].second.c_str());
    const std::string outputName = outputDirectory + "/" + std::to_string(fileIndex) + "_" + baseName;
    TFile output(outputName.c_str(), "RECREATE");
    writeAccumulator(accumulator, output);
    output.Close();

    std::cout << "Wrote perfile accumulator (status=" << status << ") to " << outputName << std::endl;
    if (status != 0) gSystem->Exit(status);
}

void DeriveRecoilCorrectionsMerge(std::string year) {
    if (year == "Run2") {
        const std::vector<std::string> years = {
            "2016preVFP", "2016postVFP", "2017", "2018"
        };
        for (const std::string& runYear : years)
            DeriveRecoilCorrectionsMerge(runYear);
        return;
    }

    const std::vector<std::pair<std::string,std::string>> files = buildFileList(year);
    const std::string perfileDirectory = OUTPUT_BASE + "/" + year + "/perfile";

    Accumulator accumulator;
    int found = 0;
    int missing = 0;
    for (size_t fileIndex = 0; fileIndex < files.size(); ++fileIndex) {
        const std::string baseName = gSystem->BaseName(files[fileIndex].second.c_str());
        const std::string perfileName = perfileDirectory + "/" + std::to_string(fileIndex) + "_" + baseName;
        if (gSystem->AccessPathName(perfileName.c_str())) {
            std::cerr << "[missing] " << perfileName << std::endl;
            ++missing;
            continue;
        }
        readAndAddAccumulator(accumulator, perfileName);
        ++found;
    }

    if (found == 0) {
        std::cerr << "ERROR: no perfile accumulators found for " << year << "; no recoil payload was written" << std::endl;
        gSystem->Exit(2);
        return;
    }

    const std::string outputDirectory = OUTPUT_BASE + "/" + year;
    writePayload(year, outputDirectory, accumulator);
    writeCompatiblePayload(year, outputDirectory, accumulator);
    std::cout << "Merged " << year << " with " << found
              << " found and " << missing << " missing perfile accumulators"
              << std::endl;
}
