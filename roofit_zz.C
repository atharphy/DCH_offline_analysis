#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooDataHist.h>
#include <RooConstVar.h>
#include <RooHistPdf.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooExtendPdf.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <fstream>
#include "include/Xsections.C"
#include "filemap/FileMap.h"
static std::string SIGNAL_MASS_FILTER = "";
#include "Stack_modules/FileDiscovery.h"
#include "roofit_formatting/RooFitStyle.h"

using namespace RooFit;
static const std::string ZZ_SF_CSV = "Dependencies/wz_zz_scale_factors/zz_scale_factors.csv";

void roofit_zz(std::string year = "Run2"){
	std::string region = "CR_0tau";
	std::string inputDir = "hists/run2_hists_noFR_roccor/";

	std::map<std::string, std::vector<TFile*>> open_files;
	openInputFiles(year, ensureTrailingSlash(inputDir), "hist_", open_files);

	std::map<std::string, TH1D*> h_bkg_group;
	std::map<std::string, double> tot_uncert_quadr;
	for (auto& kv : open_files) {
		for (auto* f : kv.second) {
			TH1D* h = nullptr;
			for (const std::string& ch : {"eeee","mmmm","eemm"}) {
				TH1D* part = (TH1D*)f->Get(("h_LT_" + ch + "_" + region).c_str());
				if (!part) continue;
				if (!h) { h = (TH1D*)part->Clone(); h->SetDirectory(nullptr); }
				else h->Add(part);
			}
			if(!h) continue;
			h->Rebin(10);
			h->Sumw2();
			if (!h_bkg_group[kv.first]){
				h_bkg_group[kv.first]= dynamic_cast<TH1D*>(h);
            	tot_uncert_quadr[kv.first] = fake_uncert_squared("4lep", f->GetName())*h->Integral()*h->Integral();
    		}
	       	else {
	       		h_bkg_group[kv.first]->Add(h);
	       		tot_uncert_quadr[kv.first] += fake_uncert_squared("4lep", f->GetName())*h->Integral()*h->Integral();
	       	}
	       	cout<<kv.first<<"\t"<<tot_uncert_quadr[kv.first]<<"\t"<<h_bkg_group[kv.first]->Integral()<<endl;
		}
	}
	TH1D* h_other_bkg = (TH1D*)h_bkg_group["other"]->Clone();
	h_other_bkg->Add(h_bkg_group["DY"]);
	h_other_bkg->Add(h_bkg_group["DY10_50"]);
	h_other_bkg->Add(h_bkg_group["WZ"]);
	h_other_bkg->Add(h_bkg_group["WW"]);
	h_other_bkg->Add(h_bkg_group["VVV"]);
	h_other_bkg->Add(h_bkg_group["ttV"]);
	h_other_bkg->Add(h_bkg_group["WJ"]);
	h_other_bkg->Add(h_bkg_group["ST"]);
	h_other_bkg->Add(h_bkg_group["TTbar"]);
	h_other_bkg->Add(h_bkg_group["QCD"]);
	double other_nominal = h_other_bkg->Integral();
	double other_uncert = sqrt(tot_uncert_quadr["other"]+tot_uncert_quadr["DY"]+tot_uncert_quadr["DY10_50"]+tot_uncert_quadr["WZ"]+tot_uncert_quadr["WW"]+tot_uncert_quadr["VVV"]+tot_uncert_quadr["ttV"]+tot_uncert_quadr["WJ"]+tot_uncert_quadr["ST"]+tot_uncert_quadr["TTbar"]+tot_uncert_quadr["QCD"])/h_other_bkg->Integral();

	RooRealVar x("x", "L_{T} variable", 0, 1000);
	RooDataHist other_hist("other_hist", "Other", x, Import(*h_other_bkg));
	RooHistPdf other_pdf("other_pdf", "Other PDF", x, other_hist);

	RooRealVar other_nuis("other_nuis", "Other nuisance", 0, -5, 5);
	RooGaussian other_constraint("other_constraint", "Other constraint", other_nuis, RooConst(0.), RooConst(1.));
	RooFormulaVar other_norm_constrained("other_norm_constrained", "@0*(1 + @1*@2)", RooArgList(RooConst(other_nominal), other_nuis, other_uncert));
	RooRealVar other_norm("other_norm", "Other yield", h_other_bkg->Integral(), 0.0, 10.0 * h_other_bkg->Integral());
	RooExtendPdf other_ext("other_ext", "Other Extended PDF", other_pdf, other_norm_constrained);

	RooDataHist zz_hist("zz_hist", "ZZ", x, Import(*h_bkg_group["ZZ"]));
	RooHistPdf zz_pdf("zz_pdf", "ZZ PDF", x, zz_hist);
	RooRealVar zz_norm("zz_norm", "ZZ yield", h_bkg_group["ZZ"]->Integral(), 0.0, 10.0 * h_bkg_group["ZZ"]->Integral());
	RooExtendPdf zz_ext("zz_ext", "ZZ Extended", zz_pdf, zz_norm);

	RooAddPdf model_core("model_core", "Total Model without constraints", RooArgList(zz_ext, other_ext));
	RooProdPdf model("model", "Model with constraints", RooArgSet(model_core, other_constraint));

	RooDataHist data_obs("data_obs", "Observed Data", x, Import(*h_bkg_group["data"]));
	model.fitTo(data_obs, Extended(true), PrintLevel(-1));
	TCanvas c(("c_roofit_zz_" + year).c_str(), "", 1000, 800);
	RooPlot* frame = x.frame();

	data_obs.plotOn(frame, Name("stack_data"));
	model_core.plotOn(frame, LineColor(kBlue), Name("stack_total"));
	model_core.plotOn(frame, Components(zz_ext), LineColor(kGreen + 1), Name("stack_zz"));
	model_core.plotOn(frame, Components(other_ext), LineColor(kRed), Name("stack_other"));

	double zz_fitted = zz_norm.getVal();
	double zz_nominal = h_bkg_group["ZZ"]->Integral();
	double zz_fit_err = zz_norm.getError();
	double scale_factor_ZZ = zz_fitted / zz_nominal;
	double sf_err = zz_fit_err / zz_nominal;
	std::cout << "ZZ scale factor = " << scale_factor_ZZ <<"+-"<< sf_err<< std::endl;

	char sfBuf[64];
	snprintf(sfBuf, sizeof(sfBuf), "ZZ SF = %.3f #pm %.3f", scale_factor_ZZ, sf_err);
	std::vector<RooFitLegendEntry> legendEntries = {
		{frame->findObject("stack_zz"), "ZZ", "l"},
		{frame->findObject("stack_other"), "Other", "l"},
		{frame->findObject("stack_data"), "Data", "lep"}
	};
	formatRooFitCanvas(c, frame, "L_{T} [GeV]", year, legendEntries, sfBuf);

	gSystem->mkdir("normfits/plots", kTRUE);
	c.SaveAs(("normfits/plots/roofit_zz_" + region + "_" + year + ".png").c_str());

	gSystem->mkdir("normfits", kTRUE);
	std::map<std::string, std::string> csvRows;
	std::ifstream fin(ZZ_SF_CSV);
	std::string line;
	if (fin) { std::getline(fin, line); while (std::getline(fin, line)) { size_t p = line.find(','); if (p != std::string::npos) csvRows[line.substr(0,p)] = line.substr(p+1); } }
	fin.close();
	csvRows[year] = std::to_string(scale_factor_ZZ) + "," + std::to_string(sf_err);
	std::ofstream fout(ZZ_SF_CSV);
	fout << "year,scale_factor,error\n";
	for (auto& kv : csvRows) fout << kv.first << "," << kv.second << "\n";
	fout.close();
	gSystem->CopyFile(ZZ_SF_CSV.c_str(), "zz_scale_factors.csv", kTRUE);
}
