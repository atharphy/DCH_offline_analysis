#include "../../filemap/FileMap.h"
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <iostream>

void quick_check_numden(std::string year = "2018") {
    ROOT::EnableImplicitMT();

    TChain ch("Events");
    auto files = getFileMap(year);
    for (const std::string& group : {"DY", "DY10_50"})
        for (const auto& f : files[group]) ch.Add(f.c_str());
    std::cout << "[" << year << "] Files in chain: " << ch.GetListOfFiles()->GetEntries() << std::endl;

    ROOT::RDataFrame frame(ch);

    auto numFrame = frame.Filter("cat==42 && genPartFlav_1==1 && genPartFlav_2==1");
    auto h_numBarrel = numFrame.Filter("abs(eta_2)<1.5").Histo1D({("h_numBarrel_" + year).c_str(), "", 200, 0, 200}, "pt_2");
    auto h_numEndcap = numFrame.Filter("abs(eta_2)>=1.5 && abs(eta_2)<2.5").Histo1D({("h_numEndcap_" + year).c_str(), "", 200, 0, 200}, "pt_2");

    auto looseFrame = frame.Filter("cat==46 && genPartFlav_1==1")
                            .Define("looseMask", "lflavor==15 && gen_match==1 && TauIDj>2")
                            .Define("loosePtBarrel", "lpt[looseMask && abs(leta)<1.5]")
                            .Define("loosePtEndcap", "lpt[looseMask && abs(leta)>=1.5 && abs(leta)<2.5]");
    auto h_looseBarrel = looseFrame.Histo1D({("h_looseBarrel_" + year).c_str(), "", 200, 0, 200}, "loosePtBarrel");
    auto h_looseEndcap = looseFrame.Histo1D({("h_looseEndcap_" + year).c_str(), "", 200, 0, 200}, "loosePtEndcap");

    h_numBarrel->GetEntries();

    TH1D* hNumBarrel = (TH1D*)h_numBarrel->Clone(("Num_" + year + "_barrel").c_str());
    TH1D* hNumEndcap = (TH1D*)h_numEndcap->Clone(("Num_" + year + "_endcap").c_str());
    TH1D* hLooseBarrel = (TH1D*)h_looseBarrel->Clone(("Loose_" + year + "_barrel").c_str());
    TH1D* hLooseEndcap = (TH1D*)h_looseEndcap->Clone(("Loose_" + year + "_endcap").c_str());

    TH1D* hDenBarrel = (TH1D*)hNumBarrel->Clone(("Den_" + year + "_barrel").c_str());
    hDenBarrel->Add(hLooseBarrel);
    TH1D* hDenEndcap = (TH1D*)hNumEndcap->Clone(("Den_" + year + "_endcap").c_str());
    hDenEndcap->Add(hLooseEndcap);

    TH1D* hRatioBarrel = (TH1D*)hNumBarrel->Clone(("Ratio_" + year + "_barrel").c_str());
    hRatioBarrel->Divide(hNumBarrel, hDenBarrel, 1.0, 1.0, "B");
    TH1D* hRatioEndcap = (TH1D*)hNumEndcap->Clone(("Ratio_" + year + "_endcap").c_str());
    hRatioEndcap->Divide(hNumEndcap, hDenEndcap, 1.0, 1.0, "B");

    std::cout << "[" << year << "] barrel: num=" << hNumBarrel->GetEntries()
              << " loose=" << hLooseBarrel->GetEntries() << " den=" << hDenBarrel->GetEntries() << std::endl;
    std::cout << "[" << year << "] endcap: num=" << hNumEndcap->GetEntries()
              << " loose=" << hLooseEndcap->GetEntries() << " den=" << hDenEndcap->GetEntries() << std::endl;

    gSystem->mkdir("cut_results", kTRUE);
    const std::string outName = "cut_results/quick_check_numden_" + year + ".root";
    TFile fout(outName.c_str(), "RECREATE");
    hNumBarrel->Write();
    hNumEndcap->Write();
    hLooseBarrel->Write();
    hLooseEndcap->Write();
    hDenBarrel->Write();
    hDenEndcap->Write();
    hRatioBarrel->Write();
    hRatioEndcap->Write();
    fout.Close();

    std::cout << "[" << year << "] output: " << outName << std::endl;
}
