#ifndef FILEMAP_H
#define FILEMAP_H

#include <string>
#include <vector>
#include <map>

const std::string MCBASE   = "root://cmseos.fnal.gov//store/user/aahmad2/run2_files/run2_skims";
const std::string DATABASE = "root://cmseos.fnal.gov//store/user/aahmad2/run2_files/run2_skims_data_updated_lumi";

const std::string SIGNALBASE = "root://cmseos.fnal.gov//store/user/aahmad2/run2_files/run2_signal";

inline std::string mcPath(const std::string& year, const std::string& f) {
    return MCBASE + "/MC_" + year + "/" + f;
}

inline std::string dataPath(const std::string& year, const std::string& f) {
    return DATABASE + "/Data_" + year + "/" + f;
}

inline std::string signalPath(const std::string& year, const std::string& f) {
    return SIGNALBASE + "/MC_" + year + "/" + f;
}

inline std::vector<std::string> mcList(const std::string& year,
                                       const std::vector<std::string>& files) {
    std::vector<std::string> out;
    for (const auto& f : files) out.push_back(mcPath(year, f));
    return out;
}

inline std::vector<std::string> dataList(const std::string& year,
                                         const std::vector<std::string>& files) {
    std::vector<std::string> out;
    for (const auto& f : files) out.push_back(dataPath(year, f));
    return out;
}

inline std::vector<std::string> signalList(const std::string& year,
                                           const std::vector<std::string>& files) {
    std::vector<std::string> out;
    for (const auto& f : files) out.push_back(signalPath(year, f));
    return out;
}

std::map<std::string, std::vector<std::string>> getFileMap(std::string year) {
    std::map<std::string, std::vector<std::string>> files;

    if (year == "2016preVFP") {
        files = {
            {"DY10_50", mcList(year, {"DYJetsToLLM10to50_2016preVFP.root"})},
            {"DY",      mcList(year, {"DYJetsToLLM50_2016preVFP.root"})},
            {"WZ",      mcList(year, {"WZTo2Q2L_2016preVFP.root", "WZTo3LNu_2016preVFP.root"})},
            {"WW",      mcList(year, {"WW_2016preVFP.root", "WWTo2L2Nu_2016preVFP.root"})},
            {"VVV",     mcList(year, {"WWW_2016preVFP.root", "WZZ_2016preVFP.root", "ZZZ_2016preVFP.root"})},
            {"ttV",     mcList(year, {"ttWJets_2016preVFP.root", "ttZJets_2016preVFP.root"})},
            {"WJ",      mcList(year, {"WJetsToLNu_NLO_2016preVFP.root"})},
            {"ZZ",      mcList(year, {"ZZTo2L2Nu_2016preVFP.root", "ZZTo2Q2L_2016preVFP.root", "ZZTo4L_2016preVFP.root"})},
            {"ST",      mcList(year, {"ST_s-channel_2016preVFP.root", "ST_t-channel_antitop_2016preVFP.root", "ST_t-channel_top_2016preVFP.root", "ST_tW_antitop_2016preVFP.root", "ST_tW_top_2016preVFP.root"})},
            {"TTbar",   mcList(year, {"TTTo2L2Nu_2016preVFP.root", "TTToSemiLeptonic_2016preVFP.root", "TTToHadronic_2016preVFP.root"})},
            {"QCD",     mcList(year, {"QCD_HT50to100_2016preVFP.root", "QCD_HT100to200_2016preVFP.root", "QCD_HT200to300_2016preVFP.root", "QCD_HT300to500_2016preVFP.root", "QCD_HT500to700_2016preVFP.root", "QCD_HT700to1000_2016preVFP.root", "QCD_HT1000to1500_2016preVFP.root", "QCD_HT1500to2000_2016preVFP.root", "QCD_HT2000toInf_2016preVFP.root"})},
            {"other",   mcList(year, {"ttHToTauTau_2016preVFP.root", "ZHToMuMu_2016preVFP.root", "ZHToTauTau_2016preVFP.root", "GluGluZH_2016preVFP.root", "ttHToEE_2016preVFP.root", "ttHTo2L2Nu_2016preVFP.root", "ttHJetToNonbb_2016preVFP.root", "TWZToLL_2016preVFP.root", "HZJ_HToWWTo2L2Nu_ZTo2L_2016preVFP.root"})},
            {"signal",  signalList(year, {"HppM500_2016preVFP.root", "HppM600_2016preVFP.root", "HppM700_2016preVFP.root", "HppM800_2016preVFP.root", "HppM900_2016preVFP.root", "HppM1000_2016preVFP.root", "HppM1100_2016preVFP.root", "HppM1200_2016preVFP.root", "HppM1300_2016preVFP.root", "HppM1400_2016preVFP.root", "HppM1500_2016preVFP.root"})},
            {"data",    dataList(year, {"SingleElectronB_2016preVFP.root", "SingleElectronC_2016preVFP.root", "SingleElectronD_2016preVFP.root", "SingleElectronE_2016preVFP.root", "SingleElectronF_2016preVFP.root", "SingleMuonB_2016preVFP.root", "SingleMuonC_2016preVFP.root", "SingleMuonD_2016preVFP.root", "SingleMuonE_2016preVFP.root", "SingleMuonF_2016preVFP.root"})}
        };
    }

    else if (year == "2016postVFP") {
        files = {
            {"DY10_50", mcList(year, {"DYJetsToLLM10to50_2016postVFP.root"})},
            {"DY",      mcList(year, {"DYJetsToLLM50_2016postVFP.root"})},
            {"WZ",      mcList(year, {"WZTo2Q2L_2016postVFP.root", "WZTo3LNu_2016postVFP.root"})},
            {"WW",      mcList(year, {"WW_2016postVFP.root", "WWTo2L2Nu_2016postVFP.root"})},
            {"VVV",     mcList(year, {"WWW_2016postVFP.root", "WZZ_2016postVFP.root", "ZZZ_2016postVFP.root"})},
            {"ttV",     mcList(year, {"ttWJets_2016postVFP.root", "ttZJets_2016postVFP.root"})},
            {"WJ",      mcList(year, {"WJetsToLNu_NLO_2016postVFP.root"})},
            {"ZZ",      mcList(year, {"ZZTo2L2Nu_2016postVFP.root", "ZZTo2Q2L_2016postVFP.root", "ZZTo4L_2016postVFP.root"})},
            {"ST",      mcList(year, {"ST_s-channel_2016postVFP.root", "ST_t-channel_antitop_2016postVFP.root", "ST_t-channel_top_2016postVFP.root", "ST_tW_antitop_2016postVFP.root", "ST_tW_top_2016postVFP.root"})},
            {"TTbar",   mcList(year, {"TTTo2L2Nu_2016postVFP.root", "TTToSemiLeptonic_2016postVFP.root", "TTToHadronic_2016postVFP.root"})},
            {"QCD",     mcList(year, {"QCD_HT50to100_2016postVFP.root", "QCD_HT100to200_2016postVFP.root", "QCD_HT200to300_2016postVFP.root", "QCD_HT300to500_2016postVFP.root", "QCD_HT500to700_2016postVFP.root", "QCD_HT700to1000_2016postVFP.root", "QCD_HT1000to1500_2016postVFP.root", "QCD_HT1500to2000_2016postVFP.root", "QCD_HT2000toInf_2016postVFP.root"})},
            {"other",   mcList(year, {"ttHToTauTau_2016postVFP.root", "ZHToMuMu_2016postVFP.root", "ZHToTauTau_ext_2016postVFP.root", "GluGluZH_2016postVFP.root", "ttHToEE_2016postVFP.root", "ttHTo2L2Nu_2016postVFP.root", "ttHJetToNonbb_2016postVFP.root", "TWZToLL_2016postVFP.root", "HZJ_HToWWTo2L2Nu_ZTo2L_2016postVFP.root"})},
            {"signal",  signalList(year, {"HppM500_2016postVFP.root", "HppM600_2016postVFP.root", "HppM700_2016postVFP.root", "HppM800_2016postVFP.root", "HppM900_2016postVFP.root", "HppM1000_2016postVFP.root", "HppM1100_2016postVFP.root", "HppM1200_2016postVFP.root", "HppM1300_2016postVFP.root", "HppM1400_2016postVFP.root", "HppM1500_2016postVFP.root"})},
            {"data",    dataList(year, {"SingleElectronF_2016postVFP.root", "SingleElectronG_2016postVFP.root", "SingleElectronH_2016postVFP.root", "SingleMuonF_2016postVFP.root", "SingleMuonG_2016postVFP.root", "SingleMuonH_2016postVFP.root"})}
        };
    }

    else if (year == "2017") {
        files = {
            {"DY10_50", mcList(year, {"DYJetsToLLM10to50_2017.root"})},
            {"DY",      mcList(year, {"DYJetsToLLM50_2017.root"})},
            {"WZ",      mcList(year, {"WZTo2Q2L_2017.root", "WZTo3LNu_2017.root"})},
            {"WW",      mcList(year, {"WW_2017.root", "WWTo2L2Nu_2017.root"})},
            {"VVV",     mcList(year, {"WWW_2017.root", "WZZ_2017.root", "ZZZ_2017.root"})},
            {"ttV",     mcList(year, {"ttWJets_2017.root", "ttZJets_2017.root"})},
            {"WJ",      mcList(year, {"WJetsToLNu_NLO_2017.root"})},
            {"ZZ",      mcList(year, {"ZZTo2L2Nu_2017.root", "ZZTo2Q2L_2017.root", "ZZTo4L_2017.root"})},
            {"ST",      mcList(year, {"ST_s-channel_2017.root", "ST_t-channel_antitop_2017.root", "ST_t-channel_top_2017.root", "ST_tW_antitop_2017.root", "ST_tW_top_2017.root"})},
            {"TTbar",   mcList(year, {"TTTo2L2Nu_2017.root", "TTToSemiLeptonic_2017.root", "TTToHadronic_2017.root"})},
            {"QCD",     mcList(year, {"QCD_HT50to100_2017.root", "QCD_HT100to200_2017.root", "QCD_HT200to300_2017.root", "QCD_HT300to500_2017.root", "QCD_HT500to700_2017.root", "QCD_HT700to1000_2017.root", "QCD_HT1000to1500_2017.root", "QCD_HT1500to2000_2017.root", "QCD_HT2000toInf_2017.root"})},
            {"other",   mcList(year, {"ttHToTauTau_2017.root", "ZHToMuMu_2017.root", "ZHToTauTau_2017.root", "GluGluZH_2017.root", "ttHToEE_2017.root", "ttHTo2L2Nu_2017.root", "ttHJetToNonbb_2017.root", "TWZToLL_2017.root", "HZJ_HToWWTo2L2Nu_ZTo2L_2017.root"})},
            {"signal",  signalList(year, {"HppM500_2017.root", "HppM600_2017.root", "HppM700_2017.root", "HppM800_2017.root", "HppM900_2017.root", "HppM1000_2017.root", "HppM1100_2017.root", "HppM1200_2017.root", "HppM1300_2017.root", "HppM1400_2017.root", "HppM1500_2017.root"})},
            {"data",    dataList(year, {"SingleElectronB_2017.root", "SingleElectronC_2017.root", "SingleElectronD_2017.root", "SingleElectronE_2017.root", "SingleElectronF_2017.root", "SingleMuonB_2017.root", "SingleMuonC_2017.root", "SingleMuonD_2017.root", "SingleMuonE_2017.root", "SingleMuonF_2017.root"})}
        };
    }

    else if (year == "2018") {
        files = {
            {"DY10_50", mcList(year, {"DYJetsToLLM10to50_2018.root"})},
            {"DY",      mcList(year, {"DYJetsToLLM50_2018.root"})},
            {"WZ",      mcList(year, {"WZTo2Q2L_2018.root", "WZTo3LNu_2018.root"})},
            {"WW",      mcList(year, {"WW_2018.root", "WWTo2L2Nu_2018.root"})},
            {"VVV",     mcList(year, {"WWW_2018.root", "WZZ_2018.root", "ZZZ_2018.root"})},
            {"ttV",     mcList(year, {"ttWJets_2018.root", "ttZJets_2018.root"})},
            {"WJ",      mcList(year, {"WJetsToLNu_NLO_2018.root"})},
            {"ZZ",      mcList(year, {"ZZTo2L2Nu_2018.root", "ZZTo2Q2L_2018.root", "ZZTo4L_2018.root"})},
            {"ST",      mcList(year, {"ST_s-channel_2018.root", "ST_t-channel_antitop_2018.root", "ST_t-channel_top_2018.root", "ST_tW_antitop_2018.root", "ST_tW_top_2018.root"})},
            {"TTbar",   mcList(year, {"TTTo2L2Nu_2018.root", "TTToSemiLeptonic_2018.root", "TTToHadronic_2018.root"})},
            {"QCD",     mcList(year, {"QCD_HT50to100_2018.root", "QCD_HT100to200_2018.root", "QCD_HT200to300_2018.root", "QCD_HT300to500_2018.root", "QCD_HT500to700_2018.root", "QCD_HT700to1000_2018.root", "QCD_HT1500to2000_2018.root", "QCD_HT2000toInf_2018.root"})},
            {"other",   mcList(year, {"ttHToTauTau_2018.root", "ZHToMuMu_2018.root", "ZHToTauTau_2018.root", "GluGluZH_2018.root", "ttHToEE_2018.root", "ttHTo2L2Nu_2018.root", "ttHJetToNonbb_2018.root", "TWZToLL_2018.root", "HZJ_HToWWTo2L2Nu_ZTo2L_2018.root"})},
            {"signal",  signalList(year, {"HppM500_2018.root", "HppM600_2018.root", "HppM700_2018.root", "HppM800_2018.root", "HppM900_2018.root", "HppM1000_2018.root", "HppM1100_2018.root", "HppM1200_2018.root", "HppM1300_2018.root", "HppM1400_2018.root"})},
            {"data",    dataList(year, {"EGammaA_2018.root", "EGammaB_2018.root", "EGammaC_2018.root", "EGammaD_2018.root", "SingleMuonA_2018.root", "SingleMuonB_2018.root", "SingleMuonC_2018.root", "SingleMuonD_2018.root"})}
        };
    }

    else if (year == "Run2") {
        files = {
            {"DY10_50", {
                mcPath("2016preVFP",  "DYJetsToLLM10to50_2016preVFP.root"),
                mcPath("2016postVFP", "DYJetsToLLM10to50_2016postVFP.root"),
                mcPath("2017",        "DYJetsToLLM10to50_2017.root"),
                mcPath("2018",        "DYJetsToLLM10to50_2018.root")
            }},
            {"DY", {
                mcPath("2016preVFP",  "DYJetsToLLM50_2016preVFP.root"),
                mcPath("2016postVFP", "DYJetsToLLM50_2016postVFP.root"),
                mcPath("2017",        "DYJetsToLLM50_2017.root"),
                mcPath("2018",        "DYJetsToLLM50_2018.root")
            }},
            {"WJ", {
                mcPath("2016preVFP",  "WJetsToLNu_NLO_2016preVFP.root"),
                mcPath("2016postVFP", "WJetsToLNu_NLO_2016postVFP.root"),
                mcPath("2017",        "WJetsToLNu_NLO_2017.root"),
                mcPath("2018",        "WJetsToLNu_NLO_2018.root")
            }},
            {"WZ", {
                mcPath("2016preVFP",  "WZTo2Q2L_2016preVFP.root"),  mcPath("2016preVFP",  "WZTo3LNu_2016preVFP.root"),
                mcPath("2016postVFP", "WZTo2Q2L_2016postVFP.root"), mcPath("2016postVFP", "WZTo3LNu_2016postVFP.root"),
                mcPath("2017",        "WZTo2Q2L_2017.root"),        mcPath("2017",        "WZTo3LNu_2017.root"),
                mcPath("2018",        "WZTo2Q2L_2018.root"),        mcPath("2018",        "WZTo3LNu_2018.root")
            }},
            {"WW", {
                mcPath("2016preVFP",  "WW_2016preVFP.root"),  mcPath("2016preVFP",  "WWTo2L2Nu_2016preVFP.root"),
                mcPath("2016postVFP", "WW_2016postVFP.root"), mcPath("2016postVFP", "WWTo2L2Nu_2016postVFP.root"),
                mcPath("2017",        "WW_2017.root"),        mcPath("2017",        "WWTo2L2Nu_2017.root"),
                mcPath("2018",        "WW_2018.root"),        mcPath("2018",        "WWTo2L2Nu_2018.root")
            }},
            {"VVV", {
                mcPath("2016preVFP",  "WWW_2016preVFP.root"),  mcPath("2016preVFP",  "WZZ_2016preVFP.root"),  mcPath("2016preVFP",  "ZZZ_2016preVFP.root"),
                mcPath("2016postVFP", "WWW_2016postVFP.root"), mcPath("2016postVFP", "WZZ_2016postVFP.root"), mcPath("2016postVFP", "ZZZ_2016postVFP.root"),
                mcPath("2017",        "WWW_2017.root"),        mcPath("2017",        "WZZ_2017.root"),        mcPath("2017",        "ZZZ_2017.root"),
                mcPath("2018",        "WWW_2018.root"),        mcPath("2018",        "WZZ_2018.root"),        mcPath("2018",        "ZZZ_2018.root")
            }},
            {"ttV", {
                mcPath("2016preVFP",  "ttWJets_2016preVFP.root"),  mcPath("2016preVFP",  "ttZJets_2016preVFP.root"),
                mcPath("2016postVFP", "ttWJets_2016postVFP.root"), mcPath("2016postVFP", "ttZJets_2016postVFP.root"),
                mcPath("2017",        "ttWJets_2017.root"),        mcPath("2017",        "ttZJets_2017.root"),
                mcPath("2018",        "ttWJets_2018.root"),        mcPath("2018",        "ttZJets_2018.root")
            }},
            {"ZZ", {
                mcPath("2016preVFP",  "ZZTo2L2Nu_2016preVFP.root"),  mcPath("2016preVFP",  "ZZTo2Q2L_2016preVFP.root"),  mcPath("2016preVFP",  "ZZTo4L_2016preVFP.root"),
                mcPath("2016postVFP", "ZZTo2L2Nu_2016postVFP.root"), mcPath("2016postVFP", "ZZTo2Q2L_2016postVFP.root"), mcPath("2016postVFP", "ZZTo4L_2016postVFP.root"),
                mcPath("2017",        "ZZTo2L2Nu_2017.root"),        mcPath("2017",        "ZZTo2Q2L_2017.root"),        mcPath("2017",        "ZZTo4L_2017.root"),
                mcPath("2018",        "ZZTo2L2Nu_2018.root"),        mcPath("2018",        "ZZTo2Q2L_2018.root"),        mcPath("2018",        "ZZTo4L_2018.root")
            }},
            {"ST", {
                mcPath("2016preVFP",  "ST_s-channel_2016preVFP.root"),  mcPath("2016preVFP",  "ST_t-channel_antitop_2016preVFP.root"),  mcPath("2016preVFP",  "ST_t-channel_top_2016preVFP.root"),  mcPath("2016preVFP",  "ST_tW_antitop_2016preVFP.root"),  mcPath("2016preVFP",  "ST_tW_top_2016preVFP.root"),
                mcPath("2016postVFP", "ST_s-channel_2016postVFP.root"), mcPath("2016postVFP", "ST_t-channel_antitop_2016postVFP.root"), mcPath("2016postVFP", "ST_t-channel_top_2016postVFP.root"), mcPath("2016postVFP", "ST_tW_antitop_2016postVFP.root"), mcPath("2016postVFP", "ST_tW_top_2016postVFP.root"),
                mcPath("2017",        "ST_s-channel_2017.root"),        mcPath("2017",        "ST_t-channel_antitop_2017.root"),        mcPath("2017",        "ST_t-channel_top_2017.root"),        mcPath("2017",        "ST_tW_antitop_2017.root"),        mcPath("2017",        "ST_tW_top_2017.root"),
                mcPath("2018",        "ST_s-channel_2018.root"),        mcPath("2018",        "ST_t-channel_antitop_2018.root"),        mcPath("2018",        "ST_t-channel_top_2018.root"),        mcPath("2018",        "ST_tW_antitop_2018.root"),        mcPath("2018",        "ST_tW_top_2018.root")
            }},
            {"TTbar", {
                mcPath("2016preVFP",  "TTTo2L2Nu_2016preVFP.root"),  mcPath("2016preVFP",  "TTToSemiLeptonic_2016preVFP.root"),  mcPath("2016preVFP",  "TTToHadronic_2016preVFP.root"),
                mcPath("2016postVFP", "TTTo2L2Nu_2016postVFP.root"), mcPath("2016postVFP", "TTToSemiLeptonic_2016postVFP.root"), mcPath("2016postVFP", "TTToHadronic_2016postVFP.root"),
                mcPath("2017",        "TTTo2L2Nu_2017.root"),        mcPath("2017",        "TTToSemiLeptonic_2017.root"),        mcPath("2017",        "TTToHadronic_2017.root"),
                mcPath("2018",        "TTTo2L2Nu_2018.root"),        mcPath("2018",        "TTToSemiLeptonic_2018.root"),        mcPath("2018",        "TTToHadronic_2018.root")
            }},
            {"QCD", {
                mcPath("2016preVFP",  "QCD_HT50to100_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT100to200_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT200to300_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT300to500_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT500to700_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT700to1000_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT1000to1500_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT1500to2000_2016preVFP.root"),  mcPath("2016preVFP",  "QCD_HT2000toInf_2016preVFP.root"),
                mcPath("2016postVFP", "QCD_HT50to100_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT100to200_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT200to300_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT300to500_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT500to700_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT700to1000_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT1000to1500_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT1500to2000_2016postVFP.root"), mcPath("2016postVFP", "QCD_HT2000toInf_2016postVFP.root"),
                mcPath("2017",        "QCD_HT50to100_2017.root"),        mcPath("2017",        "QCD_HT100to200_2017.root"),        mcPath("2017",        "QCD_HT200to300_2017.root"),        mcPath("2017",        "QCD_HT300to500_2017.root"),        mcPath("2017",        "QCD_HT500to700_2017.root"),        mcPath("2017",        "QCD_HT700to1000_2017.root"),        mcPath("2017",        "QCD_HT1000to1500_2017.root"),        mcPath("2017",        "QCD_HT1500to2000_2017.root"),        mcPath("2017",        "QCD_HT2000toInf_2017.root"),
                mcPath("2018",        "QCD_HT50to100_2018.root"),        mcPath("2018",        "QCD_HT100to200_2018.root"),        mcPath("2018",        "QCD_HT200to300_2018.root"),        mcPath("2018",        "QCD_HT300to500_2018.root"),        mcPath("2018",        "QCD_HT500to700_2018.root"),        mcPath("2018",        "QCD_HT700to1000_2018.root"),        mcPath("2018",        "QCD_HT1500to2000_2018.root"),        mcPath("2018",        "QCD_HT2000toInf_2018.root")
            }},
            {"other", {
                mcPath("2016preVFP",  "ttHToTauTau_2016preVFP.root"),  mcPath("2016preVFP",  "ZHToMuMu_2016preVFP.root"),  mcPath("2016preVFP",  "ZHToTauTau_2016preVFP.root"),  mcPath("2016preVFP",  "GluGluZH_2016preVFP.root"),  mcPath("2016preVFP",  "ttHToEE_2016preVFP.root"),  mcPath("2016preVFP",  "ttHTo2L2Nu_2016preVFP.root"),  mcPath("2016preVFP",  "ttHJetToNonbb_2016preVFP.root"),  mcPath("2016preVFP",  "TWZToLL_2016preVFP.root"),  mcPath("2016preVFP",  "HZJ_HToWWTo2L2Nu_ZTo2L_2016preVFP.root"),
                mcPath("2016postVFP", "ttHToTauTau_2016postVFP.root"), mcPath("2016postVFP", "ZHToMuMu_2016postVFP.root"), mcPath("2016postVFP", "ZHToTauTau_ext_2016postVFP.root"), mcPath("2016postVFP", "GluGluZH_2016postVFP.root"), mcPath("2016postVFP", "ttHToEE_2016postVFP.root"), mcPath("2016postVFP", "ttHTo2L2Nu_2016postVFP.root"), mcPath("2016postVFP", "ttHJetToNonbb_2016postVFP.root"), mcPath("2016postVFP", "TWZToLL_2016postVFP.root"), mcPath("2016postVFP", "HZJ_HToWWTo2L2Nu_ZTo2L_2016postVFP.root"),
                mcPath("2017",        "ttHToTauTau_2017.root"),        mcPath("2017",        "ZHToMuMu_2017.root"),        mcPath("2017",        "ZHToTauTau_2017.root"),        mcPath("2017",        "GluGluZH_2017.root"),        mcPath("2017",        "ttHToEE_2017.root"),        mcPath("2017",        "ttHTo2L2Nu_2017.root"),        mcPath("2017",        "ttHJetToNonbb_2017.root"),        mcPath("2017",        "TWZToLL_2017.root"),        mcPath("2017",        "HZJ_HToWWTo2L2Nu_ZTo2L_2017.root"),
                mcPath("2018",        "ttHToTauTau_2018.root"),        mcPath("2018",        "ZHToMuMu_2018.root"),        mcPath("2018",        "ZHToTauTau_2018.root"),        mcPath("2018",        "GluGluZH_2018.root"),        mcPath("2018",        "ttHToEE_2018.root"),        mcPath("2018",        "ttHTo2L2Nu_2018.root"),        mcPath("2018",        "ttHJetToNonbb_2018.root"),        mcPath("2018",        "TWZToLL_2018.root"),        mcPath("2018",        "HZJ_HToWWTo2L2Nu_ZTo2L_2018.root")
            }},
            {"signal", {
                signalPath("2016preVFP",  "HppM500_2016preVFP.root"),  signalPath("2016preVFP",  "HppM600_2016preVFP.root"),  signalPath("2016preVFP",  "HppM700_2016preVFP.root"),  signalPath("2016preVFP",  "HppM800_2016preVFP.root"),  signalPath("2016preVFP",  "HppM900_2016preVFP.root"),  signalPath("2016preVFP",  "HppM1000_2016preVFP.root"),  signalPath("2016preVFP",  "HppM1100_2016preVFP.root"),  signalPath("2016preVFP",  "HppM1200_2016preVFP.root"),  signalPath("2016preVFP",  "HppM1300_2016preVFP.root"),  signalPath("2016preVFP",  "HppM1400_2016preVFP.root"),  signalPath("2016preVFP",  "HppM1500_2016preVFP.root"),
                signalPath("2016postVFP", "HppM500_2016postVFP.root"), signalPath("2016postVFP", "HppM600_2016postVFP.root"), signalPath("2016postVFP", "HppM700_2016postVFP.root"), signalPath("2016postVFP", "HppM800_2016postVFP.root"), signalPath("2016postVFP", "HppM900_2016postVFP.root"), signalPath("2016postVFP", "HppM1000_2016postVFP.root"), signalPath("2016postVFP", "HppM1100_2016postVFP.root"), signalPath("2016postVFP", "HppM1200_2016postVFP.root"), signalPath("2016postVFP", "HppM1300_2016postVFP.root"), signalPath("2016postVFP", "HppM1400_2016postVFP.root"), signalPath("2016postVFP", "HppM1500_2016postVFP.root"),
                signalPath("2017",        "HppM500_2017.root"),        signalPath("2017",        "HppM600_2017.root"),        signalPath("2017",        "HppM700_2017.root"),        signalPath("2017",        "HppM800_2017.root"),        signalPath("2017",        "HppM900_2017.root"),        signalPath("2017",        "HppM1000_2017.root"),        signalPath("2017",        "HppM1100_2017.root"),        signalPath("2017",        "HppM1200_2017.root"),        signalPath("2017",        "HppM1300_2017.root"),        signalPath("2017",        "HppM1400_2017.root"),        signalPath("2017",        "HppM1500_2017.root"),
                signalPath("2018",        "HppM500_2018.root"),        signalPath("2018",        "HppM600_2018.root"),        signalPath("2018",        "HppM700_2018.root"),        signalPath("2018",        "HppM800_2018.root"),        signalPath("2018",        "HppM900_2018.root"),        signalPath("2018",        "HppM1000_2018.root"),        signalPath("2018",        "HppM1100_2018.root"),        signalPath("2018",        "HppM1200_2018.root"),        signalPath("2018",        "HppM1300_2018.root"),        signalPath("2018",        "HppM1400_2018.root")
            }},
            {"data", {
                dataPath("2016preVFP",  "SingleElectronB_2016preVFP.root"),  dataPath("2016preVFP",  "SingleElectronC_2016preVFP.root"),  dataPath("2016preVFP",  "SingleElectronD_2016preVFP.root"),  dataPath("2016preVFP",  "SingleElectronE_2016preVFP.root"),  dataPath("2016preVFP",  "SingleElectronF_2016preVFP.root"),  dataPath("2016preVFP",  "SingleMuonB_2016preVFP.root"),  dataPath("2016preVFP",  "SingleMuonC_2016preVFP.root"),  dataPath("2016preVFP",  "SingleMuonD_2016preVFP.root"),  dataPath("2016preVFP",  "SingleMuonE_2016preVFP.root"),  dataPath("2016preVFP",  "SingleMuonF_2016preVFP.root"),
                dataPath("2016postVFP", "SingleElectronF_2016postVFP.root"), dataPath("2016postVFP", "SingleElectronG_2016postVFP.root"), dataPath("2016postVFP", "SingleElectronH_2016postVFP.root"), dataPath("2016postVFP", "SingleMuonF_2016postVFP.root"), dataPath("2016postVFP", "SingleMuonG_2016postVFP.root"), dataPath("2016postVFP", "SingleMuonH_2016postVFP.root"),
                dataPath("2017",        "SingleElectronB_2017.root"),        dataPath("2017",        "SingleElectronC_2017.root"),        dataPath("2017",        "SingleElectronD_2017.root"),        dataPath("2017",        "SingleElectronE_2017.root"),        dataPath("2017",        "SingleElectronF_2017.root"),        dataPath("2017",        "SingleMuonB_2017.root"),        dataPath("2017",        "SingleMuonC_2017.root"),        dataPath("2017",        "SingleMuonD_2017.root"),        dataPath("2017",        "SingleMuonE_2017.root"),        dataPath("2017",        "SingleMuonF_2017.root"),
                dataPath("2018",        "EGammaA_2018.root"),               dataPath("2018",        "EGammaB_2018.root"),               dataPath("2018",        "EGammaC_2018.root"),               dataPath("2018",        "EGammaD_2018.root"),               dataPath("2018",        "SingleMuonA_2018.root"),        dataPath("2018",        "SingleMuonB_2018.root"),        dataPath("2018",        "SingleMuonC_2018.root"),        dataPath("2018",        "SingleMuonD_2018.root")
            }}
        };
    }

    return files;
}

#endif
