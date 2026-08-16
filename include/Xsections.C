double XSec(std::string fname){
	if(fname.find("ttHTo2L2Nu") < fname.length()) return 0.5418;
	else if(fname.find("ttHToEE") < fname.length()) return 0;
	else if(fname.find("ttHToMuMu") < fname.length()) return 0.5269*0.000218;
	else if(fname.find("ttHToTauTau") < fname.length()) return 0.5269*0.0627;
	else if(fname.find("ttWJets") < fname.length()) return 0.4611;
	else if(fname.find("ttZJets") < fname.length()) return 0.95;
	else if(fname.find("WJetsToLNu_NLO") < fname.length()) return 67350;
	else if(fname.find("WJetsToLNu_HT-70To100") < fname.length()) return  1264.0*1.1421;
	else if(fname.find("WJetsToLNu_HT-100To200") < fname.length()) return 1256.0*1.1421;
	else if(fname.find("WJetsToLNu_HT-200To400") < fname.length()) return 335.5*1.1421;
	else if(fname.find("WJetsToLNu_HT-400To600") < fname.length()) return 45.25*1.1421;
	else if(fname.find("WJetsToLNu_HT-600To800") < fname.length()) return 10.97*1.1421;
	else if(fname.find("WJetsToLNu_HT-800To1200") < fname.length()) return 4.933*1.1421;
	else if(fname.find("WJetsToLNu_HT-1200To2500") < fname.length()) return 1.160*1.1421;
	else if(fname.find("WJetsToLNu_HT-2500ToInf") < fname.length()) return 0.02624*1.1421;
	else if(fname.find("WWTo2L2Nu") < fname.length() and fname.find("HZJ") > fname.length()) return 12.178;
	else if(fname.find("WWW_") < fname.length()) return 0.2086;
	else if(fname.find("WW_") < fname.length()) return 118.7;
	else if(fname.find("WZTo2Q2L") < fname.length()) return 6.565;
	else if(fname.find("WZTo3LNu") < fname.length()) return 5.257;
	else if(fname.find("WZZ_") < fname.length()) return 0.05709;
	else if(fname.find("WZ_") < fname.length()) return 27.55;
	else if(fname.find("ZHToMuMu") < fname.length()) return 0.877*0.000218;
	else if(fname.find("ZHToTauTau") < fname.length()) return 0.877*0.06256;
	else if(fname.find("ZZTo2L2Nu") < fname.length()) return 0.9738;
	else if(fname.find("ZZTo2Q2L") < fname.length()) return 3.698;
	else if(fname.find("ZZTo4L") < fname.length()) return 1.325;
	else if(fname.find("ZZZ_") < fname.length()) return 0.01476;
	else if(fname.find("GluGluZH_") < fname.length()) return 0.0616;
	else if(fname.find("DYJetsToLLM10to50") < fname.length()) return 18610;
	else if(fname.find("DYJetsToLLM50") < fname.length()) return 6077.22;
	else if(fname.find("ST_s-channel_") < fname.length()) return 10.4;
	else if(fname.find("ST_t-channel_antitop_") < fname.length()) return 80.0;
	else if(fname.find("ST_t-channel_top_") < fname.length()) return 134.2;
	else if(fname.find("ST_tW_antitop_") < fname.length()) return 39.65;
	else if(fname.find("ST_tW_top_") < fname.length()) return 39.65;
	else if(fname.find("TTTo2L2Nu_") < fname.length()) return 88.51;
	else if(fname.find("TTToSemiLeptonic_") < fname.length()) return 366.29;
	else if(fname.find("TTToHadronic_") < fname.length()) return 378.93;
	else if(fname.find("ttHJetToNonbb") < fname.length()) return 0.24111;
	else if(fname.find("TWZToLL") < fname.length()) return 0.001669;
	else if(fname.find("HZJ") < fname.length()) return 0.00177;
	else if(fname.find("HppM500") < fname.length()) return 1.717e-3;
	else if(fname.find("HppM600") < fname.length()) return 6.953e-4;
	else if(fname.find("HppM700") < fname.length()) return 3.051e-4;
	else if(fname.find("HppM800") < fname.length()) return 1.433e-4;
	else if(fname.find("HppM900") < fname.length()) return 6.976e-5;
	else if(fname.find("HppM1000") < fname.length()) return 3.519e-5;
	else if(fname.find("HppM1100") < fname.length()) return 1.824e-5;
	else if(fname.find("HppM1200") < fname.length()) return 9.636e-6;
	else if(fname.find("HppM1300") < fname.length()) return 5.228e-6;
	else if(fname.find("HppM1400") < fname.length()) return 2.857e-6;
	else if(fname.find("HppM1500") < fname.length()) return 1.581e-6;
	else if(fname.find("EGamma") < fname.length()) return 1;
	else if(fname.find("Muon") < fname.length()) return 1;
	else if(fname.find("Tau") < fname.length()) return 1;
	else if(fname.find("Single") < fname.length()) return 1;

	else if(fname.find("QCD_Pt-15To20_MuEn")   < fname.length()) return 2813000;
	else if(fname.find("QCD_Pt-20To30_MuEn")   < fname.length()) return 2538000;
	else if(fname.find("QCD_Pt-30To50_MuEn")   < fname.length()) return 1369000;
	else if(fname.find("QCD_Pt-50To80_MuEn")   < fname.length()) return 378000;
	else if(fname.find("QCD_Pt-80To120_MuEn")  < fname.length()) return 89430;
	else if(fname.find("QCD_Pt-120To170_MuEn") < fname.length()) return 21090;
	else if(fname.find("QCD_Pt-170To300_MuEn") < fname.length()) return 7001;
	else if(fname.find("QCD_Pt-300To470_MuEn") < fname.length()) return 618.9;
	else if(fname.find("QCD_Pt-470To600_MuEn") < fname.length()) return 58.54;
	else if(fname.find("QCD_Pt-600To800_MuEn") < fname.length()) return 18.05;
	else if(fname.find("QCD_Pt-800To1000_MuEn")< fname.length()) return 3.281;
	else if(fname.find("QCD_Pt-1000_MuEn")     < fname.length()) return 1.070;

    else if(fname.find("QCD_Pt-15to20_EMEn") < fname.length()) return 1324000;
    else if(fname.find("QCD_Pt-20to30_EMEn") < fname.length()) return 4896000;
    else if(fname.find("QCD_Pt-30to50_EMEn") < fname.length()) return 6447000;
    else if(fname.find("QCD_Pt-50to80_EMEn") < fname.length()) return 1988000;
    else if(fname.find("QCD_Pt-80to120_EMEn") < fname.length()) return 367500;
    else if(fname.find("QCD_Pt-120to170_EMEn") < fname.length()) return 66590;
    else if(fname.find("QCD_Pt-170to300_EMEn") < fname.length()) return 16620;
    else if(fname.find("QCD_Pt-300toInf_EMEn") < fname.length()) return 1104;

    else if(fname.find("QCD_HT50to100")    < fname.length()) return 187300000.0;
    else if(fname.find("QCD_HT100to200")   < fname.length()) return 2350000.0;
    else if(fname.find("QCD_HT200to300")   < fname.length()) return 1555000.0;
    else if(fname.find("QCD_HT300to500")   < fname.length()) return 324500.0;
    else if(fname.find("QCD_HT500to700")   < fname.length()) return 30310.0;
    else if(fname.find("QCD_HT700to1000")  < fname.length()) return 6444.0;
    else if(fname.find("QCD_HT1000to1500") < fname.length()) return 1127.0;
    else if(fname.find("QCD_HT1500to2000") < fname.length()) return 109.8;
    else if(fname.find("QCD_HT2000toInf")  < fname.length()) return 21.98;
	else{

		return 0;
	}
}

double applyXSec(TFile* ifile){
	double lumi_2016 = 35900, lumi_2017 = 41500, lumi_2018 = 58900.0;
	double xs_weight = 1.0;
	std::string fname = ifile->GetName();
	if(XSec(fname)!=1){
			TH1D* hnevts = (TH1D*)ifile->Get("hnevts");
			if(!hnevts) return xs_weight;
			if (fname.find("2016") < fname.length())
				xs_weight = lumi_2016*XSec(fname)/hnevts->Integral();
			else if (fname.find("2017") < fname.length())
				xs_weight = lumi_2017*XSec(fname)/hnevts->Integral();
			else if (fname.find("2018") < fname.length())
				xs_weight = lumi_2018*XSec(fname)/hnevts->Integral();
	}
	return xs_weight;
}

double XSec_Uncert(std::string fname){
	if(fname.find("HppM500") < fname.length()) return 0.000002732;
	else if(fname.find("HppM600") < fname.length()) return 0.000001099;
	else if(fname.find("HppM700") < fname.length()) return 4.80e-07;
	else if(fname.find("HppM800") < fname.length()) return 2.24e-07;
	else if(fname.find("HppM900") < fname.length()) return 1.09e-07;
	else if(fname.find("HppM1000") < fname.length()) return 5.45e-08;
	else if(fname.find("HppM1100") < fname.length()) return 2.81e-08;
	else if(fname.find("HppM1200") < fname.length()) return 1.48e-08;
	else if(fname.find("HppM1300") < fname.length()) return 7.99e-09;
	else if(fname.find("HppM1400") < fname.length()) return 4.35e-09;
	else if(fname.find("HppM1500") < fname.length()) return 2.40e-09;
	else if(fname.find("ttHTo") < fname.length()) return 6.96;
	else if(fname.find("ttWJets") < fname.length()) return 7.47;
	else if(fname.find("ttZJets") < fname.length()) return 8.22;
	else if(fname.find("WW_") < fname.length()) return 5.77;

	else if(fname.find("WZTo") < fname.length()) return 0.0195477*100/1.03002;
	else if(fname.find("ZHToMuMu") < fname.length()) return 4.1;
	else if(fname.find("ZHToTauTau") < fname.length()) return 4.1;

	else if(fname.find("ZZTo") < fname.length()) return 0.0352558*100/1.30291;
	else if(fname.find("GluGluZH_") < fname.length()) return 4.1;
	else if(fname.find("DYJetsToLLM10to50") < fname.length()) return 0;
	else if(fname.find("DYJetsToLLM50") < fname.length()) return 2.49;
	else if(fname.find("ST_t-channel_antitop_") < fname.length()) return 18.82;
	else if(fname.find("ST_t-channel_top_") < fname.length()) return 14.29;
	else if(fname.find("ST_tW_antitop_") < fname.length()) return 11.05;
	else if(fname.find("ST_tW_top_") < fname.length()) return 11.05;
	else if(fname.find("TTTo") < fname.length()) return 3.19;
	else if(fname.find("WWW") < fname.length()) return 27.12;
	else if(fname.find("WZZ") < fname.length()) return 40.00;

	else if(fname.find("QCD_Pt-15To20_MuEn")   < fname.length()) return 8428;
	else if(fname.find("QCD_Pt-20To30_MuEn")   < fname.length()) return 10720;
	else if(fname.find("QCD_Pt-30To50_MuEn")   < fname.length()) return 7188;
	else if(fname.find("QCD_Pt-50To80_MuEn")   < fname.length()) return 1315;
	else if(fname.find("QCD_Pt-80To120_MuEn")  < fname.length()) return 583.5;
	else if(fname.find("QCD_Pt-120To170_MuEn") < fname.length()) return 72.11;
	else if(fname.find("QCD_Pt-170To300_MuEn") < fname.length()) return 28.96;
	else if(fname.find("QCD_Pt-300To470_MuEn") < fname.length()) return 2.54;
	else if(fname.find("QCD_Pt-470To600_MuEn") < fname.length()) return 0.3101;
	else if(fname.find("QCD_Pt-600To800_MuEn") < fname.length()) return 0.08194;
	else if(fname.find("QCD_Pt-800To1000_MuEn")< fname.length()) return 0.01323;
	else if(fname.find("QCD_Pt-1000_MuEn")     < fname.length()) return 0.004323;

    else if(fname.find("QCD_Pt-15to20_EMEn") < fname.length()) return 4183;
    else if(fname.find("QCD_Pt-20to30_EMEn") < fname.length()) return 15410;
    else if(fname.find("QCD_Pt-30to50_EMEn") < fname.length()) return 19870;
    else if(fname.find("QCD_Pt-50to80_EMEn") < fname.length()) return 5961;
    else if(fname.find("QCD_Pt-80to120_EMEn") < fname.length()) return 1091;
    else if(fname.find("QCD_Pt-120to170_EMEn") < fname.length()) return 196.6;
    else if(fname.find("QCD_Pt-170to300_EMEn") < fname.length()) return 49.1;
    else if(fname.find("QCD_Pt-300toInf_EMEn") < fname.length()) return 3.282;

    else if(fname.find("QCD_HT50to100")    < fname.length()) return 520900.0;
    else if(fname.find("QCD_HT100to200")   < fname.length()) return 68630.0;
    else if(fname.find("QCD_HT200to300")   < fname.length()) return 4626.0;
    else if(fname.find("QCD_HT300to500")   < fname.length()) return 976.7;
    else if(fname.find("QCD_HT500to700")   < fname.length()) return 91.71;
    else if(fname.find("QCD_HT700to1000")  < fname.length()) return 19.56;
    else if(fname.find("QCD_HT1000to1500") < fname.length()) return 3.426;
    else if(fname.find("QCD_HT1500to2000") < fname.length()) return 0.3341;
    else if(fname.find("QCD_HT2000toInf")  < fname.length()) return 0.0669;
	else{

		return 0;
	}
}

double fake_uncert_squared(std::string lep_ch, std::string fname){
	double uncert = 0, uncert_square = 0;
	if(lep_ch == "3lep"){
		if(fname.find("DYJetsToLL") < fname.length()) uncert = 0.03;
		else if(fname.find("TTTo2L2Nu") < fname.length()) uncert = 0.03;
		else if(fname.find("TTToSemi") < fname.length()) uncert = 0.03*0.03;
		else if(fname.find("TTToHadro") < fname.length()) uncert = 0.03*0.03*0.03;
	}
	else if(lep_ch == "4lep"){
		if(fname.find("DYJetsToLL") < fname.length()) uncert = 0.03*0.03;
		else if(fname.find("WZ") < fname.length() and fname.find("TWZ") > fname.length() and fname.find("WZZ") > fname.length()) uncert = 0.03;
		else if(fname.find("TTTo2L2Nu") < fname.length()) uncert = 0.03*.03;
		else if(fname.find("TTToSemi") < fname.length()) uncert = 0.03*0.03*0.3;
		else if(fname.find("TTToHadro") < fname.length()) uncert = 0.03*0.03*0.03*0.3;
	}
	uncert_square = uncert*uncert;
	return uncert_square;
}
