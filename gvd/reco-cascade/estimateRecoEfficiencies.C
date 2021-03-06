#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRatioPlot.h"

// #include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

#include "BEvent.h"


TH1D* GetHists(TString folderName, TString histoName)
{
    // vector<TH1D*> v_histograms;
	TFile* file = new TFile(folderName);
	TH1D* h_histo = static_cast<TH1D*>(file->Get(histoName));
	// v_histograms.emplace_back(h_energyMC);

    return h_histo;
}

TRatioPlot* GetRatioPlot(TString triggerLevel, TString recoLevel, TString histoName)
{
	TH1D* h_triggerLevelhists = GetHists(triggerLevel,histoName);
	TH1D* h_recoLevelhists = GetHists(recoLevel,histoName);
	auto rp = new TRatioPlot(h_recoLevelhists,h_triggerLevelhists);
	// c1->SetTicks(0, 1);
	// rp->Draw();
	return rp;
}

int estimateRecoEfficiencies(int processID = 0)
{
	gStyle->SetPalette(kRainBow);

	TString triggerLevel = "";
	TString recoLevel = "";

	switch (processID)
	{
		case 0:
			triggerLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver2/results_i0000.reco.cascade.00000000_01000000.root";
			recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver2/results_nueatm_c3.reco.cascade.root";
			break;
		case 1:
			triggerLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver2_50kNoise/results_nueatm_c3_ver2_50kNoise.reco.cascade.root";
			// recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver2_100kNoise/results_nueatm_c3_ver2_100kNoise.reco.cascade.root";
			recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver2_200kNoise/results_nueatm_c3_ver2_200kNoise.reco.cascade.root";
			break;
		case 2:
			triggerLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver3_50kNoise/results_i0000.reco.cascade.00000000_00300000.caus.root";
			// recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver3_50kNoise/results_i0000.reco.cascade.00000000_00300000.root";
			recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver3_50kNoise/results_nueatm_cl3_ver3_50kNoise.reco.cascade.root";
			break;
		case 3:
			triggerLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/results_i0000.reco.cascade.00000000_00243000.caus.root";
			recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/results_i0000.reco.cascade.00000000_00243000.fitPos.root";
			// recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/results_i0000.reco.cascade.00000000_00243000.closeHits.root";
			// recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/results_nueatm_c3_ver4_50kNoise.reco.cascade.root";
			break;
	}

	// TString triggerLevel = "/Data/BaikalData/2019/cluster2/exp/reco/reco1.0-1162-g2262-v1.1/0000/results_i0000.reco.cascade.00002000_00004000.root";
	// TString triggerLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver2/results_i0000.reco.cascade.00000000_00367500.root";
	// TString triggerLevel = "/Data/BaikalData/2019/cluster1/exp/reco/reco1.0-1162-g2262-v1.1/0000/results_i0000.reco.cascade.root";
	// TString recoLevel = "/Data/BaikalData/2019/cluster1/exp/reco/reco1.0-1162-g2262-v1.1/0000/results_i0000.reco.cascade.00002000_00004000.root";
	// TString recoLevel = "/Data/BaikalData/mc/ANIS/nueatm_ver2/results_nueatm_c3.reco.cascade.root";
	// TString recoLevel = "/Data/BaikalData/mc/ANIS/results_i0000.reco.cascade.root";

	TCanvas* c_EfficiencyVsE = new TCanvas("c_EfficiencyVsE","EfficiencyVsE",800,600);
	TRatioPlot* rp_efficiencyVsE = GetRatioPlot(triggerLevel,recoLevel,"h_logEnergyMC");
	rp_efficiencyVsE->Draw();
	c_EfficiencyVsE->BuildLegend();

	TCanvas* c_EfficiencyVsENW = new TCanvas("c_EfficiencyVsENW","EfficiencyVsENW",800,600);
	TRatioPlot* rp_efficiencyVsENW = GetRatioPlot(triggerLevel,recoLevel,"h_logEnergyMCNW");
	rp_efficiencyVsENW->Draw();

	TCanvas* c_EfficiencyVsTheta = new TCanvas("c_EfficiencyVsTheta","EfficiencyVsTheta",800,600);
	TRatioPlot* rp_efficiencyVsTheta = GetRatioPlot(triggerLevel,recoLevel,"h_cosThetaMC");
	rp_efficiencyVsTheta->Draw();

	TCanvas* c_EfficiencyVsThetaNW = new TCanvas("c_EfficiencyVsThetaNW","EfficiencyVsThetaNW",800,600);
	TRatioPlot* rp_efficiencyVsThetaNW = GetRatioPlot(triggerLevel,recoLevel,"h_cosThetaMCNW");
	rp_efficiencyVsThetaNW->Draw();

	TCanvas* c_EfficiencyVsDistCS = new TCanvas("c_EfficiencyVsDistCS","EfficiencyVsDistCS",800,600);
	TRatioPlot* rp_efficiencyVsDistCS = GetRatioPlot(triggerLevel,recoLevel,"h_distanceCSMC");
	rp_efficiencyVsDistCS->Draw();

	TCanvas* c_EfficiencyVsDistCSNW = new TCanvas("c_EfficiencyVsDistCSNW","EfficiencyVsDistCSNW",800,600);
	TRatioPlot* rp_efficiencyVsDistCSNW = GetRatioPlot(triggerLevel,recoLevel,"h_distanceCSMCNW");
	rp_efficiencyVsDistCSNW->Draw();


	// hs_ALow->Draw("plc nostack HIST");
	// gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	return 0;
}