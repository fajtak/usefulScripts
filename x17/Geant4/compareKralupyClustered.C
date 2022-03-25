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
#include <utility>
#include <tuple>

#include "BEvent.h"


TH1D* GetHists(TString folderName, TString histoName)
{
    // vector<TH1D*> v_histograms;
	TFile* file = new TFile(folderName);
	TH1D* h_histo = static_cast<TH1D*>(file->Get(histoName));
	// v_histograms.emplace_back(h_energyMC);

    return h_histo;
}

TRatioPlot* GetRatioPlot(vector<tuple<TString,TString,double>> &filePaths, TString histoName)
{
	TH1D* h_mcHist = GetHists(get<1>(filePaths[1]),histoName);
	if (get<2>(filePaths[1]) > 0 )
		h_mcHist->Scale(1/(get<2>(filePaths[1])));
	else
	{
		if (get<2>(filePaths[1]) == -2)
			h_mcHist->Scale(1/h_mcHist->GetMaximum());
		else
			h_mcHist->Scale(1/h_mcHist->Integral());
	}
	h_mcHist->SetLineColor(kRed);
	h_mcHist->SetMarkerColor(kRed);
	// h_mcHist->Scale(1/(realDataDuration*24*3600));
	TH1D* h_dataHist = GetHists(get<1>(filePaths[0]),histoName);
	if (get<2>(filePaths[0]) > 0 )
		h_dataHist->Scale(1/(get<2>(filePaths[1])));
	else
	{
		if (get<2>(filePaths[0]) == -2)
			h_dataHist->Scale(1/h_dataHist->GetMaximum());
		else
			h_dataHist->Scale(1/h_dataHist->Integral());
	}
		// h_dataHist->Scale(1/h_dataHist->Integral());
	h_dataHist->SetLineColor(kBlack);
	h_dataHist->SetMarkerColor(kBlack);
	// h_dataHist->Scale(1/dataExpTime);
	cout << "Data: " << h_dataHist->Integral() << " MC: " << h_mcHist->Integral() << endl;
	auto rp = new TRatioPlot(h_mcHist,h_dataHist);
	// c1->SetTicks(0, 1);
	// rp->Draw();
	return rp;
}

THStack* GetHistStack(vector<tuple<TString,TString,double>> &filePaths, TString histoName)
{
	auto stack = new THStack(histoName,"");
	double scale = 0;
	for (unsigned int i = 0; i < filePaths.size(); ++i)
	{
		TH1D* h_temp = GetHists(get<1>(filePaths[i]),histoName);
		if (i == 0)
			scale = h_temp->Integral();
		// h_temp->Rebin(4);
		if (get<2>(filePaths[i]) > 0 )
			h_temp->Scale(1/(get<2>(filePaths[i])));
		else
		{
			if (get<2>(filePaths[i]) == -2)
				h_temp->Scale(1/h_temp->GetMaximum());
			else if (get<2>(filePaths[i]) == -3)
				h_temp->Scale(1/h_temp->Integral()*scale);
			else
				h_temp->Scale(1/h_temp->Integral());
		}
		h_temp->SetLineColor(i+1);
		h_temp->SetMarkerColor(i+1);
		h_temp->SetTitle(get<0>(filePaths[i]));
		stack->Add(h_temp);
		if (i == 0)
		{
			TString title = ";";
			title +=h_temp->GetXaxis()->GetTitle();
			title += ";";
			title += h_temp->GetYaxis()->GetTitle();
			if (get<2>(filePaths[i]) == -1)
				title += " (Normalized)";
			stack->SetTitle(title);
		}
	}

	return stack;
}

void DrawRatios(vector<tuple<TString,TString,double>> &filePaths,vector<TString> histNames)
{
	for (unsigned int i = 0; i < histNames.size(); ++i)
	{
		TString canvasName = histNames[i];
		canvasName.ReplaceAll("h_","c_");
		TString canvasTitle = histNames[i];
		canvasTitle.ReplaceAll("h_","");
		TCanvas* c_nHitsNorm = new TCanvas(canvasName,canvasTitle,800,600);
		TRatioPlot* rp_nHitsNorm = GetRatioPlot(filePaths,histNames[i]);
		rp_nHitsNorm->Draw();
	}
}

void DrawNormRatios(vector<tuple<TString,TString,double>> &filePaths,vector<TString> histNames)
{
	// vector<pair<TString,TString>> filePaths;
	for (unsigned int i = 0; i < filePaths.size(); ++i)
	{
		get<2>(filePaths[i]) = -1;
	}

	for (unsigned int i = 0; i < histNames.size(); ++i)
	{
		TString canvasName = histNames[i];
		canvasName.ReplaceAll("h_","c_");
		canvasName += "Norm";
		TString canvasTitle = histNames[i];
		canvasTitle.ReplaceAll("h_","");
		canvasTitle += "Norm";
		TCanvas* c_nHitsNorm = new TCanvas(canvasName,canvasTitle,800,600);
		TRatioPlot* rp_nHitsNorm = GetRatioPlot(filePaths,histNames[i]);
		rp_nHitsNorm->Draw();
	}
}

void DrawStacks(vector<tuple<TString,TString,double>> &filePaths, vector<TString> histNames)
{
	for (unsigned int i = 0; i < histNames.size(); ++i)
	{
		TString canvasName = histNames[i];
		canvasName.ReplaceAll("h_","c_");
		canvasName += "All";
		TString canvasTitle = histNames[i];
		canvasTitle.ReplaceAll("h_","");
		canvasTitle += "all";
		TCanvas* c_temp = new TCanvas(canvasName,canvasTitle,800,600);
		THStack* hs_thetaRec = GetHistStack(filePaths,histNames[i]);
		hs_thetaRec->Draw("nostack");
	}
}

void DrawNormStacks(vector<tuple<TString,TString,double>> &filePaths,vector<TString> histNames)
{
	for (unsigned int i = 0; i < filePaths.size(); ++i)
	{
		get<2>(filePaths[i]) = -1;
	}
	for (unsigned int i = 0; i < histNames.size(); ++i)
	{
		TString canvasName = histNames[i];
		canvasName.ReplaceAll("h_","c_");
		canvasName += "AllNorm";
		TString canvasTitle = histNames[i];
		canvasTitle.ReplaceAll("h_","");
		canvasTitle += "allNorm";
		TCanvas* c_temp = new TCanvas(canvasName,canvasTitle,800,600);
		THStack* hs_thetaRec = GetHistStack(filePaths,histNames[i]);
		hs_thetaRec->Draw("nostack");
		c_temp->BuildLegend();
	}
}

int compareKralupyClustered(TString histName = "", Bool_t onlyNorm = true, Bool_t onlyRatio = true)
{
	gStyle->SetPalette(kRainBow);

	vector<tuple<TString,TString,double>> filePaths;

	// filePaths.push_back(std::make_tuple("Data","/Data/x17/Kralupy/results_Run002_G3_1MeV_90Degree_0.root",1));
	// filePaths.push_back(std::make_tuple("0.08mm","/home/fajtak/work/allpix-squared/output/kralupy_pos1_step10micron_1MeV_0.08mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));
	// filePaths.push_back(std::make_tuple("0.0977mm","/home/fajtak/work/allpix-squared/output/kralupy_pos1_step10micron_1MeV_0.0977mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));
	// filePaths.push_back(std::make_tuple("0.09mm","/home/fajtak/work/allpix-squared/output/kralupy_pos1_step10micron_1MeV_0.09mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));
	// filePaths.push_back(std::make_tuple("0.07mm","/home/fajtak/work/allpix-squared/output/kralupy_pos1_step10micron_1MeV_0.07mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));
	// filePaths.push_back(std::make_tuple("0.06mm","/home/fajtak/work/allpix-squared/output/kralupy_pos1_step10micron_1MeV_0.06mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));
	// filePaths.push_back(std::make_tuple("0.06mm - noStep","/home/fajtak/work/allpix-squared/output/kralupy_pos1_noStep_1MeV_0.06mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));

	filePaths.push_back(std::make_tuple("Data","/Data/x17/Kralupy/results_Run002_G3_1MeV_90Degree_Pos2_0.root",1));
	filePaths.push_back(std::make_tuple("0.08mm","/home/fajtak/work/allpix-squared/output/kralupy_pos2_step10micron_1MeV_0.08mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));
	filePaths.push_back(std::make_tuple("0.06mm","/home/fajtak/work/allpix-squared/output/kralupy_pos2_step10micron_1MeV_0.06mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));
	filePaths.push_back(std::make_tuple("0.09mm","/home/fajtak/work/allpix-squared/output/kralupy_pos2_step10micron_1MeV_0.09mm/ClusterFiles/root/results_data_CSAEmul_noADCres.root",1));


	vector<TString> histNames;
	if (histName == "")
		histNames = vector<TString>{"h_clusterSize","h_pixX","h_pixY","h_clusterMeanX","h_clusterMeanY","h_clusterHeightKeV","h_clusterVolumeKeV","h_clusterVolumeCentroidX","h_clusterVolumeCentroidY","h_deltaToA","h_ToTKeV"};
	else
		histNames = vector<TString>{histName};

	if (!onlyNorm)
	{
		if (onlyRatio)
			DrawRatios(filePaths,histNames);
		else
			DrawStacks(filePaths,histNames);
	}else
	{
		if (onlyRatio)
			DrawNormRatios(filePaths,histNames);
		else
			DrawNormStacks(filePaths,histNames);
	}

	return 0;
}