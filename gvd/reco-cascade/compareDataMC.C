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

TRatioPlot* GetRatioPlot(TString mc, TString data, TString histoName, Double_t mcExpTime = 1, Double_t dataExpTime = 1)
{
	TH1D* h_mcHist = GetHists(mc,histoName);
	if (mcExpTime != -1 && dataExpTime != -1 )
		h_mcHist->Scale(1/(mcExpTime/365));
	else
		h_mcHist->Scale(1/h_mcHist->Integral());
	h_mcHist->SetLineColor(kRed);
	h_mcHist->SetMarkerColor(kRed);
	// h_mcHist->Scale(1/(realDataDuration*24*3600));
	TH1D* h_dataHist = GetHists(data,histoName);
	if (mcExpTime != -1 && dataExpTime != -1 )
		h_dataHist->Scale(1/(dataExpTime/365));
	else
		h_dataHist->Scale(1/h_dataHist->Integral());
	h_dataHist->SetLineColor(kBlack);
	h_dataHist->SetMarkerColor(kBlack);
	// h_dataHist->Scale(1/dataExpTime);
	cout << "Data: " << h_dataHist->Integral() << " MC: " << h_mcHist->Integral() << endl;
	auto rp = new TRatioPlot(h_mcHist,h_dataHist);
	// c1->SetTicks(0, 1);
	// rp->Draw();
	return rp;
}

THStack* GetHistStack(TString mcBkg, TString data, TString mcSig, TString histoName, Double_t mcBkgExpTime = 1, Double_t dataExpTime = 1, Double_t mcSigExpTime = 1)
{
	TH1D* h_mcSigHist = GetHists(mcSig,histoName);
	if (mcSigExpTime != -1 && dataExpTime != -1 )
		h_mcSigHist->Scale(1/(mcSigExpTime/365));
	else
		h_mcSigHist->Scale(1/h_mcSigHist->Integral());
	h_mcSigHist->SetLineColor(kBlue);
	h_mcSigHist->SetMarkerColor(kBlue);

	TH1D* h_mcBkgHist = GetHists(mcBkg,histoName);
	if (mcBkgExpTime != -1 && dataExpTime != -1 )
		h_mcBkgHist->Scale(1/(mcBkgExpTime/365));
	else
		h_mcBkgHist->Scale(1/h_mcBkgHist->Integral());
	h_mcBkgHist->SetLineColor(kRed);
	h_mcBkgHist->SetMarkerColor(kRed);

	TH1D* h_dataHist = GetHists(data,histoName);
	if (mcSigExpTime != -1 && dataExpTime != -1 )
		h_dataHist->Scale(1/(dataExpTime/365));
	else
		h_dataHist->Scale(1/h_dataHist->Integral());
	h_dataHist->SetLineColor(kBlack);
	h_dataHist->SetMarkerColor(kBlack);
	// h_dataHist->Scale(1/dataExpTime);
	// cout << "Data: " << h_dataHist->Integral() << " MC: " << h_mcHist->Integral() << endl;
	auto stack = new THStack(histoName,h_dataHist->GetTitle());
	stack->Add(h_dataHist);
	stack->Add(h_mcSigHist);
	stack->Add(h_mcBkgHist);
	return stack;
}

void DrawRatios(TString mcBkg,TString data,Double_t mcBckExpTime,Double_t dataExpTime,vector<TString> histNames)
{
	for (unsigned int i = 0; i < histNames.size(); ++i)
	{
		TString canvasName = histNames[i];
		canvasName.ReplaceAll("h_","c_");
		TString canvasTitle = histNames[i];
		canvasTitle.ReplaceAll("h_","");
		TCanvas* c_nHitsNorm = new TCanvas(canvasName,canvasTitle,800,600);
		TRatioPlot* rp_nHitsNorm = GetRatioPlot(mcBkg,data,histNames[i],mcBckExpTime,dataExpTime);
		rp_nHitsNorm->Draw();
	}
}

void DrawNormRatios(TString mcBkg,TString data,vector<TString> histNames)
{
	for (unsigned int i = 0; i < histNames.size(); ++i)
	{
		TString canvasName = histNames[i];
		canvasName.ReplaceAll("h_","c_");
		canvasName += "Norm";
		TString canvasTitle = histNames[i];
		canvasTitle.ReplaceAll("h_","");
		canvasTitle += "Norm";
		TCanvas* c_nHitsNorm = new TCanvas(canvasName,canvasTitle,800,600);
		TRatioPlot* rp_nHitsNorm = GetRatioPlot(mcBkg,data,histNames[i],-1,-1);
		rp_nHitsNorm->Draw();
	}
}

void DrawStacks(TString mcBkg,TString data,TString mcSig,Double_t mcBckExpTime,Double_t dataExpTime,Double_t mcSigExpTime,vector<TString> histNames)
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
		THStack* hs_thetaRec = GetHistStack(mcBkg,data,mcSig,histNames[i],mcBckExpTime,dataExpTime,mcSigExpTime);
		hs_thetaRec->Draw("nostack");
	}
}

void DrawNormStacks(TString mcBkg,TString data,TString mcSig,vector<TString> histNames)
{
	for (unsigned int i = 0; i < histNames.size(); ++i)
	{
		TString canvasName = histNames[i];
		canvasName.ReplaceAll("h_","c_");
		canvasName += "AllNorm";
		TString canvasTitle = histNames[i];
		canvasTitle.ReplaceAll("h_","");
		canvasTitle += "allNorm";
		TCanvas* c_temp = new TCanvas(canvasName,canvasTitle,800,600);
		THStack* hs_thetaRec = GetHistStack(mcBkg,data,mcSig,histNames[i],-1,-1,-1);
		hs_thetaRec->Draw("nostack");
	}
}

int compareDataMC(Int_t cluster, TString histName = "", Bool_t onlyNorm = true, Bool_t onlyRatio = true)
{
	gStyle->SetPalette(kRainBow);

	TString mcBkg = "";
	TString mcSig = "";
	TString data = "";
	Double_t mcBckExpTime = 0;
	Double_t mcSigExpTime = 0;
	Double_t dataExpTime = 0;

	switch (cluster)
	{
		case 3:
			mcBkg = "/Data/BaikalData/mc/simGVD/cluster3/results_muatm_sep20_c3.reco.cascade.root";
			mcSig = "/Data/BaikalData/mc/ANIS/nueatm/results_nueatm_c5.reco.cascade.root";
			data = "/Data/BaikalData/reco-cascade/2019/cluster3/reco.cascade.results.root";
			mcBckExpTime = 130.3993;
			mcSigExpTime = 109.5;
			dataExpTime = 62.8844;
			break;
		case 4:
			mcBkg = "/Data/BaikalData/mc/simGVD/cluster4/results_muatm_sep20_c4.reco.cascade.root";
			mcSig = "/Data/BaikalData/mc/ANIS/nueatm/results_nueatm_c5.reco.cascade.root";
			data = "/Data/BaikalData/reco-cascade/2019/cluster4/reco.cascade.results.root";
			mcBckExpTime = 130.3993;
			mcSigExpTime = 109.5;
			dataExpTime = 77.5543;
			break;
		case 5:
			mcBkg = "/Data/BaikalData/mc/simGVD/cluster5/results_muatm_sep20_c5.reco.cascade.root";
			mcSig = "/Data/BaikalData/mc/ANIS/nueatm/results_nueatm_c5.reco.cascade.root";
			data = "/Data/BaikalData/reco-cascade/2019/cluster5/reco.cascade.results.root";
			mcBckExpTime = 127.8483;
			mcSigExpTime = 109.5;
			dataExpTime = 61.5483;
			break;
	}

	vector<TString> histNames;
	if (histName == "")
		histNames = vector<TString>{"h_nHits","h_nHitsReco","h_nStringsReco","h_energyRec","h_thetaRec","h_cosThetaRec","h_OMID","h_chi2","h_like","h_likeHitOnly","h_distanceCS","h_posChange","h_qRatio","h_qRatioRed","h_branchRatio","h_nTrackHits","h_qTrack","h_nTrackHitsSeg","h_qTrackSeg","h_qEarly","h_logPHit","h_tRes","h_chi2Caus"};
	else
		histNames = vector<TString>{histName};

	if (!onlyNorm)
	{
		if (onlyRatio)
			DrawRatios(mcBkg,data,mcBckExpTime,dataExpTime,histNames);
		else
			DrawStacks(mcBkg,data,mcSig,mcBckExpTime,dataExpTime,mcSigExpTime,histNames);
	}else
	{
		if (onlyRatio)
			DrawNormRatios(mcBkg,data,histNames);
		else
			DrawNormStacks(mcBkg,data,mcSig,histNames);
	}

	return 0;
}