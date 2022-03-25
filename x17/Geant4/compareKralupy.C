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
		h_mcHist->Scale(1/(mcExpTime));
	else
		h_mcHist->Scale(1/h_mcHist->Integral());
	h_mcHist->SetLineColor(kRed);
	h_mcHist->SetMarkerColor(kRed);
	// h_mcHist->Scale(1/(realDataDuration*24*3600));
	TH1D* h_dataHist = GetHists(data,histoName);
	if (mcExpTime != -1 && dataExpTime != -1 )
		h_dataHist->Scale(1/(dataExpTime));
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

THStack* GetHistStack(TString mcBkg, TString data, TString mcSig, TString mcSig2, TString mcAstro, TString histoName, Double_t mcBkgExpTime = 1, Double_t dataExpTime = 1, Double_t mcSigExpTime = 1, Double_t mcSig2ExpTime = 1, Double_t mcAstroExpTime = 1)
{
	TH1D* h_mcSigHist = GetHists(mcSig,histoName);
	if (mcSigExpTime != -1 && dataExpTime != -1 )
		h_mcSigHist->Scale(1/(mcSigExpTime/365));
	else
		h_mcSigHist->Scale(1/h_mcSigHist->Integral());
	h_mcSigHist->SetLineColor(kBlue);
	h_mcSigHist->SetMarkerColor(kBlue);

	TH1D* h_mcSig2Hist = GetHists(mcSig2,histoName);
	if (mcSig2ExpTime != -1 && dataExpTime != -1 )
		h_mcSig2Hist->Scale(1/(mcSig2ExpTime/365));
	else
		h_mcSig2Hist->Scale(1/h_mcSig2Hist->Integral());
	h_mcSig2Hist->SetLineColor(kGreen);
	h_mcSig2Hist->SetMarkerColor(kGreen);

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

	TH1D* h_mcAstro = GetHists(mcAstro,histoName);
	if (mcSigExpTime != -1 && dataExpTime != -1 )
		h_mcAstro->Scale(1/(mcAstroExpTime/365));
	else
		h_mcAstro->Scale(1/h_mcAstro->Integral());
	h_mcAstro->SetLineColor(kViolet);
	h_mcAstro->SetMarkerColor(kViolet);


	// cout << "Data: " << h_dataHist->Integral() << " MC: " << h_mcHist->Integral() << endl;
	TString title = ";";
	title +=h_mcSigHist->GetXaxis()->GetTitle();
	title += ";";
	title += h_mcSigHist->GetYaxis()->GetTitle();
	if (mcSigExpTime == -1 || dataExpTime == -1 )
		title += " (Normalized)";

	auto stack = new THStack(histoName,title);
	stack->Add(h_dataHist);
	stack->Add(h_mcSigHist);
	stack->Add(h_mcSig2Hist);
	stack->Add(h_mcBkgHist);
	stack->Add(h_mcAstro);
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

void DrawStacks(TString mcBkg,TString data,TString mcSig, TString mcSig2, TString mcAstro, Double_t mcBckExpTime,Double_t dataExpTime,Double_t mcSigExpTime, Double_t mcSig2ExpTime, Double_t mcAstroExpTime, vector<TString> histNames)
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
		THStack* hs_thetaRec = GetHistStack(mcBkg,data,mcSig,mcSig2,mcAstro,histNames[i],mcBckExpTime,dataExpTime,mcSigExpTime,mcSig2ExpTime,mcAstroExpTime);
		hs_thetaRec->Draw("nostack");
	}
}

void DrawNormStacks(TString mcBkg,TString data,TString mcSig,TString mcSig2, TString mcAstro,vector<TString> histNames)
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
		THStack* hs_thetaRec = GetHistStack(mcBkg,data,mcSig,mcSig2,mcAstro,histNames[i],-1,-1,-1,-1,-1);
		hs_thetaRec->Draw("nostack");
	}
}

int compareKralupy(TString histName = "", Bool_t onlyNorm = true, Bool_t onlyRatio = true)
{
	gStyle->SetPalette(kRainBow);

	vector<TString> filePaths;

	TString noStep = "/Data/x17/Geant4/results_test_Pos1_NoStep_1M_Elec1MeV_.06mmDiverg.root";
	TString stepLimit = "/Data/x17/Geant4/results_test_Pos1_Step10micron_1M_Elec1MeV_.06mmDiverg.root";

	vector<TString> histNames;
	if (histName == "")
		histNames = vector<TString>{"h_x","h_y","h_z","h_nHits","h_eDep","h_stepLength","h_trackLength"};
	else
		histNames = vector<TString>{histName};

	if (!onlyNorm)
	{
		if (onlyRatio)
			DrawRatios(noStep,stepLimit,1,1,histNames);
		// else
			// DrawStacks(mcBkg,data,mcSig,mcSig2,mcAstro,mcBckExpTime,dataExpTime,mcSigExpTime,mcSig2ExpTime,mcAstroExpTime,histNames);
	}else
	{
		if (onlyRatio)
			DrawNormRatios(noStep,stepLimit,histNames);
		// else
			// DrawNormStacks(mcBkg,data,mcSig,mcSig2,mcAstro,histNames);
	}

	return 0;
}