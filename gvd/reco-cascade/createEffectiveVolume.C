#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

TH1D* GetHist1D(TString folderName, TString histoName)
{
    // vector<TH1D*> v_histograms;
	TFile* file = new TFile(folderName);
	TH1D* h_histo = static_cast<TH1D*>(file->Get(histoName));
	// v_histograms.emplace_back(h_energyMC);

    return h_histo;
}

TH2D* GetHist2D(TString folderName, TString histoName)
{
    // vector<TH1D*> v_histograms;
	TFile* file = new TFile(folderName);
	TH2D* h_histo = static_cast<TH2D*>(file->Get(histoName));
	// v_histograms.emplace_back(h_energyMC);

    return h_histo;
}

int createEffectiveVolume()
{
	// TString MCGenFile = "/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/effVolume_ANISoutput.root";
	TString MCGenFile = "/Data/BaikalData/mc/ANIS/nueastro_ver4_50kNoise/effVolume_ANISoutput.root";
	// TString MCDetFile = "/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/results_nueatm_c3_ver4_50kNoise.reco.cascade.root";
	TString MCDetFile = "/Data/BaikalData/mc/ANIS/nueastro_ver4_50kNoise/results_nueastro_c3_ver4_50kNoise.reco.cascade.root";
	TString MCDetFile63 = "/Data/BaikalData/mc/ANIS/nueastro_ver4_50kNoise/results_i0000.reco.cascade.root";
	TString MCDetFileDetResp = "/Data/BaikalData/mc/ANIS/nueastro_ver4_50kNoise/effVolumeDetResp_astro_nue_ver4.dat.root";

	// TH2D* h_MCGen2D = GetHist2D(MCGenFile,"h_genEventsVslogEVsCosTheta");
	TH2D* h_MCGen2D = GetHist2D(MCGenFile,"h_genEventsVslogEVsCosThetaWeighted");
	TCanvas* c_MCGen2D = new TCanvas("c_MCGen2D","MCGen2D",800,600);
	h_MCGen2D->Draw("colz");
	TH1D* h_MCGen2Dprojection = h_MCGen2D->ProjectionX("Gen_pxFull");
	TCanvas* c_MCGen2Dprojection = new TCanvas("c_MCGen2Dprojection","MCGen2DProjection",800,600);
	h_MCGen2Dprojection->Draw();
	TH1D* h_MCGen2DprojectionHalf = h_MCGen2D->ProjectionX("Gen_pxUp",11,20);
	TCanvas* c_MCGen2DprojectionHalf = new TCanvas("c_MCGen2DprojectionHalf","MCGen2DProjectionHalf",800,600);
	h_MCGen2DprojectionHalf->Draw();
	TH1D* h_MCGen2DprojectionDownHalf = h_MCGen2D->ProjectionX("Gen_pxDown",1,10);
	TCanvas* c_MCGen2DprojectionDownHalf = new TCanvas("c_MCGen2DprojectionDownHalf","MCGen2DProjectionDownHalf",800,600);
	h_MCGen2DprojectionDownHalf->Draw();
	TH1D* h_MCGen2DprojectionY = h_MCGen2D->ProjectionY("Gen_py",1,40);
	TCanvas* c_MCGen2DprojectionY = new TCanvas("c_MCGen2DprojectionY","MCGen2DProjectionY",800,600);
	h_MCGen2DprojectionY->Draw();
	// h_MCGen2Dprojection->Scale(V);
	// h_MCGen2Dprojection->SetTitle(";neutrino energy log_{10}(E [TeV]);#nu_{e} effective volume [km^{3}]");

	// TH2D* h_MCDet2D = GetHist2D(MCDetFile,"h_genEventsVslogEVsCosThetaWeighted");
	TH2D* h_MCDet2D = GetHist2D(MCDetFile,"h_genEventsVslogEVsCosThetaWeighted");
	TCanvas* c_MCDet2D = new TCanvas("c_MCDet2D","MCDet2D",800,600);
	h_MCDet2D->Draw("colz");
	TH1D* h_MCDet2Dprojection = h_MCDet2D->ProjectionX("Det_pxFull");
	TCanvas* c_MCDet2Dprojection = new TCanvas("c_MCDet2Dprojection","MCDet2DProjection",800,600);
	h_MCDet2Dprojection->Draw();
	TH1D* h_MCDet2DprojectionHalf = h_MCDet2D->ProjectionX("Det_pxUp",11,20);
	TCanvas* c_MCDet2DprojectionHalf = new TCanvas("c_MCDet2DprojectionHalf","MCDet2DProjectionHalf",800,600);
	h_MCDet2DprojectionHalf->Draw();
	TH1D* h_MCDet2DprojectionDownHalf = h_MCDet2D->ProjectionX("Det_pxDown",1,10);
	TCanvas* c_MCDet2DprojectionDownHalf = new TCanvas("c_MCDet2DprojectionDownHalf","MCDet2DProjectionDownHalf",800,600);
	h_MCDet2DprojectionDownHalf->Draw();
	TH1D* h_MCDet2DprojectionY = h_MCDet2D->ProjectionY("Det_py",1,40);
	TCanvas* c_MCDet2DprojectionY = new TCanvas("c_MCDet2DprojectionY","MCDet2DProjectionY",800,600);
	h_MCDet2DprojectionY->Draw();

	TH2D* h_MCDet632D = GetHist2D(MCDetFile63,"h_genEventsVslogEVsCosThetaWeighted");
	TCanvas* c_MCDet632D = new TCanvas("c_MCDet632D","MCDet632D",800,600);
	h_MCDet632D->Draw("colz");
	TH1D* h_MCDet632Dprojection = h_MCDet632D->ProjectionX("Det63_pxFull");
	TCanvas* c_MCDet632Dprojection = new TCanvas("c_MCDet632Dprojection","MCDet632DProjection",800,600);
	h_MCDet632Dprojection->Draw();

	TH2D* h_MCDetResp2D = GetHist2D(MCDetFileDetResp,"h_genEventsVslogEVsCosThetaWeighted");
	TCanvas* c_MCDetResp2D = new TCanvas("c_MCDetResp2D","MCDetResp2D",800,600);
	h_MCDetResp2D->Draw("colz");
	TH1D* h_MCDetResp2Dprojection = h_MCDetResp2D->ProjectionX("DetResp_pxFull");
	TCanvas* c_MCDetResp2Dprojection = new TCanvas("c_MCDetResp2Dprojection","MCDetResp2DProjection",800,600);
	h_MCDetResp2Dprojection->Draw();

	double V = 0.264; //km^3
	int nClusters = 8;

	TH2D* h_ratio = static_cast<TH2D*>(h_MCDet2D->Clone());
	h_ratio->Divide(h_MCGen2D);
	h_ratio->Scale(V*nClusters);
	TCanvas* c_ratio = new TCanvas("c_ratio","Ratio",800,600);
	h_ratio->Draw("colz");

	TH1D* h_ratioCorrect = static_cast<TH1D*>(h_MCDet2Dprojection->Clone());
	h_ratioCorrect->Divide(h_MCGen2Dprojection);
	h_ratioCorrect->Scale(V*nClusters);
	TCanvas* c_ratioCorrect = new TCanvas("c_ratioCorrect","RatioCorrect",800,600);
	h_ratioCorrect->Draw("colz");

	TH1D* h_ratioCorrectY = static_cast<TH1D*>(h_MCDet2DprojectionY->Clone());
	h_ratioCorrectY->Divide(h_MCGen2DprojectionY);
	h_ratioCorrectY->Scale(V*nClusters);
	TCanvas* c_ratioCorrectY = new TCanvas("c_ratioCorrectY","RatioCorrectY",800,600);
	h_ratioCorrectY->Draw("colz");

	TH1D* h_ratioCorrectHalf = static_cast<TH1D*>(h_MCDet2DprojectionHalf->Clone());
	h_ratioCorrectHalf->Divide(h_MCGen2DprojectionHalf);
	h_ratioCorrectHalf->Scale(V*nClusters);
	TCanvas* c_ratioCorrectHalf = new TCanvas("c_ratioCorrectHalf","RatioCorrectHalf",800,600);
	h_ratioCorrectHalf->Draw("colz");

	TH1D* h_ratioCorrectDownHalf = static_cast<TH1D*>(h_MCDet2DprojectionDownHalf->Clone());
	h_ratioCorrectDownHalf->Divide(h_MCGen2DprojectionDownHalf);
	h_ratioCorrectDownHalf->Scale(V*nClusters);
	TCanvas* c_ratioCorrectDownHalf = new TCanvas("c_ratioCorrectDownHalf","RatioCorrectDownHalf",800,600);
	h_ratioCorrectDownHalf->Draw("colz");

	TH1D* h_ratioCorrect63 = static_cast<TH1D*>(h_MCDet632Dprojection->Clone());
	h_ratioCorrect63->Divide(h_MCGen2Dprojection);
	h_ratioCorrect63->Scale(V*nClusters);
	TCanvas* c_ratio63Correct = new TCanvas("c_ratio63Correct","Ratio63Correct",800,600);
	h_ratioCorrect63->Draw("colz");

	TH1D* h_ratioCorrectDetResp = static_cast<TH1D*>(h_MCDetResp2Dprojection->Clone());
	h_ratioCorrectDetResp->Divide(h_MCGen2Dprojection);
	h_ratioCorrectDetResp->Scale(V*nClusters);
	TCanvas* c_ratioDetRespCorrect = new TCanvas("c_ratioDetRespCorrect","RatioDetRespCorrect",800,600);
	h_ratioCorrectDetResp->Draw("colz");

	TGraph* g_projection = new TGraph();
	Int_t nPoints = 0;

	TGraph* g_projectionHalf = new TGraph();
	Int_t nPointsHalf = 0;

	TGraph* g_projectionDownHalf = new TGraph();
	Int_t nPointsDownHalf = 0;

	TGraph* g_projection63 = new TGraph();
	Int_t nPoints63 = 0;

	TGraph* g_projectionDetResp = new TGraph();
	Int_t nPointsDetResp = 0;



	for (int i = 0; i < h_ratioCorrect->GetNbinsX(); ++i)
	{
		if (h_ratioCorrect->GetBinContent(i) != 0)
		{
			g_projection->SetPoint(nPoints,TMath::Power(10,h_ratioCorrect->GetBinCenter(i))*1000,h_ratioCorrect->GetBinContent(i));
			// std::cout << i << "\t" << h_ratioCorrect->GetBinCenter(i) << "\t" << h_ratioCorrect->GetBinContent(i) << std::endl;
			nPoints++;
		}
		if (h_ratioCorrect63->GetBinContent(i) != 0)
		{
			g_projection63->SetPoint(nPoints63,TMath::Power(10,h_ratioCorrect63->GetBinCenter(i))*1000,h_ratioCorrect63->GetBinContent(i));
			nPoints63++;
		}
		if (h_ratioCorrectDetResp->GetBinContent(i) != 0)
		{
			g_projectionDetResp->SetPoint(nPointsDetResp,TMath::Power(10,h_ratioCorrectDetResp->GetBinCenter(i))*1000,h_ratioCorrectDetResp->GetBinContent(i));
			nPointsDetResp++;
		}
		if (h_ratioCorrectHalf->GetBinContent(i) != 0)
		{
			g_projectionHalf->SetPoint(nPointsHalf,TMath::Power(10,h_ratioCorrectHalf->GetBinCenter(i))*1000,h_ratioCorrectHalf->GetBinContent(i));
			nPointsHalf++;
		}
		if (h_ratioCorrectDownHalf->GetBinContent(i) != 0)
		{
			g_projectionDownHalf->SetPoint(nPointsDownHalf,TMath::Power(10,h_ratioCorrectDownHalf->GetBinCenter(i))*1000,h_ratioCorrectDownHalf->GetBinContent(i));
			nPointsDownHalf++;
		}
	}


	TCanvas* c_graphProjection = new TCanvas("c_graphProjection","GraphProjection",800,600);
	g_projection->Draw("AL");
	g_projection->SetLineColor(kGreen);
	g_projection->SetLineWidth(2);
	g_projection->SetTitle("Well reco cascades, all sky");
	g_projectionHalf->Draw("SAME L");
	g_projectionHalf->SetLineColor(kOrange);
	g_projectionHalf->SetLineWidth(2);
	g_projectionHalf->SetTitle("Well reco cascades, up-going ");
	g_projectionDownHalf->Draw("SAME L");
	g_projectionDownHalf->SetLineColor(kMagenta);
	g_projectionDownHalf->SetLineWidth(2);
	g_projectionDownHalf->SetTitle("Well reco cascades, down-going ");
	g_projection63->Draw("SAME L");
	g_projection63->SetLineColor(kRed);
	g_projection63->SetLineWidth(2);
	g_projection63->SetTitle("Causality 6/3, all sky");
	g_projectionDetResp->Draw("SAME L");
	g_projectionDetResp->SetLineColor(kBlue);
	g_projectionDetResp->SetLineWidth(2);
	g_projectionDetResp->SetTitle("Detector response (4/3), all sky");
	g_projection->GetXaxis()->SetTitle("neutrino energy [GeV]");
	g_projection->GetYaxis()->SetTitle("effective volume [km^{3}]");
	gPad->SetLogx();
	gPad->SetLogy();
	g_projection->GetXaxis()->SetLimits(1000,1.0e7);
	g_projection->GetHistogram()->SetMinimum(1e-6);
	g_projection->GetHistogram()->SetMaximum(1);
	// g_projection->Draw("AC"); // draw again to apply the axis limits
	c_graphProjection->Update();
	c_graphProjection->BuildLegend(0.5,0.4,0.9,0.6);

	TLatex latex;
	latex.SetTextSize(0.035);
	latex.SetTextAlign(13);  //align at top
	latex.DrawLatex(1.5e3,5e-6,("Baikal-GVD preliminary, cascade channel #nu_{e}, " + to_string(nClusters)+" clusters").c_str());

	return 0;
}