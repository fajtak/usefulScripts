#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "BRecoCascade.h"
#include "BEventMask.h"
#include "BEvent.h"
#include "BJointHeader.h"
#include "BMCCascadeEvent.h"
#include "BRunInfo.h"

double SetTChain(TChain &recCasc, int year, int cluster)
{
	int startID = cluster!=-1?cluster:1;
	int endID = cluster!=-1?cluster+1:10;

	int startSeason = year!=-1?year:16;
	int endSeason = year!=-1?year+1:20+1;

	for (int j = startSeason; j < endSeason; j++)
	{
		for (int i = startID; i < endID; ++i)
		{
			TString filesDir = Form("/Data/BaikalData/reco-cascade/20%d/cluster%d/*.reco.cascade.root",j,i);
			// TString filesDir = Form("/Data/BaikalData/reco-cascade/20%d/cluster%d_beg/*.reco.cascade.root",j,i);
			recCasc.Add(filesDir);
		}
	}
	return recCasc.GetEntries();
}

double SetTChain(TChain &recCasc, TString fileName)
{
	recCasc.Add(fileName);
	return recCasc.GetEntries();
}

double CalculateExpositionTime(TChain &recCasc)
{
	double expositionTimeInHours = 0;
	TObjArray* fileArray = recCasc.GetListOfFiles();
	for (int i = 0; i < fileArray->GetEntries(); ++i)
	{
		TString fileName = fileArray->At(i)->GetTitle();
		fileName.ReplaceAll("reco.cascade.root","runinfo.root").ReplaceAll("/k","/runinfo/k");
		// cout << fileName.ReplaceAll("reco.cascade.root","runinfo.root").ReplaceAll("/k","/runinfo/k") << endl;
		TFile f{fileName};

		BRunInfo* runInfo = (BRunInfo*)f.Get("RunInfo");
		// cout << runInfo->GetStartTime() << endl;
		// cout << runInfo->GetStopTime() << endl;
		// cout << (runInfo->GetStopTime() - runInfo->GetStartTime())/3600.0<< endl;
		expositionTimeInHours += (runInfo->GetStopTime() - runInfo->GetStartTime())/3600.0;
	//    TTree *tree = f.Get<TTree>(chainElement->GetName());
	}
	cout << "Exposition time: " << expositionTimeInHours << "\t" << expositionTimeInHours/24.0 <<endl;
	cout << "Number of runs: " << fileArray->GetEntries() <<endl;
	return expositionTimeInHours;
}

int SetHistograms(std::map<std::string,TH1D*> &histograms, Bool_t isMCFile)
{
	TH1D* h_energyRec = new TH1D("h_energyRec","Reconstructed energy; E [TeV]; NoE [#]",20000,0,20000);
	histograms.insert(std::make_pair("h_energyRec",h_energyRec));

	TH1D* h_thetaRec = new TH1D("h_thetaRec","Reconstructed theta; #theta [deg.]; NoE [#]",180,0,180);
	histograms.insert(std::make_pair("h_thetaRec",h_thetaRec));

	TH1D* h_cosThetaRec = new TH1D("h_cosThetaRec","Reconstructed cos(theta); cos(#theta) [1]; NoE [#]",20,-1,1);
	histograms.insert(std::make_pair("h_cosThetaRec",h_cosThetaRec));

	TH1D* h_nHits = new TH1D("h_nHits","Number of reco hits; N_{recoHits} [#]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_nHits",h_nHits));

	if (isMCFile)
	{
		TH1D* h_mismatchEnergy = new TH1D("h_mismatchEnergy","Mismatch energy; E_{rec}/E_{MC}; NoE [#]",100,0,10);
		histograms.insert(std::make_pair("h_mismatchEnergy",h_mismatchEnergy));

		TH1D* h_mismatchAngle = new TH1D("h_mismatchAngle","Mismatch angle;Mismatch angle [deg.];NoE [#]",180,0,180);
		histograms.insert(std::make_pair("h_mismatchAngle",h_mismatchAngle));

		TH1D* h_mismatchPosition = new TH1D("h_mismatchPosition","Mismatch position;Mismatch position [m];NoE [#]",200,0,200);
		histograms.insert(std::make_pair("h_mismatchPosition",h_mismatchPosition));

		TH1D* h_causHitEfficiency = new TH1D("h_causHitEfficiency","Causality Hit efficiency; Eff [%]; NoE [%]",110,0,110);
		histograms.insert(std::make_pair("h_causHitEfficiency",h_causHitEfficiency));

		TH1D* h_causHitPurity = new TH1D("h_causHitPurity","Causality Hit purity; Purity [%]; NoE [%]",110,0,110);
		histograms.insert(std::make_pair("h_causHitPurity",h_causHitPurity));

		TH1D* h_hitEfficiency = new TH1D("h_hitEfficiency","Hit efficiency; Eff [%]; NoE [%]",110,0,110);
		histograms.insert(std::make_pair("h_hitEfficiency",h_hitEfficiency));

		TH1D* h_hitPurity = new TH1D("h_hitPurity","Hit purity; Purity [%]; NoE [%]",110,0,110);
		histograms.insert(std::make_pair("h_hitPurity",h_hitPurity));

		TH1D* h_noiseQ = new TH1D("h_noiseQ","Charge of noise hits; Q [p.e.]; NoE [%]",100,0,100);
		histograms.insert(std::make_pair("h_noiseQ",h_noiseQ));
	}

	return histograms.size();
}

int DrawHistograms(std::map<std::string,TH1D*> &histograms, Bool_t isMCFile)
{
	TCanvas* c_energyRec = new TCanvas("c_energyRec","EnergyRec",800,600);
	histograms["h_energyRec"]->Draw("HIST");

	TCanvas* c_thetaRec = new TCanvas("c_thetaRec","ThetaRec",800,600);
	histograms["h_thetaRec"]->Draw("HIST");

	TCanvas* c_cosThetaRec = new TCanvas("c_cosThetaRec","CosThetaRec",800,600);
	histograms["h_cosThetaRec"]->Draw("HIST");

	TCanvas* c_nHits = new TCanvas("c_nHits","NRecoHits",800,600);
	histograms["h_nHits"]->Draw("HIST");

	if (isMCFile)
	{
		TCanvas* c_mismatchEnergy = new TCanvas("c_mismatchEnergy","MismatchEnergy",800,600);
		histograms["h_mismatchEnergy"]->Draw("HIST");

		TCanvas* c_mismatchAngle = new TCanvas("c_mismatchAngle","MismatchAngle",800,600);
		histograms["h_mismatchAngle"]->Draw("HIST");

		TCanvas* c_mismatchPosition = new TCanvas("c_mismatchPosition","MismatchPosition",800,600);
		histograms["h_mismatchPosition"]->Draw("HIST");

		TCanvas* c_causHitEfficiency = new TCanvas("c_causHitEfficiency","CausalityHitEfficiency",800,600);
		histograms["h_causHitEfficiency"]->Draw("HIST");

		TCanvas* c_causHitPurity = new TCanvas("c_causHitPurity","CausalityHitPurity",800,600);
		histograms["h_causHitPurity"]->Draw("HIST");


		TCanvas* c_hitEfficiency = new TCanvas("c_hitEfficiency","HitEfficiency",800,600);
		histograms["h_hitEfficiency"]->Draw("HIST");

		TCanvas* c_hitPurity = new TCanvas("c_hitPurity","HitPurity",800,600);
		histograms["h_hitPurity"]->Draw("HIST");

		TCanvas* c_noiseQ = new TCanvas("c_noiseQ","NoiseQ",800,600);
		histograms["h_noiseQ"]->Draw("HIST");
	}

	return 0;
}

int SaveHistograms(std::map<std::string,TH1D*> &histograms)
{
	return 0;
}

void PrintCascade(BRecoCascade* cascade, BJointHeader* header)
{
	cout << std::setprecision(4) << header->GetSeason() << "\t" << header->GetCluster() << "\t" << header->GetRun() << "\t" << std::setw(10) << header->GetEventIDCC() << "\t" << cascade->GetEnergyRec() << "\t" << cascade->GetThetaRec() << "\t" << cascade->GetNHitsTFil() << "\t" << cascade->GetDistanceCS() << "\t" << cascade->GetLikelihoodHitOnly() << "\t" << cascade->GetQTotal() << endl;
}

void PrintHECascades(BRecoCascade* cascade, BJointHeader* header, double energyThreshold)
{
	// cout << string(40,'*') << " HE cascades" << string(40,'*') << endl;
	if (cascade->GetEnergyRec() > energyThreshold)
	{
		cout << "HE" << "\t";
		PrintCascade(cascade,header);
	}
}

void PrintUpGoingCascades(BRecoCascade* cascade, BJointHeader* header, double thetaThreshold)
{
	// cout << string(40,'*') << " HE cascades" << string(40,'*') << endl;
	if (cascade->GetThetaRec() < thetaThreshold)
	{
		cout << "Up" << "\t";
		PrintCascade(cascade,header);
	}
}

int analyzeWout(int year = -1, int cluster = -1, bool calcExpTime = false, TString fileName = "")
{
	TChain MCEvents("Events");

	if (fileName == "")
		SetTChain(MCEvents,year,cluster);
	else
		SetTChain(MCEvents,fileName);

	BEvent*  	  myEvent 			= NULL;
	BMCEvent*  	  myMCEvent 		= NULL;
	BEventMaskMC* myEventMaskMC 	= NULL;

	MCEvents.SetBranchAddress("BMCEvent.", &myMCEvent);
	MCEvents.SetBranchAddress("MCEventMask.", &myEventMaskMC);
	MCEvents.SetBranchAddress("BEvent.", &myEvent);

	TH1D* h_nHitsPrimaryCascade = new TH1D("h_nHitsPrimaryCascade","Number of hits from the primary cascade; N_{hits} [#]; NoE [#]",100,0,100);

	Int_t nHitsPrimaryCascade = 0;
	Int_t nHitsMuon = 0;
	for (int i = 0; i < MCEvents.GetEntries(); ++i)
	{
		MCEvents.GetEntry(i);
		nHitsPrimaryCascade = 0;
		nHitsMuon = 0;
		for (int i = 0; i < myEvent->NHits(); ++i)
		{
			for (int j = 0; j < myEventMaskMC->GetFlagN(i); ++j)
			{
				if (myEventMaskMC->GetFlag(i,j) == -999)
					nHitsMuon++;
				if (myEventMaskMC->GetFlag(i,j) == 1001)
				{
					nHitsPrimaryCascade++;
					break;
				}
				// cout << myEventMaskMC->GetFlag(i,j) << "\t";
			}
			// cout << endl;
		}
		if (nHitsMuon == 0)
			cout << nHitsPrimaryCascade << "\t" << nHitsMuon << endl;
		h_nHitsPrimaryCascade->Fill(nHitsPrimaryCascade);
		// cout << myMCEvent->GetChannelN() <<endl;
		// cout << myEvent->NHits() << endl;
		// cout << myEventMaskMC->GetEntries() << "\t" << myEventMaskMC->GetNSignalHits() << "\t" << myEventMaskMC->GetNNoiseHits() << endl;
	}

	h_nHitsPrimaryCascade->Draw();

	// DrawHistograms(v_histograms,isMCFile);
	// SaveHistograms(v_histograms);

	return 0;
}