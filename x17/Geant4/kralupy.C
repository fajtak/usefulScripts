#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"

#include <iostream>
#include <map>


int SetHistograms(std::map<TString,TH1D*> &histograms)
{
	TH1D* h_nHits = new TH1D("h_nHits","Number of hits per event;N_{hits} [#];NoE [#]",500,0,500);
	histograms.insert(std::make_pair("h_nHits",h_nHits));

	TH1D* h_eDep = new TH1D("h_eDep","Deposited energy per event;E_{dep} [keV];NoE [#]",110,0,1100);
	histograms.insert(std::make_pair("h_eDep",h_eDep));

	TH1D* h_x = new TH1D("h_x","X position;x [mm];NoE [#]",2000,-10,10);
	histograms.insert(std::make_pair("h_x",h_x));

	TH1D* h_y = new TH1D("h_y","Y position;y [mm];NoE [#]",2000,-10,10);
	histograms.insert(std::make_pair("h_y",h_y));

	TH1D* h_z = new TH1D("h_z","Z position;z [mm];NoE [#]",2000,-10,10);
	histograms.insert(std::make_pair("h_z",h_z));

	TH1D* h_stepLength = new TH1D("h_stepLength","Step length;L_{step} [mm];NoE [#]",200,0,0.2);
	histograms.insert(std::make_pair("h_stepLength",h_stepLength));

	TH1D* h_trackLength = new TH1D("h_trackLength","Track length;L_{track} [mm];NoE [#]",3000,0,30);
	histograms.insert(std::make_pair("h_trackLength",h_trackLength));

	return 0;
}

int DrawHistograms(std::map<TString,TH1D*> histograms)
{
	auto it{ histograms.cbegin() }; // declare a const iterator and assign to start of vector
    while (it != histograms.cend()) // while it hasn't reach the end
    {
        TCanvas* c_tempCanvas = new TCanvas(Form("c_%s",it->first.Data()),it->first.Data(),800,600);
        it->second->Draw();
        ++it; // and iterate to the next element
    }

	return 0;
}

int SaveHistograms(std::map<TString,TH1D*> &histograms, TString fileName)
{

	TString outputFileName = "";
	outputFileName = fileName(0,fileName.Last('/')) + "/results_" + fileName(fileName.Last('/')+1,fileName.Length());

	TFile* f_results = new TFile(outputFileName,"RECREATE");

	auto it{ histograms.cbegin() }; // declare a const iterator and assign to start of vector
	while (it != histograms.cend()) // while it hasn't reach the end
	{
		it->second->Write();
		++it; // and iterate to the next element
	}
	return 0;
}

int kralupy()
{
	// TString fileName = "/Data/x17/Geant4/test_Pos1_NoStep_1M_Elec1MeV_.06mmDiverg.root";
	TString fileName = "/Data/x17/Geant4/test_Pos1_Step10micron_1M_Elec1MeV_.06mmDiverg.root";
	TFile* inputFile = new TFile(fileName);
	TTree* tree = (TTree*)inputFile->Get("trackX17");

	Int_t eventID, totalHits;
	Double_t totalEdep,x,y,z,stepLength,trackLength;
	Int_t previousEventID = -1;
	Int_t nHits = 0;

	tree->SetBranchAddress("Event",&eventID);
	tree->SetBranchAddress("TotalHits",&totalHits);
	tree->SetBranchAddress("TotalEdep",&totalEdep);
	tree->SetBranchAddress("Xpos",&x);
	tree->SetBranchAddress("Ypos",&y);
	tree->SetBranchAddress("Zpos",&z);
	tree->SetBranchAddress("StepLength",&stepLength);
	tree->SetBranchAddress("TrackLength",&trackLength);

	std::map<TString,TH1D*> histograms;
	SetHistograms(histograms);

	for (int i = 0; i < tree->GetEntries(); ++i)
	{
		tree->GetEntry(i);
		if (totalHits != 0)
		{
			histograms["h_nHits"]->Fill(totalHits);
			histograms["h_eDep"]->Fill(totalEdep);
			histograms["h_trackLength"]->Fill(trackLength);
		}
		histograms["h_x"]->Fill(x);
		histograms["h_y"]->Fill(y);
		histograms["h_z"]->Fill(z);
		histograms["h_stepLength"]->Fill(stepLength);
	}

	DrawHistograms(histograms);
	SaveHistograms(histograms,fileName);

	return 0;
}