#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"

#include <iostream>
#include <map>


int SetHistograms(std::map<TString,TH1D*> &histograms)
{
	TH1D* h_clusterSize = new TH1D("h_clusterSize","Cluster size;cluster size [#];NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_clusterSize",h_clusterSize));

	TH1D* h_pixX = new TH1D("h_pixX","Pix X ;X_{pix} [#];NoE [#]",256,0,256);
	histograms.insert(std::make_pair("h_pixX",h_pixX));

	TH1D* h_pixY = new TH1D("h_pixY","Pix Y ;Y_{pix} [#];NoE [#]",256,0,256);
	histograms.insert(std::make_pair("h_pixY",h_pixY));

	TH1D* h_clusterMeanX = new TH1D("h_clusterMeanX","Cluster mean X ;X_{cluster} [#];NoE [#]",256,0,256);
	histograms.insert(std::make_pair("h_clusterMeanX",h_clusterMeanX));

	TH1D* h_clusterMeanY = new TH1D("h_clusterMeanY","Cluster mean Y ;Y_{cluster} [#];NoE [#]",256,0,256);
	histograms.insert(std::make_pair("h_clusterMeanY",h_clusterMeanY));

	TH1D* h_clusterHeightKeV = new TH1D("h_clusterHeightKeV","Cluster height keV ;Height_{cluster} [keV];NoE [#]",600,0,1200);
	histograms.insert(std::make_pair("h_clusterHeightKeV",h_clusterHeightKeV));

	TH1D* h_clusterVolumeKeV = new TH1D("h_clusterVolumeKeV","Cluster volume keV ;Volume_{cluster} [keV];NoE [#]",600,0,1200);
	histograms.insert(std::make_pair("h_clusterVolumeKeV",h_clusterVolumeKeV));

	TH1D* h_ToTKeV = new TH1D("h_ToTKeV","ToT ;ToT [keV];NoE [#]",600,0,1200);
	histograms.insert(std::make_pair("h_ToTKeV",h_ToTKeV));

	TH1D* h_clusterVolumeCentroidX = new TH1D("h_clusterVolumeCentroidX","Cluster volume centroid X;X_{clusterVolCentroid} [#];NoE [#]",256,0,256);
	histograms.insert(std::make_pair("h_clusterVolumeCentroidX",h_clusterVolumeCentroidX));

	TH1D* h_clusterVolumeCentroidY = new TH1D("h_clusterVolumeCentroidY","Cluster volume centroid Y;Y_{clusterVolCentroid} [#];NoE [#]",256,0,256);
	histograms.insert(std::make_pair("h_clusterVolumeCentroidY",h_clusterVolumeCentroidY));

	TH1D* h_deltaToA = new TH1D("h_deltaToA","Delta ToA ;Delta ToA [#];NoE [#]",250,0,500);
	histograms.insert(std::make_pair("h_deltaToA",h_deltaToA));

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

int kralupyClustered()
{
	// TString fileName = "/home/fajtak/work/allpix-squared/output/kralupy_noStep/ClusterFiles/root/data_CSAEmul_noADCres.root";
	// TString fileName = "/home/fajtak/work/allpix-squared/output/kralupy_step10micron/ClusterFiles/root/data_CSAEmul_noADCres.root";
	// TString fileName = "/Data/x17/Kralupy/Run002_G3_1MeV_90Degree_Pos2_0.root";
	// TString fileName = "/Data/x17/Kralupy/Run002_G3_1MeV_90Degree_0.root";
	TString fileName = "/home/fajtak/work/allpix-squared/output/kralupy_pos2_step10micron_1MeV_0.09mm/ClusterFiles/root/data_CSAEmul_noADCres.root";
	// TString fileName = "/home/fajtak/work/allpix-squared/output/kralupy_pos1_noStep_1MeV_0.09mm/ClusterFiles/root/data_CSAEmul_noADCres.root";
	TFile* inputFile = new TFile(fileName);
	TTree* tree = (TTree*)inputFile->Get("clusteredData");

	Int_t clusterSize;
	Short_t pixX[500], pixY[500];
	Float_t ToTKeV[500];
	Float_t clusterMeanX, clusterMeanY, clusterHeightKeV, clusterVolumeKeV, clusterVolumeCentroidX, clusterVolumeCentroidY;
	Double_t deltaToA;

	tree->SetBranchAddress("clstrSize",&clusterSize);
	tree->SetBranchAddress("PixX",pixX);
	tree->SetBranchAddress("PixY",pixY);
	tree->SetBranchAddress("clstrMeanX",&clusterMeanX);
	tree->SetBranchAddress("clstrMeanY",&clusterMeanY);
	tree->SetBranchAddress("clstrHeight_keV",&clusterHeightKeV);
	tree->SetBranchAddress("clstrVolume_keV",&clusterVolumeKeV);
	tree->SetBranchAddress("ToT_keV",ToTKeV);
	tree->SetBranchAddress("clstrVolCentroidX",&clusterVolumeCentroidX);
	tree->SetBranchAddress("clstrVolCentroidY",&clusterVolumeCentroidY);
	tree->SetBranchAddress("delta_ToA",&deltaToA);

	std::map<TString,TH1D*> histograms;
	SetHistograms(histograms);

	for (int i = 0; i < tree->GetEntries(); ++i)
	{
		tree->GetEntry(i);
		histograms["h_clusterSize"]->Fill(clusterSize);
		histograms["h_clusterMeanX"]->Fill(clusterMeanX);
		histograms["h_clusterMeanY"]->Fill(clusterMeanY);
		histograms["h_clusterHeightKeV"]->Fill(clusterHeightKeV);
		histograms["h_clusterVolumeKeV"]->Fill(clusterVolumeKeV);
		histograms["h_clusterVolumeCentroidX"]->Fill(clusterVolumeCentroidX);
		histograms["h_clusterVolumeCentroidY"]->Fill(clusterVolumeCentroidY);
		histograms["h_deltaToA"]->Fill(deltaToA);
		for (int j = 0; j < clusterSize; ++j)
		{
			histograms["h_pixX"]->Fill(pixX[j]);
			histograms["h_pixY"]->Fill(pixY[j]);
			histograms["h_ToTKeV"]->Fill(ToTKeV[j]);
		}
	}

	DrawHistograms(histograms);
	SaveHistograms(histograms,fileName);

	return 0;
}