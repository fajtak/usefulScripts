int showTracks()
{
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

	// std::map<TString,TH1D*> histograms;
	// SetHistograms(histograms);
	TH2D* h_timepix = new TH2D("h_timepix","Timepix;X [pixel];Y [pixel]",256,0,256,256,0,256);

	for (int i = 0; i < tree->GetEntries(); ++i)
	{
		tree->GetEntry(i);
		// histograms["h_clusterSize"]->Fill(clusterSize);
		if (clusterSize < 10 || clusterSize > 15)
			continue;
		for (int j = 0; j < clusterSize; ++j)
		{
			h_timepix->Fill(pixX[j],pixY[j]);
		}
	}

	h_timepix->Draw("colz");
	// DrawHistograms(histograms);
	// SaveHistograms(histograms,fileName);

	return 0;
}