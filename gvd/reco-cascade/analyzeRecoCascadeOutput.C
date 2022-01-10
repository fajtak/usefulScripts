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
#include "TH2D.h"
#include "TCanvas.h"
#include "TProfile.h"

#include "BRecoCascade.h"
#include "BEvent.h"
#include "BJointHeader.h"
#include "BMCCascadeEvent.h"
#include "BMCEvent.h"
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
			// TString filesDir = Form("/Data/BaikalData/reco-cascade/20%d/cluster%d/*.reco.cascade.root",j,i);
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
		fileName.ReplaceAll("reco.cascade.root","runinfo.root").ReplaceAll("/k","/runinfo/k").ReplaceAll("/i","/runinfo/i");
		// cout << fileName.ReplaceAll("reco.cascade.root","runinfo.root").ReplaceAll("/k","/runinfo/k") << endl;
		TFile f{fileName};
		if (!f.IsOpen())
		{
			cout << "Missing runInfo file for run: " << fileName << endl;
			return -1;
		}

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

int SetHistograms(std::map<std::string,TH1D*> &histograms,std::map<std::string,TH2D*> &histograms2D, Bool_t isMCFile)
{
	TH1D* h_energyRec = new TH1D("h_energyRec","Reconstructed energy; E [TeV]; NoE [#]",5000,0,10000);
	histograms.insert(std::make_pair("h_energyRec",h_energyRec));

	TH1D* h_energyRecError = new TH1D("h_energyRecError","Reconstructed energy error; #sigma_{E} [TeV]; NoE [#]",5000,0,10000);
	histograms.insert(std::make_pair("h_energyRecError",h_energyRecError));

	TH1D* h_thetaRec = new TH1D("h_thetaRec","Reconstructed theta; #theta [deg.]; NoE [#]",36,0,180);
	histograms.insert(std::make_pair("h_thetaRec",h_thetaRec));

	TH1D* h_thetaRecError = new TH1D("h_thetaRecError","Reconstructed theta error; #sigma_{#theta} [deg.]; NoE [#]",36,0,180);
	histograms.insert(std::make_pair("h_thetaRecError",h_thetaRecError));

	TH1D* h_phiRecError = new TH1D("h_phiRecError","Reconstructed phi error; #sigma_{#phi} [deg.]; NoE [#]",72,0,360);
	histograms.insert(std::make_pair("h_phiRecError",h_phiRecError));

	TH1D* h_cosThetaRec = new TH1D("h_cosThetaRec","Reconstructed cos(theta); cos(#theta) [1]; NoE [#]",20,-1,1);
	histograms.insert(std::make_pair("h_cosThetaRec",h_cosThetaRec));

	TH1D* h_nHits = new TH1D("h_nHits","Number of event hits; N_{Hits} [#]; NoE [#]",150,0,150);
	histograms.insert(std::make_pair("h_nHits",h_nHits));

	TH1D* h_nHitsCaus = new TH1D("h_nHitsCaus","Number of event hits after Causality; N_{Hits} [#]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_nHitsCaus",h_nHitsCaus));

	TH1D* h_nHitsReco = new TH1D("h_nHitsReco","Number of reco hits; N_{recoHits} [#]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_nHitsReco",h_nHitsReco));

	TH1D* h_nStringsReco = new TH1D("h_nStringsReco","Number of reco strings; N_{recoStrings} [#]; NoE [#]",10,0,10);
	histograms.insert(std::make_pair("h_nStringsReco",h_nStringsReco));

	TH1D* h_nStringsCaus = new TH1D("h_nStringsCaus","Number of strings after Causality; N_{recoStrings} [#]; NoE [#]",10,0,10);
	histograms.insert(std::make_pair("h_nStringsCaus",h_nStringsCaus));

	TH1D* h_OMID = new TH1D("h_OMID","OM IDs of reco hits; OMID [1]; NoE [#]",300,0,300);
	histograms.insert(std::make_pair("h_OMID",h_OMID));

	TH1D* h_tRes = new TH1D("h_tRes","#delta T of reco hits; #delta T [ns]; NoE [#]",100,-50,50);
	histograms.insert(std::make_pair("h_tRes",h_tRes));

	TH1D* h_logPHit = new TH1D("h_logPHit","Log_{10}(P_{hit}); Log_{10}(P_{hit}) [1]; NoE [#]",100,-10,0);
	histograms.insert(std::make_pair("h_logPHit",h_logPHit));

	TH1D* h_chi2 = new TH1D("h_chi2","#chi^{2} ; #chi^{2} [1]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_chi2",h_chi2));

	TH1D* h_chi2NCalls = new TH1D("h_chi2NCalls","#chi^{2} N_{calls} ; N_{calls} [#]; NoE [#]",200,0,200);
	histograms.insert(std::make_pair("h_chi2NCalls",h_chi2NCalls));

	TH1D* h_chi2Caus = new TH1D("h_chi2Caus","#chi^{2} after Causality ; #chi^{2} [1]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_chi2Caus",h_chi2Caus));

	TH1D* h_like = new TH1D("h_like","Likelihood; Likelihood [1]; NoE [#]",50,0,0.5);
	histograms.insert(std::make_pair("h_like",h_like));

	TH1D* h_likeNCalls = new TH1D("h_likeNCalls","Likelihood N_{calls}; N_{calls} [#]; NoE [#]",100,0,1000);
	histograms.insert(std::make_pair("h_likeNCalls",h_likeNCalls));

	TH1D* h_likeFitStatus = new TH1D("h_likeFitStatus","Likelihood fit status; Fit status [ID]; NoE [#]",6,0,6);
	histograms.insert(std::make_pair("h_likeFitStatus",h_likeFitStatus));

	TH1D* h_likeFitCovMatStatus = new TH1D("h_likeFitCovMatStatus","Likelihood fit covMat status; Fit status [ID]; NoE [#]",6,0,6);
	histograms.insert(std::make_pair("h_likeFitCovMatStatus",h_likeFitCovMatStatus));

	TH1D* h_likeHitOnly = new TH1D("h_likeHitOnly","Likelihood hit only; Likelihood_{HitOnly} [1]; NoE [#]",500,0,5);
	histograms.insert(std::make_pair("h_likeHitOnly",h_likeHitOnly));

	TH1D* h_distanceCS = new TH1D("h_distanceCS","DistanceCS; D_{CS} [m]; NoE [#]",100,0,200);
	histograms.insert(std::make_pair("h_distanceCS",h_distanceCS));

	TH1D* h_z = new TH1D("h_z","Z position; Z [m]; NoE [#]",70,0,700);
	histograms.insert(std::make_pair("h_z",h_z));

	TH1D* h_posChange = new TH1D("h_posChange","Change in the fit position; L [m]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_posChange",h_posChange));

	TH1D* h_posError = new TH1D("h_posError","Error of the fit position; #sigma_{L} [m]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_posError",h_posError));

	TH1D* h_timeError = new TH1D("h_timeError","Error of the fit time; #sigma_{T} [ns]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_timeError",h_timeError));

	TH1D* h_qTotal = new TH1D("h_qTotal","Total charge of the event; Q_{total} [p.e.]; NoE [#]",100,0,10000);
	histograms.insert(std::make_pair("h_qTotal",h_qTotal));

	TH1D* h_nCloseHits = new TH1D("h_nCloseHits","Number of close hits; N_{closeHits} [#]; NoE [#]",20,0,20);
	histograms.insert(std::make_pair("h_nCloseHits",h_nCloseHits));

	TH1D* h_qCloseHits = new TH1D("h_qCloseHits","Charge of close hits; Q_{closeHits} [p.e.]; NoE [#]",100,0,10000);
	histograms.insert(std::make_pair("h_qCloseHits",h_qCloseHits));

	TH1D* h_nTrackHits = new TH1D("h_nTrackHits","Number of track hits; N_{trackHits} [#]; NoE [#]",25,0,25);
	histograms.insert(std::make_pair("h_nTrackHits",h_nTrackHits));

	TH1D* h_qTrack = new TH1D("h_qTrack","Charge of track hits; Q_{trackHits} [p.e.]; NoE [#]",200,0,2000);
	histograms.insert(std::make_pair("h_qTrack",h_qTrack));

	TH1D* h_nTrackHitsSeg = new TH1D("h_nTrackHitsSeg","Number of track hits in segment; N_{trackHits} [#]; NoE [#]",25,0,25);
	histograms.insert(std::make_pair("h_nTrackHitsSeg",h_nTrackHitsSeg));

	TH1D* h_qTrackSeg = new TH1D("h_qTrackSeg","Charge of track hits in segment; Q_{trackHits} [p.e.]; NoE [#]",200,0,2000);
	histograms.insert(std::make_pair("h_qTrackSeg",h_qTrackSeg));

	TH1D* h_branchRatio = new TH1D("h_branchRatio","Branch ratio; BR [1]; NoE [#]",100,0,20);
	histograms.insert(std::make_pair("h_branchRatio",h_branchRatio));

	TH1D* h_qRatio = new TH1D("h_qRatio","Q ratio; QR [1]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_qRatio",h_qRatio));

	TH1D* h_qRatioRed = new TH1D("h_qRatioRed","Reduced Q ratio; QR_{red} [1]; NoE [#]",200,0,200);
	histograms.insert(std::make_pair("h_qRatioRed",h_qRatioRed));

	TH1D* h_qEarly = new TH1D("h_qEarly","Charge of early hits; Q_{early} [p.e.]; NoE [#]",200,0,2000);
	histograms.insert(std::make_pair("h_qEarly",h_qEarly));

	TH1D* h_qRecoHits = new TH1D("h_qRecoHits","Charge of reco hits; Q_{recoHits} [p.e.]; NoE [#]",100,0,10000);
	histograms.insert(std::make_pair("h_qRecoHits",h_qRecoHits));

	if (isMCFile)
	{
		TH1D* h_energyMC = new TH1D("h_energyMC","MC energy; E_{MC} [TeV]; NoE [#]",2000,0,20000);
		histograms.insert(std::make_pair("h_energyMC",h_energyMC));

		TH1D* h_energyMCNW = new TH1D("h_energyMCNW","MC energy; E_{MC} [TeV]; NoE [#]",2000,0,20000);
		histograms.insert(std::make_pair("h_energyMCNW",h_energyMCNW));

		TH1D* h_logEnergyMC = new TH1D("h_logEnergyMC","MC Log(energy); log_{10}(E_{MC} [TeV]) ; NoE [#]",60,0,6);
		histograms.insert(std::make_pair("h_logEnergyMC",h_logEnergyMC));

		TH1D* h_logEnergyMCNW = new TH1D("h_logEnergyMCNW","MC Log(energy); log_{10}(E_{MC} [TeV]) ; NoE [#]",60,0,6);
		histograms.insert(std::make_pair("h_logEnergyMCNW",h_logEnergyMCNW));

		TH1D* h_thetaMC = new TH1D("h_thetaMC","MC theta; #theta_{MC} [deg.]; NoE [#]",20,0,200);
		histograms.insert(std::make_pair("h_thetaMC",h_thetaMC));

		TH1D* h_thetaMCNW = new TH1D("h_thetaMCNW","MC theta; #theta_{MC} [deg.]; NoE [#]",20,0,200);
		histograms.insert(std::make_pair("h_thetaMCNW",h_thetaMCNW));

		TH1D* h_cosThetaMC = new TH1D("h_cosThetaMC","MC cos(#theta); cos(#theta_{MC}) [1]; NoE [#]",20,-1,1);
		histograms.insert(std::make_pair("h_cosThetaMC",h_cosThetaMC));

		TH1D* h_cosThetaMCNW = new TH1D("h_cosThetaMCNW","MC cos(#theta); cos(#theta_{MC}) [1]; NoE [#]",20,-1,1);
		histograms.insert(std::make_pair("h_cosThetaMCNW",h_cosThetaMCNW));

		TH1D* h_distanceCSMC = new TH1D("h_distanceCSMC","MC distance to center string; Distance_{CS} [m]; NoE [#]",20,0,200);
		histograms.insert(std::make_pair("h_distanceCSMC",h_distanceCSMC));

		TH1D* h_distanceCSMCNW = new TH1D("h_distanceCSMCNW","MC distance to center string; Distance_{CS} [m]; NoE [#]",20,0,200);
		histograms.insert(std::make_pair("h_distanceCSMCNW",h_distanceCSMCNW));

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

		TH2D* h_hitPurityVsE = new TH2D("h_hitPurityVsE","Hit purity vs. energy; log_{10}(E) [TeV]; Purity [%]; NoE [#]",60,0,6,110,0,110);
		histograms2D.insert(std::make_pair("h_hitPurityVsE",h_hitPurityVsE));

		TH2D* h_mismatchAngleVsE = new TH2D("h_mismatchAngleVsE","Mismatch angle vs. energy; log_{10}(E) [TeV]; Mismatch angle [deg.]; NoE [#]",60,0,6,180,0,180);
		histograms2D.insert(std::make_pair("h_mismatchAngleVsE",h_mismatchAngleVsE));

		TH2D* h_likeVsE = new TH2D("h_likeVsE","Likelihood vs. energy; log_{10}(E) [TeV]; Likelihood [1]; NoE [#]",60,0,6,50,0,5);
		histograms2D.insert(std::make_pair("h_likeVsE",h_likeVsE));

		TH2D* h_likeHitOnlyVsE = new TH2D("h_likeHitOnlyVsE","LikelihoodHitOnly vs. energy; log_{10}(E) [TeV]; Likelihood_{HitOnly} [1]; NoE [#]",60,0,6,200,0,20);
		histograms2D.insert(std::make_pair("h_likeHitOnlyVsE",h_likeHitOnlyVsE));
	}

	return histograms.size();
}

int DrawHistograms(std::map<std::string,TH1D*> &histograms,std::map<std::string,TH2D*> &histograms2D, Bool_t isMCFile)
{
	TCanvas* c_energyRec = new TCanvas("c_energyRec","EnergyRec",800,600);
	histograms["h_energyRec"]->Draw("HIST");

	TCanvas* c_thetaRec = new TCanvas("c_thetaRec","ThetaRec",800,600);
	histograms["h_thetaRec"]->Draw("HIST");

	TCanvas* c_cosThetaRec = new TCanvas("c_cosThetaRec","CosThetaRec",800,600);
	histograms["h_cosThetaRec"]->Draw("HIST");

	TCanvas* c_nHits = new TCanvas("c_nHits","NRecoHits",800,600);
	histograms["h_nHits"]->Draw("HIST");

	TCanvas* c_nHitsReco = new TCanvas("c_nHitsReco","NRecoHits",800,600);
	histograms["h_nHitsReco"]->Draw("HIST");

	TCanvas* c_nStringsReco = new TCanvas("c_nStringsReco","NRecoStrings",800,600);
	histograms["h_nStringsReco"]->Draw("HIST");

	TCanvas* c_OMID = new TCanvas("c_OMID","OMID",800,600);
	histograms["h_OMID"]->Draw("HIST");

	TCanvas* c_tRes = new TCanvas("c_tRes","TResiduals",800,600);
	histograms["h_tRes"]->Draw("HIST");

	TCanvas* c_chi2 = new TCanvas("c_chi2","Chi2",800,600);
	histograms["h_chi2"]->Draw("HIST");

	TCanvas* c_like = new TCanvas("c_like","Likelihood",800,600);
	histograms["h_like"]->Draw("HIST");

	TCanvas* c_likeHitOnly = new TCanvas("c_likeHitOnly","LikelihoodHitOnly",800,600);
	histograms["h_likeHitOnly"]->Draw("HIST");

	TCanvas* c_distanceCS = new TCanvas("c_distanceCS","DistanceCS",800,600);
	histograms["h_distanceCS"]->Draw("HIST");


	if (isMCFile)
	{
		TCanvas* c_energyMC = new TCanvas("c_energyMC","EnergyMC",800,600);
		histograms["h_energyMC"]->Draw("HIST");

		TCanvas* c_cosThetaMC = new TCanvas("c_cosThetaMC","CosThetaMC",800,600);
		histograms["h_cosThetaMC"]->Draw("HIST");

		TCanvas* c_cosThetaMCNW = new TCanvas("c_cosThetaMCNW","CosThetaMCNW",800,600);
		histograms["h_cosThetaMCNW"]->Draw("HIST");

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

		TCanvas* c_hitPurityVsE = new TCanvas("c_hitPurityVsE","HitPurityVsE",800,600);
		histograms2D["h_hitPurityVsE"]->ProfileX()->Draw("COLZ");

		TCanvas* c_mismatchAngleVsE = new TCanvas("c_mismatchAngleVsE","MismatchAngleVsE",800,600);
		histograms2D["h_mismatchAngleVsE"]->Draw("COLZ");
		// histograms2D["h_mismatchAngleVsE"]->ProfileX()->Draw("SAME");
		histograms2D["h_mismatchAngleVsE"]->QuantilesX()->Draw("SAME");

		TCanvas* c_likeVsE = new TCanvas("c_likeVsE","LikelihoodVsE",800,600);
		histograms2D["h_likeVsE"]->Draw("COLZ");
		histograms2D["h_likeVsE"]->QuantilesX()->Draw("SAME");

		TCanvas* c_likeHitOnlyVsE = new TCanvas("c_likeHitOnlyVsE","LikelihoodHitOnlyVsE",800,600);
		histograms2D["h_likeHitOnlyVsE"]->Draw("COLZ");
		histograms2D["h_likeHitOnlyVsE"]->QuantilesX()->Draw("SAME");
	}

	return 0;
}

int SaveHistograms(std::map<std::string,TH1D*> &histograms,TString fileName,int season, int cluster)
{

	TString outputFileName = "";
	if (fileName != "")
		outputFileName = fileName(0,fileName.Last('/')) + "/results_" + fileName(fileName.Last('/')+1,fileName.Length());
	else
		outputFileName = Form("/Data/BaikalData/reco-cascade/20%d/cluster%d/reco.cascade.results.root",season,cluster);
	TFile* f_results = new TFile(outputFileName,"RECREATE");

	auto it{ histograms.cbegin() }; // declare a const iterator and assign to start of vector
	while (it != histograms.cend()) // while it hasn't reach the end
	{
		it->second->Write();
		++it; // and iterate to the next element
	}

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

int analyzeRecoCascadeOutput(int year = -1, int cluster = -1, bool calcExpTime = false, TString fileName = "")
{
	TChain reconstructedCascades("Events");

	if (fileName == "")
		SetTChain(reconstructedCascades,year,cluster);
	else
		SetTChain(reconstructedCascades,fileName);

	BRecoCascade* myCascade = NULL;
	BJointHeader* myHeader  = NULL;
	BEvent*       myEvent   = NULL;

	TBranch* b_joint = reconstructedCascades.FindBranch("BJointHeader.");
	TBranch* b_reco = reconstructedCascades.FindBranch("BRecoCascade.");
	TBranch* b_event = reconstructedCascades.FindBranch("BEvent.");
    if(!b_joint || !b_reco || !b_event) {
        std::cout << "Could not find the BJointHeader or BRecoCascade branch on tree Events, cannot continue." << std::endl;
        return -1;
    }
	reconstructedCascades.SetBranchAddress("BRecoCascade.", &myCascade);
	reconstructedCascades.SetBranchAddress("BJointHeader.", &myHeader);
	reconstructedCascades.SetBranchAddress("BEvent.", &myEvent);

	Bool_t isMCFile = false;
	Bool_t isMCCascadeFile = false;
	BMCCascadeEvent* myMCCascade = NULL;
    TBranch* b_mcCascade = reconstructedCascades.FindBranch("BMCCascadeEvent.");
    if (b_mcCascade)
    {
    	isMCFile = true;
    	isMCCascadeFile = true;
    	reconstructedCascades.SetBranchAddress("BMCCascadeEvent.", &myMCCascade);
    }
    BMCEvent* myMCEvent = NULL;
    TBranch* b_mcEvent = reconstructedCascades.FindBranch("BMCEvent.");
    if (b_mcEvent)
    {
    	isMCFile = true;
    	reconstructedCascades.SetBranchAddress("BMCEvent.", &myMCEvent);
    }

	std::map<std::string, TH1D*> v_histograms;
	std::map<std::string, TH2D*> v_histograms2D;
	// vector<TH1D*> v_histograms;
	SetHistograms(v_histograms,v_histograms2D,isMCFile);
	if (calcExpTime)
		CalculateExpositionTime(reconstructedCascades);
	Float_t eventWeight = 1.0;

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);
		if (isMCFile)
			eventWeight = myCascade->GetEventWeight();
		// if (myCascade->GetFitPos().Z() > 610)
		if (myCascade->GetFitPos().Z() > (isMCFile?220.3:580))
			continue;
		if (myCascade->GetLikelihoodHitOnly() > 1.5)
			continue;
		// if (myCascade->GetEnergyRec() < 10)
			// continue;
		if (myCascade->GetDistanceCS() > 100)
			continue;
		// if (myCascade->GetDistanceCS() > 65)
		// if (myCascade->GetDistanceCS() > 65)
			// continue;
		if (myCascade->GetNHitsTFil() < 20)
			continue;

		// if (myCascade->GetNHitsTrackSeg() > 1)
			// continue;
		// if (myCascade->GetEnergyRec() < 500)
			// continue;
		// if (myCascade->GetNHitsTFil() < 50)
		// 	continue;
		v_histograms["h_energyRec"]->Fill(myCascade->GetEnergyRec(),eventWeight);
		v_histograms["h_energyRecError"]->Fill(myCascade->GetEnergyRecError(),eventWeight);
		v_histograms["h_thetaRec"]->Fill(myCascade->GetThetaRec(),eventWeight);
		v_histograms["h_thetaRecError"]->Fill(myCascade->GetThetaRecError(),eventWeight);
		v_histograms["h_phiRecError"]->Fill(myCascade->GetPhiRecError(),eventWeight);
		v_histograms["h_cosThetaRec"]->Fill(TMath::Cos(myCascade->GetThetaRec()*TMath::DegToRad()),eventWeight);
		v_histograms["h_nHits"]->Fill(myCascade->GetNHits(),eventWeight);
		v_histograms["h_nHitsCaus"]->Fill(myCascade->GetNHitsCaus(),eventWeight);
		v_histograms["h_nHitsReco"]->Fill(myCascade->GetNHitsTFil(),eventWeight);
		v_histograms["h_nStringsCaus"]->Fill(myCascade->GetNStringsCaus(),eventWeight);
		v_histograms["h_nStringsReco"]->Fill(myCascade->GetNStringsTFil(),eventWeight);
		v_histograms["h_chi2"]->Fill(myCascade->GetChi2(),eventWeight);
		v_histograms["h_chi2NCalls"]->Fill(myCascade->GetChi2NCalls(),eventWeight);
		v_histograms["h_chi2Caus"]->Fill(myCascade->GetChi2Caus(),eventWeight);
		v_histograms["h_like"]->Fill(myCascade->GetLikelihood(),eventWeight);
		v_histograms["h_likeNCalls"]->Fill(myCascade->GetLikeNCalls(),eventWeight);
		v_histograms["h_likeFitStatus"]->Fill(myCascade->GetLikeFitStatus(),eventWeight);
		v_histograms["h_likeFitCovMatStatus"]->Fill(myCascade->GetLikeFitCovMatStatus(),eventWeight);
		v_histograms["h_likeHitOnly"]->Fill(myCascade->GetLikelihoodHitOnly(),eventWeight);
		v_histograms["h_distanceCS"]->Fill(myCascade->GetDistanceCS(),eventWeight);
		if (!isMCFile)
			v_histograms["h_z"]->Fill(myCascade->GetFitPos().Z(),eventWeight);
		else
			v_histograms["h_z"]->Fill(myCascade->GetFitPos().Z()+359.7,eventWeight);
		v_histograms["h_posChange"]->Fill((myCascade->GetFitPos()-myCascade->GetInitPos()).Mag(),eventWeight);
		v_histograms["h_posError"]->Fill(myCascade->GetFitPosError().Mag(),eventWeight);
		v_histograms["h_timeError"]->Fill(myCascade->GetFitTimeError(),eventWeight);
		v_histograms["h_qTotal"]->Fill(myCascade->GetQTotal(),eventWeight);
		v_histograms["h_nCloseHits"]->Fill(myCascade->GetNCloseHits(),eventWeight);
		v_histograms["h_qCloseHits"]->Fill(myCascade->GetQCloseHits(),eventWeight);
		v_histograms["h_nTrackHits"]->Fill(myCascade->GetNHitsTrack(),eventWeight);
		v_histograms["h_nTrackHitsSeg"]->Fill(myCascade->GetNHitsTrackSeg(),eventWeight);
		v_histograms["h_qTrack"]->Fill(myCascade->GetQTrack(),eventWeight);
		v_histograms["h_qTrackSeg"]->Fill(myCascade->GetQTrackSeg(),eventWeight);
		v_histograms["h_branchRatio"]->Fill(myCascade->GetBranchRatio(),eventWeight);
		v_histograms["h_qRatio"]->Fill(myCascade->GetQRatio(),eventWeight);
		v_histograms["h_qRatioRed"]->Fill(myCascade->GetQRatioReduced(),eventWeight);
		v_histograms["h_qEarly"]->Fill(myCascade->GetQEarly(),eventWeight);
		v_histograms["h_qRecoHits"]->Fill(myCascade->GetQRecoHits(),eventWeight);

		for (int i = 0; i < myCascade->GetNImpulseIDs(); ++i)
		{
			v_histograms["h_OMID"]->Fill(myEvent->HitChannel(myCascade->GetImpulseID(i)),eventWeight);
			v_histograms["h_tRes"]->Fill(myCascade->GetTResPerHit(i),eventWeight);
			v_histograms["h_logPHit"]->Fill(myCascade->GetLikePerHit(i),eventWeight);
		}
		if (isMCCascadeFile)
		{
			v_histograms["h_energyMC"]->Fill(myMCCascade->GetShowerEnergy(),eventWeight);
			v_histograms["h_energyMCNW"]->Fill(myMCCascade->GetShowerEnergy());
			v_histograms["h_thetaMC"]->Fill(myCascade->GetThetaMC(),eventWeight);
			v_histograms["h_thetaMCNW"]->Fill(myCascade->GetThetaMC());
			v_histograms["h_distanceCSMC"]->Fill(myCascade->GetDistanceCSMC(),eventWeight);
			v_histograms["h_distanceCSMCNW"]->Fill(myCascade->GetDistanceCSMC());
			v_histograms["h_cosThetaMC"]->Fill(TMath::Cos(myCascade->GetThetaMC()*TMath::DegToRad()),eventWeight);
			v_histograms["h_cosThetaMCNW"]->Fill(TMath::Cos(myCascade->GetThetaMC()*TMath::DegToRad()));
			v_histograms["h_logEnergyMC"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()),eventWeight);
			v_histograms["h_logEnergyMCNW"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()));
			v_histograms["h_mismatchEnergy"]->Fill(myCascade->GetEnergyRec()/myMCCascade->GetShowerEnergy(),eventWeight);
			TVector3 trueDir;
			trueDir.SetMagThetaPhi(1,myMCCascade->GetNeutrinoZenith(),myMCCascade->GetNeutrinoAzimuth());
			TVector3 recDir;
			recDir.SetMagThetaPhi(1,myCascade->GetThetaRec()*TMath::DegToRad(),myCascade->GetPhiRec()*TMath::DegToRad());
			v_histograms["h_mismatchAngle"]->Fill(trueDir.Angle(recDir)*TMath::RadToDeg(),eventWeight);
			v_histograms["h_mismatchPosition"]->Fill((myCascade->GetFitPos()-myCascade->GetPosMC()).Mag(),eventWeight);
			v_histograms2D["h_mismatchAngleVsE"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()),trueDir.Angle(recDir)*TMath::RadToDeg(),eventWeight);
			v_histograms2D["h_likeVsE"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()),myCascade->GetLikelihood(),eventWeight);
			v_histograms2D["h_likeHitOnlyVsE"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()),myCascade->GetLikelihoodHitOnly(),eventWeight);
			Int_t nSignalHitsTrue = 0;
			for (int i = 0; i < myMCCascade->GetNHits(); ++i)
			{
				if (myMCCascade->GetPulse(i)->GetMagic() > 0)
					nSignalHitsTrue++;
				if (myMCCascade->GetPulse(i)->GetMagic() == 0)
					v_histograms["h_noiseQ"]->Fill(myMCCascade->GetPulse(i)->GetAmplitude());
			}
			Int_t nSignalHitsReco = 0;
			for (int i = 0; i < myCascade->GetNImpulseIDs(); ++i)
			{
				if (myMCCascade->GetPulse(myEvent->Hit(myCascade->GetImpulseID(i)).GetUniqueID())->GetMagic() > 0)
					nSignalHitsReco++;
			}
			Int_t nSignalCausHitsReco = 0;
			for (int i = 0; i < myCascade->GetNCausImpulseIDs(); ++i)
			{
				if (myMCCascade->GetPulse(myEvent->Hit(myCascade->GetCausImpulseID(i)).GetUniqueID())->GetMagic() > 0)
					nSignalCausHitsReco++;
			}
			v_histograms["h_hitEfficiency"]->Fill((double)nSignalHitsReco/nSignalHitsTrue*100,eventWeight);
			v_histograms["h_hitPurity"]->Fill((double)nSignalHitsReco/myCascade->GetNImpulseIDs()*100,eventWeight);
			v_histograms["h_causHitEfficiency"]->Fill((double)nSignalCausHitsReco/nSignalHitsTrue*100,eventWeight);
			v_histograms["h_causHitPurity"]->Fill((double)nSignalCausHitsReco/myCascade->GetNCausImpulseIDs()*100,eventWeight);
			v_histograms2D["h_hitPurityVsE"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()),(double)nSignalHitsReco/myCascade->GetNImpulseIDs()*100,eventWeight);
		}

		PrintHECascades(myCascade,myHeader,100.0);
		PrintUpGoingCascades(myCascade,myHeader,80.0);
	}

	DrawHistograms(v_histograms,v_histograms2D,isMCFile);
	SaveHistograms(v_histograms,fileName,year,cluster);

	return 0;
}