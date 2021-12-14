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

int SetHistograms(std::map<std::string,TH1D*> &histograms,std::map<std::string,TH2D*> &histograms2D, Bool_t isMCFile)
{
	TH1D* h_energyRec = new TH1D("h_energyRec","Reconstructed energy; E [TeV]; NoE [#]",2000,0,20000);
	histograms.insert(std::make_pair("h_energyRec",h_energyRec));

	TH1D* h_thetaRec = new TH1D("h_thetaRec","Reconstructed theta; #theta [deg.]; NoE [#]",180,0,180);
	histograms.insert(std::make_pair("h_thetaRec",h_thetaRec));

	TH1D* h_cosThetaRec = new TH1D("h_cosThetaRec","Reconstructed cos(theta); cos(#theta) [1]; NoE [#]",20,-1,1);
	histograms.insert(std::make_pair("h_cosThetaRec",h_cosThetaRec));

	TH1D* h_nHits = new TH1D("h_nHits","Number of reco hits; N_{recoHits} [#]; NoE [#]",100,0,100);
	histograms.insert(std::make_pair("h_nHits",h_nHits));

	TH1D* h_OMID = new TH1D("h_OMID","OM IDs of reco hits; OMID [1]; NoE [#]",300,0,300);
	histograms.insert(std::make_pair("h_OMID",h_OMID));

	TH1D* h_tRes = new TH1D("h_tRes","#delta T of reco hits; #delta T [ns]; NoE [#]",100,-50,50);
	histograms.insert(std::make_pair("h_tRes",h_tRes));


	if (isMCFile)
	{
		TH1D* h_energyMC = new TH1D("h_energyMC","MC energy; E_{MC} [TeV]; NoE [#]",2000,0,20000);
		histograms.insert(std::make_pair("h_energyMC",h_energyMC));

		TH1D* h_logEnergyMC = new TH1D("h_logEnergyMC","MC Log(energy); log_{10}(E_{MC} [TeV]) ; NoE [#]",60,0,6);
		histograms.insert(std::make_pair("h_logEnergyMC",h_logEnergyMC));

		TH1D* h_thetaMC = new TH1D("h_thetaMC","MC theta; #theta_{MC} [deg.]; NoE [#]",20,0,200);
		histograms.insert(std::make_pair("h_thetaMC",h_thetaMC));

		TH1D* h_cosThetaMC = new TH1D("h_cosThetaMC","MC cos(#theta); cos(#theta_{MC}) [1]; NoE [#]",20,-1,1);
		histograms.insert(std::make_pair("h_cosThetaMC",h_cosThetaMC));

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

	TCanvas* c_OMID = new TCanvas("c_OMID","OMID",800,600);
	histograms["h_OMID"]->Draw("HIST");

	TCanvas* c_tRes = new TCanvas("c_tRes","TResiduals",800,600);
	histograms["h_tRes"]->Draw("HIST");


	if (isMCFile)
	{
		TCanvas* c_energyMC = new TCanvas("c_energyMC","EnergyMC",800,600);
		histograms["h_energyMC"]->Draw("HIST");

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

	}

	return 0;
}

int SaveHistograms(std::map<std::string,TH1D*> &histograms,TString fileName)
{
	TString outputFileName = fileName(0,fileName.Last('/')) + "/results_" + fileName(fileName.Last('/')+1,fileName.Length());
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
	BMCCascadeEvent* myMCCascade = NULL;
    TBranch* b_mcCascade = reconstructedCascades.FindBranch("BMCCascadeEvent.");
    if (b_mcCascade)
    {
    	isMCFile = true;
    	reconstructedCascades.SetBranchAddress("BMCCascadeEvent.", &myMCCascade);
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
			eventWeight = myMCCascade->GetEventWeight();
		// if (myCascade->GetFitPos().Z() > 610)
		if (myCascade->GetFitPos().Z() > 580)
			continue;
		if (myCascade->GetLikelihoodHitOnly() > 1.5)
			continue;
		// if (myCascade->GetEnergyRec() < 10)
			// continue;
		// if (myCascade->GetDistanceCS() > 100)
		// if (myCascade->GetDistanceCS() > 65)
		// if (myCascade->GetDistanceCS() < 100)
			// continue;
		// if (myCascade->GetEnergyRec() < 500)
			// continue;
		// if (myCascade->GetNHitsTFil() < 50)
		// 	continue;
		v_histograms["h_energyRec"]->Fill(myCascade->GetEnergyRec(),eventWeight);
		v_histograms["h_thetaRec"]->Fill(myCascade->GetThetaRec(),eventWeight);
		v_histograms["h_cosThetaRec"]->Fill(TMath::Cos(myCascade->GetThetaRec()*TMath::RadToDeg()),eventWeight);
		v_histograms["h_nHits"]->Fill(myCascade->GetNHitsTFil(),eventWeight);

		for (int i = 0; i < myCascade->GetNImpulseIDs(); ++i)
		{
			v_histograms["h_OMID"]->Fill(myEvent->HitChannel(myCascade->GetImpulseID(i)),eventWeight);
			v_histograms["h_tRes"]->Fill(myCascade->GetTResPerHit(i),eventWeight);
		}
		if (isMCFile)
		{
			v_histograms["h_energyMC"]->Fill(myMCCascade->GetShowerEnergy(),eventWeight);
			v_histograms["h_thetaMC"]->Fill(myCascade->GetThetaMC(),eventWeight);
			v_histograms["h_cosThetaMC"]->Fill(TMath::Cos(myCascade->GetThetaMC()*TMath::DegToRad()),eventWeight);
			v_histograms["h_logEnergyMC"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()),eventWeight);
			v_histograms["h_mismatchEnergy"]->Fill(myCascade->GetEnergyRec()/myMCCascade->GetShowerEnergy(),eventWeight);
			TVector3 trueDir;
			trueDir.SetMagThetaPhi(1,myMCCascade->GetNeutrinoZenith(),myMCCascade->GetNeutrinoAzimuth());
			TVector3 recDir;
			recDir.SetMagThetaPhi(1,myCascade->GetThetaRec()*TMath::DegToRad(),myCascade->GetPhiRec()*TMath::DegToRad());
			v_histograms["h_mismatchAngle"]->Fill(trueDir.Angle(recDir)*TMath::RadToDeg(),eventWeight);
			v_histograms["h_mismatchPosition"]->Fill((myCascade->GetFitPos()-myCascade->GetPosMC()).Mag(),eventWeight);
			v_histograms2D["h_mismatchAngleVsE"]->Fill(TMath::Log10(myMCCascade->GetShowerEnergy()),trueDir.Angle(recDir)*TMath::RadToDeg(),eventWeight);
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
		// PrintUpGoingCascades(myCascade,myHeader,80.0);
	}

	DrawHistograms(v_histograms,v_histograms2D,isMCFile);
	SaveHistograms(v_histograms,fileName);

	return 0;
}