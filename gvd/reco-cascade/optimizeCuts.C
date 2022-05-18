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
#include "BGeomTel.h"
#include "BPhysNamespace.h"

double SetTChain(TChain &recCasc, TString fileName)
{
	recCasc.Add(fileName);
	return recCasc.GetEntries();
}

int SetHistograms(std::map<TString,TH1D*> &histograms,std::map<TString,TH2D*> &histograms2D)
{
	TH1D* h_chargePerHit = new TH1D("h_chargePerHit","Charge per hit; Q/Hit [p.e./#]; NoE [#]",1000,0,1000);
	histograms.insert(std::make_pair("h_chargePerHit",h_chargePerHit));

	TH1D* h_chi2 = new TH1D("h_chi2","Chi^2; #chi^{2} [1]; NoE [#]",1000,0,1000);
	histograms.insert(std::make_pair("h_chi2",h_chi2));

	TH1D* h_initMismatchPosition = new TH1D("h_initMismatchPosition","Mismatch position; Initial mismatch position [m]; NoE [#]",100,0,1000);
	histograms.insert(std::make_pair("h_initMismatchPosition",h_initMismatchPosition));

	TH1D* h_fitMismatchPosition = new TH1D("h_fitMismatchPosition","Mismatch position; Initial mismatch position [m]; NoE [#]",100,0,1000);
	histograms.insert(std::make_pair("h_fitMismatchPosition",h_fitMismatchPosition));

	TH1D* h_causHitPurity = new TH1D("h_causHitPurity","Causality Hit purity; Hit selection purity [%]; NoE*eventWeight [1]",110,0,110);
	histograms.insert(std::make_pair("h_causHitPurity",h_causHitPurity));

	TH1D* h_causHitEfficiency = new TH1D("h_causHitEfficiency","Causality Hit efficiency; Hit selection efficiency [%]; NoE*eventWeight [1]",110,0,110);
	histograms.insert(std::make_pair("h_causHitEfficiency",h_causHitEfficiency));

	TH1D* h_hitEfficiency = new TH1D("h_hitEfficiency","Hit efficiency; Hit selection efficiency [%]; NoE*eventWeight [1]",110,0,110);
	histograms.insert(std::make_pair("h_hitEfficiency",h_hitEfficiency));

	TH1D* h_nHitsCaus = new TH1D("h_nHitsCaus","Number of hits after causality; N_{hits}^{caus} [#]; NoE*eventWeight [1]",110,0,110);
	histograms.insert(std::make_pair("h_nHitsCaus",h_nHitsCaus));

	TH1D* h_nHitsFit = new TH1D("h_nHitsFit","Number of hits after fit; N_{hits}^{fit} [#]; NoE*eventWeight [1]",110,0,110);
	histograms.insert(std::make_pair("h_nHitsFit",h_nHitsFit));

	TH1D* h_nHitsDiff = new TH1D("h_nHitsDiff","Number of hits difference; #Delta N_{hits} [#]; NoE*eventWeight [1]",220,-110,110);
	histograms.insert(std::make_pair("h_nHitsDiff",h_nHitsDiff));

	TH1D* h_zDistanceSignalHits = new TH1D("h_zDistanceSignalHits","Distance of signal hits; #delta Z [m]; NoE [#]",1000,0,1000);
	histograms.insert(std::make_pair("h_zDistanceSignalHits",h_zDistanceSignalHits));

	TH1D* h_zDistanceNoiseHits = new TH1D("h_zDistanceNoiseHits","Distance of noise hits; #delta Z [m]; NoE [#]",1000,0,1000);
	histograms.insert(std::make_pair("h_zDistanceNoiseHits",h_zDistanceNoiseHits));

	TH2D* h_initMisPosVsChi2 = new TH2D("h_initMisPosVsChi2","Initial mismatch position vs chi2; Initial mismatch position [m]; #chi^{2} [1]; NoE [#]",100,0,1000,100,0,1000);
	histograms2D.insert(std::make_pair("h_initMisPosVsChi2",h_initMisPosVsChi2));

	TH2D* h_initMisPosVsInitPrec = new TH2D("h_initMisPosVsInitPrec","Initial mismatch position vs init precision; Initial mismatch position [m]; Init precision [1]; NoE [#]",100,0,1000,1000,0,10000);
	histograms2D.insert(std::make_pair("h_initMisPosVsInitPrec",h_initMisPosVsInitPrec));

	TH2D* h_initMisPosVsQTotal = new TH2D("h_initMisPosVsQTotal","Initial mismatch position vs totalQ; Initial mismatch position [m]; Q_{total} [p.e.]; NoE [#]",100,0,1000,1000,0,10000);
	histograms2D.insert(std::make_pair("h_initMisPosVsQTotal",h_initMisPosVsQTotal));

	TH2D* h_initMisPosVsNHitsCaus = new TH2D("h_initMisPosVsNHitsCaus","Initial mismatch position vs number of hits after causality; Initial mismatch position [m]; N_{hits}^{caus} [#]; NoE [#]",100,0,1000,100,0,100);
	histograms2D.insert(std::make_pair("h_initMisPosVsNHitsCaus",h_initMisPosVsNHitsCaus));

	TH2D* h_initMisPosVsDistanceCS = new TH2D("h_initMisPosVsDistanceCS","Initial mismatch position vs distance from the CS; Initial mismatch position [m]; Distance_{CS} [m]; NoE [#]",100,0,1000,150,0,150);
	histograms2D.insert(std::make_pair("h_initMisPosVsDistanceCS",h_initMisPosVsDistanceCS));

	TH2D* h_initMisPosVsNNoiseHits = new TH2D("h_initMisPosVsNNoiseHits","Initial mismatch position vs number of noise hits; Initial mismatch position [m]; N_{noiseHits} [#]; NoE [#]",100,0,1000,20,0,20);
	histograms2D.insert(std::make_pair("h_initMisPosVsNNoiseHits",h_initMisPosVsNNoiseHits));

	TH2D* h_initMisPosVsChi2Init = new TH2D("h_initMisPosVsChi2Init","Initial mismatch position vs initial chi2; Initial mismatch position [m]; #chi^{2}_{init} [1]; NoE [#]",100,0,1000,1000,0,10000);
	histograms2D.insert(std::make_pair("h_initMisPosVsChi2Init",h_initMisPosVsChi2Init));

	TH2D* h_misPosVsNHits = new TH2D("h_misPosVsNHits","Mismatch position vs number of hits after TFilter; Mismatch position [m]; N_{hits}^{Tfil} [#]; NoE [#]",100,0,1000,100,0,100);
	histograms2D.insert(std::make_pair("h_misPosVsNHits",h_misPosVsNHits));


	return 0;
}


int DrawHistograms(std::map<TString,TH1D*> histograms,std::map<TString,TH2D*> &histograms2D)
{
	auto it{ histograms.cbegin() }; // declare a const iterator and assign to start of vector
    while (it != histograms.cend()) // while it hasn't reach the end
    {
        TCanvas* c_tempCanvas = new TCanvas(Form("c_%s",it->first.Data()),it->first.Data(),800,600);
        it->second->Draw();
        ++it; // and iterate to the next element
    }

    auto it2D{ histograms2D.cbegin() }; // declare a const iterator and assign to start of vector
	while (it2D != histograms2D.cend()) // while it hasn't reach the end
	{

		TCanvas* c_tempCanvas = new TCanvas(Form("c_%s",it2D->first.Data()),it2D->first.Data(),800,600);
		it2D->second->Draw("colz");
		++it2D; // and iterate to the next element
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

int PrintEvent(BEvent* myEvent, BRecoCascade* myCascade, BGeomTel* myGeom, BMCCascadeEvent* myMCCascade,std::map<TString,TH2D*> &histograms2D)
{
	cout << "EvENT: " << endl;
	for (int i = 0; i < myEvent->NHits(); ++i)
	{
		cout << setprecision(4) << i << "\t" << myEvent->HitChannel(i) << "\t" << myEvent->Q(i) << "\t" << myEvent->T(i) << "\t" << myMCCascade->GetPulse(myEvent->Hit(i).GetUniqueID())->GetMagic() << "\t\t\t";
		if (i < myCascade->GetNCausImpulseIDs())
			cout << setprecision(4) << myEvent->HitChannel(myCascade->GetCausImpulseID(i)) << "\t" << myMCCascade->GetPulse(myEvent->Hit(myCascade->GetCausImpulseID(i)).GetUniqueID())->GetMagic() << "\t" << myEvent->Q(myCascade->GetCausImpulseID(i)) << "\t\t" << myEvent->T(myCascade->GetCausImpulseID(i)) << "\t\t\t";
			// cout << setprecision(4) << myEvent->HitChannel(myCascade->GetCausImpulseID(i)) << "\t" << myMCCascade->GetPulse(myEvent->Hit(myCascade->GetCausImpulseID(i)).GetUniqueID())->GetMagic() << "\t" << myEvent->Q(myCascade->GetCausImpulseID(i)) << "\t\t" << myEvent->T(myCascade->GetCausImpulseID(i)) << "\t" << myGeom->GetX(myEvent->HitChannel(myCascade->GetCausImpulseID(i))) << "\t\t" << myGeom->GetY(myEvent->HitChannel(myCascade->GetCausImpulseID(i))) << "\t\t" << myGeom->GetZ(myEvent->HitChannel(myCascade->GetCausImpulseID(i))) << "\t" << (myCascade->GetPosMC()-myGeom->GetPosition(myEvent->HitChannel(myCascade->GetCausImpulseID(i)))).Mag()/BARS::Phys::c_baikal*1e9 - myEvent->T(myCascade->GetCausImpulseID(i)) << "\t" << (myCascade->GetFitPos()-myGeom->GetPosition(myEvent->HitChannel(myCascade->GetCausImpulseID(i)))).Mag()/BARS::Phys::c_baikal*1e9 - myEvent->T(myCascade->GetCausImpulseID(i)) << "\t" << (myCascade->GetFitPos()-myGeom->GetPosition(myEvent->HitChannel(myCascade->GetCausImpulseID(i)))).Mag() << "\t" << (myCascade->GetInitPos()-myGeom->GetPosition(myEvent->HitChannel(myCascade->GetCausImpulseID(i)))).Mag()/BARS::Phys::c_baikal*1e9 - myEvent->T(myCascade->GetCausImpulseID(i)) + myCascade->GetInitTime() << "\t" << (myGeom->GetPosition(myEvent->HitChannel(myCascade->GetCausImpulseID(0)))-myGeom->GetPosition(myEvent->HitChannel(myCascade->GetCausImpulseID(i)))).Mag();
		if (i < myCascade->GetNImpulseIDs())
			cout << setprecision(4) << myEvent->HitChannel(myCascade->GetImpulseID(i)) << "\t" << myMCCascade->GetPulse(myEvent->Hit(myCascade->GetImpulseID(i)).GetUniqueID())->GetMagic() << "\t" << myEvent->Q(myCascade->GetImpulseID(i)) << "\t\t" << myEvent->T(myCascade->GetImpulseID(i)) << "\t\t\t";
		cout << endl;
	}

	Int_t nSignalHitsTrue = 0;
	Int_t nNoiseHitsTrue = 0;
	for (int i = 0; i < myMCCascade->GetNHits(); ++i)
	{
		if (myMCCascade->GetPulse(i)->GetMagic() > 0)
			nSignalHitsTrue++;
		if (myMCCascade->GetPulse(i)->GetMagic() == 0)
			nNoiseHitsTrue++;
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
	Int_t nCausHits = myCascade->GetNCausImpulseIDs();

	cout << nSignalHitsTrue << "\t" << nNoiseHitsTrue << "\t" << nSignalHitsReco << "\t" << nSignalCausHitsReco << "\t" << nCausHits << endl;

	return 0;
}

int optimizeCuts(TString fileName = "")
{
	TChain reconstructedCascades("Events");
	SetTChain(reconstructedCascades,fileName);

	BRecoCascade* myCascade = NULL;
	BJointHeader* myHeader  = NULL;
	BEvent*       myEvent   = NULL;
	BGeomTel*     myGeom    = NULL;

	TBranch* b_joint = reconstructedCascades.FindBranch("BJointHeader.");
	TBranch* b_reco  = reconstructedCascades.FindBranch("BRecoCascade.");
	TBranch* b_event = reconstructedCascades.FindBranch("BEvent.");
	TBranch* b_geom  = reconstructedCascades.FindBranch("BGeomTel.");
    if(!b_joint || !b_reco || !b_event || !b_geom) {
        std::cout << "Could not find the BJointHeader or BRecoCascade branch on tree Events, cannot continue." << std::endl;
        return -1;
    }
	reconstructedCascades.SetBranchAddress("BRecoCascade.", &myCascade);
	reconstructedCascades.SetBranchAddress("BJointHeader.", &myHeader);
	reconstructedCascades.SetBranchAddress("BEvent.", &myEvent);
	reconstructedCascades.SetBranchAddress("BGeomTel.", &myGeom);

	BMCCascadeEvent* myMCCascade = NULL;
    TBranch* b_mcCascade = reconstructedCascades.FindBranch("BMCCascadeEvent.");
    Bool_t muonBundle = false;
    if (!b_mcCascade)
    {
    	muonBundle = true;
    // 	std::cout << "Could not find the BMCCascadeEvent branch on tree Events, cannot continue." << std::endl;
    //     return -1;
    }else
    {
    	reconstructedCascades.SetBranchAddress("BMCCascadeEvent.", &myMCCascade);
    }
    BMCEvent* myMCEvent = NULL;
    TBranch* b_mcEvent = reconstructedCascades.FindBranch("BMCEvent.");
    if (b_mcEvent)
    {
    	reconstructedCascades.SetBranchAddress("BMCEvent.", &myMCEvent);
    }

	std::map<TString, TH1D*> v_histograms;
	std::map<TString, TH2D*> v_histograms2D;
	SetHistograms(v_histograms,v_histograms2D);

	Float_t eventWeight = 1.0;

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);
		eventWeight = myCascade->GetEventWeight();

		Int_t nSignalCausHitsReco = 0;
		Int_t nSignalHitsTrue = 0;
		Int_t nNoiseHitsTrue = 0;
		Int_t nSignalHitsReco = 0;

		if (!muonBundle)
		{
			// if ((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag() > 50 && myCascade->GetNHitsCaus() > 30)
			// if ((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag() < 10 && myCascade->GetChi2Init() > 100 && myCascade->GetDistanceCSMC() < 60 && myCascade->GetPosMC().Z() < 250 && myCascade->GetPosMC().Z() > -250)
			// if (myCascade->GetChi2Init() > 1000)
			// if (myCascade->GetChi2Caus() > 80)
			 // if (i == 39961 || i == 39962)
			// if (myCascade->GetNHitsTFil()-myCascade->GetNHitsCaus() < 0)
			if (true)
			{
				cout << i << "\t" << myHeader->GetEventIDCC() << endl;
				cout << myCascade->GetChi2Init() << endl;
				cout << myCascade->GetInitTime() << endl;
				cout << myCascade->GetChi2Caus() << endl;
				cout << myCascade->GetInitPrec() << endl;
				myCascade->GetPosMC().Print();
				myCascade->GetInitPos().Print();
				myCascade->GetFitPos().Print();
				PrintEvent(myEvent,myCascade,myGeom,myMCCascade,v_histograms2D);
			}

			for (int i = 0; i < myMCCascade->GetNHits(); ++i)
			{
				if (myMCCascade->GetPulse(i)->GetMagic() > 0)
					nSignalHitsTrue++;
				if (myMCCascade->GetPulse(i)->GetMagic() == 0)
					nNoiseHitsTrue++;
			}
			for (int i = 0; i < myCascade->GetNImpulseIDs(); ++i)
			{
				if (myMCCascade->GetPulse(myEvent->Hit(myCascade->GetImpulseID(i)).GetUniqueID())->GetMagic() > 0)
					nSignalHitsReco++;
			}
			for (int i = 0; i < myCascade->GetNCausImpulseIDs(); ++i)
			{
				if (myMCCascade->GetPulse(myEvent->Hit(myCascade->GetCausImpulseID(i)).GetUniqueID())->GetMagic() > 0)
				{
					nSignalCausHitsReco++;
					v_histograms["h_zDistanceSignalHits"]->Fill(TMath::Abs(myGeom->GetZ(myEvent->HitChannel(myCascade->GetCausImpulseID(i)))-myGeom->GetZ(myEvent->HitChannel(myCascade->GetCausImpulseID(0)))));
				}else
				{
					v_histograms["h_zDistanceNoiseHits"]->Fill(TMath::Abs(myGeom->GetZ(myEvent->HitChannel(myCascade->GetCausImpulseID(i)))-myGeom->GetZ(myEvent->HitChannel(myCascade->GetCausImpulseID(0)))));
				}

			}
		}

		Int_t nCausHits = myCascade->GetNCausImpulseIDs();

		v_histograms["h_chargePerHit"]->Fill(myCascade->GetQTotal()/myCascade->GetNHitsCaus(),eventWeight);
		v_histograms["h_nHitsCaus"]->Fill(myCascade->GetNHitsCaus(),eventWeight);
		v_histograms["h_nHitsFit"]->Fill(myCascade->GetNHitsTFil(),eventWeight);
		v_histograms["h_nHitsDiff"]->Fill(myCascade->GetNHitsTFil()-myCascade->GetNHitsCaus(),eventWeight);
		v_histograms["h_chi2"]->Fill(myCascade->GetChi2Caus(),eventWeight);
		v_histograms2D["h_initMisPosVsChi2"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),myCascade->GetChi2Caus(),eventWeight);
		v_histograms2D["h_initMisPosVsInitPrec"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),TMath::Sqrt(myCascade->GetInitPrec()),eventWeight);
		v_histograms2D["h_initMisPosVsQTotal"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),myCascade->GetQTotal(),eventWeight);
		v_histograms2D["h_initMisPosVsNHitsCaus"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),myCascade->GetNHitsCaus(),eventWeight);
		v_histograms2D["h_misPosVsNHits"]->Fill((myCascade->GetPosMC()-myCascade->GetFitPos()).Mag(),myCascade->GetNHitsTFil(),eventWeight);
		v_histograms2D["h_initMisPosVsDistanceCS"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),myCascade->GetDistanceCS(),eventWeight);
		v_histograms["h_initMismatchPosition"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),eventWeight);
		v_histograms2D["h_initMisPosVsChi2Init"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),myCascade->GetChi2Init(),eventWeight);
		v_histograms["h_fitMismatchPosition"]->Fill((myCascade->GetPosMC()-myCascade->GetFitPos()).Mag(),eventWeight);
		if (!muonBundle)
		{
			v_histograms2D["h_initMisPosVsNNoiseHits"]->Fill((myCascade->GetPosMC()-myCascade->GetInitPos()).Mag(),nCausHits-nSignalCausHitsReco);
			v_histograms["h_causHitPurity"]->Fill((double)nSignalCausHitsReco/myCascade->GetNCausImpulseIDs()*100,eventWeight);
			v_histograms["h_causHitEfficiency"]->Fill((double)nSignalCausHitsReco/nSignalHitsTrue*100,eventWeight);
		}
	}

	DrawHistograms(v_histograms,v_histograms2D);
	SaveHistograms(v_histograms,fileName);

	return 0;
}