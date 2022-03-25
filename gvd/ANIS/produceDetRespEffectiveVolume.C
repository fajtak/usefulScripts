#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLegend.h"

int SetHistograms(std::map<std::string,TH1D*> &histograms,std::map<std::string,TH2D*> &histograms2D)
{
    TH1D* h_xPos = new TH1D("h_xPos","Interaction X position; x [m]; NoE [#]",2000,-1000,1000);
    histograms.insert(std::make_pair("h_xPos",h_xPos));

    TH1D* h_yPos = new TH1D("h_yPos","Interaction Y position; y [m]; NoE [#]",2000,-1000,1000);
    histograms.insert(std::make_pair("h_yPos",h_yPos));

    TH1D* h_zPos = new TH1D("h_zPos","Interaction Z position; z [m]; NoE [#]",2000,-1000,1000);
    histograms.insert(std::make_pair("h_zPos",h_zPos));

    TH2D* h_xVsyPos = new TH2D("h_xVsyPos","Interaction X,Y position; x [m]; y [m]; NoE [#]",200,-1000,1000,200,-1000,1000);
    histograms2D.insert(std::make_pair("h_xVsyPos",h_xVsyPos));

    TH2D* h_xVszPos = new TH2D("h_xVszPos","Interaction X,Z position; x [m]; z [m]; NoE [#]",200,-1000,1000,200,-1000,1000);
    histograms2D.insert(std::make_pair("h_xVszPos",h_xVszPos));

    TH2D* h_yVszPos = new TH2D("h_yVszPos","Interaction Y,Z position; y [m]; z [m]; NoE [#]",200,-1000,1000,200,-1000,1000);
    histograms2D.insert(std::make_pair("h_yVszPos",h_yVszPos));

    TH1D* h_distPos = new TH1D("h_distPos","Distance interaction point and origin; D [m]; NoE [#]",1000,0,1000);
    histograms.insert(std::make_pair("h_distPos",h_distPos));

    TH1D* h_ENuIn = new TH1D("h_ENuIn","Incoming Neutrino energy; E [TeV]; NoE [#]",1000,0,10000);
    histograms.insert(std::make_pair("h_ENuIn",h_ENuIn));

    TH1D* h_ENuInWeighted = new TH1D("h_ENuInWeighted","Incoming Neutrino energy (weighted); E_{nu} [TeV]; NoE [#]",1000,0,10000);
    histograms.insert(std::make_pair("h_ENuInWeighted",h_ENuInWeighted));

    TH1D* h_ECascade = new TH1D("h_ECascade","Cascade energy; E_{cascade} [TeV]; NoE [#]",1000,0,10000);
    histograms.insert(std::make_pair("h_ECascade",h_ECascade));

    TH1D* h_ECascadeWeighted = new TH1D("h_ECascadeWeighted","Cascade energy (weighted); E_{cascade} [TeV]; NoE [#]",1000,0,10000);
    histograms.insert(std::make_pair("h_ECascadeWeighted",h_ECascadeWeighted));

    TH2D* h_ENuInVsECascade = new TH2D("h_ENuInVsECascade","Neutrino energy vs. cascade energy; E_{nu} [TeV]; E_{cascade} [TeV]; NoE [#]",1000,0,10000,1000,0,10000);
    histograms2D.insert(std::make_pair("h_ENuInVsECascade",h_ENuInVsECascade));

    TH1D* h_cosTheta = new TH1D("h_cosTheta","Cascade cos(#theta); cos(#theta) [1]; NoE [#]",200,-1,1);
    histograms.insert(std::make_pair("h_cosTheta",h_cosTheta));

    TH1D* h_cosThetaWeighted = new TH1D("h_cosThetaWeighted","Cascade cos(#theta); cos(#theta) [1]; NoE [#]",200,-1,1);
    histograms.insert(std::make_pair("h_cosThetaWeighted",h_cosThetaWeighted));

    TH1D* h_phi = new TH1D("h_phi","Cascade #phi; #phi [rad]; NoE [#]",2000,-10,10);
    histograms.insert(std::make_pair("h_phi",h_phi));

    TH1D* h_phiWeighted = new TH1D("h_phiWeighted","Cascade #phi; #phi [rad]; NoE [#]",2000,-10,10);
    histograms.insert(std::make_pair("h_phiWeighted",h_phiWeighted));

    TH1D* h_weight = new TH1D("h_weight","Event weights; weight [1]; NoE [#]",1000,0,0.001);
    histograms.insert(std::make_pair("h_weight",h_weight));

    TH1D* h_logWeight = new TH1D("h_logWeight","Event weights; log_{10}(weight) [1]; NoE [#]",100,-20,1);
    histograms.insert(std::make_pair("h_logWeight",h_logWeight));

    TH2D* h_logWeightVsDistance = new TH2D("h_logWeightVsDistance","Event weight vs. distance; D [m]; log_{10}(Weight) [1]; NoE [#]",100,0,1000,100,-20,1);
    histograms2D.insert(std::make_pair("h_logWeightVsDistance",h_logWeightVsDistance));

    TH2D* h_logWeightVslogE = new TH2D("h_logWeightVslogE","Event weight vs. energy; log_{10}(E [TeV]); log_{10}(Weight) [1]; NoE [#]",100,0,10,100,-20,1);
    histograms2D.insert(std::make_pair("h_logWeightVslogE",h_logWeightVslogE));

    TH2D* h_logWeightVsCosTheta = new TH2D("h_logWeightVsCosTheta","Event weight vs. cosTheta; cos(#theta) [1]; log_{10}(Weight) [1]; NoE [#]",200,-1,1,100,-20,1);
    histograms2D.insert(std::make_pair("h_logWeightVsCosTheta",h_logWeightVsCosTheta));

    TH2D* h_logWeightVsPhi = new TH2D("h_logWeightVsPhi","Event weight vs. #phi; #phi [rad]; log_{10}(Weight) [1]; NoE [#]",200,-10,10,100,-20,1);
    histograms2D.insert(std::make_pair("h_logWeightVsPhi",h_logWeightVsPhi));

    TH2D* h_genEventsVslogEVsCosTheta = new TH2D("h_genEventsVslogEVsCosTheta","Number of generated events vs. E vs cos(#theta); log_{10}(E [TeV]); cos(#theta) [1]",100,0,10,20,-1,1);
    histograms2D.insert(std::make_pair("h_genEventsVslogEVsCosTheta",h_genEventsVslogEVsCosTheta));

    TH2D* h_genEventsVslogEVsCosThetaWeighted = new TH2D("h_genEventsVslogEVsCosThetaWeighted","Number of generated events vs. E vs cos(#theta); log_{10}(E [TeV]); cos(#theta) [1]",100,0,10,20,-1,1);
    histograms2D.insert(std::make_pair("h_genEventsVslogEVsCosThetaWeighted",h_genEventsVslogEVsCosThetaWeighted));

    return 0;
}

int DrawHistograms(std::map<std::string,TH1D*> &histograms, std::map<std::string,TH2D*> &histograms2D)
{
    auto it{ histograms.cbegin() }; // declare a const iterator and assign to start of vector
    while (it != histograms.cend()) // while it hasn't reach the end
    {
        TCanvas* c_tempCanvas = new TCanvas(Form("c_%s",it->first.c_str()),it->first.c_str(),800,600);
        it->second->Draw();
        ++it; // and iterate to the next element
    }

    auto it2D{ histograms2D.cbegin() }; // declare a const iterator and assign to start of vector
    while (it2D != histograms2D.cend()) // while it hasn't reach the end
    {
        TCanvas* c_tempCanvas = new TCanvas(Form("c_%s",it2D->first.c_str()),it2D->first.c_str(),800,600);
        it2D->second->Draw("colz");
        if (it2D->first == "h_logWeightVsDistance" || it2D->first == "h_logWeightVslogE" || it2D->first == "h_logWeightVsCosTheta" || it2D->first == "h_logWeightVsPhi")
            it2D->second->QuantilesX()->Draw("SAME");
        ++it2D; // and iterate to the next element
    }

    return 0;
}

int SaveHistograms(std::map<std::string,TH1D*> &histograms, std::map<std::string,TH2D*> &histograms2D, TString fileName)
{
    TString outputFileName = "";
    outputFileName = fileName(0,fileName.Last('/')) + "/effVolumeDetResp_" + fileName(fileName.Last('/')+1,fileName.Length()) + ".root";
    TFile* f_results = new TFile(outputFileName,"RECREATE");

    histograms2D["h_genEventsVslogEVsCosTheta"]->Write();
    histograms2D["h_genEventsVslogEVsCosThetaWeighted"]->Write();

    return 0;
}

int ReadAndFillHistograms(TString fileName, std::map<std::string,TH1D*> &histograms, std::map<std::string,TH2D*> &histograms2D, int nEvents)
{
    int eventID,interactionType,nHits;
    float eventWeight=0,ENuZero=0,ENuIn=0,ECasc=0,V1x=0,V1y=0,V1z=0,p1x=0,p1y=0,p1z=0,cosTheta=0,phi=0,ECascade=0;
    string holder;
    ifstream readFile;

    Float_t     fPulseT;
    Float_t     fPulseQ;
    Float_t     fPulseOMID;
    Float_t     fPulseMagic;

    readFile.open(fileName);
    if(!readFile.is_open()){
        readFile.close();
        return -1;
    }

    while(readFile>>eventID){
        if (nEvents != -1 && eventID > nEvents)
            break;
        readFile>>eventWeight>>interactionType>>ENuZero>>ENuIn;
        getline(readFile,holder);
        readFile>>nHits>>cosTheta>>phi>>V1x>>V1y>>V1z>>ECascade;
        getline(readFile,holder);

        for (int i = 0; i < nHits; ++i)
        {
          readFile>>fPulseOMID >> fPulseT >> fPulseQ;
          // pulse->SetPulse(fChID[i]-1,fCharge[fChID[i]-1],fTime[fChID[i]-1],1);
        }
        getline(readFile,holder);

        TVector3 dir;
        dir.SetMagThetaPhi(1,TMath::Pi()-TMath::ACos(cosTheta),(phi+TMath::Pi() >= 2*TMath::Pi())?phi-TMath::Pi():phi+TMath::Pi());
        TVector3 pos(V1x,V1y,V1z);
        histograms["h_xPos"]->Fill(V1x);
        histograms["h_yPos"]->Fill(V1y);
        histograms["h_zPos"]->Fill(V1z);
        histograms["h_distPos"]->Fill(pos.Mag());
        histograms2D["h_xVsyPos"]->Fill(V1x,V1y);
        histograms2D["h_xVszPos"]->Fill(V1x,V1z);
        histograms2D["h_yVszPos"]->Fill(V1y,V1z);
        histograms["h_ENuIn"]->Fill(ENuIn/1000);
        histograms["h_ECascade"]->Fill(dir.Mag()/1000);
        // if (p1z > 0)
        histograms["h_ENuInWeighted"]->Fill(ENuIn/1000,eventWeight);
        histograms["h_ECascadeWeighted"]->Fill(dir.Mag()/1000,eventWeight);
        histograms2D["h_ENuInVsECascade"]->Fill(ENuIn/1000,dir.Mag()/1000);
        histograms["h_cosTheta"]->Fill(dir.CosTheta());
        histograms["h_cosThetaWeighted"]->Fill(dir.CosTheta(),eventWeight);
        histograms["h_phi"]->Fill(dir.Phi());
        histograms["h_phiWeighted"]->Fill(dir.Phi(),eventWeight);
        histograms["h_weight"]->Fill(eventWeight);
        histograms["h_logWeight"]->Fill(TMath::Log10(eventWeight));
        histograms2D["h_logWeightVsDistance"]->Fill(pos.Mag(),TMath::Log10(eventWeight));
        histograms2D["h_logWeightVslogE"]->Fill(TMath::Log10(ENuIn/1000),TMath::Log10(eventWeight));
        histograms2D["h_logWeightVsCosTheta"]->Fill(dir.CosTheta(),TMath::Log10(eventWeight));
        histograms2D["h_logWeightVsPhi"]->Fill(dir.Phi(),TMath::Log10(eventWeight));
        histograms2D["h_genEventsVslogEVsCosTheta"]->Fill(TMath::Log10(ENuIn/1000),dir.CosTheta());
        histograms2D["h_genEventsVslogEVsCosThetaWeighted"]->Fill(TMath::Log10(ENuIn/1000),dir.CosTheta(),eventWeight);

        holder.clear();
    }
    readFile.close();

    return 0;
}

int produceDetRespEffectiveVolume(TString fileName = "/Data/BaikalData/mc/ANIS/nueastro_ver4_50kNoise/astro_nue_ver4.dat",int nEvents = -1)
{
    std::map<std::string, TH1D*> v_histograms;
    std::map<std::string, TH2D*> v_histograms2D;

    SetHistograms(v_histograms,v_histograms2D);

    if (ReadAndFillHistograms(fileName,v_histograms,v_histograms2D,nEvents) < 0)
        return -1;

    DrawHistograms(v_histograms,v_histograms2D);
    SaveHistograms(v_histograms,v_histograms2D,fileName);

    return 0;
}	