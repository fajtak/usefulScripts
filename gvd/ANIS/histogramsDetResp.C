#include <iostream>
#include <fstream>
#include <iomanip>
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

using namespace std;

TH1D* h_thetaRec = new TH1D("h_thetaRec","Reconstructed theta; #theta [deg.]; NoE [#]",36,0,180);

TH1D* h_cosThetaRec = new TH1D("h_cosThetaRec","Reconstructed cos(theta); cos(#theta) [1]; NoE [#]",22,-1.1,1.1);
TH1D* h_cosThetaRecNW = new TH1D("h_cosThetaRecNW","Reconstructed cos(theta); cos(#theta) [1]; NoE [#]",22,-1.1,1.1);

TH2D* h_cosThetaMCvsZ = new TH2D("h_cosThetaMCvsZ",";cos(#theta) [1]; Z [m]",22,-1.1,1.1,100,-500,500);

TH1D* h_EMC = new TH1D("h_EMC","MC energy; E_{MC} [TeV]; NoE*eventWeight [#]",200,0,1000);

void histogramsDetResp()
{
  std::ifstream fTxtFile;
  // fTxtFile.open("/Data/BaikalData/mc/ANIS/nueatm_ver3_50kNoise/output_single_cascades_2022.dat");
  // fTxtFile.open("/Data/BaikalData/mc/ANIS/nueastro_ver4/query2.dat");
  // fTxtFile.open("/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/output_single_cascades_2022.dat");
  // fTxtFile.open("/Data/BaikalData/mc/ANIS/numuatm_ver4_50kNoise/output_single_cascades_2022.dat");
  fTxtFile.open("/Data/BaikalData/mc/ANIS/nueastro_ver4_50kNoise/astro_nue_ver4.dat");
  // fTxtFile.open("/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/astro_nue_ver4.dat");
  // fTxtFile.open("/Data/BaikalData/mc/ANIS/elec-output-single-cascades.dat");
  // fTxtFile.open("/Data/BaikalData/mc/ANIS/nueastro_ver1/output_single_cascades_2021.dat");
    Int_t       fProcessedEvents;
    Int_t       fEventID;
    Int_t       fNHits;
    Float_t     fNeutrinoEnergy;
    Float_t     fShowerEnergy;
    Float_t     fCosTheta;
    Float_t     fPhi;
    Float_t     fPosition[3];
    Float_t     fCharge[288];
    Float_t     fTime[288];
    UShort_t    fChID[288];
    Float_t     fWeight;
    Float_t     fIntProb;
    Float_t     fSPWE2;
    Float_t     fPulseT;
    Float_t     fPulseQ;
    Float_t     fPulseOMID;
    Float_t     fPulseMagic;
    double Ezero=0;
    string holder;
    //loop for reading data till the end
    while(kTRUE)
    {
      if (fTxtFile >> fEventID)
      {
        fProcessedEvents++;
      }else
      {
        break;
      }

      fTxtFile>>fWeight>>fPulseMagic>>Ezero>>fNeutrinoEnergy;
      getline(fTxtFile,holder);
      fTxtFile>>fNHits>>fCosTheta>>fPhi>>fPosition[0]>>fPosition[1]>>fPosition[2]>>fShowerEnergy;
      getline(fTxtFile,holder);

      h_cosThetaRec->Fill(-1*fCosTheta,fWeight);
      h_cosThetaRecNW->Fill(-1*fCosTheta);
      h_EMC->Fill(fShowerEnergy*1E-3);
      h_cosThetaMCvsZ->Fill(-1*fCosTheta,fPosition[2]);

      // fEvent->SetNeutrinoEnergy(fNeutrinoEnergy*1E-3);
      // fEvent->SetShowerEnergy(fShowerEnergy*1E-3);
      // fEvent->SetNeutrinoTheta(TMath::Pi()-TMath::ACos(fCosTheta));
      // fEvent->SetNeutrinoAzimuth((fPhi+TMath::Pi() >= 2*TMath::Pi())?fPhi-TMath::Pi():fPhi+TMath::Pi());
      // fEvent->SetInteractionPosition(fPosition[0],fPosition[1],fPosition[2]);
      // fEvent->SetEventWeight(fWeight);
      // fEvent->SetIntProb(fIntProb);
      // fEvent->SetSPWE2(fSPWE2);

      // fCascade->SetPosMC(TVector3(fPosition[0],fPosition[1],fPosition[2]));
      // fCascade->SetEnergyMC(fShowerEnergy*1E-3);
      // fCascade->SetThetaMC((TMath::Pi()-TMath::ACos(fCosTheta))* TMath::RadToDeg());
      // fCascade->SetPhiMC(((fPhi+TMath::Pi() >= 2*TMath::Pi())?fPhi-TMath::Pi():fPhi+TMath::Pi())* TMath::RadToDeg());
      // fCascade->SetEventFlagMC(fPulseMagic+20);
      // fCascade->SetEventWeight(fWeight);

      // TVector3 referencePosition2D = (fGeomTel->GetPosition(270)+fGeomTel->GetPosition(269))*0.5;
      // referencePosition2D.SetZ(0.0);
      // TVector3 cascadePosition2D(fPosition[0],fPosition[1],0.0);

      // fCascade->SetDistanceCSMC((referencePosition2D-cascadePosition2D).Mag());

      for (int i = 0; i < fNHits; ++i)
      {
        // BMCCascadePulse* pulse = fEvent->AddPulse(i);
        fTxtFile>>fPulseOMID >> fPulseT >> fPulseQ;
        // pulse->SetOMID(fPulseOMID-1);
        // pulse->SetMagic(fPulseMagic+20);
        // pulse->SetAmplitude(fPulseQ);
        // pulse->SetTime(fPulseT);
        // pulse->SetPulse(fChID[i]-1,fCharge[fChID[i]-1],fTime[fChID[i]-1],1);
      }
      getline(fTxtFile,holder);
    }

    TCanvas* c_cosThetaRec = new TCanvas("c_cosThetaRec","Theta",800,600);
  h_cosThetaRec->Draw();

  TCanvas* c_cosThetaRecNW = new TCanvas("c_cosThetaRecNW","ThetaNW",800,600);
  h_cosThetaRecNW->Draw();

  TCanvas* c_EMC = new TCanvas("c_EMC","EMC",800,600);
  h_EMC->Draw();

    TCanvas* c_cosThetaMCvsZ = new TCanvas("c_cosThetaMCvsZ","CosThetaMCVsZ",800,600);
    h_cosThetaMCvsZ->Draw("colz");
}