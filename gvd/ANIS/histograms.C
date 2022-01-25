#include <iostream>
#include <fstream>
#include <iomanip>
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLegend.h"

using namespace std;

TH1D* h_thetaRec = new TH1D("h_thetaRec","Reconstructed theta; #theta [deg.]; NoE [#]",36,0,180);

TH1D* h_EMC = new TH1D("h_EMC","MC energy; E_{MC} [TeV]; NoE*eventWeight [#]",200,0,1000);
TH1D* h_cosThetaMC = new TH1D("h_cosThetaMC","MC cos(theta); cos(#theta_{MC}) [1]; NoE*eventWeight [#]",22,-1.1,1.1);
TH1D* h_cosThetaRecNW = new TH1D("h_cosThetaRecNW","Reconstructed cos(theta); cos(#theta) [1]; NoE [#]",22,-1.1,1.1);

void histograms()
{
  int Nev,nFlag;
  float hScale=1,weight=0,Ezero=0,Ein=0,Enr=0,V1x=0,V1y=0,V1z=0,p1x=0,p1y=0,p1z=0;
  float cumWeightsGeV=0,cumWeightsTeV=0;
  string holder;
  ifstream readFile;

  // readFile.open("muoncascade.dat");
  // readFile.open("/Data/BaikalData/mc/ANIS/elecANISoutput");
  // readFile.open("/Data/BaikalData/mc/ANIS/nueatm_ver2/ANISoutput");
  // readFile.open("/Data/BaikalData/mc/ANIS/nueatm_ver3_50kNoise/ANISoutput");
  // readFile.open("/Data/BaikalData/mc/ANIS/nueastro_ver4/ANISoutput");
  readFile.open("/Data/BaikalData/mc/ANIS/nueatm_ver4_50kNoise/ANISoutput");
  // readFile.open("/Data/BaikalData/mc/ANIS/nueastro_ver1/ANISoutput");
  if(!readFile.is_open()){ 
  readFile.close();
  return; 
  }

  while(readFile>>Nev){
    readFile>>weight>>nFlag>>Ezero>>Ein>>V1x>>V1y>>V1z>>p1x>>p1y>>p1z;
    getline(readFile,holder);
    TVector3 dir(p1x,p1y,p1z);
    holder.clear();
    h_thetaRec->Fill(dir.Theta()*TMath::RadToDeg(),weight);
    h_cosThetaMC->Fill(dir.CosTheta(),weight);
    h_cosThetaRecNW->Fill(dir.CosTheta());
    h_EMC->Fill(Ein*1E-3,weight);
    // if(nFlag==3){
      // cumWeightsGeV=cumWeightsGeV+weight;
    // }
  }
  readFile.close();
  
  TCanvas* c_cosThetaMC = new TCanvas("c_cosThetaMC","CosThetaMC",800,600);
  h_cosThetaMC->Draw();

  TCanvas* c_cosThetaRecNW = new TCanvas("c_cosThetaRecNW","ThetaNW",800,600);
  h_cosThetaRecNW->Draw();

  TCanvas* c_EMC = new TCanvas("c_EMC","EMC",800,600);
  h_EMC->Draw();


  // cout<<endl<<" sum of the weights is "<<cumWeightsGeV<<endl;

  /*
  
  TH1F *wElec = new TH1F("Elec nu -  casc energy", " log E_e",   18, 2.5, 6.5);
  TH1F *wMiu  = new TH1F("Muon nu -  casc energy", " log E_m",   18, 2.5, 6.5);
  TH1F *wAstr = new TH1F("Astr nu -  casc energy", " log E_a",   18, 2.5, 6.5);

  TH1F *wElecGeV = new TH1F("Elec nu -  casc energy GeV ", " E_e",   8, 1e3, 4e5);
  TH1F *wElecTeV = new TH1F("Elec nu -  casc energy TeV ", " E_e",   8, 1, 4e2);
  

  readFile.open("elec.dat");
  if(!readFile.is_open()){ 
	readFile.close();
	return; 
  }

  while(readFile>>Nev){
    readFile>>weight>>nFlag>>Ezero>>Ein>>V1x>>V1y>>V1z>>p1x>>p1y>>p1z;
    getline(readFile,holder);
    holder.clear();
    Enr=TMath::Power(p1x*p1x+p1y*p1y+p1z*p1z,0.5);
    wElec->Fill(TMath::Log10(Enr),weight);
    wElecGeV->Fill(Enr,weight);
    cumWeightsGeV+=weight;
    cumWeightsTeV+=weight;
    wElecTeV->Fill(1.e-3*Enr,1e+3*weight);
  }
  readFile.close();
 
  cout<<" GeV cummulative sum = "<< cumWeightsGeV<<endl;
  cout<<" TeV cummulative sum = "<< cumWeightsTeV<<endl;
  cout<<" GeV histogram integral = "<< wElecGeV->Integral()<<endl;
  cout<<" GeV histogram integral with widths= "<< wElecGeV->Integral("width")<<endl;
  cout<<" TeV histogram integral = "<< wElecTeV->Integral()<<endl;
  cout<<" TeV histogram integral with widths= "<< wElecTeV->Integral("width")<<endl;

  TCanvas *cGeV = new TCanvas("energia","energy scale");
  wElecGeV->Draw();

  TCanvas *cTeV = new TCanvas("energiaTeV","TeV energy scale");
  wElecTeV->Draw();

  
  readFile.open("muon.dat");
  if(!readFile.is_open()){ 
  readFile.close();
  return; 
  }

  while(readFile>>Nev){
    readFile>>weight>>nFlag>>Ezero>>Ein>>V1x>>V1y>>V1z>>p1x>>p1y>>p1z;
    getline(readFile,holder);
    holder.clear();
    Enr=TMath::Power(p1x*p1x+p1y*p1y+p1z*p1z,0.5);
    wMiu->Fill(TMath::Log10(Enr),weight);
  }
  readFile.close();
  
  readFile.open("astr.dat");
  if(!readFile.is_open()){ 
  readFile.close();
  return; 
  }
  
  int nEvents=0;
  float cumSum=0.;

  while(readFile>>Nev){
    readFile>>weight>>nFlag>>Ezero>>Ein>>V1x>>V1y>>V1z>>p1x>>p1y>>p1z;
    getline(readFile,holder);
    holder.clear();
    cumSum+=weight;
    nEvents++;
    Enr=TMath::Power(p1x*p1x+p1y*p1y+p1z*p1z,0.5);
    wAstr->Fill(TMath::Log10(Enr),weight); 
  }
  readFile.close();
  
  cout<<endl<<" events number = "<< nEvents <<" sum weights = "<<cumSum<<endl;

  TCanvas *cEnergy = new TCanvas("energy","energy distribution");
  TLegend *IncomingNeuts = new TLegend(.4,.7,.6,.9);
  
  wElec->GetXaxis()->SetTitle(" log E/GeV ");
  wElec->GetYaxis()->SetTitle(" Events [a.u.]");
  wElec->SetLineWidth(4);
  wElec->Draw();
  wMiu->SetLineColor(2);
  wMiu->SetLineWidth(4);
  wMiu->SetMarkerColor(2);
  wMiu->Draw("same");
  wAstr->SetLineColor(4);
  wAstr->SetLineWidth(4);
  wAstr->SetMarkerColor(4);
  wAstr->Draw("same");
  IncomingNeuts->AddEntry(wElec," atm. #nu_{e} ");
  IncomingNeuts->AddEntry(wMiu," atm. #nu_{#mu}  ");
  IncomingNeuts->AddEntry(wAstr," astr. #nu_{e} ");
  IncomingNeuts->Draw();
  */
  return ; 
}	