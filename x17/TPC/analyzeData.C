#include "TFile.h"
#include "TTree.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TH1F.h"

#include <iostream>

const int nWaveforms = 130;
const int nSamples = 1005;
const int nPedSamples = 50;
const int nSigmasIn = 5;
const int nSigmasOut = 5;


// using namespace std;

struct Hit
{
	Int_t channel;
	Int_t timeStart;
	Int_t timeEnd;
	Double_t charge;
	Double_t amplitude;
	Double_t pedestal;
	Double_t x;
	Double_t y;
};

TMultiGraph* VisWaveforms(Int_t Waveforms[][nSamples], Int_t size)
{
	TMultiGraph* mg_new = new TMultiGraph();
	for (int i = 0; i < size; ++i)
	{
		TGraph* g_temp = new TGraph();
		for (int j = 2; j < nSamples; ++j)
		{
			g_temp->SetPoint(j,j,Waveforms[i][j]);
		}
		mg_new->Add(g_temp,"lp");
	}
	return mg_new;
}

void GetPedestal(Int_t* waveform, Double_t &ped, Double_t &rms)
{
	ped = 0;
	rms = 0;
	Double_t pedSquared = 0;
	for (int i = 2; i < nPedSamples+2; ++i)
	{
		ped += waveform[i];
		pedSquared += waveform[i]*waveform[i];
	}
	rms = TMath::Sqrt((nPedSamples*pedSquared - ped*ped)/(nPedSamples*(nPedSamples-1)));
	ped /= nPedSamples;
}

void ExtractHits(Int_t Waveforms[][nSamples],Int_t* gl_Chn, Double_t* x, Double_t* y, Int_t size, vector<Hit> &hits)
{
	hits.clear();
	Double_t ped = 0;
	Double_t rms = 0;
	Double_t q = 0;
	Double_t a = 0;
	for (int i = 0; i < size; ++i)
	{
		GetPedestal(Waveforms[i],ped,rms);
		cout << ped << "\t" << rms << endl;
		for (int j = 2; j < nSamples; ++j)
		{
			if (Waveforms[i][j] > ped + nSigmasIn*rms)
			{
				q = 0;
				a = 0;
				for (int k = j; k < nSamples; ++k)
				{
					if (Waveforms[i][k] > ped + nSigmasOut*rms)
					{
						q += Waveforms[i][k]-ped;
						if (a < Waveforms[i][k]-ped)
						{
							a = Waveforms[i][k]-ped;
						}
					}else
					{
						if (j != k-1 && gl_Chn[i] != 58 && gl_Chn[i] != 12 && gl_Chn[i] != 1)
							hits.push_back(Hit{gl_Chn[i],j,k-1,q,a,ped,x[i],y[i]});
						j = k;
						break;
					}
				}
				cout << j << "\t" << gl_Chn[i] << "\t" << x[i] << "\t" << y[i] << endl;
			}
		}
	}
}

void PrintHits(const vector<Hit> &hits)
{
	for (unsigned int i = 0; i < hits.size(); ++i)
	{
		cout << i << "\t" << hits[i].channel << "\t" << hits[i].timeStart << "\t" << hits[i].timeEnd << "\t" << hits[i].charge << "\t" << hits[i].amplitude << "\t" <<hits[i].pedestal << "\t" << hits[i].x << "\t" << hits[i].y << endl;
	}
}

void DrawHits(const vector<Hit> &hits)
{
	TCanvas* c_event = new TCanvas("c_event","Event",800,600);
	TGraph2D *dt = new TGraph2D();
	Int_t nPoints = 0;
	for (unsigned int i = 0; i < hits.size(); ++i)
	{
		// if (hits[i].amplitude > 80 && hits[i].timeEnd-hits[i].timeStart > 4)
		{
			dt->SetPoint(nPoints,hits[i].x,hits[i].y,hits[i].timeStart);
			nPoints++;
		}
	}
	gStyle->SetPalette(1);
	dt->SetMarkerStyle(20);
	dt->Draw("pcol");

	TCanvas* c_event_xy = new TCanvas("c_event_xy","Event_XY",800,600);
	TGraph *g_xy = new TGraph(dt->GetN(),dt->GetX(),dt->GetY());
	g_xy->Draw("AP");

	TCanvas* c_event_xz = new TCanvas("c_event_xz","Event_XZ",800,600);
	TGraph *g_xz = new TGraph(dt->GetN(),dt->GetX(),dt->GetZ());
	g_xz->Draw("AP");

	TCanvas* c_event_yz = new TCanvas("c_event_yz","Event_YZ",800,600);
	TGraph *g_yz = new TGraph(dt->GetN(),dt->GetY(),dt->GetZ());

	g_yz->Draw("AP");
}

void AnalyzeHits(const vector<Hit> &hits)
{
	TH1F* h_q = new TH1F("h_q","Charge;Q [channels]; NoE [#]",1000,0,100000);
	TH1F* h_a = new TH1F("h_a","Amplitude;A [channels]; NoE [#]",1200,0,1200);
	TH1F* h_t = new TH1F("h_t","Time duration;T [channels]; NoE [#]",1000,0,1000);

	for (unsigned int i = 0; i < hits.size(); ++i)
	{
		// if (hits[i].amplitude > 60)
		{
			h_q->Fill(hits[i].charge);
			h_a->Fill(hits[i].amplitude);
			h_t->Fill(hits[i].timeEnd-hits[i].timeStart);
		}
	}
	TCanvas* c_q = new TCanvas("c_q","Charge",800,600);
	h_q->Draw();
	TCanvas* c_a = new TCanvas("c_a","Amplitude",800,600);
	h_a->Draw();
	TCanvas* c_t = new TCanvas("c_t","Time",800,600);
	h_t->Draw();

}

int analyzeData(Int_t eventID = 50)
{
	TString filePath = "/Data/x17/TPC/Run42_clean.root";

	TFile* inputFile = new TFile(filePath);
	TTree* waveforms = (TTree*)inputFile->Get("Waveforms");

	Long64_t 	triggerID;
	Int_t  		BX_counter, size;
	Int_t  		Sampa[nWaveforms], Chn[nWaveforms], Pad[nWaveforms], gl_Chn[nWaveforms], Waveform[nWaveforms][nSamples];
	Double_t  	x[nWaveforms],y[nWaveforms];

	waveforms->SetBranchAddress("trgID", &triggerID);
	waveforms->SetBranchAddress("BX_counter", &BX_counter);
	waveforms->SetBranchAddress("size", &size);
	waveforms->SetBranchAddress("Sampa", Sampa);
	waveforms->SetBranchAddress("Chn", Chn);
	waveforms->SetBranchAddress("Pad", Pad);
	waveforms->SetBranchAddress("gl_Chn", gl_Chn);
	waveforms->SetBranchAddress("Waveform", Waveform);
	waveforms->SetBranchAddress("x", x);
	waveforms->SetBranchAddress("y", y);
	// waveforms->SetBranchAddress("size", );

	// cout << waveforms->GetEntries() << endl;
	for (int i = 0; i < waveforms->GetEntries(); ++i)
	{
		waveforms->GetEntry(i);
		// cout << triggerID << "\t" << BX_counter << "\t" << size << endl;
		// for (int j = 0; j < size; ++j)
		// {
		// 	cout << "\t" << Sampa[j] << "\t" << Chn[j] << "\t" << Pad[j] << "\t" << gl_Chn[j] << "\t" << x[j] << "\t" << y[j];
		// 	for (int k = 0; k < nSamples; ++k)
		// 	{
		// 		cout << "\t" << Waveform[j][k];
		// 	}
		// 	cout << endl;
		// }
	}

	vector<Hit> extractedHits;

	waveforms->GetEntry(eventID);
	TMultiGraph* showEvent = VisWaveforms(Waveform,size);
	showEvent->Draw("a pmc plc l3d");
	ExtractHits(Waveform,gl_Chn,x,y,size,extractedHits);
	PrintHits(extractedHits);
	DrawHits(extractedHits);
	AnalyzeHits(extractedHits);

	return 0;
}

Double_t Fit(Int_t size,Double_t* x,Double_t* y,Float_t* weight,Double_t &a,Double_t &b)
{
  Double_t integ=0;
  Double_t p=0,q=0,r=0,s=0;
  Double_t resid=0;
  Double_t* w = new Double_t[size];

  for(int i = 0; i<size;i++){  //normalization of weights
    integ += weight[i];
  }
  for(int i = 0; i<size;i++){
    w[i] = weight[i]/integ;

    p += w[i]*x[i]*x[i];
    q += w[i]*x[i];
    r += w[i]*y[i];
    s += w[i]*x[i]*y[i];

  }
  if(q*q-p){ // fit results are in vertical line?
    a = (q*r-s)/(q*q-p);
    b = (q*s-p*r)/(q*q-p);
    for(int i = 0; i<size;i++)
      resid += w[i]*(y[i]-a*x[i]-b)*(y[i]-a*x[i]-b);
  }
  else return -1;
  delete[] w;
  return resid;

}
