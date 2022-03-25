#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"

TH1D* h_distPosSphere = new TH1D("h_distPosSphere","Distance interaction point and origin for sphere; D [m]; NoE [#]",1000,0,1000);
TH1D* h_distPosCylinder = new TH1D("h_distPosCylinder","Distance interaction point and origin for cylinder; D [m]; NoE [#]",1000,0,1000);
TH1D* h_xPosSphere = new TH1D("h_xPosSphere","X position sphere; x [m]; NoE [#]",2000,-1000,1000);
TH1D* h_yPosSphere = new TH1D("h_yPosSphere","Y position sphere; y [m]; NoE [#]",2000,-1000,1000);
TH1D* h_zPosSphere = new TH1D("h_zPosSphere","Z position sphere; z [m]; NoE [#]",2000,-1000,1000);
TH1D* h_xPosCylinder = new TH1D("h_xPosCylinder","X position Cylinder; x [m]; NoE [#]",2000,-1000,1000);
TH1D* h_yPosCylinder = new TH1D("h_yPosCylinder","Y position Cylinder; y [m]; NoE [#]",2000,-1000,1000);
TH1D* h_zPosCylinder = new TH1D("h_zPosCylinder","Z position Cylinder; z [m]; NoE [#]",2000,-1000,1000);


int simulateVolume()
{
	TRandom3 *eventGenerator = new TRandom3();
	Double_t x,y,z;

	for (int i = 0; i < 10000000; ++i)
	{
		x = eventGenerator->Uniform(2*520)-520;
		y = eventGenerator->Uniform(2*520)-520;
		z = eventGenerator->Uniform(2*520)-520;
		TVector3 pos(x,y,z);
		if (pos.Mag() <= 520)
		{
			h_distPosSphere->Fill(pos.Mag());
			h_xPosSphere->Fill(x);
			h_yPosSphere->Fill(y);
			h_zPosSphere->Fill(z);
		}
	}

	for (int i = 0; i < 10000000; ++i)
	{
		x = eventGenerator->Uniform(2*320)-320;
		y = eventGenerator->Uniform(2*320)-320;
		z = eventGenerator->Uniform(2*410)-410;
		TVector3 pos(x,y,0);
		if (pos.Mag() <= 320)
		{
			pos.SetXYZ(x,y,z);
			h_distPosCylinder->Fill(pos.Mag());
			h_xPosCylinder->Fill(x);
			h_yPosCylinder->Fill(y);
			h_zPosCylinder->Fill(z);
		}
	}


	TCanvas* c_distPosSphere = new TCanvas("c_distPosSphere","DistancePositionSphere",800,600);
	h_distPosSphere->Draw();

	TCanvas* c_xPosSphere = new TCanvas("c_xPosSphere","XPositionSphere",800,600);
	h_xPosSphere->Draw();

	TCanvas* c_yPosSphere = new TCanvas("c_yPosSphere","YPositionSphere",800,600);
	h_yPosSphere->Draw();

	TCanvas* c_zPosSphere = new TCanvas("c_zPosSphere","ZPositionSphere",800,600);
	h_zPosSphere->Draw();

	TCanvas* c_xPosCylinder = new TCanvas("c_xPosCylinder","XPositionCylinder",800,600);
	h_xPosCylinder->Draw();

	TCanvas* c_yPosCylinder = new TCanvas("c_yPosCylinder","YPositionCylinder",800,600);
	h_yPosCylinder->Draw();

	TCanvas* c_zPosCylinder = new TCanvas("c_zPosCylinder","ZPositionCylinder",800,600);
	h_zPosCylinder->Draw();

	TCanvas* c_distPosCylinder = new TCanvas("c_distPosCylinder","DistancePositionCylinder",800,600);
	h_distPosCylinder->Draw();

	return 0;
}