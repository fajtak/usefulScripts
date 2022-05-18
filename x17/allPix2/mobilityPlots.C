#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"

// Double_t myfunction(Double_t *x, Double_t *par)
// {
//    Float_t xx =x[0];
//    Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
//    return f;
// }
// void myfunc()
// {
//    TF1 *f1 = new TF1("myfunc",myfunction,0,10,2);
//    f1->SetParameters(2,1);
//    f1->SetParNames("constant","coefficient");
//    f1->Draw();
// }
// void myfit()
// {
//    TH1F *h1=new TH1F("h1","test",100,0,10);
//    h1->FillRandom("myfunc",20000);
//    TF1 *f1 = (TF1 *)gROOT->GetFunction("myfunc");
//    f1->SetParameters(800,1);
//    h1->Fit("myfunc");
// }

Double_t jacoboniElectrons(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t vme = 1.53e9 * TMath::Power(T,-0.87);
	Double_t Ece = 1.01 * TMath::Power(T,1.55);
	Double_t Be = 2.57e-2 * TMath::Power(T,0.66);

	return (vme/Ece)*1/TMath::Power(1+TMath::Power(x[0]/Ece,Be),1/Be);
}

Double_t jacoboniHoles(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t vmh = 1.62e8 * TMath::Power(T,-0.52);
	Double_t Ech = 1.24 * TMath::Power(T,1.68);
	Double_t Bh = 0.46 * TMath::Power(T,0.17);

	return (vmh/Ech)*1/TMath::Power(1+TMath::Power(x[0]/Ech,Bh),1/Bh);
}

Double_t jacoboniHolesTemperature(Double_t* x, Double_t* par)
{
	Double_t T = x[0];
	Double_t vmh = 1.62e8 * TMath::Power(T,-0.52);
	Double_t Ech = 1.24 * TMath::Power(T,1.68);
	Double_t Bh = 0.46 * TMath::Power(T,0.17);

	return (vmh/Ech)*1/TMath::Power(1+TMath::Power(par[0]/Ech,Bh),1/Bh);
}

Double_t jacoboniPropagationTimeElectrons(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu = jacoboniElectrons(x,par);
	Double_t sensorThickness = 0.05; //cm
	Double_t velocity = mu*x[0];
	return sensorThickness/velocity*1e9;
}

Double_t jacoboniPropagationTimeHoles(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu = jacoboniHoles(x,par);
	Double_t sensorThickness = 0.05; //cm
	Double_t velocity = mu*x[0];
	return sensorThickness/velocity*1e9;
}

Double_t jacoboniPropagationTimeHolesTrue(Double_t* x, Double_t* par)
{
	// Double_t T = par[0];
	Double_t UD = par[1];
	Double_t UB = x[0]/20;
	Double_t mu = jacoboniHoles(x,par);
	// cout << mu << endl;
	Double_t sensorThickness = 0.05; //cm
	// Double_t velocity = mu*x[0];
	return -TMath::Power(sensorThickness,2)/(2*UD*mu)*TMath::Log(1 - (2*UD)/(UD+UB))*1e9;
}

Double_t hamburgElectrons(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu0e = 1530 * TMath::Power(T/300,-2.42);
	Double_t vsat = 1.03e7 * TMath::Power(T/300,-0.226);
	Double_t mueInv = 1/mu0e + x[0]/vsat;
	return 1/mueInv;
}

Double_t hamburgHoles(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu0h = 464 * TMath::Power(T/300,-2.2);
	Double_t b = 9.57e-8 * TMath::Power(T/300,-0.101);
	Double_t c = -3.31e-13;
	Double_t E0 = 2640 * TMath::Power(T/300,0.526);

	Double_t muhInv = 0;
	if (x[0] < E0)
		muhInv = 1/mu0h;
	else
		muhInv = 1/mu0h + b * (x[0] - E0) + c * TMath::Power((x[0] - E0),2);
	return 1/muhInv;
}

Double_t hamburgPropagationTimeHoles(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu = hamburgHoles(x,par);
	Double_t sensorThickness = 0.05; //cm
	Double_t velocity = mu*x[0];
	return sensorThickness/velocity*1e9;
}

Double_t hamburgHFElectrons(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu0e = 1430 * TMath::Power(T/300,-1.99);
	Double_t vsat = 1.05e7 * TMath::Power(T/300,-0.302);
	Double_t mueInv = 1/mu0e + x[0]/vsat;
	return 1/mueInv;
}

Double_t hamburgHFHoles(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu0h = 457 * TMath::Power(T/300,-2.8);
	Double_t b = 9.57e-8 * TMath::Power(T/300,-0.155);
	Double_t c = -3.24e-13;
	Double_t E0 = 2970 * TMath::Power(T/300,0.563);

	Double_t muhInv = 0;
	if (x[0] < E0)
		muhInv = 1/mu0h;
	else
		muhInv = 1/mu0h + b * (x[0] - E0) + c * TMath::Power((x[0] - E0),2);
	return 1/muhInv;
}

Double_t hamburgHFPropagationTimeHoles(Double_t* x, Double_t* par)
{
	Double_t T = par[0];
	Double_t mu = hamburgHFHoles(x,par);
	Double_t sensorThickness = 0.05; //cm
	Double_t velocity = mu*x[0];
	return sensorThickness/velocity*1e9;
}

int mobilityPlots()
{
	TCanvas* c_mobilityPlots = new TCanvas("c_mobilityPlots","MobilityPlots",800,600);
	TF1* f_jacoboniElectrons = new TF1("Jacoboni, Electrons, T = 352K",jacoboniElectrons,0,5000,1);
	f_jacoboniElectrons->SetParameter(0,352.0);
	f_jacoboniElectrons->Draw();
	f_jacoboniElectrons->SetLineColor(kBlue);
	f_jacoboniElectrons->GetYaxis()->SetRangeUser(0,2000);
	f_jacoboniElectrons->GetYaxis()->SetTitle("#mu [cm^{2}/V/s]");
	f_jacoboniElectrons->GetXaxis()->SetTitle("E [V/cm]");

	TF1* f_jacoboniElectrons_293 = new TF1("Jacoboni, Electrons, T = 293K",jacoboniElectrons,0,5000,1);
	f_jacoboniElectrons_293->SetParameter(0,293.0);
	f_jacoboniElectrons_293->Draw("same");
	f_jacoboniElectrons_293->SetLineColor(kBlue);
	f_jacoboniElectrons_293->SetLineStyle(6);

	TF1* f_jacoboniHoles = new TF1("Jacoboni, Holes, T = 352K",jacoboniHoles,0,5000,1);
	f_jacoboniHoles->SetParameter(0,352.0);
	f_jacoboniHoles->Draw("same");
	f_jacoboniHoles->SetLineColor(kRed);

	TF1* f_jacoboniHoles_293 = new TF1("Jacoboni, Holes, T = 293K",jacoboniHoles,0,5000,1);
	f_jacoboniHoles_293->SetParameter(0,293.0);
	f_jacoboniHoles_293->Draw("same");
	f_jacoboniHoles_293->SetLineColor(kRed);
	f_jacoboniHoles_293->SetLineStyle(6);

	TF1* f_hamburgElectrons = new TF1("Hamburg, Electrons, T = 352K",hamburgElectrons,0,5000,1);
	f_hamburgElectrons->SetParameter(0,352.0);
	f_hamburgElectrons->Draw("same");
	f_hamburgElectrons->SetLineColor(kGreen);
	// f_hamburgElectrons->SetLineStyle(6);

	TF1* f_hamburgHoles = new TF1("Hamburg, Holes, T = 352K",hamburgHoles,0,5000,1);
	f_hamburgHoles->SetParameter(0,352.0);
	f_hamburgHoles->Draw("same");
	f_hamburgHoles->SetLineColor(kYellow);
	// f_hamburgHoles->SetLineStyle(6);

	TF1* f_hamburgHFElectrons = new TF1("Hamburg HF, Electrons, T = 352K",hamburgHFElectrons,0,5000,1);
	f_hamburgHFElectrons->SetParameter(0,352.0);
	f_hamburgHFElectrons->Draw("same");
	f_hamburgHFElectrons->SetLineColor(kMagenta);
	// f_hamburgHFElectrons->SetLineStyle(6);

	TF1* f_hamburgHFHoles = new TF1("Hamburg HF, Holes, T = 352K",hamburgHFHoles,0,5000,1);
	f_hamburgHFHoles->SetParameter(0,352.0);
	f_hamburgHFHoles->Draw("same");
	f_hamburgHFHoles->SetLineColor(kBlack);
	// f_hamburgHFHoles->SetLineStyle(6);

	gPad->BuildLegend();

	TCanvas* c_propagationTimes = new TCanvas("c_propagationTimes","PropagationTimes",800,600);
	TF1* f_jacoboniPropagationTimeHoles = new TF1("Jacoboni, Holes, T = 352K",jacoboniPropagationTimeHoles,3000,5000,1);
	f_jacoboniPropagationTimeHoles->SetParameter(0,352.0);
	f_jacoboniPropagationTimeHoles->Draw();
	f_jacoboniPropagationTimeHoles->SetLineColor(kRed);
	f_jacoboniPropagationTimeHoles->GetYaxis()->SetRangeUser(0,100);
	f_jacoboniPropagationTimeHoles->GetYaxis()->SetTitle("#DeltaT (500 #mum) [ns]");
	f_jacoboniPropagationTimeHoles->GetXaxis()->SetTitle("E [V/cm]");

	TF1* f_jacoboniPropagationTimeHolesTrue_Ud100 = new TF1("Jacoboni, Holes, T = 352K, E linear, Ud = 100",jacoboniPropagationTimeHolesTrue,3000,5000,2);
	f_jacoboniPropagationTimeHolesTrue_Ud100->SetParameter(0,352.0); // Temperature in Kelvins
	f_jacoboniPropagationTimeHolesTrue_Ud100->SetParameter(1,100.0); // Depletion voltage
	f_jacoboniPropagationTimeHolesTrue_Ud100->Draw("same");
	f_jacoboniPropagationTimeHolesTrue_Ud100->SetLineColor(kRed);
	f_jacoboniPropagationTimeHolesTrue_Ud100->SetLineStyle(7);

	TF1* f_jacoboniPropagationTimeHolesTrue_Ud80 = new TF1("Jacoboni, Holes, T = 352K, E linear, Ud = 80",jacoboniPropagationTimeHolesTrue,3000,5000,2);
	f_jacoboniPropagationTimeHolesTrue_Ud80->SetParameter(0,352.0); // Temperature in Kelvins
	f_jacoboniPropagationTimeHolesTrue_Ud80->SetParameter(1,80.0); // Depletion voltage
	f_jacoboniPropagationTimeHolesTrue_Ud80->Draw("same");
	f_jacoboniPropagationTimeHolesTrue_Ud80->SetLineColor(kRed);
	f_jacoboniPropagationTimeHolesTrue_Ud80->SetLineStyle(8);

	TF1* f_jacoboniPropagationTimeHoles_330 = new TF1("Jacoboni, Holes, T = 330K",jacoboniPropagationTimeHoles,3000,5000,1);
	f_jacoboniPropagationTimeHoles_330->SetParameter(0,330.0);
	f_jacoboniPropagationTimeHoles_330->Draw("same");
	f_jacoboniPropagationTimeHoles_330->SetLineColor(kRed);
	f_jacoboniPropagationTimeHoles_330->SetLineStyle(9);

	TF1* f_jacoboniPropagationTimeHoles_293 = new TF1("Jacoboni, Holes, T = 293K",jacoboniPropagationTimeHoles,3000,5000,1);
	f_jacoboniPropagationTimeHoles_293->SetParameter(0,293.0);
	f_jacoboniPropagationTimeHoles_293->Draw("same");
	f_jacoboniPropagationTimeHoles_293->SetLineColor(kRed);
	f_jacoboniPropagationTimeHoles_293->SetLineStyle(6);

	TF1* f_hamburgPropagationTimeHoles = new TF1("Hamburg, Holes, T = 352K",hamburgPropagationTimeHoles,3000,5000,1);
	f_hamburgPropagationTimeHoles->SetParameter(0,352.0);
	// f_hamburgPropagationTimeHoles->Draw("same");
	f_hamburgPropagationTimeHoles->SetLineColor(kYellow);
	// f_hamburgPropagationTimeHoles->SetLineStyle(9);

	TF1* f_hamburgHFPropagationTimeHoles = new TF1("Hamburg HF, Holes, T = 352K",hamburgHFPropagationTimeHoles,3000,5000,1);
	f_hamburgHFPropagationTimeHoles->SetParameter(0,352.0);
	// f_hamburgHFPropagationTimeHoles->Draw("same");
	f_hamburgHFPropagationTimeHoles->SetLineColor(kBlack);
	// f_hamburgHFPropagationTimeHoles->SetLineStyle(9);

	TF1* f_jacoboniPropagationTimeElectrons = new TF1("Jacoboni, Electrons, T = 352K",jacoboniPropagationTimeElectrons,3000,5000,1);
	f_jacoboniPropagationTimeElectrons->SetParameter(0,352.0);
	f_jacoboniPropagationTimeElectrons->Draw("same");
	f_jacoboniPropagationTimeElectrons->SetLineColor(kBlue);
	// f_jacoboniPropagationTimeElectrons->SetLineStyle(9);

	TF1* f_jacoboniPropagationTimeElectrons_293 = new TF1("Jacoboni, Electrons, T = 293K",jacoboniPropagationTimeElectrons,3000,5000,1);
	f_jacoboniPropagationTimeElectrons_293->SetParameter(0,293.0);
	f_jacoboniPropagationTimeElectrons_293->Draw("same");
	f_jacoboniPropagationTimeElectrons_293->SetLineColor(kBlue);
	f_jacoboniPropagationTimeElectrons_293->SetLineStyle(6);

	gPad->BuildLegend();

	TCanvas* c_mobilityPlotsTemperature = new TCanvas("c_mobilityPlotsTemperature","MobilityPlotsTemperature",800,600);
	TF1* f_jacoboniHolesTemp = new TF1("Jacoboni, Holes, E = 4200V/cm",jacoboniHolesTemperature,250,500,1);
	f_jacoboniHolesTemp->SetParameter(0,4200);
	f_jacoboniHolesTemp->Draw();
	f_jacoboniHolesTemp->SetLineColor(kRed);
	// f_jacoboniHolesTemp->GetYaxis()->SetRangeUser(0,500);
	f_jacoboniHolesTemp->GetYaxis()->SetTitle("#mu [cm^{2}/V/s]");
	f_jacoboniHolesTemp->GetXaxis()->SetTitle("T [K]");

	gPad->BuildLegend();

	return 0;
}
