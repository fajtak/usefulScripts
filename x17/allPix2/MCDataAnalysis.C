#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraph2D.h>

#include <memory>

#include "/home/fajtak/work/allpix-squared/src/objects/MCParticle.hpp"
#include "/home/fajtak/work/allpix-squared/src/objects/PixelCharge.hpp"
#include "/home/fajtak/work/allpix-squared/src/objects/PixelHit.hpp"
#include "/home/fajtak/work/allpix-squared/src/objects/PropagatedCharge.hpp"
#include "/home/fajtak/work/allpix-squared/src/objects/DepositedCharge.hpp"

void MCDataAnalysis(std::string filePath, std::string detector,int eventID = 0) {
	TFile* file = new TFile(filePath.c_str());
	TTree* depositedCharge = static_cast<TTree*>(file->Get("DepositedCharge"));
	if (!depositedCharge)
	{
		std::cout << "Could not read tree DepositedCharge, cannot continue." << std::endl;
		return;
	}
	TBranch* depositedChargeBranch = depositedCharge->FindBranch(detector.c_str());
	if (!depositedChargeBranch)
	{
		std::cout << "Could not find the branch on tree DepositedCharge for the corresponding detector, cannot continue"
                  << std::endl;
        return;
	}
	// Bind the information to a predefined vector
    std::vector<allpix::DepositedCharge*> depCharge;
    depositedChargeBranch->SetObject(&depCharge);

    TTree* propagatedCharge = static_cast<TTree*>(file->Get("PropagatedCharge"));
	if (!propagatedCharge)
	{
		std::cout << "Could not read tree propagatedCharge, cannot continue." << std::endl;
		return;
	}
	TBranch* propagatedChargeBranch = propagatedCharge->FindBranch(detector.c_str());
	if (!propagatedChargeBranch)
	{
		std::cout << "Could not find the branch on tree DepositedCharge for the corresponding detector, cannot continue"
                  << std::endl;
        return;
	}
	// Bind the information to a predefined vector
    std::vector<allpix::PropagatedCharge*> propCharge;
    propagatedChargeBranch->SetObject(&propCharge);

	// Initialise reading of the PixelHit TTrees
    TTree* pixel_hit_tree = static_cast<TTree*>(file->Get("PixelHit"));
    if(!pixel_hit_tree) {
        std::cout << "Could not read tree PixelHit, cannot continue." << std::endl;
        return;
    }
    TBranch* pixel_hit_branch = pixel_hit_tree->FindBranch(detector.c_str());
    if(!pixel_hit_branch) {
        std::cout << "Could not find the branch on tree PixelHit for the corresponding detector, cannot continue."
                  << std::endl;
        return;
    }
    // Bind the information to a predefined vector
    std::vector<allpix::PixelHit*> input_hits;
    pixel_hit_branch->SetObject(&input_hits);

    TGraph2D* g_depChargePos = new TGraph2D();
    TGraph2D* g_propChargePos = new TGraph2D();
    TGraph2D* g_propChargeTime = new TGraph2D();
    TH2D* h_sensor = new TH2D("h_sensor","SENSOR; x [pixels]; y [pixels]",256,0,256,256,0,256);


    // for(int i = 0; i < depositedCharge->GetEntries(); ++i)
    for(int i = eventID; i < eventID+1; ++i)
    {
        if(i % 100 == 0)
            std::cout << "Processing event " << i << std::endl;

        // Access next event. Pushes information into input_*
        pixel_hit_tree->GetEntry(i);
        depositedCharge->GetEntry(i);
        propagatedCharge->GetEntry(i);
        int nPointsDep = 0;
        int nPointsProp = 0;

        for (unsigned int j = 0; j < depCharge.size(); ++j)
        {
        	if (depCharge[j]->getSign() < 0 )
        		continue;

        	g_depChargePos->SetPoint(nPointsDep,depCharge[j]->getLocalPosition().X()/0.055,depCharge[j]->getLocalPosition().Y()/0.055,depCharge[j]->getLocalPosition().Z());
        	nPointsDep++;
        	std::cout << j << "\t" << depCharge[j]->getLocalTime() << "\t" << depCharge[j]->getCharge() << "\t" << depCharge[j]->getSign() << "\t" << depCharge[j]->getLocalPosition().X() << "\t" << depCharge[j]->getLocalPosition().Y() << "\t" << depCharge[j]->getLocalPosition().Z() << std::endl;
        }
        for (unsigned int j = 0; j < propCharge.size(); ++j)
        {
        	if (propCharge[j]->getSign() < 0 )
        		continue;

        	g_propChargePos->SetPoint(nPointsProp,propCharge[j]->getLocalPosition().X()/0.055,propCharge[j]->getLocalPosition().Y()/0.055,propCharge[j]->getLocalPosition().Z()+(j==0?0.000001:0));
        	g_propChargeTime->SetPoint(nPointsProp,propCharge[j]->getLocalPosition().X()/0.055,propCharge[j]->getLocalPosition().Y()/0.055,propCharge[j]->getLocalTime());
        	nPointsProp++;
        	std::cout << j << "\t" << propCharge[j]->getLocalTime() << "\t" << propCharge[j]->getCharge() << "\t" << propCharge[j]->getSign() << "\t" << propCharge[j]->getLocalPosition().X() << "\t" << propCharge[j]->getLocalPosition().Y() << "\t" << propCharge[j]->getLocalPosition().Z() << std::endl;
        }
        for (unsigned int j = 0; j < input_hits.size(); j++)
        {
        	h_sensor->Fill((input_hits[j]->getIndex().X()),(input_hits[j]->getIndex().Y()),input_hits[j]->getLocalTime());
        }

        TCanvas* c_depChargePos = new TCanvas("c_depChargePos","DepChargePos",800,600);
        g_depChargePos->Draw("P");
        g_depChargePos->SetTitle("; x [pixels]; y [pixels]; z [mm]");
        TCanvas* c_propChargePos = new TCanvas("c_propChargePos","PropChargePos",800,600);
        g_propChargePos->Draw("P");
        g_propChargePos->SetTitle("; x [pixels]; y [pixels]; z [mm]");
        TCanvas* c_propChargeTime = new TCanvas("c_propChargeTime","PropChargeTime",800,600);
        g_propChargeTime->Draw("P");
        g_propChargeTime->SetTitle("; x [pixels]; y [pixels]; Local time [ns]");
        TCanvas* c_sensor = new TCanvas("c_sensor","Sensor",800,600);
        h_sensor->Draw("colz");
        h_sensor->SetTitle("; x [pixels]; y [pixels]; Local time [ns]");

    }
}

void old(TFile* file, std::string detector)
{
    // Initialise reading of the PixelHit TTrees
    TTree* pixel_hit_tree = static_cast<TTree*>(file->Get("PixelHit"));
    if(!pixel_hit_tree) {
        std::cout << "Could not read tree PixelHit, cannot continue." << std::endl;
        return;
    }
    TBranch* pixel_hit_branch = pixel_hit_tree->FindBranch(detector.c_str());
    if(!pixel_hit_branch) {
        std::cout << "Could not find the branch on tree PixelHit for the corresponding detector, cannot continue."
                  << std::endl;
        return;
    }
    // Bind the information to a predefined vector
    std::vector<allpix::PixelHit*> input_hits;
    pixel_hit_branch->SetObject(&input_hits);

    // Initialise reading of the MCParticle TTrees
    TTree* mc_particle_tree = static_cast<TTree*>(file->Get("MCParticle"));
    if(!mc_particle_tree) {
        std::cout << "Could not read tree MCParticle" << std::endl;
        return;
    }
    TBranch* mc_particle_branch = mc_particle_tree->FindBranch(detector.c_str());
    if(!mc_particle_branch) {
        std::cout << "Could not find the branch on tree MCParticle for the corresponding detector, cannot continue"
                  << std::endl;
        return;
    }
    // Bind the information to a predefined vector
    std::vector<allpix::MCParticle*> input_particles;
    mc_particle_branch->SetObject(&input_particles);

    // Initialise histograms
    // Hitmap with arbitrary field of view
    TH2D* hitmap = new TH2D("hitmap", "Hitmap; x [mm]; y [mm]; hits", 200, 0, 20, 200, 0, 20);

    // Residuals:
    // first for all hits with the mean position of all MCParticles,
    // then only using the MCParticles that are part of the PixelHit history
    TH1D* residual_x = new TH1D("residual_x", "residual x; x_{MC} - x_{hit} [mm]; hits", 200, -5, 5);
    TH1D* residual_x_related =
        new TH1D("residual_x_related", "residual X, related hits; x_{MC} - x_{hit} [mm]; hits", 200, -5, 5);
    TH1D* residual_y = new TH1D("residual_y", "residual y; y_{MC} - y_{hit} [mm]; hits", 200, -5, 5);
    TH1D* residual_y_related =
        new TH1D("residual_y_related", "residual Y, related hits; y_{MC} - y_{hit} [mm]; hits", 200, -5, 5);

    // Spectrum of the PixelHit signal
    TH1D* spectrum = new TH1D("spectrum", "PixelHit signal spectrum; signal; hits", 200, 0, 100000);

    // Iterate over all events
    for(int i = 0; i < pixel_hit_tree->GetEntries(); ++i) {
        if(i % 100 == 0) {
            std::cout << "Processing event " << i << std::endl;
        }

        // Access next event. Pushes information into input_*
        pixel_hit_tree->GetEntry(i);
        mc_particle_tree->GetEntry(i);

        // Calculate the mean position of all MCParticles in this event
        double position_mcparts_x = 0.;
        double position_mcparts_y = 0.;
        for(auto& mc_part : input_particles) {
            position_mcparts_x += mc_part->getLocalReferencePoint().x();
            position_mcparts_y += mc_part->getLocalReferencePoint().y();
        }
        position_mcparts_x /= input_particles.size();
        position_mcparts_y /= input_particles.size();

        // Iterate over all PixelHits
        for(auto& hit : input_hits) {
            // Retrieve information using Allpix Squared methods
            double position_hit_x = hit->getPixel().getLocalCenter().x();
            double position_hit_y = hit->getPixel().getLocalCenter().y();
            double charge = hit->getSignal();

            // Access history of the PixelHit
            auto parts = hit->getMCParticles();

            // Calculate the mean position of all MCParticles related to this hit
            double position_mcparts_related_x = 0.;
            double position_mcparts_related_y = 0.;
            for(auto& part : parts) {
                position_mcparts_related_x += part->getLocalReferencePoint().x();
                position_mcparts_related_y += part->getLocalReferencePoint().y();
            }
            position_mcparts_related_x /= parts.size();
            position_mcparts_related_y /= parts.size();

            // Fill histograms
            hitmap->Fill(position_hit_x, position_hit_y);

            residual_x->Fill(position_mcparts_x - position_hit_x);
            residual_x_related->Fill(position_mcparts_related_x - position_hit_x);
            residual_y->Fill(position_mcparts_y - position_hit_y);
            residual_y_related->Fill(position_mcparts_related_y - position_hit_y);

            spectrum->Fill(charge);
        }
    }

    // Draw histograms
    TCanvas* c0 = new TCanvas("c0", "Hitmap", 600, 400);
    hitmap->Draw("colz");

    TCanvas* c1 = new TCanvas("c1", "Residuals", 1200, 800);
    c1->Divide(2, 2);

    c1->cd(1);
    residual_x->Draw();
    c1->cd(2);
    residual_y->Draw();

    c1->cd(3);
    residual_x_related->Draw();
    c1->cd(4);
    residual_y_related->Draw();

    TCanvas* c2 = new TCanvas("c2", "Signal spectrum", 600, 400);
    spectrum->Draw();
}