#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"
#include "region_particle.h"
#include "runconfig.h"
#include "region_fdet.h"

using namespace clas12;

void vz() {
    clas12root::HipoChain chain; 
    chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ob_CuSn/dst/recon/018573/rec_clas_018573.evio.00045-00049.hipo");
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref(); 

    TCanvas* can1 = new TCanvas("can1", "Canvas for 6 Sectors", 1200, 800);
   can1->Divide(3, 2);
   TH1F* hvz[6];
    for(int i = 0; i < 6; ++i) {
        hvz[i] = new TH1F(Form("hvz_S%d", i+1), Form("Vz for electrons (Sector %d); Vz (cm); counts", i+1), 200, -30.0, 30.0);
    }
   while (chain.Next()) {
        auto electrons = c12->getByID(11); 

        for (const auto& elec : electrons) {
            auto phi = elec->getPhi() * TMath::RadToDeg(); 
    auto vz = elec->par()->getVz();
            int sector = -1;
if(phi >= -180 && phi < -120) sector = 0; // Sector 1
            else if(phi >= -120 && phi < -60) sector = 1; // Sector 2
            else if(phi >= -60 && phi < 0) sector = 2; // Sector 3
            else if(phi >= 0 && phi < 60) sector = 3; // Sector 4
            else if(phi >= 60 && phi < 120) sector = 4; // Sector 5
            else if(phi >= 120 && phi < 180) sector = 5; // Sector 6
           
           if(sector != -1) {
                hvz[sector]->Fill(vz);
            }
        }
    }

	double g1Mean[] = {-8.4, -7, -8, -8, -9, -8.5}; // Adjust sector 1 value as needed
    double g1Sigma[] = {8, 8, 8, 8, 8, 8}; // Adjust sector 1 value as needed
    double g1StartRange[] = {-9.5, -9, -9, -10.5, -10.2, -10}; // Adjust sector 1 value as needed
    double g1EndRange[] = {-6, -6, -6.5, -7, -7, -7}; // Adjust sector 1 value as needed
	double g1Height[] = {70, 70, 70, 70, 70, 70};


 	 double g2Mean[] = {-4.4, -3.5, -2.5, -3.5, -4, -4}; // Adjust sector 1 value as needed
    double g2Sigma[] = {8, 8, 8, 8, 8, 8}; // Adjust sector 1 value as needed
    double g2StartRange[] = {-4.5, -4, -4.5, -4.5, -5, -5}; // Adjust sector 1 value as needed
    double g2EndRange[] = {-1, -1, -1.5, -1.5, -2.5, -2}; // Adjust sector 1 value as needed
	double g2Height[] = {90, 90, 90, 90, 90, 90};



   for(int i = 0; i < 6; ++i) {
        can1->cd(i+1); // Move to the next pad

	hvz[i]->SetStats(0); 
        hvz[i]->Draw();
		
/*	if (i == 0) { // Fit only for the first sector's histogram as a demonstration
            // Define two Gaussian fits for the two peaks
            TF1 *gaus1 = new TF1("gaus1", "gaus", -9.5, -6); // Range around the first peak
            TF1 *gaus2 = new TF1("gaus2", "gaus", -4.5, -1); // Range around the second peak
	    // Set initial parameters for the Gaussians
	    gaus1->SetParameters(70, 8.46, 8); // height, mean, standard deviation
            gaus2->SetParameters(90, 4.4, 8); // Adjust height as neededi
	
	    // Fit both Gaussians to the histogram
	    hvz[i]->Fit(gaus1, "R+");
            hvz[i]->Fit(gaus2, "R+");


	can1->cd(1);
	hvz[0]->SetStats(0);
    hvz[0]->Draw();
    gaus1->SetLineColor(kRed);
    gaus2->SetLineColor(kBlue);
    gaus1->SetLineWidth(1);
    gaus2->SetLineWidth(1);
    gaus1->Draw("same");
    gaus2->Draw("same");
 

	 TLatex latex;

	    latex.SetTextSize(0.020);
    latex.SetTextAlign(13);
    latex.DrawLatexNDC(0.6, 0.85, Form("Cu: Mean = %.2f, Sigma = %.2f", gaus1->GetParameter(1), gaus1->GetParameter(2)));
    latex.DrawLatexNDC(0.6, 0.80, Form("Sn: Mean = %.2f, Sigma = %.2f", gaus2->GetParameter(1), gaus2->GetParameter(2)));

		*/


	 TF1 *gaus1 = new TF1(Form("gaus1_S%d", i+1), "gaus", g1StartRange[i], g1EndRange[i]);
        gaus1->SetParameters(g1Height[i], g1Mean[i], g1Sigma[i]);
        hvz[i]->Fit(gaus1, "R+");

        TF1 *gaus2 = new TF1(Form("gaus2_S%d", i+1), "gaus", g2StartRange[i], g2EndRange[i]);
        gaus2->SetParameters(g2Height[i], g2Mean[i], g2Sigma[i]); // Example: adjust maximum accordingly
        hvz[i]->Fit(gaus2, "R+");

        gaus1->SetLineColor(kRed);
        gaus2->SetLineColor(kBlue);
	
      	// gaus1->SetLineWidth(1);
   	 //gaus2->SetLineWidth(1);
        gaus1->Draw("same");
        gaus2->Draw("same");
	
	 // Drawing 3 sigma lines for Cu 
	 double line1_start = gaus1->GetParameter(1) - 3 * gaus1->GetParameter(2);
        double line1_end = gaus1->GetParameter(1) + 3 * gaus1->GetParameter(2);
        TLine *line1a = new TLine(line1_start, 0, line1_start, gaus1->GetMaximum());
        TLine *line1b = new TLine(line1_end, 0, line1_end, gaus1->GetMaximum());
        line1a->SetLineColor(kRed);
        line1b->SetLineColor(kRed);
        line1a->Draw("same");
        line1b->Draw("same");

	// Drawing 3 sigma lines for Sn (gaus2)
	double line2_start = gaus2->GetParameter(1) - 3 * gaus2->GetParameter(2);
        double line2_end = gaus2->GetParameter(1) + 3 * gaus2->GetParameter(2);
        TLine *line2a = new TLine(line2_start, 0, line2_start, gaus2->GetMaximum());
        TLine *line2b = new TLine(line2_end, 0, line2_end, gaus2->GetMaximum());
        line2a->SetLineColor(kBlue);
        line2b->SetLineColor(kBlue);
        line2a->Draw("same");
        line2b->Draw("same");

	
	

        TLatex latex;
        latex.SetTextSize(0.020);
        latex.SetTextAlign(13);
        latex.DrawLatexNDC(0.55, 0.85, Form("Cu(red): Mean = %.2f, Sigma = %.2f", gaus1->GetParameter(1), gaus1->GetParameter(2)));
        latex.DrawLatexNDC(0.55, 0.80, Form("Sn(blue): Mean = %.2f, Sigma = %.2f", gaus2->GetParameter(1), gaus2->GetParameter(2)));
	latex.DrawLatexNDC(0.55, 0.75, Form("Entries: %d", int(hvz[i]->GetEntries())));





	}
 	can1->Update();	

        can1->SaveAs("vzdistributionfit3sigmacut.pdf");
    }

