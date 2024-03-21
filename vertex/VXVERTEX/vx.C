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

void vx() {
    clas12root::HipoChain chain; 
    chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ob_CuSn/dst/recon/018573/rec_clas_018573.evio.00045-00049.hipo");
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref(); 

    TCanvas* can1 = new TCanvas("can1", "Canvas for 6 Sectors", 1200, 800);
   can1->Divide(3, 2);
   TH1F* hvx[6];
    for(int i = 0; i < 6; ++i) {
        hvx[i] = new TH1F(Form("hvx_S%d", i+1), Form("Vx for electrons (Sector %d); Vx (cm); counts", i+1), 200, -30.0, 30.0);
    }
   while (chain.Next()) {
        auto electrons = c12->getByID(11); 

        for (const auto& elec : electrons) {
            auto phi = elec->getPhi() * TMath::RadToDeg(); 
    auto vx = elec->par()->getVx();
            int sector = -1;
if(phi >= -180 && phi < -120) sector = 0; // Sector 1
            else if(phi >= -120 && phi < -60) sector = 1; // Sector 2
            else if(phi >= -60 && phi < 0) sector = 2; // Sector 3
            else if(phi >= 0 && phi < 60) sector = 3; // Sector 4
            else if(phi >= 60 && phi < 120) sector = 4; // Sector 5
            else if(phi >= 120 && phi < 180) sector = 5; // Sector 6
           
           if(sector != -1) {
                hvx[sector]->Fill(vx);
            }
        }
    }
   for(int i = 0; i < 6; ++i) {
        can1->cd(i+1); // Move to the next pad
	hvx[i]->SetStats(0);
        hvx[i]->Draw();
	// Fit with a single Gaussian
	TF1 *gaus = new TF1(Form("gaus_S%d", i+1), "gaus", -0.5, 0.5); // Adjust the fitting range as needed
        hvx[i]->Fit(gaus, "R+");
        gaus->SetLineColor(kRed);
        gaus->Draw("same");
	// Optionally, draw sigma lines based on the Gaussian fit
	double line_start = gaus->GetParameter(1) - 3 * gaus->GetParameter(2);
        double line_end = gaus->GetParameter(1) + 3 * gaus->GetParameter(2);
        TLine *linea = new TLine(line_start, 0, line_start, gaus->GetMaximum());
        TLine *lineb = new TLine(line_end, 0, line_end, gaus->GetMaximum());
        linea->SetLineColor(kBlue);
        lineb->SetLineColor(kBlue);
        linea->Draw("same");
        lineb->Draw("same");

	TLatex latex;
        latex.SetTextSize(0.020);
        latex.SetTextAlign(13);
        latex.DrawLatexNDC(0.55, 0.85, Form("Mean = %.2f, Sigma = %.2f", gaus->GetParameter(1), gaus->GetParameter(2)));
        latex.DrawLatexNDC(0.55, 0.80, Form("Entries: %d", int(hvx[i]->GetEntries())));
     }
	can1->Update();

        can1->SaveAs("vxdistribution3sigmacut.pdf");
    }


