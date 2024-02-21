#include <TFile.h>
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
#include <TF1.h>
#include "clas12reader.h"
#include "HipoChain.h"
#include "region_particle.h"
#include "runconfig.h"
#include "region_fdet.h"

using namespace clas12;

void vertex2fit() {
    clas12root::HipoChain chain; 
    chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ob_CuSn/dst/recon/018573/rec_clas_018573.evio.00045-00049.hipo");
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref(); 

    TCanvas* can1 = new TCanvas("can1", "Canvas for 6 Sectors", 1200, 800);
    can1->Divide(3, 2); // Arrange the canvas to fit 6 histograms

    TH1F* hvz[6];
    TF1* fitFunc[6];
    for(int i = 0; i < 6; ++i) {
        hvz[i] = new TH1F(Form("hvz_S%d", i+1), Form("Vz for electrons (Sector %d); Vz (cm); counts", i+1), 200, -15.0, 15.0);
        fitFunc[i] = new TF1(Form("fitFunc_S%d", i+1), "gaus(0)+gaus(3)", -15.0, 15.0); // Two Gaussian peaks
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
   
    for(int i = 0; i < 6; ++i) {
        can1->cd(i+1); // Move to the next pad
        hvz[i]->Draw();

       fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), 0.0, 1.0, hvz[i]->GetMaximum(), 0.0, 1.0);
      
       hvz[i]->Fit(Form("fitFunc_S%d", i+1), "R");

      fitFunc[i]->Draw("same");
    }
 can1->SaveAs("SectorsVzDistribution.pdf");
}
