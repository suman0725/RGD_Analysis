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

void vertex5fit() {
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
        fitFunc[i] = new TF1(Form("fitFunc_S%d", i+1), "gaus", -15.0, 15.0); // Use a single Gaussian peak for fitting
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

    // Process data as before
    for(int i = 0; i < 6; ++i) {
        can1->cd(i+1);
    // Customize initial parameters based on sector
    switch(i) {
            case 0: // Sector 1
                fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), 10.0, 1.5);
                break;
	    case 1: // Sector 2
                fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), -1.0, 1.5);
                break;
            case 2: // Sector 3
                fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), 0.0, 2.0);
                break;
            case 3: // Sector 4
                fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), 1.0, 0.5);
                break;
            case 4: // Sector 5
                fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), 2.0, 1.2);
                break;
            case 5: // Sector 6
                fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), 3.0, 0.8);
                break;
            default:
               // Default parameters, just in case
               fitFunc[i]->SetParameters(hvz[i]->GetMaximum(), 0.0, 1.0);
                break;
        }

        hvz[i]->Draw();
        hvz[i]->Fit(Form("fitFunc_S%d", i+1), "R");
        fitFunc[i]->Draw("same");
    }
    can1->SaveAs("SectorsVzDistributionSinglePeak.pdf");
}                                                                                                                  
