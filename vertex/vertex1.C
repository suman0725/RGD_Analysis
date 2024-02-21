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
#include "clas12reader.h"
#include "HipoChain.h"
#include "region_particle.h"
#include "runconfig.h"
#include "region_fdet.h"

using namespace clas12;

void vertex1() {

    clas12root::HipoChain chain; 
    chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ob_CuSn/dst/recon/018573/rec_clas_018573.evio.00045-00049.hipo");
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref(); 

    TCanvas* can1 = new TCanvas("can1", "Canvas1", 800, 600); 
    TH1F* hvz[6];
    for (int i = 0; i < 6; ++i) {
        hvz[i] = new TH1F(Form("hvz_S%d", i+1), Form("Vz for electrons (S%d); Vz (cm); counts", i+1), 200, -15.0, 15.0);
    }

    while (chain.Next()){
        auto electrons = c12->getByID(11); 

        for (const auto& elec : electrons) {
            auto phie = elec->getPhi() * TMath::RadToDeg(); // Phi in degrees
            auto vz = elec->par()->getVz();
            cout << phie << endl;   
            int sector;
            if (phie >= -180 && phie < -120) sector = 0;
            else if (phie >= -120 && phie < -60) sector = 1;
            else if (phie >= -60 && phie < 0) sector = 2;
            else if (phie >= 0 && phie < 60) sector = 3;
            else if (phie >= 60 && phie < 120) sector = 4;
            else if (phie >= 120 && phie < 180) sector = 5;

            hvz[sector]->Fill(vz);
        }
    }

    for (int i = 0; i < 6; ++i) {
        can1->cd(i+1);
        hvz[i]->Draw();
    }
}
