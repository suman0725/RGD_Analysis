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

void vertex() {

    clas12root::HipoChain chain; 
    chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ob_CuSn/dst/recon/018573/rec_clas_018573.evio.00045-00049.hipo");
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref(); 

    TCanvas* can1 = new TCanvas("can1", "Canvas1", 800, 600); 
    TH1F* hvz = new TH1F("hvz_S4", "Vz for electrons (S4); Vz (cm); counts", 200, -15.0, 15.0);
   // TH2F* hvzphie = new TH2F("hvzphie", "Vz for electrons (S4); Vz (cm); counts", 200, -15.0, 15.0);

    while (chain.Next()){
        auto electrons = c12->getByID(11); 

        for (const auto& elec : electrons) {
            auto phie = elec->getPhi() * TMath::RadToDeg();
    //        cout<<phie<<endl;
            auto vz = elec->par()->getVz();
            if (phie >= 0.0 && phie < 60.0) {
                hvz->Fill(vz);
            }

        }
    }

    can1->cd();
    hvz->Draw();
}

