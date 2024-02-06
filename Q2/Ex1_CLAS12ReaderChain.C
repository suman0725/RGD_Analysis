#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"

using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

void Ex1_CLAS12ReaderChain() {
    auto db = TDatabasePDG::Instance();
    TLorentzVector beam(0, 0, 10.6, 10.6);
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass());
    TLorentzVector el; 
    auto* hQ2 = new TH1F("Q2", "Q^2 distribution", 200, 0, 10);
    clas12root::HipoChain chain;
    // Add your HIPo files
    chain.Add("/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/018451/rec_clas_018451.evio.01215-01219.hipo");
    chain.db()->turnOffQADB();
    //         // chain.Add("/path/to/your/data/file2.hipo");
    //             // etc.
    auto config_c12 = chain.GetC12Reader();
    // Configure your reader as needed
    auto& c12 = chain.C12ref();

    while (chain.Next()) {
        // Event processing code here
        auto electrons = c12->getByID(11); // Get electrons
        if (electrons.size() == 1) { // Assuming one electron per event for simplicity
            SetLorentzVector(el, electrons[0]);
            // Calculate Q^2
            TLorentzVector q = beam - el; // 4-momentum transfer
            double Q2 = -q.Mag2();
            hQ2->Fill(Q2);
         }
  
   }
   // Post-processing: Save histograms, etc.
   TCanvas* can = new TCanvas("can", "Q^2 Distribution", 800, 600);
   can->cd(1);
   hQ2->Draw();
 // Other operations like saving the histogram to a file can go here
  }

