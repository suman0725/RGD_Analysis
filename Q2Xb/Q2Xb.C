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
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"

using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

void Q2Xb() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.6; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Assuming the beam is along the z-axis
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el; 
    auto* hQ2 = new TH1F("Q2", "Q^2 distribution", 200, 0, 10);
    auto* hxB = new TH1F("xB", "Bjorken x distribution", 100, 0, 1); // Histogram for Bjorken x
    auto* hQ2vsxB = new TH2F("Q2vsxB", "Q^2 vs Bjorken x;Bjorken x;Q^2", 100, 0, 1, 200, 0, 10); // 2D histogram for Q^2 vs x_B
    clas12root::HipoChain chain;
    // Add your HIPo files
    chain.Add("/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/018451/rec_clas_018451.evio.01215-01219.hipo");
    chain.db()->turnOffQADB();
    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();

    while (chain.Next()) {
        auto electrons = c12->getByID(11); // Get electrons
         if (electrons.size() == 1) { // Assuming one electron per event for simplicity
            SetLorentzVector(el, electrons[0]);
            TLorentzVector q = beam - el; // 4-momentum transfer
            double Q2 = -q.Mag2();
            double nu = beam.Energy() - el.Energy();
            double xB = Q2 / (2 * target.M() * nu); // Calculate Bjorken x
            hQ2->Fill(Q2);
            hxB->Fill(xB);
            hQ2vsxB->Fill(xB, Q2); // Fill 2D histogram
         }
    }

// Post-processing: Draw histograms
   TCanvas* can = new TCanvas("can", "Distributions", 1200, 400);
    can->Divide(3,1);
    can->cd(1);
    hQ2->Draw();
    can->cd(2);
    hxB->Draw();
    can->cd(3);
    hQ2vsxB->Draw("COLZ"); // Draw 2D histogram
// Optionally, save the histograms to a file
//     // TFile outFile("results.root", "RECREATE");
//         // hQ2->Write();
//             // hxB->Write();
//                 // hQ2vsxB->Write();
//                     // outFile.Close();
                                    }
