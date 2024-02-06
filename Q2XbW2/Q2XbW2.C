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

void Q2XbW2() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.6; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Assuming the beam is along the z-axis
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el;
    TLorentzVector pip;
    auto* hQ2 = new TH1F("Q2", "Q^2 distribution", 200, 0, 10);
    auto* hxB = new TH1F("xB", "Bjorken x distribution", 100, 0, 1);
    auto* hQ2vsxB = new TH2F("Q2vsxB", "Q^2 vs Bjorken x;Bjorken x;Q^2", 100, 0, 1, 200, 0, 10);
    auto* hW2 = new TH1F("W2", "W^2 distribution", 200, 0, 20); // Histogram for W^2
    auto* hW = new TH1F("W", "W distribution", 200, 0, 8); // Histogram for W
    auto* hmiss = new TH1F("Mx","Missing mass distribution",200, 0,10);
    auto* hPt2 = new TH1F("Pt2", "Transverse Momentum Squared of Pi+", 100, 0, 1); // Histogram for Pt^2

    //auto* hmiss2 = new TH1F("Mx2","missing mass squared distribution",200,0,2);  
    clas12root::HipoChain chain;
     // Add your HIPo files
     chain.Add("/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/018451/rec_clas_018451.evio.01215-01219.hipo");
    chain.db()->turnOffQADB();
    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();
    while (chain.Next()) {
        auto electrons = c12->getByID(11); // Get electrons
        auto pips = c12->getByID(211);
        if (electrons.size() == 1 && pips.size() ==1) { // Assuming one electron per event for simplicity
            SetLorentzVector(el, electrons[0]);
            SetLorentzVector(pip, pips[0]); 
            TLorentzVector q = beam - el; // 4-momentum transfer
            // Calculate the missing system assuming only an electron and a pi+ are detected
            TLorentzVector miss = beam + target - el - pip; // Adjust according to detected particles
            
            double Q2 = -q.Mag2();
            double nu = beam.Energy() - el.Energy();
            double xB = Q2 / (2 * target.M() * nu); // Calculate Bjorken x
            double W2 = target.M()*target.M() + 2*target.M()*nu - Q2; // Calculate W^2
            double W = sqrt(W2); // Calculate W
            double y = nu / beamEnergy;
            double Pt2 = pip.Px()*pip.Px() + pip.Py()*pip.Py();


            hmiss->Fill(miss.M());
            hQ2->Fill(Q2);
            hxB->Fill(xB);
            hQ2vsxB->Fill(xB, Q2);
            hW2->Fill(W2); // Fill W^2 histogram
            hW->Fill(W); // Fill W histogram
            hPt2->Fill(Pt2); // Fill Pt^2 histogram
            if (Q2 > 1.5 && y> 0.25 && y < 0.85 && W > 2){
            hPt2->Fill(Pt2); // Draw Pt^2 histogram
             }
            
         }
    }

    TCanvas* can = new TCanvas("can", "Distributions", 2000, 400);
    can->Divide(4,2);
    can->cd(1);
    hQ2->Draw();
    can->cd(2);
    hxB->Draw();
    can->cd(3);
    hQ2vsxB->Draw("COLZ");
    can->cd(4);
    hW2->Draw(); // Draw W^2 histogram
    can->cd(5);
    hW->Draw(); // Draw W histogram
    can->cd(6);
    hmiss->Draw();
    can->cd(7);

    can->Print("can.pdf");
  
}
