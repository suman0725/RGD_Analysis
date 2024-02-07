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

/*double region_particle::getTheta() {
    _parts->setEntry(_pentry);
    double x=getPx();
    double y=getPy();
    double z=getPz();
    return  x == 0.0 && y == 0.0 && z == 0.0 ? 0.0
      : atan2(sqrt(x*x+y*y),z);
  }
  double region_particle::getPhi() {
    _parts->setEntry(_pentry);
    double x=getPx();
    double y=getPy();
    return atan2(y,x);
  }*/

void  negpart(){

    clas12root::HipoChain chain;
    chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v0_ib_LD2/dst/recon/018321/rec_clas_018321.evio.00085-00089.hipo"); 
    //chain.Add("/volatile/clas12/rg-d/production/trig/v3/trig/recon/018497/rec_clas_018497.evio.*.hipo");
    //chain.Add("/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/018451/rec_clas_018451.evio.01215-01219.hipo");
    chain.db()->turnOffQADB();

    //Get a c12 object for configuring 
    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();
 
    TH1F* negMomentum = new TH1F("negMomentum", "Momentum for Negative Particles; Momentum (GeV/c); Counts", 200, 0, 10);
    TH1F* hVz = new TH1F("vz", "Vertex (Vz) for Negative Particles; Vz (cm); Counts", 200, -15, 15);
    TH1F* hPhi = new TH1F("phi", "Phi for Negative Particles; Phi (degrees); Counts", 200, -180, 180);
    TH1F* hTheta = new TH1F("theta", "Theta for Negative Particles; Theta (degrees); Counts", 200, 0, 180);
    TH1F* hTrackMultiplicity = new TH1F("trackMultiplicity", "Track Multiplicity for Negative Particles;Track Multiplicity;", 15, 0, 15);
    TH2F* phiVmom = new TH2F("phiVmom", "Phi vs Momentum for Negative Particles; Momentum (GeV/c); Phi (degrees)", 200, 0, 10, 200, -180, 180);
    TH2F* vzVmom = new TH2F("vzVmom", "Vz vs Momentum for Negative Particles; Momentum (GeV/c); Vz (cm)", 200, 0, 10, 200, -15, 15);
    TH2F* vzVphi = new TH2F("vzVphi", "Vz vs Phi for Negative Particles; Phi (degrees); Vz (cm)", 200, -180, 180, 200, -15, 15);
    TH2F* betaVmom = new TH2F("betaVmom", "Beta vs Momentum for Negative Particles; Momentum (GeV/c); Beta", 200, 0, 10, 200, 0, 1);

 

    while (chain.Next()){

         auto negs = c12->getByCharge(-1);
         int negTrackCount = negs.size();
         hTrackMultiplicity->Fill(negTrackCount);

         for (const auto& neg : negs) {
             auto px = neg->par()->getPx();
             auto py = neg->par()->getPy();
             auto pz = neg->par()->getPz();
             auto momentum = TMath::Sqrt(px * px + py * py + pz * pz);
             auto phie = neg->getPhi() * TMath::RadToDeg();
             auto thetae = neg->getTheta() * TMath::RadToDeg();
             auto vz = neg->par()->getVz();
             auto beta = neg->par()->getBeta();
             
             negMomentum->Fill(momentum);
             hPhi->Fill(phie);
             hTheta->Fill(thetae);
             vzVmom->Fill(momentum, vz);
             phiVmom->Fill(momentum, phie);
             hVz->Fill(vz);
             vzVphi->Fill(phie, vz);
             betaVmom->Fill(momentum, beta);
          }
    } 
    
    /*TCanvas* can = new TCanvas("can", "Negative Particles Analysis", 1200, 800);
    can->Divide(5, 2);
    can->cd(1); hVz->Draw();
    can->cd(2); negMomentum->Draw();
    can->cd(3); vzVphi->Draw("COLZ");
    can->cd(4); hPhi->Draw();
    can->cd(5); hTheta->Draw();
    can->cd(6); hTrackMultiplicity->Draw();
    can->cd(7); phiVmom->Draw("COLZ");
    can->cd(8); vzVmom->Draw("COLZ");
    can->cd(9); betaVmom->Draw("COLZ");
    
    can->Print("negative_particles_analysis_inbending_LD2_018321.pdf");*/
    TCanvas* can1 = new TCanvas("can1", "Canvas 1", 800, 600);
    can1->Divide(2, 1); // Adjust as necessary
    can1->cd(1);
    negMomentum->Draw();
    can1->cd(2);
    hVz->Draw();
    can1->Print("negative_particles_results.pdf(", "pdf");

    TCanvas* can2 = new TCanvas("can2", "Canvas 2", 800, 600);
    can2->Divide(2, 1); // Adjust as necessary
    can2->cd(1);
    phiVmom->Draw("COLZ");
    can2->cd(2);
    vzVphi->Draw("COLZ");
    can2->Print("negative_particles_results.pdf", "pdf");
    
    TCanvas* can3 = new TCanvas("can3", "Canvas 3", 800, 600);
    can3->Divide(2, 1); 
    can3->cd(1);
    betaVmom->Draw("COLZ");
    can3->cd(2); 
    hPhi->Draw(); 
    can3->Print("negative_particles_results.pdf", "pdf");
    
    TCanvas* can4 = new TCanvas("can4", "Canvas 4", 800, 600);
    can4->Divide(2, 1); 
    can4->cd(1);
    hTheta->Draw();
    can4->cd(2); 
    hTrackMultiplicity->Draw(); 
    can4->Print("negative_particles_results.pdf", "pdf");
   
    TCanvas* can5 = new TCanvas("can5", "Canvas 5", 800, 600);
    can5->Divide(1, 1); // This canvas has only one plot
    can5->cd(1);
    vzVmom->Draw("COLZ");
    can5->Print("negative_particles_results.pdf", "pdf");

    can1->Print("negative_particles_results.pdf)", "pdf");




}

