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

void neutralpart() {
    clas12root::HipoChain chain;
    chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v0_ib_LD2/dst/recon/018321/rec_clas_018321.evio.00085-00089.hipo");
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref();

    TH1F* neutralMomentum = new TH1F("neutralMomentum", "Momentum for Neutral Particles; Momentum (GeV/c); Counts", 200, 0, 10);
    TH1F* hVz = new TH1F("vz", "Vertex (Vz) for Neutral Particles; Vz (cm); Counts", 200, -15, 15);
    TH1F* hPhi = new TH1F("phi", "Phi for Neutral Particles; Phi (degrees); Counts", 200, -180, 180);
    TH1F* hTheta = new TH1F("theta", "Theta for Neutral Particles; Theta (degrees); Counts", 200, 0, 180);
    TH1F* hTrackMultiplicity = new TH1F("trackMultiplicity", "Track Multiplicity for Neutral Particles; Track Multiplicity;", 50, 0, 50);
    TH2F* phiVmom = new TH2F("phiVmom", "Phi vs Momentum for Neutral Particles; Momentum (GeV/c); Phi (degrees)", 200, 0, 10, 200, -180, 180);
    TH2F* vzVmom = new TH2F("vzVmom", "Vz vs Momentum for Neutral Particles; Momentum (GeV/c); Vz (cm)", 200, 0, 10, 200, -15, 15);
    TH2F* vzVphi = new TH2F("vzVphi", "Vz vs Phi for Neutral Particles; Phi (degrees); Vz (cm)", 200, -180, 180, 200, -15, 15);
    TH2F* betaVmom = new TH2F("betaVmom", "Beta vs Momentum for Neutral Particles; Momentum (GeV/c); Beta", 200, 0, 10, 200, 0, 1);

    while (chain.Next()) {
        auto neutrals = c12->getByCharge(0);
        int neutralTrackCount = neutrals.size();
        hTrackMultiplicity->Fill(neutralTrackCount);

        for (const auto& neutral : neutrals) {
            auto px = neutral->par()->getPx();
            auto py = neutral->par()->getPy();
            auto pz = neutral->par()->getPz();
            auto momentum = TMath::Sqrt(px * px + py * py + pz * pz);
            auto phie = neutral->getPhi() * TMath::RadToDeg();
            auto thetae = neutral->getTheta() * TMath::RadToDeg();
            auto vz = neutral->par()->getVz();
            auto beta = neutral->par()->getBeta();

            neutralMomentum->Fill(momentum);
            hPhi->Fill(phie);
            hTheta->Fill(thetae);
            vzVmom->Fill(momentum, vz);
            phiVmom->Fill(momentum, phie);
            hVz->Fill(vz);
            vzVphi->Fill(phie, vz);
            betaVmom->Fill(momentum, beta);
        }
    }

    TCanvas* can1 = new TCanvas("can1", "Canvas 1", 800, 600);
    can1->Divide(2, 1);
    can1->cd(1); neutralMomentum->Draw();
    can1->cd(2); hVz->Draw();
    can1->Print("neutral_particles_results.pdf(", "pdf");

    TCanvas* can2 = new TCanvas("can2", "Canvas 2", 800, 600);
    can2->Divide(2, 1);
    can2->cd(1); phiVmom->Draw("COLZ");
    can2->cd(2); vzVphi->Draw("COLZ");
    can2->Print("neutral_particles_results.pdf", "pdf");

    TCanvas* can3 = new TCanvas("can3", "Canvas 3", 800, 600);
    can3->Divide(2, 1);
    can3->cd(1); betaVmom->Draw("COLZ");
    can3->cd(2); hPhi->Draw();
    can3->Print("neutral_particles_results.pdf", "pdf");

    TCanvas* can4 = new TCanvas("can4", "Canvas 4", 800, 600);
    can4->Divide(2, 1);
    can4->cd(1); hTheta->Draw();
    can4->cd(2); hTrackMultiplicity->Draw();
    can4->Print("neutral_particles_results.pdf", "pdf");

    TCanvas* can5 = new TCanvas("can5", "Canvas 5", 800, 600);
    can5->Divide(1, 1);
    can5->cd(1); vzVmom->Draw("COLZ");
    can5->Print("neutral_particles_results.pdf", "pdf");

    can1->Print("neutral_particles_results.pdf)", "pdf");
}



