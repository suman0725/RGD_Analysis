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
//chain.Add("/volatile/clas12/rg-d/production/trig/v3/trig/recon/018497/rec_clas_018497.evio.*.hipo");
chain.Add("/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/018451/rec_clas_018451.evio.01215-01219.hipo");
chain.db()->turnOffQADB();

//Get a c12 object for configuring 
auto config_c12 = chain.GetC12Reader();
auto& c12 = chain.C12ref();
 
TH1F* negMomentum = new TH1F("negMomentum", "", 200, 0, 10);
negMomentum->GetXaxis()->SetTitle("Momentum (GeV/c)"); 
negMomentum->GetYaxis()->SetTitle("counts)"); 

 TH2F* phiVmom = new TH2F("phiVmom", "PhiVSMomentumforneg",200,0,11,200,-180,180); 
phiVmom->GetXaxis()->SetTitle("Momentum (GeV/c)");
phiVmom->GetYaxis()->SetTitle("#phi");

/*TH2F* vzVmom = new TH2F("vzVmom", "VzVSMomentumforneg",200,0,11,200,-15,15); 
vzVmom->GetXaxis()->SetTitle("Momentum (GeV/c)");
vzVmom->GetYaxis()->SetTitle("Vz (cm)");*/


TH1F* vzVmom = new TH1F("vz", "Vzforneg",200,-15,15); 
vzVmom->GetXaxis()->SetTitle("Vz (cm)");
vzVmom->GetYaxis()->SetTitle("counts");

TH2F* vzVphi = new TH2F("vzVphi", "VzVSPhiforneg",200,-180,180,200,-15,15); 
vzVphi->GetXaxis()->SetTitle("Vz (cm)");
vzVphi->GetYaxis()->SetTitle("#phi");

TH2F* htccxVhtccy= new TH2F("htccxVhtccy", "HtccxVSHtccy(position)forneg",200,-100,100,200,-100,100); 
htccxVhtccy->GetXaxis()->SetTitle("X");
htccxVhtccy->GetYaxis()->SetTitle("Y");


int counter = 0 ; 

while (chain.Next()){

auto negs = c12->getByCharge(-1);
   /*int torus; 
   std::ofstream outFile("torus.txt");
   outFile << "Electrons are " << torus << " polarity"<<endl; 
   torus =c12->runconfig()->getTorus(); 
   outFile.close(); */


 for (const auto& neg : negs) {
   
   auto px = neg->par()->getPx();
   auto py = neg->par()->getPy();
   auto pz = neg->par()->getPz();
   auto momentum = TMath::Sqrt(px * px + py * py + pz * pz);
   auto phie = neg->getPhi();
   auto vz = neg->par()->getVz();
   auto htccx = neg->che(HTCC)->getX(); 
   auto htccy = neg->che(HTCC)->getY(); 


   negMomentum->Fill(momentum);
   //negMomentum->SetTitle("Momentum (GeV); Count");
   phiVmom->Fill(momentum,phie*TMath::RadToDeg());
   //vzVmom->Fill(vz,momentum);
   vzVmom->Fill(vz);
   vzVphi-> Fill(phie*TMath::RadToDeg(),vz);
   htccxVhtccy->Fill(htccx,htccy);

   
  
}
     counter++;
} 
     TCanvas* can=new TCanvas();
  can -> Divide(3,2);
 can->cd(1);
  vzVmom->DrawCopy("colz");
can->cd(2);
  vzVphi->DrawCopy("colz"); 
can->cd(3);
phiVmom->DrawCopy("colz"); 
can->cd(4);
negMomentum->DrawCopy("colz"); 
//can->cd(5);
//htccxVhtccy->DrawCopy("colz");

}

