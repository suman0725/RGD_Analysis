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

void  electrsec(){

clas12root::HipoChain chain; 
//chain.Add("/volatile/clas12/rg-d/production/trig/v3/trig/recon/018497/rec_clas_018497.evio.*.hipo");
chain.Add("/volatile/clas12/rg-d/production/trig/v3/trig/recon/018497/rec_clas_018497.evio.00385-00389.hipo");
chain.db()->turnOffQADB();

//Get a c12 object for configuring 
auto config_c12 = chain.GetC12Reader();
auto& c12 = chain.C12ref();
 
TH1F* hElectronMomentum = new TH1F("hElectronMomentum", "", 200, 0, 10);
hElectronMomentum->GetXaxis()->SetTitle("Momentum (GeV/c)"); 
hElectronMomentum->GetYaxis()->SetTitle("counts)"); 

 TH2F* phiVmom = new TH2F("phiVmom", "PhiVSMomentumforelectron",200,0,11,200,-180,180); 
phiVmom->GetXaxis()->SetTitle("Momentum (GeV/c)");
phiVmom->GetYaxis()->SetTitle("#phi");

TH2F* vzVmom = new TH2F("vzVmom", "VzVSMomentumforelectron",200,0,11,200,-15,15); 
vzVmom->GetXaxis()->SetTitle("Momentum (GeV/c)");
vzVmom->GetYaxis()->SetTitle("Vz (cm)");

TH2F* vzVphi = new TH2F("vzVphi", "VzVSPhiforelectron",200,-180,180,200,-15,15); 
vzVphi->GetXaxis()->SetTitle("Vz (cm)");
vzVphi->GetYaxis()->SetTitle("#phi");

TH2F* htccxVhtccy= new TH2F("htccxVhtccy", "HtccxVSHtccy(position)forelectron",200,-100,100,200,-100,100); 
htccxVhtccy->GetXaxis()->SetTitle("X");
htccxVhtccy->GetYaxis()->SetTitle("Y");


int counter = 0 ; 

while (chain.Next()){

auto electrons = c12->getByID(11);
   /*int torus; 
   std::ofstream outFile("torus.txt");
   outFile << "Electrons are " << torus << " polarity"<<endl; 
   torus =c12->runconfig()->getTorus(); 
   outFile.close(); */


 for (const auto& electron : electrons) {
   
   auto px = electron->par()->getPx();
   auto py = electron->par()->getPy();
   auto pz = electron->par()->getPz();
   auto momentum = TMath::Sqrt(px * px + py * py + pz * pz);
   auto phie = electron->getPhi();
   auto vz = electron->par()->getVz();
   auto htccx = electron->che(HTCC)->getX(); 
   auto htccy = electron->che(HTCC)->getY(); 


   hElectronMomentum->Fill(momentum);
   //hElectronMomentum->SetTitle("Momentum (GeV); Count");
   phiVmom->Fill(momentum,phie*TMath::RadToDeg());
   vzVmom->Fill(vz,momentum);
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
hElectronMomentum->DrawCopy("colz"); 
//can->cd(5);
//htccxVhtccy->DrawCopy("colz");

}

