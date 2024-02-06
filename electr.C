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

void  electr(){

clas12root::HipoChain chain; 
//chain.Add(" /volatile/clas12/rg-d/production/trig/v3/trig/recon/018497/rec_clas_018497.evio.00385-00389.hipo");
//chain.Add("/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/018451/rec_clas_018451.evio.*.hipo");
chain.Add("/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/018451/rec_clas_018451.evio.01215-01219.hipo");
chain.db()->turnOffQADB();

//Get a c12 object for configuring 
auto config_c12 = chain.GetC12Reader();
auto& c12 = chain.C12ref();
TFile* outputFile = new TFile("electron.root", "RECREATE");
 
TH1F* hElectronMomentum = new TH1F("Momentum", "", 200, 0, 10);
hElectronMomentum->GetXaxis()->SetTitle("Momentum (GeV/c)"); 
hElectronMomentum->GetYaxis()->SetTitle("counts"); 

 TH2F* phiVmom = new TH2F("#phi Vs Momentum", "",200,0,11,200,-180,180); 
phiVmom->GetXaxis()->SetTitle("Momentum (GeV/c)");
phiVmom->GetYaxis()->SetTitle("#phi");

//TH2F* vzVmom = new TH2F("Momentum Vs Vz", "",200,0,11,200,-15,10); 
TH1F* vze = new TH1F("Vz", "",200,-15,10); 
vze->GetXaxis()->SetTitle("Momentum (GeV/c)");
vze->GetYaxis()->SetTitle("counts");

TH2F* vzVphi = new TH2F( "#phi Vs Vz", "",200,-180,180,200,-15,15); 
vzVphi->GetXaxis()->SetTitle("Vz (cm)");
vzVphi->GetYaxis()->SetTitle("#phi");

TH2F* htccxVhtccy= new TH2F("", "",200,-100,100,200,-100,100); 
htccxVhtccy->GetXaxis()->SetTitle("X");
htccxVhtccy->GetYaxis()->SetTitle("Y");

TH1F* phi= new TH1F("#phi", "",200,-150,150); 
phi->GetXaxis()->SetTitle("#phi");
phi->GetYaxis()->SetTitle("counts");

TH1F* theta = new TH1F("#theta", "",200, 0,50); 
theta->GetXaxis()->SetTitle("#theta");
theta->GetYaxis()->SetTitle("counts");


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
   auto thet = electron->getTheta(); 
   auto vz = electron->par()->getVz();
   auto htccx = electron->che(HTCC)->getX(); 
   auto htccy = electron->che(HTCC)->getY(); 


   hElectronMomentum->Fill(momentum);
   //hElectronMomentum->SetTitle("Momentum (GeV); Count");
   phiVmom->Fill(momentum,phie*TMath::RadToDeg());
   vze->Fill(vz);
   vzVphi-> Fill(phie*TMath::RadToDeg(),vz);
   htccxVhtccy->Fill(htccx,htccy);

  phi->Fill(phie*TMath::RadToDeg()); 
  theta->Fill(thet*TMath::RadToDeg()); 
  

   
  
}
counter++;
} 
 TCanvas* can = new TCanvas();    
// TCanvas* can=new TCanvas("a","a", 800,600);
   //can->SetTitle("My Can");
  can -> Divide(3,2);
 //can->cd(1);
 // vzVmom->DrawCopy("");
 can->cd(1);
  vze->DrawCopy("colz");
can->cd(2);
  vzVphi->DrawCopy("colz"); 
can->cd(3);
phiVmom->DrawCopy("colz"); 
can->cd(4);
hElectronMomentum->DrawCopy("colz"); 
//can->cd(5);
//htccxVhtccy->DrawCopy("colz");

 can->cd(5);
  phi->DrawCopy("colz");
 can->cd(6);
  theta->DrawCopy("colz");

outputFile->Close();
    delete outputFile;
}


