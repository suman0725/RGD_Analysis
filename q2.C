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

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}


void q2(){

   auto db=TDatabasePDG::Instance();
   TLorentzVector beam(0,0,10.532,10.532);

   auto* hq2 = new TH1F("hq2", "hq2", 200, 1, 10); 

   clas12root:: HipoChain chain; 
   chain.Add("/volatile/clas12/rg-d/production/trig/v3/trig/recon/018497/rec_clas_018497.evio.00385-00389.hipo"); 
   auto config_c12=chain.GetC12Reader(); 
   
   auto& c12=chain.C12ref(); 
   
   int counter=0; 
    
   while (chain.Next()){ 

   auto electrons=c12->getByID(11); 
   auto gamma=c12->getByID(22); 
   




}

