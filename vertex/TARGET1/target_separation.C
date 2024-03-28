#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TPaveText.h>
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

void target_separation(){

	clas12root::HipoChain chain;
	chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/018564/rec_clas_018564.*.hipo");
	//chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/018564/rec_clas_018564.evio.00320-00324.hipo");
	//chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ib_CuSn/dst/recon/018348/rec_clas_018348.evio.*.hipo");
	chain.db()->turnOffQADB(); 
	auto& c12 = chain.C12ref();

	// Create three canvases for vx, vy, and vz distributions
	TCanvas* canVx = new TCanvas("canVx", "vx distribution", 1200, 800);
	canVx->Divide(3, 2); // Divide into 3 columns and 2 rows
  
	TCanvas* canVy = new TCanvas("canVy", "vy distribution", 1200, 800);
	canVy->Divide(3, 2); 

	TCanvas* canVz = new TCanvas("canVz", "vz distribution", 1200, 800);
	canVz->Divide(3, 2); 


	// Initialize histograms for vx, vy, and vz for all 6 sectors
	TH1F* hvx[6];
	TH1F* hvy[6];
	TH1F* hvz[6];
		 
	// Sector-specific bounds for vx and vy based on your calculations

	for(int i = 0; i < 6; ++i) {

		if ( i == 1 || i == 4) {
   			hvx[i] = new TH1F(Form("hvx_S%d", i+1), Form("Vx for electrons (Sector %d); Vx (cm); Counts", i+1), 200, -3.0, 3.0);
		} else if (i == 2){
   			hvx[i] = new TH1F(Form("hvx_S%d", i+1), Form("Vx for electrons (Sector %d); Vx (cm); Counts", i+1), 200, -1, 1);
		} 		 else  {
   			hvx[i] = new TH1F(Form("hvx_S%d", i+1), Form("Vx for electrons (Sector %d); Vx (cm); Counts", i+1), 200, -1.5, 1.5);
		  }


		if ( i == 0){    	
		hvy[i] = new TH1F(Form("hvy_S%d", i+1), Form("Vy for electrons (Sector %d); Vy (cm); Counts", i+1), 200, -3.0, 3.0);
		} else if ( i == 1 || i == 4){
		     hvy[i] = new TH1F(Form("hvy_S%d", i+1), Form("Vy for electrons (Sector %d); Vy (cm); Counts", i+1), 200, -1.5, 1.5);
		  } else {
		      hvy[i] = new TH1F(Form("hvy_S%d", i+1), Form("Vy for electrons (Sector %d); Vy (cm); Counts", i+1), 200, -2.5, 2.5);
		    }
	
    		hvz[i] = new TH1F(Form("hvz_S%d", i+1), Form("Vz for electrons (Sector %d); Vz (cm); Counts", i+1), 200, -22.0, 18.0);
			
	}


		
	
	while (chain.Next()) {
		auto electrons = c12->getByID(11); 
		double maxMomentum = -1.0;
		clas12::region_part_ptr  highestMomentum =  nullptr;
		
		for (auto& electron : electrons) {

			if (electron -> getRegion() == FD) {  
				double momentum = electron ->par()->getP();
				if (momentum > maxMomentum) {  
					maxMomentum = momentum; 
					highestMomentum = electron; 
			}

			}



		}
	

		if (highestMomentum != nullptr) {
            		//std::cout << "Highest momentum of an electron in FD for this event: " << maxMomentum << " GeV/c" << std::endl;
			auto phi = highestMomentum->getPhi() * TMath::RadToDeg();
        		auto vx = highestMomentum->par()->getVx();
        		auto vy = highestMomentum->par()->getVy();
        		auto vz = highestMomentum->par()->getVz();

			int sector = -1;
                 	if(phi >= -180 && phi < -120) sector = 0; // Sector 1
       			else if(phi >= -120 && phi < -60) sector = 1; // Sector 2
        		else if(phi >= -60 && phi < 0) sector = 2;  // Sector 3
        		else if(phi >= 0 && phi < 60) sector = 3; // Sector 4
        		else if(phi >= 60 && phi < 120) sector = 4; // Sector 5
        		else if(phi >= 120 && phi < 180) sector = 5; // Sector 6
			
			if(sector != -1) {
            			hvx[sector]->Fill(vx);
            			hvy[sector]->Fill(vy);
            		
         		  }


        	}

		

	}
	const int numSectors = 6; // Assuming you have 6 sectors
    	double fitStartRange[numSectors] = {-1.5, -3, -1, -1.5, -2, -1.5};
     	double fitEndRange[numSectors] = {1.5, 3, 1, 1.5, 2, 1.5};

    	double fitStartRangevy[numSectors] = {-1.5, -3, -2.5, -3, -2, -1.5};
     	double fitEndRangevy[numSectors] = {1.5, 3, 2.5, 3, 2, 1.5};

    	double fitStartRangevzCu[numSectors] = {-10.5, -10, -10, -10.5, -10.5, -10.5};
    	double fitEndRangevzCu[numSectors] = {-4.5, -4.5, -4.5, -5, -5, -5};

     	double fitStartRangevzSn[numSectors] = {1.5, 3, 2.5, 3, 2, 1.5};
     	double fitEndRangevzSn[numSectors] = {1.5, 3, 2.5, 3, 2, 1.5};
	
	double vx_cuts[6][2]; // For storing lower and upper bounds of the 3-sigma cut for vx
	double vy_cuts[6][2]; 	
 
	for(int i = 0; i < 6; ++i) {

	 canVx->cd(i+1);
	 hvx[i]->SetStats(0); 
	 hvx[i]->Draw();

	 TF1* fitFunc = new TF1(Form("fitFunc_S%d", i+1), "gaus(0)+pol3(3)", fitStartRange[i], fitEndRange[i]);	
	 fitFunc->SetParameters(1, 0, 1, 0, 1, 1);
         hvx[i]->Fit(fitFunc, "R");

         fitFunc->SetLineColor(kBlue);
         fitFunc->SetLineWidth(1);
         fitFunc->Draw("same");

	 double mean = fitFunc->GetParameter(1);
         double sigma = fitFunc->GetParameter(2);
          vx_cuts[i][0] = mean - 3 * sigma;
          vx_cuts[i][1] = mean + 3 * sigma;

	 TLine *line1 = new TLine(vx_cuts[i][0], 0, vx_cuts[i][0], fitFunc->GetMaximum());
         TLine *line2 = new TLine(vx_cuts[i][1], 0, vx_cuts[i][1], fitFunc->GetMaximum());
         line1->SetLineColor(kRed);
         line2->SetLineColor(kRed);
         line1->Draw("same");
         line2->Draw("same");


	 TLatex latex;
         latex.SetTextSize(0.025);
         latex.SetTextAlign(13); // Align top-left
         latex.DrawLatexNDC(0.15, 0.85, Form("Mean = %.2f", mean));
         latex.DrawLatexNDC(0.15, 0.80, Form("Sigma = %.2f", sigma));
         latex.DrawLatexNDC(0.15, 0.75, Form("3-sigma cut: [%.2f, %.2f]", vx_cuts[i][0], vx_cuts[i][1]));
         latex.DrawLatexNDC(0.15, 0.70, Form("p0 = %.2f", fitFunc->GetParameter(3))); // Coeff of pol1
         latex.DrawLatexNDC(0.15, 0.65, Form("p1 = %.2f", fitFunc->GetParameter(4))); // Coeff of pol1
	 latex.DrawLatexNDC(0.15, 0.60, Form("p2 = %.2f", fitFunc->GetParameter(5))); // Coeff of pol2
	 latex.DrawLatexNDC(0.15, 0.55, Form("p3 = %.2f", fitFunc->GetParameter(6))); // Cubic term of polynomial


	
	 canVy->cd(i+1);
	// hvy[i]->SetStats(0); 
	 hvy[i]->Draw();
		
	 TF1* fitFuncvy = new TF1(Form("fitFunc_S%d", i+1), "gaus(0)+pol3(3)", fitStartRangevy[i], fitEndRangevy[i]);	
	 fitFuncvy->SetParameters(1, 0, 1, 1, 0, 1, 1); 
         hvy[i]->Fit(fitFuncvy, "R+");
         
	 fitFuncvy->SetLineColor(kRed);
         fitFuncvy->SetLineWidth(1);
         fitFuncvy->Draw("same");

	 double meanvy = fitFuncvy->GetParameter(1);
         double sigmavy = fitFuncvy->GetParameter(2);
         vy_cuts[i][0]= meanvy - 3 * sigmavy;
          vy_cuts[i][1]= meanvy + 3 * sigmavy;
	 
	 TLine *line1vy = new TLine(vy_cuts[i][0], 0, vy_cuts[i][0], fitFuncvy->GetMaximum());
         TLine *line2vy = new TLine(vy_cuts[i][1], 0, vy_cuts[i][1], fitFuncvy->GetMaximum());
         line1vy->SetLineColor(kBlue);
         line2vy->SetLineColor(kBlue);
         line1vy->Draw("same");
         line2vy->Draw("same");

	 TLatex latexvy;
         latexvy.SetTextSize(0.025);
         latexvy.SetTextAlign(13); // Align top-left
         latexvy.DrawLatexNDC(0.15, 0.85, Form("Meanvy = %.2f", meanvy));
         latexvy.DrawLatexNDC(0.15, 0.80, Form("Sigmavy = %.2f", sigmavy));
         latexvy.DrawLatexNDC(0.15, 0.75, Form("3-sigma cut: [%.2f, %.2f]", vy_cuts[i][0],vy_cuts[i][1] ));
         latexvy.DrawLatexNDC(0.15, 0.70, Form("p0 = %.2f", fitFuncvy->GetParameter(3))); // Coeff of pol1
         latexvy.DrawLatexNDC(0.15, 0.65, Form("p1 = %.2f", fitFuncvy->GetParameter(4))); // Coeff of pol1
	 latexvy.DrawLatexNDC(0.15, 0.60, Form("p2 = %.2f", fitFuncvy->GetParameter(5))); // Coeff of pol2
	 latexvy.DrawLatexNDC(0.15, 0.55, Form("p3 = %.2f", fitFuncvy->GetParameter(6))); // Cubic term of polynomial
	
	
	std::ofstream cutsFile("cuts_018564.txt");
	for (int i = 0; i < 6; ++i) {
    	cutsFile << vx_cuts[i][0] << " " << vx_cuts[i][1] << " " << vy_cuts[i][0] << " " << vy_cuts[i][1] << std::endl;
	}
	cutsFile.close();
       
	TFile *outfile = new TFile("Vx_Vy_electron_018564.root", "RECREATE");
	for(int i = 0; i < 6; ++i) {
    	hvx[i]->Write();
    	hvy[i]->Write();
	}
	outfile->Close(); 
	canVx->Print("VxVy_018564.pdf(","pdf");
	canVy->Print("VxVy_018564.pdf)","pdf"); 
	}
	/*clas12root::HipoChain chain1;
	//chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/018564/rec_clas_018564.*.hipo");
	chain1.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/018564/rec_clas_018564.evio.00320-00324.hipo");
	//chain.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ib_CuSn/dst/recon/018348/rec_clas_018348.evio.*.hipo");
	chain1.db()->turnOffQADB(); 
	auto& c12s = chain.C12ref();



	while (chain1.Next()) {
		auto electrons1 = c12s->getByID(11); 
		double maxMomentum1 = -1.0;
		clas12::region_part_ptr  highestMomentum1 =  nullptr;
		
		for (auto& electron1 : electrons1) {

			if (electron1 -> getRegion() == FD) {  
				double momentum1 = electron1 ->par()->getP();
				if (momentum1 > maxMomentum1) {  
					maxMomentum1 = momentum1; 
					highestMomentum1 = electron1; 
			}

			}



		}
	

		if (highestMomentum1 != nullptr) {
            		//std::cout << "Highest momentum of an electron in FD for this event: " << maxMomentum1 << " GeV/c" << std::endl;
			auto phi1 = highestMomentum1->getPhi() * TMath::RadToDeg();
        		auto vx1 = highestMomentum1->par()->getVx();
        		auto vy1 = highestMomentum1->par()->getVy();
        		auto vz1 = highestMomentum1->par()->getVz();

			int sector = -1;
                 	if(phi1 >= -180 && phi1 < -120) sector = 0; // Sector 1
       			else if(phi1 >= -120 && phi1 < -60) sector = 1; // Sector 2
        		else if(phi1 >= -60 && phi1 < 0) sector = 2;  // Sector 3
        		else if(phi1 >= 0 && phi1 < 60) sector = 3; // Sector 4
        		else if(phi1 >= 60 && phi1 < 120) sector = 4; // Sector 5
        		else if(phi1 >= 120 && phi1 < 180) sector = 5; // Sector 6
                       		   hvz[sector]->Fill(vz1);
			
			if(sector != -1) {

    				if( vx1 > vx_cuts[sector][0] && vx1 < vx_cuts[sector][1] && vy1 > vy_cuts[sector][0] && vy1 < vy_cuts[sector][1]) {
                       		   hvz[sector]->Fill(vz1);
				}
                       		   hvz[sector]->Fill(vz1);
                                                        }


        	}

		

	}

	for (int i = 0; i < 6 ; ++i){

	 canVz->cd(i+1); 
	 // hvz[i]->SetStats(0); 
	 hvz[i]->Draw();
	 
	 TF1* fitFuncvzCu = new TF1(Form("fitFunc_S%d", i+1), "gaus(0)+pol3(3)", fitStartRangevzCu[i], fitEndRangevzCu[i]);	
	 fitFuncvzCu->SetParameters(1, 1, 1, 0, 1, 1);
         hvz[i]->Fit(fitFuncvzCu, "R");

         fitFuncvzCu->SetLineColor(kBlue);
         fitFuncvzCu->SetLineWidth(1);
         fitFuncvzCu->Draw("same");

	 double meanvzCu = fitFuncvzCu->GetParameter(1);
         double sigmavzCu = fitFuncvzCu->GetParameter(2);
         double line_startvzCu = meanvzCu - 3 * sigmavzCu;
         double line_endvzCu = meanvzCu + 3 * sigmavzCu;

	 TLine *line1vzCu = new TLine(line_startvzCu, 0, line_startvzCu, fitFuncvzCu->GetMaximum());
         TLine *line2vzCu = new TLine(line_endvzCu, 0, line_endvzCu, fitFuncvzCu->GetMaximum());
         line1vzCu->SetLineColor(kRed);
         line2vzCu->SetLineColor(kRed);
         line1vzCu->Draw("same");
         line2vzCu->Draw("same");


	 TLatex latexvzCu;
         latexvzCu.SetTextSize(0.025);
         latexvzCu.SetTextAlign(13); // Align top-left
         latexvzCu.DrawLatexNDC(0.15, 0.85, Form("Mean = %.2f", meanvzCu));
         latexvzCu.DrawLatexNDC(0.15, 0.80, Form("Sigma = %.2f", sigmavzCu));
         latexvzCu.DrawLatexNDC(0.15, 0.75, Form("3-sigma cut: [%.2f, %.2f]", line_startvzCu, line_endvzCu));
         latexvzCu.DrawLatexNDC(0.15, 0.70, Form("p0 = %.2f", fitFuncvzCu->GetParameter(3))); // Coeff of pol1
         latexvzCu.DrawLatexNDC(0.15, 0.65, Form("p1 = %.2f", fitFuncvzCu->GetParameter(4))); // Coeff of pol1
	 latexvzCu.DrawLatexNDC(0.15, 0.60, Form("p2 = %.2f", fitFuncvzCu->GetParameter(5))); // Coeff of pol2
	 latexvzCu.DrawLatexNDC(0.15, 0.55, Form("p3 = %.2f", fitFuncvzCu->GetParameter(6))); // Cubic term of polynomial

	
	}*/





}
