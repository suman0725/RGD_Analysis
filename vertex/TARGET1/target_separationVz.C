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

void target_separationVz(){   


	std::vector<std::pair<double, double>> vx_cuts(6), vy_cuts(6);
    	std::ifstream cutsFile("cuts_018564.txt");
    	std::string line;
   	int sector = 0;

    	while (std::getline(cutsFile, line) && sector < 6) {
        std::istringstream iss(line);
        double vx_low, vx_high, vy_low, vy_high;
        if (!(iss >> vx_low >> vx_high >> vy_low >> vy_high)) { break; } // error

        vx_cuts[sector] = std::make_pair(vx_low, vx_high);
        vy_cuts[sector] = std::make_pair(vy_low, vy_high);
        ++sector;
    	}
 
	for(int i = 0; i < 6; ++i) {
        std::cout << "Sector: " << (i+1)
                  << " vx_cut: (" << vx_cuts[i].first << ", " << vx_cuts[i].second << ")"
                  << " vy_cut: (" << vy_cuts[i].first << ", " << vy_cuts[i].second << ")"
                  << std::endl;
       }

	clas12root::HipoChain chain1;
	//chain1.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/018564/rec_clas_018564.*.hipo");
	//chain1.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon/018564/rec_clas_018564.evio.00320-00324.hipo");
	chain1.Add("/lustre19/expphy/volatile/clas12/rg-d/production/prod/v2_ib_CuSn/dst/recon/018348/rec_clas_018348.evio.*.hipo");
	chain1.db()->turnOffQADB(); 
	auto& c12s = chain1.C12ref();
	

	TCanvas* canVz = new TCanvas("canVz", "vz distribution", 1200, 800);
	canVz->Divide(3, 2);  
	TH1F* hvz[6];

	for (int i = 0; i < 6; ++i) {   
		hvz[i] = new TH1F(Form("hvz_S%d", i+1), Form("Vz for electrons (Sector %d); Vz (cm); Counts", i+1), 200, -22.0, 18.0);	

	}	


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
	
			if(sector != -1) {

			if( vx1 > vx_cuts[sector].first && vx1 < vx_cuts[sector].second && vy1 > vy_cuts[sector].first && vy1 < vy_cuts[sector].second) {
                       		   hvz[sector]->Fill(vz1);
				
                       		   hvz[sector]->Fill(vz1);
                                 }


           		}

		}

	}
	const int numSectors = 6;
	double fitStartRangevzCu[numSectors] = {-10, -10, -10.5, -10, -10, -10.5};
    	double fitEndRangevzCu[numSectors] = {-6, -4.5, -6, -6, -6.7, -6};

     	double fitStartRangevzSn[numSectors] = {-5,-5,-5,-6,-5,-6};
     	double fitEndRangevzSn[numSectors] = {0, 0, 0, -2, -2.4, -2};	
	 
	for (int i = 0; i < 6; ++i) {
        canVz->cd(i + 1);

         hvz[i]->SetStats(0);
                hvz[i]->Draw();

	
                  TF1* fitFuncvzCu = new TF1(Form("fitFunc_S%d", i + 1), "gaus(0)+pol3(3)", fitStartRangevzCu[i], fitEndRangevzCu[i]);
                                 fitFuncvzCu->SetParameters(1,0, 1, 1, 0, 1, 1);
                         hvz[i]->Fit(fitFuncvzCu, "R");

		fitFuncvzCu->SetLineColor(kBlue);
            fitFuncvzCu->SetLineWidth(1);
            fitFuncvzCu->Draw("same");

            double meanvzCu = fitFuncvzCu->GetParameter(1);
            double sigmavzCu = fitFuncvzCu->GetParameter(2);
            double line_startvzCu = meanvzCu - 3 * sigmavzCu;
            double line_endvzCu = meanvzCu + 3 * sigmavzCu;

            TLine* line1vzCu = new TLine(line_startvzCu, 0, line_startvzCu, fitFuncvzCu->GetMaximum());
            TLine* line2vzCu = new TLine(line_endvzCu, 0, line_endvzCu, fitFuncvzCu->GetMaximum());
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
	
	   
                  TF1* fitFuncvzSn = new TF1(Form("fitFunc_S%d", i + 1), "gaus(0)+pol3(3)", fitStartRangevzSn[i], fitEndRangevzSn[i]);
                                 fitFuncvzSn->SetParameters(1,0, 1, 1, 0, 1, 1);
                         hvz[i]->Fit(fitFuncvzSn, "R");

		fitFuncvzSn->SetLineColor(kGreen);
            fitFuncvzSn->SetLineWidth(1);
            fitFuncvzSn->Draw("same");

            double meanvzSn = fitFuncvzSn->GetParameter(1);
            double sigmavzSn = fitFuncvzSn->GetParameter(2);
            double line_startvzSn = meanvzSn - 3 * sigmavzSn;
            double line_endvzSn = meanvzSn + 3 * sigmavzSn;

            TLine* line1vzSn = new TLine(line_startvzSn, 0, line_startvzSn, fitFuncvzSn->GetMaximum());
            TLine* line2vzSn = new TLine(line_endvzSn, 0, line_endvzSn, fitFuncvzSn->GetMaximum());
            line1vzSn->SetLineColor(kBlack);
            line2vzSn->SetLineColor(kBlack);
            line1vzSn->Draw("same");
            line2vzSn->Draw("same");

            TLatex latexvzSn;
            latexvzSn.SetTextSize(0.025);
	    latexvzSn.SetTextAlign(33); // Correct comment to Align top-right
 	   latexvzSn.DrawLatexNDC(0.85, 0.85, Form("Mean = %.2f", meanvzSn)); // Adjusted for top-right
	   latexvzSn.DrawLatexNDC(0.85, 0.80, Form("Sigma = %.2f", sigmavzSn)); // Adjusted for top-right
	   latexvzSn.DrawLatexNDC(0.85, 0.75, Form("3-sigma cut: [%.2f, %.2f]", line_startvzSn, line_endvzSn)); // Adjusted for top-right
	   latexvzSn.DrawLatexNDC(0.85, 0.70, Form("p0 = %.2f", fitFuncvzSn->GetParameter(3))); // Adjusted for top-right
 	   latexvzSn.DrawLatexNDC(0.85, 0.65, Form("p1 = %.2f", fitFuncvzSn->GetParameter(4))); // Adjusted for top-right
	   latexvzSn.DrawLatexNDC(0.85, 0.60, Form("p2 = %.2f", fitFuncvzSn->GetParameter(5))); // Adjusted for top-right
 	   latexvzSn.DrawLatexNDC(0.85, 0.55, Form("p3 = %.2f", fitFuncvzSn->GetParameter(6))); // Adjusted for top-right
	canVz->Print("Vz_018564.pdf","pdf"); 
    }


