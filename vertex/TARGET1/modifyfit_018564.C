#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iostream>
#include <TH1F.h>

	void modifyfit_018564(){    




	TFile *file = TFile::Open("Vx_Vy_electron_018564.root", "UPDATE");
        if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }
	TH1F *h = (TH1F*)file->Get("hvx_S1"); // Use your histogram's name
    if (!h) {
        std::cerr << "Histogram not found." << std::endl;
        file->Close();
        delete file;
        return;
    }


	TF1 *fitFunc = h->GetFunction("fitFunc_S1"); // Use your function's name
        if (!fitFunc) {
        std::cerr << "Fit function not found." << std::endl;
        file->Close();
        delete file;
        return;
    }

	TH1F* hvx[6];        
	fitFunc->SetParameter(0,100); 
	fitFunc->SetParameter(1,0.3);
	fitFunc->SetRange(-1,1);

	TCanvas *canVx = new TCanvas("canVx", "Canvas Title", 1200, 6800);
	hvx[0]->Draw();
	fitFunc->Draw("same");
	canVx->Update();




}
