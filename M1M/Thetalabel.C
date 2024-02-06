void Thetalabel(){
TFile *file = new TFile("PhiTheta.root", "READ");
TList *histList = dynamic_cast<TList*>(file->Get("HipoHists"));
if (histList == nullptr) {
    std::cerr << "Failed to retrieve the TList from the file." << std::endl;
    return 1;
}

TH1 *histogram = dynamic_cast<TH1*>(histList->At(1));
histogram->Rebin(0.5);
//histogram->Smooth(3);
//histogram.SetRangeX(0, 40)
histogram->GetXaxis()->SetRangeUser(0, 40);
//histogram->GetXaxis()->SetRangeUser(-30, 30);
// Change the X-axis label
 histogram->GetXaxis()->SetTitle("#theta (angle in degrees)");
 histogram->GetYaxis()->SetTitle("count");

//
// // Change the Y-axis label
 histogram->GetYaxis()->SetTitle("#Theta (angle in degrees)");
 //histogram->GetYaxis()->SetTitle("#beta");
//
TCanvas *canvas = new TCanvas("canvas", "Histogram", 1200, 800);
histogram->Draw("COLZ");
canvas->Update();
canvas->Modified();
    canvas->SaveAs("Theta.pdf");
file->Close();
}
