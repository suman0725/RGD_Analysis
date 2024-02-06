void label(){
TFile *file = new TFile("PhiVz.root", "READ");
TList *histList = dynamic_cast<TList*>(file->Get("HipoHists"));
if (histList == nullptr) {
    std::cerr << "Failed to retrieve the TList from the file." << std::endl;
    return 1;
}

TH1 *histogram = dynamic_cast<TH1*>(histList->At(0));
//histogram.SetRangeX(0, 60)
// Change the X-axis label
 //histogram->GetXaxis()->SetTitle("#theta (angle in degrees)");
 histogram->GetXaxis()->SetTitle("V_z (cm)");

//
// // Change the Y-axis label
 histogram->GetYaxis()->SetTitle("#phi (angle in degrees)");
 //histogram->GetYaxis()->SetTitle("#beta");
//
TCanvas *canvas = new TCanvas("canvas", "Histogram", 800, 600);
histogram->Draw("COLZ");
canvas->Update();
canvas->Modified();
    canvas->SaveAs("PhiVz.pdf");
file->Close();
}
