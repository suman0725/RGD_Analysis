void Momlabel(){
TFile *file = new TFile("Mom.root", "READ");
TList *histList = dynamic_cast<TList*>(file->Get("HipoHists"));
if (histList == nullptr) {
std::cerr << "Failed to retrieve the TList from the file." << std::endl;
return;
}

TCanvas *canvas = new TCanvas("canvas", "Histograms", 1200, 800);
canvas->Divide(2, 3);  // Create a 2x3 grid for plotting

for (int i = 0; i < 6; ++i) {
    TH1 *histogram = dynamic_cast<TH1*>(histList->At(i));
    if (histogram != nullptr) {
        canvas->cd(i + 1);
        histogram->GetXaxis()->SetTitle("P (Gev)");
        histogram->GetYaxis()->SetTitle("counts");
        histogram->Draw("COLZ");
    }
}
canvas->SaveAs("Mom_All.pdf");
file->Close();
}
