{
hists.Hist2D("PBANK.Vz:P.Phi*TMath::RadToDeg()",200,-20,20,200,-180,180,"P.Pid==11")->Draw("colz");

hists.Save("PhiVz.root");
}
