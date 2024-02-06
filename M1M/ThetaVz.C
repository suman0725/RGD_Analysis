{
hists.Hist2D("PBANK.Vz:P.Theta*TMath::RadToDeg()",200,-30,20,200,0,180,"P.Pid==11")->Draw("colz");

hists.Save("ThetaVz.root");
}
