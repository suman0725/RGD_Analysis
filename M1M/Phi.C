{
hists.Hist1D("P.Phi*TMath::RadToDeg()",200,-180,180,"P.Pid==11");
hists.Hist1D("P.Theta*TMath::RadToDeg()",200,-180,180,"P.Pid==11")->Draw("(2x1)");

hists.Save("Phi&Theta.root");
}
