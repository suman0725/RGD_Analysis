{
hists.Hist1D("PBANK.Vz",200,-20,20,"P.Pid==11&&P.Sector==1");
hists.Hist1D("PBANK.Vz",200,-20,20,"P.Pid==11&&P.Sector==2");
hists.Hist1D("PBANK.Vz",200,-20,20,"P.Pid==11&&P.Sector==3");
hists.Hist1D("PBANK.Vz",200,-20,20,"P.Pid==11&&P.Sector==4");
hists.Hist1D("PBANK.Vz",200,-20,20,"P.Pid==11&&P.Sector==5");
hists.Hist1D("PBANK.Vz",200,-20,20,"P.Pid==11&&P.Sector==6")->Draw("(2x3)colz");

hists.Save("Vz.root");
}
