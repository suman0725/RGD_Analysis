{
hists.Hist1D("P.P",200,1,11,"P.Pid==11&&P.Sector==1");
hists.Hist1D("P.P",200,1,11,"P.Pid==11&&P.Sector==2");
hists.Hist1D("P.P",200,1,11,"P.Pid==11&&P.Sector==3");
hists.Hist1D("P.P",200,1,11,"P.Pid==11&&P.Sector==4");
hists.Hist1D("P.P",200,1,11,"P.Pid==11&&P.Sector==5");
hists.Hist1D("P.P",200,1,11,"P.Pid==11&&P.Sector==6")->Draw("(2x3)colz");

hists.Save("Mom.root");
}
