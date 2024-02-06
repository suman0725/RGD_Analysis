
void PhiThetalable() {
    TFile *file = new TFile("PhiTheta.root", "READ");

// Create two TLists, one for each histogram
 TList *histList_phi = dynamic_cast<TList*>(file->Get("HipoHists_phi"));
    TList *histList_theta = dynamic_cast<TList*>(file->Get("HipoHists_theta"));
    
