#include <TFile.h>
#include <TTree.h>
#include "functions.h"
Int_t Centrality(int gRefMult )
{
  int centrality;
  // int centFull[10]={5, 9, 16, 26, 41, 60, 86, 119, 142, 195};
  int centFull[10]={6, 12, 28, 36, 56, 84, 120, 167, 198, 290};
  if (gRefMult>=centFull[8]) centrality=8;
  else if (gRefMult>=centFull[7] && gRefMult<centFull[8] ) centrality=7;
  else if (gRefMult>=centFull[6] && gRefMult<centFull[7] ) centrality=6;
  else if (gRefMult>=centFull[5] && gRefMult<centFull[6] ) centrality=5;
  else if (gRefMult>=centFull[4] && gRefMult<centFull[5] ) centrality=4;
  else if (gRefMult>=centFull[3] && gRefMult<centFull[4] ) centrality=3;
  else if (gRefMult>=centFull[2] && gRefMult<centFull[3] ) centrality=2;
  else if (gRefMult>=centFull[1] && gRefMult<centFull[2] ) centrality=1;
  else if (gRefMult>=centFull[0] && gRefMult<centFull[1] ) centrality=0;
  else centrality = 9;

  return centrality;
}
void readtree(TString mInputlist, TString outfilename) {
  TString treename="EmbeddingTracks";
  TChain* tree = new TChain(treename.Data());

  // Get the TTree from the file
  // TTree* tree = dynamic_cast<TTree*>(file->Get("singleparticle_tree"));
  if (mInputlist.Contains(".root"))
  {
    tree->Add(mInputlist.Data());
  }
  else
  {
    int nfile = 0;
    char tmp[2000];
    ifstream readlists;
    readlists.open(mInputlist.Data());
    while (readlists.good()){
      readlists.getline(tmp,2000);
      TFile *ftmp = new TFile(tmp);
      if (!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
        cout<<"Could not open this file: "<< tmp  <<endl;
        continue;
      }
      else {
        if(nfile%1==0) cout<<"read in "<<nfile<<"th file: "<< tmp <<endl;
        tree->Add(tmp);
        nfile++;
        ftmp->Close();
      }
    }
  }
  // Define the variables to store the branch data
  Float_t runid, eventid, gref, ref, vpdVz, ispion, iskaon, isproton, isdeuteron, istriton, isHe3, cent, mcGeantId, mcHitsTpc, mcIdVx, rcId, rcCharge, rcNfit, rcNposs, rcNdedx, rcNsigPi, rcNsigK;
  Float_t vtx, vty, vtz, mcPt, mcP, mcEta, mcY, mcPhi, rcPt, rcP, rcEta, rcPhi, rcY, rcDedx, rcDca, rcDcaXy, rcDcaZ;
  Float_t map0, map1, map2;

  // Set the branch addresses
  tree->SetBranchAddress("runid", &runid);
  tree->SetBranchAddress("eventid", &eventid);
  tree->SetBranchAddress("gref", &gref);
  tree->SetBranchAddress("ref", &ref);
  tree->SetBranchAddress("vtx", &vtx);
  tree->SetBranchAddress("vty", &vty);
  tree->SetBranchAddress("vtz", &vtz);
  tree->SetBranchAddress("vpdVz", &vpdVz);
  tree->SetBranchAddress("ispion", &ispion);
  tree->SetBranchAddress("iskaon", &iskaon);
  tree->SetBranchAddress("isproton", &isproton);
  tree->SetBranchAddress("isdeuteron", &isdeuteron);
  tree->SetBranchAddress("istriton", &istriton);
  tree->SetBranchAddress("isHe3", &isHe3);
  tree->SetBranchAddress("cent", &cent);
  tree->SetBranchAddress("mcPt", &mcPt);
  tree->SetBranchAddress("mcP", &mcP);
  tree->SetBranchAddress("mcEta", &mcEta);
  tree->SetBranchAddress("mcY", &mcY);
  tree->SetBranchAddress("mcPhi", &mcPhi);
  tree->SetBranchAddress("mcGeantId", &mcGeantId);
  tree->SetBranchAddress("mcHitsTpc", &mcHitsTpc);
  tree->SetBranchAddress("mcIdVx", &mcIdVx);
  tree->SetBranchAddress("rcId", &rcId);
  tree->SetBranchAddress("rcPt", &rcPt);
  tree->SetBranchAddress("rcP", &rcP);
  tree->SetBranchAddress("rcEta", &rcEta);
  tree->SetBranchAddress("rcPhi", &rcPhi);
  tree->SetBranchAddress("rcY", &rcY);
  tree->SetBranchAddress("rcCharge", &rcCharge);
  tree->SetBranchAddress("rcNfit", &rcNfit);
  tree->SetBranchAddress("rcNposs", &rcNposs);
  tree->SetBranchAddress("rcNdedx", &rcNdedx);
  tree->SetBranchAddress("rcDedx", &rcDedx);
  tree->SetBranchAddress("rcNsigPi", &rcNsigPi);
  tree->SetBranchAddress("rcNsigK", &rcNsigK);
  tree->SetBranchAddress("rcDca", &rcDca);

  bookHists();

  // Loop over the entries in the tree
  for (Long64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry) {
    // Get the current entry in the tree
    tree->GetEntry(iEntry);
    if (iEntry%1000000==0) cout << iEntry/(1.0*tree->GetEntries())<< endl;
    // int cent = getCentrality(countrefmult);
    int particleid=-999;
    particleid = ConvertGeantToPDG(mcGeantId);
    double weight =1; 

    if (rcId>0)  { 
      bool passquality = false;
      passquality = rcDca<3 && rcNfit>15 && rcNfit/(rcNposs*1.)>0.52;

      if (rcNfit>20 && rcNdedx>10 )  hrc->Fill(rcPt,rcEta,cent,weight);
      if (passquality) {
        // if (rcNfit>20 && rcNdedx/(1.*rcNposs>0.4))  hrc->Fill(rcPt,rcEta,cent,weight);
        fillQaHistograms(particleid,rcPt,rcEta,rcNfit,rcDca,cent,weight);
        if ( rcNfit>20)hphipt->Fill(rcPt,rcPhi,weight);
        if (rcNdedx>10 && rcNfit>20) hphipt_ndedx->Fill(rcPt,rcPhi,weight);
      }
    }
   
    hmc->Fill(mcPt,mcEta,cent,weight); 

  }

  // Close the ROOT file
  Write(outfilename);
}

