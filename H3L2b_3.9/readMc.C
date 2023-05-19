#include "tree.h"
#include "Hists.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TString.h"

#include <fstream>

// void readMc(TString mInputlist="Lambda_tree_mc.root", int const mode = 1,   TString outfile="fout_Lambda.root", int const mcState=1, int centL=0, int centH=8)
void readMc(TString mInputlist="quasiMC4.root", int const mode = 0,   TString outfile="fMC_H3L_0050.root", int const mcState=1, int centL=3, int centH=8)
{
  bool fillQAplots=1;
  TStopwatch time;
  time.Start();
  int snn = 390;
  double ycm;
  if(snn==3) //target
  {
    ycm = -1.045;
  }
  if(snn==390)
  {
    ycm=-1.37;
  }
  else{
    ycm = -999.;
  }
  double const mass_ht = 2.99131;
  double const ht_width = 0.005;
  double const mass_hl = 3.9239;
  double const hl_width = 0.005;
  double const mass_ld = 1.11568;
  double const mass_p = 0.93827;
  double const mass_pi = 0.13957;
  double const mass_d = 1.8756;


  //weighting function for MC
  if (mcState!=0){
    cout <<"starting book mc functions!" << endl;
    double par[4][3]={ {0.261677, 53.1192, 2.99131}, {0.261677, 53.1192, 2.99131}, {0.181327, 14641, 2.99131}, {0.181327, 14641, 2.99131}};
    double par_0_80[4][3]={ {0.261677, 53.1192, 2.99131}, {0.261677, 53.1192, 2.99131}, {0.181327, 14641, 2.99131}, {0.181327, 14641, 2.99131}};
    double par_10_40[4][3]={ {0.315791, 1.27673, 2.99131}, { 0.315791, 1.27673, 2.99131}, { 0.167319, 24701.3, 2.99131}, {  0.167319, 24701.3, 2.99131}};
    // TH1F* hpad = new TH1F("hpad","hpad",1, 0,10);
    // hpad->GetXaxis()->SetRangeUser(1e-10,1);
    // hpad->Draw();
    for (int irap=0;irap<4;irap++) {
      fH3Ldydpt[irap] = new TF1(Form("fH3Ldydpt[%d]",irap), "2*TMath::Pi()*[1]*x*exp(-(sqrt([2]*[2]+x*x))/[0])", 0,5);
      fH3Ldydpt[irap]->SetParameters(par[irap]);
      if (centL==4 && centH==6) 
      {  
        fH3Ldydpt[irap]->SetParameters(par_10_40[irap]);
      }
    }

    ////////uniform  mc pt y distribution/////////
    TFile* fgpt_0;
    if (mode==0 && mcState==1)
    {
      fgpt_0= new TFile("fMC_H3L_wt.root","READ");
      g_pt_fine_in = (TH2F*)fgpt_0->Get("hPhase")->Clone("g_pt_fine_in");
      gJam[0] = new TGraph("model/H3L2b_0_10cent.csv");
      gJam[1] = new TGraph("model/H3L2b_10_50cent.csv");
    }

    cout <<"finish book mc weighting functions." << endl;
  }
  ///////////////////////end of MC weighting////////////////

  cout <<"start read tree!" << endl;
  TString treename;
  treename = "htriton_mc_tree";
  TChain htriton_tree(treename.Data()); 

  TH1F* hvtx  = new TH1F("hvtx",    "Vz;Vz(cm);Counts",100,195,205);
  TH1F* hvtxgood  = new TH1F("hvtxgood","Vz;Vz(cm);Counts",100,195,205);
  TH1F* hrefmult  = new TH1F("hrefmult_tot", "refmult; hrefmult; N_{evt}", 600,0,600);
  TH1F* hcent = new TH1F("hcent", "cent; centrality; N_{evt}", 9,-0.5,8.5);
  TH1F* hcentwt = new TH1F("hcentwt", "centwt; centrality; N_{evt}", 9,-0.5,8.5);
  TH2F* hvtx_xy = new TH2F("hvtx_xy",  "Vx;Vx(cm);Vy;Vy(cm)",250,-5,5,250,-5,5);
  // hrefmult->SetDirectory(0);

  if (mInputlist.Contains(".root"))
  {
    htriton_tree.Add(mInputlist.Data());
    TFile *ftmp = new TFile(mInputlist.Data());
    TH1F* htmp3 = (TH1F*)ftmp->Get("hrefmult");
    hrefmult->Add(htmp3);
    // ftmp->Close();
    TH1F* htmp = (TH1F*)ftmp->Get("hCent");
    TH1F* htmp2 = (TH1F*)ftmp->Get("hCentWt");
    hcent->Add(htmp);
    hcentwt->Add(htmp2);
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
        htriton_tree.Add(tmp);
        TH1F* htmp = (TH1F*)ftmp->Get("hCent");
        TH1F* htmp2 = (TH1F*)ftmp->Get("hCentWt");
        hcent->Add(htmp);
        hcentwt->Add(htmp2);
        TH1F* htmp3 = (TH1F*)ftmp->Get("hrefmult");
        hrefmult->Add(htmp3);
        TH1F* htmp4 = (TH1F*)ftmp->Get("hvtx");
        hvtx->Add(htmp4);
        TH1F* htmp5 = (TH1F*)ftmp->Get("hvtxgood");
        hvtxgood->Add(htmp5);
        TH2F* htmp6 = (TH2F*)ftmp->Get("hvtx_xy");
        hvtx_xy->Add(htmp6);
        nfile++;
        ftmp->Close();
      }
    }
  }


  htriton_tree.SetBranchAddress("bmcparticleid",&bparticleid);

  htriton_tree.SetBranchAddress("reweight", &reweight);
  htriton_tree.SetBranchAddress("cent9", &cent9);

  htriton_tree.SetBranchAddress("bmcrawpx",&bmcpx);
  htriton_tree.SetBranchAddress("bmcrawpy",&bmcpy);   
  htriton_tree.SetBranchAddress("bmcrawpz",&bmcpz);
  htriton_tree.SetBranchAddress("bmcrawl",&bmcl);
  htriton_tree.SetBranchAddress("bmcrawpl",&bmcpl);
  // do not stored 
  // htriton_tree.SetBranchAddress("b0mcrawpx",&b0mcpx);
  // htriton_tree.SetBranchAddress("b0mcrawpy",&b0mcpy);                 
  // htriton_tree.SetBranchAddress("b0mcrawpz",&b0mcpz);
  //
  // htriton_tree.SetBranchAddress("b1mcrawpx",&b1mcpx);
  // htriton_tree.SetBranchAddress("b1mcrawpy",&b1mcpy);                 
  // htriton_tree.SetBranchAddress("b1mcrawpz",&b1mcpz); 
  //
  TH2F* hptH3Lmass  = new TH2F("hptH3Lmass","hptH3Lmass;p_{T};H3L mass",100,0,5,200,2.95,3.05);
  hptH3Lmass->Sumw2();
  TH3F* hH3LMassPtY= new TH3F("hH3LMassPtY","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  hH3LMassPtY->Sumw2();
  TH3F* hH3LMassPtCent= new TH3F("hH3LMassPtCent","hH3LMassPtCent;p_{T};H3L mass;Centrality",100,0,5,200,2.95,3.05, 9, -0.5,8.5);
  hH3LMassPtCent->Sumw2();

  TH1F* hdecayLMc = new TH1F("hdecayLMc", "hdecayLMc", 100, 0, 100);
  TH1F* hdecayLMcWt = new TH1F("hdecayLMcWt", "hdecayLMc", 100, 0, 100);
  TH1F* hplMcWt= new TH1F("hplMcWt", "hplMcWt", 100, 0, 100);

  TH2F* hPhase = new TH2F("hPhase", "hPhase;y;pt", 1200,-2.5,1, 1000,0,5 );
  TH2F* hPhase_wt = new TH2F("hPhase_wt", "hPhase_wt;y;pt", 120,-2.5,1, 100,0,5 );
  ////////////////////////////////////////////////////////////////

  Long64_t n_Entries = htriton_tree.GetEntries();
  cout <<"start process "<< n_Entries<<" events" << endl;

  for (int i=0;i<n_Entries;i++)
  {
    htriton_tree.GetEntry(i); 
    if (i%100000000==0) cout <<"read "<<i<<" events!" << endl;
    // if (cent9<=2 || cent9>8) continue; //remove 60-80%
    if (cent9<centL || cent9>centH) continue; //remove 60-80%
    // cout <<"test1" << endl;
    // cout<< cent9<<" "<<bparticlemass << endl;
    if (mode==0) {
      // if (bparticlemass<2.95 || bparticlemass > 3.05) continue;
      if (bparticleid!=3004) continue;
    }
    else if (mode==1) 
    {
      if (bparticlemass<3.86 || bparticlemass > 3.96) continue;
    }

    double ptweight = 1; //reserved
    double rapweight = 1;

    double mcweight = 1;
    if (mcState!=0) {
      TLorentzVector mcptc;
      if (mode==0 && mcState==1 ) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ht);
      else if (mode==1 && mcState==1) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_hl);
      double bmcrap = -1* (mcptc.Rapidity() - ycm);
      double mcMotherPt = mcptc.Pt();
      double mccounts = g_pt_fine_in->GetBinContent(g_pt_fine_in->GetXaxis()->FindBin(bmcrap),g_pt_fine_in->GetYaxis()->FindBin(mcMotherPt));
      if (mccounts>0) mcweight = 1./mccounts; 
      if (cent9>=7) rapweight = gJam[0]->Eval(fabs(bmcrap));
      else if (cent9<7) rapweight = gJam[1]->Eval(fabs(bmcrap));
      if (rapweight<0) rapweight=0;
      if (bmcrap>=0 && bmcrap<=0.5) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
      if (bmcrap<0 && bmcrap>=-0.25) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
      if (bmcrap<-0.25 && bmcrap>=-0.5) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
      if (bmcrap<-0.5 && bmcrap>=-0.75) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
      if (bmcrap<-0.75) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
    }

    // cout <<"test2" << endl;
    double gweight=1; //centrality weight

    // bmcpl = bmcl/(( sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz)/mass_ht));
    // double wt_l = TMath::Exp(-1*bmcpl*1e2/2.998*( 1./223. - 1./263.)); 
    // cout <<bmcpl <<" "<<bmcl/(( sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz)/mass_ht))<< endl;
    // wt_l = 1; 
    // double weight = ptweight*rapweight*mcweight*gweight*wt_l;
    // double weight = 1;
    double weight = ptweight*rapweight*mcweight*gweight;
    double weight2 = ptweight*rapweight*mcweight*gweight;

    //compare H3L or  background
    if (mode==0  ) {
      TLorentzVector H3L;
      H3L.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ht);
      double H3LpT = sqrt(bmcpx*bmcpx+bmcpy*bmcpy);
      double H3Ly = -1*(H3L.Rapidity()-ycm);
      hptH3Lmass->Fill(H3LpT, mass_ht, weight); 
      hH3LMassPtY->Fill(H3LpT, mass_ht, H3Ly, weight);
      if (H3LpT<1.5 && H3LpT>1.) hdecayLMc->Fill( bmcl, weight2);
      if (H3LpT<1.5 && H3LpT>1.) hdecayLMcWt->Fill(bmcl, weight);
      hplMcWt->Fill(bmcpl, weight);
      if (fabs(H3Ly)<1) hH3LMassPtCent->Fill(H3LpT, mass_ht, cent9, weight);
      // hPhase->Fill( H3Ly, H3LpT);
      hPhase->Fill( H3L.Rapidity(), H3LpT);
      hPhase_wt->Fill( H3Ly, H3LpT, weight);

      // cout <<H3LpT<<" "<< weight << endl;
    }
  }

  TFile* fout = new TFile(outfile.Data(),"recreate");

  ///////////////////
  if (mode==0)
  {
    // hH3LptPionPt->Write();
    // hH3LptProtonPt->Write();
    // hH3LPPiMassPt->Write();
    hH3LMassPtY->Write();
    hH3LMassPtCent->Write();
    hPhase->Write();
    hPhase_wt->Write();
    hdecayLMcWt->Write();
    hdecayLMc->Write();
    hplMcWt->Write();
  }

  // hrefmult->Write();
  hcentwt->Write();
  hcent->Write();

  fout->Close();
  time.Stop();
  time.Print();
}
