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

#define DEBUG 0

void readtree(TString mInputlist="Lambda_tree_mc.root", int const mode = 1,   TString outfile="fout_Lambda.root", int const mcState=1, int centLow=3, int centHigh=8)
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
    double par[4][3]={{0.261677, 53.1192, 2.99131}, {0.261677, 53.1192, 2.99131}, {0.181327, 14641, 2.99131}, {0.181327, 14641, 2.99131}};
    for (int irap=0;irap<4;irap++) {
      fH3Ldydpt[irap] = new TF1(Form("fH3Ldydpt[%d]",irap), "2*TMath::Pi()*[1]*x*exp(-(sqrt([2]*[2]+x*x))/[0])", 0,5);
      fH3Ldydpt[irap]->SetParameters(par[irap]);
      // fH3Ldydpt[irap]->Draw("same");
      // fH3Ldydpt[irap]->SetLineColor(irap+1);
    }
    // return;

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
  treename = "htriton_tree";
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

  htriton_tree.SetBranchAddress("brunid", &brunid);
  htriton_tree.SetBranchAddress("beventid", &beventid);
  htriton_tree.SetBranchAddress("bVz", &bVz);
  htriton_tree.SetBranchAddress("bVzerr", &bVzerr);
  htriton_tree.SetBranchAddress("brefmult",&brefmult);
  htriton_tree.SetBranchAddress("btofmult",&btofmult);
  htriton_tree.SetBranchAddress("bparticleid",&bparticleid);
  htriton_tree.SetBranchAddress("bparticlemass",&bparticlemass);
  htriton_tree.SetBranchAddress("bpx",&bpx);
  htriton_tree.SetBranchAddress("bpy",&bpy);
  htriton_tree.SetBranchAddress("bpz",&bpz);
  htriton_tree.SetBranchAddress("chi2primary_he", &chi2primary_he);
  htriton_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
  htriton_tree.SetBranchAddress("ht_chi2topo", &ht_chi2topo);
  htriton_tree.SetBranchAddress("ht_chi2ndf", &ht_chi2ndf);
  htriton_tree.SetBranchAddress("ht_ldl", &ht_ldl);
  htriton_tree.SetBranchAddress("ht_l", &ht_l);
  htriton_tree.SetBranchAddress("ht_dl", &ht_dl);
  htriton_tree.SetBranchAddress("dca_he",&dca_he);
  htriton_tree.SetBranchAddress("dca_pi",&dca_pi);

  htriton_tree.SetBranchAddress("px_pi",&px_pi);
  htriton_tree.SetBranchAddress("py_pi",&py_pi);
  htriton_tree.SetBranchAddress("pz_pi",&pz_pi);
  htriton_tree.SetBranchAddress("px_he",&px_he);
  htriton_tree.SetBranchAddress("py_he",&py_he);
  htriton_tree.SetBranchAddress("pz_he",&pz_he);
  htriton_tree.SetBranchAddress("dedx_he",&dedx_he);
  htriton_tree.SetBranchAddress("dedx_pi",&dedx_pi);
  htriton_tree.SetBranchAddress("nhitsdedx_pi",&nhitsdedx_pi);
  htriton_tree.SetBranchAddress("nhitsdedx_he",&nhitsdedx_he);
  htriton_tree.SetBranchAddress("ht_dl", &ht_dl);
  // htriton_tree.SetBranchAddress("dca_he",&dca_he);
  // htriton_tree.SetBranchAddress("dca_pi",&dca_pi);

  htriton_tree.SetBranchAddress("px_pi",&px_pi);
  htriton_tree.SetBranchAddress("py_pi",&py_pi);
  htriton_tree.SetBranchAddress("pz_pi",&pz_pi);
  htriton_tree.SetBranchAddress("px_he",&px_he);
  htriton_tree.SetBranchAddress("py_he",&py_he);
  htriton_tree.SetBranchAddress("pz_he",&pz_he);
  htriton_tree.SetBranchAddress("dedx_he",&dedx_he);
  htriton_tree.SetBranchAddress("dedx_pi",&dedx_pi);
  htriton_tree.SetBranchAddress("nhitsdedx_pi",&nhitsdedx_pi);
  htriton_tree.SetBranchAddress("nhitsdedx_he",&nhitsdedx_he);
  htriton_tree.SetBranchAddress("nhits_pi",&nhits_pi);
  htriton_tree.SetBranchAddress("nhits_he",&nhits_he);
  htriton_tree.SetBranchAddress("ht_bdfvtx",&ht_bdfvtx);
  htriton_tree.SetBranchAddress("ht_bdfvtx2",&ht_bdfvtx2);
  htriton_tree.SetBranchAddress("ht_lifetime",&ht_lifetime);
  htriton_tree.SetBranchAddress("countrefmult",&countrefmult);
  htriton_tree.SetBranchAddress("reweight", &reweight);
  htriton_tree.SetBranchAddress("cent9", &cent9);
  htriton_tree.SetBranchAddress("bismc", &bismc);
  htriton_tree.SetBranchAddress("bpl",&bpl);
  if(mcState==1){
    // htriton_tree.SetBranchAddress("bVzerr",&bVzerr);
    htriton_tree.SetBranchAddress("bmcpx",&bmcpx);
    htriton_tree.SetBranchAddress("bmcpy",&bmcpy);
    htriton_tree.SetBranchAddress("bmcpz",&bmcpz);
    htriton_tree.SetBranchAddress("bmcl",&bmcl);
    htriton_tree.SetBranchAddress("bmcpl",&bmcpl);
  }
  htriton_tree.SetBranchAddress("FXTMult2",&FXTMult2);

  cout <<"finish read tree" << endl;

  TH1F* hcharge= new TH1F("hcharge","hcharge",2,-1.5,1.5);

  double masslo, masshi;
  if (mode==0) {masslo=2.95; masshi=3.05;}
  if (mode==1) {masslo=3.86; masshi=3.96;}

  double pth, ptl;


  TH2F* hptH3Lmass  = new TH2F("hptH3Lmass","hptH3Lmass;p_{T};H3L mass",100,0,5,200,masslo,masshi);
  hptH3Lmass->Sumw2();
  TH3F* hH3LMassPtY= new TH3F("hH3LMassPtY","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,masslo,masshi, 300, -2,1.);
  hH3LMassPtY->Sumw2();

  TH3F* hH3LMassPtCent= new TH3F("hH3LMassPtCent","hH3LMassPtCent;p_{T};H3L mass;Centrality",100,0,5,200,masslo,masshi, 9, -0.5,8.5);
  hH3LMassPtCent->Sumw2();

  //topological variable for H3L
  //signal region
  TH3F* h3H3L_lSig= new TH3F("h3H3L_lSig","hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 300, -2., 1);
  h3H3L_lSig->Sumw2();

  TH3F* h3H3L_ldlSig = new TH3F("h3H3L_ldlSig","hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,300, -2., 1);
  h3H3L_ldlSig->Sumw2();


  TH3F* h3H3L_chi2ndfSig= new TH3F("h3H3L_chi2ndfSig","hptH3L_chi2ndf;p_{T};chi2ndf;y",100,0,5,100,0,10, 300, -2., 1);
  h3H3L_chi2ndfSig->Sumw2();
  TH3F* h3H3L_chi2topoSig= new TH3F("h3H3L_chi2topoSig","hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,100,0,10, 300, -2., 1);
  h3H3L_chi2topoSig->Sumw2();

  TH3F* hHeDcaPtY= new TH3F("hHeDcaPtY", "hHeDcaPtY;p_{T} (GeV/c);Dca;y",100, 0, 5,  200, 0, 5, 300, -2,1);
  hHeDcaPtY->Sumw2();
  TH3F* hPiDcaPtY= new TH3F("hPiDcaPtY", "hPiDcaPtY;p_{T} (GeV/c);Dca;y", 100, 0, 5, 200, 0, 5, 300, -2,1);
  hPiDcaPtY->Sumw2();

  TH3F* hHechi2primPtY= new TH3F("hHechi2primPtY", "hHechi2primPtY;p_{T} (GeV/c);chi2prim;y",100, 0, 5,  200, 0, 200, 300, -2,1);
  hHechi2primPtY->Sumw2();
  TH3F* hPichi2primPtY= new TH3F("hPichi2primPtY", "hPichi2primPtY;p_{T} (GeV/c);chi2prim;y", 100, 0, 5, 200, 0, 200, 300, -2,1);
  hPichi2primPtY->Sumw2();
  TH3F* hHenHitsPtY = new TH3F("hHenHitsPtY", "hHenHitsPtY;p_{T} (GeV/c);nHits;#eta",100, 0, 5,  100, 0, 100, 300, -2,1);
  hHenHitsPtY->Sumw2();
  TH3F* hPinHitsPtY = new TH3F("hPinHitsPtY", "hPinHitsPtY;p_{T} (GeV/c);nHits;#eta", 100, 0, 5, 100, 0, 100, 300, -2,1);
  hPinHitsPtY->Sumw2();

  TH3F* h3H3L_lMass = new TH3F("h3H3L_lMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  h3H3L_lMass->Sumw2();
  TH3F* h3H3L_ldlMass = new TH3F("h3H3L_ldlMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  h3H3L_ldlMass->Sumw2();
  TH3F* h3H3L_chi2ndfMass = new TH3F("h3H3L_chi2ndfMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  h3H3L_chi2ndfMass->Sumw2();
  TH3F* h3H3L_chi2topoMass = new TH3F("h3H3L_chi2topoMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  h3H3L_chi2topoMass->Sumw2();
  TH3F* hHeDcaMass = new TH3F("hHeDcaMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  hHeDcaMass->Sumw2();
  TH3F* hPiDcaMass = new TH3F("hPiDcaMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  hPiDcaMass->Sumw2();
  TH3F* hHenHitsMass = new TH3F("hHenHitsMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  hHenHitsMass->Sumw2();
  TH3F* hPinHitsMass = new TH3F("hPinHitsMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  hPinHitsMass->Sumw2();
  TH3F* hHechi2primMass = new TH3F("hHechi2primMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  hHechi2primMass->Sumw2();
  TH3F* hPichi2primMass = new TH3F("hPichi2primMass","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
  hPichi2primMass->Sumw2();

  int const nsys=16;
  TH3F* hH3LMassPtYSysCut[nsys];
  for (int i=0;i<nsys;i++)
  {
    hH3LMassPtYSysCut[i] = new TH3F(Form("hH3LMassPtYSysCut%d",i),"hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 300, -2.,1);
    hH3LMassPtYSysCut[i]->Sumw2();
  }
  ////////////////////////////////////////////////////////////////

  Long64_t n_Entries = htriton_tree.GetEntries();
  cout <<"start process "<< n_Entries<<" events" << endl;

  for (int i=0;i<n_Entries;i++)
  {
    htriton_tree.GetEntry(i); 
    if (i%100000==0) cout <<"read "<<i<<" events!" << endl;
    if ( !(bismc == mcState)  ) continue;
    if (cent9<centLow || cent9>centHigh) continue; // since no official centraltiy definition, so remove cut

    // if (DEBUG) 
    //   cout <<bparticlemass<<" "<<weight<<endl;

    double ptweight = 1; //reserved
    double rapweight = 1;
    double mcweight = 1;
    if (mcState!=0) {
      TLorentzVector mcptc;
      if (mode==0 && mcState==1 ) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ht);
      else if (mode==1 && mcState==1) {
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_hl);
      }

      //mcMotherPt =sqrt(bmcpx*bmcpx+bmcpy*bmcpy) ;
      double bmcrap = -1*(mcptc.Rapidity() - ycm);
      double mcMotherPt = mcptc.Pt();
      // rapweight = t_quadr->Eval(bmcrap);
      if (cent9>=7) rapweight = gJam[0]->Eval(fabs(bmcrap));
      else if (cent9<7) rapweight = gJam[1]->Eval(fabs(bmcrap));
      double mccounts=0;

      if (mode==0 && mcState==1) 
        mccounts = g_pt_fine_in->GetBinContent(g_pt_fine_in->GetXaxis()->FindBin(bmcrap),g_pt_fine_in->GetYaxis()->FindBin(mcMotherPt));
      else if (mode==1 && mcState==1)
      {  
        // To Be added
      }
      if (mccounts>0) mcweight = 1./mccounts;

      // if (mode==0 && mcState ==1) ptweight = bolt1->Eval(mcMotherPt);
      if (mode==0 && mcState ==1) {
        if (bmcrap>0 && bmcrap<=0.5) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
        if (bmcrap<0 && bmcrap>=-0.25) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
        if (bmcrap<-0.25 && bmcrap>=-0.5) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
        if (bmcrap<-0.5 && bmcrap>=-0.75) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
        if (bmcrap<-0.75) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
      }
    }
    // cout << mcMotherPt <<" "<<mcweight<<" "<<ptweight<<" "<<rapweight<<endl;
    double gweight=1; //centrality weight 

    // double wt_l = TMath::Exp(-1*bmcpl*1e2/2.998*( 1./223. - 1./263.)); 
    // double weight = ptweight*rapweight*mcweight*gweight*wt_l;
    // double weight = 1;
    double weight = ptweight*rapweight*mcweight*gweight;
    if (mcState==0) weight = gweight; 
    if (mode==0) {
      if (bparticlemass<2.95 || bparticlemass > 3.05) continue;
      if (bparticleid!=3004) continue;
    }
    else if (mode==1) 
    {
      if (bparticlemass<3.86 || bparticlemass > 3.96) continue;
    }

    if (DEBUG) 
      cout <<bparticlemass<<" "<< weight<< endl;

    hcharge->Fill(bparticleid>0?1:-1);
    if (bparticleid<0) continue; //currently only look at particle

    //compare H3L or  background
    if (mode==0 || mode==1 ) {
      TLorentzVector H3L;
      H3L.SetXYZM(bpx, bpy, bpz, bparticlemass );
      double H3LpT = sqrt(bpx*bpx+bpy*bpy);
      double H3Ly = -1*(H3L.Rapidity() - ycm);
      double pi_pt = sqrt(px_pi*px_pi+py_pi*py_pi);
      double he_pt = sqrt(px_he*px_he+py_he*py_he);
      double p_pi = sqrt(px_pi*px_pi+py_pi*py_pi+pz_pi*pz_pi);
      double p_He = sqrt(px_he*px_he+py_he*py_he+pz_he*pz_he);

      bool dcaCut = pi_pt>0.1 && he_pt>0.1;
      bool nhitscut = nhits_he>20 && nhits_pi>20;
      //default cut
      double lcut=1, ldlcut=3, chi2topocut = 4., chi2ndfcut=3.5, chi2_hecut=0, chi2_picut=10, pHecutLow=0, pHecut=50, ppicut=50, bplcut=0;
      bool PIDcut = (fabs(p_He)<pHecut && fabs(p_pi)<ppicut);
      if (H3Ly<-0.5) chi2_hecut=5;
      if (centLow>=7) { chi2topocut = 4.5; chi2_picut=7;} 

      bool passTopoCuts = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && ht_chi2topo>0 &&  fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && ht_chi2ndf>0  && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhitscut;

      int const nmax = 20;
      bool passTopoCutsSys[nmax];
      // double lcut1 = 1, lcut2 = 1, ldlcut1 = 3, ldlcut2 = 3; // very loose cuts
      double  chi2topocut1 = 3.5, chi2topocut2 = 4.5;
      if (centLow>=7) { chi2topocut1 = 4; chi2topocut2 = 4.5;} 
      double  chi2ndfcut1 = 2.5, chi2ndfcut2 = 4;
      double  chi2_picut1= 7, chi2_picut2 = 13;
       if (centLow>=7) { chi2_picut1 = 4; chi2_picut2 = 10;}
      double  nhitscut1= 15, nhitscut2 = 25;

      int ncuts=0; 
      // passTopoCutsSys[ncuts++] = ht_l >lcut1 && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >lcut2 && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut1 && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut2 && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut1  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhitscut;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut2  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut1 && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut2 && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut1 && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut2 && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && ht_chi2topo>0 &&  fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && ht_chi2ndf>0  && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhits_he>nhitscut1 && nhits_pi>nhitscut1;
      passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && ht_chi2topo>0 &&  fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && ht_chi2ndf>0  && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhits_he>nhitscut2 && nhits_pi>nhitscut2;

      // //below is for tune cuts
      // passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && ht_chi2topo>0 &&  fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && ht_chi2ndf>0  && chi2primary_pi>chi2_picut && chi2primary_he>2 && bpl>bplcut && dcaCut && nhits_he>20 && nhits_pi>20;
      // passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && ht_chi2topo>0 &&  fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && ht_chi2ndf>0  && chi2primary_pi>chi2_picut && chi2primary_he>5 && bpl>bplcut && dcaCut && nhits_he>20 && nhits_pi>20;
      // passTopoCutsSys[ncuts++] = ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && ht_chi2topo>0 &&  fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && ht_chi2ndf>0  && chi2primary_pi>chi2_picut && chi2primary_he>8 && bpl>bplcut && dcaCut && nhits_he>20 && nhits_pi>20;
      // passTopoCutsSys[ncuts++] = ht_l >1 && ht_ldl>3 && ht_chi2topo<3  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<3.5 && chi2primary_pi>0 && chi2primary_he>0 && bpl>0 && dcaCut  && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >1 && ht_ldl>3 && ht_chi2topo<3  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<3.5 && chi2primary_pi>3 && chi2primary_he>0 && bpl>0 && dcaCut  && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >1 && ht_ldl>3 && ht_chi2topo<2.5  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<3.5 && chi2primary_pi>5 && chi2primary_he>0 && bpl>0 && dcaCut  && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >1 && ht_ldl>3 && ht_chi2topo<2.5  && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<3.5 && chi2primary_pi>3 && chi2primary_he>0 && bpl>0 && dcaCut  && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >1 && ht_ldl>3 && ht_chi2topo<3 && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<3.5 && chi2primary_pi>5 && chi2primary_he>0 && bpl>0 && dcaCut  && nhitscut;
      // passTopoCutsSys[ncuts++] = ht_l >1 && ht_ldl>3 && ht_chi2topo<3 && fabs(p_He)>pHecutLow && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<2.5 && chi2primary_pi>40 && chi2primary_he>10  && dcaCut  && nhitscut;

      // topo qa cuts
      bool qaldlcuts= ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhitscut;
      bool qalcuts=ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      bool qachi2topocuts=ht_l >lcut && ht_ldl>ldlcut  && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      bool qachi2ndfcuts=ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      bool qapichi2primcuts=ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      bool qapinhitscuts=ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhits_he>20; 
      bool qapiDcacuts=ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut  && nhitscut;
      bool qaHechi2primcuts=ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>0&& bpl>bplcut && dcaCut  && nhitscut;
      bool qaHenhitscuts=ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut &&nhits_pi>20;
      bool qaHeDcacuts= ht_l >lcut && ht_ldl>ldlcut && ht_chi2topo<chi2topocut && fabs(p_He)<pHecut && fabs(p_pi)<ppicut && ht_chi2ndf<chi2ndfcut && chi2primary_pi>chi2_picut && chi2primary_he>chi2_hecut && bpl>bplcut && dcaCut && nhitscut;

      if (passTopoCuts) {
        hH3LMassPtY->Fill(H3LpT, bparticlemass, H3Ly, weight);
        hH3LMassPtCent->Fill(H3LpT, bparticlemass, cent9, weight);
      }
      //add topo qa 
      if (bparticlemass>2.988 && bparticlemass<2.996){
  
        if (qalcuts)
        {
          h3H3L_lSig->Fill(H3LpT, ht_l, H3Ly, weight);
        }  
        if (qaldlcuts)
        {
          h3H3L_ldlSig->Fill(H3LpT, ht_ldl, H3Ly, weight); 
        }
        if (qachi2ndfcuts)
        {
          h3H3L_chi2ndfSig->Fill(H3LpT, ht_chi2ndf, H3Ly, weight);
        }
        if (qachi2topocuts)
        {
          h3H3L_chi2topoSig->Fill(H3LpT, ht_chi2topo, H3Ly, weight);
        }  
  
        if (qapichi2primcuts)
        {
          hPichi2primPtY->Fill(H3LpT, chi2primary_pi, H3Ly, weight);
        }
        if (qaHechi2primcuts)
        {
          hHechi2primPtY->Fill(H3LpT, chi2primary_he, H3Ly, weight);
        }
        if (qapiDcacuts)
        {
          hPiDcaPtY->Fill(H3LpT, dca_pi, H3Ly, weight);
        }
        if (qaHeDcacuts)
        {
          hHeDcaPtY->Fill(H3LpT, dca_he, H3Ly, weight);
        }
        if (qapinhitscuts)
        {
          hPinHitsPtY->Fill(H3LpT, nhits_pi, H3Ly, weight);
        }
        if (qaHenhitscuts)
        {
          hHenHitsPtY->Fill(H3LpT, nhits_he, H3Ly, weight);
        }
      }
      if (bparticlemass>2.95 && bparticlemass<3.04){
  
        if (qalcuts)
        {
          h3H3L_lMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }  
        if (qaldlcuts)
        {
          h3H3L_ldlMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
        if (qachi2ndfcuts)
        {
          h3H3L_chi2ndfMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
        if (qachi2topocuts)
        {
          h3H3L_chi2topoMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }  
  
        if (qapichi2primcuts)
        {
          hPichi2primMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
        if (qaHechi2primcuts)
        {
          hHechi2primMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
        if (qapiDcacuts)
        {
          hPichi2primMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
        if (qaHeDcacuts)
        {
          hHeDcaMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
        if (qapinhitscuts)
        {
          hPinHitsMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
        if (qaHenhitscuts)
        {
          hHenHitsMass->Fill(H3LpT, bparticlemass, H3Ly, weight);
        }
      }

      for (int isys=0;isys<nsys;isys++)
      {
        if (passTopoCutsSys[isys]) 
          hH3LMassPtYSysCut[isys]->Fill(H3LpT, bparticlemass, H3Ly, weight);
      }
    }
  }

  TFile* fout = new TFile(outfile.Data(),"recreate");

  hcharge->Write();

  if (mode==0)
  {
    if (fillQAplots){ 
      h3H3L_lSig->Write();
      h3H3L_ldlSig->Write();
      h3H3L_chi2ndfSig->Write();
      h3H3L_chi2topoSig->Write();
      hHeDcaPtY->Write();
      hPiDcaPtY->Write();
      hHechi2primPtY->Write();
      hPichi2primPtY->Write();
      hHenHitsPtY->Write();
      hPinHitsPtY->Write();

      h3H3L_lMass->Write();
      h3H3L_ldlMass->Write();
      h3H3L_chi2ndfMass->Write();
      h3H3L_chi2topoMass->Write();
      hHeDcaMass->Write();
      hPiDcaMass->Write();
      hHechi2primMass->Write();
      hPichi2primMass->Write();
      hHenHitsMass->Write();
      hPinHitsMass->Write();
    }

    hH3LMassPtY->Write();

    for (int i=0;i<nsys;i++) {
      hH3LMassPtYSysCut[i]->Write(); 
    }
    hH3LMassPtCent->Write();
    // hH3LMassPtY_5_40->Write();
  }

  hrefmult->Write();
  hcentwt->Write();
  hcent->Write();
  hvtx_xy->Write();
  hvtx->Write();
  hvtxgood->Write();

  fout->Close();
  time.Stop();
  time.Print();
}
