#include "sPhenixStyle.h"
void setHistStyle(TH1* h, int color, int markerstyle, double size , int mode = 0)
{
  h->SetLineColor(color);
  if (mode==0) //hist and graph
  {
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerstyle);
    h->SetMarkerSize(size);
  }
  else if (mode==1) //TF1
  {
    h->SetLineStyle(markerstyle);
    h->SetLineWidth(size);
  }
}
void Norm(TH1* h){
  if (h->Integral()!=0)
    h->Scale(1./h->Integral());
}
void drawplots(TString filename1,TString histname1, int charge1, TString legname1, 
    TString filename2, TString histname2, int charge2, TString legname2,TString particlename,TString outname)
{
  SetsPhenixStyle();
  // map<int,int> pididx={{211,0},{-211,1},{321,2},{-321,3},{2212,4},{-2212,5}};
  // TString particlename[]={"#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}"};
  // TString histname[]={"hpion","hpion","hkaon","hkaon","hproton","hproton"};
  // int idx = pididx[pid];

  // TString centname[9]={"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  const int ncent = 5;
  // TString centname[ncent]={"0-80%","60-80%","40-60%","20-40%","0-20%"};
  // int centH[ncent]={8,1,3,5,8};
  // int centL[ncent]={0,0,2,4,6};
  TString centname[ncent]={"0-80%","0-20%","20-40%","40-60%","60-80%"};
  int centH[ncent]={8,8,5,3,1};
  int centL[ncent]={0,6,4,2,0};

  TCanvas* c = new TCanvas("c","c",800,600);
  TPDF* pdf = new TPDF(Form("plots_%s_gDca.pdf",outname.Data()));
  pdf->Off();

  TFile* f1 = new TFile(filename1.Data()); 
  THnSparseF* hn1 = (THnSparseF*)f1->Get(histname1.Data()); 
  // hn1->SetDirectory(0);
  f1->Close();

  TFile* f2 = new TFile(filename2.Data()); 
  THnSparseF* hn2 = (THnSparseF*)f2->Get(histname2.Data()); 
  // hn2->SetDirectory(0);
  f2->Close();

  int a=3;int b=2;
  c->Divide(a,b);
  int changepage = a*b;
  int ipad=0;

  //projection dca
  for (int ic=0;ic<ncent;ic++) 
  {
    hn1->GetAxis(2)->SetRangeUser(0,3);
    hn2->GetAxis(2)->SetRangeUser(0,3);

    hn1->GetAxis(3)->SetRangeUser(20,100);
    hn2->GetAxis(3)->SetRangeUser(20,100);
    hn1->GetAxis(4)->SetRangeUser(centL[ic],centH[ic]);
    hn2->GetAxis(4)->SetRangeUser(centL[ic],centH[ic]);
    hn1->GetAxis(5)->SetRangeUser(charge1,charge1);
    hn2->GetAxis(5)->SetRangeUser(charge2,charge2);

    TH3F* h3dca1 = (TH3F*)hn1->Projection(0,1,2,"E");
    TH3F* h3dca2 = (TH3F*)hn2->Projection(0,1,2,"E");

    for (int ieta=0;ieta<20;ieta++) 
    {
      for (int ip=0;ip<15;ip++)
      { 
        c->cd(ipad+1);
        double pth=0.1*(ip+1)+0;
        double ptl=0.1*(ip)+0;
        double etal=0.1*ieta-2;
        double etah=0.1*(ieta+1)-2;

        if (etal>0. || etah<-1.6 || ptl>2.5) continue;
        // hn1->GetAxis(0)->SetRangeUser(ptl,pth);
        // hn1->GetAxis(1)->SetRangeUser(etal,etah);
        // hn1->GetAxis(4)->SetRangeUser(ic,ic);
        // hn1->GetAxis(5)->SetRangeUser(charge1,charge1);
        // TH1F* h1 = (TH1F*)hn1->Projection(2,"E");
        //
        // hn2->GetAxis(0)->SetRangeUser(ptl,pth);
        // hn2->GetAxis(1)->SetRangeUser(etal,etah);
        // hn2->GetAxis(4)->SetRangeUser(ic,ic);
        // hn2->GetAxis(5)->SetRangeUser(charge2,charge2);
        // TH1F* h2 = (TH1F*)hn2->Projection(2,"E");

        TH1D* h1 = (TH1D*)h3dca1->ProjectionZ("h1",h3dca1->GetXaxis()->FindBin(ptl),h3dca1->GetXaxis()->FindBin(pth),
            h3dca1->GetYaxis()->FindBin(etal),h3dca1->GetYaxis()->FindBin(etah));
        TH1D* h2 = (TH1D*)h3dca2->ProjectionZ("h2",h3dca2->GetXaxis()->FindBin(ptl),h3dca2->GetXaxis()->FindBin(pth),
            h3dca2->GetYaxis()->FindBin(etal),h3dca2->GetYaxis()->FindBin(etah));

        setHistStyle(h1, kRed-4, kOpenCircle, 0.8, 0); 
        setHistStyle(h2, kBlue, kFullCircle, 0.8, 0); 
        h1->GetXaxis()->SetTitle("gDca (cm)");
        h1->GetYaxis()->SetTitle("Arb. Unit");
        Norm(h1);Norm(h2);

        TH1D* hpad = new TH1D("h","h;gDca (cm);Arb. Unit",1,0,3);
        // hpad->GetYaxis()->SetRangeUser(0,0.2);
        double max2 = h2->GetMaximum();
        double max1 = h1->GetMaximum();

        hpad->GetYaxis()->SetRangeUser(0,(max2>max1?max2:max1)*1.45);
        hpad->DrawCopy("c");
        h2->DrawCopy("p same");
        h1->DrawCopy("p same");
        // gPad->SetLogy();

        TLegend* leg = new TLegend(0.75,0.75,0.9,0.92);
        leg->AddEntry(h1,legname1.Data(),"pe");
        leg->AddEntry(h2,legname2.Data(),"le");
        leg->Draw();

        drawLatex(0.4,0.86,particlename.Data(),0.055);
        drawLatex(0.2,0.86,centname[ic].Data(),0.055);
        drawLatex(0.2,0.8,Form("%0.2f<#eta<%0.2f",etal,etah),0.055);
        drawLatex(0.2,0.74,Form("%0.2f<p_{T}<%0.2f GeV/c",ptl,pth),0.055);

        ipad++;
        if (ipad==changepage) {c->cd();addpdf(pdf);c->Clear();ipad=0;c->Divide(a,b);} 
        // delete h1;delete h2;
        delete hpad;
      }
    }  
    // delete h3dca1;delete h3dca2;
  }
  c->cd();
  // drawLatex("");
  addpdf(pdf);
  pdf->On();
  pdf->Close();
  pdf = new TPDF(Form("plots_%s_nhits.pdf",outname.Data()));
  pdf->Off();

  ipad=0;
  gPad->SetLogy(0);
  //projection nhits
  for (int ic=0;ic<ncent;ic++) 
  {
    // hn1->GetAxis(4)->SetRangeUser(ic,ic);
    hn1->GetAxis(2)->SetRangeUser(0,3);
    hn2->GetAxis(2)->SetRangeUser(0,3);
    hn1->GetAxis(3)->SetRangeUser(0,100);
    hn2->GetAxis(3)->SetRangeUser(0,100);

    hn1->GetAxis(4)->SetRangeUser(centL[ic],centH[ic]);
    hn2->GetAxis(4)->SetRangeUser(centL[ic],centH[ic]);
    hn1->GetAxis(5)->SetRangeUser(charge1,charge1);
    hn2->GetAxis(5)->SetRangeUser(charge2,charge2);

    TH3F* h3hits1 = (TH3F*)hn1->Projection(0,1,3,"E");
    TH3F* h3hits2 = (TH3F*)hn2->Projection(0,1,3,"E");

    for (int ie=0;ie<20;ie++) 
    {
      for (int ip=0;ip<15;ip++)
      { 
        c->cd(ipad+1);
        double pth=0.1*(ip+1)+0;
        double ptl=0.1*(ip)+0;
        double etah=0.1*ie-2;
        double etal=0.1*(ie+1)-2;
        if (etah>0. || etah<=-2. || ptl>=2.5 || ptl<0.15) continue;
        // hn1->GetAxis(0)->SetRangeUser(ptl,pth);
        // hn1->GetAxis(1)->SetRangeUser(etal,etah);
        // hn1->GetAxis(4)->SetRangeUser(ic,ic);
        // hn1->GetAxis(5)->SetRangeUser(charge1,charge1);
        // TH1F* h1 = (TH1F*)hn1->Projection(3,"E");
        //
        // hn2->GetAxis(0)->SetRangeUser(ptl,pth);
        // hn2->GetAxis(1)->SetRangeUser(etal,etah);
        // hn2->GetAxis(4)->SetRangeUser(ic,ic);
        // hn2->GetAxis(5)->SetRangeUser(charge2,charge2);
        // TH1F* h2 = (TH1F*)hn2->Projection(3,"E");

        TH1D* h1 = (TH1D*)h3hits1->ProjectionZ("h1",h3hits1->GetXaxis()->FindBin(ptl+1e-6),h3hits1->GetXaxis()->FindBin(pth-1e-6),
            h3hits1->GetYaxis()->FindBin(etal+1e-6),h3hits1->GetYaxis()->FindBin(etah-1e-6));
        TH1D* h2 = (TH1D*)h3hits2->ProjectionZ("h2",h3hits2->GetXaxis()->FindBin(ptl+1e-6),h3hits2->GetXaxis()->FindBin(pth-1e-6),
            h3hits2->GetYaxis()->FindBin(etal+1e-6),h3hits2->GetYaxis()->FindBin(etah-1e-6));


        setHistStyle(h1, kRed-4, kOpenCircle, 0.8, 0); 
        setHistStyle(h2, kBlue, kFullCircle, 0.8, 0); 
        h1->GetXaxis()->SetTitle("nHitsFit");
        h1->GetYaxis()->SetTitle("Arb. Unit");
        Norm(h1);Norm(h2);
        TH1D* hpad = new TH1D("h","h;nHitsFit;Arb. Unit",1,0,120);
        // hpad->GetYaxis()->SetRangeUser(0,0.12);
        double max2 = h2->GetMaximum();
        double max1 = h1->GetMaximum();

        hpad->GetYaxis()->SetRangeUser(0,(max2>max1?max2:max1)*1.45);
        hpad->DrawCopy("c");
        h2->DrawCopy("p same");
        h1->DrawCopy("p same");

        TLegend* leg = new TLegend(0.75,0.75,0.9,0.92);
        leg->AddEntry(h1,legname1.Data(),"pe");
        leg->AddEntry(h2,legname2.Data(),"pe");
        leg->Draw();

        drawLatex(0.4,0.86,particlename.Data(),0.055);
        drawLatex(0.2,0.86,centname[ic].Data(),0.055);
        drawLatex(0.2,0.8,Form("%0.2f<#eta<%0.2f",etal,etah),0.055);
        drawLatex(0.2,0.74,Form("%0.2f<p_{T}<%0.2f GeV/c",ptl,pth),0.055);

        ipad++;
        if (ipad==changepage) {c->cd();addpdf(pdf);c->Clear();c->Divide(a,b);ipad=0;} 
        delete hpad;
      }
    }  
  }

  c->cd();
  // drawLatex("");
  addpdf(pdf);
  pdf->On();
  pdf->Close();
}
void calEfficiency(TString filename, TString particle, TString outname)
{
  SetsPhenixStyle();
  TFile* f = new TFile(filename.Data());
  TH3F* h3rc = (TH3F*)f->Get("hrc");
  TH3F* h3mc = (TH3F*)f->Get("hmc");

  int const ncent = 5;
  int centH[ncent]={8,1,3,5,8};
  int centL[ncent]={0,0,2,4,6};
  TString centname[ncent]={"0-80%","60-80%","40-60%","20-40%","0-20%"};
  int color[ncent]={kBlack,kRed,kOrange,kGreen,kBlue};

  TCanvas* c = new TCanvas("c","c",800,600);
  TPDF* pdf = new TPDF(Form("plots_%s.pdf",outname.Data()));
  pdf->Off();

  int a=3;int b=2;
  c->Divide(a,b);
  int changepage = a*b;
  int ipad=0;

  for (int ieta=0;ieta<24;ieta++)
  {
    c->cd(ipad+1);
    double etal=0.1*ieta-2.4;
    double etah=0.1*(ieta+1)-2.4;
    if (etal>=0.) continue; 

    TH1F* hpad = new TH1F("h","h;p_{T} (GeV/c);Efficiency",1,0,0.5);
    hpad->GetYaxis()->SetRangeUser(0,1.2); 
    hpad->DrawCopy();
    gPad->SetGridy();

    TLegend* leg;
    if (etal>=-2.2 && etah<=-0.100001) leg = new TLegend(0.65,0.2,0.9,0.5); 
    else leg = new TLegend(0.65,0.57,0.9,0.83);
    TH1F* heff[ncent];
    for (int ic=0;ic<ncent;ic++){
      TH1F* h1rc = (TH1F*)h3rc->ProjectionX("h3rc", h3rc->GetYaxis()->FindBin(etal+1e-6),  h3rc->GetYaxis()->FindBin(etah-1e-6), centL[ic]+1, centH[ic]+1);
      TH1F* h1mc = (TH1F*)h3mc->ProjectionX("h3mc", h3mc->GetYaxis()->FindBin(etal+1e-6),  h3mc->GetYaxis()->FindBin(etah-1e-6), centL[ic]+1, centH[ic]+1);

      heff[ic] = (TH1F*)h1rc->Clone("heff");
      heff[ic]->Divide(h1mc);
      setHistStyle(heff[ic],color[ic], kOpenCircle, 0.8, 0);
      heff[ic]->DrawCopy("same"); 
      leg->AddEntry(heff[ic],centname[ic].Data(),"pe");
    }
    leg->Draw();
    drawLatex(0.65,0.86,particle.Data(),0.055);
    // drawLatex(0.2,0.88,centname[ic].Data(),0.055);
    drawLatex(0.2,0.86,Form("%0.2f<#eta<%0.2f",etal,etah),0.055);
    // drawLatex(0.2,0.8,Form("%0.2f<p_{T}<%0.2f GeV/c",ptl,pth),0.055);

    ipad++;
    if (ipad==changepage) {c->cd();addpdf(pdf);c->Clear();ipad=0;c->Divide(a,b);} 
    // delete h1;delete h2;
    delete hpad;

  }

  c->cd();
  addpdf(pdf);
  drawLatex(0.2,0.8,"nHitsdEdx>10 && gDca <3 && ",0.055);
  // drawLatex(0.2,0.6,"nHitsFit>20 && nHitsFit/nHitsPoss>0.52",0.055);
  drawLatex(0.2,0.6,"nHitsFit>20 ",0.055);
  pdf->On();
  pdf->Close();

}
void drawQa()
{
  // // drawplots("data_test.root","hpion",1,"Data","out_pip.root","hpion",1,"MC","#pi^{+}","3.2GeVpip");
  // drawplots("data_test.root","hpion",-1,"Data","out_pim.root","hpion",-1,"MC","#pi^{-}","3.2GeVpim");
  // // drawplots("data.root","hproton",1,"Data","out_p.root","hproton",1,"MC","p","3.2GeVp");
  // // drawplots("data.root","hproton",-1,"Data","out_ap.root","hproton",-1,"MC","#bar{p}","3.2GeVap");
  // drawplots("data_test.root","hkaon",1,"Data","out_kp.root","hkaon",1,"MC","K^{+}","3.2GeVkp");
  // drawplots("data_test.root","hkaon",-1,"Data","out_km.root","hkaon",-1,"MC","K^{-}","3.2GeVkm");
  // //
  // drawplots("data.root","hpion",1,"Data","out_pip.root","hpion",1,"MC","#pi^{+}","3.2GeVpip");
  // drawplots("data.root","hpion",-1,"Data","out_pim.root","hpion",-1,"MC","#pi^{-}","3.2GeVpim");
  // drawplots("data.root","hproton",1,"Data","out_p.root","hproton",1,"MC","p","3.2GeVp");
  // drawplots("data.root","hproton",-1,"Data","out_ap.root","hproton",-1,"MC","#bar{p}","3.2GeVap");
  // drawplots("data.root","hkaon",1,"Data","out_kp.root","hkaon",1,"MC","K^{+}","3.2GeVkp");
  // drawplots("data.root","hkaon",-1,"Data","out_km.root","hkaon",-1,"MC","K^{-}","3.2GeVkm");
  //
  // calEfficiency("out_pip.root", "#pi^{+}", "eff_pip_ndedx");
  // calEfficiency("out_pim.root", "#pi^{-}", "eff_pim_ndedx");
  // calEfficiency("out_kp.root", "K^{+}", "eff_kp_ndedx");
  // calEfficiency("out_km.root", "K^{-}", "eff_km_ndedx");
  // calEfficiency("out_p.root", "p", "eff_p_ndedx");
  // calEfficiency("out_ap.root", "#bar{p}", "eff_ap_ndedx");
  calEfficiency("out_pip.root", "#pi^{+}", "eff_pip");
  calEfficiency("out_pim.root", "#pi^{-}", "eff_pim");
  calEfficiency("out_kp.root", "K^{+}", "eff_kp");
  calEfficiency("out_km.root", "K^{-}", "eff_km");
  calEfficiency("out_p.root", "p", "eff_p");
  calEfficiency("out_ap.root", "#bar{p}", "eff_ap");


}
