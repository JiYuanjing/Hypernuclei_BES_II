#include "sPhenixStyle.h"
#include "blaskwave.h"
// #include "style.h"

int npages=0;
bool savepng=false;
int ntmpfitplots=0;
int nsignalpdfpages=0;
bool savesignalpdf=false;

void graphscale(TGraphErrors* g, double scale)
{
  double x, y, err;
  for (int i=0;i<g->GetN();i++)
  {
    g->GetPoint( i, x, y);
    err = g->GetErrorY(i);
    g->SetPoint( i, x, y*scale);
    g->SetPointError(i, 0, err*scale);
  }
}
void graphscale(TGraph* g, double scale)
{
  double x, y, err;
  for (int i=0;i<g->GetN();i++)
  {
    g->GetPoint( i, x, y);
    // err = g->GetErrorY(i);
    g->SetPoint( i, x, y*scale);
    // g->SetPointError(i, 0, err*scale);
  }
}
TGraphErrors* combinechannels(TGraphErrors* g2b, TGraphErrors* g3b, TString gname)
{
  double x3,x2, y2,y3, err2,err3;
  double x, y,err;
  TGraphErrors* gcomb = (TGraphErrors*)g3b->Clone(gname.Data());
  for (int i=0;i<g3b->GetN()-1;i++)
  {
      g2b->GetPoint( i, x2, y2);
      g3b->GetPoint( i, x3, y3);
      cout << "y3: "<<x3<< " y2:"<<x2<< endl;
      err2 = g2b->GetErrorY(i);
      err3 = g3b->GetErrorY(i);
      x = x3;
      y = (y2/err2/err2+y3/err3/err3)/(1./err3/err3+1./err2/err2);
      err = sqrt(1./(1./err3/err3+1./err2/err2));
      gcomb->SetPoint( i, x, y);
      gcomb->SetPointError(i, 0, err);
  }
  return gcomb;
}
TGraphErrors* combinechannels_sys(TGraphErrors* g2b, TGraphErrors* g3b, TGraphErrors* gstat, TString gname)
{
  double x3,x2, y2,y3, err2,err3;
  double x, y,err;
  TGraphErrors* gcomb = (TGraphErrors*)g3b->Clone(gname.Data());
  for (int i=0;i<g3b->GetN()-1;i++)
  {
      g2b->GetPoint( i, x2, y2);
      g3b->GetPoint( i, x3, y3);
      gstat->GetPoint( i, x, y);
      cout << "y3: "<<x3<< " y2:"<<x2<< endl;
      err2 = g2b->GetErrorY(i);
      err3 = g3b->GetErrorY(i);
      // y = (y2/err2/err2+y3/err3/err3)/(1./err3/err3+1./err2/err2);
      err = sqrt(1./(1./err3/err3+1./err2/err2));
      gcomb->SetPoint( i, x, y);
      gcomb->SetPointError(i, 0, err);
  }
  return gcomb;
}

bool Norm(TH1F* h)
{
  double scale = 1./h->Integral();
  h->Scale(fabs(scale/h->GetBinWidth(1)));
  // h->Scale(fabs(1./h->GetMaximum()));
  h->SetDirectory(0);
  h->GetYaxis()->SetTitle("Arb. Unit");
  if (scale>0) return true;
  else return false;
}
double BrErr(double y2b, double y2berr, double y3b, double y3berr)
{
   return sqrt( pow((y2b/(y3b+y2b)/(y3b+y2b))*y3berr, 2) +  pow((y3b/(y3b+y2b)/(y3b+y2b))*y2berr, 2));
}
template <typename Hist>
void setHistStyle(Hist* h, int color, int markerstyle, double size , int mode = 0)
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
TH1F* reBinHist(double xmid,double xlow, double xhigh, TH1F* h, int n)
{
  int const nMax = 2000;
  double xedge[nMax];
  int nbin=0;
  xedge[0]=xlow;
  double rebin = n;
  while (xedge[nbin]<=xmid)
  {
    nbin++;
    xedge[nbin]=xedge[nbin-1]+h->GetBinWidth(1)*rebin;
  }
  int turningpoint = nbin+1;
  rebin = n*4;
  while (xedge[nbin]<xhigh)
  {
    nbin++;
    xedge[nbin]=xedge[nbin-1]+h->GetBinWidth(1)*rebin;
    // h->SetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*1.5),h->GetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*1.5))*0.5);
    // h->SetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*0.5),h->GetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*0.5))*0.5);
  }
  // return (TH1F*)h->Rebin(nbin,h->GetName(),xedge);
  h = (TH1F*)h->Rebin(nbin,h->GetName(),xedge);
  for (int ibin=0;ibin<turningpoint;ibin++)
  {
    h->SetBinContent(ibin,h->GetBinContent(ibin)/(1.0*n));
    h->SetBinError(ibin,h->GetBinError(ibin)/sqrt(1.0*n));
  }
  for (int ibin=turningpoint;ibin<h->GetNbinsX();ibin++)
  {
    h->SetBinContent(ibin,h->GetBinContent(ibin)/(1.0*rebin));
    h->SetBinError(ibin,h->GetBinError(ibin)/sqrt(1.0*rebin));
  }
  return h;
}
Double_t linebk(Double_t* x, Double_t* par)
{
  if (x[0]>2.987 && x[0]<2.997) TF1::RejectPoint();
  return par[0]+x[0]*par[1];
}

void drawCompareDataEmbedding(TString name, double scale, TFile* fH3L, TFile* fH3Lbk, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="masshistname")
{
  cout <<fH3L->GetName()<<" "<<xTitle.Data() << endl;
  int  const bins=1;
  double ptedge[bins+1]={0,4};
  double yedge[bins+1]={-1.3,0.1};

  fH3L->cd();
  TH3F* h2H3L = (TH3F*)fH3L->Get((name).Data())->Clone("h3sig");
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  TH3F* h2H3LM = (TH3F*)fH3L->Get((text).Data())->Clone("h3sigm");
  h2H3LM->Sumw2();
  h2H3LM->SetDirectory(0);
  fH3Lbk->cd();
  TH3F* h2H3Lbk = (TH3F*)fH3Lbk->Get((name).Data())->Clone("h3bk");
  h2H3Lbk->Sumw2();
  h2H3Lbk->SetDirectory(0);
  TH3F* h2H3LbkM = (TH3F*)fH3Lbk->Get((text).Data())->Clone("h3bkm");
  h2H3LbkM->Sumw2();
  h2H3LbkM->SetDirectory(0);
  // double sc = scale;
  // do the scaling with the same cuts
  double lowpt=ptedge[0], highpt=ptedge[bins], lowy=yedge[0], highy=yedge[bins];
  TH1F* hsigm = (TH1F*)h2H3LM->ProjectionY("hsigm", h2H3LM->GetXaxis()->FindBin(lowpt+1e-6),  h2H3LM->GetXaxis()->FindBin(highpt-1e-6), h2H3LM->GetZaxis()->FindBin(lowy), h2H3LM->GetZaxis()->FindBin(highy));
  hsigm->SetDirectory(0);
  TH1F* hbkm = (TH1F*)h2H3LbkM->ProjectionY("hbkm", h2H3LbkM->GetXaxis()->FindBin(lowpt+1e-6),  h2H3LbkM->GetXaxis()->FindBin(highpt-1e-6), h2H3LbkM->GetZaxis()->FindBin(lowy), h2H3LbkM->GetZaxis()->FindBin(highy));
  hbkm->SetDirectory(0);
  //scale
  double ml=3.005,mh=3.02;
  double sig_sb =  hsigm->Integral(hsigm->GetXaxis()->FindBin(ml),  hsigm->GetXaxis()->FindBin(mh));
  double bk_sb =   hbkm->Integral(hbkm->GetXaxis()->FindBin(ml),  hbkm->GetXaxis()->FindBin(mh));
  double sc = sig_sb/bk_sb;
  h2H3Lbk->Scale(sc);

  h2H3L->Add(h2H3Lbk, -1);
  cout <<"scale: "<<sc<< endl;
  // check whether the scale is correct
  // setHistStyle(hbkm,kRed, kOpenCircle, 1.2, 0);
  // setHistStyle(hsigm,kBlack, kFullCircle, 1.2, 0);
  // hbkm->Scale(sc);
  // hsigm->Draw();
  // hbkm->Draw("p same");
  // addpdf(pdf);

  fLa->cd();
  TH3F* h2La = (TH3F*)fLa->Get((name).Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  h2La->Scale(sc);
  c->Clear();

  c->Divide(1,1);

  cout <<"start compare... " <<endl;
  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);

    TH1F* h1H3= (TH1F*)h2H3L->ProjectionY(Form("h1H3%d",i), h2H3L->GetXaxis()->FindBin(ptedge[i]),h2H3L->GetXaxis()->FindBin(ptedge[i+1]), h2H3L->GetZaxis()->FindBin(yedge[i]),h2H3L->GetZaxis()->FindBin(yedge[i+1]));
    TH1F* h1La= (TH1F*)h2La->ProjectionY(Form("h1La%d",i), h2La->GetXaxis()->FindBin(ptedge[i]),h2La->GetXaxis()->FindBin(ptedge[i+1]), h2La->GetZaxis()->FindBin(yedge[i]),h2La->GetZaxis()->FindBin(yedge[i+1]));

    h1La->SetLineColor(kRed);
    h1La->SetMarkerColor(kRed);
    h1La->SetMarkerStyle(kOpenCircle);
    h1La->SetMarkerSize(1);
    h1H3->SetLineColor(kBlue);
    h1H3->SetMarkerColor(kBlue);
    h1H3->SetMarkerSize(1);
    h1H3->Rebin(4);
    if (rebin>1) {
      h1H3->Rebin(rebin);
      h1La->Rebin(rebin);
    }
    if (!Norm(h1La)) cout<<"warning:"<< fLa->GetName()<< " "<<xTitle.Data()<<" is nagetive!" <<endl;
    if (!Norm(h1H3)) cout<<"warning:"<< fH3L->GetName()<< " "<<xTitle.Data()<<" is nagetive!" <<endl;
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.4*h1H3->GetMaximum());
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    drawLatex( 0.18, 0.88, Form("%0.1f<p_{T}<%0.1f GeV/c", ptedge[i],ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.75,0.75,0.9,0.9);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
    // p12->cd();
    // TH1F* hratio = (TH1F*)h1H3->Clone("hratio");
    // hratio->Divide(h1La);
    // hratio->Draw();
    // hratio->GetYaxis()->SetRangeUser(hratio->GetMinimum()*0.5,hratio->GetMaximum()*2);
    // hratio->GetYaxis()->SetNdivisions(206);
    // hratio->GetYaxis()->SetTitleOffset(0.8);
    // hratio->GetYaxis()->SetTitleSize(0.08);
    // hratio->GetYaxis()->SetTitle(Form("%s/%s  ", legtitle1.Data(),legtitle2.Data()));
    // hratio->GetXaxis()->SetTitleSize(0.08);
    // hratio->GetXaxis()->SetTitleOffset(0.8);
    // drawLine( hratio->GetXaxis()->GetXmin(), 1., hratio->GetXaxis()->GetXmax(), 1, 1.5, 9, 1);
  }
  cout << "done!"<< endl;
  c->cd();
  addpdf(pdf);
}
void qaplots()
{
  SetsPhenixStyle();
  TFile* fsig = new TFile("fout_H3L_data_KF_0_80.root"); 
  TFile* fbk = new TFile("fout_H3L_data_RT_0_80.root"); 
  TFile* fmc = new TFile("fout_H3L_data_MC_RC_0080.root"); 
  TString histname = "hH3LMassPtY";
  double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0.;

  TCanvas* c = new TCanvas();
  TPDF* pdf = new TPDF("qa.pdf");
  pdf->Off();

  // TH3F* h2sig = (TH3F*)fsig->Get(histname.Data())->Clone("hptH3Lmass_sig");
  // h2sig->SetDirectory(0);
  // TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  // hsig->SetDirectory(0);
  //
  // TH3F* h2bk = (TH3F*)fbk->Get(histname.Data())->Clone("hptH3Lmass_bk");
  // h2bk->SetDirectory(0);
  // TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  // hbk->SetDirectory(0);
  //
  // //scale
  // double ml=3.005,mh=3.02;
  // double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(ml),  hsig->GetXaxis()->FindBin(mh));
  // double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(ml),  hbk->GetXaxis()->FindBin(mh));
  // double scale = sig_sb/bk_sb;
  double scale = 1./4.7;
  //
  c->Draw();
  drawLatex(0.1,0.9,"H3L->{}^{3}He#pi",0.055);
  drawLatex(0.1,0.82,"0<p_{T}<4 GeV/c",0.055);
  drawLatex(0.1,0.74,"-1.3<y<0.1",0.055);
  drawLatex(0.1,0.66,"0-80%",0.055);
  drawLatex(0.1,0.58,"Run 2020 FXT Au+Au 3.9 GeV",0.055);
  drawLatex(0.1,0.5,"Data samples selected within",0.055);
  drawLatex(0.1,0.42,"2.988<M(H3L)<2.996 GeV/c",0.055);
  addpdf(pdf);
  c->Clear();

  drawCompareDataEmbedding("h3H3L_lSig", scale, fsig, fbk, fmc, c, pdf, "l (cm)" , "p", 2,  "H3L Data", "MC", "h3H3L_lMass");
  drawCompareDataEmbedding("h3H3L_ldlSig", scale, fsig, fbk, fmc, c, pdf, "l/dl" , "p", 1,  "H3L Data", "MC", "h3H3L_ldlMass");
  drawCompareDataEmbedding("h3H3L_chi2ndfSig", scale, fsig, fbk, fmc, c, pdf, "#chi^{2}_{NDF}" , "p", 1,  "H3L Data", "MC", "h3H3L_chi2ndfMass");
  drawCompareDataEmbedding("h3H3L_chi2topoSig", scale, fsig, fbk, fmc, c, pdf, "#chi^{2}_{Topo}" , "p", 1,  "H3L Data", "MC", "h3H3L_chi2topoMass");
  drawCompareDataEmbedding("hHenHitsPtY", scale, fsig, fbk, fmc, c, pdf, "{}^{3}He nHitsFit" , "p", 1,  "H3L Data", "MC", "hHenHitsMass");
  drawCompareDataEmbedding("hPinHitsPtY", scale, fsig, fbk, fmc, c, pdf, "#pi nHitsFit" , "p", 1,  "H3L Data", "MC", "hPinHitsMass");
  drawCompareDataEmbedding("hHechi2primPtY", scale, fsig, fbk, fmc, c, pdf, "{}^{3}He #chi^{2}_{prim}" , "p", 1,  "H3L Data", "MC", "hHechi2primMass");
  drawCompareDataEmbedding("hPichi2primPtY", scale, fsig, fbk, fmc, c, pdf, "#pi #chi^{2}_{prim}" , "p", 1,  "H3L Data", "MC", "hPichi2primMass");

  pdf->On();
  pdf->Close();
}

void drawData()
{
  // YieldSysCuts(1); 
  // YieldSysCuts(0); 
  qaplots();
}

