#include "sPhenixStyle.h"
#include "blaskwave.h"
// #include "style.h"

int npages=0;
bool savepng=false;
int ntmpfitplots=0;
int nsignalpdfpages=0;
bool savesignalpdf=false;


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
double fityield2(  TString histname, double lowpt, double highpt, double lowy, double highy, double& err, TFile* f1, TFile* f2, TCanvas* c, TPDF* pdf, int centL=8, int centH=9,int mode=1)
{ 
  c->cd();
  int centname[10]={80,70,60,50,40,30,20,10,5,0};
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  // double nEvents_se = hcent_se->Integral(4,9);
  double nEvents_se = hcent_se->Integral(centL,centH);

  TH3F* h2sig = (TH3F*)f1->Get(Form("%s", histname.Data()))->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy+1e-6), h2sig->GetZaxis()->FindBin(highy-1e-6));
  hsig->SetDirectory(0);

  TH3F* h2bk = (TH3F*)f2->Get(Form("%s", histname.Data()))->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy+1e-6), h2bk->GetZaxis()->FindBin(highy-1e-6));
  hbk->SetDirectory(0);

  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.02),  hsig->GetXaxis()->FindBin(3.04));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.02),  hbk->GetXaxis()->FindBin(3.04));

  double scale = sig_sb/bk_sb;
  hbk->Scale(scale);
  hsig->SetMarkerColor(kRed);
  hsig->SetLineColor(kRed);
  hsig->GetXaxis()->SetRangeUser(2.97,3.04);
  // hsig->GetXaxis()->SetRangeUser(2.97,3.02);
  hsig->Draw();
  hsig->GetXaxis()->SetTitle("Mass({}^{3}He#pi) (GeV/c^{2})");
  hsig->GetYaxis()->SetTitle(Form("Counts per %0.1f MeV",hsig->GetBinWidth(1)*1000));
  hbk->Draw("same");
  TLegend* leg = new TLegend( 0.2, 0.78 ,0.33,0.9 );
  leg->AddEntry(hbk, "ME","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->Draw();

  drawLatex( 0.18,0.7,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.18,0.61,Form("%0.2f<y<%0.2f",lowy, highy ), 0.055);
  drawLatex( 0.18,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  drawLatex( 0.18,0.47,Form("0-80%s", "%"), 0.055);
  addpdf(pdf); 

  if (savepng) gPad->SaveAs(Form("plots/png/%d.png", npages++));
  if (savesignalpdf) gPad->SaveAs(Form("plots/signalpdf/%d.pdf", nsignalpdfpages++));

  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  hsig_bk->Rebin();
  hsig_bk->Draw();
  hsig_bk->GetYaxis()->SetTitle(Form("Counts per %0.1f MeV",hsig_bk->GetBinWidth(1)*1000));

  TF1* fit= new TF1("fit" ,"gausn(0)+pol1(3)", 2.97,3.05 );
  TF1* fitgaus = new TF1("fitgaus","gausn(0)", 2.97,3.05 );
  // TF1* resfit = new TF1("resfit","pol1", 2.95,3.05 );
  TF1* resfit = new TF1("resfit",linebk, 2.95,3.05,2 );
  
  hsig_bk->GetXaxis()->SetRangeUser(2.97,3.04);
  hsig_bk->Fit(resfit,"RN0");
  fit->SetLineColor(kRed);
  double yield_bc = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.98), hsig_bk->GetXaxis()->FindBin(3.));
  // double para[5]={yield_bc*hsig_bk->GetBinWidth(1)/sqrt(2*3.1415), 2.991, 0.0015,  resfit->GetParameter(0), resfit->GetParameter(1)};
  double para[7]={yield_bc*hsig_bk->GetBinWidth(1), 2.992, 0.0014,  resfit->GetParameter(0), resfit->GetParameter(1),0,0};
  fit->SetParameters(para);
  double lowx=2.97 ,highx =3.02;
  hsig_bk->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_bk->Draw("same");
  hsig_bk->Fit(fit,"RN0");
  fit->GetParameters(para);

  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );

  fit->SetParameters(para);

  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");

  setHistStyle(fit, kMagenta, 5, 2.5 ,1);
  fit->GetParameters(para);
  fit->Draw("same");
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;
  hsig_bk->GetXaxis()->SetRangeUser( 2.97, 3.02);
  fitgaus->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));

  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  double yield_me = (fit->GetParameter(0))/hsig_bk->GetBinWidth(1);
  // double yield_counts = hsig_bk->IntegralAndError(hsig_bk->GetXaxis()->FindBin(2.987), hsig_bk->GetXaxis()->FindBin(2.998), err, "width")-resfit->Integral( 2.987, 2.998)/hsig_bk->GetBinWidth(1);
  double yield_counts = hsig_bk->IntegralAndError(hsig_bk->GetXaxis()->FindBin(2.985), hsig_bk->GetXaxis()->FindBin(3.), err, "")-resfit->Integral( 2.985, 3.)/hsig_bk->GetBinWidth(1);
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double significance = yield_counts/sqrt(sp_counts);
  // double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  // err = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.987), hsig_bk->GetXaxis()->FindBin(2.998))-resfit->Integral( 2.987, 2.998)/hsig_bk->GetBinWidth(1);
 
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_counts), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.1f", significance), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.4f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.2f<y<%0.2f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  drawLatex( 0.62,0.47,Form("%d-%d%s", centname[centH],centname[centL-1],"%"), 0.055);
  //
  TLegend* leg2 = new TLegend(0.6,0.9,0.9,0.7);
  leg2->AddEntry(hsig_bk,"SE-ME","pe");
  leg2->AddEntry(fit,"total","l");
  leg2->Draw();

  addpdf(pdf);
  if (savepng) gPad->SaveAs(Form("plots/png/%d.png", npages++));
  if (savesignalpdf) gPad->SaveAs(Form("plots/signalpdf/%d.pdf", nsignalpdfpages++));

  return (significance>2)?yield_counts:0;
}
void drawMixDataScanSys_yield(int icut, TString histname, TString topohistname, TString pdfname, double topocut, double& Br, double & error, int pfitmode,int centL=1, int centH=9)
{
  double highpt = 10, lowpt = 0, lowy=-1.5, highy = 0.;
  cout << histname.Data()<<endl;
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  TPDF* pdf = new TPDF(pdfname);
  pdf->Off();
  gStyle->SetPalette(1);

  TFile *f1 = TFile::Open("fout_H3L_data_KF_0_80.root"); 
  TH3F* h2sig = (TH3F*)f1->Get(histname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  // TH3F* h2sig->Project3D("xz");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral( centL,centH);

  TFile *f2 = TFile::Open("fout_H3L_data_RT_0_80.root"); 
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  //scale
  double ml=3.005,mh=3.02;
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(ml),  hsig->GetXaxis()->FindBin(mh));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(ml),  hbk->GetXaxis()->FindBin(mh));
  double scale = sig_sb/bk_sb;
  cout<<"ME scale: " <<1./scale << endl;
  hbk->Scale(scale);
  hsig->Rebin();
  hbk->Rebin();

  hsig->Draw();
  hsig->GetXaxis()->SetTitle("Mass({}^{3}He#pi) (GeV/c^{2})");
  hsig->GetYaxis()->SetTitle("Counts");
  hsig->GetYaxis()->SetRangeUser(-0.1*hsig->GetMaximum(), hsig->GetMaximum()*1.1);
  
  setHistStyle(hbk, kRed, kOpenCircle, 0.5);
  hbk->Draw("same");
  
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "RT", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();

  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-80%s", "%"), 0.055);
  // drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  c->cd();
  // calculate the significance
  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.2);
  hsig_bk->Rebin();
  hsig_bk->Draw();

  // TF1* fit = new TF1("fit" ,"gaus(0)+pol1(3)", 2.97,3.02 );
  TF1* fit = new TF1("fit" ,"gausn(0)+pol1(3)", 2.97,3.02 );
  TF1* resfit = new TF1("resfit" ,"pol1", 2.95,3.05 );
  hsig_bk->GetXaxis()->SetRangeUser(2.97,2.985);
  hsig_bk->Fit(resfit,"R");
  fit->SetLineColor(kRed);
  double yield_bc = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.98), hsig_bk->GetXaxis()->FindBin(3.));
  // double para[5]={yield_bc*hsig_bk->GetBinWidth(1)/sqrt(2*3.1415), 2.991, 0.0015,  resfit->GetParameter(0), resfit->GetParameter(1)};
  double para[5]={yield_bc*hsig_bk->GetBinWidth(1), 2.991, 0.0014,  resfit->GetParameter(0), resfit->GetParameter(1)};
  fit->SetParameters(para);
  double lowx=2.97 ,highx =3.02;
  hsig_bk->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_bk->Draw("same");
  hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );


  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  // double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1)*fit->GetParameter(2)*sqrt(2*3.1415);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-3*sigma), hsig_bk->GetXaxis()->FindBin(mean+3*sigma))-resfit->Integral(mean-3*sigma,mean+3*sigma)/hsig_bk->GetBinWidth(1);
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-3*sigma), hbk->GetXaxis()->FindBin(mean+3*sigma));
  double significance = yield_counts/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->Draw();
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_me), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.0f", significance), 0.055);
  drawLatex( 0.2,0.47,Form("S/#DeltaS=%0.0f (ME)", s_me), 0.055);
  // drawLatex( 0.2,0.4,Form("S/#DeltaS=%0.0f (RT)", s_rt), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.4f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.62,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.62,0.47,Form("0-50%s", "%"), 0.055);
  // drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  // drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_bk, "SE-RT","pl");
  leg2->Draw();

  addpdf(pdf);

  //draw phase space

  // pdf->On();
  // pdf->Close();
  // return;

  TFile* fMc = TFile::Open("fout_H3L_data_MC_RC_0080.root");

  int const nybins = 2;

  double edge[nybins][9]={
    // { 0.6,1, 1.4, 1.8, 2.2, 3.2}, 
    { 0.5, 1.5, 2., 2.5, 3.5},
    // { 0.5,1.5, 2, 2.5, 3.5},
    { 1,1.5, 2, 3.5}
  };
  TH1F* hPhase[nybins];
  TH1F* hPhaseCor[nybins];

  // double ybin[nybins+1]={ 0, -0.2, -0.4, -0.6, -0.8};
  // double ybin[nybins+1]={ 0, -0.25, -0.5, -0.75, -1};
  // double ybin[nybins+1]={ 0, -0.5, -1};
  double ybin[nybins+1]={ 0, -0.5, -1};
  double ymid[nybins],yerr[nybins];
  for (int iy=0;iy<nybins;iy++)
  {
     ymid[iy] = 0.5*(ybin[iy]+ybin[iy+1]); 
     yerr[iy] = 0.5*fabs(ybin[iy]-ybin[iy+1]); 
  }

  int nptbins[nybins]={ 4, 4};
  for (int iy=0;iy<nybins;iy++){
    hPhase[iy] = new TH1F(Form("hPhase%d", iy), Form("hPhase%d;pt", iy), 300, 0.5, 3.5);
    hPhase[iy]=(TH1F*)hPhase[iy]->Rebin( nptbins[iy], Form("hPhase%d", iy), edge[iy]);
  }
  TFile* fMcH3L = new TFile("fout_H3L_data_MC_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  TFile* fRcH3L = new TFile("fout_H3L_data_MC_RC_0080.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get(histname.Data())->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(2);
  h2temp->RebinX(2);
  h2Eff->RebinY(2);
  h2Eff->RebinX(2);
  h2Eff->Divide(h2temp);
  // h2Eff->Draw("col text");
  h2Eff->Draw("col");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  h2Eff->GetYaxis()->SetTitle("Efficiency");
  addpdf(pdf);

  TH1F* heff[nybins];
  for (int ij=0;ij<nybins;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( nptbins[ij], Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin(nptbins[ij] , Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double dy=0.5;
  double dpt,pt;

  for (int ij=0; ij<nybins;ij++){
    double npt =nptbins[ij];
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield;
      gPad->SetLogy(0);
      dy = fabs(ybin[ij+1]-ybin[ij]);
      // if (icut!=14) yield = fityield2(histname, edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf, 1, 9); // normal mode
      yield = fityield2(histname, edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf, 1, 9); // normal mode
      // if (icut==14) yield = fityield2(histname, edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf, 1, 9); // normal mode
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );

      // cout << yield << " "<< yield/dy/dpt/pt/2./3.1415926/nEvents_se<<endl;
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    gPad->SetLogy();
    addpdf(pdf); 
    gPad->SetLogy(1);

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;

      double yield_cor = y3b/eff;
      double err = sqrt(y3berr*y3berr + efferr*efferr)*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);

      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << " eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");

    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

  }
  cout <<"finish "<< icut << " sys calculation" << endl;

  TFile * fout = new TFile(Form("outfile/fYield_0010_sys_%d.root", icut),"recreate");
  fout->cd();
  for (int i=0;i<nybins;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    heff[i]->Write();
  }
  h2Eff->Write();
  fout->Close();

  pdf->On();
  pdf->Close();
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
  // TFile* fsig = new TFile("fout_H3L_data_KF_0_10.root"); 
  // TFile* fbk = new TFile("fout_H3L_data_RT_0_10.root"); 
  // TFile* fmc = new TFile("fout_H3L_data_MC_RC_0_10.root"); 
  TString histname = "hH3LMassPtY";
  double highpt = 4, lowpt = 0, lowy=-1., highy = 0.;

  TCanvas* c = new TCanvas();
  TPDF* pdf = new TPDF("qa_0_10.pdf");
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
  // drawLatex(0.1,0.66,"0-80%",0.055);
  drawLatex(0.1,0.66,"0-10%",0.055);
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
void YieldSysCuts(int runscandata=1)
{
  // int runscandata=1;
  int const cent=0; // cent=1->0-10, =2->10-50, =0->0-50%
  TString centname[3]={"0080", "0010", "1040"};
  TString centname2[3]={"0-80", "0-10", "10-40"};
  SetsPhenixStyle();
  TCanvas* cc= new TCanvas();
  TPDF* pdf = new TPDF(Form("plots/finalYield_%s.pdf", centname[cent].Data()));
  pdf->On();

  int const nsys = 10;
  double p[nsys+1],perr[nsys+1];
  for (int i=0;i<nsys;i++)
  { 
    double ndfcut = 4;
    if (i==6) ndfcut=4;
    else if (i==7) ndfcut=4;
    if (runscandata) { 
     drawMixDataScanSys_yield(i, Form("hH3LMassPtYSysCut%d",i ), Form("h3H3L_chi2ndfSysCut%d",i), Form("plots/syscuts_%s_ndfcut%d.pdf", centname[cent].Data(), i),  ndfcut, p[i], perr[i], 1); 
    }
  }
  if (runscandata) {
     npages=0;
     savepng=true;    
     savesignalpdf=true;    
     drawMixDataScanSys_yield(nsys, "hH3LMassPtY", "h3H3L_chi2ndf", Form("plots/syscuts_%s_ndfcut%d.pdf", centname[cent].Data(), nsys),  3.5, p[nsys], perr[nsys], 1); 
     savepng=false;
     savesignalpdf=false;
  }
  if (runscandata) return; 

  double cuts[20]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

  // drawMixDataScanSys(4, Form("hH3LMassPtYSysCut%d",4 ), Form("h3H3L_chi2topoSysCut%d",4), Form("syscuts_0050_topocut%d.pdf", i),  topocut, p[i], perr[i]); 
  TString cutsname[nsys]={"chi2topo<3","chi2topo<5","chi2ndf<3","chi2ndf<5","pi chi2prim>3","pi chi2prim>40"};
  TFile* f[nsys];
  TH1F* h[5][nsys+1]; 

  double br[nsys+1], brerr[nsys+1], y3brelerr[nsys+1];
  int color[20]={ 
                      kGray+2, kMagenta,  kSpring+4, kGreen+2, kBlue+1, kOrange+5,
                      kBlack, kCyan+2, kViolet-1, kPink+10, kTeal+2, kAzure+2,
                      kYellow-3, kMagenta+1 ,kRed
                    };
  // TLegend* leg = new TLegend(0.2,0.2,0.5,0.5);
  int const ny = 3;
  double yedge[ny+1]={0, -0.25,-0.5,-0.75};
  // double yedge[ny+1]={0, -0.2,-0.4, -0.6, -0.8};
  double syserr[ny][10], pt[ny][10], err[ny][10], yield[ny][10], toterr[ny][10];
  int nbins[ny];
  int nycol[5]={kRed, kGreen+2, kBlue, kMagenta, kGray+2};
  int marker[5]={kFullCircle, kFullSquare, kFullStar, 27};
  TGraphErrors* gYield[ny];
  TGraphErrors* gYieldSys[ny];
  TGraphErrors* gYieldToterr[ny];

  TCanvas* cyield= new TCanvas("cyield","cyield",800,700);
  cyield->cd();
  gPad->SetLogy();
  // TH1F* hpad = new TH1F("hpad", "hpad;p_{T} (GeV/c); B.r.#timesd^{2}N/2#pidp_{T}dy", 10, 1, 3.);
  TH1F* hpad = new TH1F("hpad", "hpad;p_{T} (GeV/c); d^{2}N/2#pip_{T}dp_{T}dy (c^{2}/GeV^{2})", 10, 1, 4);
  hpad->GetYaxis()->SetRangeUser(6e-10, 5e-2);
  hpad->GetYaxis()->SetNdivisions(504);
  hpad->GetXaxis()->SetNdivisions(406);
  hpad->Draw("XY");

  double scale[]={1, 0.1,0.01,0.001,0.0001};
  double _pmass = 2.99089;
  
  TF1 *bolt[ny];
	TF1 *ptbolt[ny];
	TF1 *bolthi[ny];
	TF1 *ptbolthi[ny];

	TF1 *mtexp[ny];
	TF1 *ptmtexp[ny];
	TF1 *mtexphi[ny];
	TF1 *mtexplo[ny];
	// //pt times pt corresponding to mt
	//
	TF1 *ptmptexp[ny];
	TF1 *ptptmptexp[ny];

  TF1 *pt3exp[ny];
  TF1 *ptpt3exp[ny];

  TF1 *pt2exp[5];
  TF1 *ptpt2exp[5];
  
  TF1* bw[4];
  TF1* ptbw[4];
	
  for (int i=0;i<ny;i++){
    bolt[i] = new TF1(Form("bolt[%d]",i), "[1]*sqrt([2]*[2]+x*x)*exp(-sqrt([2]*[2]+x*x)/[0])", 0,4);
    bolt[i]->FixParameter(2,_pmass);
    bolt[i]->SetParameters(0.2,5*pow(10, 4-i),3);
    ptbolt[i] = new TF1(Form("ptbolt[%d]",i), "(2*TMath::Pi())*[1]*x*sqrt([2]*[2]+x*x)*exp(-sqrt([2]*[2]+x*x)/[0])", 0,4);
    ptbolt[i]->SetParameters(bolt[i]->GetParameter(0),bolt[i]->GetParameter(1)*2*TMath::Pi(),bolt[i]->GetParameter(2));
    ptbolt[i]->FixParameter(2,_pmass);
    ptbolt[i]->SetLineColor(kBlue);
    ptbolt[i]->SetLineStyle(kBlue);
    bolt[i]->SetLineColor(kBlue);

    mtexp[i] = new TF1(Form("mtexp[%d]",i), "[1]*exp(-(sqrt([2]*[2]+x*x))/[0])", 0,4);
    ptmtexp[i] = new TF1(Form("ptmtexp[%d]",i), "2*TMath::Pi()*[1]*x*exp(-(sqrt([2]*[2]+x*x))/[0])", 0,4);
    mtexp[i]->SetParameters(0.2, 5e0, _pmass );
    ptmtexp[i]->SetParameters(0.2, 5e0, _pmass );
    mtexp[i]->FixParameter(2, _pmass);
    ptmtexp[i]->FixParameter(2, _pmass);
    mtexp[i]->SetLineColor(kMagenta+2);
    mtexp[i]->SetLineStyle(2);
    ptmtexp[i]->SetLineColor(kMagenta+2);

    // ptmptexp[i] = new TF1(Form("ptmptexp[%d]",i), "x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,4);
    // ptptmptexp[i] = new TF1(Form("ptptmptexp[%d]",i), "2*TMath::Pi()*x*x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,4);
    // ptmptexp[i]->SetParameters(mtexp[i]->GetParameter(0),mtexp[i]->GetParameter(1),_pmass);
    // ptptmptexp[i]->SetParameters(mtexp[i]->GetParameter(0),mtexp[i]->GetParameter(1),_pmass);
    // ptmptexp[i]->FixParameter(2, _pmass);
    // ptptmptexp[i]->FixParameter(2, _pmass);
    // ptmptexp[i]->SetLineColor(kGreen+2);
    // ptptmptexp[i]->SetLineColor(kGreen+2);

    pt2exp[i] = new TF1(Form("pt2exp%d",i), "[1]*exp(-x*x/[0])", 0,4);
    double pt2exppara[4][2]={{3.60820e+00, 2.00623e-04}, {2.51772e+00,4.61740e-05}, {2.32190e+00 ,3.79600e-06},
                             {2.31222e+00 ,2.20509e-07}};
    double pt2exppara2[4][2]={{2.44728e+00,7.81707e-05},{2.33004e+00,8.63919e-06},{1.70539e+00,2.27190e-06},{1.56005e+00,3.61407e-07}};

	  if (cent==1 || cent ==0)
    {
      pt2exp[i]->SetParameters(pt2exppara[i]);
    } 
    else if (cent==2)
    {
      pt2exp[i]->SetParameters(pt2exppara2[i]);
    }
    ptpt2exp[i] = new TF1(Form("ptpt2exp%d",i), "2*TMath::Pi()*[1]*x*exp(-x*x/[0])", 0,4);
	  ptpt2exp[i]->SetParameters(pt2exp[i]->GetParameter(0),pt2exp[i]->GetParameter(1));
    pt2exp[i]->SetLineColor(kGreen+2);
    pt2exp[i]->SetLineStyle(5);
    ptpt2exp[i]->SetLineColor(kGreen+2);

    pt3exp[i] = new TF1(Form("pt3exp%d",i), "[1]*exp(-pow(x,1.5)/[0])", 0,4);
    // pt3exp[i] = new TF1(Form("pt3exp%d",i), "[1]*exp(-pow(x,3)/[0])", 0,4);
    ptpt3exp[i] = new TF1(Form("ptpt3exp%d",i), "2*TMath::Pi()*[1]*x*exp(-pow(x,1.5)/[0])", 0,4);
    // ptpt3exp[i] = new TF1(Form("ptpt3exp%d",i), "2*TMath::Pi()*[1]*x*exp(-pow(x,3)/[0])", 0,4);
    pt3exp[i]->SetParameters(0.2,5e3*pow(10,-i));
    pt3exp[i]->SetLineColor(kOrange+7);
    pt3exp[i]->SetLineStyle(3);
    ptpt3exp[i]->SetLineColor(kOrange+7);

    if(cent==1 || cent==0){
      bw[0]=GetBGBWdNdpt(_pmass,-1.92877e+05 , 3.61565e-01 , 8.72736e+04, 1.63704e+01, "bw[0]");
      bw[1]=GetBGBWdNdpt(_pmass,-1.76682e+05, 2.68207e-01, 8.72736e+04, 1.80140e+02, "bw[1]");
      bw[2]=GetBGBWdNdpt(_pmass,-2.25520e+05, 2.51423e-01, 8.72736e+04, 3.96031e+01, "bw[2]");
      bw[3]=GetBGBWdNdpt(_pmass,-1.53363e+05, 2.56562e-01, 8.72736e+04, 1.61763e+00, "bw[3]");
    }
    else if (cent==2){
      bw[0]=GetBGBWdNdpt(_pmass,-1.64556e+05, 2.53835e-01, 8.72736e+04, 7.84797e+02, "bw[0]");
      bw[1]=GetBGBWdNdpt(_pmass,-1.90876e+05, 2.48533e-01, 8.72736e+04, 1.12982e+02, "bw[1]");
      bw[2]=GetBGBWdNdpt(_pmass,-2.09422e+05, 1.88065e-01, 8.72736e+04, 5.36434e+03, "bw[2]");
      bw[3]=GetBGBWdNdpt(_pmass,-1.75081e+05, 1.75796e-01, 8.72736e+04, 3.70492e+03, "bw[3]");
    }
    ptbw[i] = GetBGBWdNdptTimesPt(_pmass, bw[i]->GetParameter(1), bw[i]->GetParameter(2), bw[i]->GetParameter(3), bw[i]->GetParameter(4)*pow(10,i),Form("ptbw[%d]",i));
    ptbw[i]->SetRange(0,4);
    bw[i]->SetRange(0,4);
    ptbw[i]->FixParameter(0, _pmass);
    bw[i]->FixParameter(0, _pmass);
    // bw[i]->SetLineColor(kMagenta+2);
    // bw[i]->SetLineStyle(4);
    // ptbw[i]->SetLineColor(kMagenta+2);
	}

  double const R3=0.3;

  TLegend* leg = new TLegend(0.18,0.18,0.4,0.4);
  leg->SetTextSize(0.05);
  for (int ic=0; ic<ny; ic++) 
  {
    for (int i=0;i<nsys+1;i++)
    {
       f[i] = new TFile( Form("outfile/fYield_%s_sys_%d.root", centname[cent].Data(), i));
       h[ic][i] = (TH1F*)f[i]->Get(Form("hPhaseCor_%d", ic))->Clone(Form("htest_%d_%d", ic, i));
       h[ic][i]->SetDirectory(0);
       h[ic][i]->SetLineColor(color[i]);
       h[ic][i]->SetMarkerColor(color[i]);
       h[ic][i]->Scale(scale[ic]);
       f[i]->Close();
    }
  
    nbins[ic]=h[ic][nsys]->GetNbinsX();
    
    for (int ip=0;ip<h[ic][nsys]->GetNbinsX();ip++)
    {
       pt[ic][ip]=h[ic][nsys]->GetBinCenter(ip+1);
       err[ic][ip]=h[ic][nsys]->GetBinError(ip+1);
       yield[ic][ip]=h[ic][nsys]->GetBinContent(ip+1);
       syserr[ic][ip]=0;

       double sys_topo=0, sys_bincount=0, sys_track=0;

       for (int i=0;i<12;i++)
       {
          double tmp=h[ic][i]->GetBinContent(ip+1);
          syserr[ic][ip]+=pow(tmp-yield[ic][ip], 2)*0.5; // each sys has 2 cuts
          sys_topo+=pow(tmp-yield[ic][ip], 2)*0.5;
          // cout << "sys topo "<< i<<" "<< fabs(tmp-yield[ic][ip])/yield[ic][ip]<< endl;
       }
       {
          double tmp1=h[ic][12]->GetBinContent(ip+1);
          double tmp2=h[ic][13]->GetBinContent(ip+1);
          double delta1 = fabs(tmp1-yield[ic][ip]);
          double delta2 = fabs(tmp2-yield[ic][ip]);
          // double delta = delta1>delta2?delta1:delta2;
          double delta = delta2;
          syserr[ic][ip]+=delta*delta;
          sys_track=delta;
       }

       {
         double tmp=h[ic][14]->GetBinContent(ip+1);
         double delta_bincounting = fabs(tmp - yield[ic][ip]);
         syserr[ic][ip]+=delta_bincounting*delta_bincounting;
         sys_bincount = delta_bincounting;
       }
       cout <<"ybin=" <<ic<<" pt="<<pt[ic][ip] << " sys bin count: "<< sys_bincount/yield[ic][ip]<<" "
       << "sys_topo: "<<sqrt(sys_topo)/yield[ic][ip] <<" "
        << "sys_track: "<<sys_track/yield[ic][ip] <<endl;

       syserr[ic][ip] = sqrt(syserr[ic][ip]);
       toterr[ic][ip] = sqrt(syserr[ic][ip]*syserr[ic][ip]+err[ic][ip]*err[ic][ip]);
       cout<<yield[ic][ip]<<" "<< err[ic][ip]<<" rel: "<< syserr[ic][ip]/yield[ic][ip]<<endl;
    }

    gYieldSys[ic] = new TGraphErrors(nbins[ic], pt[ic], yield[ic], 0, syserr[ic]);
    gYieldSys[ic]->SetName(Form("gYieldsys_%s_%d", centname[cent].Data(), ic));
    gYieldToterr[ic] = new TGraphErrors(nbins[ic], pt[ic], yield[ic], 0, toterr[ic]);
    gYieldToterr[ic]->SetName(Form("gYieldtoterr_%s_%d", centname[cent].Data(), ic));
    
    gYield[ic] = new TGraphErrors(nbins[ic], pt[ic], yield[ic], 0, err[ic]);
    gYield[ic]->SetName(Form("gYield%d",ic));
   
    gYield[ic]->SetLineColor(nycol[ic]); 
    setHistStyle(gYield[ic],nycol[ic], kFullCircle, 1.5, 0);
    setHistStyle(gYieldSys[ic],nycol[ic], kFullCircle, 1.5, 0);
    // graphscale(gYield[ic], 1./((1-R3)*0.98*2./3.));
    // graphscale(gYieldSys[ic], 1./((1-R3)*0.98*2./3.));

    drawGraphWithSys(gYield[ic], gYieldSys[ic], nycol[ic], kFullCircle, 1.5, 1, 0,0.03);
    gYield[ic]->Fit(mtexp[ic],"N");
    mtexp[ic]->DrawCopy("same");

    if (ic==2) leg->AddEntry(gYield[ic], Form("y=(%0.1f,%0.2f)#times10^{%d}",yedge[ic],yedge[ic+1], -1*ic),"pe");
    if (ic==1) leg->AddEntry(gYield[ic], Form("y=(%0.2f,%0.1f)#times10^{%d}",yedge[ic],yedge[ic+1], -1*ic),"pe");
    if (ic==0) leg->AddEntry(gYield[ic], Form("y=(%0.0f,%0.2f)",yedge[ic],yedge[ic+1]),"pe");
  }
  leg->Draw();

  TGraphErrors* gempty = (TGraphErrors*)gYield[0]->Clone("gempty");
  gempty->SetMarkerColor(0);
  gempty->SetLineColor(0);

  TLegend* l2b = new TLegend( 0.65, 0.18, 0.85, 0.38);
  l2b->SetTextSize(0.05);
  l2b->SetNColumns(3);
  // l2b->AddEntry(g2b[0]," ","p");
  // // l2b->AddEntry();
  // l2b->AddEntry(g2b[1],"","p");
  // // l2b->AddEntry(g2b[2],"{}^{3}_{#Lambda}H#rightarrow ^{3}He#pi","pe");
  // l2b->AddEntry(gempty,"{}^{3}_{#Lambda}H#rightarrow ^{3}He#pi","p");
  // l2b->Draw();
  // TLegend* l3b = new TLegend( 0.2, 0.8, 0.5, 0.9);
  // l3b->SetTextSize(0.055);
  // l3b->SetNColumns(3);
  l2b->AddEntry(gYieldSys[0]," ","p");
  l2b->AddEntry(gYieldSys[1]," ","p");
  // l2b->AddEntry(g2b[1],"","pe");
  // l2b->AddEntry(g2b[2],"{}^{3}_{#Lambda}H#rightarrow{}^{3}He#pi","pe");
  l2b->AddEntry(gYieldSys[2],"{}^{3}_{#Lambda}H#rightarrowp#pid","p");
  l2b->Draw();
  TLegend* lf = new TLegend( 0.65, 0.75, 0.85, 0.85);
  lf->SetTextSize(0.05);
  lf->AddEntry(mtexp[0], "m_{T}-exp. fit", "l");
  lf->Draw("same");
  // drawSTAR(0.6,0.85,0.05);
  drawLatex(0.2,0.85,"Au+Au #sqrt{s_{NN}} = 3 GeV",0.05);
  // drawLatex(0.2,0.85,"Au+Au 3 GeV",0.05);
  drawLatex(0.25, 0.42, Form("%s%s", centname2[cent].Data(),"%"), 0.05);
  addpdf(pdf);
  gPad->SaveAs(Form("plots/spectra_%d.pdf", cent));
  cc->cd();
  gPad->SetLogy(0);
  for (int ic=0;ic<ny;ic++)
  {
    cout << "parameters:  "<<  ic <<" " <<mtexp[ic]->GetParameter(0)<<", "<< mtexp[ic]->GetParameter(1)/scale[ic]<< ", "<<  mtexp[ic]->GetParameter(2)<< endl;
  }
  // return;

  //fit dndy yield
  int const nmodel=5;
  double dndy[ny][nsys+1], dndym[ny], dndysys_model[ny],dndysys[ny],dndystat[ny]; //lastbin is default
  double dndymd[ny][nsys+1][5]; //lastbin is default
  double dndypt4[ny][nsys+1], dndympt4[ny], dndysys_modelpt4[ny],dndysyspt4[ny],dndystatpt4[ny]; //lastbin is default
  double dndymdpt4[ny][nsys+1][5]; //lastbin is default
  double my[ny];
  double par[10];

  TH1F* hpad2 = new TH1F("hpad2", "hpad2;p_{T} (GeV/c); B.r.#timesd^{2}N/2#pidp_{T}dy", 10, 0, 3.5);
  hpad2->GetYaxis()->SetRangeUser(5e-11, 5e-2);
  hpad2->Draw();
  gPad->SetLogy(1);
  for (int ic=0;ic<ny;ic++){
    my[ic]=(yedge[ic]+yedge[ic+1])*0.5;
    // cc->Clear();
    // cc->Divide(4,4);
    for (int is=0;is<nsys+1;is++)
    // for (int is=nsys;is<nsys+1;is++)
    {
      // cout <<is<<" "<<ic <<endl;
      // cc->cd(is+1);
      hpad2->GetYaxis()->SetRangeUser(1e-7*pow(10, -1*ic), 1e-2*pow(10, -1*ic));
      hpad2->Draw();
      gPad->SetLogy(1);

      h[ic][is]->Draw("same");
      // g2b[ic]->Draw("p same");
      h[ic][is]->Fit(bolt[ic], "B");
      // g2b[ic]->Fit(bolt[ic], "B");
      bolt[ic]->GetParameters(par);
      ptbolt[ic]->SetParameters(par[0], par[1], _pmass);
      bolt[ic]->Draw("same");

      if (is<nsys) h[ic][is]->Fit(mtexp[ic], "B");
      else if (is==nsys) {
        TFitResultPtr r = h[ic][is]->Fit(mtexp[ic], "BS");
        dndystat[ic] = ptmtexp[ic]->IntegralError(0,10,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() )/scale[ic];
        dndystatpt4[ic] = ptmtexp[ic]->IntegralError(0.4*3,10,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() )/scale[ic];
      }
      // g2b[ic]->Fit(mtexp[ic], "B");
      mtexp[ic]->GetParameters(par);
      ptmtexp[ic]->SetParameters(par[0], par[1], _pmass);
      mtexp[ic]->Draw("same");

      h[ic][is]->Fit(pt2exp[ic], "B");
      // g2b[ic]->Fit(pt2exp[ic], "B");
      pt2exp[ic]->GetParameters(par);
      ptpt2exp[ic]->SetParameters(par[0], par[1], _pmass);
      pt2exp[ic]->Draw("same");

      h[ic][is]->Fit(pt3exp[ic], "B");
      // g2b[ic]->Fit(pt3exp[ic], "B");
      pt3exp[ic]->GetParameters(par);
      ptpt3exp[ic]->SetParameters(par[0], par[1]);
      pt3exp[ic]->Draw("same");

      h[ic][is]->Fit(bw[ic],"B"); 
      // g2b[ic]->Fit(bw[ic],"B"); 
      bw[ic]->GetParameters(par);
      ptbw[ic]->SetParameters(_pmass,bw[ic]->GetParameter(1),bw[ic]->GetParameter(2),bw[ic]->GetParameter(3),bw[ic]->GetParameter(4)*2*TMath::Pi());
      bw[ic]->Draw("same");

      //integral dndydpt in measured region
      dndymd[ic][is][0] = ptbolt[ic]->Integral( 0, 10)/scale[ic];
      dndymd[ic][is][1] = ptmtexp[ic]->Integral( 0, 10)/scale[ic];
      dndymd[ic][is][2] = ptpt2exp[ic]->Integral( 0, 10)/scale[ic];
      dndymd[ic][is][3] = ptpt3exp[ic]->Integral( 0, 10)/scale[ic];
      dndymd[ic][is][4] = ptbw[ic]->Integral( 0, 10)/scale[ic];

      dndymdpt4[ic][is][0] = ptbolt[ic]->Integral( 0.4*3, 10)/scale[ic];
      dndymdpt4[ic][is][1] = ptmtexp[ic]->Integral( 0.4*3, 10)/scale[ic];
      dndymdpt4[ic][is][2] = ptpt2exp[ic]->Integral( 0.4*3, 10)/scale[ic];
      dndymdpt4[ic][is][3] = ptpt3exp[ic]->Integral( 0.4*3, 10)/scale[ic];
      dndymdpt4[ic][is][4] = ptbw[ic]->Integral( 0.4*3, 10)/scale[ic];

      
      dndy[ic][is]=0;
      dndypt4[ic][is]=0;
      // for (int im=0;im<nmodel;im++)
      // {
      //    dndy[ic][is]+=dndymd[ic][is][im]; 
      // }
      // dndy[ic][is]/=(1.*nmodel); //mean as default value, the 0.5*(max-min) as uncertainty
      dndy[ic][is]=dndymd[ic][is][1]; //mt as default?
      dndypt4[ic][is]=dndymdpt4[ic][is][1]; //mt as default?
      if (is<nsys) drawLatex( 0.2,0.2, Form("sys-%d %0.2f<y<%0.2f", is, yedge[ic],yedge[ic+1]), 0.055);
      else drawLatex( 0.2,0.2, Form("%0.2f<y<%0.2f", yedge[ic],yedge[ic+1]), 0.055);
      drawLatex( 0.2,0.25, Form("dn/dy=%f#times10^{-3}", dndy[ic][is]*1e3), 0.055);
      drawLatex( 0.2,0.35, Form("%s%s", centname2[cent].Data(), "%"), 0.055);
      TLegend* lf = new TLegend(0.7,0.65,0.9,0.9);
      lf->AddEntry(mtexp[ic],"m_{T} exp" ,"l" );
      lf->AddEntry(bolt[ic],"Boltzmann" ,"l" );
      lf->AddEntry(bw[ic],"Blaskwave" ,"l" );
      lf->AddEntry(pt3exp[ic],"p_{T}^{1.5}","l" );
      lf->AddEntry(pt2exp[ic],"p_{T} Gaus","l" );
      lf->Draw();
    }

    cc->cd();  
    addpdf(pdf);
    dndysys[ic]=0;
    dndysyspt4[ic]=0;
    // for (int is=0;is<nsys;is++)
    // {
    //    dndysys[ic]+=pow( dndy[ic][is]-dndy[ic][nsys], 2)*0.5;
    //    dndysyspt4[ic]+=pow( dndypt4[ic][is]-dndypt4[ic][nsys], 2)*0.5;
    // }
    double dndysys_topo=0, dndysys_track=0, dndysys_counts=0;
     for (int is=0;is<nsys;is++)
     {
        dndysys[ic]+=pow( dndy[ic][is]-dndy[ic][nsys], 2)*0.5;
        dndysyspt4[ic]+=pow( dndypt4[ic][is]-dndypt4[ic][nsys], 2)*0.5;
        dndysys_topo+=pow( dndypt4[ic][is]-dndypt4[ic][nsys], 2)*0.5;
     }
     // {
     //    dndysys[ic]+=pow( dndy[ic][13]-dndy[ic][nsys], 2);
     //    dndysyspt4[ic]+=pow( dndypt4[ic][13]-dndypt4[ic][nsys], 2);
     //    dndysys_track+=pow( dndypt4[ic][13]-dndypt4[ic][nsys], 2);
     // }
     // {
     //    dndysys[ic]+=pow( dndy[ic][14]-dndy[ic][nsys], 2);
     //    dndysyspt4[ic]+=pow( dndypt4[ic][14]-dndypt4[ic][nsys], 2);
     //    dndysys_counts+=pow( dndypt4[ic][14]-dndypt4[ic][nsys], 2);
     // }
    double tmplo=0, tmphi=0;
    tmplo = min( min( min( min( dndymd[ic][nsys][0], dndymd[ic][nsys][1]), dndymd[ic][nsys][2]), dndymd[ic][nsys][3]),  dndymd[ic][nsys][4]);
    tmphi = max( max( max( max( dndymd[ic][nsys][0], dndymd[ic][nsys][1]), dndymd[ic][nsys][2]), dndymd[ic][nsys][3]),  dndymd[ic][nsys][4]);
    double tmplopt4 = min( min( min( min( dndymdpt4[ic][nsys][0], dndymdpt4[ic][nsys][1]), dndymdpt4[ic][nsys][2]), dndymdpt4[ic][nsys][3]),  dndymdpt4[ic][nsys][4]);
    double tmphipt4 = max( max( max( max( dndymdpt4[ic][nsys][0], dndymdpt4[ic][nsys][1]), dndymdpt4[ic][nsys][2]), dndymdpt4[ic][nsys][3]),  dndymdpt4[ic][nsys][4]);

    dndysys[ic]=sqrt(dndysys[ic] + pow( 0.5*(tmphi-tmplo), 2));
    dndysyspt4[ic]=sqrt(dndysyspt4[ic] + pow( 0.5*(tmphipt4-tmplopt4), 2));
    // dndysys[ic]=sqrt(dndysys[ic]);
    dndym[ic]=dndy[ic][nsys];
    dndympt4[ic]=dndypt4[ic][nsys];
    dndysys_model[ic]=sqrt(pow( 0.5*(tmphi-tmplo), 2));
    dndysys_modelpt4[ic]=sqrt(pow( 0.5*(tmphipt4-tmplopt4), 2));
     cout << "yields at pT/A>0.4: topo: y="<< ic<<" "<<sqrt(dndysys_topo)/dndypt4[ic][nsys]<< " track: "<<sqrt(dndysys_track)/dndypt4[ic][nsys]<< " counts: "<<sqrt(dndysys_counts)/dndypt4[ic][nsys]<<" model: "<<dndysys_modelpt4[ic]/dndypt4[ic][nsys] <<" total: "<< dndysyspt4[ic]/dndypt4[ic][nsys] <<endl;
  }
  for (int ic=0;ic<ny;ic++)
  {
    cout << "parameters:  "<<  my[ic]<<" " <<ptmtexp[ic]->GetParameter(0)<<", "<< ptmtexp[ic]->GetParameter(1)/scale[ic]<< ", "<<  ptmtexp[ic]->GetParameter(2)<< endl;
  }

  gPad->SetLogy(0);
  TH1F* hpad3 = new TH1F("hpad3", "hpad3;Rapidity;B.R.#timesdN/dy", 4, -1, 0);
  hpad3->Draw();
  hpad3->GetYaxis()->SetMaxDigits(2);
  hpad3->GetYaxis()->SetTitleOffset(0.8);
  hpad3->GetXaxis()->SetTitleOffset(1);
  // hpad3->GetYaxis()->SetRangeUser(0,3e-2);
  hpad3->GetYaxis()->SetRangeUser(0,40e-3);
  hpad3->GetYaxis()->SetTitleSize(0.06);
  hpad3->GetXaxis()->SetTitleSize(0.06);
  hpad3->GetYaxis()->SetLabelSize(0.06);
  hpad3->GetXaxis()->SetLabelSize(0.06);
  hpad3->GetXaxis()->SetNdivisions(206);
  hpad3->GetYaxis()->SetNdivisions(206);

  // TGraphErrors* gdndystat = new TGraphErrors(ny-1, my+1, dndym+1, 0, dndystat+1); 
  // TGraphErrors* gdndysys = new TGraphErrors(ny-1, my+1, dndym+1, 0, dndysys+1); 
  TGraphErrors* gdndystat = new TGraphErrors(ny, my, dndym, 0, dndystat); 
  TGraphErrors* gdndysys = new TGraphErrors(ny, my, dndym, 0, dndysys); 
  gdndysys->SetName("gdndysys_3b"); 
  gdndystat->SetName("gdndystat_3b"); 

  TGraphErrors* gdndystatpt4 = new TGraphErrors(ny, my, dndympt4, 0, dndystatpt4); 
  TGraphErrors* gdndysyspt4 = new TGraphErrors(ny, my, dndympt4, 0, dndysyspt4); 
  gdndysyspt4->SetName("gdndysys_3b_pt4"); 
  gdndystatpt4->SetName("gdndystat_3b_pt4"); 

  // drawGraphWithSys(gdndystat, gdndysys, kPink+3, kFullCircle, 1.5, 1, 0);
  drawGraphWithSys(gdndystat, gdndysys, kPink+2, kFullCircle, 1.5, 1, 0);
  //
  // TGraphErrors* gdndysys_md = new TGraphErrors(ny-1, my+1, dndym+1, 0, dndysys_model+1); 
  TGraphErrors* gdndysys_md = new TGraphErrors(ny, my, dndym, 0, dndysys_model); 
  drawGraphWithSys(gdndystat, gdndysys_md, kPink+3, kFullCircle, 1.5, 1, 0);

  //10-50
  // double my2b[2]={-0.125, -0.375};
  // double val2b[2]={ 0.000491959, 0.000883598};
  // double stat2b[2]={ 0.000103048, 0.000136774};
  // double sys2b[2]={ 0.000225092, 0.000335826};
  double yshift = 0.02;
  // double yshift = 0.0;
  double my2b[2]={-0.125+yshift, -0.375+yshift};
  double val2b[2]={ 0.00282647, 0.00348034};
  double stat2b[2]={ 0.000523963, 0.000960732};
  double sys2b[2]={ 0.000605459, 0.00105419};
  TGraphErrors* g2bdndystat = new TGraphErrors( 2, my2b, val2b, 0, stat2b);
  TGraphErrors* g2bdndysys = new TGraphErrors( 2, my2b, val2b, 0, sys2b);
  g2bdndystat->SetName("gdndystat_2b");
  g2bdndysys->SetName("gdndysys_2b");
  graphscale( g2bdndysys, (1-R3)/R3);  
  graphscale( g2bdndystat, (1-R3)/R3);  
  drawGraphWithSys(g2bdndystat, g2bdndysys, kBlue+3, kFullCircle, 1.5, 1, 0);
  drawLatex( 0.2, 0.2, Form("%s%s", centname2[cent].Data(),"%"),0.055);

  TLegend* ldndy = new TLegend( 0.2,0.7,0.5,0.9);
  ldndy->AddEntry( g2bdndysys, " ^{3}_{#Lambda}H#rightarrow ^{3}He#pi#times(1-R_{3})/R_{3}", "pfe");
  ldndy->AddEntry( gdndysys_md, " ^{3}_{#Lambda}H#rightarrow dp#pi(sys.: model only)","pfe");
  ldndy->AddEntry( gdndysys, " ^{3}_{#Lambda}H#rightarrow dp#pi(sys.: model+Topo.)", "pfe");
  // ldndy->AddEntry( gdndysys, " ^{3}_{#Lambda}H#rightarrow dp#pi", "pfe");
  ldndy->Draw();

  addpdf(pdf);

  // drawLatex( );
  graphscale(g2bdndysys,  1./((1-R3)*0.9822*2./3.));
  graphscale(g2bdndystat, 1./((1-R3)*0.9822*2./3.));
  graphscale( gdndysys, 1./((1-R3)*0.9822*2./3.));
  graphscale( gdndystat, 1./((1-R3)*0.9822*2./3.));
  graphscale( gdndysyspt4, 1./((1-R3)*0.9822*2./3.));
  graphscale( gdndystatpt4, 1./((1-R3)*0.9822*2./3.));

  TGraphErrors* gcombdndystat  = (TGraphErrors*)combinechannels(g2bdndystat, gdndystat,"gcombdndystat");
  TGraphErrors* gcombdndysys  = (TGraphErrors*)combinechannels_sys(g2bdndysys, gdndysys, gcombdndystat, "gcombdndysys");

  hpad3->GetYaxis()->SetTitle("dN/dy");
  hpad3->Draw();
  drawGraphWithSys(g2bdndystat, g2bdndysys, kBlue+3, kFullCircle, 1.5, 1, 0);
  drawGraphWithSys(gcombdndystat, gcombdndysys, kGreen+2, kFullCircle, 1.5, 1, 0);
  drawGraphWithSys(gdndystat, gdndysys, kPink+2, kFullCircle, 1.5, 1, 0);
  TGraph* gjam0_10 = new TGraph("../model/H3L2b_0_10cent.csv");
  gjam0_10->SetName("gjam0_10");
  graphscale( gjam0_10, 1.*0.001/((R3)*0.98*2./3.));
  gjam0_10->SetLineStyle(4);
  gjam0_10->SetLineWidth(2);
  gjam0_10->Draw("l same");

  TLegend* ldndycom = new TLegend( 0.2,0.7,0.5,0.9);
  ldndycom->SetTextSize(0.06);
  ldndycom->AddEntry( g2bdndysys, "{}^{3}_{#Lambda}H#rightarrow ^{3}He#pi", "pfe");
  ldndycom->AddEntry( gdndysys, " ^{3}_{#Lambda}H#rightarrow dp#pi","pfe");
  // ldndycom->AddEntry( gdndysys_md, " ^{3}_{#Lambda}H#rightarrow dp#pi(sys.: model only)","pfe");
  // ldndycom->AddEntry( gdndysys, " ^{3}_{#Lambda}H#rightarrow dp#pi(sys.: model+Topo.)", "pfe");
  ldndycom->AddEntry( gcombdndysys, "comb. ^{3}_{#Lambda}H#rightarrow ^{3}He#pi and ^{3}_{#Lambda}H#rightarrow dp#pi", "pfe");
  ldndycom->AddEntry( gjam0_10, "coales.(JAM)", "l");
  ldndycom->Draw();
  // drawSTAR( 0.6, 0.2,0.06);
  drawLatex( 0.7, 0.28, "0-10%", 0.06);
  drawLatex( 0.2, 0.2, "Au+Au #sqrt{s_{NN}} = 3 GeV", 0.06);
  gPad->SaveAs("dNdy_0_10.pdf");

  pdf->On();
  pdf->Close();

  TFile* fresult = new TFile(Form("fsysYield_H3L_cent%s.root", centname[cent].Data()), "recreate");
  for (int iy=0;iy<ny;iy++)
  {
    gYield[iy]->Write();
    gYieldSys[iy]->Write();
    gYieldToterr[iy]->Write();
    mtexp[iy]->SetParameter(1, mtexp[iy]->GetParameter(1)*1./((1-R3)*0.98*2./3.)); //times Br
    mtexp[iy]->Write();
    h[iy][nsys]->Write(); // before times BR
  }
  gdndysys->Write();
  gdndystat->Write();
  gdndystatpt4->Write();
  gdndysyspt4->Write();
  gjam0_10->Write();
  g2bdndysys->Write();
  g2bdndystat->Write();
  gcombdndysys->Write();
  gcombdndystat->Write();

  fresult->Close();
  
  // leg->Draw();
}
void drawEfficiency()
{
  SetsPhenixStyle();
  gStyle->SetPalette(1);
  // TFile* f = TFile::Open("outfile/fout_0010_sys_15.root");
  TFile* f = TFile::Open("outfile/fYield_0010_sys_15.root");
  TH2F* h2d = (TH2F*)f->Get("h2Eff");
  h2d->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  gStyle->SetPaintTextFormat("1.4f");
  h2d->Draw("col text");
  drawLatex(0.2,0.2,"0-10%",0.055);
  gPad->SaveAs("plots/eff_0_10.pdf");
  // return;
  int col[]={kViolet+1,kMagenta,kRed,kOrange+1,kYellow+1,kGreen+2,kBlue,kCyan+1,kViolet-1,kPink-1,kBlack};
  TLegend* leg = new TLegend(0.2,0.2,0.6,0.6); 
  TCanvas* c =new TCanvas("c","c",1000,400);
  c->Divide(3,1);
  TH1F* hpad = new TH1F("hpad","hpad;p_{T} (GeV/c);Efficiency",1,0.5,3.5);
  int color[3]={kRed,kGreen+2,kBlue};
  TH1F* heff[3];
  for (int i=0;i<3;i++) {
    c->cd(i+1);
    hpad->GetYaxis()->SetRangeUser(0,0.16);
    hpad->Draw("c");
    heff[i] = (TH1F*)f->Get(Form("heff%d",i));
    setHistStyle(heff[i], color[i], kFullCircle, 1.5);
    heff[i]->Draw("L same");
    drawLatex(0.2,0.78,Form("%0.2f<y<%0.2f",-1*0.25*(i+1), -1*0.25*(i)),0.055);
    drawLatex(0.2,0.88,"Au+Au 3 GeV",0.055);
    drawLatex(0.3,0.83,"0-10%",0.055);
  }
  // c->Divide(3,4);
  // for (int i=0;i<h2d->GetXaxis()->GetNbins();i++)
  // {
  //   c->cd(i+1);
  //   TH1F* h = (TH1F*)h2d->ProjectionY("h2d",i+1,i+1);
  //   h->SetLineColor(col[i]);
  //   h->SetLineWidth(2);
  //   h->SetMarkerSize(1.);
  //   h->SetMarkerStyle(kOpenCircle);
  //   h->SetMarkerColor(col[i]);
  //   // h->GetYaxis()->SetRangeUser(0.,h->GetMaximum()*1.2);
  //   // if (i>0) h->DrawCopy("p c same");
  //   // if (i==0) h->DrawCopy("p c");
  //   h->GetYaxis()->SetTitle("Efficiency");
  //   h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //   h->DrawCopy("C");
  //   // leg->AddEntry(h, Form("%0.1f<y<%0.1f",h2d->GetXaxis()->GetBinLowEdge(i+1), h2d->GetXaxis()->GetBinLowEdge(i+1)+h2d->GetXaxis()->GetBinWidth(1)), "leg");
  //   drawLatex(0.2,0.88,Form("%0.1f<y<%0.1f",h2d->GetXaxis()->GetBinLowEdge(i+1), h2d->GetXaxis()->GetBinLowEdge(i+1)+h2d->GetXaxis()->GetBinWidth(1)),0.055);
  //   drawLatex(0.2,0.8,"0-50%",0.055);
  //   drawLatex(0.2,0.72,"{}_{#Lambda}^{3}H#rightarrow dp#pi^{-}",0.055);
  //   // drawLatex(0.2,0.9,Form("%0.1f<p_{T}<%0.1f GeV/c",h2d->GetYaxis()->GetBinLowEdge(i+1), h2d->GetYaxis()->GetBinLowEdge(i+1)+h2d->GetYaxis()->GetBinWidth(1+1)),0.05);
  // }
  // c->Divide(3,2);
  // for (int i=1;i<h2d->GetYaxis()->GetNbins()-1;i++)
  // {
  //   c->cd(i);
  //   TH1F* h = (TH1F*)h2d->ProjectionX("h2d",i+1,i+1);
  //   h->SetLineColor(col[i]);
  //   h->SetLineWidth(2);
  //   h->SetMarkerSize(1.);
  //   h->SetMarkerStyle(kOpenCircle);
  //   h->SetMarkerColor(col[i]);
  //   // h->GetYaxis()->SetRangeUser(0.,h->GetMaximum()*1.2);
  //   // if (i>0) h->DrawCopy("p c same");
  //   // if (i==0) h->DrawCopy("p c");
  //   h->GetYaxis()->SetTitle("Efficiency");
  //   h->DrawCopy("C");
  //   // leg->AddEntry(h, Form("%0.1f<y<%0.1f",h2d->GetXaxis()->GetBinLowEdge(i+1), h2d->GetXaxis()->GetBinLowEdge(i+1)+h2d->GetXaxis()->GetBinWidth(1)), "leg");
  //   // drawLatex(0.2,0.2,Form("%0.1f<y<%0.1f",h2d->GetXaxis()->GetBinLowEdge(i+1), h2d->GetXaxis()->GetBinLowEdge(i+1)+h2d->GetXaxis()->GetBinWidth(1)),0.045);
  //   drawLatex(0.2,0.9,Form("%0.1f<p_{T}<%0.1f GeV/c",h2d->GetYaxis()->GetBinLowEdge(i+1), h2d->GetYaxis()->GetBinLowEdge(i+1)+h2d->GetYaxis()->GetBinWidth(1+1)),0.05);
  // }
  // leg->Draw();
  c->cd();
  gPad->SaveAs("plots/eff_0_10_small.pdf");
}
void drawData()
{
  YieldSysCuts(1); 
  YieldSysCuts(0); 
  // qaplots();
}

