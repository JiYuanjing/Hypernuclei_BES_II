//define histograms

const float m_pi = 0.13957;
const float m_p = 0.938272;
const float m_K = 0.493677;

THnSparseF* hpion;
THnSparseF* hproton;
THnSparseF* hkaon;

TH3F* hrc;
TH3F* hmc;
TH2F* hphipt;
TH2F* hphipt_ndedx;

void bookHists()
{
  // {pt,eta,dca,nhits,centnumber,charge};
  const int n = 6;
  Int_t bins[n] = {100, 120, 100, 100, 9, 2}; 
  Double_t xmin[n] = {0., -2.4, 0, 0, -0.5, -1};
  Double_t xmax[n] = {5, 2.4, 5, 100, 8.5, 1};
 
  hpion = new THnSparseF("hpion","hpion;pT;eta;dca;nhits;cent;charge",n,bins,xmin,xmax);
  hkaon = new THnSparseF("hkaon","hkaon;pT;eta;dca;nhits;cent;charge",n,bins,xmin,xmax);
  hproton = new THnSparseF("hproton","hproton;pT;eta;dca;nhits;cent;charge",n,bins,xmin,xmax);

  hphipt = new TH2F("hphipt","hphipt;pt;phi",1000,0,2.5,1200,-3.14,3.14);
  hphipt_ndedx = new TH2F("hphipt_ndedx","hphipt;pt;phi",1000,0,2.5,1200,-3.14,3.14);
  hrc = new TH3F("hrc","hrc;pt;eta;cent",100,0,0.5,120,-2.4,2.4,9,-0.5,8.5);
  hmc = new TH3F("hmc","hmc;pt;eta;cent",100,0,0.5,120,-2.4,2.4,9,-0.5,8.5);
}
void fillQaHistograms(int bparticleid,float pt,float eta,int bnhits,float bdca,int cent,double weight)
{
  // float pt;float eta;int charge;
  double array[6]={pt,eta,bdca,bnhits*1.,cent*1.,bparticleid>0?1.:-1.};
  if (abs(bparticleid)==211)  
  {
    hpion->Fill(array,weight); 
  }
  if (abs(bparticleid)==321)  
  {
    hkaon->Fill(array,weight);
  }
  if (abs(bparticleid)==2212)  
  {
    hproton->Fill(array,weight);
  }
}

// void fillQaHistograms(int bparticleid,float bpx,float bpy,float bpz,int bnhits,float bdca,int cent)
// {
//   // float pt;float eta;int charge;
//   TLorentzVector particle;
//   if (abs(bparticleid==211))  
//   {
//     particle.SetXYZM(px,py,pz,m_pi);
//     float array[n]={particle.Perp(),particle.Eta(),bdca,bnhits,cent,bparticleid>0?1:-1};
//     hpion->Fill(array); 
//   }
//   if (abs(bparticleid==321))  
//   {
//     particle.SetXYZM(px,py,pz,m_K);
//     float array[n]={particle.Perp(),particle.Eta(),bdca,bnhits,cent,bparticleid>0?1:-1};
//     hkaon->Fill(array);
//   }
//   if (abs(bparticleid==2212))  
//   {
//     particle.SetXYZM(px,py,pz,m_p);
//     float array[n]={particle.Perp(),particle.Eta(),bdca,bnhits,cent,bparticleid>0?1:-1};
//     hproton->Fill(array);
//   }
// }

int ConvertGeantToPDG(int geantid)
{
   if (geantid==8) return 211;
   if (geantid==9) return -211;
   if (geantid==11) return 321;
   if (geantid==12) return -321;
   if (geantid==14) return 2212;
   if (geantid==15) return -2212;
   return -999;
}

void Write(TString outfilename)
{
  TFile* fout = new TFile(outfilename, "recreate");
  hpion->Write();
  hkaon->Write();
  hproton->Write();
  hrc->Write();
  hmc->Write();
  hphipt->Write();
  hphipt_ndedx->Write();
  fout->Close();
}
