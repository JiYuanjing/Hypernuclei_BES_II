#include "TF1.h"
#include "TMath.h"
TF1 *fLastFunc;
Double_t IntegrandBG(const double * x, const double* p){
  // integrand for boltzman-gibbs blast wave
     // x[0] -> r (radius)
     // p[0] -> mass
     // p[1] -> pT (transverse momentum)
     // p[2] -> beta_max (surface velocity)
     // p[3] -> T (freezout temperature)
     // p[4] -> n (velocity profile)


  double x0 = x[0];

  double mass     = p[0];
  double pT       = p[1];
  double beta_max = p[2];
  double temp     = p[3];
  Double_t n      = p[4];

  // Keep beta within reasonable limits
  Double_t beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;

  double mT      = TMath::Sqrt(mass*mass+pT*pT);

  double rho0   = TMath::ATanH(beta);
  double arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  double arg01 = mT*TMath::CosH(rho0)/temp;
  double f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);

  //  printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", x0, pT, beta_max, temp, n, mT, beta, rho0, arg00, arg01);

  return f0;
}


Double_t StaticBGdNdPt(const double * x, const double* p) {

  // implementation of BGBW (1/pt dNdpt)

  double pT = x[0];;


  double mass    = p[0];
  double beta    = p[1];
  double temp    = p[2];
  double n       = p[3];
  double norm    = p[4];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", IntegrandBG, 0, 1, 5);

  fIntBG->SetParameters(mass, pT, beta, temp,n);
  double result = fIntBG->Integral(0,1);
  //  printf ("[%4.4f], Int :%f\n", pT, result);
  return result*norm;//*1e30;;

}


Double_t StaticBGdNdPtTimesPt(const double * x, const double* p) {
  // BGBW dNdpt implementation
  return x[0]*StaticBGdNdPt(x,p);
}


TF1 *GetBGBWdNdpt(Double_t mass, Double_t beta, Double_t temp,
                              Double_t n, Double_t norm, const char * name){

  // BGBW 1/pt dNdpt

  fLastFunc = new TF1 (name, StaticBGdNdPt, 0.0, 10, 5);
  fLastFunc->SetParameters(mass,beta,temp,n,norm);
  fLastFunc->FixParameter(0,mass);

  fLastFunc->SetParNames("mass", "#beta", "T", "n", "norm");
  fLastFunc->SetLineWidth(1);
  return fLastFunc;

}

TF1 *GetBGBWdNdptTimesPt(Double_t mass, Double_t beta, Double_t temp,
                              Double_t n, Double_t norm, const char * name){


  fLastFunc = new TF1 (name, StaticBGdNdPtTimesPt, 0.0, 10, 5);
  fLastFunc->SetParameters(mass,beta,temp,n,norm);
  fLastFunc->FixParameter(0,mass);

  fLastFunc->SetParNames("mass", "#beta", "T", "n", "norm");
  fLastFunc->SetLineWidth(1);
  return fLastFunc;

}
