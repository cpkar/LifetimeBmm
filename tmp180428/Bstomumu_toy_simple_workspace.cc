#include <TROOT.h>
#include "TMath.h"
#include <iostream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <TFile.h>

#include <TTree.h>

#include "RooNumIntConfig.h"
#include <RooRandom.h>

#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "RooHist.h"

#include "RooExponential.h"
#include "RooGaussModel.h"
#include "TCanvas.h"
#include "RooLognormal.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "RooGenericPdf.h"
#include "RooTruthModel.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"
#include <TRandom3.h>
#include "RooRandom.h"
#include "RooPoisson.h"
#include "RooGamma.h"
#include "RooBernstein.h"
#include "RooGaussian.h"
#include "iostream"
#include "fstream"
#include "RooFitResult.h"
#include "RooCBShape.h"
#include "RooConstVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooWorkspace.h"

bool plot0=false;
bool plot1=false;
bool simu=false;
using namespace std;
using namespace RooFit;

double bdyield(RooLognormal bdconstraintch0, RooRealVar *nbd1){
  RooDataSet *tmp_evt1BF;
  tmp_evt1BF = bdconstraintch0.generate(RooArgSet(*nbd1),1);
  nbd1->setVal(tmp_evt1BF->get(0)->getRealValue(nbd1->GetName()));
  delete tmp_evt1BF;
  return nbd1->getVal();
}


double semiyield(RooLognormal semiconstraintch0, RooRealVar *nsemi){
  RooDataSet *tmp_evt1BF;
  tmp_evt1BF = semiconstraintch0.generate(RooArgSet(*nsemi),1);
  nsemi->setVal(tmp_evt1BF->get(0)->getRealValue(nsemi->GetName()));
  delete tmp_evt1BF;
  return nsemi->getVal();
}

void definemass(RooWorkspace* wspace, RooWorkspace* Wsig, RooWorkspace* Wsemi, int era, int channel, double l0, double l1, double l2,double l3, double l4){

  RooRealVar* mass = wspace->var("mass");

  string sigera[4]={"16BF","16GH","12","11"};
  string semiera[4]={"16BFs01","16GHs01","12s01","11s01"};

  //  Wsig->Print();
  //Wsig->var(Form("mean%d_%s_%d",channel,sigera[era].c_str(),channel))->Print();
  // Wsemi->var(Form("Mean_semi_20%s_%d",semiera[era].c_str(),channel))->Print();
  RooRealVar B1(Form("B1_%d_ch%d",era,channel), "B_1comb", l0, 0.1, 0.8);
  RooFormulaVar B2(Form("B2_%d_ch%d",era,channel), "B_2comb", "1.-@0", RooArgList(B1));
  RooBernstein mass_comb(Form("mass_comb_%d_ch%d",era,channel), "mass_comb", *mass, RooArgSet(B1,B2));

  RooRealVar meanBs(Form("meanBs_%d_ch%d",era,channel), "mean",Wsig->var(Form("mean%d_%s_%d",channel,sigera[era].c_str(),channel))->getVal());
  RooRealVar sigmacbbs(Form("sigmacbbs_%d_ch%d",era,channel),"m1",Wsig->var(Form("sigma%d_%s_%d",channel,sigera[era].c_str(),channel))->getVal());
  RooRealVar tailcbbs(Form("tailcbbs_%d_ch%d",era,channel),"",Wsig->var(Form("tail%d_%s_%d",channel,sigera[era].c_str(),channel))->getVal());
  RooRealVar powcbbs(Form("powcbbs_%d_ch%d",era,channel),"",Wsig->var(Form("pow%d_%s_%d",channel,sigera[era].c_str(),channel))->getVal());
  RooCBShape cbsbs(Form("cbsbs_%d_ch%d",era,channel),"Signal Lineshape",*mass,meanBs,sigmacbbs,tailcbbs,powcbbs);

  RooRealVar meanBd(Form("meanBd__%d_ch%d",era,channel), "mean",l1);
  RooRealVar sigmacbbd(Form("sigmacbbd_%d_ch%d",era,channel),"m1",l2);
  RooRealVar tailcbbd(Form("tailcbbd_%d_ch%d",era,channel),"",l3);
  RooRealVar powcbbd(Form("powcbbd_%d_ch%d",era,channel),"",l4);
  RooCBShape cbsbd(Form("cbsbd_%d_ch%d",era,channel),"Signal Lineshape",*mass,meanBd,sigmacbbd,tailcbbd,powcbbd);

  RooRealVar meansemi(Form("meansemi_%d_ch%d",era,channel), "mean",Wsemi->var(Form("Mean_semi_20%s_%d",semiera[era].c_str(),channel))->getVal());
  RooRealVar sigmasemi(Form("sigmasemi_%d_ch%d",era,channel),"sigma1 of Gaussian",Wsemi->var(Form("Sigma_semi_20%s_%d",semiera[era].c_str(),channel))->getVal());
  RooGaussian gausiansemi(Form("gausiansemi_%d_ch%d",era,channel),"gaus1", *mass, meansemi, sigmasemi);

  wspace->import(mass_comb);
  wspace->import(cbsbs);
  wspace->import(cbsbd);
  wspace->import(gausiansemi);

}

void definelifetimeerr(RooWorkspace* wspace,RooWorkspace* Wsig, RooWorkspace* Wsemi, int era, int channel, double cer2, double cer3, double cer4, double cer5){

  RooRealVar* trecoe = wspace->var("trecoe");
  RooRealVar Gmean("Gmean","Gmean",0);
  string sigera[4]={"16BF","16GH","12","11"};
  string semiera[4]={"16BFs01","16GHs01","12s01","11s01"};
  
  RooRealVar GmEsig01 (Form("GmEsig01_%d_ch%d",era,channel),"",Wsig->var(Form("mE%d_%s_%d",channel,sigera[era].c_str(),channel))->getVal()) ;
  RooRealVar GEsig01(Form("GEsig01_%d_ch%d",era,channel),"",Wsig->var(Form("wE%d_%s_%d",channel,sigera[era].c_str(),channel))->getVal());
  RooRealVar sigdlE01(Form("sigdlE01_%d_ch%d",era,channel),"c#tauE (cm)",Wsig->var(Form("dlE%d_%s_%d",channel,sigera[era].c_str(),channel))->getVal());
  RooGaussModel REsig01(Form("REsig01_%d_ch%d",era,channel),"Resolution",*trecoe,GmEsig01, GEsig01);
  RooDecay SigE01(Form("SigE01_%d_ch%d",era,channel),"decay model",*trecoe,sigdlE01,REsig01,RooDecay::SingleSided);

  RooRealVar Galphabkg01(Form("Galphabkg01_%d_ch%d",era,channel),"",cer2) ;
  RooRealVar Gthetabkg01(Form("Gthetabkg01_%d_ch%d",era,channel),"",cer3) ;
  RooGamma Gammabkg01(Form("Gammabkg01_%d_ch%d",era,channel),"",*trecoe,Galphabkg01,Gthetabkg01,Gmean);

  RooRealVar GmEsemi01 (Form("GmEsemi01_%d_ch%d",era,channel),"",Wsemi->var(Form("mE_20%s_%d",semiera[era].c_str(),channel))->getVal()) ;
  RooRealVar GEsemi01(Form("GEsemi01_%d_ch%d",era,channel),"",Wsemi->var(Form("wE_20%s_%d",semiera[era].c_str(),channel))->getVal());
  RooRealVar semidlE01(Form("semidlE01_%d_ch%d",era,channel),"c#tauE (cm)",Wsemi->var(Form("dlE_20%s_%d",semiera[era].c_str(),channel))->getVal());
  RooGaussModel REsemi01(Form("REsemi01_%d_ch%d",era,channel),"Resolution",*trecoe,GmEsemi01, GEsemi01);
  RooDecay SemiE01(Form("SemiE01_%d_ch%d",era,channel),"decay model",*trecoe,semidlE01,REsemi01,RooDecay::SingleSided);

  RooRealVar Galphapeak01(Form("Galphapeak01_%d_ch%d",era,channel),"",cer4);
  RooRealVar Gthetapeak01(Form("Gthetapeak01_%d_ch%d",era,channel),"",cer5);
  RooGamma Gammapeak01(Form("Gammapeak01_%d_ch%d",era,channel),"",*trecoe,Galphapeak01,Gthetapeak01,Gmean);

  wspace->import(Gammabkg01);
  wspace->import(SemiE01);
  wspace->import(SigE01);
  wspace->import(Gammapeak01);
}

void lifetime(RooWorkspace* wspace, RooWorkspace* Wsemi,int era, int channel)
{
  RooRealVar* treco = wspace->var("treco");
  RooRealVar* trecoe = wspace->var("trecoe");
  string sigera[4]={"16BF","16GH","12","11"};
  string semiera[4]={"16BFs01","16GHs01","12s01","11s01"};
  
  RooGaussModel gm2(Form("gm2_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);
  RooGaussModel gm3(Form("gm3_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);
  RooGaussModel gm4(Form("gm4_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);

  RooRealVar  TauBd(Form("TauBd_%d_ch%d",era,channel),"",1.42);
  RooDecay Bkg_Ctaubd (Form("Bkg_Ctaubd_%d_ch%d",era,channel),"",*treco,TauBd,gm2,RooDecay::SingleSided);

  RooRealVar  Tausemi(Form("Tausemi_%d_ch%d",era,channel),"",Wsemi->var(Form("tp1_20%s_%d",semiera[era].c_str(),channel))->getVal());
  RooDecay semi_Ctau(Form("semi_Ctau_%d_ch%d",era,channel),"",*treco,Tausemi,gm4,RooDecay::SingleSided);

  wspace->import(semi_Ctau);
  wspace->import(Bkg_Ctaubd);

  RooGaussModel gm1(Form("gm_s1_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);

  RooRealVar  TauSig3(Form("TauSig3_%d_ch%d",era,channel),"",1.7,0.0,10.0);
  RooDecay sig_Ctau3(Form("sig_Ctau3_%d_ch%d",era,channel),"",*treco,TauSig3,gm1,RooDecay::SingleSided);

  //wspace->import(sig_Ctau3);
  RooFormulaVar*  RooEffForch0 = new RooFormulaVar("RooEffForch0","0.00349-0.000287*treco-0.00739/(1+exp(treco*1.023))", RooArgSet(*treco));
  RooFormulaVar*  RooEffFor0BF = new RooFormulaVar("RooEffFor0BF","0.00482-0.000431*treco-0.0109/(1+exp(treco*0.847))", RooArgSet(*treco));
  RooFormulaVar*  RooEffFor1 = new RooFormulaVar("RooEffFor1","0.00737-0.000698*treco-0.0184/(1+exp(treco*0.859))", RooArgSet(*treco));
  RooFormulaVar*  RooEffFor1BF = new RooFormulaVar("RooEffFor1BF","0.00665-0.000531*treco-0.02186/(1+exp(treco*1.2577))", RooArgSet(*treco));
  RooFormulaVar*  RooEffForch012 = new RooFormulaVar("RooEffForch012","0.003+0.0000701*treco-0.0123/(1+exp(treco*1.482))", RooArgSet(*treco));
  RooFormulaVar*  RooEffForch011 = new RooFormulaVar("RooEffForch011","0.00213+0.000015*treco-0.01457/(1+exp(treco*1.871))", RooArgSet(*treco));
  RooFormulaVar*  RooEffForch112 = new RooFormulaVar("RooEffForch112","0.00209-0.000055*treco-0.0077/(1+exp(treco*1.197))", RooArgSet(*treco));
  RooFormulaVar*  RooEffForch111 = new RooFormulaVar("RooEffForch111","0.000896+0.0000023*treco-0.00316/(1+exp(treco*1.0559))", RooArgSet(*treco));

  // wspace->import(*RooEffForch0);
  // wspace->import(*RooEffFor0BF);
  // wspace->import(*RooEffFor1);
  // wspace->import(*RooEffFor1BF);
  // wspace->import(*RooEffForch012);
  // wspace->import(*RooEffForch112);
  // wspace->import(*RooEffForch011);
  // wspace->import(*RooEffForch111);
    
  if(era==1 && channel==1){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffFor1);wspace->import(CtEffSig);}
  else if(era==0 && channel==1){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffFor1BF);wspace->import(CtEffSig);}
  else if(era==0 && channel==0){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffFor0BF);wspace->import(CtEffSig);}
  else if(era==1 && channel==0){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffForch0);wspace->import(CtEffSig);}
  else if(era==2 && channel==0){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffForch012);wspace->import(CtEffSig);}
  else if(era==2 && channel==1){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffForch112);wspace->import(CtEffSig);}
  else if(era==3 && channel==0){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffForch011);wspace->import(CtEffSig);}
  else if(era==3 && channel==1){ RooEffProd CtEffSig(Form("CtEffSig_%d_ch%d",era,channel),"",sig_Ctau3,*RooEffForch111);wspace->import(CtEffSig);}
  
  
}

void lifetimebkg1(RooWorkspace* wspace, int era, int channel, double bkg3)
{
  RooRealVar* treco = wspace->var("treco");
  RooRealVar* trecoe = wspace->var("trecoe");

  RooGaussModel gm1(Form("gm1_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);

  RooRealVar  TauBkg3(Form("TauBkg3_%d_ch%d",era,channel),"",bkg3,0.0,10.0);
  RooDecay Bkg_Ctau3(Form("Bkg_Ctau3_%d_ch%d",era,channel),"",*treco,TauBkg3,gm1,RooDecay::SingleSided);

  wspace->import(Bkg_Ctau3);

}

void lifetimebkg2(RooWorkspace* wspace, int era, int channel, double bkg1, double bkg2, double fg1)
{
  RooRealVar* treco = wspace->var("treco");
  RooRealVar* trecoe = wspace->var("trecoe");

  RooGaussModel gm1(Form("gm5_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);

  RooRealVar  TauBkg1(Form("TauBkg1_%d_ch%d",era,channel),"",bkg1,0.0,10.0);
  RooDecay Bkg_Ctau1(Form("Bkg_Ctau1_%d_ch%d",era,channel),"",*treco,TauBkg1,gm1,RooDecay::SingleSided);

  RooRealVar  TauBkg2(Form("TauBkg2_%d_ch%d",era,channel),"",bkg2,0.0,10.0);
  RooDecay Bkg_Ctau2(Form("Bkg_Ctau2_%d_ch%d",era,channel),"",*treco,TauBkg2,gm1,RooDecay::SingleSided);
  RooRealVar frg1(Form("frg1_%d_ch%d",era,channel),"",fg1,0.0,1.0);

  wspace->import(frg1);
  wspace->import(Bkg_Ctau2);
  wspace->import(Bkg_Ctau1);
}                           

void yield(RooWorkspace* wspace, int era, int channel, double signal, double bkg, double semi, double peak){

  RooRealVar nbkg(Form("nbkg_%d_ch%d",era,channel),"Background fraction",bkg,0,280);
  RooRealVar Nsig(Form("Nsig_%d_ch%d",era,channel),"signal fraction",signal,0,70);
  RooRealVar nbd1(Form("nbd1_%d_ch%d",era,channel),"Bd",peak,0,20);
  RooRealVar nsemi(Form("nsemi_%d_ch%d",era,channel),"semi",semi,0,520);

  wspace->import(nsemi);
  wspace->import(Nsig);
  wspace->import(nbd1);
  wspace->import(nbkg);

}

void TotalPdfdefinition(RooWorkspace* wspace, int era, int channel)
{
  
  RooAbsPdf* CtEffSig= wspace->pdf(Form("CtEffSig_%d_ch%d",era,channel));
  
  RooRealVar* Nsig=wspace->var(Form("Nsig_%d_ch%d",era,channel));
  RooRealVar* nbkg=wspace->var(Form("nbkg_%d_ch%d",era,channel));
  RooRealVar* nbd1=wspace->var(Form("nbd1_%d_ch%d",era,channel));
  RooRealVar* nsemi=wspace->var(Form("nsemi_%d_ch%d",era,channel));

  RooAbsPdf* cbsbs = wspace->pdf(Form("cbsbs_%d_ch%d",era,channel));
  RooAbsPdf* cbsbd = wspace->pdf(Form("cbsbd_%d_ch%d",era,channel));
  RooAbsPdf* mass_comb = wspace->pdf(Form("mass_comb_%d_ch%d",era,channel));
  RooAbsPdf* gausiansemi = wspace->pdf(Form("gausiansemi_%d_ch%d",era,channel));

  RooAbsPdf* CtEffSig1 = wspace->pdf(Form("CtEffSig_%d_ch%d",era,channel));
  RooAbsPdf* Bkg_Ctaubd = wspace->pdf(Form("Bkg_Ctaubd_%d_ch%d",era,channel));
  RooAbsPdf* Bkg_Ctau3 = wspace->pdf(Form("Bkg_Ctau3_%d_ch%d",era,channel));
  RooAbsPdf* semi_Ctau = wspace->pdf(Form("semi_Ctau_%d_ch%d",era,channel));

  RooAbsPdf* Bkg_Ctau12D = wspace->pdf(Form("Bkg_Ctau1_%d_ch%d",era,channel));
  RooAbsPdf* Bkg_Ctau22D = wspace->pdf(Form("Bkg_Ctau2_%d_ch%d",era,channel));

  RooRealVar* frg12D=wspace->var(Form("frg1_%d_ch%d",era,channel));

  RooAbsPdf* Gammapeak01 = wspace->pdf(Form("Gammapeak01_%d_ch%d",era,channel));
  RooAbsPdf* Gammabkg01 = wspace->pdf(Form("Gammabkg01_%d_ch%d",era,channel));
  RooAbsPdf* SemiE01 = wspace->pdf(Form("SemiE01_%d_ch%d",era,channel));
  RooAbsPdf*  SigE01 = wspace->pdf(Form("SigE01_%d_ch%d",era,channel));
  RooProdPdf BsPdf (Form("BsPdf_%d_ch%d",era,channel),"",RooArgList(*cbsbs,*CtEffSig));
  RooProdPdf  BdPdf (Form("BdPdf_%d_ch%d",era,channel),"",RooArgList(*cbsbd,*Bkg_Ctaubd) );
  RooProdPdf  semiPdf (Form("semiPdf_%d_ch%d",era,channel),"",RooArgList(*gausiansemi,*semi_Ctau) );

  RooProdPdf  BsPdf3D (Form("BsPdf3D_%d_ch%d",era,channel),"",RooArgList(*cbsbs,*CtEffSig,*SigE01));
  RooProdPdf  BdPdf3D (Form("BdPdf3D_%d_ch%d",era,channel),"",RooArgList(*cbsbd,*Bkg_Ctaubd,*Gammapeak01) );
  RooProdPdf  semiPdf3D (Form("semiPdf3D_%d_ch%d",era,channel),"",RooArgList(*gausiansemi,*semi_Ctau,*SemiE01) );

  RooProdPdf  BkgPdf (Form("BkgPdf_%d_ch%d",era,channel),"",RooArgList(*mass_comb,*Bkg_Ctau3) );
  RooProdPdf  BkgPdf3D (Form("BkgPdf3D_%d_ch%d",era,channel),"",RooArgList(*mass_comb,*Bkg_Ctau3,*Gammabkg01) );
  
  RooAddPdf TotPdf3D(Form("TotPdf3D_%d_ch%d",era,channel)," Signal + Bkg Pdf",RooArgList(BsPdf3D,BdPdf3D,BkgPdf3D,semiPdf3D),RooArgList(*Nsig,*nbd1,*nbkg,*nsemi));
  wspace->import(TotPdf3D);
  
  RooAddPdf TotPdf(Form("TotPdf_%d_ch%d",era,channel)," Signal + Bkg Pdf",RooArgList(BsPdf,BdPdf,BkgPdf,semiPdf),RooArgList(*Nsig,*nbd1,*nbkg,*nsemi));
  wspace->import(TotPdf);
}
void TotalPdfdefinition2(RooWorkspace* wspace, int era, int channel, double bkg1, double bkg2, double fg1)
{
  RooRealVar* Nsig=wspace->var(Form("Nsig_%d_ch%d",era,channel));
  RooRealVar* nbkg=wspace->var(Form("nbkg_%d_ch%d",era,channel));
  RooRealVar* nbd1=wspace->var(Form("nbd1_%d_ch%d",era,channel));
  RooRealVar* nsemi=wspace->var(Form("nsemi_%d_ch%d",era,channel));

  RooAbsPdf* cbsbs = wspace->pdf(Form("cbsbs_%d_ch%d",era,channel));
  RooAbsPdf* cbsbd = wspace->pdf(Form("cbsbd_%d_ch%d",era,channel));
  RooAbsPdf* mass_comb = wspace->pdf(Form("mass_comb_%d_ch%d",era,channel));
  RooAbsPdf* gausiansemi = wspace->pdf(Form("gausiansemi_%d_ch%d",era,channel));

  RooAbsPdf* CtEffSig1 = wspace->pdf(Form("CtEffSig_%d_ch%d",era,channel));
  RooAbsPdf* Bkg_Ctaubd = wspace->pdf(Form("Bkg_Ctaubd_%d_ch%d",era,channel));
  RooAbsPdf* semi_Ctau = wspace->pdf(Form("semi_Ctau_%d_ch%d",era,channel));

  RooAbsPdf* Gammapeak01 = wspace->pdf(Form("Gammapeak01_%d_ch%d",era,channel));
  RooAbsPdf* Gammabkg01 = wspace->pdf(Form("Gammabkg01_%d_ch%d",era,channel));
  RooAbsPdf* Gammasemi01 = wspace->pdf(Form("SemiE01_%d_ch%d",era,channel));  

  RooProdPdf BsPdf (Form("BsPdf_%d_ch%d",era,channel),"",RooArgList(*cbsbs,*CtEffSig1));
  RooProdPdf  BdPdf (Form("BdPdf_%d_ch%d",era,channel),"",RooArgList(*cbsbd,*Bkg_Ctaubd) );
  RooProdPdf  semiPdf (Form("semiPdf_%d_ch%d",era,channel),"",RooArgList(*gausiansemi,*semi_Ctau) );

  RooProdPdf  BdPdf3D (Form("BdPdf3D_%d_ch%d",era,channel),"",RooArgList(*cbsbd,*Bkg_Ctaubd,*Gammapeak01) );
  RooProdPdf  semiPdf3D (Form("semiPdf3D_%d_ch%d",era,channel),"",RooArgList(*gausiansemi,*semi_Ctau,*Gammasemi01) );

  RooRealVar* treco = wspace->var("treco");
  RooRealVar* trecoe = wspace->var("trecoe");

  RooGaussModel gm1(Form("gm5_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);
  RooGaussModel gm2(Form("gm6_%d_ch%d",era,channel),"gm", *treco,RooConst(0),RooConst(1),*trecoe);
  RooRealVar  TauBkg1(Form("TauBkg1_%d_ch%d",era,channel),"",bkg1,0.0,10.0);
  RooDecay Bkg_Ctau1(Form("Bkg_Ctau1_%d_ch%d",era,channel),"",*treco,TauBkg1,gm1,RooDecay::SingleSided);

  RooRealVar  TauBkg2(Form("TauBkg2_%d_ch%d",era,channel),"",bkg2,0.0,10.0);
  RooDecay Bkg_Ctau2(Form("Bkg_Ctau2_%d_ch%d",era,channel),"",*treco,TauBkg2,gm1,RooDecay::SingleSided);
  RooRealVar frg1(Form("frg1_%d_ch%d",era,channel),"",fg1,0.0,1.0);
  // 2nd time
  RooRealVar  TauBkg5(Form("TauBkg5_%d_ch%d",era,channel),"",bkg1,0.0,10.0);
  RooDecay Bkg_Ctau5(Form("Bkg_Ctau5_%d_ch%d",era,channel),"",*treco,TauBkg5,gm2,RooDecay::SingleSided);

  RooRealVar  TauBkg6(Form("TauBkg6_%d_ch%d",era,channel),"",bkg2,0.0,10.0);
  RooDecay Bkg_Ctau6(Form("Bkg_Ctau6_%d_ch%d",era,channel),"",*treco,TauBkg6,gm2,RooDecay::SingleSided);
  RooRealVar frg4(Form("frg4_%d_ch%d",era,channel),"",fg1,0.0,1.0);

  RooAbsPdf* SigE01 = wspace->pdf(Form("SigE01_%d_ch%d",era,channel));
  
  RooProdPdf  BkgPdf1 (Form("BkgPdf1_%d_ch%d",era,channel),"",RooArgList(*mass_comb,Bkg_Ctau5) );
  RooProdPdf  BkgPdf2 (Form("BkgPdf2_%d_ch%d",era,channel),"",RooArgList(*mass_comb,Bkg_Ctau6) );
  RooAddPdf   BkgPdf (Form("BkgPdf_%d_ch%d",era,channel),"",RooArgSet(BkgPdf1,BkgPdf2),frg4); 
  
  RooProdPdf  BsPdf3D (Form("BsPdf3D_%d_ch%d",era,channel),"",RooArgList(*cbsbs,*CtEffSig1,*SigE01));
  RooProdPdf  BkgPdf3D1 (Form("BkgPdf3D1_%d_ch%d",era,channel),"",RooArgList(*mass_comb,Bkg_Ctau1,*Gammabkg01) );
  RooProdPdf  BkgPdf3D2 (Form("BkgPdf3D2_%d_ch%d",era,channel),"",RooArgList(*mass_comb,Bkg_Ctau2,*Gammabkg01) );
  RooAddPdf   BkgPdf3D (Form("BkgPdf3D_%d_ch%d",era,channel),"",RooArgSet(BkgPdf3D1,BkgPdf3D2),frg1);
  
  RooAddPdf TotPdf3D(Form("TotPdf3D_%d_ch%d",era,channel)," Signal + Bkg Pdf",RooArgList(BsPdf3D,BdPdf3D,BkgPdf3D,semiPdf3D),RooArgList(*Nsig,*nbd1,*nbkg,*nsemi));
  wspace->import(TotPdf3D);
  RooAddPdf TotPdf(Form("TotPdf_%d_ch%d",era,channel)," Signal + Bkg Pdf",RooArgList(BsPdf,BdPdf,BkgPdf,semiPdf),RooArgList(*Nsig,*nbd1,*nbkg,*nsemi));
  wspace->import(TotPdf);
  
}

void fit_chan(RooWorkspace* wspace,int era, int channel,double bd,double semi){

  RooAbsPdf* TotPdf3D=wspace->pdf(Form("TotPdf3D_%d_ch%d",era,channel));
  RooAbsPdf* TotPdf=wspace->pdf(Form("TotPdf_%d_ch%d",era,channel));

  RooRealVar* nbd1=wspace->var(Form("nbd1_%d_ch%d",era,channel));
  RooRealVar* nsemi=wspace->var(Form("nsemi_%d_ch%d",era,channel));

  RooLognormal bdconstraint(Form("bdconstraint_%d_ch%d",era,channel),"bd constraint",*nbd1,RooConst(bd),RooConst(bd*0.25)) ;
  RooLognormal semiconstraint(Form("semiconstraint_%d_ch%d",era,channel),"semi constraint",*nsemi,RooConst(semi),RooConst(semi*0.15)) ;
  RooArgSet ext_constrch0(bdconstraint,semiconstraint);
  
  RooRealVar* treco = wspace->var("treco");
  RooRealVar* trecoe = wspace->var("trecoe");
  RooRealVar* mass = wspace->var("mass");

  RooAbsReal* Tauch0=wspace->var(Form("TauSig3_%d_ch%d",era,channel));  
  bdyield(bdconstraint,nbd1);
  // cout<<"nbd1 "<<nbd1ch0->getVal()<<endl;
  // cout<<"semi val "<<nsemich0->getVal()<<endl;
  semiyield(semiconstraint,nsemi);
  // cout<<"nsemi1 "<<nsemich0->getVal()<<endl;

  RooDataSet* datach0=TotPdf3D->generate(RooArgSet(*mass,*treco,*trecoe),Extended(kTRUE));
  cout<<"datach0 entry "<<datach0->sumEntries()<<endl;
  //RooFitResult* fit3dch0=TotPdf->fitTo(*datach0,ExternalConstraints(ext_constrch0),Extended(kTRUE),Minos(RooArgSet(*Tauch0)),Save(kTRUE),ConditionalObservables(*trecoe),Strategy(2));
  RooFitResult* fit3dch0=TotPdf->fitTo(*datach0,Extended(kTRUE),Minos(RooArgSet(*Tauch0)),Save(kTRUE),ConditionalObservables(*trecoe),ExternalConstraints(ext_constrch0));
  //fit3dch0->Print("v");
  wspace->import(*Tauch0);//,Rename("TauSig3_%d_ch%",era,channel));
  wspace->import(*datach0,Rename(Form("data_%d_%d",era,channel)));
  
  //wspace->Print();
  delete datach0;
}

void cons_par(RooWorkspace* wspace,int era, int channel){

  RooRealVar* Nsig=wspace->var(Form("Nsig_%d_ch%d",era,channel));
  RooRealVar* nbkg=wspace->var(Form("nbkg_%d_ch%d",era,channel));
  RooRealVar* nbd1=wspace->var(Form("nbd1_%d_ch%d",era,channel));
  RooRealVar* Nsemi=wspace->var(Form("nsemi_%d_ch%d",era,channel));
  
  Nsig->setConstant(kTRUE);
  Nsemi->setConstant(kTRUE);
  nbkg->setConstant(kTRUE);
  nbd1->setConstant(kTRUE);
  RooRealVar* B1=wspace->var(Form("B1_%d_ch%d",era,channel));
  B1->setConstant(kTRUE);

  RooRealVar* TauBkg3=wspace->var(Form("TauBkg3_%d_ch%d",era,channel));
  TauBkg3->setConstant(kTRUE);
}

void cons_par2(RooWorkspace* wspace,int era, int channel){

  RooRealVar* Nsig=wspace->var(Form("Nsig_%d_ch%d",era,channel));
  RooRealVar* nbkg=wspace->var(Form("nbkg_%d_ch%d",era,channel));
  RooRealVar* nbd1=wspace->var(Form("nbd1_%d_ch%d",era,channel));
  RooRealVar* Nsemi=wspace->var(Form("nsemi_%d_ch%d",era,channel));

  Nsig->setConstant(kTRUE);
  Nsemi->setConstant(kTRUE);
  nbkg->setConstant(kTRUE);
  nbd1->setConstant(kTRUE);
  RooRealVar* B1=wspace->var(Form("B1_%d_ch%d",era,channel));
  B1->setConstant(kTRUE);

  RooRealVar* TauBkg5=wspace->var(Form("TauBkg5_%d_ch%d",era,channel));
  RooRealVar* TauBkg6=wspace->var(Form("TauBkg6_%d_ch%d",era,channel));
  RooRealVar* frg4=wspace->var(Form("frg4_%d_ch%d",era,channel));

  TauBkg5->setConstant(kTRUE);
  TauBkg6->setConstant(kTRUE);
  frg4->setConstant(kTRUE);
  
}

void free_par(RooWorkspace* wspace,int era, int channel){

  RooRealVar* Nsig=wspace->var(Form("Nsig_%d_ch%d",era,channel));
  RooRealVar* nbkg=wspace->var(Form("nbkg_%d_ch%d",era,channel));
  RooRealVar* nbd1=wspace->var(Form("nbd1_%d_ch%d",era,channel));
  RooRealVar* Nsemi=wspace->var(Form("nsemi_%d_ch%d",era,channel));

  Nsig->setConstant(kFALSE);
  Nsemi->setConstant(kFALSE);
  nbkg->setConstant(kFALSE);
  nbd1->setConstant(kFALSE);
  RooRealVar* B1=wspace->var(Form("B1_%d_ch%d",era,channel));
  B1->setConstant(kFALSE);
  RooRealVar* TauBkg3=wspace->var(Form("TauBkg3_%d_ch%d",era,channel));
  TauBkg3->setConstant(kFALSE);
  
}
void free_par2(RooWorkspace* wspace,int era, int channel){

  RooRealVar* Nsig=wspace->var(Form("Nsig_%d_ch%d",era,channel));
  RooRealVar* nbkg=wspace->var(Form("nbkg_%d_ch%d",era,channel));
  RooRealVar* nbd1=wspace->var(Form("nbd1_%d_ch%d",era,channel));
  RooRealVar* Nsemi=wspace->var(Form("nsemi_%d_ch%d",era,channel));

  Nsig->setConstant(kFALSE);
  Nsemi->setConstant(kFALSE);
  nbkg->setConstant(kFALSE);
  nbd1->setConstant(kFALSE);
  RooRealVar* B1=wspace->var(Form("B1_%d_ch%d",era,channel));
  B1->setConstant(kFALSE);

  RooRealVar* TauBkg5=wspace->var(Form("TauBkg5_%d_ch%d",era,channel));
  RooRealVar* TauBkg6=wspace->var(Form("TauBkg6_%d_ch%d",era,channel));
  RooRealVar* frg4=wspace->var(Form("frg4_%d_ch%d",era,channel));

  TauBkg5->setConstant(kFALSE);
  TauBkg6->setConstant(kFALSE);
  frg4->setConstant(kFALSE);
}

void produceplot(RooWorkspace* wspace,int era, int channel){

  RooRealVar* treco = wspace->var("treco");
  RooRealVar* trecoe = wspace->var("trecoe");
  RooRealVar* mass = wspace->var("mass");

  RooAbsPdf* TotPdf=wspace->pdf(Form("TotPdf_%d_ch%d",era,channel));
  RooDataSet* data=(RooDataSet*)wspace->data(Form("data_%d_%d",era,channel));
  RooPlot* mframe  = mass->frame(Title("B_{s} mass distribution"),Bins(50));
  data->plotOn(mframe);
  TotPdf->plotOn(mframe);
  TotPdf->plotOn(mframe, Components(*wspace->pdf(Form("mass_comb_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kGreen));
  TotPdf->plotOn(mframe, Components(*wspace->pdf(Form("cbsbd_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kCyan));
  TotPdf->plotOn(mframe, Components(*wspace->pdf(Form("cbsbs_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kRed));
  TotPdf->plotOn(mframe, Components(*wspace->pdf(Form("gausiansemi_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kMagenta));

  RooPlot* tframe  = treco->frame(Title("Decay time"),Bins(40));
  data->plotOn(tframe);
  TotPdf->plotOn(tframe);
  TotPdf->plotOn(tframe, Components(*wspace->pdf(Form("CtEffSig_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kRed));
  TotPdf->plotOn(tframe, Components(*wspace->pdf(Form("Bkg_Ctau3_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kGreen));
  TotPdf->plotOn(tframe, Components(*wspace->pdf(Form("Bkg_Ctaubd_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kCyan));
  TotPdf->plotOn(tframe, Components(*wspace->pdf(Form("semi_Ctau_%d_ch%d",era,channel))), LineStyle(kDashed),LineColor(kMagenta));
  TCanvas* tmm=new TCanvas("tmm","",600,600);
  tmm->Divide(1,2);
  tmm->cd(1);
  mframe->Draw();
  tmm->cd(2);
  tframe->Draw();
  tmm->SaveAs("2.pdf");

}




void simul_fit(RooWorkspace* wspace){
  
  RooRealVar* treco = wspace->var("treco");
  RooRealVar* trecoe = wspace->var("trecoe");
  RooRealVar* mass = wspace->var("mass");
  
  RooDataSet* data0BF=(RooDataSet*)wspace->data("data_0_0");
  RooDataSet* data1BF=(RooDataSet*)wspace->data("data_0_1");
  RooDataSet* data0GH=(RooDataSet*)wspace->data("data_1_0");
  RooDataSet* data1GH=(RooDataSet*)wspace->data("data_1_1");
  RooDataSet* data012=(RooDataSet*)wspace->data("data_2_0");
  RooDataSet* data112=(RooDataSet*)wspace->data("data_2_1");
  RooDataSet* data011=(RooDataSet*)wspace->data("data_3_0");
  RooDataSet* data111=(RooDataSet*)wspace->data("data_3_1");

  RooRealVar* t_wh=new RooRealVar("t_wh","whole",1.70, 0.0, 10.0) ;
  RooGaussModel gm42("gm42","gm", *treco,RooConst(0),RooConst(1),*trecoe);
  
  RooDecay* decaytime8=new RooDecay("decaytime8","decay",*treco,*t_wh,gm42,RooDecay::SingleSided) ;
  RooEffProd* CtSigchs012 =new RooEffProd("CtSigchs012","",*decaytime8,*wspace->function("RooEffForch012"));
  RooEffProd* CtSigchs112 =new RooEffProd("CtSigchs112","",*decaytime8,*wspace->function("RooEffForch112"));
  RooEffProd* CtSigchs011 =new RooEffProd("CtSigchs011","",*decaytime8,*wspace->function("RooEffForch011"));
  RooEffProd* CtSigchs111 =new RooEffProd("CtSigchs111","",*decaytime8,*wspace->function("RooEffForch111"));
  RooEffProd* CtSigchs0BF =new RooEffProd("CtSigchs0BF","",*decaytime8,*wspace->function("RooEffFor0BF"));
  RooEffProd* CtSigchs1BF =new RooEffProd("CtSigchs1BF","",*decaytime8,*wspace->function("RooEffFor1BF"));
  RooEffProd* CtSigchs0GH =new RooEffProd("CtSigchs0GH","",*decaytime8,*wspace->function("RooEffForch0"));
  RooEffProd* CtSigchs1GH =new RooEffProd("CtSigchs1GH","",*decaytime8,*wspace->function("RooEffFor1"));
  
  
  RooProdPdf* SigTotallch011=new RooProdPdf("SigTotallch011","",RooArgList(*wspace->pdf("cbsbs_3_ch0"),*CtSigchs011));
  RooProdPdf* SigTotallch111=new RooProdPdf("SigTotallch111","",RooArgList(*wspace->pdf("cbsbs_3_ch1"),*CtSigchs111));
  RooProdPdf* SigTotallch012=new RooProdPdf("SigTotallch012","",RooArgList(*wspace->pdf("cbsbs_2_ch0"),*CtSigchs012));
  RooProdPdf* SigTotallch112=new RooProdPdf("SigTotallch112","",RooArgList(*wspace->pdf("cbsbs_2_ch1"),*CtSigchs112));

  RooProdPdf* SigTotallch0GH=new RooProdPdf("SigTotallch0GH","",RooArgList(*wspace->pdf("cbsbs_1_ch0"),*CtSigchs0GH));
  RooProdPdf* SigTotallch1GH=new RooProdPdf("SigTotallch1GH","",RooArgList(*wspace->pdf("cbsbs_1_ch1"),*CtSigchs1GH));
  RooProdPdf* SigTotallch0BF=new RooProdPdf("SigTotallch0BF","",RooArgList(*wspace->pdf("cbsbs_0_ch0"),*CtSigchs0BF));
  RooProdPdf* SigTotallch1BF=new RooProdPdf("SigTotallch1BF","",RooArgList(*wspace->pdf("cbsbs_0_ch1"),*CtSigchs1BF));


  RooAddPdf* TotAllchan011=new RooAddPdf("TotAllchan011"," Signal + Bkg Pdf",RooArgList(*SigTotallch011,*wspace->pdf("BdPdf_3_ch0"),*wspace->pdf("BkgPdf_3_ch0"),*wspace->pdf("semiPdf_3_ch0")),RooArgList(*wspace->var("Nsig_3_ch0"),*wspace->var("nbd1_3_ch0"),*wspace->var("nbkg_3_ch0"),*wspace->var("nsemi_3_ch0")));
  RooAddPdf* TotAllchan111=new RooAddPdf("TotAllchan111"," Signal + Bkg Pdf",RooArgList(*SigTotallch111,*wspace->pdf("BdPdf_3_ch1"),*wspace->pdf("BkgPdf_3_ch1"),*wspace->pdf("semiPdf_3_ch1")),RooArgList(*wspace->var("Nsig_3_ch1"),*wspace->var("nbd1_3_ch1"),*wspace->var("nbkg_3_ch1"),*wspace->var("nsemi_3_ch1")));
  
  RooAddPdf* TotAllchan012=new RooAddPdf("TotAllchan012"," Signal + Bkg Pdf",RooArgList(*SigTotallch012,*wspace->pdf("BdPdf_2_ch0"),*wspace->pdf("BkgPdf_2_ch0"),*wspace->pdf("semiPdf_2_ch0")),RooArgList(*wspace->var("Nsig_2_ch0"),*wspace->var("nbd1_2_ch0"),*wspace->var("nbkg_2_ch0"),*wspace->var("nsemi_2_ch0")));
  RooAddPdf* TotAllchan112=new RooAddPdf("TotAllchan112"," Signal + Bkg Pdf",RooArgList(*SigTotallch112,*wspace->pdf("BdPdf_2_ch1"),*wspace->pdf("BkgPdf_2_ch1"),*wspace->pdf("semiPdf_2_ch1")),RooArgList(*wspace->var("Nsig_2_ch1"),*wspace->var("nbd1_2_ch1"),*wspace->var("nbkg_2_ch1"),*wspace->var("nsemi_2_ch1")));

  RooAddPdf* TotAllchan0GH=new RooAddPdf("TotAllchan0GH"," Signal + Bkg Pdf",RooArgList(*SigTotallch0GH,*wspace->pdf("BdPdf_1_ch0"),*wspace->pdf("BkgPdf_1_ch0"),*wspace->pdf("semiPdf_1_ch0")),RooArgList(*wspace->var("Nsig_1_ch0"),*wspace->var("nbd1_1_ch0"),*wspace->var("nbkg_1_ch0"),*wspace->var("nsemi_1_ch0")));
  RooAddPdf* TotAllchan1GH=new RooAddPdf("TotAllchan1GH"," Signal + Bkg Pdf",RooArgList(*SigTotallch1GH,*wspace->pdf("BdPdf_1_ch1"),*wspace->pdf("BkgPdf_1_ch1"),*wspace->pdf("semiPdf_1_ch1")),RooArgList(*wspace->var("Nsig_1_ch1"),*wspace->var("nbd1_1_ch1"),*wspace->var("nbkg_1_ch1"),*wspace->var("nsemi_1_ch1")));
  RooAddPdf* TotAllchan0BF=new RooAddPdf("TotAllchan0BF"," Signal + Bkg Pdf",RooArgList(*SigTotallch0BF,*wspace->pdf("BdPdf_0_ch0"),*wspace->pdf("BkgPdf_0_ch0"),*wspace->pdf("semiPdf_0_ch0")),RooArgList(*wspace->var("Nsig_0_ch0"),*wspace->var("nbd1_0_ch0"),*wspace->var("nbkg_0_ch0"),*wspace->var("nsemi_0_ch0")));
  RooAddPdf* TotAllchan1BF=new RooAddPdf("TotAllchan1BF"," Signal + Bkg Pdf",RooArgList(*SigTotallch1BF,*wspace->pdf("BdPdf_0_ch1"),*wspace->pdf("BkgPdf_0_ch1"),*wspace->pdf("semiPdf_0_ch1")),RooArgList(*wspace->var("Nsig_0_ch1"),*wspace->var("nbd1_0_ch1"),*wspace->var("nbkg_0_ch1"),*wspace->var("nsemi_0_ch1")));

 
  RooCategory sampdfall("sampdfall","sampdf all") ;
  sampdfall.defineType("PdfChan0BF",1) ;
  sampdfall.defineType("PdfChan1BF",2) ;
  sampdfall.defineType("PdfChan0",3) ;
  sampdfall.defineType("PdfChan1",4) ;
  sampdfall.defineType("PdfChan012",5) ;
  sampdfall.defineType("PdfChan112",6) ;
  sampdfall.defineType("PdfChan011",7) ;
  sampdfall.defineType("PdfChan111",8) ;

  sampdfall.setIndex(1);
  data0BF->addColumn(sampdfall);
  sampdfall.setIndex(2);
  data1BF->addColumn(sampdfall);
  sampdfall.setIndex(3);
  data0GH->addColumn(sampdfall);
  sampdfall.setIndex(4);
  data1GH->addColumn(sampdfall);
  sampdfall.setIndex(5);
  data012->addColumn(sampdfall);
  sampdfall.setIndex(6);
  data112->addColumn(sampdfall);
  sampdfall.setIndex(7);
  data011->addColumn(sampdfall);
  sampdfall.setIndex(8);
  data111->addColumn(sampdfall);
  RooDataSet* combData2Dall=new RooDataSet("combData2Dall","combined data",RooArgSet(*mass,*treco,*trecoe,sampdfall));
  combData2Dall->append(*data0BF);
  combData2Dall->append(*data1BF);
  combData2Dall->append(*data0GH);
  combData2Dall->append(*data1GH);
  combData2Dall->append(*data012);
  combData2Dall->append(*data112);
  combData2Dall->append(*data011);
  combData2Dall->append(*data111);

  RooSimultaneous simulPdfall("simulPdfall","simultaneous pdf",sampdfall) ;
  simulPdfall.addPdf(*TotAllchan0BF,"PdfChan0BF") ;
  simulPdfall.addPdf(*TotAllchan1BF,"PdfChan1BF") ;
  simulPdfall.addPdf(*TotAllchan0GH,"PdfChan0") ;
  simulPdfall.addPdf(*TotAllchan1GH,"PdfChan1") ;
  simulPdfall.addPdf(*TotAllchan012,"PdfChan012") ;
  simulPdfall.addPdf(*TotAllchan112,"PdfChan112") ;
  simulPdfall.addPdf(*TotAllchan011,"PdfChan011") ;
  simulPdfall.addPdf(*TotAllchan111,"PdfChan111") ;

  RooFitResult* fitsimall=  simulPdfall.fitTo(*combData2Dall,Save(true),ConditionalObservables(*trecoe),Extended(true),Minos(RooArgSet(*t_wh)));
  //fitsimall->Print("v");
  wspace->import(*t_wh);
  delete combData2Dall;

}


void Bstomumu_toy_simple_workspace(){
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(1,false);
  RooMsgService::instance().getStream(1).removeTopic(Integration) ;
  RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
  RooMsgService::instance().getStream(1).removeTopic(Fitting) ;
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
  RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
  RooMsgService::instance().getStream(1).removeTopic(Eval) ;
  RooMsgService::instance().Print() ;
  double taumeansimall,tauerrorsimall,pullsimall;
  double taumean0,tauerror0,pull0;
  double taumean1,tauerror1,pull1;
  double taumean2,tauerror2,pull2;
  double taumean3,tauerror3,pull3;
  double taumean4,tauerror4,pull4;
  double taumean5,tauerror5,pull5;
  

  TFile *file = new TFile("ToyMC_BsToMuMu_all_workspace_plot.root","RECREATE");
  TTree *tree0 = new TTree("tree0","tree0");
  file->cd();
  tree0->Branch("taumean",&taumean0);
  tree0->Branch("tauerror",&tauerror0);
  tree0->Branch("pull",&pull0);

  TTree *tree1 = new TTree("tree1","tree1");
  file->cd();
  tree1->Branch("taumean",&taumean1);
  tree1->Branch("tauerror",&tauerror1);
  tree1->Branch("pull",&pull1);

  TTree *tree2 = new TTree("tree2","tree2");
  file->cd();
  tree2->Branch("taumean",&taumean2);
  tree2->Branch("tauerror",&tauerror2);
  tree2->Branch("pull",&pull2);

  TTree *tree3 = new TTree("tree3","tree3");
  file->cd();
  tree3->Branch("taumean",&taumean3);
  tree3->Branch("tauerror",&tauerror3);
  tree3->Branch("pull",&pull3);

  TTree *tree4 = new TTree("tree4","tree4");
  file->cd();
  tree4->Branch("taumean",&taumean4);
  tree4->Branch("tauerror",&tauerror4);
  tree4->Branch("pull",&pull4);
  
  TTree *tree5 = new TTree("tree5","tree5");
  file->cd();
  tree5->Branch("taumean",&taumeansimall);
  tree5->Branch("tauerror",&tauerrorsimall);
  tree5->Branch("pull",&pullsimall);
  
  RooRandom::randomGenerator()->SetSeed(1234);
  RooRealVar* mass=new RooRealVar("mass","M_{#mu^{+}#mu^{-}}[ GeV/c^{2}]",4.9,5.9);
  RooRealVar* treco=new RooRealVar("treco","t_{reco}[ps]",1.2,10) ;
  RooRealVar* trecoe=new RooRealVar("trecoe","terr[ps]",0.003,0.26);
  for(int i=0;i<1000;i++){
    Int_t ObjectNumber=TProcessID::GetObjectCount();
    RooWorkspace *wspace = new RooWorkspace("wspace");
    
    wspace->import(*mass);
    wspace->import(*treco);
    wspace->import(*trecoe);
  
    TFile *sigpdfMC = new TFile("workspcae_signal_BFGH.root");
    RooWorkspace* Wsig=(RooWorkspace*)sigpdfMC->Get("wspace");
    TFile *semipdfMC = new TFile("workspcae_semileptonic_BFGH.root");
    RooWorkspace* Wsemi=(RooWorkspace*)semipdfMC->Get("wspace");
    
    double bdinput[8]={0.372,0.24,1.309,0.564,0.607,1.195,0.727,1.045};
    double semiinput[8]={2.027,1.0,9.292,11.328,7.456,5.139,12.997,18.19};
    double Bsinput[8]={3.02,1.785,10.271,4.509,4.541,8.626,4.989,7.501};
    double combinput[8]={6.667,6.667,33.333,17.78,22.22,93.33,31.11,51.11};
  
    
    
    //2016BF channel 0
    
    definemass(wspace, Wsig, Wsemi,0,0, 0.521, 5.277, 0.034, 1.998, 1.191);
    definelifetimeerr(wspace,Wsig,Wsemi,0,0, 15.236, 0.006,12.4,0.004);
    lifetimebkg1(wspace,0, 0,0.585 );
    lifetime(wspace,Wsemi,0, 0);
    yield(wspace,0, 0, Bsinput[4],combinput[4],semiinput[4],bdinput[4]);
    TotalPdfdefinition(wspace,0, 0);
    fit_chan(wspace,0,0,bdinput[4],semiinput[4]);
    cons_par(wspace,0,0);
    //produceplot(wspace,0,0);

    RooRealVar* t_BF=wspace->var("TauSig3_0_ch0");
    cout<<"toy " <<i<<"\t"<<t_BF->getVal()<<"\t"<<t_BF->getErrorHi()<<"\t low "<<t_BF->getErrorLo()<<endl;
    taumean0=t_BF->getVal();
    tauerror0=t_BF->getError();
    double higherror0=fabs(t_BF->getErrorHi());
    double lowerror0=fabs(t_BF->getErrorLo());
    double meanlife0=t_BF->getVal();
    if(higherror0==0)higherror0=t_BF->getError();
    if(lowerror0==0)lowerror0=t_BF->getError();
    if(meanlife0 > 1.70){
      pull0=(meanlife0-1.70)/lowerror0;
      
    }
    else if(meanlife0 < 1.70){
      pull0=(meanlife0-1.70)/higherror0;
    }
    tree0->Fill();

    //2016BF channel 1
    definemass(wspace, Wsig, Wsemi,0,1, 0.51, 5.277, 0.048, 1.818, 1.728);
    definelifetimeerr(wspace,Wsig,Wsemi,0,1, 10.4, 0.01,12.4,0.004);
    lifetime(wspace,Wsemi,0,1);
    yield(wspace,0,1, Bsinput[5],combinput[5],semiinput[5],bdinput[5]);
    TotalPdfdefinition2(wspace,0,1,0.62,4.23,0.75);
    fit_chan(wspace,0,1,bdinput[5],semiinput[5]);
    cons_par2(wspace,0,1);
    RooRealVar* t_BF1=wspace->var("TauSig3_0_ch1");
    cout<<"toy " <<i<<"\t"<<t_BF1->getVal()<<"\t"<<t_BF1->getErrorHi()<<"\t low "<<t_BF1->getErrorLo()<<endl;
    taumean1=t_BF1->getVal();
    tauerror1=t_BF1->getError();
    double higherror1=fabs(t_BF1->getErrorHi());
    double lowerror1=fabs(t_BF1->getErrorLo());
    double meanlife1=t_BF1->getVal();
    if(higherror1==0)higherror1=t_BF1->getError();
    if(lowerror1==0)lowerror1=t_BF1->getError();
    if(meanlife1 > 1.70){
      pull1=(meanlife1-1.70)/lowerror1;

    }
    else if(meanlife1 < 1.70){
      pull1=(meanlife1-1.70)/higherror1;
    }
    tree1->Fill();
    
    // 2016GH channel 0
    definemass(wspace, Wsig, Wsemi,1,0, 0.516, 5.277, 0.034, 1.998, 1.191);
    definelifetimeerr(wspace,Wsig,Wsemi,1,0, 20.8, 0.007,12.4,0.004);
    lifetimebkg1(wspace,1,0,0.658 );
    lifetime(wspace,Wsemi,1,0);
    yield(wspace,1,0,  Bsinput[6],combinput[6],semiinput[6],bdinput[6]);
    TotalPdfdefinition(wspace,1,0);
    fit_chan(wspace,1,0,bdinput[6],semiinput[6]);
    cons_par(wspace,1,0);
    RooRealVar* t_GH=wspace->var("TauSig3_1_ch0");
    cout<<"toy " <<i<<"\t"<<t_GH->getVal()<<"\t"<<t_GH->getErrorHi()<<"\t low "<<t_GH->getErrorLo()<<endl;
    taumean2=t_GH->getVal();
    tauerror2=t_GH->getError();
    double higherror2=fabs(t_GH->getErrorHi());
    double lowerror2=fabs(t_GH->getErrorLo());
    double meanlife2=t_GH->getVal();
    if(higherror2==0)higherror2=t_GH->getError();
    if(lowerror2==0)lowerror2=t_GH->getError();
    if(meanlife2 > 1.70){
      pull2=(meanlife2-1.70)/lowerror2;
      
    }
    else if(meanlife2 < 1.70){
      pull2=(meanlife2-1.70)/higherror2;
    }
    tree2->Fill();
    
    
    //2016GH channel1
    definemass(wspace, Wsig, Wsemi,1,1, 0.46, 5.277, 0.048, 1.818, 1.728);
    definelifetimeerr(wspace,Wsig,Wsemi,1,1, 10.98, 0.007,12.4,0.004);
    lifetime(wspace,Wsemi,1,1);
    yield(wspace,1,1,  Bsinput[7],combinput[7],semiinput[7],bdinput[7]);
    TotalPdfdefinition2(wspace,1,1,0.334,2.527,0.664);
    fit_chan(wspace,1,1,bdinput[7],semiinput[7]);
    cons_par2(wspace,1,1);
    RooRealVar* t_GH1=wspace->var("TauSig3_1_ch1");
    cout<<"toy " <<i<<"\t"<<t_GH1->getVal()<<"\t"<<t_GH1->getErrorHi()<<"\t low "<<t_GH1->getErrorLo()<<endl;
    taumean3=t_GH1->getVal();
    tauerror3=t_GH1->getError();
    double higherror3=fabs(t_GH1->getErrorHi());
    double lowerror3=fabs(t_GH1->getErrorLo());
    double meanlife3=t_GH1->getVal();
    if(higherror3==0)higherror3=t_GH1->getError();
    if(lowerror3==0)lowerror3=t_GH1->getError();
    if(meanlife3 > 1.70){
    pull3=(meanlife3-1.70)/lowerror3;
    
    }
    else if(meanlife3 < 1.70){
      pull3=(meanlife3-1.70)/higherror3;
    }
    tree3->Fill();
    
    //2012 channel0
    definemass(wspace, Wsig, Wsemi,2,0, 0.515, 5.277, 0.04, 1.908, 1.429);
    definelifetimeerr(wspace,Wsig,Wsemi,2,0, 13.4, 0.007,12.4,0.004);
    lifetimebkg1(wspace,2, 0,0.429 );
    lifetime(wspace,Wsemi,2, 0);
    yield(wspace,2,0, Bsinput[2],combinput[2],semiinput[2],bdinput[2]);
    TotalPdfdefinition(wspace,2, 0);
    fit_chan(wspace,2,0,bdinput[2],semiinput[2]);
    cons_par(wspace,2,0);
    RooRealVar* t_12=wspace->var("TauSig3_2_ch0");
    cout<<"toy " <<i<<"\t"<<t_12->getVal()<<"\t"<<t_12->getErrorHi()<<"\t low "<<t_12->getErrorLo()<<endl;
    taumean4=t_12->getVal();
    tauerror4=t_12->getError();
    double higherror4=fabs(t_12->getErrorHi());
    double lowerror4=fabs(t_12->getErrorLo());
    double meanlife4=t_12->getVal();
    if(higherror4==0)higherror4=t_12->getError();
    if(lowerror4==0)lowerror4=t_12->getError();
    if(meanlife4 > 1.70){
      pull4=(meanlife4-1.70)/lowerror4;
      
    }
    else if(meanlife4 < 1.70){
      pull4=(meanlife4-1.70)/higherror4;
    }
    tree4->Fill();

    
    //2012 channel 1
    definemass(wspace, Wsig, Wsemi,2,1, 0.41, 5.277, 0.062, 1.918, 1.91);
    definelifetimeerr(wspace,Wsig,Wsemi,2,1, 11.54, 0.009,12.4,0.004);
    lifetimebkg1(wspace,2,1,0.622 );
    lifetime(wspace,Wsemi,2,1);
    yield(wspace,2,1, Bsinput[3],combinput[3],semiinput[3],bdinput[3]);
    TotalPdfdefinition(wspace,2,1);
    fit_chan(wspace,2,1,bdinput[3],semiinput[3]);
    cons_par(wspace,2,1);
    
    // 2011 channel 0
    definemass(wspace, Wsig, Wsemi,3,0, 0.531, 5.277, 0.040, 2.139, 0.991);
    definelifetimeerr(wspace,Wsig,Wsemi,3,0, 18.3, 0.01,12.4,0.004);
    lifetimebkg1(wspace,3,0,0.517 );
    lifetime(wspace,Wsemi,3,0);
    yield(wspace,3,0, Bsinput[0],combinput[0],semiinput[0],bdinput[0]);
    TotalPdfdefinition(wspace,3,0);
    fit_chan(wspace,3,0,bdinput[0],semiinput[0]);
    cons_par(wspace,3,0);
  //2011 channel1
    
    definemass(wspace, Wsig, Wsemi,3,1, 0.41, 5.277, 0.060, 1.919, 3.371);
    definelifetimeerr(wspace,Wsig,Wsemi,3,1, 11.4, 0.007,12.4,0.004);
    lifetimebkg1(wspace,3,1,0.697);
    lifetime(wspace,Wsemi,3,1);
    yield(wspace,3,1, Bsinput[1],combinput[1],semiinput[1],bdinput[1]);
    TotalPdfdefinition(wspace,3,1);
    fit_chan(wspace,3,1,bdinput[1],semiinput[1]);
    cons_par(wspace,3,1);
    //simultaneous fitting all channel
    simul_fit(wspace);
    RooRealVar* t_final=wspace->var("t_wh");
    cout<<"toy " <<i<<"\t"<<t_final->getVal()<<"\t"<<t_final->getErrorHi()<<"\t low "<<t_final->getErrorLo()<<endl;
    taumeansimall=t_final->getVal();
    tauerrorsimall=t_final->getError();
    double higherror13all=fabs(t_final->getErrorHi());
    double lowerror13all=fabs(t_final->getErrorLo());
    double meanlife13all=t_final->getVal();
    if(higherror13all==0)higherror13all=t_final->getError();
    if(lowerror13all==0)lowerror13all=t_final->getError();
    if(meanlife13all > 1.70){
      pullsimall=(meanlife13all-1.70)/lowerror13all;
      
    }
    else if(meanlife13all < 1.70){
      pullsimall=(meanlife13all-1.70)/higherror13all;
    }
    tree5->Fill();
    //*/
    
    free_par(wspace,0,0);
    free_par2(wspace,0,1);
    free_par(wspace,1,0);
    free_par2(wspace,1,1);
    free_par(wspace,2,0);
    free_par(wspace,2,1);
    free_par(wspace,3,0);
    free_par(wspace,3,1);

  sigpdfMC->Close();
  semipdfMC->Close();
  TProcessID::SetObjectCount(ObjectNumber);  
  }
  file->Write("",TObject::kOverwrite);
}

