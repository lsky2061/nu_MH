#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TLegend.h"
#include <iostream>
#include "TMath.h"
#include "TCanvas.h"
#include "stdio.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TLine.h"

void plotsensitivity_comb(int mctype = 1011, double mupmumscale = 1.0){
  //TFile* f = new TFile(); 
  // TString basestring = "/minos/app/corwin/atnu_plots/";
  TString basestring = "/minos/data/users/corwin/masshierachy/";
  TString filestring ="";
  TString filestringup = "";

  TString eventtype = "";
  
  if(mctype == 1){
    filestring = "atnu_mc_upmu_nuetonumu_reco.root";
    eventtype = "UpMu";
  } else if(mctype == 0) {
    filestring = "atnu_mc_cv.root";
  } else if(mctype == 2) {
    filestring = "atnu_mc_2flv_cv_nuetonumu_MCtruth.root";
    eventtype = "CV2flv";
  } else if(mctype == 3) {
    filestring = "atnu_mc_2flav_upmu.root";
  } else if(mctype == 1111) {
    filestring = "atnu_mc_cv_nuetonumu_reco_dmsq2.43E-3.root";
    eventtype = "CVreco";
  } else if(mctype == 1211){
    filestring = "atnu_mc_upmu_nuetonumu_reco_dmsq2.43E-3.root";
    eventtype = "UpMu";
  } else if(mctype == 1311){
    filestring = "atnu_mc_cv_nuetonumu_reco_dmsq2.43E-3.root";
    filestringup = "atnu_mc_upmu_nuetonumu_reco_dmsq2.43E-3.root";
  } else if(mctype == 5) {
    filestring = "atnu_mc_cv_nuetonumu_MCtruth.root";
    eventtype = "CV";
  } else if(mctype == 6){
    filestring = "atnu_mc_upmu_nuetonumu_MCtruth.root";
    eventtype = "UpMu";
  } else if(mctype == 1201){
    filestring = "atnu_mc_upmu_nuetonumu_reco_sinsq230.92.root";
    eventtype = "UpMu";
  } else if (mctype == 1301){
    filestring = "atnu_mc_cv_nuetonumu_reco_sinsq230.92_dmsq2.43E-3.root";
    filestringup = "atnu_mc_upmu_nuetonumu_reco_sinsq230.92_dmsq2.43E-3.root"; 
  } else if(mctype == 1101){
    filestring = "atnu_mc_cv_nuetonumu_reco_sinsq230.92_dmsq2.43E-3.root";
    eventtype = "CV";
    //wxyz
    //w = 1 for reco 2 for true
    //x = 1 for CV, 2 for UpMu, 3 for Combined
    //y = sinsq23 (0 = downward fluc, 1 = nominal, 2 = upward fluc)
    //z = deltamsq (Same as above)
  } else if(mctype == 1102) {
    filestring = "atnu_mc_cv_nuetonumu_reco_sinsq230.92_dmsq2.65E-3.root";
    eventtype = "CV";
  } else if(mctype == 1202) {
    filestring = "atnu_mc_upmu_nuetonumu_reco_sinsq230.92_dmsq2.65E-3.root";
    eventtype = "UpMu";
  } else if(mctype == 1302) {
    filestring = "atnu_mc_cv_nuetonumu_reco_sinsq230.92_dmsq2.65E-3.root";
    filestringup = "atnu_mc_upmu_nuetonumu_reco_sinsq230.92_dmsq2.65E-3.root";
  } else if(mctype == 1100){
    filestring = "atnu_mc_cv_nuetonumu_reco_sinsq230.92_dmsq2.22E-3.root";
    eventtype = "CV";
  } else if(mctype == 1200){
    filestring = "atnu_mc_upmu_nuetonumu_reco_sinsq230.92_dmsq2.22E-3.root";
    eventtype = "UpMu";
  } else if(mctype == 1300){
    filestring = "atnu_mc_cv_nuetonumu_reco_sinsq230.92_dmsq2.22E-3.root";
    filestringup = "atnu_mc_upmu_nuetonumu_reco_sinsq230.92_dmsq2.22E-3.root";
  } else if(mctype == 1110){
    filestring = "atnu_mc_cv_nuetonumu_reco_dmsq2.22E-3.root";
    eventtype = "CV";
  } else if(mctype == 1210){
    filestring = "atnu_mc_upmu_nuetonumu_reco_dmsq2.22E-3.root";
    eventtype = "CV";
  } else if(mctype == 1310){
    filestring = "atnu_mc_cv_nuetonumu_reco_dmsq2.22E-3.root";
    filestringup = "atnu_mc_upmu_nuetonumu_reco_dmsq2.22E-3.root";
  } else if(mctype == 1112){
    filestring = "atnu_mc_cv_nuetonumu_reco_dmsq2.65E-3.root";
    eventtype = "CV";
  } else if(mctype == 1212){
    filestring = "atnu_mc_upmu_nuetonumu_reco_dmsq2.65E-3.root";
    eventtype = "UpMu";
  } else if(mctype == 1312){
    filestring = "atnu_mc_cv_nuetonumu_reco_dmsq2.65E-3.root";
    filestringup = "atnu_mc_upmu_nuetonumu_reco_dmsq2.65E-3.root";
  } else {
    cout<<"Unrecoginzed MC Type -- Ending"<<endl;
    return;
  }
 

  TString fullstring = basestring + filestring;

  TFile *f = new TFile(fullstring); //Needs to be a pointer in order for Draw() to work.  See http://root.cern.ch/phpBB2/viewtopic.php?t=9977&sid=5a33c04cb02ba67fc0ed6a81970c2134
  TFile *fup;
  
  //For Combined Analysis
  int tmp;
  int typecountmax = 0;
 

  if(mctype >=2000){ //true 
    tmp = mctype - 2000;
  } else if (mctype >= 1000 && mctype < 2000){ //reco
    tmp = mctype - 1000;
  }

  if(tmp >=300){ 
    fup = new TFile(basestring + filestringup);
    eventtype = "CO";
    typecountmax = 1;
  }


  Double_t pi = TMath::Pi();
  
  cout<<"Pincer 1.0"<<endl;  //For Comined analysis, all of the binning must be the same between the two files, so we can just get this info from one of them.
  TVector3* theta13vec = (TVector3*)f->Get("theta13vec"); //bins, min, max
  TVector3* deltavec = (TVector3*)f->Get("deltavec");
  
  int theta13bins = theta13vec->X();
  double theta13min = theta13vec->Y();
  double theta13max = theta13vec->Z();

  int CPdeltabins = deltavec->X();
  double CPdeltamin = deltavec->Y();
  double CPdeltamax = deltavec->Z();
  cout<<"Pincer 2.0"<<endl;
  TH2D *hchi2sensnorm = new TH2D("hchi2sensnorm","Chi2 dist Normal Hierarchy",
				 theta13bins,theta13min,theta13max,
				 CPdeltabins, CPdeltamin,CPdeltamax);

  TH2D *hchi2sensinv = new TH2D("hchi2sensnorm","Chi2 dist Inverse Hierarchy",
				theta13bins,theta13min,theta13max,
				CPdeltabins, CPdeltamin,CPdeltamax);

  TH2D *hCLnorm = new TH2D("hCLnorm","90% CL for Excluding Inverse Hierarchy",
			   theta13bins,theta13min,theta13max,
			   CPdeltabins, CPdeltamin,CPdeltamax);

  TH2D *hCLinv = new TH2D("hCLinv","90% CL for Excluding Normal Hierarchy",
			   theta13bins,theta13min,theta13max,
			   CPdeltabins, CPdeltamin,CPdeltamax);

  TH1D *hCLnorm1D = new TH1D("hCLnorm1D","",
			     theta13bins,theta13min,theta13max);

  TH1D *hCLinv1D = new TH1D("hCLinv1D","",
			     theta13bins,theta13min,theta13max);

  double tempchi2 = 0.0;
  cout<<"Pincer 3.0"<<endl;
  for(int ith = 0; ith<theta13bins; ith++){ 
    for(int jcp = 0; jcp<CPdeltabins; jcp++){
    //int jcp = 0;

    double sinsq2theta13 = ith*(theta13max/theta13bins);
    double CPdelta = jcp*(CPdeltamax/CPdeltabins);

    char theta13str[10] = {"declared!"};
    char CPdeltastr[10] = {"declared!"};

    int ts = sprintf(theta13str,"%.3f",sinsq2theta13);
    TString theta13tstring = theta13str;

    ts = sprintf(CPdeltastr,"%.3f",CPdelta);
    TString CPdeltastring = CPdeltastr;
  
    //TH2D *hpnorm = (TH2D*)f.Get("costhenup_numu_norm_"+theta13tstring+"_"+CPdeltastring);
    //TH2D *hpinvert = (TH2D*)f.Get("costhenup_numu_invert_"+theta13tstring+"_"+CPdeltastring);
    
    //TH2D *hmnorm = (TH2D*)f.Get("costhenum_numu_norm_"+theta13tstring+"_"+CPdeltastring);
    //TH2D *hminvert = (TH2D*)f.Get("costhenum_numu_invert_"+theta13tstring+"_"+CPdeltastring);
    
    //  TH2D *hpnorm = (TH2D*)f->Get("costhenup_mc_norm_"+theta13tstring+"_"+CPdeltastring);
    //TH2D *hpinvert = (TH2D*)f->Get("costhenup_mc_invert_"+theta13tstring+"_"+CPdeltastring);
    
    //TH2D *hmnorm = (TH2D*)f->Get("costhenum_mc_norm_"+theta13tstring+"_"+CPdeltastring);
    //TH2D *hminvert = (TH2D*)f->Get("costhenum_mc_invert_"+theta13tstring+"_"+CPdeltastring);
    TH2D *hpnorm =  new TH2D();
    TH2D *hpinvert = new TH2D();
    
    TH2D *hmnorm = new TH2D();
    TH2D *hminvert = new TH2D();    


    /*
      hnormp[i]->Rebin2D(10,10);
      hnormm[i]->Rebin2D(10,10);
      hinvtp[i]->Rebin2D(10,10);
      hinvtm[i]->Rebin2D(10,10);
    */

    //hpnorm->Rebin2D(5,5);
    // hpinvert->Rebin2D(5,5);
    //hmnorm->Rebin2D(5,5);
    //hminvert->Rebin2D(5,5);
    
    //For a given value of delta and theta_13, integrate over histograms once each for inverted and normal hierarchy.
    
    //Calculate likelihood of one assuming "data" for the opposite scenario.
    
    double chi2 = 0;
    double epsilon = 0.0001;
    

    for(int h = 0; h<=typecountmax; h++){
      if(h==0){ //CV
	hpnorm = (TH2D*)f->Get("costhenup_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hpinvert = (TH2D*)f->Get("costhenup_mc_invert_"+theta13tstring+"_"+CPdeltastring);
	
	hmnorm = (TH2D*)f->Get("costhenum_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hminvert = (TH2D*)f->Get("costhenum_mc_invert_"+theta13tstring+"_"+CPdeltastring);
      } else {
	hpnorm = (TH2D*)fup->Get("costhenup_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hpinvert = (TH2D*)fup->Get("costhenup_mc_invert_"+theta13tstring+"_"+CPdeltastring);
	
	hmnorm = (TH2D*)fup->Get("costhenum_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hminvert = (TH2D*)fup->Get("costhenum_mc_invert_"+theta13tstring+"_"+CPdeltastring);
      }
    
    //Expect NH, get IH positive mu
      for (int j = 1; j<=hpnorm->GetNbinsX(); j++){
	for (int k = 1; k<=hpnorm->GetNbinsY(); k++){
	  double Nexp = hpnorm->GetBinContent(j,k);
	  double Nobs = hpinvert->GetBinContent(j,k);
	  // Nexp *= mupmumscale;
	  Nobs *= mupmumscale;  //If our observed is different from our expected ratio
	  if (Nexp>0&&Nobs>0){ 
	    //chi2 += 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	    tempchi2 = 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	  } else { 
	    //chi2 += 2*(Nexp-Nobs);
	    tempchi2 = 2*(Nexp-Nobs);
	  }else if(Nexp == 0){
	    tempchi2=0;
	  }
	  chi2 += tempchi2;
	  //  if(sinsq2theta13 == 0.0 && CPdelta==0.0&& tempchi2>1E-5) cout<<h<<"-mu ("<<j<<","<<k<<") ==> + "<<tempchi2<<"-->"<<"chi2 = "<<chi2<<endl;
	  //if(sinsq2theta13 == 0.0 && CPdelta==0.0 && tempchi2<1E-5) cout<<h<<"+mu ("<<j<<","<<k<<") ==> "<<"chi2 = "<<chi2<<endl;
	}
      }
      
      //Expect NH, get IH negative mu
      for (int j = 1; j<=hmnorm->GetNbinsX(); j++){
	for (int k = 1; k<=hmnorm->GetNbinsY(); k++){
	  double Nexp = hmnorm->GetBinContent(j,k);
	  double Nobs = hminvert->GetBinContent(j,k);
	  if (Nexp>0&&Nobs>0){ 
	    //chi2 += 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	    tempchi2 = 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	  } else if(Nobs == 0) { 
	    //chi2 += 2*(Nexp-Nobs);
	    tempchi2 = 2*(Nexp-Nobs);
	  } else if(Nexp == 0){
	    tempchi2=0;
	  }
	  chi2 += tempchi2;
     
	  //  if(sinsq2theta13 == 0.0 && CPdelta==0.0&& tempchi2>1E-5) cout<<h<<"-mu ("<<j<<","<<k<<") ==> + "<<tempchi2<<"-->"<<"chi2 = "<<chi2<<endl;
	}
      }
    }
    //cout<<"chi2 for expecting normal, getting inverted is "<<chi2<<endl;
    //double epsilon = 0.001;
    hchi2sensinv->Fill(sinsq2theta13+epsilon,CPdelta+epsilon,chi2);
    if(chi2>=2.706){ // 90%CL for 1 free parameter
      hCLinv->Fill(sinsq2theta13+epsilon,CPdelta+epsilon);
    }
    if(CPdelta == 0) hCLinv1D->Fill(sinsq2theta13+epsilon,chi2);
    

    chi2 = 0;
    for(int h = 0; h<=typecountmax; h++){
      if(h==0){ //CV
	hpnorm = (TH2D*)f->Get("costhenup_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hpinvert = (TH2D*)f->Get("costhenup_mc_invert_"+theta13tstring+"_"+CPdeltastring);
	
	hmnorm = (TH2D*)f->Get("costhenum_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hminvert = (TH2D*)f->Get("costhenum_mc_invert_"+theta13tstring+"_"+CPdeltastring);
      } else {
	hpnorm = (TH2D*)fup->Get("costhenup_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hpinvert = (TH2D*)fup->Get("costhenup_mc_invert_"+theta13tstring+"_"+CPdeltastring);
	
	hmnorm = (TH2D*)fup->Get("costhenum_mc_norm_"+theta13tstring+"_"+CPdeltastring);
	hminvert = (TH2D*)fup->Get("costhenum_mc_invert_"+theta13tstring+"_"+CPdeltastring);
      }
      
      //Expect IH, get NH positive mu
      for (int j = 1; j<=hpinvert->GetNbinsX(); j++){
	for (int k = 1; k<=hpinvert->GetNbinsY(); k++){
	  double Nexp = hpinvert->GetBinContent(j,k);
	  double Nobs = hpnorm->GetBinContent(j,k);
	  //Nexp *= mupmumscale;
	  Nobs *= mupmumscale;
	  //if (Nexp>0&&Nobs>0) chi2 += 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	  //else chi2 += 2*(Nexp-Nobs);
	  if (Nexp>0&&Nobs>0){ 
	    chi2 += 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	  } else { 
	    chi2 += 2*(Nexp-Nobs);
	  }
	}
      }
      
      
      //Expect IH, get NH negative mu
      for (int j = 1; j<=hminvert->GetNbinsX(); j++){
	for (int k = 1; k<=hminvert->GetNbinsY(); k++){
	  double Nexp = hminvert->GetBinContent(j,k);
	  double Nobs = hmnorm->GetBinContent(j,k);
	  //	if (Nexp>0&&Nobs>0) chi2 += 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	  //else chi2 += 2*(Nexp-Nobs);
	  if (Nexp>0&&Nobs>0){ 
	    chi2 += 2*(Nexp-Nobs+Nobs*log(Nobs/Nexp));
	  } else { 
	    chi2 += 2*(Nexp-Nobs);
	  }
	}
      }
    }
    //cout<<"chi2 for expecting inverted, getting normal is "<<chi2<<endl;
    
    //double epsilon = 0.001;
    hchi2sensnorm->Fill(sinsq2theta13+epsilon,CPdelta+epsilon,chi2);
    if(chi2>=2.706){ // 90%CL for 1 free parameter
      hCLnorm->Fill(sinsq2theta13+epsilon,CPdelta+epsilon);
    }
    if(CPdelta == 0) hCLnorm1D->Fill(sinsq2theta13+epsilon,chi2);

    hpnorm->Delete();
    hpinvert->Delete();
    hmnorm->Delete();
    hminvert->Delete();

    }
  }
  
   cout<<"Entries "<<hchi2sensnorm->GetEntries()<<endl;
   
   gStyle->SetPalette(1);

   TCanvas* can1 = new TCanvas("can1","can1",600,800);
   can1->Divide(1,2);
   can1->cd(1);

   hchi2sensnorm->Draw("COLZ"); 

   can1->cd(2);
   hchi2sensinv->Draw("COLZ");


   TCanvas* canCL = new TCanvas("canCL","90% CL",600,800);
   canCL->Divide(1,2);
   
   canCL->cd(1);
   hCLnorm->Draw("COLZ");

   canCL->cd(2);
   hCLinv->Draw("COLZ");

   TCanvas* can1D = new TCanvas("can1D","chi^2",600,400);
   //can1D->Divide(1,2);

   TLine* CL90line = new TLine(theta13min,2.703,theta13max,2.703);
   
   // can1D->cd(1);
   hCLnorm1D->SetLineColor(2);
   hCLnorm1D->SetLineStyle(2);
   hCLnorm1D->SetXTitle("sin^{2}(2 #theta_{13})");
   hCLnorm1D->SetYTitle("#Delta #chi^{2}");
   hCLnorm1D->SetTitle(filestring);

   hCLnorm1D->Draw("C");

   // can1D->cd(2);
   
   hCLinv1D->SetTitle(filestring);
   hCLinv1D->Draw("C SAME");

   CL90line->SetLineColor(3);
   CL90line->Draw("SAME");

   TLegend* leg = new TLegend(0.7,0.82,1.0,1.0);
   cout<<"Line Color Should be 2, is "<<hCLnorm1D->GetLineColor()<<endl;

   leg->AddEntry("hCLnorm1D","Real NH","L");
   leg->AddEntry("hCLinv1D","Real IN","L");
   leg->AddEntry(CL90line,"90% CL", "L");

   leg->Draw();
  
   
  
}

