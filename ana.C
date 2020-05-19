#define ana_cxx
#include "ana.h"
#include "Riostream.h"
#include "THists.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

using namespace std;

//double rewei(double enu, double costh, double *sigma1, double *sigma2, TH2D *h1, TH2D *h2);


double flux[4][40][20];//[4] = type of simulaiton: nue,nuebar,numu,numubar;
//40 energy bins
//20 costh bins


double rewei(double enu, double costh, int inu);

double func(double *x, double *par){

  //Function that calculates oscillation probabilities for neutrinos traveling throuhg the layers of the Earth's interior 
  
  //const double pi = TMath::Pi();

  //x[0]: neutrino energy (GeV)
  //x[1]: cos(nadir angle)
  
  //par[0]: theta13
  //par[1]: theta23
  //par[2]: theta12
  //par[3]: dmsq32
  //par[4]: dmsq21
  //par[5]: inu0
  //par[6]: inu1
  //par[7]: inmatter

  if (x[1]<=0) return 1;

  double Re = 6371;
  double Rc = 3485.7;

  double rhom = 4.5;
  double rhoc = 11.5;
  
  if (par[5]<0&&par[6]<0){
    rhom = -1*rhom;
    rhoc = -1*rhoc;
  }
  
  int inu0 = int(TMath::Abs(par[5])+0.5);
  int inu1 = int(TMath::Abs(par[6])+0.5);

  rhom = par[7]*rhom;
  rhoc = par[7]*rhoc;

  double Yem = 0.49;
  double Yec = 0.467;

  double Lm = 0;
  double Lc = 0;
  
  if (1-x[1]*x[1]<(Rc*Rc)/(Re*Re)) {

    Lm = Re*(x[1]-sqrt(Rc*Rc/(Re*Re)-(1-x[1]*x[1])));
    Lc = 2*Re*sqrt(Rc*Rc/(Re*Re)-(1-x[1]*x[1]));
  }
  else{
    Lm = Re*x[1];
  }
  //cout<<Lm<<" "<<Lc<<endl;
  //for (int i = 0; i<8; i++) cout<<i<<" "<<par[i]<<endl;

  

  double xm = 1.52648e-4*rhom*Yem*x[0]/par[3];
  double xc = 1.52648e-4*rhoc*Yec*x[0]/par[3];
  //cout<<par[7]<<" "<<rhom<<" "<<Yem<<" "<<x[0]<<" "<<par[3]<<endl;

  double tempm = sqrt(pow(cos(2*par[0])-xm,2)+pow(sin(2*par[0]),2));
  double tempc = sqrt(pow(cos(2*par[0])-xc,2)+pow(sin(2*par[0]),2));
  
  double thetam = 0.5*acos((cos(2*par[0])-xm)/tempm);
  double thetac = 0.5*acos((cos(2*par[0])-xc)/tempc);
  
  //reproduce simona's result
//  if (cos(2*par[0])-xm<0) thetam = 0.5*acos((xm-cos(2*par[0]))/tempm) + pi/2;
//  if (cos(2*par[0])-xc<0) thetac = 0.5*acos((xc-cos(2*par[0]))/tempc) + pi/2;

  double dmsqm = par[3]*tempm;
  double dmsqc = par[3]*tempc;

  double phim = 1.267*dmsqm*Lm/x[0];
  double phic = 1.267*dmsqc*Lc/x[0];

  //cout<<xm<<" "<<xc<<" "<<phim<<" "<<phic<<endl;

  double t1 = cos(2*phim)*cos(phic)-cos(2*(thetac-thetam))*sin(2*phim)*sin(phic);
  double t2 = cos(2*thetam)*(sin(phic)*cos(2*phim)*cos(2*(thetac-thetam))+cos(phic)*sin(2*phim))-sin(phic)*sin(2*thetam)*sin(2*(thetac-thetam));
  double P2 = 1-t1*t1-t2*t2;

  double phi = 1.267*par[3]*(2*Lm+Lc)/x[0];
  phi += 3.8679e-4*rhom*Yem*Lm;
  phi += 0.5*3.8679e-4*rhoc*Yec*Lc;

  double Pemu = pow(sin(par[1]),2)*P2;
  double Petau = pow(cos(par[1]),2)*P2;
  double Pee = 1-P2;
  double Pmutau = pow(sin(par[1]),2)*pow(cos(par[1]),2)*(2-P2-2*(t1*cos(phi)-t2*sin(phi)));
  //double Pmutau = pow(sin(2*par[1]),2)*pow(cos(par[0]),4)*pow(sin(1.267*par[3]*(2*Lm+Lc)/x[0]),2);
  //cout<<Pemu<<" "<<Petau<<" "<<Pee<<" "<<Pmutau<<endl;
  //cout<<pow(sin(par[1]),2)<<" "<<pow(cos(par[1]),2)<<" "<<P2<<" "<<t1<<" "<<t2<<" "<<phi<<endl;

  //Return the probabiliy of oscillation Pee = electron neutrino to electron neurtino, etc.
  if (inu0==12&&inu1==12) return Pee;
  if (inu0==12&&inu1==14) return Pemu;
  if (inu0==12&&inu1==16) return Petau;
  if (inu0==14&&inu1==12) return Pemu;
  if (inu0==14&&inu1==16) return Pmutau;
  if (inu0==14&&inu1==14) return 1-Pemu-Pmutau;
  return 0;
}
  
void ana::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ana.C
//      Root > ana t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //Amount of data in kiloton-years
   double exp_data = 24.6; //updated according to MINOS-doc-6302-v3

   //Amount of simulation in kiloton-years
   double exp_atm  = 7132;
   double exp_cos  = 170./365*5.4;

   double R = 6378;//radius of Earth in im
   double l = 20;

 
   TFile f(Form("atnu_mc/atnu_mc_%.2f.root",th13),"recreate"); //Simulation output file
 
   const int ntot = 6;
   char tag[ntot][100] = {"data","mc","cosmic","numu","nue","nc"}; //Whether the file is data or a type of simulation
   const int nosc = 5;
   char osci[nosc][100] = {"noosc","norm","invert","emu","emurew"}; //Which oscillation model is used (we are comparing "norm" and "invert" to the data.

   //Declare arrays of histograms for storing processed data and simulation
   TH1D *enu[ntot][nosc];
   TH1D *logle[ntot][nosc];
   TH1D *costh[ntot][nosc];
   TH1D *qpeqp[ntot][nosc];
   TH2D *costhenup[ntot][nosc];
   TH2D *costhenum[ntot][nosc];

   TH1D *ntracks[ntot][nosc];
   TH1D *nshowers[ntot][nosc];
   TH1D *vtxx[ntot][nosc];
   TH1D *vtxy[ntot][nosc];
   TH1D *vtxz[ntot][nosc];
   TH1D *trkdirx[ntot][nosc];
   TH1D *trkdiry[ntot][nosc];
   TH1D *trkdirz[ntot][nosc];
   TH1D *trkemcharge[ntot][nosc];
   TH1D *trkplanes[ntot][nosc];
   TH1D *shwplanes[ntot][nosc];
   TH1D *trackextension[ntot][nosc];


   THists *thetaep[2][3][8][2]; //tag, osci, k, l; Angle vs. Energy for mu+
   THists *thetaem[2][3][8][2]; //Same as above for mu-

   //For the following, h = histogram,
   //  enu = reconstructed neutrino energy,
   //  costh = reconstructed angle of neutrino with zenith, which tells us how much earth the neutrino has traversed,
   //  nue = electron neutrinos
   // numu = muon neutrinos
   //  rew = reweighted 

   TH1D *henu_nue = new TH1D("henu_nue","henu_nue",100,0,10);
   TH1D *henu_numu = new TH1D("henu_numu","henu_numu",100,0,10);
   TH1D *henu_nue_rew = new TH1D("henu_nue_rew","henu_nue_rew",100,0,10);

   TH1D *hcosth_nue = new TH1D("hcosth_nue","hcosth_nue",100,-1,1);
   TH1D *hcosth_numu = new TH1D("hcosth_numu","hcosth_numu",100,-1,1);
   TH1D *hcosth_nue_rew = new TH1D("hcosth_nue_rew","hcosth_nue_rew",100,-1,1);

   //set up histograms to properly propagate uncertainty
   henu_nue->Sumw2();
   henu_numu->Sumw2();
   henu_nue_rew->Sumw2();
   hcosth_nue->Sumw2();
   hcosth_numu->Sumw2();
   hcosth_nue_rew->Sumw2();

   char title[1000];
   //Set titles for the arrays of histograms for different data, simulation, and oscillation models
   for (int i = 0; i<ntot; i++){
     for (int j = 0 ; j<nosc; j++){
       sprintf(title,"enu_%s_%s",tag[i],osci[j]);
       enu[i][j] = new TH1D(title,title,12,-1,2);
       enu[i][j]->Sumw2();
       sprintf(title,"logle_%s_%s",tag[i],osci[j]);
       logle[i][j] = new TH1D(title,title,24,-1,5);
       logle[i][j]->Sumw2();
       sprintf(title,"costh_%s_%s",tag[i],osci[j]);
       costh[i][j] = new TH1D(title,title,12,-1,1);
       costh[i][j]->Sumw2();
       sprintf(title,"qpeqp_%s_%s",tag[i],osci[j]);
       qpeqp[i][j] = new TH1D(title,title,20,-10,10);
       qpeqp[i][j]->Sumw2();
       sprintf(title,"costhenup_%s_%s",tag[i],osci[j]);
       costhenup[i][j] = new TH2D(title,title,120,-1,1,90,-1,2);
       costhenup[i][j]->Sumw2();       
       sprintf(title,"costhenum_%s_%s",tag[i],osci[j]);
       costhenum[i][j] = new TH2D(title,title,120,-1,1,90,-1,2);
       costhenum[i][j]->Sumw2();     

       sprintf(title,"ntracks_%s_%s",tag[i],osci[j]);
       ntracks[i][j] = new TH1D(title,title,5,0,5);
       ntracks[i][j]->Sumw2();     
       sprintf(title,"nshowers_%s_%s",tag[i],osci[j]);
       nshowers[i][j] = new TH1D(title,title,5,0,5);
       nshowers[i][j]->Sumw2();     
       sprintf(title,"vtxx_%s_%s",tag[i],osci[j]);
       vtxx[i][j] = new TH1D(title,title,20,-5,5);
       vtxx[i][j]->Sumw2();     
       sprintf(title,"vtxy_%s_%s",tag[i],osci[j]);
       vtxy[i][j] = new TH1D(title,title,20,-5,5);
       vtxy[i][j]->Sumw2();     
       sprintf(title,"vtxz_%s_%s",tag[i],osci[j]);
       vtxz[i][j] = new TH1D(title,title,30,0,30);
       vtxz[i][j]->Sumw2();     
       sprintf(title,"trkdirx_%s_%s",tag[i],osci[j]);
       trkdirx[i][j] = new TH1D(title,title,20,-1,1);
       trkdirx[i][j]->Sumw2();     
       sprintf(title,"trkdiry_%s_%s",tag[i],osci[j]);
       trkdiry[i][j] = new TH1D(title,title,20,-1,1);
       trkdiry[i][j]->Sumw2();     
       sprintf(title,"trkdirz_%s_%s",tag[i],osci[j]);
       trkdirz[i][j] = new TH1D(title,title,20,-1,1);
       trkdirz[i][j]->Sumw2();     
       sprintf(title,"trkemcharge_%s_%s",tag[i],osci[j]);
       trkemcharge[i][j] = new TH1D(title,title,3,-1,2);
       trkemcharge[i][j]->Sumw2();     
       sprintf(title,"trkplanes_%s_%s",tag[i],osci[j]);
       trkplanes[i][j] = new TH1D(title,title,30,0,150);
       trkplanes[i][j]->Sumw2();     
       sprintf(title,"shwplanes_%s_%s",tag[i],osci[j]);
       shwplanes[i][j] = new TH1D(title,title,30,0,30);
       shwplanes[i][j]->Sumw2();     
       sprintf(title,"trackextension_%s_%s",tag[i],osci[j]);
       trackextension[i][j] = new TH1D(title,title,30,-20,130);
       trackextension[i][j]->Sumw2();     
     }
   }

   //Name the 2D histograms for Energy vs. Angle with zenith
   for (int i = 0; i<2; i++){
     for (int j = 0; j<3; j++){
       for (int k = 0; k<8; k++){
	 for (int l = 0; l<2; l++){
	   sprintf(title,"thetaep_%s_%s_%d_%d",tag[i],osci[j],k,l);
	   thetaep[i][j][k][l] = new THists(title);
	   sprintf(title,"thetaem_%s_%s_%d_%d",tag[i],osci[j],k,l);
	   thetaem[i][j][k][l] = new THists(title);
	 }
       }
     }
   }
   //Declare variables and read in the data and simulation
   double mc_costheta;
   double mc_mucostheta;
   double mc_shwcostheta;
   double mc_pshw;
   double mc_pmu;
   double reco_shwcostheta;
   double reco_pshw;
   double reco_pmu;
   double reco_enu;
   double reco_costheta;
   double reco_nucostheta;
   double reco_emu;
   double reco_eshw;
   TTree *t1 = new TTree("atnu","atmospheric tree");
   t1->Branch("simflag",&simflag,"simflag/I"); //1:data 4:mc
   t1->Branch("mc_inu",&mc_inu,"mc_inu/I"); //12:nue 14 numu 16 nutau
   t1->Branch("mc_iact",&mc_iact,"mc_iact/I"); //1: CC 0: NC
   t1->Branch("mc_enu",&mc_enu,"mc_enu/D"); 
   t1->Branch("mc_costheta",&mc_costheta,"mc_costheta/D");
   t1->Branch("mc_mucostheta",&mc_mucostheta,"mc_mucostheta/D");
   t1->Branch("mc_shwcostheta",&mc_shwcostheta,"mc_shwcostheta/D");
   t1->Branch("mc_pshw",&mc_pshw,"mc_pshw/D");
   t1->Branch("mc_pmu",&mc_pmu,"mc_pmu/D");
   t1->Branch("mc_emu",&mc_emu,"mc_emu/D");
   t1->Branch("mc_ehad",&mc_ehad,"mc_ehad/D"); 
   t1->Branch("reco_enu",&reco_enu,"reco_enu/D");
   t1->Branch("reco_costheta",&reco_costheta,"reco_costheta/D");
   t1->Branch("reco_nucostheta",&reco_nucostheta,"reco_nucostheta/D");
   t1->Branch("reco_shwcostheta",&reco_shwcostheta,"reco_shwcostheta/D");
   t1->Branch("reco_emu",&reco_emu,"reco_emu/D");
   t1->Branch("reco_pmu",&reco_pmu,"reco_pmu/D");
   t1->Branch("reco_eshw",&reco_eshw,"reco_eshw/D");
   t1->Branch("reco_pshw",&reco_pshw,"reco_pshw/D");
   t1->Branch("evt_gooddirection",&evt_gooddirection,"evt_gooddirection/I");
   t1->Branch("evt_positivecharge",&evt_positivecharge,"evt_positivecharge/I");
   t1->Branch("evt_negativecharge",&evt_negativecharge,"evt_negativecharge/I");


   int nent = fChain->GetEntries();//The real number of entries (events in data or simulation)

   Long64_t nentries = fChain->GetEntriesFast();//A large number that lets us loop over all entries quickly

   Long64_t nbytes = 0, nb = 0;

   //Main data processing loop. Once for each neutrino event
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<5; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000000==0) cout<<jentry<<"/"<<nent<<endl;//Progress check every million entries
      // if (Cut(ientry) < 0) continue;

      double weight = 1.; //1 for data. Weight simulations down to be the same size as the data set 

      
      int itype = -1;

      //We use and combine several different simulations, this detects which kind (or if we have data) and weights it accordingly. Also records which type it is for analysis further on
      // 0 = data
      // 1 is not used
      // 2 = cosmic rays
      // 3 = muon neutrino simulation (mc = Monte Carlo)
      // 4 = electron neutrino simulation
      // 5 = Neutral current (unknown neutrino type)

      if (simflag==1){//data
	itype = 0;
      }
      else {
	if (mc_inu==0){//cosmic mc
	  itype = 2;
	  weight = exp_data/exp_cos;
	}
	else if (mc_iact==0){//nc mc
	  itype = 5;
	  weight = exp_data/exp_atm;
	}
	else if (abs(mc_inu)==14){//numu mc
	  itype = 3;
	  weight = exp_data/exp_atm;
	}
	else if (abs(mc_inu)==12){//nue mc
	  itype = 4;
	  weight = exp_data/exp_atm;
	}
      }

      //calculate the angle relative to zenith (exact if it is MC)
      if (itype==4){ 
	henu_nue->Fill(mc_enu);
	hcosth_nue->Fill(-mc_pnuy/sqrt(mc_pnux*mc_pnux+mc_pnuy*mc_pnuy+mc_pnuz*mc_pnuz));
      }
      if (itype==3){
	henu_numu->Fill(mc_enu);
	henu_nue_rew->Fill(mc_enu,rewei(mc_enu,-mc_pnuy/sqrt(mc_pnux*mc_pnux+mc_pnuy*mc_pnuy+mc_pnuz*mc_pnuz),mc_inu));
	hcosth_numu->Fill(-mc_pnuy/sqrt(mc_pnux*mc_pnux+mc_pnuy*mc_pnuy+mc_pnuz*mc_pnuz));
	hcosth_nue_rew->Fill(-mc_pnuy/sqrt(mc_pnux*mc_pnux+mc_pnuy*mc_pnuy+mc_pnuz*mc_pnuz),rewei(mc_enu,-mc_pnuy/sqrt(mc_pnux*mc_pnux+mc_pnuy*mc_pnuy+mc_pnuz*mc_pnuz),mc_inu));
      }

      bool contained_vertex = false;

      //We are combining two separate data sets. Contained vertex events occur when a neutrino interacts inside the detector. Upward going muon events occur when a neutrino interacts outside of the detector and produces a muon that goes (upward) into the detector.  We treat these two differently. fc = fully contained, pc = partially contained.

      if (index==0 && evt_goodevent==1 && evt_atmosnumu==1 && evt_goodenergy==1 && evt_gooddirection==1 && evt_goodcharge==1 && (evt_fc==1||evt_pc==1))
	contained_vertex = true;

      if (itype == -1) continue;
      if (contained_vertex == false) continue;
      
      double oscwei[5] = {1.,1.,1.,1.,1.};

      //reset variables to zero for each time through the loop

      mc_costheta = 0;
      mc_mucostheta = 0;
      mc_shwcostheta = 0;
      mc_pshw = 0;

      //calculate reconstructed momentum, angle, and energy for data
      reco_enu = fabs(evt_recoemu)+evt_recoeshwdwgt;
      reco_pmu = 0;
      if (fabs(evt_recoemu)>0.10566){
	reco_pmu = sqrt(pow(evt_recoemu,2)-pow(0.10566,2));
      }
      //reco_costheta = -trk_vtxdiry;
      reco_costheta = -evt_trkdiry;
      if (evt_shwreco==0){
	reco_nucostheta = -evt_trkdiry;
      }
      else{
	double b = 2*reco_pmu*(evt_trkdirx*evt_shwdirx+evt_trkdiry*evt_shwdiry+evt_trkdirz*evt_shwdirz);
	double c = pow(reco_pmu,2)-pow(fabs(evt_recoemu)+evt_recoeshwdwgt,2);
	reco_pshw = (sqrt(b*b-4*c)-b)/2;
	//cout<<pshw<<" "<<evt_recoeshwdwgt<<" "<<sqrt(pow(mc_pnux-mc_pmux,2)+pow(mc_pnuy-mc_pmuy,2)+pow(mc_pnuz-mc_pmuz,2))<<" "<<mc_emu<<" "<<evt_recoemu<<endl;
	reco_nucostheta = -(reco_pmu*evt_trkdiry+reco_pshw*evt_shwdiry)/
	  sqrt(pow(reco_pmu*evt_trkdirx+reco_pshw*evt_shwdirx,2)+pow(reco_pmu*evt_trkdiry+reco_pshw*evt_shwdiry,2)+pow(reco_pmu*evt_trkdirz+reco_pshw*evt_shwdirz,2));	
	//reco_nucostheta = -(reco_emu*evt_trkdiry+reco_pshw*evt_shwdiry)/
	//sqrt(pow(reco_emu*evt_trkdirx+reco_pshw*evt_shwdirx,2)+pow(reco_emu*evt_trkdiry+reco_pshw*evt_shwdiry,2)+pow(reco_emu*evt_trkdirz+reco_pshw*evt_shwdirz,2));	
      }
      reco_shwcostheta = -evt_shwdiry;
      //reco_emu = sqrt(evt_recoemu*evt_recoemu+0.10566*0.10566);
      reco_emu = fabs(evt_recoemu);
      reco_eshw = evt_recoeshwdwgt;

      //Some information we only have in the case of simulation, so we make the calculations using that
      if (itype>0){//mc
	double x[2];
	double par[8];
	x[0] = mc_enu;
	if (mc_pnux*mc_pnux+mc_pnuy*mc_pnuy+mc_pnuz*mc_pnuz>0){
	  x[1] = mc_pnuy/sqrt(mc_pnux*mc_pnux+mc_pnuy*mc_pnuy+mc_pnuz*mc_pnuz);
	  mc_costheta = x[1];
	}
	if (mc_pmux*mc_pmux+mc_pmuy*mc_pmuy+mc_pmuz*mc_pmuz>0){
	  mc_mucostheta = mc_pmuy/sqrt(mc_pmux*mc_pmux+mc_pmuy*mc_pmuy+mc_pmuz*mc_pmuz);
	}
	if (mc_phadx*mc_phadx+mc_phady*mc_phady+mc_phadz*mc_phadz>0){
	  mc_shwcostheta = mc_phady/sqrt(mc_phadx*mc_phadx+mc_phady*mc_phady+mc_phadz*mc_phadz);
	}
	else x[1] = 0;
	mc_pshw = sqrt(pow(mc_pnux-mc_pmux,2)+pow(mc_pnuy-mc_pmuy,2)+pow(mc_pnuz-mc_pmuz,2));
	mc_pmu = sqrt(pow(mc_pmux,2)+pow(mc_pmuy,2)+pow(mc_pmuz,2));
	//par[0] = 0.5*asin(sqrt(0.10));
	//par[0] = 0.;
	if (th13){
	  par[0] = 0.5*asin(sqrt(th13));
	}
	else par[0] = 0.;

	//Apply oscillation parameters to neturinos traveling through Earth
	par[1] = asin(sqrt(0.50));
	par[2] = atan(sqrt(0.46));
	par[3] = 2.43e-3;
	par[4] = 8e-5;
	par[5] = mc_inu;
	par[6] = mc_inu;
	par[7] = 1;
	if (mc_inu!=0&&mc_iact==1){
	  oscwei[1] = func(x,par);
	  if (TMath::Abs(mc_inu)==14){//add nue->numu oscillations
	    if (mc_inu==14){
	      par[5] = 12;
	    }
	    else if(mc_inu==-14){
	      par[5] = -12;
	    }
	    oscwei[1] += func(x,par)*rewei(mc_enu,-mc_costheta,mc_inu);
	    oscwei[3] = func(x,par);
	    oscwei[4] = func(x,par)*rewei(mc_enu,-mc_costheta,mc_inu);
	    par[5] = mc_inu;
	  }
	}
	par[3] = -2.43e-3;
	//Apply oscillation parameters to neturinos traveling through Earth
	if (mc_inu!=0&&mc_iact==1){
	  oscwei[2] = func(x,par);
	  if (TMath::Abs(mc_inu)==14){
	    if (mc_inu==14){
	      par[5] = 12;
	    }
	    else if(mc_inu==-14){
	      par[5] = -12;
	    }
	    oscwei[2] += func(x,par)*rewei(mc_enu,-mc_costheta,mc_inu);
	    par[5] = mc_inu;
	  }
	}
      }

      //if (itype==2) cout<<reco_enu<<" "<<weight<<" "<<oscwei[0]<<endl;

      for (int iosc = 0; iosc<5; iosc++){//iosc

	enu[itype][iosc]->Fill(log10(reco_enu),weight*oscwei[iosc]);
	if (itype>0){
	  enu[1][iosc]->Fill(log10(reco_enu),weight*oscwei[iosc]);
	}
	ntracks[itype][iosc]->Fill(evt_ntracks,weight*oscwei[iosc]);
	if (itype>0) ntracks[1][iosc]->Fill(evt_ntracks,weight*oscwei[iosc]);
	nshowers[itype][iosc]->Fill(evt_nshowers,weight*oscwei[iosc]);
	if (itype>0) nshowers[1][iosc]->Fill(evt_nshowers,weight*oscwei[iosc]);
	vtxx[itype][iosc]->Fill(evt_vtxx,weight*oscwei[iosc]);
	if (itype>0) vtxx[1][iosc]->Fill(evt_vtxx,weight*oscwei[iosc]);
	vtxy[itype][iosc]->Fill(evt_vtxy,weight*oscwei[iosc]);
	if (itype>0) vtxy[1][iosc]->Fill(evt_vtxy,weight*oscwei[iosc]);
	vtxz[itype][iosc]->Fill(evt_vtxz,weight*oscwei[iosc]);
	if (itype>0) vtxz[1][iosc]->Fill(evt_vtxz,weight*oscwei[iosc]);
	trackextension[itype][iosc]->Fill(evt_trackextension,weight*oscwei[iosc]);
	if (itype>0) trackextension[1][iosc]->Fill(evt_trackextension,weight*oscwei[iosc]);

	if (evt_trkreco){
	  trkdirx[itype][iosc]->Fill(evt_trkdirx,weight*oscwei[iosc]);
	  if (itype>0) trkdirx[1][iosc]->Fill(evt_trkdirx,weight*oscwei[iosc]);
	  trkdiry[itype][iosc]->Fill(evt_trkdiry,weight*oscwei[iosc]);
	  if (itype>0) trkdiry[1][iosc]->Fill(evt_trkdiry,weight*oscwei[iosc]);
	  trkdirz[itype][iosc]->Fill(evt_trkdirz,weight*oscwei[iosc]);
	  if (itype>0) trkdirz[1][iosc]->Fill(evt_trkdirz,weight*oscwei[iosc]);
	  trkemcharge[itype][iosc]->Fill(evt_trkemcharge,weight*oscwei[iosc]);
	  if (itype>0) trkemcharge[1][iosc]->Fill(evt_trkemcharge,weight*oscwei[iosc]);
	  trkplanes[itype][iosc]->Fill(evt_trkplanes,weight*oscwei[iosc]);
	  if (itype>0) trkplanes[1][iosc]->Fill(evt_trkplanes,weight*oscwei[iosc]);
	}
	if (evt_shwreco){
	  shwplanes[itype][iosc]->Fill(evt_shwplanes,weight*oscwei[iosc]);
	  if (itype>0) shwplanes[1][iosc]->Fill(evt_shwplanes,weight*oscwei[iosc]);
	}
	if (evt_gooddirection==1){
	  
	  costh[itype][iosc]->Fill(-trk_vtxdiry,weight*oscwei[iosc]);
	  if (itype>0){
	    costh[1][iosc]->Fill(-trk_vtxdiry,weight*oscwei[iosc]);
	  }
	  
	  double L = sqrt((R+l)*(R+l)-R*R*(1-trk_vtxdiry*trk_vtxdiry))+R*trk_vtxdiry;
	  
	  logle[itype][iosc]->Fill(log10(L/(reco_enu)),weight*oscwei[iosc]);
	  if (itype>0){
	    logle[1][iosc]->Fill(log10(L/(reco_enu)),weight*oscwei[iosc]);
	  }
	  if (evt_positivecharge==1){
	    costhenup[itype][iosc]->Fill(-trk_vtxdiry,log10(reco_enu),weight*oscwei[iosc]);
	    //cout<<"pos"<<endl;
	    //thetaep[itype][iosc]->Fill(reco_enu,TMath::ACos(-trk_vtxdiry),weight*oscwei[iosc]);
	    if (itype>0){
	      costhenup[1][iosc]->Fill(-trk_vtxdiry,log10(reco_enu),weight*oscwei[iosc]);
	      if (iosc<3){
		for (int is = 0; is<8; is++){
		  for (int ir = 0; ir<2; ir++){
		    if (jentry%(is+1)==ir){
		      thetaep[1][iosc][is][ir]->Fill(reco_enu,TMath::ACos(-trk_vtxdiry),weight*oscwei[iosc]*(is+1));
		    }
		  }
		}
	      }

	    }
	  }
	  if (evt_negativecharge==1){
	    costhenum[itype][iosc]->Fill(-trk_vtxdiry,log10(reco_enu),weight*oscwei[iosc]);
	    //cout<<"neg"<<endl;
	    //thetaem[itype][iosc]->Fill(reco_enu,TMath::ACos(-trk_vtxdiry),weight*oscwei[iosc]);
	    //cout<<"***********"<<endl;
	    if (itype>0){
	      costhenum[1][iosc]->Fill(-trk_vtxdiry,log10(reco_enu),weight*oscwei[iosc]);
	      if (iosc<3){
		for (int is = 0; is<8; is++){
		  for (int ir = 0; ir<2; ir++){
		    if (jentry%(is+1)==ir){
		      thetaem[1][iosc][is][ir]->Fill(reco_enu,TMath::ACos(-trk_vtxdiry),weight*oscwei[iosc]*(is+1));
		    }
		  }
		}
	      }
	    }
	  }
	      
	}
      }
      t1->Fill();
   }
   t1->Write();
   f.Write();
//   for (int i = 0; i<ntot; i++){
//     for (int j = 0 ; j<nosc; j++){
//       if (i==1&&j<=2){
//	 thetaep[i][j]->SaveObj();
//	 thetaem[i][j]->SaveObj();
//	 thetaep2[i][j]->SaveObj();
//	 thetaem2[i][j]->SaveObj();
//	 thetaep3[i][j]->SaveObj();
//	 thetaem3[i][j]->SaveObj();
//	 thetaep4[i][j]->SaveObj();
//	 thetaem4[i][j]->SaveObj();
//	 thetaep5[i][j]->SaveObj();
//	 thetaem5[i][j]->SaveObj();
//	 thetaep6[i][j]->SaveObj();
//	 thetaem6[i][j]->SaveObj();
//	 thetaep7[i][j]->SaveObj();
//	 thetaem7[i][j]->SaveObj();
//	 thetaep8[i][j]->SaveObj();
//	 thetaem8[i][j]->SaveObj();
//       }
//     }
//   }

   for (int i = 1; i<2; i++){
     for (int j = 0; j<3; j++){
       for (int k = 0; k<8; k++){
	 for (int l = 0; l<2; l++){
	   thetaep[i][j][k][l]->SaveObj();
	   thetaem[i][j][k][l]->SaveObj();
	 }
       }
     }
   }
   f.Close();
}

//double rewei(double enu, double costh, double *sigma1, double *sigma2, TH2D *h1, TH2D *h2){
//
//  double wei = 0;
//  
//  if (enu<0||enu>=30) return 0;
//  
//  int bin = int(enu/30*1000);
//  
//  double binc1 = sigma1[bin];
//  double binc2 = sigma2[bin];
//
//  if (binc1>0&&binc2>0) wei = binc2/binc1;
//  else return 0;
//
//  binc1 = h1->GetBinContent(h1->FindBin(enu, costh));
//  binc2 = h2->GetBinContent(h2->FindBin(enu, costh));
//
//  if (binc1>0&&binc2>0) {
//    wei *= binc1/binc2;
//    return wei;
//  }
//  else return 0;
//
//}

void ana::testrew(double enu, double costh, int inu){
  cout<<rewei(enu,costh,inu)<<endl;
}

double rewei(double enu, double costh, int inu){

  static bool first = true;
  static double aenu[40];
  static double acosth[20];

  if (first){//read flux files
    ifstream in;
    for (int i = 0; i<4; i++){
      if (i==0) in.open("barrflux/fmax20_0401z.sou_nue");
      if (i==1) in.open("barrflux/fmax20_0401z.sou_nbe");
      if (i==2) in.open("barrflux/fmax20_0401z.sou_num");
      if (i==3) in.open("barrflux/fmax20_0401z.sou_nbm");
      char tmp[100];
      double x0,x1,x2,x3,x4;
      int count = 0;
      while (1){
	if (count==0){
	  in.getline(tmp,100);
	}
	else {
	  in>>x0>>x1>>x2>>x3>>x4;
	}
	if (!in.good()) break;
	if (count){
	  flux[i][(count-1)%40][(count-1)/40] = x2;
	  if (count<=40) aenu[count-1] = x0;
	  if ((count-1)%40==0) acosth[(count-1)/40] = x1;
	}
	count++;
      }
      in.clear();
      in.close();
    }
    first = false;
  }

//  for (int i = 0; i<40; i++) cout<<aenu[i]<<endl;
//  for (int i = 0; i<20; i++) cout<<acosth[i]<<endl;
  int ienu = -1;
  int icosth = -1;
  double min = 99999;
  for (int i = 0; i<40; i++){
    if (TMath::Abs(aenu[i]-enu)<min){
      ienu = i;
      min = TMath::Abs(aenu[i]-enu);
    }
  }
  min = 99999;
  for (int i = 0; i<20; i++){
    if (TMath::Abs(acosth[i]-costh)<min){
      icosth = i;
      min = TMath::Abs(acosth[i]-costh);
    }
  } 

  double f0 = flux[0][ienu][icosth]; //nue
  double f1 = flux[1][ienu][icosth]; //nuebar
  double f2 = flux[2][ienu][icosth]; //numu
  double f3 = flux[3][ienu][icosth]; //numubar

  double weight = 0;

  if (inu==14&&f2>0) weight = f0/f2;
  if (inu==-14&&f3>0) weight = f1/f3;
  //cout<<ienu<<" "<<icosth<<" "<<f0<<" "<<f1<<" "<<f2<<" "<<f3<<endl;
  return weight;

}
