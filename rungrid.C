void rungrid(bool upmu = false, double sinsq2theta23 = 1.0, double deltamsq23 = 2.43E-3, TString inputpath, TString outputpath)
{
  //Load compiled libraries from elsewhere
  gSystem->Load("libAtNuOscillation.so");

  //Compile and load main analysis code
  gROOT->ProcessLine(".L"+inputpath+"/ana.C+");
  gSystem->Load("ana_C.so");
  
  cout<<"rungrid - Pincer 2.0"<<endl;
  TChain * chain = new TChain("ntuple","");
  //

  //Which data set are we using (this projects combines two data sets)
  if(upmu){
    chain->Add(inputpath+"ntpsummary.mc.upmu.0000.root/ntuple");
  } else {
    chain->Add(inputpath+"ntpsummary.mc.atmos.0000.root/ntuple");
  }
  //
  ana t(chain); // Load data into the analysis code
  cout<<"rungrid - Pincer 3.0"<<endl;
  // fupmu = upmu;
  //fsinsq2theta23 = sinsq2theta23;
  //fdeltamsq23 = deltamsq23;

  t.SetParams(upmu,sinsq2theta23,deltamsq23,outputpath,true,inputpath); // Set parameters for analyis. 
  cout<<"rungrid - Pincer 4.0"<<endl;
  t.Loop();// Call the analysis code
  
}
	     
