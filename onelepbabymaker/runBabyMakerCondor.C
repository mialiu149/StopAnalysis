#include "looper.h"
#include "TChain.h"
#include "TString.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "config.h"

using namespace std;

int main(int argc, char **argv){
  
  //
  // Input sanitation
  //
  if(argc<2){
    cout<<" runBabyMaker takes five arguments: ./runBabyMakerCondor inputfile nevents outfield samplelist" << endl;
    cout<<" Need to provide at least outputfield infile, optional: nevents=-1"<<endl;
    return 0;
  }

  TString infile(argv[1]);
  int max_events = -1;
  std::string outfileid("output");
  std::string dirpath = "./";
  if (argc > 2) max_events = atoi(argv[2]);
  std::cout << "set max number of events to: " << max_events << std::endl;
  if (argc > 3) outfileid = argv[3];
  if (argc > 4) dirpath = argv[4];

  bool isFastsim = bool(infile.Contains("FSPremix") || infile.Contains("FastAsympt25ns") || infile.Contains("Spring16Fast") || infile.Contains("TChiWH"));
  bool isBadMiniAodV1 = bool(infile.Contains("V07-04-12_miniaodv1_FS"));
  bool isDataFromFileName = bool(infile.Contains("run2_data")||infile.Contains("Run2016"));
  bool isSignalFromFileName = bool(infile.Contains("SMS"))||bool(infile.Contains("TChiWH"));

  std::string dataperiod="2016C";
  if (infile.Contains("2016B")) dataperiod="2016B";
  if (infile.Contains("2016D")) dataperiod="2016D";
  if (infile.Contains("2016E")) dataperiod="2016E";
  if (infile.Contains("2016F")) dataperiod="2016F";
  if (infile.Contains("2016H")) dataperiod="2016H";
  if (infile.Contains("2016G")) dataperiod="2016G";
  cout<<__LINE__<<":"<<dataperiod<<endl;
  if (isDataFromFileName) cout << "running on DATA, based on file name: " << infile<<endl;
  else if(isSignalFromFileName) cout << "running on SIGNAL, based on file name: " << infile<<endl;
  else {
    cout << "running on MC, based on file name: " << infile<<endl;
  } 

  sampleconfig sampleConfig;
  sampleConfig.isdata = isDataFromFileName;
  sampleConfig.isfastsim = isFastsim;
  sampleConfig.issignal = isSignalFromFileName;

  //
  // Initialize looper
  //
  babyMaker *mylooper = new babyMaker();

  //
  // Skim Parameters 
  //
  bool noskim = true;
  int nVtx              = 1;        
  float met             = 50;      if (noskim) met=-1;
  int  nGoodLeptons      = 1;      if (noskim) nGoodLeptons = -1;
  float goodLep_el_pt   = 20.0;    
  float goodLep_el_eta  = 1.4442;
  float goodLep_mu_pt   = 20.0;
  float goodLep_mu_eta  = 2.4;

  float looseLep_el_pt  = 10.0;
  float looseLep_el_eta = 2.4;
  float looseLep_mu_pt  = 10.0;
  float looseLep_mu_eta = 2.4;

  float vetoLep_el_pt   = 5.0;
  float vetoLep_el_eta  = 2.4;
  float vetoLep_mu_pt   = 5.0;
  float vetoLep_mu_eta  = 2.4;

  int   nJets           = 2;      if(noskim)  nJets=-1000000;
  float jet_pt          = 30.0;  
  float jet_eta         = 2.4;

  int nBJets            = 0;      if(noskim)  nBJets = -1000000; 

  bool applyJECfromFile = true;
  int  JES_central_up_down = 0;  //0 cetranl, 1 up, -1 down;
  bool applyBtagSFs = true; 
  bool applyLeptonSFs = true;
  bool applyVetoLeptonSFs = true;
  bool apply2ndlepVeto =  false;

  float jet_ak8_pt      = 100.0;
  float jet_ak8_eta     = 2.4;

  int   nphs            = 0;
  float phs_pt          = 20.0;
  float phs_eta         = 2.4;

  bool filltaus_  =  true;
  bool filltracks_  =  true;
  bool fillZll_  =  false;
  bool fillPhoton_  =  false;
  bool fillMETfilt_  =  false;
  bool fill2ndlep_  =  false;
  bool fillExtraEvtVar_  =  false;

  bool fillAK4EF_  =  false;
  bool fillAK4_Other_  =  false;
  bool fillOverleps_  =  false;
  bool fillAK4Synch_  =  false;
  bool fillElID_  =  false;
  bool fillIso_  =  false;
  bool fillLepSynch_  =  false;
  // Input sanitation
  if( !( goodLep_mu_pt>looseLep_mu_pt && looseLep_mu_pt>vetoLep_mu_pt) ){
    cout << "   Problem with muon pT hierachy for good, loose, and veto pT!" << endl;
    cout << "     Exiting..." << endl;
    return 0;
  }

  if( !(goodLep_el_pt>looseLep_el_pt && looseLep_el_pt>vetoLep_el_pt) ){
    cout << "   Problem with electron pT hierarchy for good, loose, and veto pT!" << endl;
    cout << "     Exiting..." << endl;
    return 0;
  }
  // Set Skim Variables
  mylooper->setSkimVariables(isDataFromFileName,isSignalFromFileName,nVtx, met, nGoodLeptons, goodLep_el_pt,  goodLep_el_eta,  goodLep_mu_pt,  goodLep_mu_eta, looseLep_el_pt, looseLep_el_eta, looseLep_mu_pt, looseLep_mu_eta, vetoLep_el_pt, vetoLep_el_eta, vetoLep_mu_pt, vetoLep_mu_eta, apply2ndlepVeto,nJets, jet_pt, jet_eta, jet_ak8_pt, jet_ak8_eta, nBJets,  nphs, phs_pt, phs_eta, applyJECfromFile, JES_central_up_down, applyLeptonSFs, applyVetoLeptonSFs, applyBtagSFs, isFastsim, filltaus_,  filltracks_,  fillZll_,  fillPhoton_, fillMETfilt_,  fill2ndlep_,  fillExtraEvtVar_,  fillAK4EF_,  fillAK4_Other_,  fillOverleps_,  fillAK4Synch_,  fillElID_,  fillIso_,  fillLepSynch_);
  // Intialize TChain, load samples
  TChain *sample = new TChain("Events");
  sample->Add(infile.Data());
  if (sample->GetEntries() == 0) std::cout << "WARNING: no entries in chain. filename was: " << infile << std::endl;
  // Run Looper
  mylooper->looper(sample,outfileid,max_events,dirpath,dataperiod);
  return 0;
}
