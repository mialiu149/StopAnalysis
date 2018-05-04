#include "looper.h"
#include "TChain.h"
#include "TString.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

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
  bool isDataFromFileName = bool(infile.Contains("run2_data")||infile.Contains("Run201"));
  bool isSignalFromFileName = bool(infile.Contains("SMS"))||bool(infile.Contains("TChiWH"));

  std::string dataperiod="2016C";
  if (infile.Contains("2016B")) dataperiod="2016B";
  if (infile.Contains("2016D")) dataperiod="2016D";
  if (infile.Contains("2016E")) dataperiod="2016E";
  if (infile.Contains("2016F")) dataperiod="2016F";
  if (infile.Contains("2016H")) dataperiod="2016H";
  if (infile.Contains("2016G")) dataperiod="2016G";
  
  if (infile.Contains("2017B")) dataperiod="2017B";
  if (infile.Contains("2017C")) dataperiod="2017C";
  if (infile.Contains("2017D")) dataperiod="2017D";
  if (infile.Contains("2017E")) dataperiod="2017E";
  if (infile.Contains("2017F")) dataperiod="2017F";
  

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
  int  nVtx              = 1;        
  float met             = 50;      
  if (noskim) met=-1;

  lepconfig lepConfig;
  if (noskim) lepConfig.nlep = -1;
  lepConfig.dosf = true;  
  lepConfig.dolepveto = false;  
  jetconfig jetConfig;
  jetConfig.njet = 2;
  if(noskim)  jetConfig.njet =-1000000;
  jetConfig.dojec = true;
  jetConfig.dobtagsf = true;
  int nBJets            = 0;      if(noskim)  nBJets = -1000000; 

  gammaconfig gammaConfig;

  fillextra fillextraConfig; 
  fillextraConfig.filltaus  =  true;
  fillextraConfig.filltracks  =  true;
  fillextraConfig.filltracks  =  true;

  TChain *sample = new TChain("Events");
  sample->Add(infile.Data());
  if (sample->GetEntries() == 0) std::cout << "WARNING: no entries in chain. filename was: " << infile << std::endl;
  // Run Looper
  mylooper->looper(sample,outfileid,max_events,dirpath,dataperiod,sampleConfig,lepConfig,jetConfig,gammaConfig,fillextraConfig,nVtx,met);
  cout<<"finished"<<endl;
  return 0;
}
