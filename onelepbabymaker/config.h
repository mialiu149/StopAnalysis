// I don't like this function, can initialize as the property of the class....
// write a struct to separate flags here.
// isdata, issignal, isfastsim,
// lepton cuts, nlep, leppt, eta, applySF
// jet met cuts, njets,jetpt, eta, met, dobtagSF
// fillextrabranch flags.
struct sampleconfig{
  bool isdata;
  bool issignal;
  bool isfastsim;
  sampleconfig(){
   isdata = false;
   issignal = false;
   isfastsim = false;
  }
};

struct lepconfig {
  int nlep;
  float ptcut;
  float etacut_el;
  float etacut_mu;
  float ptcut_loose;
  float etacut_loose;
  float ptcut_veto;
  float etacut_veto;
  bool  dosf;
  bool  dolepveto;
  lepconfig(){
     nlep = 1;
     ptcut = 20;
     etacut_el = 2.4;
     etacut_mu = 2.4;
     ptcut_loose= 10;
     etacut_loose = 2.4;
     ptcut_veto = 5;
     etacut_veto = 2.4;
     dosf = false;
     dolepveto = false;
   }
}; 

struct jetconfig {
  int njet, nfatjet;
  float ptcut;
  float etacut;
  float ptcut_ak8;
  float etacut_ak8;
  bool  dojec;
  bool  dojecunc;
  bool  dobtagsf;
  jetconfig(){
     njet = 2;//2
     ptcut = 20;//30
     etacut = 2.4;
     nfatjet = 1;
     ptcut_ak8 = 200;
     etacut_ak8 = 2.4;
     dojec = false;
     dojecunc = false;
     dobtagsf = false;
  } 
};

struct gammaconfig {
  float ngamma;
  float ptcut;
  float etacut;
  gammaconfig(){
    ngamma= 0;
    ptcut = 20;
    etacut = 2.4;
   }
};

struct fillextra{
  bool filltaus;
  bool filltracks;
  bool fillZll;
  bool fillPhoton;
  bool fillMETfilt;
  bool fill2ndlep;
  bool fillExtraEvtVar;
  bool fillAK4EF;
  bool fillAK4_Other;
  bool fillOverleps;
  bool fillAK4Synch;
  bool fillElID;
  bool fillIso;
  bool fillLepSynch;

  fillextra(){
   filltaus = false;
   filltracks = false;
   fillZll = false;
   fillPhoton = false;
   fillMETfilt = false;
   fill2ndlep = false;
   fillExtraEvtVar = false;
   fillAK4EF = false;
   fillAK4_Other = false;
   fillOverleps = false;
   fillAK4Synch = false;
   fillElID = false;
   fillIso = false;
   fillLepSynch = false;
  }
};
    // Constructor/destructor
