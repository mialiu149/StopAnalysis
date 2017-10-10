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
};

struct lepconfig {
  unsigned int nlep;
  float ptcut;
  float etacut;
  bool  dosf;
}; 

struct jetconfig {
  unsigned int njet;
  float ptcut;
  float etacut;
  bool  dojec;
  bool  dobtagsf;
};

struct fillextra {
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
};
