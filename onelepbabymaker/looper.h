#include "EventTree.h"
#include "GenParticleTree.h"
#include "LeptonTree.h"
#include "PhotonTree.h"
#include "JetTree.h"
#include "TauTree.h"
#include "IsoTracksTree.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "Math/VectorUtil.h"
#include "TChain.h"
#include "Math/LorentzVector.h"

// typedefs
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

#ifndef LOOPER_H
#define LOOPER_H

#pragma GCC diagnostic ignored "-Wwrite-strings"

using namespace std;

class babyMaker {

  public:
    // Constructor/destructor
    babyMaker ();
    babyMaker (const std::string &prefix);
    virtual ~babyMaker (){ delete BabyFile; }

    void MakeBabyNtuple(std::string output_name="output");
    void InitBabyNtuple();
    int looper(TChain* chain, std::string output_name="output", int nEvents = -1, std::string path = "./");
    std::string babypath;
    std::string output;
    
    // Variables for baby skim
    int skim_nvtx;
    float skim_met;

    int   skim_nGoodLep;
    float skim_goodLep_el_pt;
    float skim_goodLep_el_eta;
    float skim_goodLep_mu_pt;
    float skim_goodLep_mu_eta;
  
    float skim_looseLep_el_pt;
    float skim_looseLep_el_eta;
    float skim_looseLep_mu_pt;
    float skim_looseLep_mu_eta;
    
    float skim_vetoLep_el_pt;
    float skim_vetoLep_el_eta;
    float skim_vetoLep_mu_pt;
    float skim_vetoLep_mu_eta;
  
    int   skim_nJets;
    float skim_jet_pt;
    float skim_jet_eta;

    float skim_jet_ak8_pt;
    float skim_jet_ak8_eta;

    int   skim_nPhotons;
    float skim_ph_pt;
    float skim_ph_eta;

    bool applyJECfromFile;
    int JES_type;
    bool skim_isDataFromFileName;
    bool skim_isSignalFromFileName;
    int  skim_nBJets;
    bool skim_2ndlepveto;

    bool skim_applyBtagSFs;
    bool skim_applyLeptonSFs; 
    bool skim_applyVetoLeptonSFs; 
    bool skim_isFastsim;

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

    void setSkimVariables(bool isDataFromFileName, bool isSignalFromFileName,int nvtx, float met, int nGoodLep, float goodLep_el_pt, float goodLep_el_eta, float goodLep_mu_pt, float goodLep_mu_eta, float looseLep_el_pt, float looseLep_el_eta, float looseLep_mu_pt, float looseLep_mu_eta, float vetoLep_el_pt, float vetoLep_el_eta, float vetoLep_mu_pt, float vetoLep_mu_eta, bool apply2ndlepveto, int njets, float jet_pt, float jet_eta, float jet_ak8_pt, float jet_ak8_eta, int nbjets, int nphs, float phs_pt, float phs_eta, bool applyJEC, int JES_type_central_up_down, bool applyLeptonSFs, bool applyVetoLeptonSFs, bool applyBtagSFs, bool isFastsim,bool filltaus_, bool filltracks_, bool fillZll_, bool fillPhoton_,bool fillMETfilt_, bool fill2ndlep_, bool fillExtraEvtVar_, bool fillAK4EF_, bool fillAK4_Other_, bool fillOverleps_, bool fillAK4Synch_, bool fillElID_, bool fillIso_, bool fillLepSynch_);


  protected:
    TFile* BabyFile;
    TFile* histFile;
    TTree* BabyTree;
    TH1D*  histcounter;
  private:

    // Tree Branches
    EventTree StopEvt;
    LeptonTree lep1;
    LeptonTree lep2;
    PhotonTree ph;
    JetTree jets;
    JetTree jets_jup;
    JetTree jets_jdown;
    TauTree Taus;
    IsoTracksTree Tracks;  
    //GenParticleTree gen_els;
    //GenParticleTree gen_mus;
    //GenParticleTree gen_taus;
    GenParticleTree gen_leps;
    GenParticleTree gen_nus;
    //GenParticleTree gen_nuels;
    //GenParticleTree gen_numus;
    //GenParticleTree gen_nutaus;
    GenParticleTree gen_tops;
    //GenParticleTree gen_bs;
    //GenParticleTree gen_cs;
    GenParticleTree gen_qs;
    //GenParticleTree gen_glus;
    GenParticleTree gen_bosons;
    //GenParticleTree gen_ws;
    //GenParticleTree gen_zs;
    //GenParticleTree gen_phs;
    //GenParticleTree gen_hs;
    GenParticleTree gen_susy;
    
  // for btag SFs
  TH2D* h_btag_eff_b;
  TH2D* h_btag_eff_c;
  TH2D* h_btag_eff_udsg;
  TH2D* h_btag_eff_b_loose;
  TH2D* h_btag_eff_c_loose;
  TH2D* h_btag_eff_udsg_loose;
  
  TH2D* h_btag_eff_b_fastsim;
  TH2D* h_btag_eff_c_fastsim;
  TH2D* h_btag_eff_udsg_fastsim;
  TH2D* h_btag_eff_b_fastsim_loose;
  TH2D* h_btag_eff_c_fastsim_loose;
  TH2D* h_btag_eff_udsg_fastsim_loose;
  
};

struct val_err_t { float value; float error; };

#endif
