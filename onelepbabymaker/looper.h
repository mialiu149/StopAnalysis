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
#include "config.h"
// typedefs
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

#ifndef LOOPER_H
#define LOOPER_H

#pragma GCC diagnostic ignored "-Wwrite-strings"

using namespace std;

class babyMaker {

  public:
    babyMaker ();
    babyMaker (const std::string prefix);
    virtual ~babyMaker (){ delete BabyFile; }

    void MakeBabyNtuple(std::string output_name,fillextra fillextraConfig);
    void InitBabyNtuple();
    int looper(TChain* chain, std::string output_name, int nEvents, std::string path, std::string dataperiod, struct sampleconfig,lepconfig lepConfig,jetconfig jetConfig,  gammaconfig gammaConfig, fillextra fillextraConfig,int skim_nvtx, float skim_met);
    std::string babypath;
    std::string output;
    
    // Variables for baby skim

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
    GenParticleTree gen_leps;
    GenParticleTree gen_nus;
    GenParticleTree gen_tops;
    GenParticleTree gen_qs;
    GenParticleTree gen_bosons;
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
  
  TFile *f_el_SF;
  // Fullsim Electron file
  TFile *f_el_SF_tracking;
  // Fullsim Muon files
  TFile *f_mu_SF_id;
  TFile *f_mu_SF_iso;
  TFile *f_mu_SF_ip;
  TFile *f_mu_SF_tracking;
  TFile *f_mu_SF_veto_id;
  TFile *f_mu_SF_veto_iso;
  TFile *f_mu_SF_veto_ip;
  // Fullsim/Fastsim Electron files
  TFile *f_el_FS_ID;
  TFile *f_el_FS_Iso;
  TFile *f_el_veto_FS_ID;
  TFile *f_el_veto_FS_Iso;

  // Fullsim/Fastsim Muon files
  TFile *f_mu_FS_ID;
  TFile *f_mu_FS_Iso;
  TFile *f_mu_FS_Ip;
  TFile *f_mu_veto_FS_ID;
  TFile *f_mu_veto_FS_Iso;
  TFile *f_mu_veto_FS_Ip;
  TFile *f_vetoLep_eff;

};

//====================//
// Utility Structures //
//====================//
struct val_err_t { float value; float error; };

struct Lepton{
        int id;
        int idx;
        LorentzVector p4;
        //Lepton(id, idx, p4) {id = id; idx = idx; p4 = p4;}
};

struct sortbypt{
  bool operator () (const pair<int, LorentzVector> &v1, const pair<int,LorentzVector> &v2){
    return v1.second.pt() > v2.second.pt();
  }
};


struct sortLepbypt{
  bool operator () (const Lepton &lep1, const Lepton &lep2){
    return lep1.p4.pt() > lep2.p4.pt();
  }
};

struct sortP4byPt {
  bool operator () (const LorentzVector &lv1, const LorentzVector &lv2) { return lv1.pt() > lv2.pt(); }
};
#endif
