#include <string.h>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <cmath>
#include <utility>

#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TGraphAsymmErrors.h"
#include "CMS3.h"
#include "MCSelections.h"
#include "StopSelections.h"
#include "TriggerSelections.h"
#include "../analysisutils/mt2w.cc"
#include "../analysisutils/mt2w_bisect.cpp"
#include "../analysisutils/chi2.cc"
#include "../analysisutils/JetUtil.cc"
#include "../analysisutils/topness.cc"
#include "../analysisutils/MT2_implementations.cc"
#include "../analysisutils/histTools.cc"
#include "JetCorrector.h"
#include "IsoTrackVeto.h"
#include "PhotonSelections.h"
#include "MuonSelections.h" //93991
  
#include "IsolationTools.h" //93991
  
#include "MetSelections.h"
#include "looper.h"
#include "goodrun.h"
#include "dorky.cc"

bool debug = false;
bool saveHighPtPFcands = true;
typedef ROOT::Math::LorentzVector < ROOT::Math::PxPyPzE4D < float > > LorentzVector;

using namespace std;
using namespace tas;
using namespace duplicate_removal;

//==============//
// object trees //
//==============//

babyMaker::babyMaker() {
  StopEvt = EventTree();
  lep1 = LeptonTree("lep1_");
  lep2 = LeptonTree("lep2_");
  ph = PhotonTree("ph_");
  jets = JetTree();
  jets_jup = JetTree("jup_");
  jets_jdown = JetTree("jdown_");
  Taus = TauTree();
  Tracks = IsoTracksTree();
  gen_leps = GenParticleTree("leps_");
  gen_nus = GenParticleTree("nus_");
  gen_qs = GenParticleTree("qs_");
  gen_bosons = GenParticleTree("bosons_");
  gen_susy = GenParticleTree("susy_");
}

void babyMaker::MakeBabyNtuple(std::string output_name, fillextra fillextraConfig) {
  std::string str = babypath + output + ".root";
  const char * cstr = str.c_str();
  BabyFile = new TFile(cstr, "RECREATE");
  BabyTree = new TTree("t", "Stop2015 Baby Ntuple");
  StopEvt.SetBranches(BabyTree);
  lep1.SetBranches(BabyTree);
  lep2.SetBranches(BabyTree);
  ph.SetBranches(BabyTree);
  jets.SetAK4Branches(BabyTree);
  jets.SetAK8Branches(BabyTree);
  jets_jup.SetAK4Branches(BabyTree);
  jets_jdown.SetAK4Branches(BabyTree);
  gen_leps.SetBranches(BabyTree);
  gen_nus.SetBranches(BabyTree);
  gen_qs.SetBranches(BabyTree);
  gen_bosons.SetBranches(BabyTree);
  gen_susy.SetBranches(BabyTree);

  //optional, now still computes these variables, does it make sense not to do it when it's false, are we going to save some cpu?
  if (fillextraConfig.filltaus) Taus.SetBranches(BabyTree);
  if (fillextraConfig.filltracks) Tracks.SetBranches(BabyTree);
  if (fillextraConfig.fillZll) StopEvt.SetZllBranches(BabyTree);
  if (fillextraConfig.fillPhoton) StopEvt.SetPhotonBranches(BabyTree);
  if (fillextraConfig.fillMETfilt) StopEvt.SetMETFilterBranches(BabyTree);
  if (fillextraConfig.fill2ndlep) StopEvt.SetSecondLepBranches(BabyTree);
  if (fillextraConfig.fillExtraEvtVar) StopEvt.SetExtraVariablesBranches(BabyTree);
  if (fillextraConfig.fillAK4EF) jets.SetAK4Branches_EF(BabyTree);
  if (fillextraConfig.fillAK4_Other) jets.SetAK4Branches_Other(BabyTree);
  if (fillextraConfig.fillOverleps) jets.SetAK4Branches_Overleps(BabyTree);
  if (fillextraConfig.fillAK4Synch) jets.SetAK4Branches_SynchTools(BabyTree);
  if (fillextraConfig.fillElID) lep1.SetBranches_electronID(BabyTree);
  if (fillextraConfig.fillIso) lep1.SetBranches_Iso(BabyTree);
  if (fillextraConfig.fillLepSynch) lep1.SetBranches_SynchTools(BabyTree);
  if (fillextraConfig.fillElID) lep2.SetBranches_electronID(BabyTree);
  if (fillextraConfig.fillIso) lep2.SetBranches_Iso(BabyTree);
  if (fillextraConfig.fillLepSynch) lep2.SetBranches_SynchTools(BabyTree);
}

void babyMaker::InitBabyNtuple() {
    StopEvt.Reset();
    lep1.Reset();
    lep2.Reset();
    jets.Reset();
    jets_jup.Reset();
    jets_jdown.Reset();
    ph.Reset();
    Taus.Reset();
    Tracks.Reset();
    gen_leps.Reset();
    gen_nus.Reset();
    gen_qs.Reset();
    gen_bosons.Reset();
    gen_susy.Reset();
  }
  //================//
int babyMaker::looper(TChain * chain, std::string output_name, int nEvents, std::string path, std::string dataperiod, sampleconfig sampleConfig, lepconfig lepConfig, jetconfig jetConfig, gammaconfig gammaConfig, fillextra fillextraConfig, int skim_nvtx, float skim_met) {
  // Set output file path
  babypath = path;
  output = output_name;
  // Benchmark
  TBenchmark * bmark = new TBenchmark();
  bmark -> Start("benchmark");
  //Set up loop over chain
  unsigned int nEvents_processed = 0;
  unsigned int nEvents_pass_skim_nVtx = 0;
  unsigned int nEvents_pass_skim_met = 0;
  unsigned int nEvents_pass_skim_nGoodLeps = 0;
  unsigned int nEvents_pass_skim_nGoodJets = 0;
  unsigned int nEvents_pass_skim_nBJets = 0;
  unsigned int nEvents_pass_skim_2ndlepVeto = 0;
  unsigned int nEventsToDo = chain -> GetEntries();
  int jet_overlep1_idx = -9999;
  int jet_overlep2_idx = -9999;

  if (nEvents >= 0) nEventsToDo = nEvents;
  TObjArray * listOfFiles = chain -> GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile * currentFile = 0;
  int JES_type = 0;
  // Lepton MC reco efficiency for veto lep IDs
  // Final scale factor histograms for selected leptons
  TH2D * h_el_SF = NULL;
  TH2D * h_mu_SF = NULL;
  TH2D * h_el_SF_tracking = NULL;
  // Final scale factor histograms for veto leptons
  TH2D * h_el_SF_veto = NULL;
  TH2D * h_mu_SF_veto = NULL;
  TH1D * h_mu_SF_tracking = NULL;
  // Final scale factor histograms for fastim/fullsim for selected leptons
  TH2D * h_el_FS = NULL;
  TH2D * h_mu_FS = NULL;
  // Final scale factor histograms for fastim/fullsim for veto leptons
  TH2D * h_el_veto_FS = NULL;
  TH2D * h_mu_veto_FS = NULL;
  // Final scale factor histograms for lost leptons
  TH2D * h_el_vetoLepEff = NULL;
  TH2D * h_mu_vetoLepEff = NULL;
  // Matching requirement for gen/reco leptons
  double matched_dr = 0.1;

  if ((lepConfig.dosf || lepConfig.dosf) && !sampleConfig.isdata) {
    cout << "  Grabbing lepton scale factors " << endl;
    // Get final fullsim, selected el, sfs
    TH2D * h_el_SF_iso = NULL;
    TH2D * h_el_SF_veto_id = NULL;
    TH2D * h_el_SF_veto_iso = NULL;
    getHist(h_el_SF, "$pwd/lepsf/analysis2016_36p46fb/scaleFactors.root", "GsfElectronToCutBasedSpring15M");
    getHist(h_el_SF_iso, "$pwd/lepsf/analysis2016_36p46fb/scaleFactors.root", "MVAVLooseElectronToMini");
    getHist(h_el_SF_tracking, "$pwd/lepsf/analysis2016_36p46fb/egammaEffi.txt_EGM2D.root", "EGamma_SF2D");
    getHist(h_el_SF_veto_id, "$pwd/lepsf/analysis2016_36p46fb/scaleFactors.root", "GsfElectronToCutBasedSpring15V");
    getHist(h_el_SF_veto_iso, "$pwd/lepsf/analysis2016_36p46fb/scaleFactors.root", "MVAVLooseElectronToMini2");
    h_el_SF -> Multiply(h_el_SF_iso);
    h_el_SF = (TH2D * ) h_el_SF -> Clone("h_el_SF"); // fix the name
    h_el_SF_veto = (TH2D * ) h_el_SF_veto_id -> Clone("h_el_SF_veto");
    h_el_SF_veto -> Multiply(h_el_SF_veto_iso);
    // Fastsim/Fullsim el files
    TH2D * h_el_FS_ID = NULL;
    TH2D * h_el_FS_Iso = NULL;
    TH2D * h_el_veto_FS_ID = NULL;
    TH2D * h_el_veto_FS_Iso = NULL;
    getHist(h_el_FS_ID, "$pwd/lepsf/analysis2016_36p46fb/sf_el_mediumCB.root", "histo2D");
    getHist(h_el_FS_Iso, "$pwd/lepsf/analysis2016_36p46fb/sf_el_mini01.root", "histo2D");
    getHist(h_el_veto_FS_ID, "$pwd/lepsf/analysis2016_36p46fb/sf_el_vetoCB.root", "histo2D");
    getHist(h_el_veto_FS_Iso, "$pwd/lepsf/analysis2016_36p46fb/sf_el_mini02.root", "histo2D");
    h_el_FS = (TH2D * ) h_el_FS_ID -> Clone("h_el_FS");
    h_el_FS -> Multiply(h_el_FS_Iso);
    h_el_veto_FS = (TH2D * ) h_el_veto_FS_ID -> Clone("h_el_FS");
    h_el_veto_FS -> Multiply(h_el_veto_FS_Iso);

    getHist(h_el_vetoLepEff, "lepsf/analysis2016_36p46fb/lepeff__moriond17__ttbar_powheg_pythia8_25ns.root", "h2_lepEff_vetoSel_Eff_el");
    getHist(h_mu_vetoLepEff, "lepsf/analysis2016_36p46fb/lepeff__moriond17__ttbar_powheg_pythia8_25ns.root", "h2_lepEff_vetoSel_Eff_mu");
    // Muon files
    TH2D * h_mu_SF_id = NULL;
    TH2D * h_mu_SF_iso = NULL;
    TH2D * h_mu_SF_ip = NULL;
    TH2D * h_mu_SF_veto_id = NULL;
    TH2D * h_mu_SF_veto_iso = NULL;
    TH2D * h_mu_SF_veto_ip = NULL;
    TH2D * h_mu_FS_ID = NULL;
    TH2D * h_mu_FS_Iso = NULL;
    TH2D * h_mu_FS_Ip = NULL;
    TH2D * h_mu_veto_FS_ID = NULL;
    TH2D * h_mu_veto_FS_Iso = NULL;
    TH2D * h_mu_veto_FS_Ip = NULL;
    f_mu_SF_tracking = new TFile("lepsf/analysis2016_12p9fb/muons_tracking_sf.root", "read");
    getHist(h_mu_SF_id, "lepsf/analysis2016_36p46fb/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root", "SF");
    getHist(h_mu_SF_iso, "lepsf/analysis2016_36p46fb/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root", "SF");
    getHist(h_mu_SF_ip, "lepsf/analysis2016_36p46fb/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root", "SF");
    getHist(h_mu_SF_veto_id, "lepsf/analysis2016_36p46fb/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root", "SF");
    getHist(h_mu_SF_veto_iso, "lepsf/analysis2016_36p46fb/TnP_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root", "SF");
    getHist(h_mu_SF_veto_ip, "lepsf/analysis2016_36p46fb/TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root", "SF");
    getHist(h_mu_FS_ID, "lepsf/analysis2016_36p46fb/sf_mu_mediumID.root", "histo2D");
    getHist(h_mu_FS_Iso, "lepsf/analysis2016_36p46fb/sf_mu_mediumID_mini02.root", "histo2D");
    getHist(h_mu_FS_Ip, "lepsf/analysis2016_36p46fb/sf_mu_mediumID_tightIP2D.root", "histo2D");
    getHist(h_mu_veto_FS_ID, "lepsf/analysis2016_36p46fb/sf_mu_looseID.root", "histo2D");
    getHist(h_mu_veto_FS_Iso, "lepsf/analysis2016_36p46fb/sf_mu_looseID_mini02.root", "histo2D");
    getHist(h_mu_veto_FS_Ip, "lepsf/analysis2016_36p46fb/sf_mu_mediumID_looseIP2D.root", "histo2D");

    // Grab mu histos
    TGraphAsymmErrors * h_mu_SF_tracking_temp = (TGraphAsymmErrors * ) f_mu_SF_tracking -> Get("ratio_eta");
    // Double unc. on selected muon id sfs, since not for our exact wp
    doubleSysError(h_mu_SF_id);
    doubleSysError(h_mu_SF_iso);
    doubleSysError(h_mu_SF_ip);
    h_mu_SF = (TH2D * ) h_mu_SF_id -> Clone("h_mu_SF");
    h_mu_SF -> Multiply(h_mu_SF_iso);
    h_mu_SF -> Multiply(h_mu_SF_ip);

    // Get final fullsim, selected muon, tracking sfs, convert TGraphErrors
    int nX = h_mu_SF_tracking_temp -> GetN();
    Double_t * x_val = h_mu_SF_tracking_temp -> GetX();
    Double_t * y_val = h_mu_SF_tracking_temp -> GetY();
    Double_t * y_err_up = h_mu_SF_tracking_temp -> GetEYhigh();
    Double_t * y_err_low = h_mu_SF_tracking_temp -> GetEYhigh();
    h_mu_SF_tracking = new TH1D("h_mu_SF_tracking", "h_mu_SF_tracking", nX - 1, x_val);
    for (int i = 0; i < nX; i++) {
      h_mu_SF_tracking -> SetBinContent(i + 1, y_val[i]);
      h_mu_SF_tracking -> SetBinError(i + 1, std::max(y_err_up[i], y_err_low[i]));
    }
    // Double unc. on veto muon ip sfs, since not for our exact wp
    doubleSysError(h_mu_SF_veto_ip);
    h_mu_SF_veto = (TH2D * ) h_mu_SF_veto_id -> Clone("h_mu_SF_veto");
    h_mu_SF_veto -> Multiply(h_mu_SF_veto_iso);
    h_mu_SF_veto -> Multiply(h_mu_SF_veto_ip);
    // Double unc. on selected muon FS id sfs, since not for our exact wp
    doubleSysError(h_mu_FS_Iso);
    doubleSysError(h_mu_FS_Ip);
    h_mu_FS = (TH2D * ) h_mu_FS_ID -> Clone("h_mu_FS");
    h_mu_FS -> Multiply(h_mu_FS_Iso);
    h_mu_FS -> Multiply(h_mu_FS_Ip);
    // Double unc. on selected muon FS ip sfs, since not for our exact wp
    doubleSysError(h_mu_veto_FS_Ip);
    h_mu_veto_FS = (TH2D * ) h_mu_veto_FS_ID -> Clone("h_mu_veto_FS");
    h_mu_veto_FS -> Multiply(h_mu_veto_FS_Iso);
    h_mu_veto_FS -> Multiply(h_mu_veto_FS_Ip);
  }

  TFile * fxsec;
  TH1D * hxsec;
  if (sampleConfig.issignal) {
    fxsec = new TFile("xsec_susy_13tev.root", "READ");
    if (fxsec -> IsZombie()) {
      std::cout << "Somehow xsec_stop_13TeV.root is corrupted. Exit..." << std::endl;
      exit(0);
    }
    hxsec = (TH1D * ) fxsec -> Get("h_xsec_c1n2");
  }

  TFile * pileupfile;
  TH1D * hPU;
  TH1D * hPUup;
  TH1D * hPUdown;
  if (!sampleConfig.isdata) {
    pileupfile = new TFile("puWeights_2016data_36p6fbinv.root", "READ");
    if (pileupfile -> IsZombie()) {
      std::cout << "Somehow puWeights_2016data_36p6fbinv.root is corrupted. Exit..." << std::endl;
      exit(0);
    }
    hPU = (TH1D * ) pileupfile -> Get("puWeight");
    hPUup = (TH1D * ) pileupfile -> Get("puWeightUp");
    hPUdown = (TH1D * ) pileupfile -> Get("puWeightDown");
  }
  //////////////////////////
  ////counter histograms////
  //////////////////////////
  string labels[] = {
    "nominal,muR=1 muF=1",
    "muR=1 muF=2",
    "muR=1 muF=0.5",
    "muR=1 muF=0.5",
    "muR=2 muF=1",
    "muR=2 muF=2",
    "muR=2 muF=0.5",
    "muR=0.5 muF=1",
    "muR=0.5 muF=2",
    "muR=0.5 muF=0.5",
    "pdf_up",
    "pdf_down",
    "pdf_alphas_var_1",
    "pdf_alphas_var_2",
    "weight_btagsf",
    "weight_btagsf_heavy_UP",
    "weight_btagsf_light_UP",
    "weight_btagsf_heavy_DN",
    "weight_btagsf_light_DN",
    "weight_ISR_nominal",
    "weight_ISR_up",
    "weight_ISR_down",
    "NEvents",
    "weight_btagsf_fastsim_UP",
    "weight_btagsf_fastsim_DN",
    "weight_ISRnjets_nominal",
    "weight_ISRnjets_up",
    "weight_ISRnjets_down",
    "weight_lepSF",
    "weight_lepSF_up",
    "weight_lepSF_down",
    "weight_vetoLepSF",
    "weight_vetoLepSF_up",
    "weight_vetoLepSF_down",
    "weight_lepSF_fastSim",
    "weight_lepSF_fastSim_up",
    "weight_lepSF_fastSim_down"
  };
  TH1D * counterhist = new TH1D("h_counter", "h_counter", 36, 0.5, 36.5);
  setcounterLabel(counterhist, labels);
  TH3D * counterhistSig;
  TH2F * histNEvts; //count #evts per signal point
  if (sampleConfig.issignal) { //create histos only for signals
    counterhistSig = new TH3D("h_counterSMS", "h_counterSMS", 37, 99, 1024, 19, -1, 474, 35, 0.5, 35.5); //15000 bins!
    setcounterLabel3D(counterhistSig, labels);
    histNEvts = new TH2F("histNEvts", "h_histNEvts", 37, 99, 1024, 19, -1, 474); //x=mStop, y=mLSP
    histNEvts -> Sumw2();
  }
  // output Ntuple  and branches
  MakeBabyNtuple(output_name + ".root", fillextraConfig);
  InitBabyNtuple();
  // Set JSON file
  const char * json_file = "json_files/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);
  // set JEC files
  std::vector < std::string > jetcorr_filenames_pfL1FastJetL2L3;
  FactorizedJetCorrector * jet_corrector_pfL1FastJetL2L3(0);
  JetCorrectionUncertainty * jetcorr_uncertainty(0);
  JetCorrectionUncertainty * jetcorr_uncertainty_sys(0);
  jetcorr_filenames_pfL1FastJetL2L3.clear();
  std::string jetcorr_uncertainty_filename;
  char * jecpath;
  jecpath = getenv("TOOLSPATH");
  // files for RunIISpring15 MC
  // From Nick Amin:
  //I put the latest JECs into CORE (2017 branch) [1]https://github.com/cmstas/CORE/tree/2017/Tools/jetcorr/data/run2_25ns
  //Fall17_*_V6_{DATA,MC}
  if (sampleConfig.isdata) {
    if (dataperiod.find("2016B") != std::string::npos ||
      dataperiod.find("2016C") != std::string::npos ||
      dataperiod.find("2016D") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
      jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }

    if (dataperiod.find("2016E") != std::string::npos ||
      dataperiod.find("2016F") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
      jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }

    if (dataperiod.find("2016G") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
      jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }

    if (dataperiod.find("2016H") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
      jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }
    //Adding New...
    
    if (dataperiod.find("2017B") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
                   jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }
    if (dataperiod.find("2017C") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
                   jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }
    if (dataperiod.find("2017D") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
                   jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }
    if (dataperiod.find("2017E") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
                   jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }
    if (dataperiod.find("2017F") != std::string::npos) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L1FastJet_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L2Relative_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L3Absolute_AK4PFchs.txt", jecpath));
      jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L2L3Residual_AK4PFchs.txt", jecpath));
                   jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_Uncertainty_AK4PFchs.txt", jecpath);
    }
    
  
//  } else if (sampleConfig.issignal) {
//    jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Spring16_FastSimV1_L1FastJet_AK4PFchs.txt", jecpath));
//    jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Spring16_FastSimV1_L2Relative_AK4PFchs.txt", jecpath));
//    jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Spring16_FastSimV1_L3Absolute_AK4PFchs.txt", jecpath));
//                 jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Spring16_FastSimV1_Uncertainty_AK4PFchs.txt", jecpath);
  } else {
    jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK4PFchs.txt", jecpath));
    jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK4PFchs.txt", jecpath));
    jetcorr_filenames_pfL1FastJetL2L3.push_back(Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK4PFchs.txt", jecpath));
                 jetcorr_uncertainty_filename = Form("%s/jetcorr/data/run2_25ns/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt", jecpath);
  }

  cout << "applying JEC from following files:" << endl;
  for (unsigned int ifile = 0; ifile < jetcorr_filenames_pfL1FastJetL2L3.size(); ++ifile) {
    cout << "   " << jetcorr_filenames_pfL1FastJetL2L3.at(ifile) << endl;
  }

  jet_corrector_pfL1FastJetL2L3 = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  if (!sampleConfig.isdata && jetConfig.dojecunc != 0) {
    cout << "applying JEC uncertainties with weight " << jetConfig.dojecunc << " from file: " << endl << "   " << jetcorr_uncertainty_filename << endl;
    jetcorr_uncertainty = new JetCorrectionUncertainty(jetcorr_uncertainty_filename);
  }
  jetcorr_uncertainty_sys = new JetCorrectionUncertainty(jetcorr_uncertainty_filename);

  if (jetConfig.dobtagsf) {
    jets.InitBtagSFTool(sampleConfig.isfastsim);
    jets_jup.InitBtagSFTool(sampleConfig.isfastsim);
    jets_jdown.InitBtagSFTool(sampleConfig.isfastsim);
  }
  // Loop over the trees

  while ((currentFile = (TFile * ) fileIter.Next())) {
    // Get File Content
    TFile file(currentFile -> GetTitle());
    TTree * tree = (TTree * ) file.Get("Events");
    cms3.Init(tree);
    TString thisfilename = file.GetName();
    cout << "file name is " << file.GetName() << endl;

    bool isbadrawMET = false;
    if (thisfilename.Contains("V07-04-12_miniaodv1_FS")) {
      cout << "This file seems to have a badrawMET, thus MET needs to be recalculated" << endl;
      isbadrawMET = true;
    }
    //
    // Loop over Events 
    //
    unsigned int nEventsTree = tree -> GetEntriesFast();
    cout << "Total number of events:" << nEventsTree << endl;

    for (unsigned int evt = 0; evt < nEventsTree; evt++) {
      if (nEvents_processed >= nEventsToDo) break;
      cms3.GetEntry(evt); // Get Event Content
      nEvents_processed++;
      CMS3::progress(nEvents_processed, nEventsToDo); // Progress bar
      InitBabyNtuple(); // Intialize Baby NTuple Branches
      // calculate sum of weights and save them in a hisogram.
      float pdf_weight_up = 1;
      float pdf_weight_down = 1;
      float sum_of_weights = 0;
      float average_of_weights = 0;
      counterhist -> Fill(22, 1.);
      if (!evt_isRealData()) {
        //error on pdf replicas
        if (genweights().size() > 109) {
          for (int ipdf = 9; ipdf < 109; ipdf++) {
            average_of_weights += cms3.genweights().at(ipdf);
          } // average of weights
          average_of_weights = average_of_weights / 100;
          for (int ipdf = 9; ipdf < 109; ipdf++) {
            sum_of_weights += (cms3.genweights().at(ipdf) - average_of_weights) * (cms3.genweights().at(ipdf) - average_of_weights);
          } //std of weights.     
          pdf_weight_up = (average_of_weights + sqrt(sum_of_weights / 99));
          pdf_weight_down = (average_of_weights - sqrt(sum_of_weights / 99));
          StopEvt.pdf_up_weight = pdf_weight_up;
          StopEvt.pdf_down_weight = pdf_weight_down;
          counterhist -> Fill(1, genweights()[0]);
          counterhist -> Fill(2, genweights()[1]);
          counterhist -> Fill(3, genweights()[2]);
          counterhist -> Fill(4, genweights()[3]);
          counterhist -> Fill(5, genweights()[4]);
          counterhist -> Fill(6, genweights()[5]);
          counterhist -> Fill(7, genweights()[6]);
          counterhist -> Fill(8, genweights()[7]);
          counterhist -> Fill(9, genweights()[8]);
          counterhist -> Fill(10, pdf_weight_up);
          counterhist -> Fill(11, pdf_weight_down);
          counterhist -> Fill(12, genweights()[109]); // α_s variation. 
          counterhist -> Fill(13, genweights()[110]); // α_s variation. 
        }
      }
      //////////////////////////////////////////
      //              good run list           //
      //////////////////////////////////////////
      //  if( evt_isRealData() && !goodrun(evt_run(), evt_lumiBlock()) ) continue;
      if (evt_isRealData()) {
        DorkyEventIdentifier id(evt_run(), evt_event(), evt_lumiBlock());
        if (is_duplicate(id)) continue;
      }

      if (debug) std::cout << "[babymaker::looper]: filling event vars LINE:" << __LINE__ << std::endl;
      StopEvt.FillCommon(file.GetName()); // Fill Event Variables 
      if (debug) std::cout << "[babymaker::looper]: filling event vars completed LINE:" << __LINE__ << std::endl;

      if (!StopEvt.is_data) {
        StopEvt.weight_PU = hPU -> GetBinContent(hPU -> FindBin(StopEvt.pu_ntrue));
        StopEvt.weight_PUup = hPUup -> GetBinContent(hPUup -> FindBin(StopEvt.pu_ntrue));
        StopEvt.weight_PUdown = hPUdown -> GetBinContent(hPUdown -> FindBin(StopEvt.pu_ntrue));
      }

      if (sampleConfig.issignal) {
        //get susy particle masses from sparms
        //for (unsigned int nsparm = 0; nsparm < sparm_names().size(); ++nsparm) {
        //  if (sparm_names().at(nsparm).Contains("mCh")) StopEvt.mass_chargino = sparm_values().at(nsparm);
        //  StopEvt.mass_stop = StopEvt.mass_chargino; // modification for whmet
        //  if (sparm_names().at(nsparm).Contains("mLSP")) StopEvt.mass_lsp = sparm_values().at(nsparm);
        //  if (sparm_names().at(nsparm).Contains("mGl")) StopEvt.mass_gluino = sparm_values().at(nsparm);
        //}
        if (genps_weight() > 0) histNEvts -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 1);
        else if (genps_weight() < 0) histNEvts -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, -1);
        StopEvt.xsec = hxsec -> GetBinContent(hxsec -> FindBin(StopEvt.mass_stop));
        StopEvt.xsec_uncert = hxsec -> GetBinError(hxsec -> FindBin(StopEvt.mass_stop));

        float SMSpdf_weight_up = 1;
        float SMSpdf_weight_down = 1;
        float SMSsum_of_weights = 0;
        float SMSaverage_of_weights = 0;
        //error on pdf replicas
        //fastsim has first genweights bin being ==1
        if (genweights().size() > 111) { //fix segfault
          for (int ipdf = 10; ipdf < 110; ipdf++) {
            SMSaverage_of_weights += cms3.genweights().at(ipdf);
          } // average of weights
          SMSaverage_of_weights = average_of_weights / 100.;
          for (int ipdf = 10; ipdf < 110; ipdf++) {
            SMSsum_of_weights += pow(cms3.genweights().at(ipdf) - SMSaverage_of_weights, 2);
          } //std of weights.
          SMSpdf_weight_up = (average_of_weights + sqrt(SMSsum_of_weights / 99.));
          SMSpdf_weight_down = (average_of_weights - sqrt(SMSsum_of_weights / 99.));
          StopEvt.pdf_up_weight = SMSpdf_weight_up; //overwrite here, although it should not matter
          StopEvt.pdf_down_weight = SMSpdf_weight_down; //overwrite here, although it should not matter
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 1, genweights()[1]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 2, genweights()[2]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 3, genweights()[3]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 4, genweights()[4]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 5, genweights()[5]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 6, genweights()[6]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 7, genweights()[7]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 8, genweights()[8]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 9, genweights()[9]);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 10, SMSpdf_weight_up);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 11, SMSpdf_weight_down);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 12, genweights()[110]); // α_s variation. 
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 13, genweights()[111]); // α_s variation.
        }
      } // is signal
      //
      // Gen Information
      //
      if (debug) std::cout << "[babymaker::looper]: filling gen particles vars  LINE:" << __LINE__ << std::endl;

      //ttbar counters using neutrinos:
      int n_nutaufromt(0), n_nuelfromt(0), n_numufromt(0), nLepsHardProcess(0), nNusHardProcess(0);
      bool ee0lep(false), ee1lep(false), ge2lep(false), zToNuNu(false);
      bool isstopevent(false), istopevent(false), isWevent(false), isZevent(false);

      TString thisFile = chain -> GetFile() -> GetName();
      vector < LorentzVector > genpnotISR;
      genpnotISR.clear();

      if (thisFile.Contains("DYJets") || thisFile.Contains("ZJets") || thisFile.Contains("ZZ") || thisFile.Contains("WZ") || thisFile.Contains("TTZ") || thisFile.Contains("tZq"))
        isZevent = true;
      if (thisFile.Contains("WJets") || thisFile.Contains("WW") || thisFile.Contains("WZ") || thisFile.Contains("TTW")) isWevent = true;
      if (thisFile.Contains("TTJets") || thisFile.Contains("TTTo2L2Nu") || thisFile.Contains("TT_") || thisFile.Contains("ST_") || thisFile.Contains("tZq") ||
        thisFile.Contains("TTW") || thisFile.Contains("TTZ")) istopevent = true;
      if (sampleConfig.issignal || thisFile.Contains("T2tt") || thisFile.Contains("T2tb") || thisFile.Contains("T2bW") || thisFile.Contains("TChWH"))
        isstopevent = true;
      //gen particles, put all these in gen particle tree. sooooo messy.
      if (!evt_isRealData()) {
        for (unsigned int genx = 0; genx < genps_p4().size(); genx++) {
          if (genps_id().at(genx) == genps_id_simplemother().at(genx) && !genps_isLastCopy().at(genx)) continue;
          if (genps_isHardProcess().at(genx) ||
            genps_fromHardProcessDecayed().at(genx) ||
            genps_fromHardProcessFinalState().at(genx) ||
            genps_isLastCopy().at(genx) ||
            genps_status().at(genx) == 1) {
            if (abs(genps_id().at(genx)) == pdg_el || abs(genps_id().at(genx)) == pdg_mu || abs(genps_id().at(genx)) == pdg_tau || abs(genps_id().at(genx)) == pdg_b || abs(genps_id().at(genx)) == pdg_c || abs(genps_id().at(genx)) == pdg_s || abs(genps_id().at(genx)) == pdg_d || abs(genps_id().at(genx)) == pdg_u) { //if have lepton or quark --> those could be matched to a jet
              int mother_idx = genps_idx_mother().at(genx);
              int mother_id = -999;
              if (mother_idx >= 0) mother_id = abs(genps_id().at(mother_idx));
              if (mother_id == pdg_t || mother_id == pdg_W || mother_id == 1000006 || mother_id == 2000006 || mother_id == pdg_chi_1plus1 || mother_id == pdg_chi_2plus1) {
                genpnotISR.push_back(genps_p4().at(genx));
              }
            }
            if (abs(genps_id().at(genx)) == pdg_el || abs(genps_id().at(genx)) == pdg_mu || abs(genps_id().at(genx)) == pdg_tau) gen_leps.FillCommon(genx);
            if (abs(genps_id().at(genx)) == pdg_numu || abs(genps_id().at(genx)) == pdg_nue || abs(genps_id().at(genx)) == pdg_nutau) gen_nus.FillCommon(genx);
            if (abs(genps_id().at(genx)) == pdg_t || abs(genps_id().at(genx)) == pdg_b || abs(genps_id().at(genx)) == pdg_c || abs(genps_id().at(genx)) == pdg_s || abs(genps_id().at(genx)) == pdg_d || abs(genps_id().at(genx)) == pdg_u) gen_qs.FillCommon(genx);
            if (abs(genps_id().at(genx)) == pdg_W || abs(genps_id().at(genx)) == pdg_Z || (abs(genps_id().at(genx)) == pdg_ph &&
                genps_p4().at(genx).Pt() > 5.0) || abs(genps_id().at(genx)) == pdg_h) gen_bosons.FillCommon(genx);

            //SUSY particles
            if (abs(genps_id().at(genx)) >= 1000000 && abs(genps_id().at(genx)) <= 1000040) gen_susy.FillCommon(genx);

            if (abs(genps_id_mother().at(genx)) == pdg_W && abs(genps_id().at(genx)) == pdg_nue && genps_status().at(genx) == 1 && abs(genps_id_mother().at(genps_idx_mother().at(genx))) == pdg_t) n_nuelfromt++;

            if (abs(genps_id_mother().at(genx)) == pdg_W && abs(genps_id().at(genx)) == pdg_numu && genps_status().at(genx) == 1 && abs(genps_id_mother().at(genps_idx_mother().at(genx))) == pdg_t) n_numufromt++;

            if (abs(genps_id_mother().at(genx)) == pdg_W && abs(genps_id().at(genx)) == pdg_nutau && genps_status().at(genx) == 1 && abs(genps_id_mother().at(genps_idx_mother().at(genx))) == pdg_t) n_nutaufromt++;
            // count leptons and neutrinos
            if (abs(genps_id().at(genx)) == pdg_el && genps_fromHardProcessFinalState().at(genx) && genps_isLastCopy().at(genx)) nLepsHardProcess++;
            if (abs(genps_id().at(genx)) == pdg_mu && genps_fromHardProcessFinalState().at(genx) && genps_isLastCopy().at(genx)) nLepsHardProcess++;
            if (abs(genps_id().at(genx)) == pdg_tau && genps_fromHardProcessDecayed().at(genx) && genps_isLastCopy().at(genx)) nLepsHardProcess++;

            if (abs(genps_id().at(genx)) == pdg_nue && genps_fromHardProcessFinalState().at(genx) && genps_isLastCopy().at(genx)) nNusHardProcess++;
            if (abs(genps_id().at(genx)) == pdg_numu && genps_fromHardProcessFinalState().at(genx) && genps_isLastCopy().at(genx)) nNusHardProcess++;
            if (abs(genps_id().at(genx)) == pdg_nutau && genps_fromHardProcessFinalState().at(genx) && genps_isLastCopy().at(genx)) nNusHardProcess++;

          }
        }
        //calculate system pt
        LorentzVector hardsystem;
        hardsystem.SetPxPyPzE(0., 0., 0., 0.);
        for (unsigned int genx = 0; genx < gen_susy.id.size(); genx++) {
          if (isstopevent) {
            if ((abs(gen_susy.id.at(genx)) == pdg_chi_1plus1 || abs(gen_susy.id.at(genx)) == pdg_chi_2neutral) && gen_susy.isLastCopy.at(genx)) hardsystem += gen_susy.p4.at(genx);
          }
          if (istopevent) {
            if (abs(gen_qs.id.at(genx)) == pdg_t && gen_qs.isLastCopy.at(genx)) hardsystem += gen_qs.p4.at(genx);
          }
          if (isWevent) {
            if (abs(gen_bosons.id.at(genx)) == pdg_W && abs(gen_bosons.id.at(genx)) != pdg_t && gen_bosons.isLastCopy.at(genx)) hardsystem += +gen_bosons.p4.at(genx);
          }
          if (isZevent) {
            if (abs(gen_bosons.id.at(genx)) == pdg_Z && abs(gen_bosons.id.at(genx)) != pdg_t && gen_bosons.isLastCopy.at(genx)) hardsystem += gen_bosons.p4.at(genx);
          }
        }
        StopEvt.hardgenpt = hardsystem.Pt();
        StopEvt.weight_ISR = 1.;
        //note - these weights don't contain the renormalization, hardcoded
        if (StopEvt.hardgenpt > 600.) {
          StopEvt.weight_ISR = .783;
        } else if (StopEvt.hardgenpt > 400.) {
          StopEvt.weight_ISR = .912;
        } else if (StopEvt.hardgenpt > 300.) {
          StopEvt.weight_ISR = 1.;
        } else if (StopEvt.hardgenpt > 200.) {
          StopEvt.weight_ISR = 1.057;
        } else if (StopEvt.hardgenpt > 150.) {
          StopEvt.weight_ISR = 1.150;
        } else if (StopEvt.hardgenpt > 100.) {
          StopEvt.weight_ISR = 1.179;
        } else if (StopEvt.hardgenpt > 50.) {
          StopEvt.weight_ISR = 1.052;
        } else {
          StopEvt.weight_ISR = 1.;
        }
        float isrweight_unc = fabs(StopEvt.weight_ISR - 1);
        StopEvt.weight_ISRup = 1. + isrweight_unc;
        StopEvt.weight_ISRdown = 1. - isrweight_unc;
        counterhist -> Fill(19, StopEvt.weight_ISR);
        counterhist -> Fill(20, StopEvt.weight_ISRup);
        counterhist -> Fill(21, StopEvt.weight_ISRdown);

        if (sampleConfig.issignal) {
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 19, StopEvt.weight_ISR);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 20, StopEvt.weight_ISRup);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 21, StopEvt.weight_ISRdown);
        }
      } //no data

      // Gen lepton counting and event classification
      StopEvt.genlepsfromtop = n_nuelfromt + n_numufromt + n_nutaufromt;
      gen_leps.gen_nfromt = n_nuelfromt + n_numufromt + n_nutaufromt;
      StopEvt.genLepsHardProcess = nLepsHardProcess;
      StopEvt.genNusHardProcess = nNusHardProcess;

      if (nLepsHardProcess == 0) ee0lep = true;
      if (nLepsHardProcess == 1) ee1lep = true;
      if (nLepsHardProcess >= 2) ge2lep = true;

      if (thisFile.Contains("DYJets") ||
        thisFile.Contains("ZJets") ||
        thisFile.Contains("ZZ")) {
        if (nNusHardProcess >= 2) zToNuNu = true;
      }
      if (thisFile.Contains("WZ") ||
        thisFile.Contains("TTZ") ||
        thisFile.Contains("tZq")) {
        if (nNusHardProcess - nLepsHardProcess >= 2) zToNuNu = true;
      }

      if (zToNuNu) {
        StopEvt.isZtoNuNu = 1;
        StopEvt.is0lep = 0;
        StopEvt.is1lep = 0;
        StopEvt.is2lep = 0;
      } else if (ee0lep) {
        StopEvt.isZtoNuNu = 0;
        StopEvt.is0lep = 1;
        StopEvt.is1lep = 0;
        StopEvt.is2lep = 0;
      } else if (ee1lep) {
        StopEvt.isZtoNuNu = 0;
        StopEvt.is0lep = 0;
        StopEvt.is1lep = 1;
        StopEvt.is2lep = 0;
      } else if (ge2lep) {
        StopEvt.isZtoNuNu = 0;
        StopEvt.is0lep = 0;
        StopEvt.is1lep = 0;
        StopEvt.is2lep = 1;
      }

      if ((ee1lep) && ((StopEvt.genLepsHardProcess - StopEvt.genlepsfromtop) > 0)) StopEvt.is1lepFromW = 1;
      else StopEvt.is1lepFromW = 0;

      if ((ee1lep) && ((StopEvt.genLepsHardProcess - StopEvt.genlepsfromtop) == 0)) StopEvt.is1lepFromTop = 1;
      else StopEvt.is1lepFromTop = 0;

      // nVertex Cut
      if (StopEvt.nvtxs < skim_nvtx) continue;
      if (StopEvt.firstGoodVtxIdx != 0) continue; //really check that first vertex is good
      nEvents_pass_skim_nVtx++;
      //recalculate type 1 corrected met using new JEC
      if (jetConfig.dojec) {
        pair < float, float > newmet;
        pair < float, float > newmet_jup;
        pair < float, float > newmet_jdown;
        if (isbadrawMET) newmet = getT1CHSMET_fromMINIAOD(jet_corrector_pfL1FastJetL2L3, NULL, 0, isbadrawMET);
        else newmet = getT1CHSMET_fromMINIAOD(jet_corrector_pfL1FastJetL2L3);
        newmet_jup = getT1CHSMET_fromMINIAOD(jet_corrector_pfL1FastJetL2L3, jetcorr_uncertainty_sys, true, isbadrawMET);
        newmet_jdown = getT1CHSMET_fromMINIAOD(jet_corrector_pfL1FastJetL2L3, jetcorr_uncertainty_sys, false, isbadrawMET);
        StopEvt.pfmet = newmet.first;
        StopEvt.pfmet_phi = newmet.second;
        StopEvt.pfmet_jup = newmet_jup.first;
        StopEvt.pfmet_phi_jup = newmet_jup.second;
        StopEvt.pfmet_jdown = newmet_jdown.first;
        StopEvt.pfmet_phi_jdown = newmet_jdown.second;
        if (TMath::IsNaN(StopEvt.pfmet) || (!TMath::Finite(StopEvt.pfmet)) || StopEvt.pfmet > 14000.) continue;
      } else {
        StopEvt.pfmet = evt_pfmet();
        StopEvt.pfmet_phi = evt_pfmetPhi();
      }
      // pfmet over calomet filter
      StopEvt.filt_pfovercalomet = !(StopEvt.calomet > 0 && StopEvt.pfmet / StopEvt.calomet > 5);
      /* if(evt_isRealData() && StopEvt.cms3tag.find("CMS4") == 0){
           StopEvt.pfmet = evt_muegclean_pfmet();
           StopEvt.pfmet_phi = evt_muegclean_pfmetPhi();
       }
       */
      //
      //Lepton Variables
      //
      if (debug) std::cout << "[babymaker::looper]: selecting leptons" << std::endl;

      int nGoodLeptons = 0; //stored lep1,lep2
      int nVetoLeptons = 0;

      vector < Lepton > GoodLeps;
      vector < Lepton > LooseLeps;
      vector < Lepton > VetoLeps;
      vector < Lepton > AllLeps;
      vector < unsigned int > idx_alloverlapjets;
      vector < unsigned int > idx_alloverlapjets_jup;
      vector < unsigned int > idx_alloverlapjets_jdown;

      // Electrons
      for (unsigned int eidx = 0; eidx < els_p4().size(); eidx++) {
        if (!PassElectronVetoSelections(eidx, lepConfig.ptcut_veto, lepConfig.etacut_veto)) continue;
        Lepton mylepton;
        mylepton.id = -11 * els_charge().at(eidx);
        mylepton.idx = eidx;
        mylepton.p4 = els_p4().at(eidx);
        int overlapping_jet = getOverlappingJetIndex(mylepton.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, jetConfig.dojec, jetcorr_uncertainty, JES_type, sampleConfig.isfastsim); //don't care about jid
        if (overlapping_jet >= 0) idx_alloverlapjets.push_back(overlapping_jet); //overlap removal for all jets w.r.t. all leptons
        overlapping_jet = getOverlappingJetIndex(mylepton.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, 1, sampleConfig.isfastsim); //don't care about jid
        if (overlapping_jet >= 0) idx_alloverlapjets_jup.push_back(overlapping_jet); //overlap removal for all jets w.r.t. all leptons
        overlapping_jet = getOverlappingJetIndex(mylepton.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, -1, sampleConfig.isfastsim); //don't care about jid
        if (overlapping_jet >= 0) idx_alloverlapjets_jdown.push_back(overlapping_jet); //overlap removal for all jets w.r.t. all leptons

        if ((PassElectronVetoSelections(eidx, lepConfig.ptcut_veto, lepConfig.etacut_veto)) &&
          (!PassElectronPreSelections(eidx, lepConfig.ptcut_loose, lepConfig.etacut_loose))) VetoLeps.push_back(mylepton);
        if ((PassElectronPreSelections(eidx, lepConfig.ptcut_loose, lepConfig.etacut_loose)) &&
          (!PassElectronPreSelections(eidx, lepConfig.ptcut, lepConfig.etacut_el))) LooseLeps.push_back(mylepton);
        if (PassElectronPreSelections(eidx, lepConfig.ptcut, lepConfig.etacut_el)) GoodLeps.push_back(mylepton);
      }

      // Muons
      for (unsigned int midx = 0; midx < mus_p4().size(); midx++) {
        if (!PassMuonVetoSelections(midx, lepConfig.ptcut_veto, lepConfig.etacut_veto)) continue;
        Lepton mylepton;
        mylepton.id = -13 * mus_charge().at(midx);
        mylepton.idx = midx;
        mylepton.p4 = mus_p4().at(midx);
        int overlapping_jet = getOverlappingJetIndex(mylepton.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, jetConfig.dojec, jetcorr_uncertainty, JES_type, sampleConfig.isfastsim); //don't care about jid
        if (overlapping_jet >= 0) idx_alloverlapjets.push_back(overlapping_jet); //overlap removal for all jets w.r.t. all leptons
        overlapping_jet = getOverlappingJetIndex(mylepton.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, 1, sampleConfig.isfastsim); //don't care about jid
        if (overlapping_jet >= 0) idx_alloverlapjets_jup.push_back(overlapping_jet); //overlap removal for all jets w.r.t. all leptons
        overlapping_jet = getOverlappingJetIndex(mylepton.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, -1, sampleConfig.isfastsim); //don't care about jid
        if (overlapping_jet >= 0) idx_alloverlapjets_jdown.push_back(overlapping_jet); //overlap removal for all jets w.r.t. all leptons	
        if ((PassMuonVetoSelections(midx, lepConfig.ptcut_veto, lepConfig.etacut_veto)) &&
          (!PassMuonPreSelections(midx, lepConfig.ptcut_loose, lepConfig.etacut_loose))) VetoLeps.push_back(mylepton);
        if ((PassMuonPreSelections(midx, lepConfig.ptcut_loose, lepConfig.etacut_loose)) &&
          (!PassMuonPreSelections(midx, lepConfig.ptcut, lepConfig.etacut_mu))) LooseLeps.push_back(mylepton);
        if (PassMuonPreSelections(midx, lepConfig.ptcut, lepConfig.etacut_mu)) GoodLeps.push_back(mylepton);

      }
      //sort by pt
      sort(GoodLeps.begin(), GoodLeps.end(), sortLepbypt());
      sort(LooseLeps.begin(), LooseLeps.end(), sortLepbypt());
      sort(VetoLeps.begin(), VetoLeps.end(), sortLepbypt());
      // remove overlaps
      if (GoodLeps.size() > 0) {
        for (unsigned int lep = 1; lep < GoodLeps.size(); lep++) {
          if (ROOT::Math::VectorUtil::DeltaR(GoodLeps.at(0).p4, GoodLeps.at(lep).p4) < 0.01) GoodLeps.erase(GoodLeps.begin() + lep);
        }
        for (unsigned int lep = 0; lep < LooseLeps.size(); lep++) {
          if (ROOT::Math::VectorUtil::DeltaR(GoodLeps.at(0).p4, LooseLeps.at(lep).p4) < 0.01) LooseLeps.erase(LooseLeps.begin() + lep);
        }
        for (unsigned int lep = 0; lep < VetoLeps.size(); lep++) {
          if (ROOT::Math::VectorUtil::DeltaR(GoodLeps.at(0).p4, VetoLeps.at(lep).p4) < 0.01) VetoLeps.erase(VetoLeps.begin() + lep);
        }
      }
      nGoodLeptons = GoodLeps.size();
      int nLooseLeptons = GoodLeps.size() + LooseLeps.size(); //use for Zll
      nVetoLeptons = GoodLeps.size() + LooseLeps.size() + VetoLeps.size();

      StopEvt.ngoodleps = nGoodLeptons;
      StopEvt.nlooseleps = nLooseLeptons; // these are used for Zll
      StopEvt.nvetoleps = nVetoLeptons;

      if (debug) std::cout << "[babymaker::looper]: filling lepton variables LINE:" << __LINE__ << std::endl;
      AllLeps.clear();
      AllLeps.insert(AllLeps.end(), GoodLeps.begin(), GoodLeps.end());
      AllLeps.insert(AllLeps.end(), LooseLeps.begin(), LooseLeps.end());
      AllLeps.insert(AllLeps.end(), VetoLeps.begin(), VetoLeps.end());
      if (AllLeps.size() > 0) lep1.FillCommon(AllLeps.at(0).id, AllLeps.at(0).idx);
      if (AllLeps.size() > 1) lep2.FillCommon(AllLeps.at(1).id, AllLeps.at(1).idx);
      if (debug) std::cout << "[babymaker::looper]: filling lepton variables LINE:" << __LINE__ << std::endl;

      // Lepton SFs
      float lepSF_pt_cutoff(99.999), lepSF_pt_min(10.001);
      double lepSF(1.0), lepSF_Up(1.0), lepSF_Dn(1.0);
      double vetoLepSF(1.0), vetoLepSF_Up(1.0), vetoLepSF_Dn(1.0);
      float lepSF_FS_pt_cutoff(199.999), lepSF_FS_pt_min(10.001);
      double lepSF_FS(1.0), lepSF_FS_Up(1.0), lepSF_FS_Dn(1.0);

      // Lep1 SF
      if (lepConfig.dosf && !StopEvt.is_data && nVetoLeptons > 0) {
        if (abs(lep1.pdgid) == pdg_el) {
          vector < float > updownerr = getupdownerr(h_el_SF, (float) lep1.p4.Pt(), (float) lep1.p4.Eta(), lepSF_pt_cutoff, lepSF_pt_min, 2.39, true);
          vector < float > updownerr_track = getupdownerr(h_el_SF_tracking, (float) lep1.p4.Pt(), (float) lep1.p4.Eta(), 199, 21, 2.39, false);
          lepSF = updownerr.at(0);
          lepSF_Up = updownerr.at(1);
          lepSF_Dn = updownerr.at(2);
          lepSF *= updownerr_track.at(0);
          lepSF_Up *= updownerr_track.at(1);
          lepSF_Dn *= updownerr_track.at(2);
          if (sampleConfig.isfastsim) {
            vector < float > updownerr_fs = getupdownerr(h_el_FS, (float) lep1.p4.Pt(), (float) lep1.p4.Eta(), lepSF_FS_pt_cutoff, lepSF_FS_pt_min, 2.39, true);
            lepSF *= updownerr_fs.at(0);
            lepSF_Up *= updownerr_fs.at(1);
            lepSF_Dn *= updownerr_fs.at(2);
          }
        }
        if (abs(lep1.pdgid) == pdg_mu) {
          vector < float > updownerr_mu = getupdownerr(h_mu_SF, (float) lep1.p4.Pt(), (float) lep1.p4.Eta(), lepSF_pt_cutoff, lepSF_pt_min, 2.1, true);
          lepSF = updownerr_mu.at(0);
          lepSF_Up = updownerr_mu.at(1);
          lepSF_Dn = updownerr_mu.at(2);
          int binX = h_mu_SF_tracking -> GetXaxis() -> FindBin(std::max(-2.1, (double) lep1.p4.Eta()));
          lepSF *= h_mu_SF_tracking -> GetBinContent(binX);
          lepSF_Up *= (h_mu_SF_tracking -> GetBinContent(binX) + h_mu_SF_tracking -> GetBinError(binX));
          lepSF_Dn *= (h_mu_SF_tracking -> GetBinContent(binX) - h_mu_SF_tracking -> GetBinError(binX));
          if (sampleConfig.isfastsim) {
            vector < float > updownerr_fs_mu = getupdownerr(h_mu_FS, (float) lep1.p4.Pt(), (float) lep1.p4.Eta(), lepSF_FS_pt_cutoff, lepSF_FS_pt_min, 2.1, true);
            lepSF *= updownerr_fs_mu.at(0);
            lepSF_Up *= updownerr_fs_mu.at(1);
            lepSF_Dn *= updownerr_fs_mu.at(2);
          }
        }
      }
      // Lep2 SF
      if (lepConfig.dosf && !StopEvt.is_data && nVetoLeptons > 1) {
        if (abs(lep2.pdgid) == pdg_el) {
          if (nGoodLeptons > 1) {
            vector < float > updownerr_lep2 = getupdownerr(h_el_SF, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_pt_cutoff, lepSF_pt_min, 2.39, true);
            lepSF *= updownerr_lep2.at(0);
            lepSF_Up *= updownerr_lep2.at(1);
            lepSF_Dn *= updownerr_lep2.at(2);
            if (sampleConfig.isfastsim) {
              vector < float > updownerr_fs_lep2 = getupdownerr(h_el_FS, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_FS_pt_cutoff, lepSF_FS_pt_min, 2.39, true);
              lepSF *= updownerr_fs_lep2.at(0);
              lepSF_Up *= updownerr_fs_lep2.at(1);
              lepSF_Dn *= updownerr_fs_lep2.at(2);
            }
          } // end if 2 good electrons
          else {
            vector < float > updownerr_lep2_veto = getupdownerr(h_el_SF_veto, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_pt_cutoff, lepSF_pt_min, 2.39, true);
            lepSF *= updownerr_lep2_veto.at(0);
            lepSF_Up *= updownerr_lep2_veto.at(1);
            lepSF_Dn *= updownerr_lep2_veto.at(2);
            if (sampleConfig.isfastsim) {
              vector < float > updownerr_fs_veto = getupdownerr(h_el_veto_FS, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_FS_pt_cutoff, lepSF_FS_pt_min, 2.39, true);
              lepSF *= updownerr_fs_veto.at(0);
              lepSF_Up *= updownerr_fs_veto.at(1);
              lepSF_Dn *= updownerr_fs_veto.at(2);
            }
          }
          vector < float > updownerr_track_lep2 = getupdownerr(h_el_SF_tracking, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), 199, 21, 2.39, false);
          lepSF *= updownerr_track_lep2.at(0);
          lepSF_Up *= updownerr_track_lep2.at(1);
          lepSF_Dn *= updownerr_track_lep2.at(2);
        } // end if 2nd lep if el
        if (abs(lep2.pdgid) == pdg_mu) {
          if (nGoodLeptons > 1) {
            vector < float > updownerr_lep2_mu = getupdownerr(h_mu_SF, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_pt_cutoff, lepSF_pt_min, 2.2, true);
            lepSF *= updownerr_lep2_mu.at(0);
            lepSF_Up *= updownerr_lep2_mu.at(1);
            lepSF_Dn *= updownerr_lep2_mu.at(2);
            if (sampleConfig.isfastsim) {
              vector < float > updownerr_fs_lep2_mu = getupdownerr(h_mu_FS, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_FS_pt_cutoff, lepSF_FS_pt_min, 2.2, true);
              lepSF *= updownerr_fs_lep2_mu.at(0);
              lepSF_Up *= updownerr_fs_lep2_mu.at(1);
              lepSF_Dn *= updownerr_fs_lep2_mu.at(2);
            }
          } // end if 2 good leptons
          else {
            vector < float > updownerr_veto_mu = getupdownerr(h_mu_SF_veto, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_pt_cutoff, lepSF_pt_min, 2.2, true);
            lepSF *= updownerr_veto_mu.at(0);
            lepSF_Up *= updownerr_veto_mu.at(1);
            lepSF_Dn *= updownerr_veto_mu.at(2);
            if (sampleConfig.isfastsim) {
              vector < float > updownerr_veto_mu_fs = getupdownerr(h_mu_veto_FS, (float) lep2.p4.Pt(), (float) lep2.p4.Eta(), lepSF_pt_cutoff, lepSF_pt_min, 2.2, true);
              lepSF *= updownerr_veto_mu_fs.at(0);
              lepSF_Up *= updownerr_veto_mu_fs.at(1);
              lepSF_Dn *= updownerr_veto_mu_fs.at(2);
            }
          }
        } // end if 2nd lep is mu
        int binX = h_mu_SF_tracking -> GetXaxis() -> FindBin(std::max(-2.23, (double) lep2.p4.Eta()));
        lepSF *= h_mu_SF_tracking -> GetBinContent(binX);
        lepSF_Up *= (h_mu_SF_tracking -> GetBinContent(binX) + h_mu_SF_tracking -> GetBinError(binX));
        lepSF_Dn *= (h_mu_SF_tracking -> GetBinContent(binX) - h_mu_SF_tracking -> GetBinError(binX));
      } // end if 2nd lepton reco lepton exists

      // If only 1 reco lepton, and is2lep event, then find lost gen lepton
      if (lepConfig.dosf && !StopEvt.is_data && nVetoLeptons == 1 && StopEvt.is2lep) {
        for (int iGen = 0; iGen < (int) gen_leps.p4.size(); iGen++) {
          if (abs(gen_leps.id.at(iGen)) != 11 && abs(gen_leps.id.at(iGen)) != 13) continue;
          if (!gen_leps.fromHardProcessFinalState.at(iGen)) continue;
          if (!gen_leps.isLastCopy.at(iGen)) continue;
          if (ROOT::Math::VectorUtil::DeltaR(gen_leps.p4.at(iGen), lep1.p4) < matched_dr) continue;
          if (gen_leps.p4.at(iGen).Pt() < 5 || fabs(gen_leps.p4.at(iGen).Eta()) > 2.4) continue;

          TH2D * h_vetoLep_eff = NULL;
          if (abs(gen_leps.id.at(iGen)) == 11) h_vetoLep_eff = h_el_vetoLepEff;
          if (abs(gen_leps.id.at(iGen)) == 13) h_vetoLep_eff = h_mu_vetoLepEff;

          int binX_eff = h_vetoLep_eff -> GetXaxis() -> FindBin(std::max(std::min(lepSF_pt_cutoff, (float) gen_leps.p4.at(iGen).Pt()), lepSF_pt_min));
          int binY_eff = h_vetoLep_eff -> GetYaxis() -> FindBin(fabs(gen_leps.p4.at(iGen).Eta()));
          double vetoEff = h_vetoLep_eff -> GetBinContent(binX_eff, binY_eff);

          TH2D * h_lep_sf = NULL;
          if (abs(gen_leps.id.at(iGen)) == 11) h_lep_sf = h_el_SF_veto;
          if (abs(gen_leps.id.at(iGen)) == 13) h_lep_sf = h_mu_SF_veto;

          int binX_sf = h_lep_sf -> GetXaxis() -> FindBin(std::max(std::min(lepSF_pt_cutoff, (float) gen_leps.p4.at(iGen).Pt()), lepSF_pt_min));
          int binY_sf = h_lep_sf -> GetYaxis() -> FindBin(fabs(gen_leps.p4.at(iGen).Eta()));

          double vetoLepSF_temp = h_lep_sf -> GetBinContent(binX_sf, binY_sf);
          double vetoLepSF_temp_Up = vetoLepSF_temp + h_lep_sf -> GetBinError(binX_sf, binY_sf);
          double vetoLepSF_temp_Dn = vetoLepSF_temp - h_lep_sf -> GetBinError(binX_sf, binY_sf);

          if (vetoEff == 1.0) {
            vetoLepSF = 1.0;
            vetoLepSF_Up = 1.0;
            vetoLepSF_Dn = 1.0;
          } else {
            vetoLepSF = (1 - (vetoEff * vetoLepSF_temp)) / (1 - vetoEff);
            vetoLepSF_Up = (1 - (vetoEff * vetoLepSF_temp_Up)) / (1 - vetoEff);
            vetoLepSF_Dn = (1 - (vetoEff * vetoLepSF_temp_Dn)) / (1 - vetoEff);
          }

          break; // break after finding 2nd hard gen lepton

        } // end loop over gen leptons

      } // end if finding gen lost lepton for vetoEff SF

      StopEvt.weight_lepSF = lepSF;
      StopEvt.weight_lepSF_up = lepSF_Up;
      StopEvt.weight_lepSF_down = lepSF_Dn;

      StopEvt.weight_vetoLepSF = vetoLepSF;
      StopEvt.weight_vetoLepSF_up = vetoLepSF_Up;
      StopEvt.weight_vetoLepSF_down = vetoLepSF_Dn;

      StopEvt.weight_lepSF_fastSim = lepSF_FS;
      StopEvt.weight_lepSF_fastSim_up = lepSF_FS_Up;
      StopEvt.weight_lepSF_fastSim_down = lepSF_FS_Dn;
      if (debug) std::cout << "[babymaker::looper]: filling got lepton sf variables LINE:" << __LINE__ << std::endl;
      // save the sum of weights for normalization offline to n-babies.
      if (!evt_isRealData() && lepConfig.dosf) {
        counterhist -> Fill(28, StopEvt.weight_lepSF);
        counterhist -> Fill(29, StopEvt.weight_lepSF_up);
        counterhist -> Fill(30, StopEvt.weight_lepSF_down);
      }
      if (!evt_isRealData() && lepConfig.dosf) {
        counterhist -> Fill(31, StopEvt.weight_vetoLepSF);
        counterhist -> Fill(32, StopEvt.weight_vetoLepSF_up);
        counterhist -> Fill(33, StopEvt.weight_vetoLepSF_down);
      }
      if (!evt_isRealData() && lepConfig.dosf && sampleConfig.isfastsim) {
        counterhist -> Fill(34, StopEvt.weight_lepSF_fastSim);
        counterhist -> Fill(35, StopEvt.weight_lepSF_fastSim_up);
        counterhist -> Fill(36, StopEvt.weight_lepSF_fastSim_down);
      }

      // Signal
      if (sampleConfig.issignal && !evt_isRealData() && lepConfig.dosf) {
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 27, StopEvt.weight_lepSF);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 28, StopEvt.weight_lepSF_up);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 29, StopEvt.weight_lepSF_down);
      }
      if (sampleConfig.issignal && !evt_isRealData() && lepConfig.dosf) {
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 30, StopEvt.weight_vetoLepSF);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 31, StopEvt.weight_vetoLepSF_up);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 32, StopEvt.weight_vetoLepSF_down);
      }
      if (sampleConfig.issignal && !evt_isRealData() && lepConfig.dosf && sampleConfig.isfastsim) {
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 33, StopEvt.weight_lepSF_fastSim);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 34, StopEvt.weight_lepSF_fastSim_up);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 35, StopEvt.weight_lepSF_fastSim_down);
      }
      //
      // Jet Selection
      //
      float btagprob_data(1), btagprob_mc(1), btagprob_heavy_UP(1), btagprob_heavy_DN(1), btagprob_light_UP(1), btagprob_light_DN(1), btagprob_FS_UP(1), btagprob_FS_DN(1);
      float btagprob_data_jup(1), btagprob_mc_jup(1), btagprob_heavy_UP_jup(1), btagprob_heavy_DN_jup(1), btagprob_light_UP_jup(1), btagprob_light_DN_jup(1), btagprob_FS_UP_jup(1), btagprob_FS_DN_jup(1);
      float btagprob_data_jdown(1), btagprob_mc_jdown(1), btagprob_heavy_UP_jdown(1), btagprob_heavy_DN_jdown(1), btagprob_light_UP_jdown(1), btagprob_light_DN_jdown(1), btagprob_FS_UP_jdown(1), btagprob_FS_DN_jdown(1);
      //std::cout << "[babymaker::looper]: filling jets vars" << std::endl;         
      // Get the jets overlapping with the selected leptons
      if (debug) std::cout << "[babymaker::looper]: about to fill jet variables variables LINE:" << __LINE__ << std::endl;
      if (pfjets_p4().size() > 0) {
        jet_overlep1_idx = -9999;
        jet_overlep2_idx = -9999;
        if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
        if (nVetoLeptons > 0) jet_overlep1_idx = getOverlappingJetIndex(lep1.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, jetConfig.dojec, jetcorr_uncertainty, JES_type, sampleConfig.isfastsim); //don't care about jid
        if (nVetoLeptons > 1) jet_overlep2_idx = getOverlappingJetIndex(lep2.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, jetConfig.dojec, jetcorr_uncertainty, JES_type, sampleConfig.isfastsim); //don't care about jid

        if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
        // Jets and b-tag variables feeding the index for the jet overlapping the selected leptons
        jets.SetJetSelection("ak4", jetConfig.ptcut, jetConfig.etacut, true); //save only jets passing jid
        jets.SetJetSelection("ak8", jetConfig.ptcut_ak8, jetConfig.etacut_ak8, true); //save only jets passing jid
        if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
        jets.FillCommon(idx_alloverlapjets, jet_corrector_pfL1FastJetL2L3, btagprob_data, btagprob_mc, btagprob_heavy_UP, btagprob_heavy_DN, btagprob_light_UP, btagprob_light_DN, btagprob_FS_UP, btagprob_FS_DN, jet_overlep1_idx, jet_overlep2_idx, jetConfig.dojec, jetcorr_uncertainty, JES_type, jetConfig.dobtagsf, sampleConfig.isfastsim);

        if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
        //JEC up
        jet_overlep1_idx = -9999;
        jet_overlep2_idx = -9999;
        if (nVetoLeptons > 0) jet_overlep1_idx = getOverlappingJetIndex(lep1.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, 1, sampleConfig.isfastsim); //don't care about jid
        if (nVetoLeptons > 1) jet_overlep2_idx = getOverlappingJetIndex(lep2.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, 1, sampleConfig.isfastsim); //don't care about jid

        // Jets and b-tag variables feeding the index for the jet overlapping the selected leptons
        jets_jup.SetJetSelection("ak4", jetConfig.ptcut, jetConfig.etacut, true); //save only jets passing jid
        jets_jup.SetJetSelection("ak8", jetConfig.ptcut_ak8, jetConfig.etacut_ak8, true); //save only jets passing jid
        jets_jup.FillCommon(idx_alloverlapjets_jup, jet_corrector_pfL1FastJetL2L3, btagprob_data_jup, btagprob_mc_jup, btagprob_heavy_UP_jup, btagprob_heavy_DN_jup, btagprob_light_UP_jup, btagprob_light_DN_jup, btagprob_FS_UP_jup, btagprob_FS_DN_jup, jet_overlep1_idx, jet_overlep2_idx, true, jetcorr_uncertainty_sys, 1, false, sampleConfig.isfastsim);

        //JEC down
        jet_overlep1_idx = -9999;
        jet_overlep2_idx = -9999;
        if (nVetoLeptons > 0) jet_overlep1_idx = getOverlappingJetIndex(lep1.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, -1, sampleConfig.isfastsim); //don't care about jid
        if (nVetoLeptons > 1) jet_overlep2_idx = getOverlappingJetIndex(lep2.p4, pfjets_p4(), 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, true, jetcorr_uncertainty_sys, -1, sampleConfig.isfastsim); //don't care about jid

        // Jets and b-tag variables feeding the index for the jet overlapping the selected leptons
        jets_jdown.SetJetSelection("ak4", jetConfig.ptcut, jetConfig.etacut, true); //save only jets passing jid
        jets_jdown.SetJetSelection("ak8", jetConfig.ptcut_ak8, jetConfig.etacut_ak8, true); //save only jets passing jid
        jets_jdown.FillCommon(idx_alloverlapjets_jdown, jet_corrector_pfL1FastJetL2L3, btagprob_data_jdown, btagprob_mc_jdown, btagprob_heavy_UP_jdown, btagprob_heavy_DN_jdown, btagprob_light_UP_jdown, btagprob_light_DN_jdown, btagprob_FS_UP_jdown, btagprob_FS_DN_jdown, jet_overlep1_idx, jet_overlep2_idx, true, jetcorr_uncertainty_sys, -1, false, sampleConfig.isfastsim);
      }

      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;

      if (!evt_isRealData()) {
        int NISRjets = 0;
        int nonISRjets = 0;
        for (unsigned int jix = 0; jix < jets.ak4pfjets_p4.size(); ++jix) {
          bool ismatched = false;
          for (unsigned int gpix = 0; gpix < genpnotISR.size(); ++gpix) {
            if (dRbetweenVectors(jets.ak4pfjets_p4[jix], genpnotISR[gpix]) <= 0.3) {
              ismatched = true;
              break;
            }
          }
          if (!ismatched) ++NISRjets;
          else ++nonISRjets;
        }
        if (debug) std::cout << "[babymaker::looper]: filling ISRnjets variables LINE:" << __LINE__ << std::endl;
        // moriond 2017
        if (NISRjets == 0) {
          StopEvt.weight_ISRnjets = 1.000;
          StopEvt.weight_ISRnjets_UP = 1.0000;
          StopEvt.weight_ISRnjets_DN = 1.0000;
        } else if (NISRjets == 1) {
          StopEvt.weight_ISRnjets = 0.920;
          StopEvt.weight_ISRnjets_UP = 0.960;
          StopEvt.weight_ISRnjets_DN = 0.880;
        } else if (NISRjets == 2) {
          StopEvt.weight_ISRnjets = 0.821;
          StopEvt.weight_ISRnjets_UP = 0.911;
          StopEvt.weight_ISRnjets_DN = 0.731;
        } else if (NISRjets == 3) {
          StopEvt.weight_ISRnjets = 0.715;
          StopEvt.weight_ISRnjets_UP = 0.858;
          StopEvt.weight_ISRnjets_DN = 0.572;
        } else if (NISRjets == 4) {
          StopEvt.weight_ISRnjets = 0.662;
          StopEvt.weight_ISRnjets_UP = 0.832;
          StopEvt.weight_ISRnjets_DN = 0.492;
        } else if (NISRjets == 5) {
          StopEvt.weight_ISRnjets = 0.561;
          StopEvt.weight_ISRnjets_UP = 0.782;
          StopEvt.weight_ISRnjets_DN = 0.340;
        } else {
          StopEvt.weight_ISRnjets = 0.511;
          StopEvt.weight_ISRnjets_UP = 0.769;
          StopEvt.weight_ISRnjets_DN = 0.253;
        }
        StopEvt.NISRjets = NISRjets;
        StopEvt.NnonISRjets = nonISRjets;

        counterhist -> Fill(25, StopEvt.weight_ISRnjets);
        counterhist -> Fill(26, StopEvt.weight_ISRnjets_UP);
        counterhist -> Fill(27, StopEvt.weight_ISRnjets_DN);

        if (sampleConfig.issignal) {
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 24, StopEvt.weight_ISRnjets);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 25, StopEvt.weight_ISRnjets_UP);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 26, StopEvt.weight_ISRnjets_DN);
        }
      }
      // SAVE B TAGGING SF 
      if (!evt_isRealData() && jetConfig.dobtagsf) {
        StopEvt.weight_btagsf = btagprob_data / btagprob_mc;
        StopEvt.weight_btagsf_heavy_UP = btagprob_heavy_UP / btagprob_mc;
        StopEvt.weight_btagsf_light_UP = btagprob_light_UP / btagprob_mc;
        StopEvt.weight_btagsf_heavy_DN = btagprob_heavy_DN / btagprob_mc;
        StopEvt.weight_btagsf_light_DN = btagprob_light_DN / btagprob_mc;
        if (sampleConfig.isfastsim) {
          StopEvt.weight_btagsf_fastsim_UP = btagprob_FS_UP / btagprob_mc;
          StopEvt.weight_btagsf_fastsim_DN = btagprob_FS_DN / btagprob_mc;
        }
      }
      // save the sum of weights for normalization offline to n-babies.
      if (!evt_isRealData() && jetConfig.dobtagsf) {
        counterhist -> Fill(14, StopEvt.weight_btagsf);
        counterhist -> Fill(15, StopEvt.weight_btagsf_heavy_UP);
        counterhist -> Fill(16, StopEvt.weight_btagsf_light_UP);
        counterhist -> Fill(17, StopEvt.weight_btagsf_heavy_DN);
        counterhist -> Fill(18, StopEvt.weight_btagsf_light_DN);
        if (sampleConfig.isfastsim) {
          counterhist -> Fill(23, StopEvt.weight_btagsf_fastsim_UP);
          counterhist -> Fill(24, StopEvt.weight_btagsf_fastsim_DN);
        }
      }
      if (sampleConfig.issignal && !evt_isRealData() && jetConfig.dobtagsf) {
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 14, StopEvt.weight_btagsf);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 15, StopEvt.weight_btagsf_heavy_UP);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 16, StopEvt.weight_btagsf_light_UP);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 17, StopEvt.weight_btagsf_heavy_DN);
        counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 18, StopEvt.weight_btagsf_light_DN);
        if (sampleConfig.isfastsim) {
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 22, StopEvt.weight_btagsf_fastsim_UP);
          counterhistSig -> Fill(StopEvt.mass_stop, StopEvt.mass_lsp, 23, StopEvt.weight_btagsf_fastsim_DN);
        }
      }
      // apply skim
      if (debug) std::cout << "[babymaker::looper]: going to apply skims LINE:" << __LINE__ << std::endl;
      if (nGoodLeptons < lepConfig.nlep) continue;
      nEvents_pass_skim_nGoodLeps++;
      //if (!(jets.ngoodjets >= jetConfig.njet) && !(jets_jup.ngoodjets >= jetConfig.njet) && !(jets_jdown.ngoodjets >= jetConfig.njet)) continue;
      bool pass1jet = (jets.ak8GoodPFJets >= jetConfig.nfatjet); //Sicong: Add 1 jet region
      bool pass2jet = (jets.ngoodjets >= jetConfig.njet) || (jets_jup.ngoodjets >= jetConfig.njet) || (jets_jdown.ngoodjets >= jetConfig.njet);
      if (!pass1jet && !pass2jet) continue;
      //bool pass_old2jet = false;//Sicong: Specifically for getting 1 jet region, s.t. we could estimate the background yield 
      //if (jets.ngoodjets >= 2){
      //    if (jets.ak4pfjets_p4[1].Pt()>30&&jets.ak4pfjets_p4[0].Pt()>30) pass_old2jet = true;
      //}
      //if (!pass1jet || pass_old2jet) continue;
      //cout<<"Got 1 jet."<<endl; 
      nEvents_pass_skim_nGoodJets++;
      int skim_nBJets = 1;
      if (!(pass1jet) && !(jets.ngoodbtags >= skim_nBJets) && !(jets_jup.ngoodbtags_jup >= skim_nBJets) && !(jets_jdown.ngoodbtags_jdown >= skim_nBJets)) continue;
      nEvents_pass_skim_nBJets++; //
      // fastsim filter for bad jets//
      if (debug) std::cout << "[babymaker::looper]: fastsim filter LINE:" << __LINE__ << std::endl;
      bool isbadmuonjet = false;
      for (unsigned int jdx = 0; jdx < jets.ak4pfjets_p4.size(); ++jdx) {
        jets.dphi_ak4pfjet_met.push_back(getdphi(jets.ak4pfjets_p4[jdx].Phi(), StopEvt.pfmet_phi));
        if (jets.ak4pfjets_p4[jdx].Pt() > 200 && jets.dphi_ak4pfjet_met[jdx] > (TMath::Pi() - 0.4) && jets.ak4pfjets_muf[jdx] > 0.5) isbadmuonjet = true;
      }
      StopEvt.filt_jetWithBadMuon = !isbadmuonjet;
      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
      isbadmuonjet = false;
      for (unsigned int jdx = 0; jdx < jets_jup.ak4pfjets_p4.size(); ++jdx) {
        jets_jup.dphi_ak4pfjet_met.push_back(getdphi(jets_jup.ak4pfjets_p4[jdx].Phi(), StopEvt.pfmet_phi_jup));
        if (jets_jup.ak4pfjets_p4[jdx].Pt() > 200 && jets_jup.dphi_ak4pfjet_met[jdx] > (TMath::Pi() - 0.4) && jets_jup.ak4pfjets_muf[jdx] > 0.5) isbadmuonjet = true;
      }
      StopEvt.filt_jetWithBadMuon_jup = !isbadmuonjet;
      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
      isbadmuonjet = false;
      for (unsigned int jdx = 0; jdx < jets_jdown.ak4pfjets_p4.size(); ++jdx) {
        jets_jdown.dphi_ak4pfjet_met.push_back(getdphi(jets_jdown.ak4pfjets_p4[jdx].Phi(), StopEvt.pfmet_phi_jdown));
        if (jets_jdown.ak4pfjets_p4[jdx].Pt() > 200 && jets_jdown.dphi_ak4pfjet_met[jdx] > (TMath::Pi() - 0.4) && jets_jdown.ak4pfjets_muf[jdx] > 0.5) isbadmuonjet = true;
      }
      StopEvt.filt_jetWithBadMuon_jdown = !isbadmuonjet;

      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
      bool fastsimfilt = false;
      for (unsigned int jix = 0; jix < jets.ak4pfjets_p4.size(); ++jix) {
        if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
        if (jets.ak4pfjets_p4[jix].Pt() < 30) continue;
        if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
        if (fabs(jets.ak4pfjets_p4[jix].Eta()) > 2.4) continue;
        if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
        bool isgenmatch = false;
        for (unsigned int gix = 0; gix < jets.ak4genjets_p4.size(); ++gix) {
          if (dRbetweenVectors(jets.ak4genjets_p4[gix], jets.ak4pfjets_p4[jix]) < 0.3) {
            isgenmatch = true;
            break;
          }
        }
        if (isgenmatch) continue;
        if (jets.ak4pfjets_chf[jix] > 0.1) continue; {
          fastsimfilt = true;
          break;
        }
      }
      StopEvt.filt_fastsimjets = !fastsimfilt;
      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;

      // FastSim filter, JESup Jets
      bool fastsimfilt_jup = false;
      for (unsigned int jix = 0; jix < jets_jup.ak4pfjets_p4.size(); ++jix) {
        if (jets_jup.ak4pfjets_p4[jix].Pt() < 30) continue;
        if (fabs(jets_jup.ak4pfjets_p4[jix].Eta()) > 2.4) continue;
        bool isgenmatch = false;
        for (unsigned int gix = 0; gix < jets_jup.ak4genjets_p4.size(); ++gix) {
          if (dRbetweenVectors(jets_jup.ak4genjets_p4[gix], jets_jup.ak4pfjets_p4[jix]) < 0.3) {
            isgenmatch = true;
            break;
          }
        }
        if (isgenmatch) continue;
        if (jets_jup.ak4pfjets_chf[jix] > 0.1) continue;
        fastsimfilt_jup = true;
        break;
      }
      StopEvt.filt_fastsimjets_jup = !fastsimfilt_jup;
      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;

      // FastSim filter, JESdn Jets
      bool fastsimfilt_jdown = false;
      for (unsigned int jix = 0; jix < jets_jdown.ak4pfjets_p4.size(); ++jix) {
        if (jets_jdown.ak4pfjets_p4[jix].Pt() < 30) continue;
        if (fabs(jets_jdown.ak4pfjets_p4[jix].Eta()) > 2.4) continue;
        bool isgenmatch = false;
        for (unsigned int gix = 0; gix < jets_jdown.ak4genjets_p4.size(); ++gix) {
          if (dRbetweenVectors(jets_jdown.ak4genjets_p4[gix], jets_jdown.ak4pfjets_p4[jix]) < 0.3) {
            isgenmatch = true;
            break;
          }
        }
        if (isgenmatch) continue;
        if (jets_jdown.ak4pfjets_chf[jix] > 0.1) continue;
        fastsimfilt_jdown = true;
        break;
      }
      StopEvt.filt_fastsimjets_jdown = !fastsimfilt_jdown;
      //
      // Photon Selection
      //
      if (debug) std::cout << "[babymaker::looper]: Photon Selection LINE:" << __LINE__ << std::endl;
      ph.SetPhotonSelection(gammaConfig.ptcut, gammaConfig.etacut);
      ph.FillCommon();
      StopEvt.nPhotons = ph.p4.size();
      if (StopEvt.nPhotons < gammaConfig.ngamma) continue;
      int leadph = -1; //use this in case we have a wide photon selection (like loose id), but want to use specific photon
      for (unsigned int i = 0; i < ph.p4.size(); ++i) {
        int overlapping_jet = getOverlappingJetIndex(ph.p4.at(i), jets.ak4pfjets_p4, 0.4, jetConfig.ptcut, jetConfig.etacut, false, jet_corrector_pfL1FastJetL2L3, jetConfig.dojec, jetcorr_uncertainty, JES_type, sampleConfig.isfastsim);
        ph.overlapJetId.at(i) = overlapping_jet;
        if (leadph != -1 && ph.p4.at(i).Pt() < ph.p4.at(leadph).Pt()) continue;
        if (StopEvt.ngoodleps > 0 && ROOT::Math::VectorUtil::DeltaR(ph.p4.at(i), lep1.p4) < 0.2) continue;
        if (StopEvt.ngoodleps > 1 && ROOT::Math::VectorUtil::DeltaR(ph.p4.at(i), lep2.p4) < 0.2) continue;
        leadph = i;
      }
      StopEvt.ph_selectedidx = leadph;

      // Event Variables

      if (debug) std::cout << "[babymaker::looper]: calculate event level variables LINE:" << __LINE__ << std::endl;
      // MET & Leptons
      if (nVetoLeptons > 0) StopEvt.mt_met_lep = calculateMt(lep1.p4, StopEvt.pfmet, StopEvt.pfmet_phi);
      if (nVetoLeptons > 0) StopEvt.mt_met_lep_jup = calculateMt(lep1.p4, StopEvt.pfmet_jup, StopEvt.pfmet_phi_jup);
      if (nVetoLeptons > 0) StopEvt.mt_met_lep_jdown = calculateMt(lep1.p4, StopEvt.pfmet_jdown, StopEvt.pfmet_phi_jdown);
      if (nVetoLeptons > 1) StopEvt.mt_met_lep2 = calculateMt(lep2.p4, StopEvt.pfmet, StopEvt.pfmet_phi);
      if (nVetoLeptons > 0) StopEvt.dphi_Wlep = DPhi_W_lep(StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4);
      if (jets.ak4pfjets_p4.size() > 0) StopEvt.MET_over_sqrtHT = StopEvt.pfmet / TMath::Sqrt(jets.ak4_HT);

      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
      StopEvt.ak4pfjets_rho = evt_fixgridfastjet_all_rho();
      vector < int > jetIndexSortedCSV;
      if (sampleConfig.isfastsim) jetIndexSortedCSV = JetUtil::JetIndexCSVsorted(jets.ak4pfjets_CSV, jets.ak4pfjets_p4, jets.ak4pfjets_loose_pfid, jetConfig.ptcut, jetConfig.etacut, false);
      else jetIndexSortedCSV = JetUtil::JetIndexCSVsorted(jets.ak4pfjets_CSV, jets.ak4pfjets_p4, jets.ak4pfjets_loose_pfid, jetConfig.ptcut, jetConfig.etacut, true);
      vector < LorentzVector > mybjets;
      vector < LorentzVector > myaddjets;
      for (unsigned int idx = 0; idx < jetIndexSortedCSV.size(); ++idx) {
        if (jets.ak4pfjets_passMEDbtag.at(jetIndexSortedCSV[idx]) == true) mybjets.push_back(jets.ak4pfjets_p4.at(jetIndexSortedCSV[idx]));
        else if (mybjets.size() <= 1 && (mybjets.size() + myaddjets.size()) < 3) myaddjets.push_back(jets.ak4pfjets_p4.at(jetIndexSortedCSV[idx]));
      }

      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
      vector < int > jetIndexSortedCSV_jup;
      if (sampleConfig.isfastsim) jetIndexSortedCSV_jup = JetUtil::JetIndexCSVsorted(jets_jup.ak4pfjets_CSV, jets_jup.ak4pfjets_p4, jets_jup.ak4pfjets_loose_pfid, jetConfig.ptcut, jetConfig.etacut, false);
      else jetIndexSortedCSV_jup = JetUtil::JetIndexCSVsorted(jets_jup.ak4pfjets_CSV, jets_jup.ak4pfjets_p4, jets_jup.ak4pfjets_loose_pfid, jetConfig.ptcut, jetConfig.etacut, true);
      vector < LorentzVector > mybjets_jup;
      vector < LorentzVector > myaddjets_jup;
      for (unsigned int idx = 0; idx < jetIndexSortedCSV_jup.size(); ++idx) {
        if (jets_jup.ak4pfjets_passMEDbtag.at(jetIndexSortedCSV_jup[idx]) == true) mybjets_jup.push_back(jets_jup.ak4pfjets_p4.at(jetIndexSortedCSV_jup[idx]));
        else if (mybjets_jup.size() <= 1 && (mybjets_jup.size() + myaddjets_jup.size()) < 3) myaddjets_jup.push_back(jets_jup.ak4pfjets_p4.at(jetIndexSortedCSV_jup[idx]));
      }

      if (debug) std::cout << "[babymaker::looper]: filling jet variables LINE:" << __LINE__ << std::endl;
      vector < int > jetIndexSortedCSV_jdown;
      if (sampleConfig.isfastsim) jetIndexSortedCSV_jdown = JetUtil::JetIndexCSVsorted(jets_jdown.ak4pfjets_CSV, jets_jdown.ak4pfjets_p4, jets_jdown.ak4pfjets_loose_pfid, jetConfig.ptcut, jetConfig.etacut, false);
      else jetIndexSortedCSV_jdown = JetUtil::JetIndexCSVsorted(jets_jdown.ak4pfjets_CSV, jets_jdown.ak4pfjets_p4, jets_jdown.ak4pfjets_loose_pfid, jetConfig.ptcut, jetConfig.etacut, true);
      vector < LorentzVector > mybjets_jdown;
      vector < LorentzVector > myaddjets_jdown;
      for (unsigned int idx = 0; idx < jetIndexSortedCSV_jdown.size(); ++idx) {
        if (jets_jdown.ak4pfjets_passMEDbtag.at(jetIndexSortedCSV_jdown[idx]) == true) mybjets_jdown.push_back(jets_jdown.ak4pfjets_p4.at(jetIndexSortedCSV_jdown[idx]));
        else if (mybjets_jdown.size() <= 1 && (mybjets_jdown.size() + myaddjets_jdown.size()) < 3) myaddjets_jdown.push_back(jets_jdown.ak4pfjets_p4.at(jetIndexSortedCSV_jdown[idx]));
      }
      // the following variables need jets to be calculated. 
      vector < float > dummy_sigma;
      dummy_sigma.clear(); //move outside of if-clause to be able to copy for photon selection
      for (size_t idx = 0; idx < jets.ak4pfjets_p4.size(); ++idx) {
        dummy_sigma.push_back(0.1);
      }

      if (jets.ak4pfjets_p4.size() > 1) {
        // Hadronic Chi2
        StopEvt.hadronic_top_chi2 = calculateChi2(jets.ak4pfjets_p4, dummy_sigma, jets.ak4pfjets_passMEDbtag);
        StopEvt.mbb = jets.ak4_mbb;
        StopEvt.mct = jets.ak4_mct;
        // Jets & MET
        StopEvt.mindphi_met_j1_j2 = getMinDphi(StopEvt.pfmet_phi, jets.ak4pfjets_p4.at(0), jets.ak4pfjets_p4.at(1));
        // DR(lep, leadB) with medium discriminator
        if (nVetoLeptons > 0) {
          StopEvt.dR_lep_leadb = dRbetweenVectors(jets.ak4pfjets_leadMEDbjet_p4, lep1.p4);
          StopEvt.MT2W = CalcMT2W_(mybjets, myaddjets, lep1.p4, StopEvt.pfmet, StopEvt.pfmet_phi);
          StopEvt.topness = CalcTopness_(0, StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, mybjets, myaddjets);
          StopEvt.MT2_lb_b_mass = CalcMT2_lb_b_(StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, mybjets, myaddjets, 0, true);
          StopEvt.MT2_lb_b = CalcMT2_lb_b_(StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, mybjets, myaddjets, 0, false);
          StopEvt.topnessMod = CalcTopness_(1, StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, mybjets, myaddjets);
          StopEvt.MT2_lb_bqq_mass = CalcMT2_lb_bqq_(StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, true);
          StopEvt.MT2_lb_bqq = CalcMT2_lb_bqq_(StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, false);
        }
        if (nVetoLeptons > 1) {
          StopEvt.dR_lep2_leadb = dRbetweenVectors(jets.ak4pfjets_leadMEDbjet_p4, lep2.p4);
          StopEvt.MT2W_lep2 = CalcMT2W_(mybjets, myaddjets, lep2.p4, StopEvt.pfmet, StopEvt.pfmet_phi);
          StopEvt.topness_lep2 = CalcTopness_(0, StopEvt.pfmet, StopEvt.pfmet_phi, lep2.p4, mybjets, myaddjets);
          StopEvt.topnessMod_lep2 = CalcTopness_(1, StopEvt.pfmet, StopEvt.pfmet_phi, lep2.p4, mybjets, myaddjets);
          StopEvt.MT2_lb_b_mass_lep2 = CalcMT2_lb_b_(StopEvt.pfmet, StopEvt.pfmet_phi, lep2.p4, mybjets, myaddjets, 0, true);
          StopEvt.MT2_lb_b_lep2 = CalcMT2_lb_b_(StopEvt.pfmet, StopEvt.pfmet_phi, lep2.p4, mybjets, myaddjets, 0, false);
          StopEvt.MT2_lb_bqq_mass_lep2 = CalcMT2_lb_bqq_(StopEvt.pfmet, StopEvt.pfmet_phi, lep2.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, true);
          StopEvt.MT2_lb_bqq_lep2 = CalcMT2_lb_bqq_(StopEvt.pfmet, StopEvt.pfmet_phi, lep2.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, false);
          StopEvt.MT2_l_l = CalcMT2_(StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, lep2.p4, false, 0);
        }
      }
      if (jets_jup.ak4pfjets_p4.size() > 1) {
        StopEvt.mindphi_met_j1_j2_jup = getMinDphi(StopEvt.pfmet_phi_jup, jets_jup.ak4pfjets_p4.at(0), jets_jup.ak4pfjets_p4.at(1));
        StopEvt.mbb_jup = jets_jup.ak4_mbb;
        StopEvt.mct_jup = jets_jup.ak4_mct;
        if (nVetoLeptons > 0) StopEvt.MT2W_jup = CalcMT2W_(mybjets_jup, myaddjets_jup, lep1.p4, StopEvt.pfmet_jup, StopEvt.pfmet_phi_jup);
        if (nVetoLeptons > 0) StopEvt.topnessMod_jup = CalcTopness_(1, StopEvt.pfmet_jup, StopEvt.pfmet_phi_jup, lep1.p4, mybjets_jup, myaddjets_jup);
      }
      if (jets_jdown.ak4pfjets_p4.size() > 1) {
        StopEvt.mindphi_met_j1_j2_jdown = getMinDphi(StopEvt.pfmet_phi_jdown, jets_jdown.ak4pfjets_p4.at(0), jets_jdown.ak4pfjets_p4.at(1));
        StopEvt.mbb_jdown = jets_jdown.ak4_mbb;
        StopEvt.mct_jdown = jets_jdown.ak4_mct;
        if (nVetoLeptons > 0) StopEvt.MT2W_jdown = CalcMT2W_(mybjets_jdown, myaddjets_jdown, lep1.p4, StopEvt.pfmet_jdown, StopEvt.pfmet_phi_jdown);
        if (nVetoLeptons > 0) StopEvt.topnessMod_jdown = CalcTopness_(1, StopEvt.pfmet_jdown, StopEvt.pfmet_phi_jdown, lep1.p4, mybjets_jdown, myaddjets_jdown);
      }

      vector < pair < float, int > > rankminDR;
      vector < pair < float, int > > rankmaxDPhi;
      vector < pair < float, int > > rankminDR_lep2;
      vector < pair < float, int > > rankmaxDPhi_lep2;
      for (unsigned int idx = 0; idx < jets.ak4pfjets_p4.size(); ++idx) {
        if (nVetoLeptons == 0) continue;
        pair < float, int > mypair;
        mypair.second = idx;
        mypair.first = getdphi(jets.ak4pfjets_p4.at(idx).Phi(), lep1.p4.Phi());
        rankmaxDPhi.push_back(mypair);
        mypair.first = dRbetweenVectors(jets.ak4pfjets_p4.at(idx), lep1.p4);
        rankminDR.push_back(mypair);
        if (nVetoLeptons <= 1) continue;
        mypair.first = getdphi(jets.ak4pfjets_p4.at(idx).Phi(), lep2.p4.Phi());
        rankmaxDPhi_lep2.push_back(mypair);
        mypair.first = dRbetweenVectors(jets.ak4pfjets_p4.at(idx), lep2.p4);
        rankminDR_lep2.push_back(mypair);
      }
      sort(rankminDR.begin(), rankminDR.end(), CompareIndexValueSmallest);
      sort(rankminDR_lep2.begin(), rankminDR_lep2.end(), CompareIndexValueSmallest);
      sort(rankmaxDPhi.begin(), rankmaxDPhi.end(), CompareIndexValueGreatest);
      sort(rankmaxDPhi_lep2.begin(), rankmaxDPhi_lep2.end(), CompareIndexValueGreatest);

      if (jets.ak4pfjets_p4.size() > 0) {
        for (unsigned int idx = 0; idx < rankminDR.size(); ++idx) {
          if (nVetoLeptons == 0) continue;
          if (!(jets.ak4pfjets_passMEDbtag.at(rankminDR[idx].second))) continue;
          StopEvt.Mlb_closestb = (jets.ak4pfjets_p4.at(rankminDR[idx].second) + lep1.p4).M();
          break;
        }
        for (unsigned int idx = 0; idx < rankminDR_lep2.size(); ++idx) {
          if (nVetoLeptons <= 1) continue;
          if (!(jets.ak4pfjets_passMEDbtag.at(rankminDR_lep2[idx].second))) continue;
          StopEvt.Mlb_closestb_lep2 = (jets.ak4pfjets_p4.at(rankminDR_lep2[idx].second) + lep2.p4).M();
          break;
        }
        if (nVetoLeptons > 0) StopEvt.Mlb_lead_bdiscr = (jets.ak4pfjets_p4.at(jetIndexSortedCSV[0]) + lep1.p4).M();
        if (nVetoLeptons > 1) StopEvt.Mlb_lead_bdiscr_lep2 = (jets.ak4pfjets_p4.at(jetIndexSortedCSV[0]) + lep2.p4).M();
        if (rankmaxDPhi.size() >= 3) {
          StopEvt.Mjjj = (jets.ak4pfjets_p4.at(rankmaxDPhi[0].second) + jets.ak4pfjets_p4.at(rankmaxDPhi[1].second) + jets.ak4pfjets_p4.at(rankmaxDPhi[2].second)).M();
        }
        if (rankmaxDPhi_lep2.size() >= 3) {
          StopEvt.Mjjj_lep2 = (jets.ak4pfjets_p4.at(rankmaxDPhi_lep2[0].second) + jets.ak4pfjets_p4.at(rankmaxDPhi_lep2[1].second) + jets.ak4pfjets_p4.at(rankmaxDPhi_lep2[2].second)).M();
        }

      } // end if >0 jets      

      if (fillextraConfig.fillZll) {
        //
        // Zll Event Variables
        //
        //  first find a Zll
        //  fill only for 2 or three lepton events//Zl2 will have always idx 1(2) for 2l(3l) events
        //  if four lepton events test only leading three leptons
        //  Zll must be always OS, then prefer OF, then prefer Zmass
        //  Zll needs to go before ph, as we recalculate myaddjets, mybjets
        if (nLooseLeptons >= 2) {
          int Zl1 = -1;
          if (nLooseLeptons == 2 && AllLeps[0].id * AllLeps[1].id < 0) Zl1 = 0;
          else if (nLooseLeptons >= 3) {
            if (nGoodLeptons == 1 && AllLeps[1].id * AllLeps[2].id < 0) Zl1 = 1; //0 is taken by selection lepton
            else if (nGoodLeptons >= 2) { //selection lepton can be either 0(lep1) or 1(lep2)
              if ((AllLeps[0].id * AllLeps[2].id) < 0 && (AllLeps[1].id * AllLeps[2].id) > 0) Zl1 = 0;
              else if ((AllLeps[0].id * AllLeps[2].id) > 0 && (AllLeps[1].id * AllLeps[2].id) < 0) Zl1 = 1;
              else if ((AllLeps[0].id * AllLeps[2].id) < 0 && (AllLeps[1].id * AllLeps[2].id) < 0) {
                if ((abs(AllLeps[0].id) == abs(AllLeps[2].id)) && !(abs(AllLeps[1].id) == abs(AllLeps[2].id))) Zl1 = 0;
                else if (!(abs(AllLeps[0].id) == abs(AllLeps[2].id)) && (abs(AllLeps[1].id) == abs(AllLeps[2].id))) Zl1 = 1;
                else { //Z will be SF/OF for either combination, decide on Zmass
                  if (fabs((AllLeps[0].p4 + AllLeps[2].p4).M() - 91.) > fabs((AllLeps[1].p4 + AllLeps[2].p4).M() - 91.)) Zl1 = 1;
                  else Zl1 = 0;
                }
              }
            }
          } //end of Zl1 selection
          if (Zl1 >= 0) { //found a Zll candidate
            int Zl2 = -1;
            if (nLooseLeptons == 2) Zl2 = 1;
            else Zl2 = 2;
            StopEvt.Zll_idl1 = AllLeps[Zl1].id;
            StopEvt.Zll_idl2 = AllLeps[Zl2].id;
            StopEvt.Zll_p4l1 = AllLeps[Zl1].p4;
            StopEvt.Zll_p4l2 = AllLeps[Zl2].p4;
            StopEvt.Zll_OS = (AllLeps[Zl1].id * AllLeps[Zl2].id < 0);
            StopEvt.Zll_SF = (abs(AllLeps[Zl1].id) == abs(AllLeps[Zl2].id));
            StopEvt.Zll_isZmass = (fabs((AllLeps[Zl1].p4 + AllLeps[Zl2].p4).M() - 91.) < 15);
            StopEvt.Zll_M = (AllLeps[Zl1].p4 + AllLeps[Zl2].p4).M();
            StopEvt.Zll_p4 = (AllLeps[Zl1].p4 + AllLeps[Zl2].p4);
            if (nLooseLeptons == 2) StopEvt.Zll_selLep = 1; //selection lepton (not Zll) is lep1
            else if (Zl1 == 0) StopEvt.Zll_selLep = 2; //selection lepton (not Zll) is lep2
            else StopEvt.Zll_selLep = 1; //selection lepton (not Zll) is lep1
            //Zll_selLep decides if one has to use Mlb_closestb or Mlb_closestb_lep2 for variables without MET.
            //it also decides if we recalculate mt_met_lep using lep1 or lep2, etc.
            double Zllmetpx = (StopEvt.pfmet * cos(StopEvt.pfmet_phi)) + (StopEvt.Zll_p4.Px());
            double Zllmetpy = (StopEvt.pfmet * sin(StopEvt.pfmet_phi)) + (StopEvt.Zll_p4.Py());
            StopEvt.Zll_met = sqrt(pow(Zllmetpx, 2) + pow(Zllmetpy, 2));
            StopEvt.Zll_met_phi = atan2(Zllmetpy, Zllmetpx);
            if (jets.ak4pfjets_p4.size() > 1) StopEvt.Zll_mindphi_met_j1_j2 = getMinDphi(StopEvt.Zll_met_phi, jets.ak4pfjets_p4.at(0), jets.ak4pfjets_p4.at(1));
            if (nVetoLeptons > 2) StopEvt.Zll_MT2_l_l = CalcMT2_(StopEvt.pfmet, StopEvt.pfmet_phi, AllLeps[Zl1].p4, AllLeps[Zl2].p4, false, 0);
            if (StopEvt.Zll_selLep == 1) {
              StopEvt.Zll_mt_met_lep = calculateMt(lep1.p4, StopEvt.Zll_met, StopEvt.Zll_met_phi);
              StopEvt.Zll_dphi_Wlep = DPhi_W_lep(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep1.p4);
              if (jets.ak4pfjets_p4.size() > 1) {
                StopEvt.Zll_MT2W = CalcMT2W_(mybjets, myaddjets, lep1.p4, StopEvt.Zll_met, StopEvt.Zll_met_phi);
                StopEvt.Zll_topness = CalcTopness_(0, StopEvt.Zll_met, StopEvt.Zll_met_phi, lep1.p4, mybjets, myaddjets);
                StopEvt.Zll_topnessMod = CalcTopness_(1, StopEvt.Zll_met, StopEvt.Zll_met_phi, lep1.p4, mybjets, myaddjets);
                StopEvt.Zll_MT2_lb_b_mass = CalcMT2_lb_b_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep1.p4, mybjets, myaddjets, 0, true);
                StopEvt.Zll_MT2_lb_b = CalcMT2_lb_b_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep1.p4, mybjets, myaddjets, 0, false);
                StopEvt.Zll_MT2_lb_bqq_mass = CalcMT2_lb_bqq_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep1.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, true);
                StopEvt.Zll_MT2_lb_bqq = CalcMT2_lb_bqq_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep1.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, false);
              }
            } else {
              StopEvt.Zll_mt_met_lep = calculateMt(lep2.p4, StopEvt.Zll_met, StopEvt.Zll_met_phi);
              StopEvt.Zll_dphi_Wlep = DPhi_W_lep(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep2.p4);
              if (jets.ak4pfjets_p4.size() > 1) {
                StopEvt.Zll_MT2W = CalcMT2W_(mybjets, myaddjets, lep2.p4, StopEvt.Zll_met, StopEvt.Zll_met_phi);
                StopEvt.Zll_topness = CalcTopness_(0, StopEvt.Zll_met, StopEvt.Zll_met_phi, lep2.p4, mybjets, myaddjets);
                StopEvt.Zll_topnessMod = CalcTopness_(1, StopEvt.Zll_met, StopEvt.Zll_met_phi, lep2.p4, mybjets, myaddjets);
                StopEvt.Zll_MT2_lb_b_mass = CalcMT2_lb_b_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep2.p4, mybjets, myaddjets, 0, true);
                StopEvt.Zll_MT2_lb_b = CalcMT2_lb_b_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep2.p4, mybjets, myaddjets, 0, false);
                StopEvt.Zll_MT2_lb_bqq_mass = CalcMT2_lb_bqq_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep2.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, true);
                StopEvt.Zll_MT2_lb_bqq = CalcMT2_lb_bqq_(StopEvt.Zll_met, StopEvt.Zll_met_phi, lep2.p4, mybjets, myaddjets, jets.ak4pfjets_p4, 0, false);
              }
            }
          } //end of Zll filling
        } //end of Zll
      }
      // Photon Event Variables
      if (fillextraConfig.fillPhoton) {
        if (StopEvt.ph_selectedidx >= 0) {
          int oljind = ph.overlapJetId.at(StopEvt.ph_selectedidx);
          double phmetpx = (StopEvt.pfmet * cos(StopEvt.pfmet_phi)) + (ph.p4.at(StopEvt.ph_selectedidx).Px());
          double phmetpy = (StopEvt.pfmet * sin(StopEvt.pfmet_phi)) + (ph.p4.at(StopEvt.ph_selectedidx).Py());
          StopEvt.ph_met = sqrt(pow(phmetpx, 2) + pow(phmetpy, 2));
          StopEvt.ph_met_phi = atan2(phmetpy, phmetpx);
          StopEvt.ph_ngoodjets = jets.ngoodjets;
          StopEvt.ph_ngoodbtags = jets.ngoodbtags;
          if (oljind >= 0) {
            StopEvt.ph_ngoodjets = jets.ngoodjets - 1;
            if (jets.ak4pfjets_passMEDbtag.at(oljind) == true) StopEvt.ph_ngoodbtags = jets.ngoodbtags - 1;
          }
          vector < LorentzVector > jetsp4_phcleaned;
          jetsp4_phcleaned.clear(); //used
          vector < float > jetsCSV_phcleaned;
          jetsCSV_phcleaned.clear();
          vector < bool > jetsbtag_phcleaned;
          jetsbtag_phcleaned.clear();
          vector < float > dummy_sigma_phcleaned;
          dummy_sigma_phcleaned.clear(); //used
          int leadbidx = -1;
          float htph(0), htssmph(0), htosmph(0);
          for (unsigned int idx = 0; idx < jets.ak4pfjets_p4.size(); ++idx) {
            if ((int) idx == oljind) continue;
            jetsp4_phcleaned.push_back(jets.ak4pfjets_p4.at(idx));
            jetsCSV_phcleaned.push_back(jets.ak4pfjets_CSV.at(idx));
            jetsbtag_phcleaned.push_back(jets.ak4pfjets_passMEDbtag.at(idx));
            dummy_sigma_phcleaned.push_back(dummy_sigma.at(idx));
            htph += jets.ak4pfjets_pt.at(idx);
            float dPhiM = getdphi(StopEvt.ph_met_phi, jets.ak4pfjets_phi.at(idx));
            if (dPhiM < (TMath::Pi() / 2)) htssmph += jets.ak4pfjets_pt.at(idx);
            else htosmph += jets.ak4pfjets_pt.at(idx);
            if (leadbidx == -1 && jets.ak4pfjets_passMEDbtag.at(idx)) leadbidx = idx; //leading bjet photon cleaned
          }
          StopEvt.ph_HT = htph;
          StopEvt.ph_htssm = htssmph;
          StopEvt.ph_htosm = htosmph;
          StopEvt.ph_htratiom = StopEvt.ph_htssm / (StopEvt.ph_htosm + StopEvt.ph_htssm);
          mybjets.clear();
          myaddjets.clear();
          for (unsigned int idx = 0; idx < jetIndexSortedCSV.size(); ++idx) {
            if (jetIndexSortedCSV[idx] == oljind) continue;
            if (jets.ak4pfjets_passMEDbtag.at(jetIndexSortedCSV[idx]) == true) mybjets.push_back(jets.ak4pfjets_p4.at(jetIndexSortedCSV[idx]));
            else if (mybjets.size() <= 1 && (mybjets.size() + myaddjets.size()) < 3) myaddjets.push_back(jets.ak4pfjets_p4.at(jetIndexSortedCSV[idx]));
          }
          if (nVetoLeptons > 0) {
            StopEvt.ph_mt_met_lep = calculateMt(lep1.p4, StopEvt.ph_met, StopEvt.ph_met_phi);
            StopEvt.ph_dphi_Wlep = DPhi_W_lep(StopEvt.ph_met, StopEvt.ph_met_phi, lep1.p4);
          }
          if (nVetoLeptons > 1) StopEvt.ph_MT2_l_l = CalcMT2_(StopEvt.pfmet, StopEvt.pfmet_phi, lep1.p4, lep2.p4, false, 0);
          if (jetsp4_phcleaned.size() > 1) {
            if (nVetoLeptons > 0) {
              StopEvt.ph_MT2W = CalcMT2W_(mybjets, myaddjets, lep1.p4, StopEvt.ph_met, StopEvt.ph_met_phi);
              StopEvt.ph_topness = CalcTopness_(0, StopEvt.ph_met, StopEvt.ph_met_phi, lep1.p4, mybjets, myaddjets);
              StopEvt.ph_topnessMod = CalcTopness_(1, StopEvt.ph_met, StopEvt.ph_met_phi, lep1.p4, mybjets, myaddjets);
              StopEvt.ph_MT2_lb_b_mass = CalcMT2_lb_b_(StopEvt.ph_met, StopEvt.ph_met_phi, lep1.p4, mybjets, myaddjets, 0, true);
              StopEvt.ph_MT2_lb_b = CalcMT2_lb_b_(StopEvt.ph_met, StopEvt.ph_met_phi, lep1.p4, mybjets, myaddjets, 0, false);
              StopEvt.ph_MT2_lb_bqq_mass = CalcMT2_lb_bqq_(StopEvt.ph_met, StopEvt.ph_met_phi, lep1.p4, mybjets, myaddjets, jetsp4_phcleaned, 0, true);
              StopEvt.ph_MT2_lb_bqq = CalcMT2_lb_bqq_(StopEvt.ph_met, StopEvt.ph_met_phi, lep1.p4, mybjets, myaddjets, jetsp4_phcleaned, 0, false);
            }
            StopEvt.ph_hadronic_top_chi2 = calculateChi2(jetsp4_phcleaned, dummy_sigma_phcleaned, jetsbtag_phcleaned);
            StopEvt.ph_mindphi_met_j1_j2 = getMinDphi(StopEvt.ph_met_phi, jetsp4_phcleaned.at(0), jetsp4_phcleaned.at(1));
          } //at least two jets
          if (jetsp4_phcleaned.size() > 0) {
            if (nVetoLeptons > 0) {
              StopEvt.ph_Mlb_lead_bdiscr = (jets.ak4pfjets_p4.at(jetIndexSortedCSV[0]) + lep1.p4).M();
              if (oljind == jetIndexSortedCSV[0]) StopEvt.ph_Mlb_lead_bdiscr = (jets.ak4pfjets_p4.at(jetIndexSortedCSV[1]) + lep1.p4).M(); //exists as index=0 doesn't count for jetsp4_phcleaned.size()
              if (leadbidx >= 0) StopEvt.ph_dR_lep_leadb = dRbetweenVectors(jets.ak4pfjets_p4.at(leadbidx), lep1.p4);
              for (unsigned int idx = 0; idx < rankminDR.size(); ++idx) {
                if (rankminDR[idx].second == oljind) continue;
                if (nVetoLeptons == 0) continue;
                if (!(jets.ak4pfjets_passMEDbtag.at(rankminDR[idx].second))) continue;
                StopEvt.ph_Mlb_closestb = (jets.ak4pfjets_p4.at(rankminDR[idx].second) + lep1.p4).M();
                break;
              }
              LorentzVector threejetsum;
              threejetsum.SetPxPyPzE(0., 0., 0., 0.);
              int threejetcounter = 0;
              for (unsigned int idx = 0; idx < rankmaxDPhi.size(); ++idx) {
                if (rankmaxDPhi[idx].second == oljind) continue;
                threejetsum = threejetsum + jets.ak4pfjets_p4.at(rankmaxDPhi[idx].second);
                ++threejetcounter;
                if (threejetcounter >= 3) break;
              }
              if (threejetcounter == 3) StopEvt.ph_Mjjj = threejetsum.M();
            } //at least one lepton
          } //at least one jet
        } //end of photon additions
      }
      //
      // Loop over Taus to check tau veto candiates.
      //
      if (debug) std::cout << "[babymaker::looper]: tau  vars. LINE:" << __LINE__ << std::endl;
      int vetotaus = 0;
      double tau_pt = 20.0;
      double tau_eta = 2.3;

      Taus.tau_IDnames = taus_pf_IDnames();

      for (unsigned int iTau = 0; iTau < taus_pf_p4().size(); iTau++) {

        if (taus_pf_p4().at(iTau).Pt() < tau_pt) continue;
        if (fabs(taus_pf_p4().at(iTau).Eta()) > tau_eta) continue;

        Taus.FillCommon(iTau, tau_pt, tau_eta);
        if (isVetoTau(iTau, lep1.p4, lep1.charge)) {
          Taus.tau_isVetoTau.push_back(true);
          vetotaus++;
        } else Taus.tau_isVetoTau.push_back(false);
      }

      if (vetotaus < 1) StopEvt.PassTauVeto = true;
      else StopEvt.PassTauVeto = false;
      Taus.ngoodtaus = vetotaus;
      //
      // Loop over IsoTracks (Charged pfLeptons and pfChargedHadrons)
      //
      if (debug) std::cout << "[babymaker::laooper]: filling isotrack vars  LINE:" << __LINE__ << std::endl;
      int vetotracks_v2(0), vetotracks_v3(0);
      if (debug) cout << "before highPtPFcands" << endl;

      for (unsigned int ipf = 0; ipf < pfcands_p4().size(); ipf++) {
        float cand_pt = cms3.pfcands_p4().at(ipf).pt();
        if (pfcands_charge().at(ipf) == 0) continue;
        if (cand_pt < 5) continue;
        if (fabs(pfcands_p4().at(ipf).eta()) > 2.4) continue;
        if (fabs(pfcands_dz().at(ipf)) > 0.1) continue;

        //remove everything that is within 0.1 of selected lead and subleading leptons
        if (nVetoLeptons > 0) {
          if (ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(ipf), lep1.p4) < 0.1) continue;
        }
        if (nVetoLeptons > 1) {
          if (ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(ipf), lep2.p4) < 0.1) continue;
        }
        Tracks.FillCommon(ipf);
        LorentzVector temp(-99.9, -99.9, -99.9, -99.9);
        // 13 TeV Track Isolation Configuration, pfLep and pfCH
        if (nVetoLeptons > 0) {
          if (isVetoTrack_v2(ipf, lep1.p4, lep1.charge)) {
            Tracks.isoTracks_isVetoTrack_v2.push_back(true);
            vetotracks_v2++;
          } else Tracks.isoTracks_isVetoTrack_v2.push_back(false);
          if (isVetoTrack_v3(ipf, lep1.p4, lep1.charge)) {
            Tracks.isoTracks_isVetoTrack_v3.push_back(true);
            vetotracks_v3++;
          } else Tracks.isoTracks_isVetoTrack_v3.push_back(false);
        } else {
          if (isVetoTrack_v2(ipf, temp, 0)) {
            Tracks.isoTracks_isVetoTrack_v2.push_back(true);
            vetotracks_v2++;
          } else Tracks.isoTracks_isVetoTrack_v2.push_back(false);
          if (isVetoTrack_v3(ipf, temp, 0)) {
            Tracks.isoTracks_isVetoTrack_v3.push_back(true);
            vetotracks_v3++;
          } else Tracks.isoTracks_isVetoTrack_v3.push_back(false);
        }
      } // end loop over pfCands
      if (vetotracks_v3 < 1) StopEvt.PassTrackVeto = true;
      else StopEvt.PassTrackVeto = false;
      if (lepConfig.dolepveto && !(StopEvt.nvetoleps != 1 && StopEvt.PassTrackVeto && StopEvt.PassTauVeto)) continue;
      nEvents_pass_skim_2ndlepVeto++;
      if (debug) std::cout << "[babymaker::looper]: updating geninfo for recoleptons LINE:" << __LINE__ << std::endl;
      // Check that we have the gen leptons matched to reco leptons
      int lep1_match_idx = -99;
      int lep2_match_idx = -99;

      double min_dr_lep1 = 999.9;
      int min_dr_lep1_idx = -99;

      double min_dr_lep2 = 999.9;
      int min_dr_lep2_idx = -99;

      if (!evt_isRealData()) {

        if (nVetoLeptons > 0) {
          for (int iGen = 0; iGen < (int) gen_leps.p4.size(); iGen++) {
            if (abs(gen_leps.id.at(iGen)) == abs(lep1.pdgid)) {
              double temp_dr = ROOT::Math::VectorUtil::DeltaR(gen_leps.p4.at(iGen), lep1.p4);
              if (temp_dr < matched_dr) {
                lep1_match_idx = gen_leps.genpsidx.at(iGen);
                break;
              } else if (temp_dr < min_dr_lep1) {
                min_dr_lep1 = temp_dr;
                min_dr_lep1_idx = gen_leps.genpsidx.at(iGen);
              }
            } // end if lep1 id matches genLep id
          } // end loop over gen leps
        }
        if (nVetoLeptons > 1) {
          for (int iGen = 0; iGen < (int) gen_leps.p4.size(); iGen++) {
            if (lep1_match_idx == iGen) continue;
            if (abs(gen_leps.id.at(iGen)) == abs(lep2.pdgid)) {
              double temp_dr = ROOT::Math::VectorUtil::DeltaR(gen_leps.p4.at(iGen), lep2.p4);
              if (temp_dr < matched_dr) {
                lep2_match_idx = gen_leps.genpsidx.at(iGen);
                break;
              } else if (temp_dr < min_dr_lep2) {
                min_dr_lep2 = temp_dr;
                min_dr_lep2_idx = gen_leps.genpsidx.at(iGen);
              }
            }
          }
        }
        // If lep1 isn't matched to a lepton already stored, then try to find another match
        if ((nVetoLeptons > 0 && lep1_match_idx < 0)) {
          for (unsigned int genx = 0; genx < genps_p4().size(); genx++) {
            if (!genps_isLastCopy().at(genx)) continue;
            if (abs(genps_id().at(genx)) == abs(lep1.pdgid)) {
              double temp_dr = ROOT::Math::VectorUtil::DeltaR(genps_p4().at(genx), lep1.p4);
              if (temp_dr < matched_dr) {
                lep1_match_idx = genx;
                gen_leps.FillCommon(genx);
                break;
              } else if (temp_dr < min_dr_lep1) {
                min_dr_lep1 = temp_dr;
                min_dr_lep1_idx = genx;
              }
            }
          }
          // if lep1 is still unmatched, fill with closest match, if possible
          if (lep1_match_idx < 0 && min_dr_lep1_idx > 0) gen_leps.FillCommon(min_dr_lep1_idx);
        }

        // If lep2 isn't matched to a lepton already stored, then try to find another match
        if ((nVetoLeptons > 1 && lep2_match_idx < 0)) {
          for (unsigned int genx = 0; genx < genps_p4().size(); genx++) {
            if (!genps_isLastCopy().at(genx)) continue;
            if (abs(genps_id().at(genx)) == abs(lep2.pdgid)) {
              double temp_dr = ROOT::Math::VectorUtil::DeltaR(genps_p4().at(genx), lep2.p4);
              if (temp_dr < matched_dr) {
                lep2_match_idx = genx;
                gen_leps.FillCommon(genx);
                break;
              } else if (temp_dr < min_dr_lep2) {
                min_dr_lep2 = temp_dr;
                min_dr_lep2_idx = genx;
              }
            }
          }
          // if lep2 is still unmatched, fill with closest match, if possible
          if (lep2_match_idx < 0 && min_dr_lep2_idx > 0) gen_leps.FillCommon(min_dr_lep2_idx);
        }
      } // end if not data      

      //-----------------------------------------------------------------------------------//
      // calculate new met variable with the 2nd lepton removed
      float new_pfmet_x = 0;
      float new_pfmet_y = 0;
      // if fail the vetos, sum up second lepton, iso track, or tau 4 momentum
      if (nVetoLeptons > 1) {
        new_pfmet_x += lep2.p4.px();
        new_pfmet_y += lep2.p4.py();
      } else if (!StopEvt.PassTrackVeto) {
        vecLorentzVector isotrk_p4s;
        for (unsigned int idx = 0; idx < Tracks.isoTracks_isVetoTrack_v3.size(); idx++) {
          if (!Tracks.isoTracks_isVetoTrack_v3.at(idx)) continue;
          isotrk_p4s.push_back(Tracks.isoTracks_p4.at(idx));
        }
        if (isotrk_p4s.size() == 0) {
          std::cout << "ERROR!!! Event fails iso track veto but no isolated track found!!!" << std::endl;
        } else {
          std::sort(isotrk_p4s.begin(), isotrk_p4s.end(), sortP4byPt());
          new_pfmet_x += isotrk_p4s.at(0).px();
          new_pfmet_y += isotrk_p4s.at(0).py();
        }
      } else if (!StopEvt.PassTauVeto) {
        vecLorentzVector tau_p4s;
        for (unsigned int idx = 0; idx < Taus.tau_isVetoTau.size(); idx++) {
          if (!Taus.tau_isVetoTau.at(idx)) continue;
          tau_p4s.push_back(Taus.tau_p4.at(idx));
        }
        if (tau_p4s.size() == 0) {
          std::cout << "ERROR!!! Event fails tau veto but no tau found!!!" << std::endl;
        } else {
          std::sort(tau_p4s.begin(), tau_p4s.end(), sortP4byPt());
          new_pfmet_x += tau_p4s.at(0).px();
          new_pfmet_y += tau_p4s.at(0).py();
        }
      }
      // add the subtract back to the met.
      float new_pfmet_x_jup = new_pfmet_x + (StopEvt.pfmet_jup * std::cos(StopEvt.pfmet_phi_jup));
      float new_pfmet_y_jup = new_pfmet_y + (StopEvt.pfmet_jup * std::sin(StopEvt.pfmet_phi_jup));
      float new_pfmet_x_jdown = new_pfmet_x + (StopEvt.pfmet_jdown * std::cos(StopEvt.pfmet_phi_jdown));
      float new_pfmet_y_jdown = new_pfmet_y + (StopEvt.pfmet_jdown * std::sin(StopEvt.pfmet_phi_jdown));
      new_pfmet_x += StopEvt.pfmet * std::cos(StopEvt.pfmet_phi);
      new_pfmet_y += StopEvt.pfmet * std::sin(StopEvt.pfmet_phi);

      StopEvt.pfmet_rl = std::sqrt(pow(new_pfmet_x, 2) + pow(new_pfmet_y, 2));
      StopEvt.pfmet_phi_rl = std::atan2(new_pfmet_y, new_pfmet_x);
      StopEvt.pfmet_rl_jup = std::sqrt(pow(new_pfmet_x_jup, 2) + pow(new_pfmet_y_jup, 2));
      StopEvt.pfmet_phi_rl_jup = std::atan2(new_pfmet_y_jup, new_pfmet_x_jup);
      StopEvt.pfmet_rl_jdown = std::sqrt(pow(new_pfmet_x_jdown, 2) + pow(new_pfmet_y_jdown, 2));
      StopEvt.pfmet_phi_rl_jdown = std::atan2(new_pfmet_y_jdown, new_pfmet_x_jdown);
      //
      // calclate other quantities with new met.
      //
      if (debug) std::cout << "[babymaker::looper]: filling reclculated mt etc LINE:" << __LINE__ << std::endl;
      if (nVetoLeptons > 0) {
        StopEvt.mt_met_lep_rl = calculateMt(lep1.p4, StopEvt.pfmet_rl, StopEvt.pfmet_phi_rl);
        StopEvt.MT2W_rl = CalcMT2W_(mybjets, myaddjets, lep1.p4, StopEvt.pfmet_rl, StopEvt.pfmet_phi_rl);
        StopEvt.topnessMod_rl = CalcTopness_(1, StopEvt.pfmet_rl, StopEvt.pfmet_phi_rl, lep1.p4, mybjets, myaddjets);
        //JEC up
        StopEvt.mt_met_lep_rl_jup = calculateMt(lep1.p4, StopEvt.pfmet_rl_jup, StopEvt.pfmet_phi_rl_jup);
        StopEvt.MT2W_rl_jup = CalcMT2W_(mybjets_jup, myaddjets_jup, lep1.p4, StopEvt.pfmet_rl_jup, StopEvt.pfmet_phi_rl_jup);
        StopEvt.topnessMod_rl_jup = CalcTopness_(1, StopEvt.pfmet_rl_jup, StopEvt.pfmet_phi_rl_jup, lep1.p4, mybjets_jup, myaddjets_jup);
        //JEC down
        StopEvt.mt_met_lep_rl_jdown = calculateMt(lep1.p4, StopEvt.pfmet_rl_jdown, StopEvt.pfmet_phi_rl_jdown);
        StopEvt.MT2W_rl_jdown = CalcMT2W_(mybjets_jdown, myaddjets_jdown, lep1.p4, StopEvt.pfmet_rl_jdown, StopEvt.pfmet_phi_rl_jdown);
        StopEvt.topnessMod_rl_jdown = CalcTopness_(1, StopEvt.pfmet_rl_jdown, StopEvt.pfmet_phi_rl_jdown, lep1.p4, mybjets_jdown, myaddjets_jdown);
      }

      if (jets.ngoodjets > 1) {
        StopEvt.mindphi_met_j1_j2_rl = getMinDphi(StopEvt.pfmet_phi_rl, jets.ak4pfjets_p4.at(0), jets.ak4pfjets_p4.at(1));
        //  StopEvt.mindphi_met_j1_j2_rl_jup = getMinDphi(StopEvt.pfmet_phi_rl_jup,jets_jup.ak4pfjets_p4.at(0),jets_jup.ak4pfjets_p4.at(1));
        //  StopEvt.mindphi_met_j1_j2_rl_jdown = getMinDphi(StopEvt.pfmet_phi_rl_jdown,jets_jdown.ak4pfjets_p4.at(0),jets_jdown.ak4pfjets_p4.at(1));
      }
      if (!(StopEvt.pfmet >= skim_met) && !(StopEvt.pfmet_rl >= skim_met) && !(StopEvt.pfmet_rl_jup >= skim_met) && !(StopEvt.pfmet_rl_jdown >= skim_met) && !(StopEvt.pfmet_jup >= skim_met) && !(StopEvt.pfmet_jdown >= skim_met)) continue;
      nEvents_pass_skim_met++;
      //----------------------------------------------//
      // Trigger 				      //
      //----------------------------------------------//
      if (debug) std::cout << "[babymaker::looper]: filling HLT vars" << std::endl;
      //////////////// 2015 Run II  //////////////////////
      if (!sampleConfig.issignal) {
        StopEvt.HLT_MET = passHLTTriggerPattern("HLT_PFMET170_NoiseCleaned_v") || passHLTTriggerPattern("HLT_PFMET170_JetIdCleaned_v") || passHLTTriggerPattern("HLT_PFMET170_HBHECleaned_v") || passHLTTriggerPattern("HLT_PFMET170_NotCleaned_v");
        StopEvt.HLT_MET100_MHT100 = passHLTTriggerPattern("HLT_PFMET100_PFMHT100_IDTight_v");
        StopEvt.HLT_SingleEl = passHLTTriggerPattern("HLT_Ele27_eta2p1_WPLoose_Gsf_v") || passHLTTriggerPattern("HLT_Ele27_eta2p1_WPTight_Gsf_v") || passHLTTriggerPattern("HLT_Ele105_CaloIdVT_GsfTrkIdT_v") || passHLTTriggerPattern("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
        //StopEvt.HLT_SingleMu = passHLTTriggerPattern("HLT_IsoMu20_v") || passHLTTriggerPattern("HLT_IsoTkMu20_v") || passHLTTriggerPattern("HLT_IsoMu22_v") || passHLTTriggerPattern("HLT_IsoTkMu22_v") || passHLTTriggerPattern("HLT_IsoMu24_v") || passHLTTriggerPattern("HLT_IsoTkMu24_v");
        StopEvt.HLT_SingleMu = passHLTTriggerPattern("HLT_IsoMu24_v") || passHLTTriggerPattern("HLT_IsoTkMu24_v") || passHLTTriggerPattern("HLT_Mu50_v") || passHLTTriggerPattern("HLT_TkMu50_v") || passHLTTriggerPattern("HLT_Mu55_v");
        StopEvt.HLT_DiEl = passHLTTriggerPattern("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
        StopEvt.HLT_DiMu = passHLTTriggerPattern("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") || passHLTTriggerPattern("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
        StopEvt.HLT_MuE = passHLTTriggerPattern("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") || passHLTTriggerPattern("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
        //photons more complicated because of prescales
        StopEvt.HLT_Photon90_CaloIdL_PFHT500 = passHLTTriggerPattern("HLT_Photon90_CaloIdL_PFHT500_v");
        StopEvt.HLT_Photon22_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon22_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon30_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon30_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon36_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon36_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon50_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon50_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon75_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon75_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon90_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon90_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon120_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon120_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon165_R9Id90_HE10_IsoM = HLT_prescale(triggerName("HLT_Photon165_R9Id90_HE10_IsoM_v"));
        StopEvt.HLT_Photon175 = passHLTTriggerPattern("HLT_Photon175_v");
        StopEvt.HLT_Photon165_HE10 = passHLTTriggerPattern("HLT_Photon165_HE10_v");
      }
      BabyTree -> Fill(); // fill tree
    }
    file.Close(); //close event loop
  } //close file loop
  // Write and Close baby file
  BabyFile -> cd();
  BabyTree -> Write();
  counterhist -> Write();
  if (sampleConfig.issignal) {
    counterhistSig -> Write();
    histNEvts -> Write();
  }
  BabyFile -> Close();
  // clean up
  if (sampleConfig.issignal) {
    fxsec -> Close();
    delete fxsec;
  }
  if (!sampleConfig.isdata) {
    pileupfile -> Close();
    delete pileupfile;
  }
  if (jetConfig.dobtagsf) jets.deleteBtagSFTool();
  // Benchmarking
  bmark -> Stop("benchmark");
  cout << endl;
  cout << "Wrote babies into file " << BabyFile -> GetName() << endl;
  cout << "-----------------------------" << endl;
  cout << "Events Processed                     " << nEvents_processed << endl;
  cout << "Events with " << skim_nvtx << " Good Vertex            " << nEvents_pass_skim_nVtx << endl;
  cout << "Events with at least " << lepConfig.nlep << " Good Lepton   " << nEvents_pass_skim_nGoodLeps << endl;
  cout << "Events with at least " << jetConfig.njet << " Good Jets     " << nEvents_pass_skim_nGoodJets << endl;
  //  cout << "Events with at least " << skim_nBJets << " Good BJets   " << nEvents_pass_skim_nBJets << endl;
  cout << "Events passing 2nd Lep Veto " << lepConfig.dolepveto << "    " << nEvents_pass_skim_2ndlepVeto << endl;
  cout << "Events with MET > " << skim_met << " GeV             " << nEvents_pass_skim_met << endl;
  cout << "-----------------------------" << endl;
  cout << "CPU  Time:   " << Form("%.01f", bmark -> GetCpuTime("benchmark")) << endl;
  cout << "Real Time:   " << Form("%.01f", bmark -> GetRealTime("benchmark")) << endl;
  cout << endl;
  delete bmark;
  return 0;
}