#include "EventTree.h"
#include "TString.h"
#include "CMS3.h"
#include "VertexSelections.h"
#include "StopSelections.h"
#include "MetSelections.h"

using namespace tas; 
EventTree::EventTree ()
{
}
 
EventTree::EventTree (const std::string &prefix)
    : prefix_(prefix)
{
}
 
void EventTree::FillCommon (const std::string &root_file_name)
{
    bool signal=false;

    std::string filestr (root_file_name);
    string signalstr ("mStop");
    string lspstr ("mLSP");//add those for testing purpose
    string smsstr ("SMS");//add those for testing purpose
    string tchwhstr("TChiWH");
    if (filestr.find(signalstr) != string::npos || filestr.find(lspstr)    != string::npos || filestr.find(smsstr)    != string::npos ||filestr.find(tchwhstr)    != string::npos) signal = true;

    run = evt_run();
    ls  = evt_lumiBlock();
    evt = evt_event();
    nvtxs = numberOfGoodVertices();
    firstGoodVtxIdx = firstGoodVertex();
 
    calomet = evt_calomet();
    calomet_phi = evt_calometPhi();

    is_data = evt_isRealData();

    // the recommended met filters //
    if(!signal){
      if(nvtxs>0) filt_met = true;
      else filt_met = false;
      filt_met = filt_met*filt_ecalTP()*filt_eeBadSc()*filt_hbheNoise()*filt_hbheNoiseIso()*filt_goodVertices()*filt_globalTightHalo2016();
      if (is_data) filt_badMuonFilter = badMuonFilterV2(); //still some problems with MC
      filt_badChargedCandidateFilter = badChargedCandidateFilterV2();
      filt_cscbeamhalo = filt_cscBeamHalo();
      filt_cscbeamhalo2015 = filt_cscBeamHalo2015();
      filt_globaltighthalo2016 = filt_globalTightHalo2016();
      filt_globalsupertighthalo2016 = filt_globalSuperTightHalo2016();
      filt_eebadsc = filt_eeBadSc();
      filt_goodvtx = filt_goodVertices(); //not working but same as our 1goodvertex requirement
      filt_ecallaser = filt_ecalLaser();
      filt_ecaltp = filt_ecalTP();
      filt_hcallaser = filt_hcalLaser();
      filt_trkfail = filt_trackingFailure();
      filt_trkPOG = filt_trkPOGFilters();
      filt_trkPOG_logerr_tmc = filt_trkPOG_logErrorTooManyClusters();
      filt_trkPOG_tmc =filt_trkPOG_manystripclus53X();
      filt_trkPOG_tms = filt_trkPOG_toomanystripclus53X();
      filt_hbhenoise = filt_hbheNoise(); // hbheNoiseFilter_25ns();
      filt_hbheisonoise = filt_hbheNoiseIso();//hbheIsoNoiseFilter();
//      filt_badmuons = filt_badMuons();
//      filt_duplicatemuons = filt_duplicateMuons();
//      filt_nobadmuons = filt_noBadMuons();
    }
    
    if (!is_data)
    {
      scale1fb = evt_scale1fb();
      xsec     = evt_xsec_incl();
      kfactor  = evt_kfactor();
//      pu_nvtxs = puInfo_nPUvertices().at(6);
      pu_ntrue = puInfo_trueNumInteractions().at(0);
      genweights = cms3.genweights();
      genweightsID = cms3.genweightsID();
       
      //if(signal){
      //	sparms_values = sparm_values();
      //	for ( auto name : sparm_names() )
      //	  sparms_names.push_back(name.Data());
      //          sparms_subProcessId     = sparm_subProcessId();	
      //}
      genmet = gen_met();
      genmet_phi = gen_metPhi();

      // calculate nupt
         LorentzVector nup4 (0.,0.,0.,0.);
         for(int iGen=0; iGen<(int)genps_p4().size(); iGen++){
              if(!(abs(genps_id().at(iGen) ) ==12 || abs(genps_id().at(iGen) ) == 14 || abs(genps_id().at(iGen) ) ==16)) continue;
                if(!genps_isHardProcess().at(iGen)) continue;
                nup4=nup4+genps_p4().at(iGen);
          }
          nupt = nup4.pt();
     
       //calculate genht
         float _genht=0;
         for (size_t gidx = 0; gidx < genps_p4().size(); gidx++){
           bool is_b_a_jet = true; //madgraph parameter called maxjetflavor that sets up to which id is the parton counted in this thing.  for 5f PDF (like the one used in WJets sample production), this is set to 5...therefore is_b_a_jet set to true
           if (genps_status().at(gidx) != 23) continue;
           int id = abs(genps_id().at(gidx));
           if ((id >= 1 && id<=4) || (id == 5 && is_b_a_jet) || (id == 21)) _genht += genps_p4().at(gidx).pt();
         }
         genht = _genht;

 
    }
    dataset = evt_dataset().at(0).Data();
    filename = root_file_name;
    cms3tag = evt_CMS3tag().at(0).Data();

    firstGoodVtxIdx = firstGoodVertex();   
 
}

void EventTree::Reset ()
{
    run = 0;
    ls = 0;
    evt = 0;

    ngoodleps      =  -9999;
    nlooseleps     =  -9999;
    nvetoleps      =  -9999;
    genlepsfromtop = -9999;

    genLepsHardProcess = -9999;
    genNusHardProcess  = -9999;
    
    is0lep    = -9999;
    is1lep    = -9999;
    is2lep    = -9999;
    isZtoNuNu = -9999;

    is1lepFromW   = -9999;
    is1lepFromTop = -9999;

    nvtxs 	=  -9999;
    firstGoodVtxIdx   = -9999;//not necessarily first vertex
    pfmet                = -9999.;
    pfmet_phi            = -9999.;
    pfmet_jup            = -9999.;
    pfmet_phi_jup        = -9999.;
    pfmet_jdown          = -9999.;
    pfmet_phi_jdown      = -9999.;
    pfmet_rl             = -9999.;
    pfmet_phi_rl         = -9999.;
    pfmet_rl_jdown             = -9999.;
    pfmet_phi_rl_jdown         = -9999.;
    pfmet_rl_jup             = -9999.;
    pfmet_phi_rl_jup         = -9999.;
    calomet              = -9999.;
    calomet_phi          = -9999.;
    scale1fb             = -9999.;
    xsec                 = -9999.;
    xsec_uncert          = -9999.;
    kfactor              = -9999.;
    pu_ntrue             = -9999.;
    dR_lep_leadb         = -9999.;
    dR_lep2_leadb        = -9999.;
    mbb                  = -9999.;
    mbb_jup              = -9999.;
    mbb_jdown            = -9999.;
    mct                  = -9999.;
    mct_jup              = -9999.;
    mct_jdown            = -9999.;
    MT2W                 = -9999.;
    MT2W_lep2            = -9999.;
    MT2W_rl              = -9999.;
    MT2W_jup             = -9999.;
    MT2W_jdown           = -9999.;
    mindphi_met_j1_j2    = -9999.;
    mindphi_met_j1_j2_rl = -9999.;
    mindphi_met_j1_j2_jup= -9999.;
    mindphi_met_j1_j2_jdown= -9999.;
    mt_met_lep           = -9999.;
    mt_met_lep2          = -9999.;
    mt_met_lep_rl        = -9999.;
    mt_met_lep_jup       = -9999.;
    mt_met_lep_jdown     = -9999.;
    hadronic_top_chi2    = -9999.; 
    is_data              = false;

    topness                = -9999.; 
    topness_lep2           = -9999.; 
    topnessMod             = -9999.; 
    topnessMod_lep2        = -9999.;
    topnessMod_rl          = -9999.;
    topnessMod_jup         = -9999.;
    topnessMod_jdown       = -9999.;
    MT2_lb_b               = -9999.; 
    MT2_lb_b_lep2          = -9999.; 
    MT2_lb_b_mass          = -9999.; 
    MT2_lb_b_mass_lep2     = -9999.; 
    MT2_lb_bqq             = -9999.; 
    MT2_lb_bqq_lep2        = -9999.; 
    MT2_lb_bqq_mass        = -9999.; 
    MT2_lb_bqq_mass_lep2   = -9999.; 
    Mlb_closestb           = -9999.; 
    Mlb_lead_bdiscr        = -9999.; 
    Mlb_closestb_lep2      = -9999.; 
    Mlb_lead_bdiscr_lep2   = -9999.; 
    Mjjj                   = -9999.; 
    Mjjj_lep2              = -9999.; 

    dphi_Wlep       = -9999.;
    MET_over_sqrtHT = -9999.;
    ak4pfjets_rho   = -9999.;

    dataset  = "";
    filename = "";
    cms3tag  = "";
    pdf_up_weight = -9999.;
    pdf_down_weight = -9999.;
    genweights.clear();
    genweightsID.clear();

    weight_btagsf = -9999;
    weight_btagsf_heavy_UP = -9999;
    weight_btagsf_light_UP = -9999;
    weight_btagsf_heavy_DN = -9999;
    weight_btagsf_light_DN = -9999;
    weight_btagsf_fastsim_UP = -9999;
    weight_btagsf_fastsim_DN = -9999;
    weight_lepSF      = -9999;
    weight_lepSF_up   = -9999;
    weight_lepSF_down = -9999;
    weight_vetoLepSF      = -9999;
    weight_vetoLepSF_up   = -9999;
    weight_vetoLepSF_down = -9999;
    weight_lepSF_fastSim      = -9999;
    weight_lepSF_fastSim_up   = -9999;
    weight_lepSF_fastSim_down = -9999;
    weight_ISR     = -9999;
    weight_ISRup   = -9999;
    weight_ISRdown = -9999;
    weight_PU      = -9999;
    weight_PUup    = -9999;
    weight_PUdown  = -9999;
    hardgenpt      = -9999.;
    weight_ISRnjets     = -9999;
    weight_ISRnjets_UP  = -9999;
    weight_ISRnjets_DN  = -9999;
    NISRjets            = -9999;
    NnonISRjets         = -9999;
    filt_fastsimjets    = false;
    filt_fastsimjets_jup   = false;
    filt_fastsimjets_jdown = false;
    
    sparms_names.clear();
    sparms_values.clear();
    sparms_subProcessId 	= -9999;
    mass_stop                   = -9999;
    mass_lsp                    = -9999;
    mass_chargino               = -9999;
    mass_gluino                 = -9999;

    genmet 	= -9999.;
    genmet_phi 	= -9999.;
    nupt = -99999.;
    genht = -9999.;
    PassTrackVeto    = false;
    PassTauVeto      = false;

    HLT_SingleMu           = -9999.; 
    HLT_SingleEl           = -9999.;
    HLT_MET                = -9999.;
    HLT_MET100_MHT100 = -9999.;
    HLT_DiEl               = -9999.;
    HLT_MuE                = -9999.;
    HLT_DiMu               = -9999.;
    HLT_SingleMuNoIso      = -9999.;
    HLT_SingleElNoIso      = -9999.;

    HLT_Photon90_CaloIdL_PFHT500 = -9999.;
    HLT_Photon22_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon30_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon36_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon50_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon75_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon90_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon120_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon165_R9Id90_HE10_IsoM = -9999.;
    HLT_Photon175 = -9999.;
    HLT_Photon165_HE10 = -9999.;

    filt_cscbeamhalo = false;
     filt_cscbeamhalo2015 = false;
     filt_globaltighthalo2016 = false;
     filt_globalsupertighthalo2016 = false;
     filt_ecallaser = false;
     filt_ecaltp = false;
     filt_eebadsc = false;
     filt_goodvtx = false;
     filt_badevents = false;
     filt_hbhenoise = false;
     filt_hbheisonoise = false;
     filt_hcallaser = false;
     filt_met = false;
     filt_trkfail = false;
     filt_trkPOG = false;
     filt_trkPOG_logerr_tmc = false;
     filt_trkPOG_tmc = false;
     filt_trkPOG_tms = false;
     filt_badMuonFilter = false;
     filt_badChargedCandidateFilter = false;
     filt_badMuonFilter = false;

     filt_jetWithBadMuon        = false;
     filt_jetWithBadMuon_jup    = false;
     filt_jetWithBadMuon_jdown  = false;
     filt_pfovercalomet         = false;
     filt_badmuons              = false;
     filt_duplicatemuons        = false;
     filt_nobadmuons            = false; 
    nPhotons             = -9999;
    ph_selectedidx       = -9999;
    ph_ngoodjets         = -9999;
    ph_ngoodbtags        = -9999;
    ph_met               = -9999.;
    ph_met_phi           = -9999.;
    ph_HT                = -9999.;
    ph_htssm             = -9999.;
    ph_htosm             = -9999.;
    ph_htratiom          = -9999.;
    ph_mt_met_lep        = -9999.;
    ph_dphi_Wlep         = -9999.;
    ph_MT2W              = -9999.;
    ph_topness           = -9999.;
    ph_topnessMod        = -9999.;
    ph_MT2_l_l           = -9999.;
    ph_MT2_lb_b_mass     = -9999.;
    ph_MT2_lb_b          = -9999.;
    ph_MT2_lb_bqq_mass   = -9999.;
    ph_MT2_lb_bqq        = -9999.;
    ph_hadronic_top_chi2 = -9999.;
    ph_mindphi_met_j1_j2 = -9999.;
    ph_Mlb_lead_bdiscr   = -9999.;
    ph_dR_lep_leadb      = -9999.;
    ph_Mlb_closestb      = -9999.;
    ph_Mjjj              = -9999.;

    Zll_idl1              = -9999;
    Zll_idl2              = -9999;
    Zll_p4l1              = LorentzVector(0,0, 0,0);
    Zll_p4l2              = LorentzVector(0,0, 0,0);
    Zll_OS                = false;
    Zll_SF                = false;
    Zll_isZmass           = false;
    Zll_M                 = -9999.;
    Zll_p4                = LorentzVector(0,0, 0,0);
    Zll_selLep            = -9999;
    Zll_met               = -9999.;
    Zll_met_phi           = -9999.;
    Zll_mindphi_met_j1_j2 = -9999.;
    Zll_mt_met_lep        = -9999.;
    Zll_dphi_Wlep         = -9999.;
    Zll_MT2W              = -9999.;
    Zll_topness           = -9999.;
    Zll_topnessMod        = -9999.;
    Zll_MT2_lb_b_mass     = -9999.;
    Zll_MT2_l_l           = -9999.;
    Zll_MT2_lb_b          = -9999.;
    Zll_MT2_lb_bqq_mass   = -9999.;
    Zll_MT2_lb_bqq        = -9999.;

}

void EventTree::SetBranches (TTree* tree)                                                                                          
{
    tree->Branch("run", &run);
    tree->Branch("ls", &ls);
    tree->Branch("evt", &evt);   
    tree->Branch("nvtxs", &nvtxs);
    tree->Branch("pu_nvtxs", &pu_nvtxs);
    tree->Branch("pfmet", &pfmet);
    tree->Branch("pfmet_phi", &pfmet_phi);
    tree->Branch("pfmet_jup", &pfmet_jup);
    tree->Branch("pfmet_phi_jup", &pfmet_phi_jup);
    tree->Branch("pfmet_jdown", &pfmet_jdown);
    tree->Branch("pfmet_phi_jdown", &pfmet_phi_jdown);
    tree->Branch("pfmet_rl", &pfmet_rl);
    tree->Branch("pfmet_phi_rl", &pfmet_phi_rl);
    tree->Branch("pfmet_rl_jup", &pfmet_rl_jup);
    tree->Branch("pfmet_phi_rl_jup", &pfmet_phi_rl_jup);
    tree->Branch("pfmet_rl_jdown", &pfmet_rl_jdown);
    tree->Branch("pfmet_phi_rl_jdown", &pfmet_phi_rl_jdown);
    tree->Branch("scale1fb", &scale1fb);
    tree->Branch("xsec", &xsec);
    tree->Branch("xsec_uncert", &xsec_uncert);
    tree->Branch("kfactor", &kfactor);
    tree->Branch("pu_ntrue", &pu_ntrue);    
    tree->Branch("ngoodleps",&ngoodleps);
    tree->Branch("nlooseleps",&nlooseleps);
    tree->Branch("nvetoleps",&nvetoleps);
    tree->Branch("is_data", &is_data);
    tree->Branch("dataset", &dataset);
    tree->Branch("filename", &filename);
    tree->Branch("filt_badmuons", &filt_badmuons);
    tree->Branch("filt_duplicatemuons", &filt_duplicatemuons);
    tree->Branch("filt_nobadmuons", &filt_nobadmuons);
    tree->Branch("cms3tag", &cms3tag);
    tree->Branch("nEvents", &nEvents);
    tree->Branch("nEvents_goodvtx", &nEvents_goodvtx);
    tree->Branch("nEvents_MET30", &nEvents_MET30);
    tree->Branch("nEvents_1goodlep", &nEvents_1goodlep);
    tree->Branch("nEvents_2goodjets", &nEvents_2goodjets);
    tree->Branch("is0lep", &is0lep);
    tree->Branch("is1lep", &is1lep);
    tree->Branch("is2lep", &is2lep);
    tree->Branch("isZtoNuNu", &isZtoNuNu);
    tree->Branch("is1lepFromW", &is1lepFromW);
    tree->Branch("is1lepFromTop", &is1lepFromTop);
    tree->Branch("mbb",&mbb);
    tree->Branch("mbb_jup",&mbb_jup);
    tree->Branch("mbb_jdown",&mbb_jdown);
    tree->Branch("mct",&mct);
    tree->Branch("mct_jup",&mct_jup);
    tree->Branch("mct_jdown",&mct_jdown);
    tree->Branch("MT2W",&MT2W);
    tree->Branch("MT2W_rl",&MT2W_rl); 
    tree->Branch("mindphi_met_j1_j2", &mindphi_met_j1_j2);
    tree->Branch("mindphi_met_j1_j2_rl", &mindphi_met_j1_j2_rl);
    tree->Branch("mt_met_lep", &mt_met_lep);
    tree->Branch("mt_met_lep_rl", &mt_met_lep_rl);
    tree->Branch("MT2W_jup",&MT2W_jup);
    tree->Branch("MT2W_rl_jup",&MT2W_rl_jup);
    tree->Branch("mindphi_met_j1_j2_jup", &mindphi_met_j1_j2_jup);
    tree->Branch("mindphi_met_j1_j2_rl_jup", &mindphi_met_j1_j2_rl_jup);
    tree->Branch("mt_met_lep_jup", &mt_met_lep_jup);
    tree->Branch("mt_met_lep_rl_jup", &mt_met_lep_rl_jup);
    tree->Branch("MT2W_jdown",&MT2W_jdown);
    tree->Branch("MT2W_rl_jdown",&MT2W_rl_jdown);
    tree->Branch("mindphi_met_j1_j2_jdown", &mindphi_met_j1_j2_jdown);
    tree->Branch("mindphi_met_j1_j2_rl_jdown", &mindphi_met_j1_j2_rl_jdown);
    tree->Branch("mt_met_lep_jdown", &mt_met_lep_jdown);
    tree->Branch("mt_met_lep_rl_jdown", &mt_met_lep_rl_jdown);
    tree->Branch("hadronic_top_chi2", &hadronic_top_chi2);
    tree->Branch("ak4pfjets_rho", &ak4pfjets_rho);
    tree->Branch("pdf_up_weight", &pdf_up_weight);
    tree->Branch("pdf_down_weight", &pdf_down_weight);
    tree->Branch("genweightsID", &genweightsID);
    tree->Branch("genweights", &genweights);
    tree->Branch("weight_btagsf", &weight_btagsf);
    tree->Branch("weight_btagsf_heavy_UP", &weight_btagsf_heavy_UP);
    tree->Branch("weight_btagsf_light_UP", &weight_btagsf_light_UP);
    tree->Branch("weight_btagsf_heavy_DN", &weight_btagsf_heavy_DN);
    tree->Branch("weight_btagsf_light_DN", &weight_btagsf_light_DN);
    tree->Branch("weight_btagsf_fastsim_UP", &weight_btagsf_fastsim_UP);
    tree->Branch("weight_btagsf_fastsim_DN", &weight_btagsf_fastsim_DN);
    tree->Branch("weight_lepSF", &weight_lepSF);
    tree->Branch("weight_lepSF_up", &weight_lepSF_up);
    tree->Branch("weight_lepSF_down", &weight_lepSF_down);
    tree->Branch("weight_vetoLepSF", &weight_vetoLepSF);
    tree->Branch("weight_vetoLepSF_up", &weight_vetoLepSF_up);
    tree->Branch("weight_vetoLepSF_down", &weight_vetoLepSF_down);
    tree->Branch("weight_lepSF_fastSim", &weight_lepSF_fastSim);
    tree->Branch("weight_lepSF_fastSim_up", &weight_lepSF_fastSim_up);
    tree->Branch("weight_lepSF_fastSim_down", &weight_lepSF_fastSim_down);
    tree->Branch("weight_ISR", &weight_ISR);
    tree->Branch("weight_ISRup", &weight_ISRup);
    tree->Branch("weight_ISRdown", &weight_ISRdown);
    tree->Branch("weight_PU", &weight_PU);
    tree->Branch("weight_PUup", &weight_PUup);
    tree->Branch("weight_PUdown", &weight_PUdown);
    tree->Branch("weight_ISRnjets", &weight_ISRnjets);
    tree->Branch("weight_ISRnjets_UP", &weight_ISRnjets_UP);
    tree->Branch("weight_ISRnjets_DN", &weight_ISRnjets_DN);
    tree->Branch("NISRjets", &NISRjets);
    tree->Branch("NnonISRjets", &NnonISRjets);
    tree->Branch("filt_fastsimjets", &filt_fastsimjets);
    tree->Branch("filt_fastsimjets_jup", &filt_fastsimjets_jup);
    tree->Branch("filt_fastsimjets_jdown", &filt_fastsimjets_jdown);
    tree->Branch("filt_jetWithBadMuon", &filt_jetWithBadMuon);
    tree->Branch("filt_jetWithBadMuon_jup", &filt_jetWithBadMuon_jup);
    tree->Branch("filt_jetWithBadMuon_jdown", &filt_jetWithBadMuon_jdown);
    tree->Branch("filt_pfovercalomet", &filt_pfovercalomet);
    tree->Branch("sparms_names", &sparms_names);
    tree->Branch("sparms_values", &sparms_values);
    tree->Branch("sparms_subProcessId", &sparms_subProcessId);
    tree->Branch("mass_lsp", &mass_lsp);
    tree->Branch("mass_chargino", &mass_chargino);
    tree->Branch("mass_stop", &mass_stop);
    tree->Branch("mass_gluino", &mass_gluino);
    tree->Branch("genmet", &genmet);
    tree->Branch("genmet_phi", &genmet_phi);
    tree->Branch("nupt", &nupt);
    tree->Branch("genht", &genht);
    tree->Branch("PassTrackVeto",&PassTrackVeto);
    tree->Branch("PassTauVeto",&PassTauVeto);
    tree->Branch("topness", &topness); 
    tree->Branch("topnessMod", &topnessMod);
    tree->Branch("topnessMod_rl", &topnessMod_rl);  
    tree->Branch("topnessMod_jup", &topnessMod_jup);
    tree->Branch("topnessMod_rl_jup", &topnessMod_rl_jup);
    tree->Branch("topnessMod_jdown", &topnessMod_jdown);
    tree->Branch("topnessMod_rl_jdown", &topnessMod_rl_jdown);
    tree->Branch("Mlb_closestb", &Mlb_closestb); 
    tree->Branch("HLT_SingleEl", &HLT_SingleEl );
    tree->Branch("HLT_SingleMu", &HLT_SingleMu );
    tree->Branch("HLT_SingleMuNoIso", &HLT_SingleMuNoIso );
    tree->Branch("HLT_SingleElNoIso", &HLT_SingleElNoIso );
    tree->Branch("HLT_MET", &HLT_MET);
    tree->Branch("HLT_MET100_MHT100", &HLT_MET100_MHT100);
    tree->Branch("HLT_DiEl", &HLT_DiEl );
    tree->Branch("HLT_DiMu", &HLT_DiMu );
    tree->Branch("HLT_MuE", &HLT_MuE);
    tree->Branch("HLT_Photon90_CaloIdL_PFHT500", &HLT_Photon90_CaloIdL_PFHT500);
    tree->Branch("HLT_Photon22_R9Id90_HE10_IsoM", &HLT_Photon22_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon30_R9Id90_HE10_IsoM", &HLT_Photon30_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon36_R9Id90_HE10_IsoM", &HLT_Photon36_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon175", &HLT_Photon175);
    tree->Branch("HLT_Photon165_HE10", &HLT_Photon165_HE10);
    tree->Branch("nPhotons",             &nPhotons);
    tree->Branch("ph_ngoodjets",         &ph_ngoodjets);
    tree->Branch("ph_ngoodbtags",        &ph_ngoodbtags);
    tree->Branch("filt_met", &filt_met);
    tree->Branch("filt_badMuonFilter", &filt_badMuonFilter);
    tree->Branch("filt_badChargedCandidateFilter", &filt_badChargedCandidateFilter);
    tree->Branch("filt_cscbeamhalo2015", &filt_cscbeamhalo2015);
    tree->Branch("filt_cscbeamhalo", &filt_cscbeamhalo);
    tree->Branch("hardgenpt", &hardgenpt);
    tree->Branch("calomet", &calomet);
    tree->Branch("calomet_phi", &calomet_phi);
}

void EventTree::SetSecondLepBranches (TTree* tree)
{
    tree->Branch("MT2W_lep2",&MT2W_lep2);
    tree->Branch("mt_met_lep2", &mt_met_lep2);
    tree->Branch("dR_lep2_leadb", &dR_lep2_leadb);
    tree->Branch("topness_lep2", &topness_lep2); 
    tree->Branch("topnessMod_lep2", &topnessMod_lep2); 
    tree->Branch("Mlb_lead_bdiscr_lep2", &Mlb_lead_bdiscr_lep2);
    tree->Branch("Mlb_closestb_lep2", &Mlb_closestb_lep2);
    tree->Branch("MT2_lb_bqq_mass_lep2", &MT2_lb_bqq_mass_lep2);
    tree->Branch("MT2_lb_bqq_lep2", &MT2_lb_bqq_lep2);
    tree->Branch("MT2_lb_b_mass_lep2", &MT2_lb_b_mass_lep2); 
    tree->Branch("MT2_lb_b_lep2", &MT2_lb_b_lep2);
}

void EventTree::SetExtraVariablesBranches (TTree* tree)
{
    tree->Branch("dR_lep_leadb", &dR_lep_leadb);
    tree->Branch("dphi_Wlep", &dphi_Wlep);
    tree->Branch("MET_over_sqrtHT", &MET_over_sqrtHT);
    tree->Branch("MT2_l_l", &MT2_l_l); 
    tree->Branch("MT2_lb_b", &MT2_lb_b); 
    tree->Branch("MT2_lb_b_mass", &MT2_lb_b_mass); 
    tree->Branch("MT2_lb_bqq", &MT2_lb_bqq); 
    tree->Branch("MT2_lb_bqq_mass", &MT2_lb_bqq_mass);
    tree->Branch("Mlb_lead_bdiscr", &Mlb_lead_bdiscr);
    tree->Branch("Mjjj", &Mjjj);
    tree->Branch("Mjjj_lep2", &Mjjj_lep2);
 //   tree->Branch("calomet", &calomet);
 //   tree->Branch("calomet_phi", &calomet_phi);  
}

void EventTree::SetMETFilterBranches (TTree* tree)
{
    tree->Branch("filt_cscbeamhalo", &filt_cscbeamhalo);
    tree->Branch("filt_cscbeamhalo2015", &filt_cscbeamhalo2015);
    tree->Branch("filt_globaltighthalo2016", &filt_globaltighthalo2016);
    tree->Branch("filt_globalsupertighthalo2016", &filt_globalsupertighthalo2016);
    tree->Branch("filt_ecallaser", &filt_ecallaser);
    tree->Branch("filt_ecaltp", &filt_ecaltp);
    tree->Branch("filt_eebadsc", &filt_eebadsc);
    tree->Branch("filt_goodvtx", &filt_goodvtx);
    tree->Branch("filt_badevents", &filt_badevents);
    tree->Branch("filt_hbhenoise", &filt_hbhenoise);
    tree->Branch("filt_hbheisonoise", &filt_hbheisonoise);
    tree->Branch("filt_hcallaser", &filt_hcallaser);
    tree->Branch("filt_trkfail", &filt_trkfail);
    tree->Branch("filt_trkPOG", &filt_trkPOG);
    tree->Branch("filt_trkPOG_logerr_tmc", &filt_trkPOG_logerr_tmc);
    tree->Branch("filt_trkPOG_tmc", &filt_trkPOG_tmc);
    tree->Branch("filt_trkPOG_tms", &filt_trkPOG_tms);
    tree->Branch("firstGoodVtxIdx", &firstGoodVtxIdx);
//    tree->Branch("filt_badChargedCandidateFilter", &filt_badChargedCandidateFilter);
}

void EventTree::SetPhotonBranches (TTree* tree)
{
    tree->Branch("ph_selectedidx",       &ph_selectedidx);
    tree->Branch("ph_met",               &ph_met);
    tree->Branch("ph_met_phi",           &ph_met_phi);
    tree->Branch("ph_HT",                &ph_HT);
    tree->Branch("ph_htssm",             &ph_htssm);
    tree->Branch("ph_htosm",             &ph_htosm);
    tree->Branch("ph_htratiom",          &ph_htratiom);
    tree->Branch("ph_mt_met_lep",        &ph_mt_met_lep);
    tree->Branch("ph_dphi_Wlep",         &ph_dphi_Wlep);
    tree->Branch("ph_MT2W",              &ph_MT2W);
    tree->Branch("ph_topness",           &ph_topness);
    tree->Branch("ph_topnessMod",        &ph_topnessMod);
    tree->Branch("ph_MT2_l_l",           &ph_MT2_l_l);
    tree->Branch("ph_MT2_lb_b_mass",     &ph_MT2_lb_b_mass);
    tree->Branch("ph_MT2_lb_b",          &ph_MT2_lb_b);
    tree->Branch("ph_MT2_lb_bqq_mass",   &ph_MT2_lb_bqq_mass);
    tree->Branch("ph_MT2_lb_bqq",        &ph_MT2_lb_bqq);
    tree->Branch("ph_hadronic_top_chi2", &ph_hadronic_top_chi2);
    tree->Branch("ph_mindphi_met_j1_j2", &ph_mindphi_met_j1_j2);
    tree->Branch("ph_Mlb_lead_bdiscr",   &ph_Mlb_lead_bdiscr);
    tree->Branch("ph_dR_lep_leadb",      &ph_dR_lep_leadb);
    tree->Branch("ph_Mlb_closestb",      &ph_Mlb_closestb);
    tree->Branch("ph_Mjjj",              &ph_Mjjj);
}

void EventTree::SetZllBranches (TTree* tree)
{
    tree->Branch("Zll_idl1", &Zll_idl1);
    tree->Branch("Zll_idl2", &Zll_idl2);
    tree->Branch("Zll_p4l1", &Zll_p4l1);
    tree->Branch("Zll_p4l2", &Zll_p4l2);
    tree->Branch("Zll_OS", &Zll_OS);
    tree->Branch("Zll_SF", &Zll_SF);
    tree->Branch("Zll_isZmass", &Zll_isZmass);
    tree->Branch("Zll_M", &Zll_M);
    tree->Branch("Zll_p4", &Zll_p4);
    tree->Branch("Zll_selLep", &Zll_selLep);
    tree->Branch("Zll_met", &Zll_met);
    tree->Branch("Zll_met_phi", &Zll_met_phi);
    tree->Branch("Zll_mindphi_met_j1_j2", &Zll_mindphi_met_j1_j2);
    tree->Branch("Zll_mt_met_lep", &Zll_mt_met_lep);
    tree->Branch("Zll_dphi_Wlep", &Zll_dphi_Wlep);
    tree->Branch("Zll_MT2W", &Zll_MT2W);
    tree->Branch("Zll_topness", &Zll_topness);
    tree->Branch("Zll_topnessMod", &Zll_topnessMod);
    tree->Branch("Zll_MT2_l_l", &Zll_MT2_l_l);
    tree->Branch("Zll_MT2_lb_b_mass", &Zll_MT2_lb_b_mass);
    tree->Branch("Zll_MT2_lb_b", &Zll_MT2_lb_b);
    tree->Branch("Zll_MT2_lb_bqq_mass", &Zll_MT2_lb_bqq_mass);
    tree->Branch("Zll_MT2_lb_bqq", &Zll_MT2_lb_bqq);
    tree->Branch("genlepsfromtop", &genlepsfromtop);
}

/*/obsolete branches///
    tree->Branch("lep_sf", &lep_sf);
    tree->Branch("btag_sf", &btag_sf);
    tree->Branch("HLT_SingleEl_eff", &HLT_SingleEl_eff);
    tree->Branch("HLT_SingleMu_eff", &HLT_SingleMu_eff); 
    tree->Branch("HLT_MET170", &HLT_MET170 );
    tree->Branch("HLT_MET120Btag", &HLT_MET120Btag );
    tree->Branch("HLT_MET120Mu5", &HLT_MET120Mu5 );
    tree->Branch("HLT_HT350MET120", &HLT_HT350MET120 );
    tree->Branch("HLT_DiEl_17_12", &HLT_DiEl_17_12 );
    tree->Branch("HLT_Mu8El17", &HLT_Mu8El17 );
    tree->Branch("HLT_Mu8El23", &HLT_Mu8El23 );
    tree->Branch("HLT_Mu17El12", &HLT_Mu17El12 );
    tree->Branch("HLT_Mu23El12", &HLT_Mu23El12 );
    tree->Branch("HLT_SingleEl23", &HLT_SingleEl23);
    tree->Branch("HLT_SingleEl27", &HLT_SingleEl27 );
    tree->Branch("HLT_SingleEl27Tight", &HLT_SingleEl27Tight );
    tree->Branch("HLT_SingleElTight", &HLT_SingleElTight );
    tree->Branch("HLT_SingleElHT200", &HLT_SingleElHT200 );
    tree->Branch("HLT_SingleMuNoEta", &HLT_SingleMuNoEta );
    tree->Branch("HLT_SingleMuNoIso", &HLT_SingleMuNoIso );
    tree->Branch("HLT_SingleMuNoIsoNoEta", &HLT_SingleMuNoIsoNoEta );
    tree->Branch("HLT_Mu6HT200MET100", &HLT_Mu6HT200MET125 );
    tree->Branch("HLT_HT350MET100", &HLT_HT350MET100);
    tree->Branch("HLT_SingleMu17", &HLT_SingleMu17);
    tree->Branch("HLT_SingleMu18", &HLT_SingleMu18);
    tree->Branch("HLT_SingleMu20", &HLT_SingleMu20);
    tree->Branch("HLT_SingleMu24", &HLT_SingleMu24);
    tree->Branch("HLT_MonoCentPFJet80_METNoMu90_MHTNoMu90_IDTight", &HLT_MonoCentPFJet80_METNoMu90_MHTNoMu90_IDTight);
    tree->Branch("HLT_MET90_MHT90_IDTight", &HLT_MET90_MHT90_IDTight);
    tree->Branch("HLT_METNoMu90_NoiseCleaned_MHTNoMu90_IDTight", &HLT_METNoMu90_NoiseCleaned_MHTNoMu90_IDTight);
    tree->Branch("HLT_Photon90_CaloIdL_PFHT500", &HLT_Photon90_CaloIdL_PFHT500);
    tree->Branch("HLT_Photon22_R9Id90_HE10_IsoM", &HLT_Photon22_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon30_R9Id90_HE10_IsoM", &HLT_Photon30_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon36_R9Id90_HE10_IsoM", &HLT_Photon36_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM);
    tree->Branch("HLT_Photon175", &HLT_Photon175);
    tree->Branch("HLT_Photon165_HE10", &HLT_Photon165_HE10);

    tree->Branch("EA_all_rho", &EA_all_rho);   
    tree->Branch("EA_allcalo_rho", &EA_allcalo_rho); 
    tree->Branch("EA_centralcalo_rho", &EA_centralcalo_rho); 
    tree->Branch("EA_centralchargedpileup_rho", &EA_centralchargedpileup_rho); 
    tree->Branch("EA_centralneutral_rho", &EA_centralneutral_rho); 
    tree->Branch("filt_eff", &filt_eff);
    tree->Branch("nlooseleps",&nlooseleps);
    tree->Branch("genlepsfromtop", &genlepsfromtop);

//dont need to save these but do need them in the event tree
    tree->Branch("genLepsHardProcess", &genLepsHardProcess);
    tree->Branch("genNusHardProcess", &genNusHardProcess);
/////

    tree->Branch("sparms_filterEfficiency", &sparms_filterEfficiency);
    tree->Branch("sparms_pdfScale", &sparms_pdfScale);
    tree->Branch("sparms_pdfWeight1", &sparms_pdfWeight1);
    tree->Branch("sparms_pdfWeight2", &sparms_pdfWeight2);
    tree->Branch("sparms_weight", &sparms_weight);

    tree->Branch("firstGoodVtxIdx", &firstGoodVtxIdx);                                                                     
    tree->Branch("firstVtx_isfake", &firstVtx_isfake);
    tree->Branch("firstVtx_ndof", &firstVtx_ndof);
    tree->Branch("firstVtx_posRho", &firstVtx_posRho);
    tree->Branch("firstVtx_posZ", &firstVtx_posZ);
    tree->Branch("firstVtx_posp4", &firstVtx_posp4);
*/
