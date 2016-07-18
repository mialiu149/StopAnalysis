#!/bin/bash
function run () {
    echo root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/${TAG}_$1/\",\"${OUTPUTDIR}/$1.root\"\)
    nohup nice -n 19 root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/${TAG}_$1/\",\"$1\",\"${OUTPUTDIR}/$1.root\"\) >& merge_logs_${TAG}/log_merge_$1.txt &
}

source settings.sh 

mkdir -p $OUTPUTDIR
chmod -R a+wrx $OUTPUTDIR

if [ ! -d "merge_logs_${TAG}" ]; then
  mkdir merge_logs_${TAG}
fi
run  tbar_tch_4f_powheg_pythia8_25ns
run  t_sch_4f_amcnlo_pythia8_25ns
run  t_tbarW_5f_powheg_pythia8_25ns
run  t_tW_5f_powheg_pythia8_25ns
run  TTWJetsToQQ_amcnlo_pythia8_25ns
run  TTWJetsToLNu_amcnlo_pythia8_25ns
run  TTZToQQ_amcnlo_pythia8_25ns
run  TTZToLLNuNu_M-10_amcnlo_pythia8_25ns
run  ttbar_singleLeptFromT_madgraph_pythia8_25ns
run  ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns
run  ttbar_diLept_madgraph_pythia8_ext1_25ns
run  ttbar_powheg_pythia8_ext4_25ns
run  WJetsToLNu_madgraph_pythia8_25ns
run  WJetsToLNu_HT100To200_madgraph_pythia8_25ns
run  WJetsToLNu_HT100To200_madgraph_pythia8_ext1_25ns
run  WJetsToLNu_HT200To400_madgraph_pythia8_25ns
run  WJetsToLNu_HT200To400_madgraph_pythia8_ext1_25ns
run  WJetsToLNu_HT400To600_madgraph_pythia8_25ns
run  WJetsToLNu_HT600To800_madgraph_pythia8_25ns
run  WJetsToLNu_HT800To1200_madgraph_pythia8_25ns
run  WJetsToLNu_HT800To1200_madgraph_pythia8_ext1_25ns
run  WJetsToLNu_HT1200To2500_madgraph_pythia8_25ns
run  WJetsToLNu_HT2500ToInf_madgraph_pythia8_25ns
run  WJetsToLNu_HT600ToInf_madgraph_pythia8_25ns
run  DYJetsToLL_m10To50_amcnlo_pythia8_25ns
run  DYJetsToLL_m50_amcnlo_pythia8_25ns
run  WW_pythia8_25ns
run  WWToLNuQQ_powheg_25ns
run  WWTo2l2Nu_powheg_25ns
run  WZ_pythia8_25ns
run  WZTo3LNu_powheg_pythia8_25ns
run  WZTo2L2Q_amcnlo_pythia8_25ns
run  WZTo1L3Nu_amcnlo_pythia8_25ns
run  WZTo1LNu2Q_amcnlo_pythia8_25ns
run  ZZ_pythia8_25ns
run  ZZTo2L2Q_amcnlo_pythia8_25ns
run  ZZTo2Q2Nu_amcnlo_pythia8_25ns
run  ZZTo2L2Nu_powheg_pythia8_25ns
run  data_double_eg_Run2016B_MINIAOD_PromptReco-v2
run  data_double_mu_Run2016B_MINIAOD_PromptReco-v2
run  data_muon_eg_Run2016B_MINIAOD_PromptReco-v2
run  data_single_electron_Run2016B_MINIAOD_PromptReco-v2
run  data_single_muon_Run2016B_MINIAOD_PromptReco-v2
run  data_met_Run2016B_MINIAOD_PromptReco-v2
run  data_double_eg_Run2016C_MINIAOD_PromptReco-v2
run  data_double_mu_Run2016C_MINIAOD_PromptReco-v2
run  data_muon_eg_Run2016C_MINIAOD_PromptReco-v2
run  data_single_electron_Run2016C_MINIAOD_PromptReco-v2
run  data_single_muon_Run2016C_MINIAOD_PromptReco-v2
run  data_met_Run2016C_MINIAOD_PromptReco-v2
run  W1JetsToLNu_NuPt-200
run  W2JetsToLNu_NuPt-200
run  W3JetsToLNu_NuPt-200
run  W4JetsToLNu_NuPt-200
run  WplusH_HToBB_WToLNu
run  WminusH_HToBB_WToLNu-ext1
run  WZZ
run  WWW
run  ZZZ
run  WWZ
run  WWG



#run data_double_eg_Run2016B_MINIAOD_PromptReco-v2
#run data_double_mu_Run2016B_MINIAOD_PromptReco-v2
#run data_muon_eg_Run2016B_MINIAOD_PromptReco-v2
#run data_single_electron_Run2016B_MINIAOD_PromptReco-v2
#run data_single_muon_Run2016B_MINIAOD_PromptReco-v2
#run data_met_Run2016B_MINIAOD_PromptReco-v2
#run SMS_t2bw_scan
#run data_double_eg_Run2016C_MINIAOD_PromptReco-v2
#run data_double_mu_Run2016C_MINIAOD_PromptReco-v2
#run data_muon_eg_Run2016C_MINIAOD_PromptReco-v2
#run data_single_electron_Run2016C_MINIAOD_PromptReco-v2
#run data_single_muon_Run2016C_MINIAOD_PromptReco-v2
#run data_met_Run2016C_MINIAOD_PromptReco-v2
#run  W1JetsToLNu_NuPt-200
#run  W2JetsToLNu_NuPt-200
#run  W3JetsToLNu_NuPt-200
#run  W4JetsToLNu_NuPt-200
#run WplusH_HToBB_WToLNu
#run WminusH_HToBB_WToLNu-ext1
#run WminusH_HToBB_WToLNu
#run WZZ
#run WWW
#run ZZZ
#run WWZ
#run WWG 
