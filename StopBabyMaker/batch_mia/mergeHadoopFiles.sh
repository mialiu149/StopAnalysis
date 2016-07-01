#!/bin/bash
function run () {
    echo root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/${TAG}_$1/\",\"${OUTPUTDIR}/$1.root\"\)
    nohup nice -n 19 root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/${TAG}_$1/\",\"${OUTPUTDIR}/$1.root\"\) >& merge_logs_${TAG}/log_merge_$1.txt &
}

source settings.sh 

mkdir -p $OUTPUTDIR
chmod -R a+wrx $OUTPUTDIR

if [ ! -d "merge_logs_${TAG}" ]; then
  mkdir merge_logs_${TAG}
fi

#run data_double_eg_Run2016B_MINIAOD_PromptReco-v2
#run data_double_mu_Run2016B_MINIAOD_PromptReco-v2
#run data_muon_eg_Run2016B_MINIAOD_PromptReco-v2
#run data_single_electron_Run2016B_MINIAOD_PromptReco-v2
#run data_single_muon_Run2016B_MINIAOD_PromptReco-v2
#run data_met_Run2016B_MINIAOD_PromptReco-v2
run  W1JetsToLNu_NuPt-200
run  W2JetsToLNu_NuPt-200
run  W3JetsToLNu_NuPt-200
run  W4JetsToLNu_NuPt-200
run WplusH_HToBB_WToLNu
run WminusH_HToBB_WToLNu-ext1
run WminusH_HToBB_WToLNu
run WZZ
run WWW
run ZZZ
run WWZ
run WWG 
