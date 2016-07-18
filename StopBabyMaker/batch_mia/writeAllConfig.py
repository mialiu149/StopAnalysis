data = ['data_double_eg_Run2016B_MINIAOD_PromptReco-v2','data_double_mu_Run2016B_MINIAOD_PromptReco-v2','data_muon_eg_Run2016B_MINIAOD_PromptReco-v2','data_single_electron_Run2016B_MINIAOD_PromptReco-v2','data_single_muon_Run2016B_MINIAOD_PromptReco-v2','data_met_Run2016B_MINIAOD_PromptReco-v2']
dataC = ['data_double_eg_Run2016C_MINIAOD_PromptReco-v2','data_double_mu_Run2016C_MINIAOD_PromptReco-v2','data_muon_eg_Run2016C_MINIAOD_PromptReco-v2','data_single_electron_Run2016C_MINIAOD_PromptReco-v2','data_single_muon_Run2016C_MINIAOD_PromptReco-v2','data_met_Run2016C_MINIAOD_PromptReco-v2']
signal = ['SMS_t2bw_scan']
qcd = [
       'QCD_Pt-20to30_EMEnriched',
       'QCD_Pt-30to50_EMEnriched-ext1',
       'QCD_Pt-50to80_EMEnriched',
       'QCD_Pt-80to120_EMEnriched',
       'QCD_Pt-120to170_EMEnriched',
       'QCD_Pt-170to300_EMEnriched',
       'QCD_Pt-300toInf_EMEnriched',
       'QCD_Pt-30to50_MuEnrichedPt5',
       'QCD_Pt-50to80_MuEnrichedPt5',
       'QCD_Pt-80to120_MuEnrichedPt5',
       'QCD_Pt-120to170_MuEnrichedPt5',
       'QCD_Pt-20toInf_MuEnrichedPt5'
]
ttbar = [ 'ttbar_singleLeptFromT_madgraph_pythia8_25ns',
          #'ttbar_singleLeptFromTbar_madgraph_pythia8_25ns',
          'ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns',
          #'ttbar_diLept_madgraph_pythia8_25ns',
          'ttbar_diLept_madgraph_pythia8_ext1_25ns',
          #'ttbarTo2L2Nu_powheg_25ns',
          'ttbar_powheg_pythia8_ext4_25ns']
singletop = ['tbar_tch_4f_powheg_pythia8_25ns',
              't_sch_4f_amcnlo_pythia8_25ns',
              't_tbarW_5f_powheg_pythia8_25ns',
              't_tW_5f_powheg_pythia8_25ns']
wjets = [
         #'W1JetsToLNu_madgraph_pythia8_25ns',
         #'W2JetsToLNu_madgraph_pythia8_25ns',
         #'W3JetsToLNu_madgraph_pythia8_25ns',
         #'W4JetsToLNu_madgraph_pythia8_25ns',
         'WJetsToLNu_madgraph_pythia8_25ns',
         'WJetsToLNu_HT100To200_madgraph_pythia8_25ns',
         'WJetsToLNu_HT100To200_madgraph_pythia8_ext1_25ns',
         'WJetsToLNu_HT200To400_madgraph_pythia8_25ns',
         'WJetsToLNu_HT200To400_madgraph_pythia8_ext1_25ns',
         'WJetsToLNu_HT400To600_madgraph_pythia8_25ns',
         'WJetsToLNu_HT600To800_madgraph_pythia8_25ns',
         'WJetsToLNu_HT800To1200_madgraph_pythia8_25ns',
         'WJetsToLNu_HT800To1200_madgraph_pythia8_ext1_25ns',
         'WJetsToLNu_HT1200To2500_madgraph_pythia8_25ns',
         'WJetsToLNu_HT2500ToInf_madgraph_pythia8_25ns',
         'WJetsToLNu_HT600ToInf_madgraph_pythia8_25ns']
zjets = [
         'DYJetsToLL_m10To50_amcnlo_pythia8_25ns',
          'DYJetsToLL_m50_amcnlo_pythia8_25ns'
    #     'DYJetsToLL_M-50_HT-100to200',
    #     'DYJetsToLL_M-50_HT-200to400',
    #     'DYJetsToLL_M-50_HT-400to600',
   #      'DYJetsToLL_M-50_HT-600toInf'
         ]
diboson = ['WW_pythia8_25ns','WWToLNuQQ_powheg_25ns','WWTo2l2Nu_powheg_25ns','WZ_pythia8_25ns','WZTo3LNu_powheg_pythia8_25ns','WZTo2L2Q_amcnlo_pythia8_25ns','WZTo1L3Nu_amcnlo_pythia8_25ns','WZTo1LNu2Q_amcnlo_pythia8_25ns','ZZ_pythia8_25ns','ZZTo2L2Q_amcnlo_pythia8_25ns','ZZTo2Q2Nu_amcnlo_pythia8_25ns','ZZTo2L2Nu_powheg_pythia8_25ns']
addon = ['WplusH_HToBB_WToLNu','WminusH_HToBB_WToLNu-ext1','WZZ','WWW','ZZZ','WWZ','WWG']
wjets_highpt = ['W1JetsToLNu_NuPt-200','W2JetsToLNu_NuPt-200','W3JetsToLNu_NuPt-200','W4JetsToLNu_NuPt-200']
ttv = ['TTWJetsToQQ_amcnlo_pythia8_25ns',
       'TTWJetsToLNu_amcnlo_pythia8_25ns',
       'TTZToQQ_amcnlo_pythia8_25ns',
       'TTZToLLNuNu_M-10_amcnlo_pythia8_25ns']

towrite = open('towrite.sh','w')
f = open('../sample_2016.dat')
towrite.write('source settings.sh\n')
towrite.write('TAG=${analysis_version}\n')
for line in f:
    #for sample in data+dataC:
#    for sample in signal:
    for sample in singletop+ttv+ttbar+wjets+zjets+diboson+data+dataC+wjets_highpt+addon:
        print 'run ',sample
        if sample in line and 'Name' in line:
           print sample
           path = next(f).split()[1]
           print './writeConfig.sh '+path + ' ${TAG}_'+sample+' '+sample
           towrite.write('./writeConfig.sh '+path + ' ${TAG}_'+sample+' '+sample+'\n')
towrite.write('mkdir -p configs_${TAG}\n')
towrite.write('mv condor_${TAG}*.cmd configs_${TAG}\n')
towrite.write('echo "#!/bin/bash" > submitAll.sh\n')
towrite.write('echo "voms-proxy-init -voms cms -valid 240:00" >> submitAll.sh\nfor file in configs_${TAG}/*.cmd\ndo\n')
towrite.write('   echo "condor_submit ${file}" >> submitAll.sh\n')
towrite.write('done\n')
towrite.write('chmod +x submitAll.sh\n')
towrite.write('echo "[writeAllConfig] wrote submit script submitAll.sh"')
towrite.close()


