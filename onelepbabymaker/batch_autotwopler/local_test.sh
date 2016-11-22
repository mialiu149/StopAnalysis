#
# Environment
#
#if [ "$CMSSW_VER" == "CMSSW_8_0_5_patch1" ]; then CMSSW_VER=CMSSW_8_0_5; fi
FILENAME="/hadoop/cms/store/group/snt/run2_25ns_80MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/V08-00-05/merged_ntuple_1.root"
ANALYSIS=onelepbabymaker
BABY_TAG=${BABY_VERSION}
NEVENTS=1000
OUTPUT_NAMES=output.root,skim.root
IMERGED=$(echo $FILENAME | sed 's/.*merged_ntuple_\([0-9]\+\)\.root/\1/')

echo filename: $FILENAME
echo analysis: $ANALYSIS
echo baby_tag: $BABY_TAG
echo shortname: $SHORTNAME
echo nevents: $NEVENTS
echo output_names: $OUTPUT_NAMES
echo exe_args: $EXE_ARGS

cd $MAKER_PATH/NtupleTools/AutoTwopler/

echo "[wrapper] setting env"

setup_CMSSW_80X

echo "finished setting up env"

echo $PWD
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD

# Split Output files line into an array
OUTPUT_NAMES=(`echo $OUTPUT_NAMES | sed s/,/" "/g`)
echo "OUTPUT_NAMES=${OUTPUT_NAMES[*]}"

# Split Extra Args into an array 
EXE_ARGS=(`echo $EXE_ARGS | sed s/,/" "/g`)
echo "EXE_ARGS=${EXE_ARGS[*]}"

#
# Place for user code
#
cd package.$BABY_TAG/$ANALYSIS/

echo "ready to run:"
echo $PWD
ls

if [ ! -z ${EXE_ARGS[1]} ]; then
  ISFASTSIM=${EXE_ARGS[1]}
fi

source setupforcondor.sh    
echo $PWD
echo " Running BabyMaker:"
./runBabyMakerCondor ${FILENAME} ${NEVENTS}

# Format output for gfal transfer
mv ${SAMPLE_NAME}.root ${OUTPUT_NAMES[0]}
root -l -b -q skimBaby.C++'("'${OUTPUT_NAMES[0]}'", "'${OUTPUT_NAMES[1]}'")'
