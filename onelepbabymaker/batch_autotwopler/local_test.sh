#
# Environment
#
#if [ "$CMSSW_VER" == "CMSSW_8_0_5_patch1" ]; then CMSSW_VER=CMSSW_8_0_5; fi

DATASET=/SingleElectron/Run2016B-PromptReco-v2/MINIAOD
FILENAME=output.root
ANALYSIS=onelepbabymaker
BABY_TAG=${BABY_VERSION}
SHORTNAME=$7
NEVENTS=1000
OUTPUT_NAMES=output.root,skim.root
EXE_ARGS=${10}
IMERGED=$(echo $FILENAME | sed 's/.*merged_ntuple_\([0-9]\+\)\.root/\1/')

echo dataset: $DATASET
echo filename: $FILENAME
echo analysis: $ANALYSIS
echo baby_tag: $BABY_TAG
echo shortname: $SHORTNAME
echo nevents: $NEVENTS
echo output_names: $OUTPUT_NAMES
echo exe_args: $EXE_ARGS

cd $MAKER_PATH/NtupleTools/AutoTwopler/

#python batch.py

echo "[wrapper] setting env"

setup_CMSSW_80X

echo "finished setting up env"

echo $PWD
# Need this for .so files
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

SAMPLE_NAME=data_single_electron_Run2016B_MINIAOD_PromptReco-v2re_test
ISFASTSIM=0
if [ ! -z ${EXE_ARGS[1]} ]; then
  ISFASTSIM=${EXE_ARGS[1]}
fi

source setupforcondor.sh    
echo $PWD
echo " Running BabyMaker:"
#echo "    ./runBabyMaker $SAMPLE_NAME $NEVENTS $IMERGED ./ sample_2016.dat $ISFASTSIM"
echo "    ./runBabyMaker $SAMPLE_NAME $NEVENTS"
#./runBabyMaker $SAMPLE_NAME $NEVENTS $IMERGED ./ sample_2016.dat $ISFASTSIM
./runBabyMaker $SAMPLE_NAME $NEVENTS 

# Format output for gfal transfer
mv ${SAMPLE_NAME}.root ${OUTPUT_NAMES[0]}
root -l -b -q skimBaby.C++'("'${OUTPUT_NAMES[0]}'", "'${OUTPUT_NAMES[1]}'")'
