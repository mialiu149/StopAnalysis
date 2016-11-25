setup_CMSSW_80X
#! /bin/bash
#
# CONFIGURATION VARS
#
ANALYSIS_NAME=onelepbabymaker

BABY_VERSION=36.2.v7
#BABY_VERSION=test.v12

# do not modify this TARBALL_NAME
TARBALL_NAME=package.$BABY_VERSION

INSTRUCTIONS_FILE=instructions_tmp.txt
export INSTRUCTIONS_FILE=instructions_tmp.txt
EXECUTABLE_NAME=condor_executable.sh

BATCH_DIR=`pwd`

MAKER_NAME=$ANALYSIS_NAME
MAKER_PATH=/home/users/mliu/WHAnalysis/
MAKER_DIR=$MAKER_PATH/$MAKER_NAME/

CORE_NAME=CORE
CORE_PATH=/home/users/mliu/
CORE_DIR=$CORE_PATH/$CORE_NAME


cd $MAKER_PATH

#
# Checkout NtupleTools from git
#
if [ ! -d NtupleTools/ ]; then
    git clone ssh://git@github.com/cmstas/NtupleTools
fi

cd NtupleTools/AutoTwopler/

. setup.sh 

export basedir=`pwd`
if [ ! -e ~/public_html/.htaccess ]; then
    echo "[setup] don't have .htaccess file. copying one to ~/public_html/ for you"
    cp dashboard/htaccess ~/public_html/.htaccess
    chmod 755 ~/public_html/.htaccess
fi
export DASHBOARD=$(python -c "import utils; utils.make_dashboard()")
chmod -R 755 "$HOME/public_html/$DASHBASE/"
echo "[setup] dashboard is at: $DASHBOARD"
#sed -i 's/campaign = "80X"/campaign = "80X_stopBabyMaker"/g' params.py

cd ../../NtupleTools/condorMergingTools/libC/
make -j 20
cd -
#
# Copy over relevant scripts for batch
#
cp -r $BATCH_DIR/$INSTRUCTIONS_FILE ./
sed -i s/ANALYSIS/$ANALYSIS_NAME/g $INSTRUCTIONS_FILE
sed -i s/VERSION/$BABY_VERSION/g $INSTRUCTIONS_FILE
sed -i s/TARBALL/$TARBALL_NAME/g $INSTRUCTIONS_FILE
sed -i s/EXECUTABLE/$EXECUTABLE_NAME/g $INSTRUCTIONS_FILE

cp -r $BATCH_DIR/$EXECUTABLE_NAME ./

cp -r $BATCH_DIR/batch.py ./
cp -r ../condorMergingTools/libC/sweepRoot ./
cp -r $BATCH_DIR/sweepRoot.sh ./
cp -r $BATCH_DIR/mergeHadoopFiles.C ./
cp -r $BATCH_DIR/mergeScript.sh ./


#
# copy and compile tarball for babyMaker
#
if [ ! -d $TARBALL_NAME ]; then

    mkdir $TARBALL_NAME/
    mkdir $TARBALL_NAME/$CORE_NAME/
    mkdir $TARBALL_NAME/$MAKER_NAME/
    
    CONDOR_DIR_NAME=`pwd`/$TARBALL_NAME
    
    cp -r $MAKER_DIR/*.cc $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/*.h $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/Makefile $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/setupforcondor.sh $CONDOR_DIR_NAME/$MAKER_NAME/
    
    # modify makefile so it knows where CORE is on condor node, 
    #  relative to the untarred directory location.
    sed -i '/^COREPATH/d' $CONDOR_DIR_NAME/$MAKER_NAME/Makefile
    sed -i '1i COREPATH = ../CORE/' $CONDOR_DIR_NAME/$MAKER_NAME/Makefile
    
    cp -r $MAKER_DIR/xsec_susy_13tev.root $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/puWeights_2015data_2p2fbinv.root $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/*.dat $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/*.C $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/stop_variables/ $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/json_files/ $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/jecfiles/ $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/btagsf/ $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $MAKER_DIR/lepsf/ $CONDOR_DIR_NAME/$MAKER_NAME/

    cp -r $BATCH_DIR/sweepRoot $CONDOR_DIR_NAME/$MAKER_NAME/
    cp -r $BATCH_DIR/sweepRoot.sh $CONDOR_DIR_NAME/$MAKER_NAME/

    cp -r $CORE_DIR/* $CONDOR_DIR_NAME/$CORE_NAME/
    
    echo ""
    echo "  All files copied, compiling code..."
    echo ""
    
    cd $CONDOR_DIR_NAME/$CORE_NAME/
    make clean
    make -j 30
    
    cd ../$MAKER_NAME/stop_variables/
    make clean
    make -j 30
    
    cd ../
    make clean
    make -j 30

    echo ""
    echo "  Compiling code complete, creating tarball..."
    echo ""
    cd ../../
    tar czf $TARBALL_NAME.tar.gz $TARBALL_NAME/
fi

#
# Exectue NtupleTools/Autotwopler/setup.sh
#
#. setup.sh


#
# Move back to batch_autotwopler
#
cd $BATCH_DIR


#
# Print details about job submission to screen
#
echo "[setup] Merged Babies and Log Files Written to: "
echo "[setup]   /nfs-7/userdata/$USER/tupler_babies/"
echo "[setup] Unmerged Babies Written to: "
echo "[setup]   /hadoop/cms/store/user/$USER/AutoTwopler_babies/"
echo "[setup] To launch jobs do: "
echo "[setup]  . makeBabies.sh"


#
# Flag for successful environment setup
#
CONDOR_ENV_READY=true
