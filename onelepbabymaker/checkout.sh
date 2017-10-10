#!/bin/bash
setup_CMSSW_80X
#dirname="WHAnalysis_2017dev"
#git clone git@github.com:mialiu149/WHAnalysis.git ${dirname}
#cd ${dirname}
git clone git@github.com:cmstas/CORE.git
cd CORE
git checkout cms4
make -j 30
cd ../onelepbabymaker/
. setup.sh
