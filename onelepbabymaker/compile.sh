source setup.sh
cd ../analysisutils/
make -j 30
cd ../onelepbabymaker
make clean
make -j 30
