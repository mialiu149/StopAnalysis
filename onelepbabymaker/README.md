#WH Analysis

Please check if CORE and Software from common TAS repository are up-to-date,
and that the directory matches the one in the makefile

to compile:
```
source compile.sh
```

runBabyMaker takes four arguments: ./runBabyMaker sample_name nevents file_name outpath sampledat

The sampledat contains the directory to the CMS3 ntuples.

Need to provide at least sample_name; nevents=-1 (-1=all events), file_name ,output=/nfs-7/userdata/stopRun2/, sampledat=sample.dat  by default

samplenames are hardcoded and can be found in sample.dat
