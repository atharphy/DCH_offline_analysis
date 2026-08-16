#!/bin/bash
YEAR=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/a/atahmad/updated_analysis/CMSSW_13_0_10/src
eval `scramv1 runtime -sh`
cd Offline_framework/offline_analysis

root -l -b -q "YieldPlots.C+(\"${YEAR}\")"
