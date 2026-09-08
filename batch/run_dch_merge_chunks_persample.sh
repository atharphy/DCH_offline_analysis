#!/bin/bash
CHUNKDIRBASE=$1
TARGETDIRBASE=$2
YEAR=$3
SAMPLE=$4
PREFIX=${5:-hist}

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/a/atahmad/updated_analysis/CMSSW_13_0_10/src
eval `scramv1 runtime -sh`

cd Updated_offline_framework

root -l -b -q -e ".L MergeChunkHists.C+" -e "MergeChunkHists(\"${CHUNKDIRBASE}/${YEAR}\",\"${TARGETDIRBASE}/${YEAR}\",{\"${PREFIX}_${SAMPLE}.root\"})"
