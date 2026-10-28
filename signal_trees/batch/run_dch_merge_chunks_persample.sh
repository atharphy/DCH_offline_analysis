#!/bin/bash
CHUNKDIRBASE=$1
TARGETDIRBASE=$2
YEAR=$3
SAMPLE=$4
PREFIX=${5:-hist}

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$REPO_DIR"
eval `scramv1 runtime -sh`

root -l -b -q -e ".L signal_trees/MergeChunkHists.C+" -e "MergeChunkHists(\"${CHUNKDIRBASE}/${YEAR}\",\"${TARGETDIRBASE}/${YEAR}\",{\"${PREFIX}_${SAMPLE}.root\"})"
