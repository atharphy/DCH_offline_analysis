#!/bin/bash
set -o pipefail
SCRIPT=$1
YEAR=$2
IDX=$3
CHUNKIDX=$4

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/a/atahmad/updated_analysis/CMSSW_13_0_10/src
eval `scramv1 runtime -sh`

CORR_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-correctionlib/2.1.0-6dc02863165bc2126b8299c5b63785af/lib/python3.9/site-packages/correctionlib
export ROOT_INCLUDE_PATH="${CORR_BASE}/include:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="${CORR_BASE}/lib:${LD_LIBRARY_PATH:-}"

cd Updated_offline_framework

# Hard cap on captured output -- keeps AFS log storage bounded regardless
# of how chatty a given run turns out to be (6550 jobs' worth of unbounded
# logs is how the AFS quota filled up last time).
#
# ${SCRIPT}_Chunk lives inside ${SCRIPT}.C for some drivers (DCH_tauFR.C)
# but in its own ${SCRIPT}_Chunk.C file for others (DCH_tight.C) -- load
# whichever file actually defines it, or every job fails immediately with
# "undeclared identifier" without ever touching an input file.
if [ -f "${SCRIPT}_Chunk.C" ]; then
    root -l -b -q -e ".L ${SCRIPT}_Chunk.C+" -e "${SCRIPT}_Chunk(\"${YEAR}\",${IDX},${CHUNKIDX})" 2>&1 | tail -c 200000
else
    root -l -b -q -e ".L ${SCRIPT}.C+" -e "${SCRIPT}_Chunk(\"${YEAR}\",${IDX},${CHUNKIDX})" 2>&1 | tail -c 200000
fi
