#!/bin/bash
set -o pipefail
YEAR=$1
IDX=$2
BATCHSTART=$3
NUMINBATCH=$4

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$REPO_DIR"
eval `scramv1 runtime -sh`
export X509_USER_PROXY=${X509_USER_PROXY:-$HOME/.grid_proxy/x509up.pem}

CORR_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-correctionlib/2.1.0-6dc02863165bc2126b8299c5b63785af/lib/python3.9/site-packages/correctionlib
export ROOT_INCLUDE_PATH="${CORR_BASE}/include:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="${CORR_BASE}/lib:${LD_LIBRARY_PATH:-}"

export DCH_FILEMAP=qcdht
root -l -b -q -e ".L signal_trees/DCH_tight_TreeOnly_ChunkBatch.C+" \
    -e "DCH_tight_TreeOnly_ChunkBatch(\"${YEAR}\",${IDX},${BATCHSTART},${NUMINBATCH})" \
    2>&1 | tail -c 200000
