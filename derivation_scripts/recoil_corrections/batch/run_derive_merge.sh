#!/bin/bash
YEAR=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/a/atahmad/updated_analysis/CMSSW_13_0_10/src
eval `scramv1 runtime -sh`

CORR_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-correctionlib/2.1.0-6dc02863165bc2126b8299c5b63785af/lib/python3.9/site-packages/correctionlib
export ROOT_INCLUDE_PATH="${CORR_BASE}/include:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="${CORR_BASE}/lib:${LD_LIBRARY_PATH:-}"

cd Updated_offline_framework

root -l -b -q -e '.L recoil_studies/DeriveRecoilCorrections.C+' -e "DeriveRecoilCorrectionsMerge(\"${YEAR}\")"
