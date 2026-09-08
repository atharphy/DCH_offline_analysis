#!/bin/bash
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/a/atahmad/updated_analysis/CMSSW_13_0_10/src/Updated_offline_framework
eval `scramv1 runtime -sh`

BASE=/eos/user/a/atahmad/DCH_offline_analysis/new_hists
LOGDIR=/tmp/merge_parallel_logs
mkdir -p "$LOGDIR"
rm -f "$LOGDIR"/*.log

for YEAR in 2016preVFP 2016postVFP 2017 2018; do
    CHUNKDIR="$BASE/chunks/run2_tauFR/$YEAR"
    TARGETDIR="$BASE/run2_tauFR/$YEAR"
    mapfile -t BASES < <(ls "$CHUNKDIR"/*.root 2>/dev/null | sed -E 's|.*/||; s/_chunk[0-9]+\.root$/.root/' | sort -u)
    N=${#BASES[@]}
    if [ "$N" -eq 0 ]; then echo "[$YEAR] nothing left to merge"; continue; fi
    HALF=$(( (N + 1) / 2 ))

    for PART in 0 1; do
        START=$(( PART * HALF ))
        SLICE=("${BASES[@]:$START:$HALF}")
        if [ "${#SLICE[@]}" -eq 0 ]; then continue; fi
        # Build a ROOT C++ vector<string> literal for this worker's slice.
        LIST=""
        for b in "${SLICE[@]}"; do LIST="${LIST}\"${b}\","; done
        LOG="$LOGDIR/${YEAR}_part${PART}.log"
        echo "[$YEAR part$PART] ${#SLICE[@]} sample(s) -> $LOG"
        root -l -b -q "MergeChunkHists.C(\"$CHUNKDIR\", \"$TARGETDIR\", {${LIST}})" > "$LOG" 2>&1 &
    done
done

wait
echo "ALL_PARALLEL_MERGES_DONE"
