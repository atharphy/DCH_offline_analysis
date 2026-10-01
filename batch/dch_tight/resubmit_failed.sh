#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")"

CHUNKLIST="chunklist.txt"
RETRYLIST="retry_chunklist.txt"
LOGDIR="logs"

if [ ! -f "$CHUNKLIST" ]; then
  echo "error: $CHUNKLIST not found in $(pwd)." >&2
  echo "Regenerate it from the repo root with:" >&2
  echo "  root -l -b -q batch/dch_tight/gen_chunklist.C" >&2
  exit 1
fi

: > "$RETRYLIST"
total=0
failed=0
while IFS=',' read -r year idx chunkidx; do
  [ -z "$year" ] && continue
  total=$((total + 1))
  logfile="$LOGDIR/dch_tight_chunk_${year}_${idx}_${chunkidx}.log"
  errfile="$LOGDIR/dch_tight_chunk_${year}_${idx}_${chunkidx}.err"
  ok=0
  if [ -f "$logfile" ] && grep -q "return value 0" "$logfile"; then
    if [ ! -s "$errfile" ] || ! grep -qiE "error|segmentation|abort|exception" "$errfile"; then
      ok=1
    fi
  fi
  if [ "$ok" -ne 1 ]; then
    echo "${year},${idx},${chunkidx}" >> "$RETRYLIST"
    failed=$((failed + 1))
  fi
done < "$CHUNKLIST"

echo "Checked ${total} dch_tight chunk jobs, ${failed} need resubmission -> ${RETRYLIST}"

if [ "${1:-}" = "--submit" ]; then
  if [ "$failed" -eq 0 ]; then
    echo "Nothing to resubmit."
  else
    condor_submit dch_tight_chunk_retry.sub
  fi
fi
