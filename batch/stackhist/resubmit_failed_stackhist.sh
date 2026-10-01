#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:-}"
SUBMIT_FLAG="${2:-}"
FILELIST="filelist_stackhist.txt"
RETRYLIST="run2_only.txt"
LOGDIR="logs"

case "$MODE" in
  stackhist)
    PREFIX="stackhist"
    RETRYSUB="stackhist_run2_retry.sub"
    ;;
  stackhist_multiplicity)
    PREFIX="stackhist_multiplicity"
    RETRYSUB="stackhist_multiplicity_run2_retry.sub"
    ;;
  *)
    echo "usage: $0 {stackhist|stackhist_multiplicity} [--submit]" >&2
    exit 1
    ;;
esac

if [ ! -f "$FILELIST" ]; then
  echo "error: $FILELIST not found in $(pwd)." >&2
  exit 1
fi

: > "$RETRYLIST"
total=0
failed=0
while IFS=',' read -r year first; do
  [ -z "$year" ] && continue
  total=$((total + 1))
  logfile="$LOGDIR/${PREFIX}_${year}_${first}.log"
  errfile="$LOGDIR/${PREFIX}_${year}_${first}.err"
  ok=0
  if [ -f "$logfile" ] && grep -q "return value 0" "$logfile"; then
    if [ ! -s "$errfile" ] || ! grep -qiE "error|segmentation|abort|exception" "$errfile"; then
      ok=1
    fi
  fi
  if [ "$ok" -ne 1 ]; then
    echo "${year},${first}" >> "$RETRYLIST"
    failed=$((failed + 1))
  fi
done < "$FILELIST"

echo "Checked ${total} ${MODE} jobs, ${failed} need resubmission -> ${RETRYLIST}"

if [ "$SUBMIT_FLAG" = "--submit" ]; then
  if [ "$failed" -eq 0 ]; then
    echo "Nothing to resubmit."
  else
    condor_submit "$RETRYSUB"
  fi
fi
