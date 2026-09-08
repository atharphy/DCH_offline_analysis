#!/usr/bin/env python3
"""merge_qcd_fakerate_skims.py: merges the per-input-file condor outputs from
Online_framework/DCH/qcd_fakerate_skim.py into one file per (year, sample
group), ready for BuildQCDFakeRates.C.

The skim trees are flat and unweighted-structure (no histograms to sum
specially), so plain `hadd` is the right tool -- no custom merge logic
needed. The Runs tree (MC only, one row per input file holding that file's
own genEventSumw/genEventSumw2) is also just concatenated by hadd; downstream
code sums over all rows to get the per-sample total, same convention as
central NanoAOD's own Runs tree.

Grouping:
  - Data: all PDs (EGamma/SingleElectron + SingleMuon) merge into one
    Data_<year>.root -- rows are already tagged by `cat`, and PD routing was
    only ever about which trigger each row's e/mu channel came from, not
    about keeping outputs apart.
  - MC: QCD_HT / QCD_MuEnriched / QCD_EMEnriched stay in SEPARATE merged
    files (QCD_HT_<year>.root, QCD_MuEnriched_<year>.root,
    QCD_EMEnriched_<year>.root) -- different filter efficiencies/
    normalizations, per gen_fakerate_skims.py's own docstring. Do not merge
    across these three.

Run: python3 merge_qcd_fakerate_skims.py -year 2018
     python3 merge_qcd_fakerate_skims.py -all
"""
import argparse
import glob
import os
import subprocess

SKIM_ROOT = "/eos/user/a/atahmad/DCH_offline_analysis/QCD_FakeRate_skims"
MERGED_ROOT = "/eos/user/a/atahmad/DCH_offline_analysis/QCD_FakeRate_merged"
YEARS = ["2016preVFP", "2016postVFP", "2017", "2018"]
MC_CATEGORIES = ["QCD_HT", "QCD_MuEnriched", "QCD_EMEnriched"]


def hadd(out_path, in_files, label):
    if not in_files:
        print(f"[WARNING] no input files for {label}, skipping")
        return
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if os.path.exists(out_path):
        os.remove(out_path)
    cmd = ["hadd", "-f", out_path] + in_files
    print(f">>> merging {len(in_files)} files -> {out_path}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        print(f"[ERROR] hadd failed for {label}:\n{result.stdout[-2000:]}")
    else:
        print(f"[OK] {label} -> {out_path}")


def merge_year(year):
    data_files = sorted(glob.glob(f"{SKIM_ROOT}/{year}/Data/*/*.root"))
    hadd(f"{MERGED_ROOT}/{year}/Data_{year}.root", data_files, f"Data {year}")

    for category in MC_CATEGORIES:
        mc_files = sorted(glob.glob(f"{SKIM_ROOT}/{year}/MC/{category}/*/*.root"))
        hadd(f"{MERGED_ROOT}/{year}/{category}_{year}.root", mc_files, f"{category} {year}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-year", choices=YEARS)
    p.add_argument("-all", action="store_true")
    args = p.parse_args()

    if args.all:
        for y in YEARS:
            merge_year(y)
        return

    if not args.year:
        p.error("specify -year or -all")
    merge_year(args.year)


if __name__ == "__main__":
    main()
