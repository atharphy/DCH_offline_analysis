# DCH Offline Analysis

CMS Run 2 offline analysis framework for a Double Charged Higgs decaying into multilepton final
states (electron, muon, hadronic tau) search. This repository contains the event-selection and
histogramming stage (`DCH_tauFR.C`, `DCH_tight.C`), the jet-fake-rate estimation used by the
loose-to-tight (fake-factor) method, the correction/systematic modules applied to both data and
simulation, and the plotting stage (stacked-histogram and yield-table production).

It is a curated copy of a larger working area: only the scripts, modules, and result files actually
needed to run the two main event-selection macros and the plotting scripts are included. Comments
have been stripped from all `.C`/`.h` files; this README is the documentation instead.

> Looking for the previous, shorter version of this document? See `README_v1_backup.md`.

## Table of contents

1. [Layout](#layout)
2. [What is *not* included](#what-is-not-included)
3. [Prerequisites](#prerequisites)
4. [Complete step-by-step workflow (condor)](#complete-step-by-step-workflow-condor)
5. [Running interactively (no condor)](#running-interactively-no-condor)
6. [Corrections and systematics applied](#corrections-and-systematics-applied)
7. [Adding a new systematic](#adding-a-new-systematic)
8. [Full configuration reference](#full-configuration-reference)

## Layout

```
DCH_tauFR.C, DCH_tight.C          event selection + histogram production (loose-tau-fake-rate
                                   and tight-baseline variants of the same analysis)
DCH_tauFR_SF.C                    DCH_tauFR.C variant that applies the jet->tau fake rate as a
                                   data/MC scale factor (Dependencies/fake_rates/tau_fake_rates/
                                   DY_tau_fake_rate_SF_2D.root) instead of the fake-factor reweight
MergeChunkHists.C                 merges per-chunk histogram output from the condor chunking mode
                                   (see batch/*_chunk.sub) back into one file per year
Stackhist.C                       stacked background/data comparison plots
Stackhist_multiplicity.C          same, split by lepton/jet multiplicity
YieldPlots.C                      per-category event-yield tables and bar charts
roofit_wz.C, roofit_zz.C          RooFit-based WZ/ZZ normalization scale-factor extraction

DCH_modules/                      corrections, object/event weighting, systematics plumbing
                                   used by DCH_tauFR.C and DCH_tight.C
Stack_modules/                    shared code for Stackhist.C / Stackhist_multiplicity.C
yield_modules/                    shared code for YieldPlots.C
roofit_formatting/                shared plot-style header for roofit_wz.C / roofit_zz.C

include/                          low-level helpers: branch reading (MyBranch.C), 4-vector /
                                   category utilities (Kinematics.C/.h, cat.h), MC cross
                                   sections (Xsections.C)
filemap/                          FileMap.h: per-year, per-process input file lists (data, MC, and
                                   the H++ signal samples), all read over xrootd from the FNAL LPC
                                   EOS store (root://cmseos.fnal.gov//store/user/aahmad2/run2_files/...)
json_files/                       correctionlib payloads for the official MET phi-modulation fix
roccor/                           CMS RoccoR muon momentum scale/smearing correction (code + UL
                                   calibration tables)
HTT-utilities/RecoilCorrections/  CMS-HTT MET recoil-correction and MET-systematics classes,
                                   bundled with our re-derived recoil calibration payloads
batch/                            HTCondor submission files for all scripts above, including the
                                   chunked (per-file-split) DCH_tauFR.C/DCH_tight.C submission mode
systematics_plots/                standalone plotting suite for lepton/tau ID+iso/trigger SF maps,
                                   RoccoR, tau ES, and Z pT weights -- not read by DCH_tauFR.C/
                                   DCH_tight.C, run independently to produce review plots
sf_uncertainty/                   cross-check tool (TestSFUncertainty.C) comparing the on-the-fly
                                   correctionlib SF-uncertainty recompute against precomputed skim
                                   branches, plus the correctionlib JSONs/ROOT files it reads

derivation_scripts/                scripts that produced everything under Dependencies/ (kept
                                    separate from the main pipeline: none of these are called by
                                    DCH_tauFR.C/DCH_tight.C at runtime, they only need to be rerun
                                    if a correction/fake-rate/systematic needs updating). Unlike the
                                    main pipeline, these are run with THEIR OWN subfolder as the
                                    working directory, not the repo root -- see Running below.
  etau_fake_rate/                  electron->tau fake rate (MC-based DY closure)
  tau_fr_systematics/              jet->tau fake-rate systematic-variation histograms
  recoil_corrections/              re-derivation of the HTT-utilities recoil payload
  zpt_reweighting/                 Z pT reweighting derivation
  qcd_fake_rates/                  jet->electron/muon/tau fake-rate measurement from a QCD-enriched
                                    control region (Tight-Loose matrix method, AN-19-111/X53-AN2016
                                    style). Not yet wired into DCH_tauFR.C/DCH_tight.C -- see note
                                    below. The condor skim producer itself (qcd_fakerate_skim.py,
                                    gen_fakerate_skims.py) lives in the skim-production repo, not
                                    here, and reads directly from NanoAOD rather than from anything
                                    in this repository.

Dependencies/                      result files read at runtime by DCH_tauFR.C/DCH_tight.C
  fake_rates/tau_fake_rates/       jet->tau fake rate: the final 2D pT/|eta| histogram DCH_tauFR.C
                                    reads, plus rates/ (the 6 per-process inputs
                                    tau_fr_systematics/FakeRateSources.h consolidates), plus
                                    DY_tau_fake_rate_SF_2D.root (DCH_tauFR_SF.C's SF variant)
  fake_rates/qcd_fake_rates/       QCD-region-derived jet->e/mu/tau fake-rate histograms (per year,
                                    per flavor, per Data/MC category) from derivation_scripts/
                                    qcd_fake_rates/BuildQCDFakeRates.C -- not yet read by any
                                    driver script, kept here for review/iteration
  etau_fake_rates/cut_results/     electron->tau fake rate
  systematics/results/             jet->tau fake-rate systematic variations (5 sources x 4 years
                                    + Run2 combination each)
  zpt_weights/                     Z pT reweighting histograms, per year
  wz_zz_scale_factors/             WZ and ZZ process normalization scale factors (CSV) -- also
                                    read back in by roofit_wz.C/roofit_zz.C themselves (see below)
```

## What is *not* included

- The exploratory/legacy derivation history for the jet->tau fake rate (many superseded binning
  and background-composition attempts) is not included -- only the final
  `DY_data_tau_fake_rate_2D.root` that `DCH_tauFR.C` actually reads. If that measurement needs to
  be redone from scratch, it was not a clean, single-script pipeline in the source area this was
  copied from.
- The jet->electron / jet->muon / jet->tau QCD-region fake-rate measurement
  (`derivation_scripts/qcd_fake_rates/`, `Dependencies/fake_rates/qcd_fake_rates/`) is included as
  of this update, but is still a standalone, still-being-validated effort: it is not read by
  `DCH_tauFR.C` or `DCH_tight.C` at runtime yet. The ~10k-job condor skim that feeds it
  (`qcd_fakerate_skim.py`/`gen_fakerate_skims.py`) is also not included here -- it lives in the
  skim-production repository and reads NanoAOD directly, producing the per-year/per-sample skims
  that `derivation_scripts/qcd_fake_rates/merge_qcd_fakerate_skims.py` then merges and
  `BuildQCDFakeRates.C` turns into the fake-rate histograms checked in here.
- `stacking_scripts/` (an older, pre-refactor version of the stacking code that does not use
  `Stack_modules/`) was left out in favor of the current `Stackhist.C` / `Stackhist_multiplicity.C`.
- A stray pre-edit backup snapshot of `DCH_tauFR.C`/`DCH_tight.C` (predating the chunking support
  and the precomputed-SF-branch reads) was found alongside the working area and left out as stale
  scratch, not a real dependency.

## Prerequisites

- A CMSSW release with ROOT, RooFit, and `correctionlib` available (this was built and run against
  `CMSSW_13_0_10` on lxplus, using the `correctionlib` build under
  `/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-correctionlib/...` -- see `batch/run_dch.sh` for
  the exact `ROOT_INCLUDE_PATH`/`LD_LIBRARY_PATH` additions this requires).
- The `HTT-utilities/RecoilCorrections` package checked out and built as a sibling CMSSW package,
  i.e. at `$CMSSW_BASE/src/HTT-utilities/RecoilCorrections`, built with `scram b`. The `.h`/`.cc`
  sources bundled in this repo under `HTT-utilities/` are a copy for reference and to keep the
  repository self-contained for review; they are **not** picked up automatically by `#include
  "HTT-utilities/RecoilCorrections/interface/..."` unless this repository itself is cloned into
  `$CMSSW_BASE/src/` and `HTT-utilities/RecoilCorrections` is placed (or symlinked) directly under
  `$CMSSW_BASE/src/HTT-utilities/RecoilCorrections`. The upstream package (with the full stock data
  payload) is at `https://github.com/CMS-HTT/RecoilCorrections`; this repo's copy replaces the
  `TypeI-PFMet_<year>.root` recoil calibration files with our own re-derivation (see
  `derivation_scripts/recoil_corrections/`) but keeps HTT's stock `PFMEtSys_*.root` /
  `MEtSys_2017.root` MET-systematics payloads.
- An X.509 grid proxy for reading input files over xrootd (`root://...`), referenced in the
  `.sub` files via `x509userproxy`.

Every path in every script here is relative -- nothing is hardcoded to a personal `/eos/` or `/afs/`
location, except `filemap/FileMap.h`'s `MCBASE`/`DATABASE`/`SIGNALBASE` (shared xrootd endpoints,
not a personal path, so left as-is). Two different working-directory conventions apply:

- **Main pipeline** (everything at the repo root, plus `batch/`): run with the repository root as
  the current working directory (`Dependencies/...`, `hists/...`, `json_files/...` etc. are all
  relative to there) -- exactly as they were originally (`cd Offline_framework/offline_analysis`
  before invoking `root`). All `batch/*.sub`/`run_*.sh` files already `cd` there for you.
- **`derivation_scripts/*/`**: run with that script's own subfolder as the working directory (e.g.
  `cd derivation_scripts/tau_fr_systematics` before invoking `root`) -- they reach the shared
  `Dependencies/`, `filemap/`, `include/`, `json_files/` via `../../`, and write their own
  intermediate output into a local subfolder (e.g. `results/`, `output/`) next to themselves.

## Complete step-by-step workflow (condor)

This walks through the entire pipeline exactly as it's meant to be run at scale: event selection
via condor, then plotting via condor, including the one manual step people most often trip on --
**the plotting scripts don't take an input directory as a command-line argument.** `INPUT_DIR` and
`OUTPUT_DIR` are `static std::string` constants compiled into `Stackhist.C` /
`Stackhist_multiplicity.C` / `YieldPlots.C`. If you point `DCH_tauFR.C`/`DCH_tight.C` at a new
output location (a new `outDirBase`, e.g. because you changed a correction and don't want to
overwrite the old histograms), the plotting scripts will keep reading the *old* location until you
edit their `INPUT_DIR` to match and recompile. Step 5 below is exactly that.

### Step 0 -- one-time setup

```
git clone git@github.com:atharphy/DCH_offline_analysis.git $CMSSW_BASE/src/Offline_framework/offline_analysis
cd $CMSSW_BASE/src/Offline_framework/offline_analysis
```

Place (or symlink) `HTT-utilities/RecoilCorrections` at `$CMSSW_BASE/src/HTT-utilities/RecoilCorrections`
and build it (`scram b` from `$CMSSW_BASE/src`) -- see [Prerequisites](#prerequisites). Get an X.509
grid proxy and update `x509userproxy` in every `batch/*.sub` file to point at it (and update
`request_memory`/`+JobFlavour` to taste). If this repo lives under a different `$CMSSW_BASE` or a
different username than the one it was built under, update the CMSSW path baked into every
`batch/run_*.sh` (`cd /afs/.../CMSSW_13_0_10/src`) to match.

### Step 1 -- submit event selection (`DCH_tauFR.C`)

```
cd batch
condor_submit dch_tauFR.sub
```

This queues 4 jobs (one per year: `2016preVFP`, `2016postVFP`, `2017`, `2018` -- see the `queue year
in (...)` line in `dch_tauFR.sub`), each running `root -l -b -q 'DCH_tauFR.C+("<year>",0,-1,8)'`
with 8 CPUs / 16 GB (`request_cpus`/`request_memory` in the `.sub` file; `8` is also the
`nProc` argument passed to `DCH_tauFR.C`, i.e. how many worker processes `ROOT::TProcessExecutor`
spreads the input files across -- keep these two in sync if you change one). Output lands under
`DCH_tauFR.C`'s `outDirBase` (a relative path near the top of the file, currently
`hists/run2_hists_tauFR_etau_roccor`), one subfolder per year:
`hists/run2_hists_tauFR_etau_roccor/<year>/hist_<process>.root`.

**Alternative: one job per input file.** `dch_tauFR_perfile.sub` queues one job per `(year, idx)`
pair from `filelist.txt` (240 jobs total -- regenerate this list with
`root -l -b -q 'batch/MakeFileList.C+'` if the input file counts in `filemap/FileMap.h` change),
each single-threaded (`request_cpus = 1`, `nProc=1`, `nFilesToRun=1`). This is finer-grained and
more condor-friendly for a busy pool, at the cost of many more, smaller jobs; use whichever suits
the batch system's current load. Both variants write into the exact same `outDirBase` layout.

### Step 2 -- submit event selection (`DCH_tight.C`)

```
condor_submit dch_tight.sub          # or dch_tight_perfile.sub
```

Same mechanics as step 1, but for the tight-baseline (no fake-rate) selection. Output goes to
`DCH_tight.C`'s own `outDirBase` (currently `hists/run2_hists_noFR_reccor`) -- deliberately a
*different* directory from step 1's, since the two are different selections and several plotting
scripts (`Stackhist_multiplicity.C`, `YieldPlots.C`, `roofit_wz.C`, `roofit_zz.C`) read from this
one specifically.

### Step 3 -- check the jobs

```
condor_q $USER
```

Once everything shows `0 idle, 0 running` (or is gone from the queue), check `batch/logs/*.err` for
any job that didn't produce output, and confirm the expected `.root` files exist:

```
ls hists/run2_hists_tauFR_etau_roccor/2018/
ls hists/run2_hists_noFR_reccor/2018/
```

### Step 4 -- (only if you changed `outDirBase`) update the plotting scripts' `INPUT_DIR`/`OUTPUT_DIR`

If you ran steps 1-2 straight out of the box, `INPUT_DIR` in the plotting scripts already matches
and you can skip to step 5. If you changed `outDirBase` in `DCH_tauFR.C`/`DCH_tight.C` (e.g. to
`hists/run2_hists_tauFR_myNewCorrection`) to avoid overwriting a previous production, update the
matching `INPUT_DIR` (and, if you want the new plots kept separate from old ones, `OUTPUT_DIR`) at
the top of the plotting script(s) that read it, then force ACLiC to recompile:

| Script | reads (`INPUT_DIR`) | writes (`OUTPUT_DIR`) |
|---|---|---|
| `Stackhist.C` | `DCH_tauFR.C`'s `outDirBase` | `plots/run2_plots_tauFR_etau_roccor/` |
| `Stackhist_multiplicity.C` | `DCH_tight.C`'s `outDirBase` | `multiplicity_plots/run2_plots_noFR_roccor/` |
| `YieldPlots.C` | `DCH_tauFR.C`'s `outDirBase` | `yield_plot/run2_plots_tauFR_roccor_updated_SF/` |
| `roofit_wz.C` / `roofit_zz.C` | `DCH_tight.C`'s `outDirBase` (local `inputDir` inside the function) | `normfits/` |

```
sed -i 's|hists/run2_hists_tauFR_etau_roccor/|hists/run2_hists_tauFR_myNewCorrection/|' Stackhist.C
rm -f Stackhist_C.so Stackhist_C.d Stackhist_C_ACLiC_dict_rdict.pcm   # ACLiC's stale-header cache
root -l -b -q -e 'gSystem->CompileMacro("Stackhist.C","fkO")'        # confirm it compiles before submitting
```

(ACLiC normally notices when the `.C` file itself changed and recompiles automatically on the next
`+`-suffixed run; deleting the cached `.so`/`.d`/`.pcm` first is just a guaranteed-clean way to force
it, useful if you're not sure a change was picked up.) Repeat for whichever of the four scripts in
the table above needs to point at the new location.

### Step 5 -- submit the stacked-plot jobs

```
condor_submit stackhist.sub               # reads DCH_tauFR.C's output
condor_submit stackhist_multiplicity.sub  # reads DCH_tight.C's output
```

Each queues 40 jobs from `batch/filelist_stackhist.txt` -- 8 chunks of 4 variables (there are 32
plotted variables total, see `Stack_modules/StackConfig.h`'s `allVariables`) times 5 periods
(`2016preVFP`, `2016postVFP`, `2017`, `2018`, `Run2`). Each job runs
`root -l -b -q 'Stackhist.C+("<year>",<first>,4)'` (or `Stackhist_multiplicity.C`). If you add or
remove plotted variables so the total is no longer a multiple of 4, either adjust the `4` in
`stackhist.sub`'s `arguments` line to a divisor of the new total, or regenerate
`filelist_stackhist.txt` by hand (`year,first` pairs, `first` stepping by whatever chunk size you
pick, one block of periods `2016preVFP/2016postVFP/2017/2018/Run2` per chunk-size sweep).

**Final plots land in:**
```
plots/run2_plots_tauFR_etau_roccor/<year>/              # Stackhist.C
multiplicity_plots/run2_plots_noFR_roccor/<year>/       # Stackhist_multiplicity.C
```
(`OUTPUT_DIR` + `/<year>/`, per `Stackhist.C`'s `const std::string outdir = ensureTrailingSlash(OUTPUT_DIR) + year + "/";`)
-- one `.png`/`.pdf` per (variable, region, final-state-group) combination.

### Step 6 -- submit the yield-table jobs

```
condor_submit yieldplots.sub
```

5 jobs (4 years + `Run2`), each `root -l -b -q 'YieldPlots.C+("<year>")'`. **Final plots/tables land
in** `yield_plot/run2_plots_tauFR_roccor_updated_SF/<year>/` (`YieldPlots.C`'s `OUTPUT_DIR` +
`/<year>`).

### Step 7 -- WZ/ZZ normalization fits (not batched -- quick enough to run interactively)

```
root -l -b -q 'roofit_wz.C+("Run2")'
root -l -b -q 'roofit_zz.C+("Run2")'
```

Reads `DCH_tight.C`'s output (`inputDir` inside the function, currently
`hists/run2_hists_noFR_roccor/` -- update this the same way as step 4 if `DCH_tight.C`'s
`outDirBase` changed), fits the process normalization in a control region, prints the plot to
`normfits/plots/roofit_wz_<region>_<year>.png`, and updates
`Dependencies/wz_zz_scale_factors/{wz,zz}_scale_factors.csv` in place (read-modify-write -- it reads
the existing CSV, replaces/adds this year's row, and rewrites the file, so previous years' fitted
values survive re-running for a single year).

### Recap: where everything ends up

```
hists/run2_hists_tauFR_etau_roccor/<year>/     <- DCH_tauFR.C
hists/run2_hists_noFR_reccor/<year>/           <- DCH_tight.C
plots/run2_plots_tauFR_etau_roccor/<year>/     <- Stackhist.C            (final plots)
multiplicity_plots/run2_plots_noFR_roccor/<year>/ <- Stackhist_multiplicity.C (final plots)
yield_plot/run2_plots_tauFR_roccor_updated_SF/<year>/ <- YieldPlots.C    (final plots/tables)
normfits/plots/                                <- roofit_wz.C / roofit_zz.C (final plots)
Dependencies/wz_zz_scale_factors/*.csv         <- roofit_wz.C / roofit_zz.C (updated in place)
```

All of these are relative to the repository root and get created automatically (`gSystem->mkdir(...,
kTRUE)`) the first time each script runs -- nothing needs to be created by hand.

## Running interactively (no condor)

Same scripts, same environment (source `/cvmfs/cms.cern.ch/cmsset_default.sh`, `cmsenv`, plus the
`correctionlib` exports for `DCH_tauFR.C`/`DCH_tight.C` -- see `batch/run_dch.sh` for the exact
three lines), just invoked directly instead of via `condor_submit`:

```
root -l -b -q 'DCH_tauFR.C+("2018", 0, -1, 8)'      # year, firstEntry, nEntries(-1=all), nWorkers
root -l -b -q 'DCH_tight.C+("2018", 0, -1, 8)'
root -l -b -q 'Stackhist.C+("2018", 0, -1)'          # year, firstVar, nVars(-1=all)
root -l -b -q 'Stackhist_multiplicity.C+("2018", 0, -1)'
root -l -b -q 'YieldPlots.C+("2018")'                # also accepts "Run2"
root -l -b -q 'roofit_wz.C+("Run2")'
root -l -b -q 'roofit_zz.C+("Run2")'
```

`gSystem->CompileMacro("Script.C","fkO")` from an interactive `root` session compiles without
running `main()` -- useful as a syntax-only check.

### Derivation scripts (only if a correction/fake-rate/systematic input needs regenerating)

Run from inside the script's own subfolder (see the working-directory note above):

```
cd derivation_scripts/tau_fr_systematics
root -l -b -q 'ComputeFRSystematic_DYMC.C+()'    # also QcdMc, TtMc, WJetsData, WJetsMc
cd ../etau_fake_rate
root -l -b -q 'DeriveEtauFakeRates.C+("2018")'
root -l -b -q 'RebinAndPlotFakeRates.C+()'
cd ../recoil_corrections
root -l -b -q 'DeriveRecoilCorrections.C+("2018")'
cd ../zpt_reweighting
root -l -b -q 'DeriveZPtWeights.C+("2018")'
```

Each writes into a local subfolder next to itself (`results/`, `output/`, `cut_results/`) rather
than directly into `Dependencies/`; copy the relevant output file(s) over once you're satisfied
with the result, mirroring the existing `Dependencies/` layout.

## Corrections and systematics applied

- **MET**: official phi-modulation correction (`DCH_modules/METCorrections.h`, `json_files/`),
  MET recoil correction for DY/W+jets (`DCH_modules/RecoilCorrections.h`,
  `HTT-utilities/RecoilCorrections/`), MET recoil systematics (`DCH_modules/MEtSysWrapper.h`).
- **Muons**: Rochester momentum scale/smearing correction, official method
  (`DCH_modules/RoccoRCorrections.h`, `roccor/`).
- **Z pT reweighting**: `DCH_modules/ZPtReweight.h`, weights in `Dependencies/zpt_weights/`.
- **Fake rates** (loose-to-tight / fake-factor method): jet->tau
  (`DCH_modules/TauFakeRate.h`, `Dependencies/fake_rates/`), electron->tau
  (`DCH_modules/EtauFakeRate.h`, `Dependencies/etau_fake_rates/`), each with a corresponding
  statistical/source systematic variation set applied through `DCH_modules/TauFRSystematics.h`
  and `Dependencies/systematics/results/`.
- **Truth matching / object weighting**: `DCH_modules/TruthMatching.h`, `DCH_modules/EventWeights.h`,
  `DCH_modules/ObjectAccessors.h`, `DCH_modules/PairBuilder.h` apply generator weight, pileup,
  L1-prefiring, trigger, ID/ISO scale factors already stored as branches in the input ntuples.
  `DCH_modules/FlatXsecSystematic.h` provides the flat cross-section-normalization systematic.
- **WZ/ZZ normalization**: `roofit_wz.C`/`roofit_zz.C` fit the process yield directly from data in
  a control region; the fitted scale factors are stored in `Dependencies/wz_zz_scale_factors/` (not
  currently read back into `DCH_tauFR.C`/`DCH_tight.C` automatically -- apply by hand via
  `include/Xsections.C` or a per-process weight if/when needed).

All systematic sources known to the stacking stage are enumerated in
`Stack_modules/SystematicSources.h`; the yield-table equivalent is
`yield_modules/YieldSystematics.h`.

## Adding a new systematic (e.g. a new scale-factor up/down variation)

Every event-weight systematic already in this framework (pileup, L1-prefiring, Z pT, flat
cross-section, fake-rate source/stat variations, ...) follows the same three-step pattern; a new
one -- for example, an up/down variation from an alternate SF branch not yet used -- plugs in the
same way without touching the plotting code at all:

**1. Produce the varied weight in `DCH_tauFR.C`/`DCH_tight.C` and push it onto `weightVariants`.**
`generateAndFillConfigs(...)` (in `DCH_tauFR.C`) already takes a `const vector<WeightSystematic>&
weightVariants` and, for every event, loops over it once per configuration/final state:
```cpp
for (const auto& variant : weightVariants) {
    const double wVariant = variant.weight * wFR * objectSF;
    fillConfiguration(cfg, wVariant, ..., variant.suffix);
}
```
(`WeightSystematic` is just `{std::string suffix; double weight;}`, defined in
`DCH_modules/SystematicPlan.h`.) This loop is generic -- it doesn't know or care what's in
`weightVariants`, so a new source only needs to get added to that vector before this point in the
event loop. The existing sources show the pattern (`DCH_modules/EventWeights.h`):
```cpp
inline void appendPileupWeightVariants(double evtWeight, std::vector<WeightSystematic>& weightVariants) {
    if (weightPUtruejson <= 0.0) return;
    weightVariants.push_back({"_puUp", evtWeight / weightPUtruejson * weightPUtruejson_up});
    weightVariants.push_back({"_puDown", evtWeight / weightPUtruejson * weightPUtruejson_down});
}
```
i.e. divide the nominal weight by the nominal factor and multiply back in the varied one, so you
don't have to recompute the full weight from scratch. Write an analogous `appendMyNewSFVariants(...)`
in `DCH_modules/EventWeights.h` (if it's a per-event or per-object SF read the same way the nominal
SF is, e.g. via a new accessor in `DCH_modules/ObjectAccessors.h` alongside `idSFByIndex`/
`isoSFByIndex` if the alternate value lives in a branch like `IDSF_1_up`/`IDSF_1_down`), then call it
next to the existing `appendPileupWeightVariants(evtWeight, weightVariants);` /
`appendL1PrefiringWeightVariants(evtWeight, weightVariants);` calls in `DCH_tauFR.C` (and `DCH_tight.C`
if the baseline selection should carry it too). Only MC needs this -- data isn't reweighted, so
there's nothing to branch on for `isData`.

**2. Register the new suffix pair so the plotting stage picks it up automatically.** Add one line to
`kAllKnownSystSources` in `Stack_modules/SystematicSources.h` (and to `kAllSystSources` in the same
file if it should be included in the drawn total-uncertainty band, not just tracked -- the two lists
exist because a few fake-rate sub-sources are computed but not currently drawn), and the matching
list in `yield_modules/YieldSystematics.h`:
```cpp
{"myNewSF", "_myNewSFUp", "_myNewSFDown"},
```
`Stackhist.C`/`Stackhist_multiplicity.C`/`YieldPlots.C` iterate these lists generically to build the
envelope and the systematic breakdown table -- no other code change is needed for the new source to
show up in the output.

**3. Recompile and rerun.** `DCH_tauFR.C`/`DCH_tight.C` first (produces the new `_myNewSFUp`/
`_myNewSFDown`-suffixed histograms alongside the existing ones, per year), then `Stackhist.C`/
`Stackhist_multiplicity.C`/`YieldPlots.C` (folds the new source into the band/table). ACLiC only
recompiles a `.C` file whose source actually changed, so this is safe to do incrementally --
`_C.so`/`.d`/`_ACLiC_dict_rdict.pcm` next to any file you touch are ACLiC's cache; delete them if a
rebuild seems stale (e.g. after editing a `.h` a `.C` includes, which ACLiC doesn't always notice).

A fake-rate-specific systematic (a new alternate fake-rate *source*, e.g. a different background
composition or binning) follows a slightly different path: it's a `FRSystSource` (reader + suffix +
flavor flag), appended to `frSourceVariants` where `DCH_tauFR.C` builds the fake-rate readers, mirroring
how `frDyMc`/`frWJets`/`frWJetsMc`/`frQcdMc`/`frTtMc` are already wired via
`DCH_modules/TauFRSystematics.h` and `Dependencies/systematics/results/<source>/`. Step 2 (registering
in both `SystematicSources.h` files) is identical either way.

## Full configuration reference

Every user-facing toggle, threshold, list, and lookup table in the framework, by file. All of these
are plain constants near the top of their file -- change the value, delete that file's ACLiC cache
(`*_C.so`/`.d`/`_ACLiC_dict_rdict.pcm`) if you're not relying on an automatic recompile, and rerun.

### `DCH_modules/CommonConfig.h` -- top-level correction toggles (read by both `DCH_tauFR.C` and `DCH_tight.C`)

| Name | Default | Effect |
|---|---|---|
| `APPLY_OFFICIAL_MET_CORRECTION` | `true` | Apply the official MET phi-modulation correction (`METCorrections.h`) |
| `APPLY_ZPT_REWEIGHTING` | `true` | Apply Z pT reweighting to DY MC (`ZPtReweight.h`) |
| `APPLY_DY_RECOIL_CORRECTION` | `true` | Apply MET recoil correction to DY MC |
| `APPLY_WJ_RECOIL_CORRECTION` | `false` | Apply MET recoil correction to W+jets MC |
| `USE_CUSTOM_RECOIL_CORRECTIONS` | `true` | Use the re-derived per-year `TypeI-PFMet_<year>.root` payloads instead of HTT's stock Run2016BtoH-combined ones (`RecoilCorrections.h`) |
| `NlepMax` | `4` | Max leptons per event the framework's fixed-size arrays (`h_pt[NlepMax]`, etc.) support -- raising this needs matching branches (`pt_5`, ...) in the input ntuple, not just this constant |
| `finalStates` | 45 category strings | Every lepton-flavor final state the framework recognizes (2- through 4-lepton) |

### `DCH_modules/RoccoRCorrections.h`

| Name | Default | Effect |
|---|---|---|
| `APPLY_ROCCOR_DATA` | `true` | Apply Rochester correction to data muons |
| `APPLY_ROCCOR_MC` | `true` | Apply Rochester correction (spread/scale) to MC muons |

### `DCH_modules/FlatXsecSystematic.h`

`kFlatXsecUncertainties`: a `(filename-prefix, fractional uncertainty)` list -- the flat
cross-section-normalization systematic (`_xsecUp`/`_xsecDown`) applied per process, matched by
input filename prefix (`baseName.rfind(prefix, 0) == 0`). Add a new process by appending
`{"MyNewSample", 0.05}` (5% uncertainty); any file whose name doesn't match any entry gets `0.0`
(no `_xsec` variant produced for it).

### `DCH_modules/HistUtils.h` -- raw per-lepton/per-event histogram binning

`S_mZ1`, `S_mZ2`, `S_mH1`, `S_mH2`, `S_zPt`, `S_met`, `S_metphi`, `S_LT`, and the per-lepton
`S_pt[NlepMax]`, `S_eta[NlepMax]`, `S_phi[NlepMax]`, `S_d0[NlepMax]`, `S_dZ[NlepMax]`,
`S_iso[NlepMax]` -- each a `{name, title, bins, xmin, xmax}` `HistSpec`. This is the binning
`DCH_tauFR.C`/`DCH_tight.C` actually fill at production time (deliberately fine -- e.g. 1000 bins,
0-1000 GeV for `mZ1`); `Stack_modules/StackConfig.h`'s `binning` map (below) rebins these down for
plotting. Change a range here only if events are landing outside it (they'd silently fall in the
under/overflow bin and not show up).

### `DCH_modules/TauFRSystematics.h`

- `kTauFRSystematicsBase` (`"Dependencies/systematics/results"`): where the 5 fake-rate systematic
  sources are read from.
- `kTauFRSystematicSources`: the 5 `(directory-name, suffix)` pairs -- `dy_mc`/`DyMc`,
  `wjets_data`/`WJets`, `wjets_mc`/`WJetsMc`, `qcd_mc`/`QcdMc`, `tt_mc`/`TtMc`. Adding a 6th source
  needs a matching subfolder under `Dependencies/systematics/results/`, a new entry here, a new
  `FRSystSource` pushed onto `frSourceVariants` in `DCH_tauFR.C` (see "Adding a new systematic"
  above), and a new line in both `SystematicSources.h` files.

### `include/Xsections.C`

`XSec(filename)`: a long `if/else if` chain matching MC filename substrings to their cross section
in pb (data files return `1`; a file matching nothing falls through -- check the bottom of the
function for the default). To add a new MC sample, add a new `else if (fname.find("MySample") <
fname.length()) return <xsec_in_pb>;` line, and add the file to the relevant process group in
`filemap/FileMap.h`'s `getFileMap(year)`.

### `Stack_modules/StackConfig.h` -- plotting-stage configuration (`Stackhist.C`/`Stackhist_multiplicity.C`)

| Name | Purpose |
|---|---|
| `year` | Set at runtime from the script's `inYear` argument; used by `Labels.h`'s luminosity label |
| `SIGNAL_MASS_FILTER` | If non-empty, only this H++ mass point (e.g. `"HppM500"`) is drawn as the signal overlay |
| `binning` | `variable -> {nbins, xmin, xmax}` (or a `variable:region -> ...` override, e.g. `"mZ1:SR"`) rebin map applied on top of the fine production binning from `HistUtils.h` |
| `allVariables` | The 32 variables plotted -- this list's length determines how `batch/filelist_stackhist.txt` should be chunked (currently groups of 4) |
| `finalStates` | Lepton-flavor final states considered for grouping/labeling |
| `regions` | All 23 analysis regions (`DYCR_*`, `CR_*`, `VR_*`, `SR_*`) the stacking code recognizes by histogram-name suffix |

Also in `Stackhist.C`/`Stackhist_multiplicity.C` themselves (top of file, before the `#include
"Stack_modules/..."` lines -- these have to be defined before the headers that use them):

| Name | Default (`Stackhist.C`) | Default (`Stackhist_multiplicity.C`) | Effect |
|---|---|---|---|
| `DRAW_BANDS` | `true` | `true` | Draw the total-systematic uncertainty band on the stack |
| `USE_LOG_Y` | `false` | `false` | Log-scale y-axis |
| `EVENTS_PER_BIN_WIDTH` | `false` | `false` | Divide bin content by bin width (for variable-width binning) |
| `blind_SR` | `true` | *(not defined -- always unblinded)* | Blind the data points in `SR_*`-suffixed regions |
| `INPUT_DIR` | `hists/run2_hists_tauFR_etau_roccor/` | `hists/run2_hists_noFR_roccor/` | Where to read `DCH_tauFR.C`/`DCH_tight.C` output from |
| `OUTPUT_DIR` | `plots/run2_plots_tauFR_etau_roccor/` | `multiplicity_plots/run2_plots_noFR_roccor/` | Where plots are written |
| `fill_colors` (inside the `Stackhist(...)`/`Stackhist_multiplicity(...)` function body) | 14-process color map | same | ROOT color index per background process, `data`, and `signal` |

The background stacking order itself (`bkgOrder`, bottom-to-top) is in `Stack_modules/StackDraw.h`:
`DY, DY10_50, WZ, WW, ZZ, TTbar, ttV, WJ, ST, VVV, QCD, other, signal`.

### `Stack_modules/SystematicSources.h` / `yield_modules/YieldSystematics.h`

See [Adding a new systematic](#adding-a-new-systematic) above -- `kAllKnownSystSources` /
`kAllSystSources` are the full registry and the drawn-envelope subset, respectively.

### `yield_modules/YieldConfig.h` -- yield-table configuration (`YieldPlots.C`)

| Name | Purpose |
|---|---|
| `kFinalStates` | Same 45 final states as `CommonConfig.h`/`StackConfig.h` (kept as an independent copy here) |
| `kProcesses` | The 13 process groups summed into the yield table, in order |
| `kColors` | ROOT color index per process (mirrors `Stackhist.C`'s `fill_colors`) |
| `kLabels` | Display label per process (currently identity except `data` -> `"Data"`) |
| `kBinSuffixes` / `kBinLabels` | The 7 category bins (`0tau`...`3lep2tau`) and their axis labels for the yield bar chart |
| `yearsFor(year)` | Expands `"Run2"` -> all 4 years, `"2016"` -> `preVFP`+`postVFP`, else passes through |
| `periodLabel(year)` | Integrated-luminosity label per period (keep in sync with `Stack_modules/Labels.h`'s `getLumiLabel` if luminosities are revised) |

Also in `YieldPlots.C` itself (top of file): `DRAW_BANDS` (`false`), `USE_LOG_Y` (`true`),
`blind_SR` (`true`), `INPUT_DIR` (`hists/run2_hists_tauFR_roccor/`), `OUTPUT_DIR`
(`yield_plot/run2_plots_tauFR_roccor_updated_SF/`) -- same meaning as the `Stackhist.C` table above.

### `roofit_wz.C` / `roofit_zz.C`

Top of file: `SIGNAL_MASS_FILTER` (empty by default -- no signal overlay in the fit plot),
`ZZ_SF_CSV`/`WZ_SF_CSV` (`Dependencies/wz_zz_scale_factors/{zz,wz}_scale_factors.csv` -- read *and*
rewritten by the script). Inside the `roofit_wz(year)`/`roofit_zz(year)` function body: `region`
(the control region fit against -- `"CR_3lep0tau"` for WZ, `"CR_0tau"` for ZZ) and `inputDir`
(`hists/run2_hists_noFR_roccor/`, i.e. `DCH_tight.C`'s output).

### `batch/*.sub` / `batch/run_*.sh`

Per-job resources (`request_cpus`, `request_memory`, `request_disk`, `+JobFlavour`) and the
environment setup (CMSSW path, `x509userproxy`, `correctionlib` CVMFS path) -- see
[Step 0](#step-0----one-time-setup) above for what needs updating when running under a different
account/CMSSW area. `batch/MakeFileList.C` regenerates `filelist.txt` from `filemap/FileMap.h`'s
actual per-year file counts if the input file lists change; `filelist_stackhist.txt` is a plain
`year,first` list and needs updating by hand if the chunk size or variable count changes (see
[Step 5](#step-5----submit-the-stacked-plot-jobs)).
