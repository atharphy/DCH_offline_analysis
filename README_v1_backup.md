# DCH Offline Analysis

CMS Run 2 offline analysis framework for a di-Higgs / multi-lepton (electron, muon, hadronic tau)
search. This repository contains the event-selection and histogramming stage (`DCH_tauFR.C`,
`DCH_tight.C`), the jet-fake-rate estimation used by the loose-to-tight (fake-factor) method, the
correction/systematic modules applied to both data and simulation, and the plotting stage
(stacked-histogram and yield-table production).

It is a curated copy of a larger working area: only the scripts, modules, and result files actually
needed to run the two main event-selection macros and the plotting scripts are included. Comments
have been stripped from all `.C`/`.h` files; this README is the documentation instead.

## Layout

```
DCH_tauFR.C, DCH_tight.C          event selection + histogram production (loose-tau-fake-rate
                                   and tight-baseline variants of the same analysis)
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
batch/                            HTCondor submission files for all five scripts above

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

Dependencies/                      result files read at runtime by DCH_tauFR.C/DCH_tight.C
  fake_rates/tau_fake_rates/       jet->tau fake rate: the final 2D pT/|eta| histogram DCH_tauFR.C
                                    reads, plus rates/ (the 6 per-process inputs
                                    tau_fr_systematics/FakeRateSources.h consolidates)
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
- The newer jet->electron / jet->muon / jet->tau fake-rate measurement (a separate, still-unwired,
  still-being-validated effort) is intentionally excluded. It is not read by `DCH_tauFR.C` or
  `DCH_tight.C` and is not part of this analysis yet.
- `stacking_scripts/` (an older, pre-refactor version of the stacking code that does not use
  `Stack_modules/`) was left out in favor of the current `Stackhist.C` / `Stackhist_multiplicity.C`.

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
  before invoking `root`).
- **`derivation_scripts/*/`**: run with that script's own subfolder as the working directory (e.g.
  `cd derivation_scripts/tau_fr_systematics` before invoking `root`) -- they reach the shared
  `Dependencies/`, `filemap/`, `include/`, `json_files/` via `../../`, and write their own
  intermediate output into a local subfolder (e.g. `results/`, `output/`) next to themselves.

## Compiling and running the full pipeline

Every script here is a ROOT macro compiled with ACLiC -- the trailing `+` on `root -l -b -q
'Script.C+(...)'` is what triggers compilation (a `Script_C.so`/`.d`/`_ACLiC_dict_rdict.pcm` get
built next to the source the first time; ACLiC only recompiles when the source changes, so
subsequent runs of the same script can drop the `+` if you prefer, though leaving it on is harmless
and safest). `gSystem->CompileMacro("Script.C","fkO")` from an interactive `root` session is
equivalent and useful for a syntax-only check without running `main()`.

### 0. Environment (once per shell)

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_BASE/src        # the CMSSW area this repo is checked out under
eval `scramv1 runtime -sh`
cd Offline_framework/offline_analysis   # repository root -- see the working-directory note above
```

`DCH_tauFR.C` and `DCH_tight.C` additionally need `correctionlib` on the include/library path
(everything else in this repo does not):

```
CORR_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-correctionlib/2.1.0-6dc02863165bc2126b8299c5b63785af/lib/python3.9/site-packages/correctionlib
export ROOT_INCLUDE_PATH="${CORR_BASE}/include:${ROOT_INCLUDE_PATH:-}"
export LD_LIBRARY_PATH="${CORR_BASE}/lib:${LD_LIBRARY_PATH:-}"
```

(`batch/run_dch.sh` sets these same two exports for the batch jobs -- copy the version number from
there if the `correctionlib` CVMFS build has moved on since this was written.)

### 1. Event selection -> histograms

```
root -l -b -q 'DCH_tauFR.C+("2018", 0, -1, 8)'      # year, firstEntry, nEntries(-1=all), nWorkers
root -l -b -q 'DCH_tight.C+("2018", 0, -1, 8)'
```

Repeat per year (`2016preVFP`, `2016postVFP`, `2017`, `2018`). `DCH_tauFR.C` produces the
loose-to-tight fake-rate-weighted histograms (jet->tau, jet->electron, electron->tau fake
contributions applied as event weights, plus the full systematic-variant set -- see "Adding a new
systematic" below); `DCH_tight.C` produces the tight-selection-only baseline without the fake-rate
machinery. Both read `DCH_modules/CommonConfig.h` for the top-level correction toggles
(`APPLY_OFFICIAL_MET_CORRECTION`, `APPLY_ZPT_REWEIGHTING`, `APPLY_DY_RECOIL_CORRECTION`,
`APPLY_WJ_RECOIL_CORRECTION`, `USE_CUSTOM_RECOIL_CORRECTIONS`) and write output under `outDirBase`
(a relative path, `hists/...`, near the top of each file -- created under the repository root at
runtime; point it elsewhere, e.g. your own EOS workspace, if you'd rather not keep multi-GB
histogram output inside the repo checkout).

### 2. Plotting (reads step 1's output)

```
root -l -b -q 'Stackhist.C+("2018", 0, -1)'          # year, firstVar, nVars(-1=all)
root -l -b -q 'Stackhist_multiplicity.C+("2018", 0, -1)'
root -l -b -q 'YieldPlots.C+("2018")'                # also accepts "Run2"
```

These read `DCH_tauFR.C`/`DCH_tight.C`'s `hists/...` output back in via their own `INPUT_DIR` and
write plots to `OUTPUT_DIR` (also relative, also overridable) -- `INPUT_DIR` must point at whichever
of the two step-1 outputs the plot is meant to compare against data, so keep it in sync with
whichever `outDirBase` you actually ran.

### 3. Normalization fits (independent of steps 1-2, reads `DCH_tight.C`'s output directly)

```
root -l -b -q 'roofit_wz.C+("Run2")'
root -l -b -q 'roofit_zz.C+("Run2")'
```

### Batch (HTCondor)

```
cd batch
condor_submit dch_tauFR.sub          # 4 jobs, one per year, 8 cores/16GB each
condor_submit dch_tight.sub
condor_submit stackhist.sub          # chunked over (year, variable-range) pairs, see filelist_stackhist.txt
condor_submit stackhist_multiplicity.sub
condor_submit yieldplots.sub         # 5 jobs: 4 years + Run2
```

`batch/MakeFileList.C` regenerates `filelist.txt`/`filelist_stackhist.txt` from `filemap/FileMap.h`'s
actual per-year file counts if the input file lists change.

Every `.sub`/`run_*.sh` file has the CMSSW area path, grid-proxy path, and `correctionlib` CVMFS
path hardcoded to the environment this was built in -- update those three before submitting from a
different account or CMSSW area.

### Derivation scripts (only if a correction/fake-rate/systematic input needs regenerating)

Run from inside the script's own subfolder (see the working-directory note above), same
environment as step 0:

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
