# DCH Offline Analysis

Offline framework for the doubly-charged Higgs (H±±, "DCH") multi-lepton
search: tight-selection histogram/tree production from NanoAOD-derived
skims, condor batch processing, and downstream stacking.

This snapshot covers the **tight-selection histogram pipeline**
(`DCH_tight*.C` and everything it depends on). The tau-fake-rate
**derivation** pipeline (`DCH_tauFR.C`, `DCH_tauFR_SF.C`) is intentionally
not included yet — it will be added separately. Modules that the tight
pipeline uses to *apply* an already-derived fake rate (`TauFakeRate.h`,
`TauFakeSFReader.h`, `TauFRSystematics.h`) are included, since they're a
runtime dependency of `DCH_tight.C` itself. `EtauFakeRate.h` and
`etau_fake_rates/` were left out for the same reason — nothing in this
snapshot calls them.

## Requirements

- A CMSSW area with ROOT 6.26+ (tested under `CMSSW_13_0_10`, `el9_amd64_gcc11`,
  on lxplus). If you don't have one yet:
  ```bash
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  cmsrel CMSSW_13_0_10          # or: scram project CMSSW CMSSW_13_0_10
  cd CMSSW_13_0_10/src
  cmsenv                        # equivalent to `eval `scramv1 runtime -sh``
  git clone <this-repo-url> DCH_offline_analysis
  cd DCH_offline_analysis
  ```
- **Build the vendored `HTT-utilities` package once** (needed before running
  anything that uses recoil corrections, i.e. essentially every driver):
  ```bash
  cd $CMSSW_BASE/src
  scram b -j4
  ```
  This isn't optional bookkeeping: `DCH_tight*.C` load recoil corrections via
  `R__LOAD_LIBRARY(libHTT-utilitiesRecoilCorrections.so)`, a *compiled*
  library that CMSSW's build system discovers from `HTT-utilities/`'s
  `BuildFile.xml` — ACLiC alone (`.C+`) never builds it. `scram b` finds and
  builds `BuildFile.xml`-bearing directories anywhere under `$CMSSW_BASE/src`
  recursively, including nested inside this repo, so no special path setup is
  needed beyond running it once from `$CMSSW_BASE/src` after cloning. Rebuild
  (`scram b -j4`) any time you pull changes to `HTT-utilities/`.
- [`correctionlib`](https://github.com/cms-nanoAOD/correctionlib), via CVMFS:
  ```bash
  CORR_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-correctionlib/2.1.0-6dc02863165bc2126b8299c5b63785af/lib/python3.9/site-packages/correctionlib
  export ROOT_INCLUDE_PATH="${CORR_BASE}/include:${ROOT_INCLUDE_PATH:-}"
  export LD_LIBRARY_PATH="${CORR_BASE}/lib:${LD_LIBRARY_PATH:-}"
  ```
  Every `root -l -b -q` invocation below assumes `cmsenv` has been sourced,
  `correctionlib` is exported as above, `HTT-utilities` has been built, and
  the shell's current directory is the repository root.
- A valid grid proxy, for reading input skims over `xrootd` and for condor
  job submission:
  ```bash
  voms-proxy-init --voms cms --valid 168:00
  ```
- HTCondor (lxplus) for batch production.

Everything else the pipeline needs at runtime — recoil-correction payloads,
Rochester muon-correction tables, Z-pT/tau-FR/QCD-FR weight files, golden
JSONs — is vendored in this repository and resolved with plain paths
relative to the repository root (never `$CMSSW_BASE`-relative, never a
hardcoded personal `/eos` or `/afs` path), so no separate checkout of
external packages is required beyond `correctionlib` itself, and cloning
this repo anywhere under `$CMSSW_BASE/src` is sufficient — no sibling
packages need to be checked out alongside it. The only paths that are
*meant* to stay external are input skim locations (`filemap/FileMap.h`) and
each driver's own output directory (`DCH_tight*.C`, `Stackhist*.C` — see
"Where output directories are set" below).

## Directory map

```
DCH_tight.C                 driver: run over full files (or a subrange), 1..N processes
DCH_tight_Chunk.C           driver: process one 500k-event chunk of one file (unbatched)
DCH_tight_ChunkBatch.C      driver: process a batch of consecutive chunks, merged on local scratch (production path)
MergeChunkHists.C           merges per-chunk hist_/masstree_ files back into one file per sample
Stackhist.C                 downstream: data/MC stack plots from the merged per-sample histograms
Stackhist_multiplicity.C    downstream: same, split/binned by lepton multiplicity
roofit_wz.C, roofit_zz.C    downstream: RooFit-based WZ/ZZ background fits
parallel_merge.sh           shell helper for running several MergeChunkHists.C jobs in parallel

DCH_modules/                 the tight-selection pipeline itself (see table below)
Mass_reco_modules/           DCH signal-shape mass reconstruction (mDCH1/mDCH2), 3-/4-lepton events
include/                     file map input, gen-level kinematics helpers, branch loader, cross sections
filemap/                     FileMap.h: year -> {process -> [skim files]}
batch/                       condor submission files + helper scripts (see "Running the pipeline")
systematics_plots/           standalone ROOT macros that plot each systematic's SF curve
Stack_modules/                helpers used by Stackhist.C / Stackhist_multiplicity.C
sf_uncertainty/               standalone SF-uncertainty testing code + correctionlib JSONs (not a runtime dependency)
recoil_studies/                standalone macros that derived the recoil-correction parameters
systematics/results/           cached intermediate results from various systematic derivations

roccor/                      Rochester muon-momentum correction library (vendored, third-party)
HTT-utilities/                 recoil-correction library (vendored, third-party -- see note below)
fake_rates/                    tau/QCD fake-rate ROOT files + the scripts that built them
json_files/                     golden JSON (lumi mask) files, gzipped, per year
zpt_weights/                    Z-pT reweighting histograms, per year (vendored)
normfits/                       WZ/ZZ normalization scale factors used by roofit_wz.C/roofit_zz.C (vendored)
```

**Note on `HTT-utilities/`**: the upstream `RecoilCorrector`/`MEtSys`
constructors (`src/RecoilCorrector.cc`, `src/MEtSys.cc`) hardcode their
payload path as `$CMSSW_BASE/src/<fileName>`, which assumes `HTT-utilities`
is checked out as its own top-level sibling package directly under
`$CMSSW_BASE/src/` -- not nested inside another repo. Since it's vendored
*inside* this repo instead, both constructors were patched here to just use
the path they're given as-is (repo-relative, matching every other module).
This is the one piece of vendored third-party code in this repo that isn't
byte-for-byte upstream; everything else is copied verbatim.

## DCH_modules — what each file does and how it links to DCH_tight.C

`DCH_tight.C` / `DCH_tight_Chunk.C` / `DCH_tight_ChunkBatch.C` all funnel
into a single entry point, **`ProcessTightFile()`** in
`DCH_modules/TightFileProcessor.h` — that's the one function that ties
every module below together into the actual per-event loop. Everything
else in `DCH_modules/` is either a helper `TightFileProcessor.h` calls
directly, or a helper one of those helpers calls. The table below is
ordered the same way `TightFileProcessor.h` actually uses them, not
alphabetically, so it doubles as a walk-through of the event loop.

| File | Role |
|---|---|
| `FileJob.h` | Plain struct describing one processing job: year, sample name, file name, output dir(s), luminosity, and (for chunked jobs) the entry range + chunk index. This is what gets passed into `ProcessTightFile`. |
| `FileSetup.h` | `setupFileContext()`: opens the input file/tree, figures out the sample's cross-section weight and whether it's Data/DY/WJets/TTbar, loads the Z-pT weight histogram, recoil corrector, MET-systematics reader and official-MET-correction payload for the right year, and bundles it all into a `FileContext`. |
| `BranchSetup.h` | `enableBranches()`: turns on only the NanoAOD branches actually read (`SetBranchStatus`) for read-speed; called from `FileSetup.h`. |
| `CommonConfig.h` | Shared constants: max lepton multiplicity (`NlepMax`), the list of valid final-state category strings (`finalStates`), and **boolean toggles for whole correction sources** — see "Turning correction sources on/off" below. |
| `EventSelection.h` | `passTightEvent()`, `passChargeTopology()`, `hasDuplicateTightObjects()`, `getLT()` — the tight preselection cuts applied per event, before any weight/kinematics work happens. |
| `include/Kinematics.C`/`.h`, `ObjectAccessors.h` | Low-level accessors: `numberToCat()` (int category code -> string like `"eemm"`), index-based getters (`ptByIndex`, `etaByIndex`, ...) into the flat `pt_1..pt_4`-style branches, OS/SS pairing helpers (`processPairs`, `removeOverlap`), and the Zwin/Zveto classifier. |
| `PairBuilder.h` | `buildPairs()`: builds the merged OS/SS lepton-pair lists for a category string, on top of `processPairs`/`removeOverlap`. |
| `RoccoRCorrections.h` | Loads the Rochester correction tables (`roccor/`), applies the nominal muon-momentum correction/smearing, and provides `shiftMuonPts()`/`restorePts()` for the `_roccorUp`/`_roccorDown` shape systematic. |
| `METCorrections.h` | `loadOfficialMETCorrections()` / `applyOfficialMETCorrection()`: correctionlib-based Type-1 MET correction, gated by `APPLY_OFFICIAL_MET_CORRECTION`. |
| `GenBosonMomentum.h` | `getGenBosonMomentum()`: gen-level vector-boson (full + visible) momentum — the input the recoil correction needs, for DY/WJets MC. |
| `RecoilCorrections.h` | `loadRecoilCorrector()` / `applyRecoilCorrection()`: wraps the vendored `HTT-utilities` recoil correction library. |
| `MEtSysWrapper.h` | `loadMEtSys()` / `computeRecoilSystShift()`: the recoil-response/-resolution systematic shifts layered on top of the nominal recoil correction. |
| `EventWeights.h` | `computeNominalWeight()` (gen weight × object SFs × trigger SF × top-pT weight) and `appendStandardWeightVariants()` (builds every *weight-only* systematic variant — see the systematics table below). This is where most systematics are wired in. |
| `TopPtReweighting.h` | The CMS Top PAG NNLO-NLO top-pT reweighting SF and its Up/Down bracket, applied only to simulated ttbar. See the header itself for the formula and its source. |
| `TauFakeRate.h`, `TauFakeSFReader.h`, `TauFRSystematics.h` | Read already-derived tau fake-rate files (`fake_rates/tau_fake_rates/`) and evaluate the nominal fake-rate weight and its systematic-variation readers. (Not called by the nominal tight-selection path directly today — kept here as the application-side counterpart of the not-yet-migrated `DCH_tauFR.C`.) |
| `FlatXsecSystematic.h` | `getFlatXsecUncertainty()`: a hardcoded per-process flat cross-section uncertainty table, used for `_xsecUp`/`_xsecDown`. |
| `GenZMatching.h` | `matchGenZ()`: gen-level Z matching used only to pick the right generated Z pT for the Z-pT reweighting. |
| `ZPtReweight.h` | `loadZPtWeights()` / `getZPtWeight()`: reads and applies the Z-pT reweighting histogram + its uncertainty, gated by `APPLY_ZPT_REWEIGHTING` and DY-only. |
| `PairBuilder.h` (again) + `EventKinematics.h` | `buildEventKinematics()`: computes the per-event kinematic bundle (mZ1/mZ2/mH1/mH2, MET, LT, Z pT, ...); `collectHistKeys()` works out which region/category histogram keys this event falls into. |
| `SystematicPlan.h` | The `WeightSystematic` / `MetSystematic` structs that carry a systematic's name suffix + shifted weight (or shifted MET px/py) around, plus `cartesianToPolar()`. |
| `TauESSystematic.h` | `computeTauESRelShift()` / `shiftTauESPts()`: tau energy-scale relative pT shifts, for the `_tauESUp`/`_tauESDown` shape systematic. |
| `HistUtils.h` | Generic per-key `TH1D` get-or-create/fill/write helpers (`getHist`, `fillOne`, `writeHists`) shared by every histogram category — this is the low-level primitive `HistFilling.h` is built on. |
| `HistFilling.h` | The `HistBundle` struct (one `unordered_map<string,TH1D*>` per variable, keyed by region+systematic suffix) and `fillEventHistograms()`, which fills the nominal + every weight-variant + every MET-variant histogram for a given event in one call. |

## Mass_reco_modules — signal-shape mass reconstruction

Filled only for 3- and 4-lepton events (`catstr.size() == 3 || 4`), gated
inside `TightFileProcessor.h`, and only when the driver set
`job.treeOutDir` (currently: `DCH_tight_Chunk.C` and, indirectly through
its own per-batch local jobs, `DCH_tight_ChunkBatch.C`).

| File | Role |
|---|---|
| `SignalShapeReco.h` | `computeSignalShapeMasses()`: runs both neutrino-pz solvers below and keeps whichever same-sign-pair combination gives the smaller `\|mDCH1 - mDCH2\|`, returning the two DCH candidate masses. |
| `SignalShapeSelection.h` | `catToNumber()` (reverse of `numberToCat`), `buildSignalShapePairs()` (OSSF/SS pairs for this selection, reusing `processPairs`/`removeOverlap`), `computeSignalShapeVariables()` (the ~30-variable bundle: mll1/mll2/mDCH1/mDCH2/ST/mZ1-4/mT.../dR.../cutflow), and the tau-multiplicity region classifier (0tau/1tau/2tau/3tau, 3lep0tau/3lep1tau/3lep2tau) with each region's n-1 cutflow thresholds. |
| `SignalShapeHistograms.h` | `SignalShapeHistBundle`: a `TTree` (`Events_<region>`, lazily created per region on first fill — see "Why a Data chunk might have no 0tau/1tau/... tree" below) plus its full histogram set, per region and per systematic suffix. `fillSignalShapeVariant()` fills one variant into the right region's tree+histograms. |
| `massreco/rjm_test.C` (`get_legs_comb`), `massreco/rjm_guru.C` (`get_legs_tau`) | The two independent neutrino-pz solvers the reconstruction is built on. Self-contained — only depend on the `3nu_*.h`/`fir_alphas_*.h` headers alongside them, not on anything else in this repo. |

### Why a Data chunk might have no `Events_0tau`/`1tau`/`2tau`/`3tau` tree

`SignalShapeHistBundle::getTree()` only inserts a region into its internal
map the first time an event actually lands in that region, and `writeAll()`
only writes a tree if it has `GetEntries() > 0`. So if zero events in a
given chunk satisfy "3- or 4-lepton, ≥1 OSSF pair, passes that region's n-1
cutflow", that region's tree simply won't exist in that chunk's output file
at all — not a bug, just no events. This is common for Data (background-
dominated) chunks and essentially never happens for signal MC, which is
enriched in exactly this topology.

## Turning correction sources on/off

`DCH_modules/CommonConfig.h` has one boolean per whole correction source:

```cpp
const bool APPLY_OFFICIAL_MET_CORRECTION = true;
const bool APPLY_ZPT_REWEIGHTING         = true;
const bool APPLY_DY_RECOIL_CORRECTION    = true;
const bool APPLY_WJ_RECOIL_CORRECTION    = true;
const bool USE_CUSTOM_RECOIL_CORRECTIONS = true;
const bool APPLY_ROCCOR_CORRECTION       = true;
```

Flip one to `false` to disable that whole correction (and, for the ones
with an associated systematic, its Up/Down variant disappears too, since
the variant-building code is gated on the same flag or on the correction
having actually run). There's no flag for top-pT reweighting or trigger/
object SFs — those are unconditional for MC (top-pT is additionally gated
on `isTTbar` inside `TopPtReweighting.h`/`EventWeights.h`).

## Systematics reference: what exists, and how to add or remove one

Every systematic ends up as either (a) an entry in the `weightVariants`
vector, (b) an entry in the `metVariants` vector, or (c) a separate
re-histogramming block in `TightFileProcessor.h`'s event loop, depending on
whether it only rescales the event weight, only shifts MET, or actually
moves lepton momenta (and therefore every derived kinematic quantity).

**Category (a) — weight-only systematics** (built in
`EventWeights.h::appendStandardWeightVariants`, plus one exception built
directly in `TightFileProcessor.h`):

| Suffix pair | Source |
|---|---|
| `_puUp` / `_puDown` | pileup reweighting (`weightPUtruejson_up/down` branches) |
| `_l1PrefireUp` / `_l1PrefireDown` | L1 prefiring weight |
| `_xsecUp` / `_xsecDown` | flat cross-section uncertainty (`FlatXsecSystematic.h`) |
| `_eRecoUp/Down`, `_eIdIsoUp/Down` | electron reco/ID-iso SF |
| `_muIdUp/Down`, `_muIsoUp/Down` | muon ID/iso SF |
| `_tauVsEleUp/Down`, `_tauVsMuUp/Down`, `_tauVsJetUp/Down` | tau discriminator SFs |
| `_eTrigUp/Down`, `_muTrigUp/Down` | trigger SF |
| `_topPtUp` / `_topPtDown` | top-pT reweighting (ttbar MC only) |
| `_zptUp` / `_zptDown` | Z-pT reweighting (DY MC only; built inline in `TightFileProcessor.h`, not in `EventWeights.h`, right after the Z-pT nominal weight is applied) |

**To add a new weight-only systematic**: follow the same "divide out the
nominal factor already baked into `weight`, multiply in the shifted one"
pattern as every entry above, and `push_back` a `{"_suffixUp", shiftedWeight}`
/ `{"_suffixDown", ...}` pair into `weightVariants` — either inside
`appendStandardWeightVariants()` (if it should apply everywhere, like
pileup) or inline in `TightFileProcessor.h`'s event loop (if it's
conditional on the sample, like `_zptUp`/`_zptDown` is on `ctx.isDY`).
Nothing else needs to change — `HistFilling.h::fillEventHistograms()`
iterates whatever's in `weightVariants` automatically, and (for 3-/4-lepton
events) the signal-shape tree loop in `TightFileProcessor.h` (`for (const
auto& variant : weightVariants) fillSignalShapeVariant(...)`) picks it up
too, reusing the nominal reconstruction since weight-only systematics don't
touch lepton momenta.

**To remove one**: delete its `push_back` call(s). No other file needs to
change.

**Category (b) — MET-only shape systematics** (recoil response/resolution,
built inline in `TightFileProcessor.h` from `MEtSysWrapper.h`'s
`computeRecoilSystShift()`):

| Suffix pair | Source |
|---|---|
| `_recoilResponseUp` / `_recoilResponseDown` | recoil-correction response uncertainty |
| `_recoilResolutionUp` / `_recoilResolutionDown` | recoil-correction resolution uncertainty |

Only pushed into `metVariants` when `hasRecoilVariation` is true (DY/WJets
MC with a valid gen boson and a loaded `MEtSys` for that year). Each entry
only re-fills the MET/METphi histograms for every region
(`HistFilling.h::fillMetOnly`), not the full kinematic bundle — cheaper
than a shape systematic, appropriate since only MET itself is shifted, not
lepton momenta. To add another MET-only systematic, push a
`{"_suffix", shiftedPx, shiftedPy}` into `metVariants`.

**Category (c) — shape (lepton-momentum) systematics** (each needs its own
explicit re-histogramming block, since every derived kinematic quantity
must be recomputed):

| Suffix pair | Source |
|---|---|
| `_roccorUp` / `_roccorDown` | Rochester muon-momentum correction uncertainty (only events containing a muon) |
| `_tauESUp` / `_tauESDown` | tau energy-scale uncertainty (only events containing a tau) |

The pattern (see `TightFileProcessor.h` around the `catstr.find('m')` /
`catstr.find('t')` blocks) is: save the original `pt_1..pt_4`, apply the
shift, rebuild `EventKinematics`, fill histograms + (for 3-/4-lepton
events) recompute and fill the signal-shape variant under the shifted
kinematics, then restore the original pT before moving to the next
systematic. **To add a new shape systematic**, copy that four-step pattern
(shift -> rebuild kinematics -> fill -> restore) with your own shift
function. Note the comment in the code: gate the block on the category
string containing the relevant flavor (structural), *not* on whether a
shift was actually computed for a specific lepton slot — gating on
per-event availability instead silently drops the event from the Up/Down
histograms entirely (this exact bug collapsed the tauES/roccor bands to
near-zero before it was caught).

## Where output directories are set

Every driver hardcodes its own EOS output path near the top of the
function body — there's no central config for this, so all four locations
need to be changed together when starting a new production tag:

| File | Line(s) | What it sets |
|---|---|---|
| `DCH_tight.C` | `outDirBase` (~line 54) | histogram output for unchunked/full-file runs |
| `DCH_tight_Chunk.C` | `outDir` (~line 43), `treeOutDir` (~line 45) | per-chunk histogram + mass-tree output (unbatched) |
| `DCH_tight_ChunkBatch.C` | `eosHistChunkDir` / `eosTreeChunkDir` (~lines 60-61) | per-chunk-batch histogram + mass-tree output on EOS (production path) |
| `batch/dch_tight_merge_persample.sub` | `arguments =` line | chunk-source dir and merged-target dir passed to `MergeChunkHists.C` |

In the current snapshot these are all pinned to
`new_hists/.../run2_noFR_metphi_zpt_recoil_roccor_toppt_syst_v2` (chunk
output) and `new_hists/.../run2_noFR_metphi_zpt_recoil_roccor_toppt_syst`
(merged output, no `_v2` — a naming mismatch inherited as-is from the
existing production, not something this repo update tried to "fix"). When
starting a genuinely new production (e.g. after a code change that
invalidates old output), bump the tag in all four places above to a new
name so old and new output can never be mixed on disk.

## Running the pipeline

### 1. Compile check

```bash
root -l -b -q -e '.L DCH_tight_ChunkBatch.C+'
```

If a header changed but the `.so` doesn't visibly rebuild, force it:

```bash
rm -f DCH_tight_ChunkBatch_C.so DCH_tight_ChunkBatch_C.d DCH_tight_ChunkBatch_C_ACLiC_dict_rdict.pcm
root -l -b -q -e '.L DCH_tight_ChunkBatch.C+'
```

ACLiC's header-dependency tracking on a shared, concurrently-accessed AFS
checkout is not fully reliable — a stale `.so` can silently keep running
old code even after the header on disk changed. **Never edit a file this
pipeline includes while a condor batch is in flight**: a job that
(re)compiles mid-batch picks up whatever the header looks like at that
exact moment, silently producing a batch with mixed code versions across
its jobs. Let a submitted batch finish (or hold it) before editing
anything, then force a clean recompile before resubmitting.

### 2. Running over a single sample / single file (troubleshooting)

Interactively, for fast iteration while debugging:

```bash
# one file, whole file, single process
root -l -b -q 'DCH_tight.C+("2018",0,1,1,"TTbar")'
```

Arguments: `(year, firstFile, nFilesToRun, nProc, processFilter, skipData)`.
`processFilter` must match a key in `FileMap.h`'s per-year map (e.g.
`"TTbar"`, `"signal"`, `"data"`) — leave it empty to run every sample.
`nFilesToRun = -1` runs all remaining files from `firstFile` onward;
`nProc > 1` uses `ROOT::TProcessExecutor` to fan out across files.

For a single, small, fast-to-inspect chunk (useful when you only want to
check the output schema, not run a whole file):

```bash
root -l -b -q 'DCH_tight_Chunk.C+("2018",0,0)'   # year, file index, chunk index 0
```

then inspect with `root -l new_hists/chunks/.../hist_<sample>_chunk0.root`
and `.ls`. To find a file's index for a given sample, check the order
`getFileMap(year)` iterates in `filemap/FileMap.h` (it's just the vector
order for that sample's key).

### 3. Condor batch production (the actual production path)

Large samples are too big to run as one file per condor job (a single file
can be tens of millions of events), so production goes through
`DCH_tight_ChunkBatch.C`, which processes a *batch* of consecutive
500k-event chunks in one job, merging them locally on the job's own
scratch disk before uploading a single merged-batch file to EOS — this
keeps the peak EOS footprint far below what staging every raw chunk
individually would cost (staging ~6550 raw chunks peaked at >100GB before
the final merge shrank it back to <20GB, against a ~265GB EOS quota).

```bash
# 1) list every (year, file-index, chunk-index) chunk that exists
root -l -b -q batch/gen_chunklist.C          # writes batch/chunklist.txt
# (needs a valid grid proxy -- it opens every file over xrootd just for GetEntries())

# 2) group consecutive chunks of the same file into batches of 15
python3 batch/gen_batched_chunklist.py       # writes batch/batched_chunklist.txt

# 3) submit
cd batch
condor_submit dch_tight_chunk_batch.sub
```

Each job runs
`run_dch_chunkbatch.sh DCH_tight_ChunkBatch <year> <idx> <batchstart> <numinbatch>`,
which itself does `.L DCH_tight_ChunkBatch.C+` fresh on the shared AFS
checkout (see the compile-check warning above). Batch size (15 chunks/job)
is set by `BATCH_SIZE` in `batch/gen_batched_chunklist.py`; change it there
and regenerate if you want bigger/smaller jobs (bigger = fewer condor jobs
but longer wall-time and more risk per job; smaller = the opposite).

Chunk-batch output lands in
`new_hists/chunks/<tag>/<year>/hist_<sample>_chunk<N>.root` and the
matching `new_trees/chunks/<tag>/<year>/masstree_<sample>_<year>_chunk<N>.root`
(`<N>` here is the batch's own starting chunk index, per the naming
convention `MergeChunkHists.C` already expects — no change needed there
even though each "chunk" file is really a merged batch of 15).

Monitor with `condor_q <cluster>`; jobs commonly go on hold transiently
from AFS token/permission errors (`HoldReason` containing "errno 27" or
"errno 5") — `condor_release <cluster>` clears those. If a job holds
repeatedly, check `batch/logs/*.err` for that job first.

### 4. Merging chunks into per-sample files

`MergeChunkHists.C` is generic: it groups any `<base>_chunk<N>.root` files
in a source directory by `<base>`, sums their histograms, concatenates
(via `TChain::CloneTree`) any `TTree`s they contain, writes one merged file
per base into the target directory, and deletes the consumed chunk files on
success. It doesn't care whether `<base>` starts with `hist_` or
`masstree_` — but **the histogram output (`new_hists/`) and the
signal-shape tree output (`new_trees/`) live in separate directory trees**
(see "Where output directories are set" above), so they need two separate
merge passes, not one:

```bash
cd batch
condor_submit dch_tight_merge_persample.sub          # new_hists/  (hist_<sample>_<year>.root)
condor_submit dch_tight_merge_persample_trees.sub     # new_trees/  (masstree_<sample>_<year>.root)
```

Both read the same `batch/mergelist_persample.txt` (one `year,sample` line
per merge job — regenerate it from whatever samples/years you actually
processed if it's out of date) and both run
`run_dch_merge_chunks_persample.sh`, which takes the chunk-source dir,
merged-target dir, year, sample, and a file-prefix (`hist` or `masstree`)
and calls `MergeChunkHists(chunkDir, targetDir, {"<prefix>_<sample>.root"})`.
A sample/year with no signal-shape chunks at all (e.g. a background sample
that never produced a 3-/4-lepton event) simply merges zero files in the
tree pass — not an error, nothing to do.

Once a merge job succeeds, its consumed chunk files are gone — so
re-running it for a sample that already succeeded finds no chunks left and
does nothing (safe to re-submit, e.g. after fixing a held job elsewhere).

For a one-off/interactive merge of a single sample:

```bash
root -l -b -q 'MergeChunkHists.C("new_hists/chunks/<tag>/2018", "new_hists/<merged-tag>/2018", {"hist_TTTo2L2Nu_2018.root"})'
root -l -b -q 'MergeChunkHists.C("new_trees/chunks/<tag>/2018", "new_trees/<merged-tag>/2018", {"masstree_TTTo2L2Nu_2018.root"})'
```

(third argument is the list of *source* file base names, i.e. what the
unchunked run would have produced — `MergeChunkHists.C` looks for
`<base minus .root>_chunk*.root` under the source dir and writes `<base>`
under the target dir.)

### 5. Downstream stacking

Interactively, for one year:

```bash
root -l -b -q 'Stackhist.C+("2018")'
root -l -b -q 'Stackhist_multiplicity.C+("2018")'
```

reads the merged `hist_<sample>_<year>.root` files and produces data/MC
stack plots per region (and, for the `_multiplicity` variant, split by
lepton multiplicity).

Both scripts loop over many histogram variables, which is slow enough to be
worth condor-parallelizing the same way production is: `batch/filelist_stackhist.txt`
has one `year,first` line per variable index, and

```bash
cd batch
condor_submit stackhist.sub                # -> Stackhist.C
condor_submit stackhist_multiplicity.sub   # -> Stackhist_multiplicity.C
```

runs one job per `(year, variable index)` via `run_stackhist.sh`, which does
`.L <Stackhist|Stackhist_multiplicity>.C+` and calls
`<Script>("<year>",<first>,1)` (i.e. exactly one variable per job). Both
submission files read the *input* histogram directory from `Stackhist.C`'s/
`Stackhist_multiplicity.C`'s own hardcoded `INPUT_DIR` (see "Where output
directories are set" above) — update that constant, not the `.sub` file, to
point at a different production tag.

## Notes

- `DCH_tauFR.C` / `DCH_tauFR_SF.C` (tau fake-rate *derivation*) are not
  part of this snapshot yet.
- `sf_uncertainty/`, `recoil_studies/`, `systematics/results/` are
  derivation/validation code and cached results, not runtime dependencies
  of `DCH_tight*.C` — kept here for provenance, not required to run the
  pipeline itself.
- `systematics_plots/` macros (`PlotZPtWeight.C`, `PlotTopPtWeight.C`,
  `PlotRoccorSF.C`, `PlotTauSF.C`, `PlotRecoilSF.C`, `PlotLeptonSF.C`) are
  standalone — each reads the same correction source `DCH_modules/`
  applies and plots its SF (+ Up/Down band) vs. the relevant kinematic
  variable, in a shared house style (`SFPlotStyle.h`). Useful as a
  cross-check that a correction module reads the right payload before
  trusting it in a full production run.
