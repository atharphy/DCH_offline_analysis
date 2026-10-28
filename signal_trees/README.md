# signal_trees

Signal-region tree production for the 3lep2tau / DCH analysis: same selection, corrections and systematic
variants as the tight histogram production, but writing one flat tree per region (signal-shape variables, including
`met` and `metphi`) instead of histograms. No fake-rate code is part of this folder.

Everything here is self-contained and relative to the repository root, except the skim locations in
`filemap/FileMap.h`, which point to the remote storage the skims live on.

## Layout

Outputs do not live in this folder: they go to `output/{hists,trees,plots}/` at the repository root (ignored by git).
The tree chain writes to `output/trees/`; `hists/` and `plots/` are there for the downstream steps.

| Path | Content |
|---|---|
| `DCH_tight_TreeOnly_ChunkBatch.C` | driver: one file, a batch of 500k-event chunks, merged locally and copied to `<repo>/output/trees/chunks/` |
| `DCH_modules/` | event selection, weights, corrections and systematic plan used by the chain (`TightTreeOnlyFileProcessor.h` is the processor) |
| `Mass_reco_modules/` | signal-shape selection, histogram/tree booking and the mass reconstruction (`massreco/`) |
| `include/`, `roccor/` | branch declarations, kinematics, cross-sections, Rochester code |
| `filemap/` | `FileMap.h` (skim lists) and `FileMapHTQCD.h` (HT-binned QCD instead of FlatPt) |
| `data/zpt_studies/<year>/` | Z pT weights (DY) |
| `batch/` | HTCondor submit files and run scripts |

Payloads the chain loads from elsewhere in the repository: `charge_flip/DD_flip_rates/`, `roccor/*.txt` (repository root),
`json_files/met_<year>.json.gz`, `HTT-utilities/RecoilCorrections/data/`.

## Run

From the repository root (inside a CMSSW_13_0_10 `src` area, with `correctionlib` on the include and library path as in
`batch/run_*.sh`):

```
root -l -b -q -e '.L signal_trees/DCH_tight_TreeOnly_ChunkBatch.C+' \
     -e 'DCH_tight_TreeOnly_ChunkBatch("2018",0,0,15)'
```

arguments: year, file index in the file map, first chunk, number of chunks in the batch.
`DCH_FILEMAP=qcdht` switches the QCD samples to the HT-binned set and writes to its own folder.

Batch: submit from `signal_trees/batch/` (the submit files take the repository root from `$ENV(PWD)/../..`):

```
condor_submit dch_tight_treeonly_chunkbatch.sub        # chunk batches, list in batched_chunklist.txt
condor_submit dch_tight_treeonly_merge_persample.sub   # merge the chunks per sample, list in mergelist_treeonly_persample.txt
```

The `*_qcdht*` files do the same for the HT-binned QCD set; `gen_treeonly_qcdht_list.C` writes its chunk list.
The chunk and merge lists are produced from the file map and are not checked in.

Output: `<repo>/output/trees/chunks/<tag>/<year>/masstree_<sample>_chunk<N>.root`, merged into `<repo>/output/trees/<tag>/<year>/`.

## Notes

- The Rochester payload, recoil/MEt-sys payloads and charge-flip rates are read relative to the working directory, so
  always run from the repository root.
- The submit files use `$ENV(HOME)/.grid_proxy/x509up.pem` as the grid proxy.
