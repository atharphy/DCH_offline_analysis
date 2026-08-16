# Reconstructed Z-pT reweighting

`DeriveZPtWeights.C` derives separate `ee` and `mm` shape weights from the
already-weighted `h_zPt_<channel>_DYCR_0tau` histograms.

The numerator is background-subtracted Data. The denominator is the `DY`
sample only. Both are normalized to unit area before division, so this changes
the DY Z-pT shape without intentionally changing its inclusive normalization.
`DY10_50` is treated as a background in the Z-window control region.
For Data, `ee` reads only SingleElectron/EGamma files and `mm` reads only
SingleMuon files, preventing overlap between primary datasets.

Run from the repository root:

```bash
root -l -b -q 'zpt_studies/DeriveZPtWeights.C+("2018")'
root -l -b -q 'zpt_studies/DeriveZPtWeights.C+("Run2")'
```

The default input directory is the MET-XY-corrected histogram production. The
default output base is
`/eos/user/a/atahmad/DCH_offline_analysis/zpt_studies`. The input and output
locations can be overridden with the second and third arguments:

```bash
root -l -b -q 'zpt_studies/DeriveZPtWeights.C+("2018","/path/to/hists","/path/to/output")'
```

For each requested year, the macro writes:

- `ZPtWeights_<year>.root`, including nominal/up/down weights;
- `ZPtWeight_ee.png` and `ZPtWeight_mm.png` diagnostics;
- human-readable `ZPtWeight_ee.txt` and `ZPtWeight_mm.txt` tables.

The final bin includes all overflow events. Apply the nominal correction only
to DY simulation by multiplying the nominal event weight by the content of
`h_zpt_weight_ee` or `h_zpt_weight_mm` at the reconstructed Z pT. Values below
the first edge should use the first bin and values above the final edge should
use the last bin.

The correction uses 1 GeV bins from 0 to 20 GeV, followed by edges at 30, 40,
50, 70, 100, 150, 200, 300, and 500 GeV.
