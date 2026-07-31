# Figure rebuild

Run from the project root:

```bash
manuscripts/20260731T101031_v03/rebuild/figures/run_all_figures.sh OUTPUT_DIR
```

`OUTPUT_DIR` must be empty. The runner copies all seven assembly-local polishing
packages to temporary working directories, regenerates their figures from the
declared data inputs, normalizes the 23 manuscript PNGs into `OUTPUT_DIR`, and
requires every regenerated PNG to match the approved checksum in
`figure_rebuild_manifest.tsv`. It never writes to upstream polishing packages.

Figure 3a is the sole approved immutable-raster exception: its model-family
schematic is read from `figures/user-approved-raster-figures/`. All statistical
plots are regenerated locally. Microscopy shown in the SUM-159 package is read
from source image data by the polishing code rather than copied as a whole
manuscript subpanel.
