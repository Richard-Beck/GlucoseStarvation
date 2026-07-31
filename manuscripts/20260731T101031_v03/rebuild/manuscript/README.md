# Manuscript rebuild

From the project root, run:

```bash
manuscripts/20260731T101031_v03/rebuild/python_runner.sh \
  manuscripts/20260731T101031_v03/rebuild/manuscript/build_manuscript_html.py
manuscripts/20260731T101031_v03/rebuild/python_runner.sh \
  manuscripts/20260731T101031_v03/rebuild/manuscript/validate_manuscript.py \
  --rebuild-dir TEMPORARY_REPLAY_OUTPUT
```

The renderer reads assembly-local A/I/D and References files, five Results
sidecars, canonical Methods, the integrated legend file, the 79-row figure
manifest, and all 23 final PNGs. It produces a self-contained HTML review draft
with embedded figure data. It does not alter manuscript prose or legends.
