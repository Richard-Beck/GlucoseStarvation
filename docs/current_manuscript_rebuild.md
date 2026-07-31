# Current manuscript rebuild boundary

The cleanup keep set targets `manuscripts/20260731T101031_v03`, the latest
assembled manuscript at cleanup time. The exact repository-relative keep roots
are in `docs/current_manuscript_rebuild_keep_roots.txt`.

## Meaning of ground-up

For this repository, ground-up means rebuilding the current analytical objects,
23 figures, and manuscript render from retained primary or curated scientific
sources, the transformations that remain reproducible, and the current code.
It does not mean reacquiring microscopy data or recreating external
software/model artifacts that the current Methods package itself records as
external or unresolved.

A `terminal` row in the Methods provenance graph closes the Methods-facing
trace; it does not establish that the row is a primary-data boundary. Every
scientific terminal row therefore requires an independent disposition in
`docs/current_manuscript_rebuild_source_boundaries.tsv`. The cleanup audit
requires that review to cover the current terminal set exactly and verifies
that every locally retained source root exists and is protected by the keep
manifest.

The canonical graph is
`agent-dev/manuscript_methods/20260730_v03_methods_provenance_reconstruction/locked_provenance_table.md`.
At cleanup time its deterministic checks reported 271 graph rows, 234 unchanged
hashed local objects, no changed or missing locks, and an exact match for all 79
panel endpoints. Thirty-seven semantic collection, external, or unresolved
nodes have no stored file hash. Their concrete repository collections are kept
conservatively, including `all_raw`, the current image-processing releases,
classifier-training inputs, and the active model releases. Four external or
unresolved software/model identities remain limitations already disclosed by
the Methods package.

The extracellular-glucose source package is retained as the complete
`data/glucose` tree. Its `raw` directory contains the four XLSX workbooks as
received, its `processed` directory contains the CSV inputs consumed by the
model, and its README records that assay/count matching, reformatting, and
calibration extraction were performed manually. The processed CSVs remain
locked terminal objects in the Methods graph, while the source-boundary review
prevents that graph status from excluding their primary workbooks.
This compact source package is versioned explicitly despite the repository's
broad ignore policy for other large `data` collections.

## Rebuild and verification

From the repository root, regenerate all approved figure pixels into a new,
empty directory:

```bash
manuscripts/20260731T101031_v03/rebuild/figures/run_all_figures.sh OUTPUT_DIR
```

Render the manuscript:

```bash
manuscripts/20260731T101031_v03/rebuild/python_runner.sh \
  manuscripts/20260731T101031_v03/rebuild/manuscript/build_manuscript_html.py
```

Validate the package and a figure replay:

```bash
manuscripts/20260731T101031_v03/rebuild/python_runner.sh \
  manuscripts/20260731T101031_v03/rebuild/manuscript/validate_manuscript.py \
  --rebuild-dir OUTPUT_DIR
manuscripts/20260731T101031_v03/rebuild/python_runner.sh \
  manuscripts/20260731T101031_v03/rebuild/validate_assembly.py
```

Audit the filesystem against the keep-root policy:

```bash
python3 scripts/audit_cleanup_keep_roots.py
```

The audit must report no unexpected paths outside `junk/`, complete coverage of
the current scientific terminal review, and no missing retained source roots.
The `junk/` tree is ignored by Git and preserves each relocated path relative
to the repository root. It is intentionally local and recoverable by moving an
entry back to its original relative path.
