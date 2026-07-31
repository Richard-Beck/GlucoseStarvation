# Rebuild

Manuscript render:

```bash
manuscripts/20260731T101031_v03/rebuild/python_runner.sh \
  manuscripts/20260731T101031_v03/rebuild/manuscript/build_manuscript_html.py
```

Complete figure replay into an empty destination:

```bash
manuscripts/20260731T101031_v03/rebuild/figures/run_all_figures.sh OUTPUT_DIR
```

The figure runner copies its seven local polishing packages to temporary work
directories and does not modify upstream packages. See the nested README files
for validation commands and scope. `python_runner.sh` requires Python 3.9 or
newer, honors `MANUSCRIPT_PYTHON`, and otherwise falls back to the established
HPC `cpose` environment when the system `python3` is too old.
