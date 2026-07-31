# Polishing Gate Contract

Use a small JSON file to make validation deterministic. Paths are resolved relative to the directory where the validator is run unless `project_root` is set.

## Minimal Example

```json
{
  "wp_id": "WP4",
  "polish_root": "figures/manuscript_draft_v3/WP4/polishing",
  "feedback_file": "agent-dev/manuscript_work_packages/WP4/drafting/feedback.md",
  "draft_doc": "agent-dev/manuscript_work_packages/WP4/drafting/drafting_panels.md",
  "panel_map": "figures/manuscript_draft_v3/WP4/polishing/panel_map.csv",
  "source_pngs": [
    "agent-dev/manuscript_work_packages/WP4/drafting/refined_subpanels/F4A_parameter_effect_intervals.png"
  ],
  "approved_raster_subpanels": [
    "figures/user-approved-subpanels/model_family_schematic.png"
  ],
  "figure_script": "figures/manuscript_draft_v3/WP4/polishing/scripts/polish_figures.R",
  "subpanel_dimensions": "figures/manuscript_draft_v3/WP4/polishing/layout/subpanel_dimensions.csv",
  "layout_optimizer_script": ".agents/skills/manuscript-polishing-gate/scripts/optimize_panel_layout.R",
  "layout_optimizer_command": "scripts/agentRrunner.sh .agents/skills/manuscript-polishing-gate/scripts/optimize_panel_layout.R --input figures/manuscript_draft_v3/WP4/polishing/layout/subpanel_dimensions.csv --output-dir figures/manuscript_draft_v3/WP4/polishing/layout --target-width 7 --max-height 9.25",
  "layout_plan": "figures/manuscript_draft_v3/WP4/polishing/layout/layout_plan.csv",
  "layout_report": "figures/manuscript_draft_v3/WP4/polishing/layout/layout_report.md",
  "layout_qc_files": [
    "figures/manuscript_draft_v3/WP4/polishing/layout/Figure_4_layout_preview.png",
    "figures/manuscript_draft_v3/WP4/polishing/layout/Figure_S7_layout_preview.png"
  ],
  "expected_outputs": [
    "figures/manuscript_draft_v3/WP4/polishing/final_images/figure_4.png",
    "figures/manuscript_draft_v3/WP4/polishing/final_images/figure_s7.png"
  ],
  "legend_files": [
    "figures/manuscript_draft_v3/WP4/polishing/legend.md"
  ],
  "provenance_table": "figures/manuscript_draft_v3/WP4/polishing/provenance.csv",
  "output_manifest": "figures/manuscript_draft_v3/WP4/polishing/manifest.csv",
  "polishing_notes": "figures/manuscript_draft_v3/WP4/polishing/notes.md",
  "visual_qc_file": "figures/manuscript_draft_v3/WP4/polishing/visual_qc.md",
  "validation_report": "figures/manuscript_draft_v3/WP4/polishing/validation_report.json",
  "report_dir": "figures/manuscript_draft_v3/WP4/polishing"
}
```

## Fields

- `wp_id`: required work-package id, such as `WP2`.
- `project_root`: optional root path for resolving relative paths.
- `polish_root`: required by default; all polish-generated scripts, reports, subpanels, final images, legends, manifests, provenance, notes, and validation outputs should live under this directory.
- `require_polish_root`: optional boolean, default true. Set false only for legacy validation or an explicitly approved exception.
- `feedback_file`: required; must exist and be non-empty.
- `approval_required`: optional boolean. If true, feedback text must contain one approval marker.
- `approval_markers`: optional list of accepted approval strings. Defaults include `APPROVED`, `approved`, `good enough to proceed`, and `proceed`.
- `draft_doc`: optional but recommended; the current drafting documentation.
- `panel_map`: required by default; CSV crosswalk that defines the exact final panel set and connects each panel to approved draft identity and regeneration strategy.
- `require_panel_contract`: optional boolean, default true when regenerated subpanels are required. Set false only for legacy validation or a documented, user-approved exception.
- `source_files`: optional list of non-image source files that must exist before polishing.
- `source_pngs`: optional list of approved draft/reference PNG panels that must exist and be readable before polishing. These are visual or approval references, not final subpanel inputs.
- `approved_raster_subpanels`: optional list of exact immutable subpanel assets supplied by the user under `figures/user-approved-subpanels/`. These are the only raster image files that may be loaded as final-composite subpanel inputs, and each listed asset must correspond to an individual panel-map row. Agents must treat this directory as read-only and must not crop, trim, recolor, resample, rewrite, or otherwise edit these files.
- `figure_script`: required by default for postflight; single R script that regenerates every displayed subpanel as an audit export and assembles final composites from the same in-memory panel constructors.
- `subpanel_dimensions`: required by default for postflight; CSV dimensions table written by `figure_script`.
- `require_regenerated_subpanels`: optional boolean, default true. Set false only for a documented, user-approved exception.
- `allow_split_scripts`: optional boolean, default false. Set true only for a documented exception allowing legacy `subpanel_script` plus `assembly_script` fields instead of `figure_script`.
- `allow_raster_assembly`: deprecated legacy field. Do not use it for new polishing contracts; it is too broad. Use `approved_raster_subpanels` only for exact immutable user-supplied subpanel assets under `figures/user-approved-subpanels/`.
- `allow_in_figure_titles`: optional boolean, default false. Set true only for a documented, user-approved exception allowing figure-level titles/subtitles/captions in final PNGs.
- `layout_optimizer_script`: optional; defaults to `.agents/skills/manuscript-polishing-gate/scripts/optimize_panel_layout.R` when omitted.
- `layout_optimizer_command`: required by default for postflight; exact command used to run the layout optimizer.
- `layout_plan`: required by default for postflight; CSV layout/scaling plan written by the optimizer.
- `layout_report`: optional but recommended; Markdown summary written by the optimizer.
- `layout_qc_files`: required by default for postflight; list of optimizer preview or assembly-QC artifacts that verify panel order and placement.
- `require_layout_qc`: optional boolean, default true. Set false only for an explicitly intermediate validation pass or documented exception.
- `expected_outputs`: required for postflight; polished PNG outputs that must exist and be readable.
- `legend_files`: required by default for postflight; separate text/Markdown figure legend files.
- `require_legend_files`: optional boolean, default true. Set false only for an explicitly intermediate validation pass.
- `provenance_table`: required by default for postflight; CSV table with one row per displayed subpanel.
- `require_provenance_table`: optional boolean, default true. Set false only for an explicitly intermediate validation pass.
- `output_manifest`: recommended; CSV manifest written by polishing.
- `polishing_notes`: recommended; Markdown notes written by polishing.
- `visual_qc_file`: required by default for postflight; Markdown record of direct visual inspection of every final PNG composite.
- `require_visual_qc`: optional boolean, default true. Set false only for an explicitly intermediate validation pass or documented exception.
- `validation_report`: optional path for saved validation output.
- `report_dir`: optional report/export directory.
- `output_dir`: optional figure-output directory.
- `optional_outputs`: optional list. Missing files are warnings, not errors.
- `check_manifest_paths`: optional boolean, default true. Checks CSV manifest cells that look like local file paths.
- `require_project_map_decision`: optional boolean, default true. Warns if notes do not mention a project-map decision.

## Panel Map Columns

The panel map is the preflight contract for panel identity. The validator requires these columns when `require_panel_contract` is true:

- `figure`: polished figure id, such as `Figure 4`.
- `panel`: lowercase panel label, such as `a`.
- `approved_source`: approved draft panel image, draft documentation row, or source script path used to identify the intended panel. Semicolon-separated paths are allowed.
- `approved_generator`: script, notebook, or command that generated the approved draft panel, when available.
- `intended_content`: concise description of what the panel must show.
- `regeneration_strategy`: how polishing will regenerate the same panel identity from data/report exports.

The final `subpanel_dimensions`, `layout_plan`, and `provenance_table` panel keys must match the panel map exactly. If a panel cannot be traced to the approved draft identity, stop and ask for direction instead of substituting nearby content.

## Subpanel Dimension Table Columns

The validator requires these columns when `require_regenerated_subpanels` is true:

- `figure`: polished figure id, such as `Figure 4`.
- `panel`: lowercase panel label, such as `a`.
- `subpanel_png`: regenerated subpanel PNG path.
- `width_px`: regenerated subpanel width in pixels.
- `height_px`: regenerated subpanel height in pixels.
- `width_in`: regenerated subpanel width in inches.
- `height_in`: regenerated subpanel height in inches.

The `figure_script` should write this table after writing the audit subpanel PNGs. Pixel dimensions should match the PNG headers.

## Layout Plan Columns

The bundled optimizer writes the required layout columns:

- `figure`, `panel`, `subpanel_png`
- `x_in`, `y_in`, `width_in`, `height_in`
- `sx`, `sy`
- `x_npc`, `y_npc`, `width_npc`, `height_npc`
- `layout_width_in`, `layout_height_in`

Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `x_npc`, `y_npc`, `width_npc`, and `height_npc` directly with layout systems such as `cowplot::draw_plot()`. For custom `grid` layouts, place the viewport center at `x_npc + width_npc / 2` and `y_npc + height_npc / 2`; do not invert `y_npc`.

## Figure Script Pattern

The default polishing script is one file, usually `scripts/polish_figures.R`, run twice:

```text
scripts/agentRrunner.sh <polish_root>/scripts/polish_figures.R --phase subpanels
scripts/agentRrunner.sh .agents/skills/manuscript-polishing-gate/scripts/optimize_panel_layout.R --input <polish_root>/layout/subpanel_dimensions.csv --output-dir <polish_root>/layout --target-width 7 --max-height 9.25
scripts/agentRrunner.sh <polish_root>/scripts/polish_figures.R --phase final
```

The script should define reusable panel constructors. The subpanel phase saves audit PNGs and dimensions. The final phase rebuilds the same panel objects and composes them directly. Do not re-read polish-generated audit PNGs to build final composites. The only raster source allowed in final assembly is an exact immutable subpanel listed in `approved_raster_subpanels` under `figures/user-approved-subpanels/`; load it directly as one subpanel object and combine it with regenerated plot objects.

Before writing final composites, remove figure-level titles/subtitles/captions from the final panel objects and layout wrapper. Avoid `ggtitle()`, `labs(title=...)`, `labs(subtitle=...)`, `labs(caption=...)`, `plot_annotation(title=..., subtitle=..., caption=...)`, and figure-level `draw_label()` prose in final PNGs. Put figure titles and explanatory prose in `legend_files`, not images.

## Visual QC File

The visual QC file documents visual inspection of the rendered final PNG images themselves. File existence, dimensions, logs, CSVs, and validation output are not substitutes. It should contain one row or bullet per final PNG with checks for:

- no visible figure title/subtitle/caption/header text
- no clipping of panel labels, axis text, legends, colorbars, or plotted data
- acceptable spacing, gutters, margins, and alignment
- readable text and symbols at intended print size
- correct panel order and labels

If any check fails, revise `polish_figures.R`, rerun the needed phases, inspect again, and record the rerun. Do not pass postflight with unresolved significant visual defects unless the user explicitly accepts the caveat.

## Provenance Table Columns

The validator requires these columns when `provenance_table` is present and regenerated subpanels are required:

- `figure`: polished figure id, such as `Figure 4`.
- `panel`: lowercase panel label, such as `a`.
- `subpanel_image`: regenerated subpanel PNG used for the panel, or the exact `figures/user-approved-subpanels/` source path for an approved immutable raster subpanel.
- `generator`: script, notebook, or workflow that generated the subpanel, usually `polish_figures.R`.
- `command`: command needed to regenerate the subpanel and final artifact.
- `data_inputs`: semicolon-separated direct data/report-export inputs.
- `layout_plan`: layout/scaling plan used to place the panel.
- `output_image`: polished figure image containing the panel.
- `notes`: caveats or context-only dependencies.

Use `context_inputs` as an optional extra column for files read for interpretation but not directly plotted.
For an approved immutable raster subpanel, keep a normal provenance row and add an `approved_raster_source` column with the same `figures/user-approved-subpanels/` path. Do not use copied-artifact or whole-figure raster-assembly exceptions for new polishing work.

## Error Policy

The validator exits `1` when it reports any `ERROR`. Treat that as a blocking validation failure. It exits `0` when it reports only `WARN` or `OK` items.
