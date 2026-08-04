# Figure Polishing Contract

Use a small JSON file to make polishing validation deterministic. Paths are
resolved relative to the directory where the validator is run unless
`project_root` is set.

## Minimal Example

```json
{
  "wp_id": "WP4",
  "polish_root": "<polish_root>",
  "draft_doc": "<draft_root>/drafting_panels.md",
  "panel_map": "<polish_root>/panel_map.csv",
  "source_pngs": [
    "<draft_root>/refined_subpanels/F4A_parameter_effect_intervals.png"
  ],
  "approved_raster_subpanels": [
    "<discovered-canonical-approved-raster-directory>/model_family_schematic.png"
  ],
  "figure_script": "<polish_root>/scripts/polish_figures.R",
  "subpanel_dimensions": "<polish_root>/layout/subpanel_dimensions.csv",
  "layout_optimizer_script": "<skill_dir>/scripts/optimize_panel_layout.R",
  "layout_optimizer_command": "<r_runner> <skill_dir>/scripts/optimize_panel_layout.R --input <polish_root>/layout/subpanel_dimensions.csv --output-dir <polish_root>/layout --target-width 7 --max-height 9.25",
  "layout_plan": "<polish_root>/layout/layout_plan.csv",
  "layout_report": "<polish_root>/layout/layout_report.md",
  "layout_qc_files": [
    "<polish_root>/layout/Figure_4_layout_preview.png"
  ],
  "expected_outputs": [
    "<polish_root>/final_images/figure_4.png"
  ],
  "provenance_table": "<polish_root>/provenance.csv",
  "output_manifest": "<polish_root>/manifest.csv",
  "rebuild_manifest": "<polish_root>/figure_rebuild_manifest.tsv",
  "byte_identity_report": "<polish_root>/figure_byte_identity_report.tsv",
  "polishing_notes": "<polish_root>/notes.md",
  "visual_qc_file": "<polish_root>/visual_qc.md",
  "validation_report": "<polish_root>/validation_report.json",
  "report_dir": "<polish_root>"
}
```

## Fields

- `wp_id`: required work-package id, such as `WP2`.
- `project_root`: optional root path for resolving relative paths.
- `polish_root`: required by default; all polish-generated files should live
  here.
- `require_polish_root`: optional boolean, default true.
- `feedback_manager_context`: optional path to a short local note pointing to a
  supplied managed-feedback handoff.
- `draft_doc`: optional but recommended current drafting documentation.
- `panel_map`: required by default; CSV map defining the exact final panel set.
- `require_panel_contract`: optional boolean, default true when regenerated
  subpanels are required.
- `source_files`: optional non-image source files required before polishing.
- `source_pngs`: optional approved draft/reference PNG panels used as visual or
  approval references, not final subpanel inputs.
- `approved_raster_subpanels`: optional exact immutable subpanel assets that
  satisfy the global raster policy in `SKILL.md`. The validator checks these
  paths before allowing final-composite raster loading. Use the actual path
  discovered for the repository; the placeholder in the example is not a
  literal directory name.
- `figure_script`: required by default for postflight; single R script that
  regenerates displayed subpanels and assembles final composites.
- `subpanel_dimensions`: required by default for postflight.
- `require_regenerated_subpanels`: optional boolean, default true.
- `allow_split_scripts`: optional boolean, default false.
- `allow_in_figure_titles`: optional boolean, default false.
- `layout_optimizer_script`: optional; defaults to this skill's
  `scripts/optimize_panel_layout.R` when omitted by local convention.
- `layout_optimizer_command`: required by default for postflight.
- `layout_plan`: required by default for postflight.
- `layout_report`: optional but recommended.
- `layout_qc_files`: required by default for postflight.
- `require_layout_qc`: optional boolean, default true.
- `expected_outputs`: required for postflight.
- `legend_files`: optional existing legend or legend-source files to validate if
  the polishing package already carries them. Figure polishing should not require
  manuscript-facing legend prose; use `manuscript-legend-writing` for that.
- `require_legend_files`: optional boolean, default false. Set true only for a
  legacy or explicitly requested package-local legend requirement.
- `provenance_table`: required by default for postflight.
- `require_provenance_table`: optional boolean, default true.
- `output_manifest`: recommended.
- `rebuild_manifest`: optional but recommended package-level rebuild manifest
  following `RebuildRecord` in `references/shared_figure_objects.md`, with
  `stage = polishing`.
- `byte_identity_report`: optional derived checksum or byte-identity report for
  final polished whole-figure outputs.
- `polishing_notes`: recommended.
- `visual_qc_file`: required by default for postflight.
- `require_visual_qc`: optional boolean, default true.
- `validation_report`: optional path for saved validation output.
- `report_dir`: optional report/export directory.
- `output_dir`: optional figure-output directory.
- `optional_outputs`: optional list. Missing files are warnings, not errors.
- `check_manifest_paths`: optional boolean, default true.
- `require_project_map_decision`: optional boolean, default true.

## Panel Map Columns

The validator requires these columns when `require_panel_contract` is true:

- `figure`: polished figure id, such as `Figure 4`.
- `panel`: lowercase panel label, such as `a`.
- `approved_source`: approved draft panel image, draft documentation row, or
  source script path used to identify the intended panel.
- `approved_generator`: script, notebook, or command that generated the approved
  draft panel when available.
- `intended_content`: concise description of what the panel must show.
- `regeneration_strategy`: how polishing will regenerate the same panel identity
  from data/report exports.

The final `subpanel_dimensions`, `layout_plan`, and `provenance_table` panel keys
must match the panel map exactly.

## Subpanel Dimension Table Columns

Required when `require_regenerated_subpanels` is true:

- `figure`
- `panel`
- `subpanel_png`
- `width_px`
- `height_px`
- `width_in`
- `height_in`

Pixel dimensions should match the PNG headers.

## Layout Plan Columns

The bundled optimizer writes:

- `figure`, `panel`, `subpanel_png`
- `x_in`, `y_in`, `width_in`, `height_in`
- `sx`, `sy`
- `x_npc`, `y_npc`, `width_npc`, `height_npc`
- `layout_width_in`, `layout_height_in`

Coordinate convention: `x_in`, `y_in`, `x_npc`, and `y_npc` are lower-left panel
origins. For custom `grid` layouts, place the viewport center at
`x_npc + width_npc / 2` and `y_npc + height_npc / 2`; do not invert `y_npc`.

## Figure Script Pattern

The default polishing script is one file, usually `scripts/polish_figures.R`,
run in phases:

```text
<r_runner> <polish_root>/scripts/polish_figures.R --phase subpanels
<r_runner> <skill_dir>/scripts/optimize_panel_layout.R --input <polish_root>/layout/subpanel_dimensions.csv --output-dir <polish_root>/layout --target-width 7 --max-height 9.25
<r_runner> <polish_root>/scripts/polish_figures.R --phase final
```

The script should define reusable panel constructors. The subpanel phase saves
audit PNGs and dimensions. The final phase rebuilds the same panel objects and
composes them directly. Do not re-read polish-generated audit PNGs to build
final composites. Apply the global raster policy in `SKILL.md` for any
final-composite raster input.

Before writing final composites, remove figure-level titles, subtitles, and
captions from final panel objects and layout wrappers.

## Rebuild Outputs

`rebuild_manifest` is the polishing-stage instance of `RebuildRecord` from
`references/shared_figure_objects.md`. Use `stage = polishing` and one row per
polished whole-figure output.

`byte_identity_report` is optional derived validation detail. Do not treat it as
a primary contract object. If exact byte identity cannot be achieved, record the
issue as a validation blocker and in polishing notes.

## Visual QC File

The visual QC file documents direct visual inspection of rendered final PNGs. It
should contain one row or bullet per final PNG with checks for:

- no visible figure title, subtitle, caption, or header text
- no clipping of panel labels, axis text, legends, colorbars, or plotted data
- acceptable spacing, gutters, margins, and alignment
- readable text and symbols at intended print size
- correct panel order and labels

If any check fails, revise `polish_figures.R`, rerun the needed phases, inspect
again, and record the rerun. Do not pass postflight with unresolved significant
visual defects unless the user explicitly accepts the caveat.

## Provenance Table Columns

Required when `provenance_table` is present and regenerated subpanels are
required:

- `figure`
- `panel`
- `subpanel_image`
- `generator`
- `command`
- `data_inputs`
- `layout_plan`
- `output_image`
- `notes`

Use `context_inputs` as an optional extra column for files read for
interpretation but not directly plotted. For a raster input allowed by the
global raster policy in `SKILL.md`, add `approved_raster_source` with the same
path.

## Error Policy

The validator exits `1` when it reports any `ERROR`. Treat that as a blocking
validation failure. It exits `0` when it reports only `WARN` or `OK` items.
