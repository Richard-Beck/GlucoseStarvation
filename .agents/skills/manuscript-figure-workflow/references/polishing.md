# Figure Polishing Workflow

## Purpose

Promote approved manuscript figure draft material into polished, validated figure
outputs while preserving source provenance. Stop before manuscript-wide final
integration unless the user asks for integration.

## Polishing Rules

1. Confirm the requested work package or figure package, figure scope, and
   approved draft inputs.
2. Keep polishing under the lead agent. Optional subagents may assist with
   bounded inventory or independent QC only when permitted and useful;
   subagent availability never gates polishing.
3. Do not redo ideation or drafting. Use approved draft inputs and any
   feedback-manager handoff supplied for polishing.
4. Do not launch expensive refits, segmentation reruns, classifier retraining, or
   long jobs. If polishing reveals a scientific/input blocker, document it and
   stop.
5. Read `.agents/references/manuscript_figure_style.md` and apply it to polished
   figure images and visual QC. 
6. Run R scripts through the project-approved R runner.
7. Polishing owns final subpanel generation and layout. Use one R script as the
   polishing entrypoint; do not split the default workflow into one script that
   exports PNG panels and another that assembles those PNGs.
8. Before regenerating panels, create a panel-identity map that links every final
    figure/panel label to approved draft panel, intended content, draft
    generator, and regeneration strategy.
9. Regenerate by adapting or calling the approved draft-generation path when it
    exists. Do not substitute simplified plots from convenient summary tables
    unless the user explicitly approves the changed panel identity.
10. Subpanel PNGs are audit exports, not sources for final composites. In final
    assembly, rebuild panel plot/grob objects in memory and follow the global
    raster policy in `SKILL.md`.
11. Treat the layout optimizer's `x_npc` and `y_npc` as lower-left coordinates.
    Do not invert `y_npc` during assembly.
12. Inspect every final PNG composite visually before postflight validation.

## Optional Subagent Assignments

The lead agent owns the polishing contract,
`scripts/polish_figures.R`, final implementation decisions, validation reruns,
and final status. When delegation is permitted, useful focused assignments are:

- Panel/provenance inventory: inventory approved draft figures, source scripts,
  direct data/report-export inputs, manifests, prior-code fidelity records, and
  regeneration risks.
- Independent final PNG QC: inspect rendered final PNGs against the style
  reference, panel map, and feedback handoff; report clipping, spacing, label
  order, readability, raster defects, and panel-identity mismatches.

For multiple work packages, polish sequentially unless the user explicitly asks
for parallel polishing. 

## Output Layout

Create one `polish_root` and keep all polish-generated files under it:

```text
<manuscript_draft_or_wp>/polishing/
  contract.json
  panel_map.csv
  feedback_manager_context.md   # optional supplied managed-feedback handoff
  legend.md                    # optional existing/linkage artifact
  manifest.csv
  provenance.csv
  figure_rebuild_manifest.tsv
  notes.md
  visual_qc.md
  validation_report.json
  scripts/
    polish_figures.R
  subpanels/
  layout/
    subpanel_dimensions.csv
    layout_plan.csv
    layout_report.md
    *_layout_preview.png
  final_images/
```

## Output Object Boundary

- `panel_map.csv` is a polishing-local identity map. It defines the intended
  final panel set before regeneration and feeds `FigureIndexRow`; it is not a
  shared handoff object.
- `manifest.csv` is the polishing-local projection of `FigureIndexRow` fields
  that integration needs: source-local names, source roots, panel labels,
  polished final image paths/hashes, source fields, and subpanel audit paths.
  Polishing does not assign final manuscript numbering.
- `provenance.csv` is a validator-local projection of `FigureIndexRow` panel
  source fields. Keep its required columns aligned with
  `references/polishing_contract.md`.
- `figure_rebuild_manifest.tsv` follows `RebuildRecord` with
  `stage = polishing`.
- `visual_qc.md` records `VisualQCRecord` entries or bullet-equivalent checks
  from direct PNG inspection.
- `validation_report.json` records `ValidationSummary` for the polishing stage.
- `contract.json`, `subpanel_dimensions.csv`, `layout_plan.csv`, layout reports,
  layout previews, optional legend links, and feedback-manager context are
  polishing-local implementation or validation artifacts.

## Polished Figure Requirements

- Apply `.agents/references/manuscript_figure_style.md`.
- Strip figure-level title, subtitle, and caption calls from reused draft code,
  including `ggtitle()`, `labs(title=...)`, `labs(subtitle=...)`,
  `labs(caption=...)`, `patchwork::plot_annotation(title=..., subtitle=...,
  caption=...)`, and figure-level `cowplot::draw_label()` prose.
- Record integration-ready source, image, hash, and provenance metadata for each
  final panel using `FigureIndexRow` fields.
- Record one `RebuildRecord` row per polished whole-figure output.
- If a subpanel cannot be regenerated from data/report exports, treat it as a
  blocker. Follow the global raster policy in `SKILL.md` for user-approved
  raster inputs.
- Record existing legend or legend-source paths when available, but do not make
  legend prose a required polishing output. Do not modify legends. 

## Workflow

1. Build a polishing contract.
   - Use `references/polishing_contract.md` for JSON fields.
   - Write `panel_map.csv` first, with one row per displayed final subpanel.
   - Include the package id, feedback-manager handoff if applicable,
     `polish_root`, panel map, selected draft/reference inputs, one R figure
     script, dimension table, layout optimizer command/output, layout QC files,
     expected polished outputs, optional legend/source links, manifest,
     provenance, rebuild manifest, notes, visual QC, and validation report.
   - Save the contract inside `polish_root`.

2. Run preflight validation before editing figures.
   - Command:
     `python3 <skill_dir>/scripts/validate_polishing.py <contract.json> --phase preflight`
   - Treat `ERROR` items as blockers. Fix missing or malformed inputs before
     polishing; do not bypass validation with package-local approval fields.

3. Write the single polishing script and run its subpanel phase.
   - Write `scripts/polish_figures.R` under `polish_root`.
   - The script must define panel constructors used for audit subpanel exports
     and final composites.
   - The subpanel phase must regenerate every final displayed subpanel listed in
     `panel_map.csv` from direct data/report-export inputs and follow the
     global raster policy in `SKILL.md` for user-approved raster inputs.
   - Run:
     `<r_runner> <polish_root>/scripts/polish_figures.R --phase subpanels`

4. Optimize the figure layout before assembly.
   - Run:
     `<r_runner> <skill_dir>/scripts/optimize_panel_layout.R --input <polish_root>/layout/subpanel_dimensions.csv --output-dir <polish_root>/layout --target-width 7 --max-height 9.25`
   - Read `layout/layout_plan.csv`, `layout/layout_report.md`, and optimizer
     preview PNGs. Preview labels must be in intended reading order.
   - If the optimizer recommends nontrivial scaling, revise dimensions in
     `polish_figures.R`, rerun subpanels, and rerun the optimizer.

5. Run the final phase from the same script.
   - Run:
     `<r_runner> <polish_root>/scripts/polish_figures.R --phase final`
   - The final phase must rebuild panel plot/grob objects in memory and assemble
     them directly. Do not re-read polish-generated `subpanels/*.png` for final
     composites.
   - Use optimizer `x_npc`, `y_npc`, `width_npc`, and `height_npc` as lower-left
     coordinates.
   - Write final PNG composites under `final_images/`, compute hashes, and write
     manifest, provenance, rebuild manifest, and notes.

6. Visually inspect final composite PNGs and revise if needed.
   - Check for visible titles/subtitles/captions, clipping, overlapping elements,
     awkward spacing, unreadable text, incorrect panel order, inconsistent panel
     labels, and cropped legends/colorbars.
   - If significant issues are found, edit `polish_figures.R`, rerun needed
     phases, inspect again, and record the rerun.
   - Write `visual_qc.md` with one row or bullet per final PNG.

7. Run postflight validation after outputs are written.
   - Command:
     `python3 <skill_dir>/scripts/validate_polishing.py <contract.json> --phase postflight --write-report <validation-report.json>`
   - A nonzero exit means polishing is not complete. Fix issues and rerun unless
     the user explicitly accepts the unresolved blocker.

8. Finish with a short status listing polished figure files, manifest, notes,
    validation report, feedback-manager handoff/evidence status when applicable,
    and whether validation passed.

## Output Expectations

A completed polishing task should normally produce:

- one polish root containing all polish-generated artifacts
- a polishing-local panel map tying every final panel to approved draft identity
- one R script that regenerates displayed subpanels as audit artifacts and
  assembles final composites from the same in-memory panel constructors
- regenerated subpanel PNGs and a subpanel dimension table
- layout optimizer outputs and previews
- polished PNG composite figure(s)
- integration-ready `FigureIndexRow` metadata
- `RebuildRecord` manifest
- polishing notes
- feedback-manager context when a managed-feedback handoff was supplied
- visual QC notes from direct PNG inspection
- validator report
