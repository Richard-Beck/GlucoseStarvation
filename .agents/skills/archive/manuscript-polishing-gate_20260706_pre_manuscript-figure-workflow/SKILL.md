---
name: archived-manuscript-polishing-gate
description: "Archived historical copy of manuscript-polishing-gate, superseded by manuscript-figure-workflow on 2026-07-06. Use only when explicitly auditing or comparing the pre-promotion manuscript polishing workflow."
---

# Manuscript Polishing Gate

## Purpose

Use this skill to execute the polishing gate for a reviewed manuscript work package. The gate promotes approved draft material into polished, validated figure outputs while preserving source provenance and stopping before manuscript-wide final integration unless the user asks for integration.

## Gate Rules

1. Confirm the requested work package, figure scope, and approved draft inputs.
2. Use subagents for polishing. At the start of the task, confirm that subagent deployment is available. If subagents are not available, immediately abort the polishing attempt before editing files or running polishing code, and explicitly ask the user whether to enable/permit subagents or continue without them.
3. Do not redo ideation or drafting. Use the approved feedback and selected draft panels as inputs.
4. Do not launch expensive refits, segmentation reruns, classifier retraining, or long jobs. If polishing reveals a scientific/input blocker, document it and stop.
5. Read `docs/manuscript_plan_v3.txt` when scope or numbering is unclear.
6. Do not automatically read `docs/project_map.txt`. Read it only if the polishing pass adds or changes a maintained entrypoint, canonical output path, or workflow description. Always record a project-map decision in the polishing notes.
7. Read `.agents/references/manuscript_figure_style.md` and apply it to polished figure images, legends, and visual QC.
8. Run R scripts through `scripts/agentRrunner.sh`.
9. Polishing owns final subpanel generation and layout. Use one R script as the polishing entrypoint; do not split the default workflow into one script that exports PNG panels and another that assembles those PNGs.
10. Before regenerating panels, create a panel-identity map that links every final figure/panel label to the approved draft panel, intended content, draft generator, and regeneration strategy. The final panel set must match this map exactly.
11. Regenerate by adapting or calling the approved draft-generation path whenever it exists. Do not substitute simplified plots from convenient summary tables unless the user explicitly approves the changed panel identity.
12. Subpanel PNGs are audit exports, not sources for final composites. In the final phase, rebuild panel plot/grob objects in memory and compose those objects with `patchwork`, `cowplot`, or equivalent R tooling. The only raster exception is an individual subpanel asset explicitly placed by the user under `figures/user-approved-subpanels/`; load that asset directly as one immutable subpanel object, record it in the panel map/provenance, and do not crop, trim, recolor, resample, rewrite, or otherwise edit the raster pixels. Do not use polish-generated `subpanels/*.png`, draft PNG crops, or copied refined panels as final-composite inputs.
13. Treat the layout optimizer's `x_npc`/`y_npc` as lower-left coordinates. Do not invert `y_npc` during assembly; use it directly with `cowplot::draw_plot()` or with a `grid::viewport()` centered at `y_npc + height_npc / 2`.
14. Inspect every final PNG composite visually before postflight validation. If clipping, awkward spacing, unreadable text, overlapping elements, wrong panel order, or title/subtitle text is visible, revise `polish_figures.R` and rerun the needed phases until the issue is fixed or documented as a user-approved exception.

## Subagent Deployment

The lead agent owns the polishing contract, `panel_map.csv`, `scripts/polish_figures.R`, final implementation decisions, validation reruns, and final status. Use subagents in parallel within each work package for focused review tasks that reduce context pressure without fragmenting the executable polishing path.

Recommended subagent assignment for each work package:

- Panel/provenance inventory subagent: inventory approved draft figures, candidate source scripts, direct data/report-export inputs, existing manifests, prior-code fidelity records, and regeneration risks; return a concise panel-by-panel handoff for the lead agent to convert into `panel_map.csv` and provenance rows.
- User-feedback integration subagent: read raw feedback, approval notes, critic/reviewer reports, redrafter notes, and closeout handoffs; record every relevant user-feedback item that has been integrated into the work package, every item explicitly deferred to polishing, and any rejected or superseded feedback that must not re-enter the polished output.
- Independent final PNG QC subagent: after the lead agent generates final composites, inspect the rendered PNGs directly against the shared style reference, approved panel map, and feedback handoff; report title/caption leakage, clipping, spacing, label order, readability, raster defects, and any mismatch between final images and approved panel identity.

For requests that include multiple work packages, polish work packages sequentially. For each work package, deploy the parallel subagents above, complete that work package through postflight validation and notes, then compact context before starting the next work package. Do not run multiple work-package polishing gates concurrently unless the user explicitly overrides this sequencing.

## Output Layout

Create one `polish_root` for the gate and keep all polish-generated files under it. A typical structure is:

```text
<manuscript_draft_or_wp>/polishing/
  contract.json
  panel_map.csv
  legend.md
  manifest.csv
  provenance.csv
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

Use stable names, but keep this single-root pattern unless the user requests a different storage layout.

## Polished Figure Requirements

Apply these requirements to polished figure images:

- Apply `.agents/references/manuscript_figure_style.md` as the shared figure-formatting standard for figure images, panel labels, legends, color/annotation definitions, and visual QC.
- Strip figure-level title calls from reused draft code before final export, including `ggtitle()`, `labs(title=...)`, `labs(subtitle=...)`, `labs(caption=...)`, `patchwork::plot_annotation(title=..., subtitle=..., caption=...)`, and figure-level `cowplot::draw_label()` text.
- Treat the shared style reference as mandatory for final polished outputs unless the contract, notes, and provenance record a user-approved exception.

Apply these provenance requirements to every polished figure:

- Record each subpanel as a separate provenance row keyed by figure id and panel label.
- For each subpanel, record the regenerated subpanel image path, generating script or command, direct input data files/report exports, layout plan, upstream draft documentation or feedback, and the final polished output that contains it.
- Preserve enough command detail that another agent can rerun or audit the subpanel without guessing which workflow produced it.
- Distinguish direct inputs used to draw the panel from contextual files used only for interpretation.
- If a subpanel cannot be regenerated from data/report exports, treat it as a polishing blocker unless the user has supplied that exact immutable subpanel under `figures/user-approved-subpanels/`.
- Treat missing provenance for any displayed subpanel as a polishing blocker unless the user explicitly accepts the limitation.

Apply these requirements to separate figure legends:

- Follow the legend requirements in `.agents/references/manuscript_figure_style.md`.

## Workflow

1. Build a polishing contract.
   - Use `references/polishing_gate_contract.md` for the JSON fields.
   - Write `panel_map.csv` first, with one row per displayed final subpanel and columns documented in `references/polishing_gate_contract.md`.
   - Include the WP id, feedback file, `polish_root`, panel map, selected draft/reference inputs, one R figure script, subpanel dimension table, layout optimizer command/output, layout QC files, expected polished outputs, separate legend file(s), provenance table, manifest path, notes path, visual QC path, and validation report.
   - Save the contract inside `polish_root`.

2. Run preflight validation before editing figures.
   - Command: `python3 .agents/skills/manuscript-polishing-gate/scripts/validate_polishing_gate.py <contract.json> --phase preflight`
   - Treat `ERROR` items as blockers. Fix missing or malformed inputs before polishing, or ask the user if the blocker reflects an intentional exception.
   - Treat `WARN` items as review notes, not blockers.

3. Write the single polishing script and run its subpanel phase.
   - Write `scripts/polish_figures.R` under `polish_root`.
   - The script must define the panel constructors used for both audit subpanel exports and final composites.
   - The subpanel phase must regenerate every final displayed subpanel listed in `panel_map.csv` from direct data/report-export inputs, not from draft PNG crops or copied refined panels. The only exception is an exact user-approved raster subpanel under `figures/user-approved-subpanels/`; record its dimensions and provenance without modifying the source file.
   - Prefer adapting or calling the scripts that produced the approved draft panels. If the approved draft panel source cannot be found or rerun, stop and document a blocker instead of inventing replacement content.
   - It must write audit subpanel images under `subpanels/`.
   - It must write `layout/subpanel_dimensions.csv` with one row per displayed subpanel and columns for `figure`, `panel`, `subpanel_png`, `width_px`, `height_px`, `width_in`, and `height_in`.
   - Run it with `scripts/agentRrunner.sh <polish_root>/scripts/polish_figures.R --phase subpanels`.

4. Optimize the figure layout before assembly.
   - Run the bundled optimizer on the subpanel dimension table, for example:
     `scripts/agentRrunner.sh .agents/skills/manuscript-polishing-gate/scripts/optimize_panel_layout.R --input <polish_root>/layout/subpanel_dimensions.csv --output-dir <polish_root>/layout --target-width 7 --max-height 9.25`
   - Read `layout/layout_plan.csv`, `layout/layout_report.md`, and the optimizer preview PNGs. The preview must show panel labels in the intended reading order before assembly begins.
   - If the optimizer recommends nontrivial `sx` or `sy` changes, revise the dimensions in `polish_figures.R`, rerun the subpanel phase, and rerun the optimizer before assembling final figures.
   - Keep the optimizer outputs as polish artifacts.

5. Run the final phase from the same script.
   - Run `scripts/agentRrunner.sh <polish_root>/scripts/polish_figures.R --phase final`.
   - The final phase must rebuild the same panel plot/grob objects in memory and assemble those objects. Do not re-read `subpanels/*.png` to make final composites. For an exact user-approved raster subpanel under `figures/user-approved-subpanels/`, load only that source PNG as a single immutable subpanel grob and combine it with the regenerated plot objects.
   - Prefer `patchwork` for rectangular ggplot layouts or `cowplot::ggdraw()` plus `cowplot::draw_plot()` for custom positions. Use raw `grid` coordinate math only when necessary and document the coordinate convention in notes.
   - When using the optimizer plan, `x_npc` and `y_npc` are the lower-left corner of each panel in normalized parent coordinates. For `grid`, use viewport centers `x_npc + width_npc / 2` and `y_npc + height_npc / 2`; do not use `1 - y_npc - height_npc / 2`.
   - Remove all figure-level title/subtitle/caption code before saving final composites. For reused draft plot objects, explicitly blank `plot.title`, `plot.subtitle`, and `plot.caption` where needed and avoid `plot_annotation(title=...)` in final images.
   - Write final PNG composites under `final_images/`.
   - Use stable, non-draft filenames where possible.
   - Preserve a manifest mapping polished outputs to regenerated subpanels, layout plans, and feedback.
   - Write a panel-level provenance table before postflight validation.
   - Write separate figure legend file(s).
   - Write concise polishing notes covering inputs, panel-map decisions, commands, feedback applied, layout decisions, caveats, validation status, and the project-map decision.

6. Visually inspect final composite PNGs and revise if needed.
   - Visually inspect the rendered final PNG images themselves. Do not rely only on file existence, dimensions, logs, CSVs, or validation output.
   - Look for visible figure titles/subtitles/captions, clipping of labels or plot content, overlapping elements, awkward or excessive spacing, unreadable text, incorrect panel order, inconsistent panel labels, and cropped legends/colorbars.
   - If any significant issue is found, edit `polish_figures.R`, rerun `--phase subpanels` and the optimizer if dimensions changed, rerun `--phase final`, then inspect again.
   - Write `visual_qc.md` with one row or bullet per final PNG, recording pass/fail for title removal, clipping, spacing/layout, and readability, plus any reruns or accepted caveats.

7. Run postflight validation after outputs are written.
   - Command: `python3 .agents/skills/manuscript-polishing-gate/scripts/validate_polishing_gate.py <contract.json> --phase postflight --write-report <validation-report.json>`
   - A nonzero exit means polishing is not complete. Fix the issue and rerun validation, unless the user explicitly accepts the unresolved blocker.

8. Finish with a short status.
   - List polished figure files, manifest, notes, and validation report.
   - State whether validation passed and whether `docs/project_map.txt` was checked or updated.

## Output Expectations

A completed polishing gate should normally produce:

- a single polish root containing all polish-generated artifacts
- a panel map tying every final panel to the approved draft panel/source and intended content
- one R script that regenerates all displayed subpanels as audit artifacts and assembles final composites from the same in-memory panel constructors, except for any exact immutable subpanel assets supplied under `figures/user-approved-subpanels/`
- regenerated subpanel PNGs plus a subpanel dimension table
- layout optimizer outputs and scaling recommendations
- polished PNG composite figure(s)
- separate figure legend file(s)
- a panel-level provenance table with scripts, commands, and data dependencies
- a polished figure manifest CSV
- polishing notes in Markdown
- a visual QC Markdown file documenting final PNG inspection and any reruns
- a validation report from `validate_polishing_gate.py`

Keep this gate distinct from final integration. Do not create or renumber the entire manuscript figure set unless explicitly asked.
