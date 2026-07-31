---
name: manuscript-figure-workflow
description: "Determine and execute the next manuscript figure task from a prompt, 
plan, user feedback, data, and optional preexisting figures, legends, or manuscript text. 
Use when Codex needs to choose among and run figure ideation, draft generation/revision, 
polishing, figure-set integration, or subpanel semantic interpretation."
---

# Manuscript Figure Workflow

## Purpose

Use this skill to intake a prompt, plan, user feedback, data, and optionally one
or more preexisting figures, legends, reports, or manuscript text; determine the
next manuscript figure task or workflow; execute that task; and prepare the
resulting outputs for review or handoff. Detailed task procedures live in
optional reference files so agents load only the instructions needed for the
current turn.

## Single-Turn Workflow

Always follow this workflow within a turn:

1. Orient and determine inputs.
   - Identify the user's requested outcome, prompt/plan constraints, data inputs,
     feedback handoff if any, and any preexisting figures, legends, reports, or
     manuscript text.
   - Decide whether the task is a figure idea, draft, polish, figure-set
     integration, visual interpretation, or a small connected sequence of these.

2. Determine workflow.
   - Select the task module or modules needed for the current turn.
   - Load only the relevant reference files from `references/`.
   - Use more than one module when the user's request naturally requires a short
     connected sequence, such as polishing selected draft figures before
     figure-set integration.
   - If the request requires expensive computation, model refits, segmentation
     reruns, classifier retraining, GPU use, large RAM, or expected runtime over
     5 minutes, use the major-analysis or batch-compute workflow before figure
     work.

3. Execute workflow.
   - Follow the selected module instructions.
   - Preserve figure identity, data scope, and feedback intent unless the user or
     supplied plan clearly asks for a change.
   - If execution reveals that the selected module is wrong for the task, switch
     modules and record why in the local notes or review report.

4. Prepare outputs.
   - Write the review, figure, interpretation, manifest/index, validation, or
     handoff artifacts required by the executed module.
   - Make outputs self-contained enough for review: include paths to source
     scripts/data, visual QC notes for generated images, feedback handoff context
     when applicable, and unresolved blockers.

## Task Modules

Use these modules as task procedures, not as a fixed sequence:

- Figure ideation: greenfield figure planning, targeted revision brainstorming,
  candidate panel options, or a user-review checkpoint before plotting.
- Figure drafting: reviewable draft figures, subpanel options,
  figure-generation scripts, draft revisions from feedback, prior-code fidelity
  checks, and draft review reports.
- Figure polishing: final regeneration, layout optimization, final PNG
  composites, integration-ready source/image/hash metadata, rebuild records,
  validation, and visual QC for already selected draft material.
- Figure set integration: normalized manuscript figure sets, final image
  packaging, canonical figure index assembly, rebuild/regeneration scripts,
  semantic interpretation links, optional derived lineage audit, and packaging
  validation.
- Subpanel semantic interpretation: one fresh-context visual inspection per
  subpanel, followed by main-agent source clarification and reconciliation for
  later integration, manuscript-legend-writing, Methods, or evidence review.

## Shared Inputs

Common inputs include:

- A work-package id such as `WP3` or redraft figure-package id such as
  `figure_4_s8`, `FG3_wp4_posterior`, or similar.
- A planning source such as `docs/manuscript_plan.txt`,
  `docs/manuscript_plan_v3.txt`, or a redraft root containing `redraft_plan.md`.
- User feedback or approval context.
- Existing figures, polished outputs, draft panels, manifests, provenance tables,
  semantic interpretations, and figure legends relevant to the requested package.
- Integrated figure roots with `figure_set_manifest.csv`, final images, claim
  graph, semantic interpretation paths, or existing legend files when they
  clarify figure-package context. Use
  `manuscript-legend-writing` for renderer-facing integrated legend prose.

## Output Roots

Use the most local review root implied by the request:

- Legacy work-package workflow:
  `agent-dev/manuscript_work_packages/<WP_ID>/<task>/`
- Redraft figure-package workflow:
  `<redraft_root>/figure_generation/<package_id>/<task>/`
- Figure polishing root:
  keep all polish-generated files under one `polishing/` root selected by the
  polishing contract.
- Semantic interpretation outputs:
  use the user-specified path or a local `subpanel_interpretations/` folder near
  the relevant figure package or integration root.
- Figure-set integration outputs:
  default to a timestamped or versioned integration root, usually
  `agent-dev/manuscript_integration/<run_id>/`.

Write self-contained HTML review reports for ideation and drafting. For
polishing, write final PNGs, integration-ready metadata, rebuild records, notes,
visual QC, and validation output under the polish root. For semantic
interpretation, keep visual/source/audit material separate from manuscript-facing
interpretation text. 

## Canonical Figures Handoff Contract

When handing off completed figure sets for manuscript consumption, ensure the
Figure set integration contract is followed. Its manifest is the canonical
figure index and may require source regeneration of unmodified figures.


## Shared Figure Rules

- Read `.agents/references/manuscript_figure_style.md` before drafting or
  polishing figure images.
- Do not launch expensive refits, segmentation reruns, classifier retraining, or
  long jobs unless explicitly approved.
- Prefer existing analysis data and helpers under `data/`, `R/`, and maintained
  report exports over ad hoc reconstruction.
- For inherited panels, preserve the most recent user-reviewed panel identity
  unless feedback or the approved plan requires a targeted change or
  replacement.
- Presented draft, polished, and integrated panels must be generated by local scripts in the
  task package, not by copying old PNGs.
- Final manuscript figure composites must be assembled from regenerated
  plot/grob objects. Before embedding an immutable raster, discover the
  repository's canonical centralized directory named
  `user-approved-raster-figures`; do not assume its parent path. Exactly one
  such repository-level directory must exist. If none or more than one exists,
  explain the ambiguity to the user and block raster use.
- Determine whether the canonical directory contains the correct raster for the
  intended panel. Compare the candidate with the user-designated or most recent
  user-reviewed reference by checksum first. A checksum match establishes exact
  identity. A checksum mismatch is not by itself a blocker: directly inspect
  both images and accept the canonical candidate when their displayed content
  is a visual match. Record the compared paths, checksums, and whether identity
  was established by checksum or visual inspection. If no canonical candidate
  can be matched, explain the issue to the user and block.
- The only raster image files that may be embedded into manuscript-visible
  composites are the matched immutable assets in that discovered canonical
  directory and declared in the workflow's raster input fields. Load the
  canonical file directly and do not crop, trim, recolor, resample, rewrite, or
  replace its pixels. Treat the directory as read-only. Any other
  raster-derived displayed panel is a blocker.
- Visually inspect reviewable or final PNG outputs. File existence, dimensions,
  logs, and validation output are not substitutes for looking at the image.
- Record blockers where the requested panel cannot be regenerated, provenance is
  missing, feedback conflicts, or the output would misrepresent the data.

## Bundled Helpers

Resolve `<skill_dir>` to this installed skill folder:

```text
.agents/skills/manuscript-figure-workflow
```

- `scripts/review_report_template.R`: drafting review-report template.
- `scripts/optimize_panel_layout.R`: polishing layout optimizer.
- `scripts/validate_polishing.py`: polishing contract validator.
- `references/polishing_contract.md`: polishing contract schema and column
  conventions.

Use this installed path in generated local scripts or notes that need to call
bundled helper scripts.

## Optional Task Modules

Additional task modules are documented in `references/`. Read the module files
needed for the current turn; do not load unrelated references.

name: Figure ideation
description: Create or revise reviewable manuscript figure ideation packages,
including revision context, existing-panel disposition, candidate ideas,
and a self-contained HTML review report.
source: references/ideation.md

name: Figure drafting
description: Generate or revise reviewable draft figure material from approved
ideation, an approved redraft plan, or explicit user instruction, while
preserving directive coverage, prior-code fidelity, local regeneration,
and draft review reports.
source: references/drafting.md

name: Figure polishing
description: Promote approved draft material into regenerated,
layout-optimized, validated polished manuscript figures with panel maps,
integration-ready metadata, rebuild records, final PNG visual QC, and validator
output.
source: references/polishing.md

name: Figure set integration
description: Assemble polished figure-package outputs into a normalized,
auditable manuscript figure set with final images, a canonical figure index,
semantic interpretation links, optional derived lineage audit, rebuild scripts,
and packaging validation.
source: references/figure_set_integration.md

name: Subpanel semantic interpretation
description: Interpret one or more rendered manuscript subpanels using fresh
per-subpanel visual inspectors followed by grouped main-agent source
clarification, with explicit limits on provenance, manuscript evidence, and
visual claims.
source: references/semantic_interpretation.md
