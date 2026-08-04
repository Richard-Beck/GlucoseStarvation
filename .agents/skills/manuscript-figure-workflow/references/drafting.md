# Figure Drafting Workflow

## Purpose

Turn an approved ideation proposal, approved redraft figure package, or explicit
user instruction into reviewable draft figure material.

Approval to draft does not mean the user carefully endorsed every agent
interpretation. Drafts must match the supplied prompt, plan, and revision
context and the scientific role of each requested panel.

Read `.agents/references/manuscript_figure_style.md` before creating draft
figures or reviewing visual quality.

## Inputs

Required:

- A work-package id or redraft figure-package id.
- An approved ideation file, approved redraft plan, or explicit user instruction
  to draft.
- A planning source, usually `docs/manuscript_plan_v3.txt`,
  `docs/manuscript_plan.txt`, or a redraft root containing `redraft_plan.md`.

Optional:

- Supplied feedback-manager handoff or tracking pointer, when present.
- Existing draft outputs to diagnose or revise.
- Legacy feedback files or sibling package review notes supplied as revision
  context.

## Output Location

Use the most local review root implied by the request:

- Legacy work-package workflow:
  `<work_package_root>/<WP_ID>/drafting/`
- Redraft figure-package workflow:
  `<redraft_root>/figure_generation/<package_id>/drafting/`

Write the applicable outputs:

- `initial_subpanels/<subpanel_id>/`: source-generated draft options and notes.
- `refined_subpanels/`: selected or refined PNG subpanels.
- `final_figures/conservative/`, `final_figures/recommended/`, and optional
  `final_figures/exploratory/`: assembled draft figures plus folder-level
  `legend.md`.
- `review_report.html`: self-contained human-review report with embedded PNGs,
  directive status, recommended drafts first, legends, feedback history, review
  context, caveats, prior-code fidelity, and provenance summary.
- `drafting_panels.md`: purpose, generation, selection rationale, provenance,
  caveats, visual QC, and revision history.
- `feedback_manager_context.md`: optional concise pointer when a managed-feedback
  handoff was supplied.
- `feedback_intake.md`: required consumer ledger of prompt/feedback inputs,
  scope decisions, directive ids, assignments, and final output coverage.
- `report_manifest.csv`: required coverage manifest for every reviewable PNG.
- `prior_panel_disposition.csv`: prior panels marked preserved, targeted fix,
  moved, replaced, or dropped.
- `prior_code_fidelity.csv`: provenance check for every affected prior panel.
- `prior_code/<panel_id>/` and `code_diffs/<panel_id>.diff` for inherited panels
  with identifiable prior reviewed code.
- `not_drafted.md` if any approved option is skipped.

## Directive Intake

Record in `feedback_intake.md`:

- every prompt/plan/feedback-manager handoff source used, with enough identifier
  detail to find it again
- package scope decisions such as `draft`, `abort`, `subsume_into_<package>`,
  `defer_to_integration`, `defer_to_results_text`, or
  `blocked_pending_investigation`
- concrete user directives, each with a stable `directive_id`
- cross-version and cross-package implications
- each directive's planned output, report section, subagent, or explicit
  non-drafting disposition
- each directive's final status: `addressed`, `partially_addressed`, `blocked`,
  `dropped`, or `deferred`, plus the output path or blocker

When a managed-feedback handoff is supplied, its raw feedback takes precedence
over ideation interpretations. If feedback contradicts a candidate idea, follow
the feedback and document the disposition.

## Prior Code Fidelity Check

For every affected prior figure/subpanel, classify the drafting path before
plotting:

- `inherited_preserve`: regenerate the prior reviewed panel unchanged except for
  local-path relocation.
- `inherited_targeted_fix`: modify only to satisfy explicit feedback or prompt
  instruction.
- `inherited_move`: reuse in a different position with substantive content
  preserved.
- `inherited_replace`: intentionally replace because feedback, prompt, or plan
  requires it.
- `novel_no_prior`: no prior reviewed panel/code exists or user requested a novel
  panel.
- `blocked_prior_missing`: prior reviewed output or code should exist but cannot
  be identified.

Write `prior_code_fidelity.csv` before generating reviewable PNGs, with columns:

```text
panel_id,inheritance_class,prior_png_path,prior_code_path,copied_baseline_code_path,
active_local_script_path,diff_path,allowed_change_directive_ids,fidelity_status,blocker
```

For inherited preserve, targeted-fix, or move panels, copy prior reviewed code
into local `prior_code/<panel_id>/` before adaptation. Adapt that copy into the
active local script under `scripts/`, then generate `code_diffs/<panel_id>.diff`.
Each substantive diff must be justified by a feedback or prompt directive id, a
mechanical local-path/output-path relocation, or a blocker.

Do not source old manuscript draft scripts, polishing scripts, prior draft PNGs,
or ad hoc figure-generation folders as a shortcut. Those files may be inspected
and copied into the local package, but the local drafting package must contain
the reproducible source used for every presented inherited panel.

## Subagent Brief Contract

When using drafting subagents, give each a short written brief containing:

- source paths for ideation, feedback-manager handoff if applicable, planning
  source, and prior panels
- inheritance class for assigned panels
- directive ids and feedback excerpts governing the assigned work
- expected drafted outputs, including variants, investigation tasks, or
  abort/defer decisions
- local-script requirement and output folder
- copied-prior-code and diff requirements for inherited panels
- request for directive-by-directive coverage mapping outputs, skipped outputs,
  blockers, and caveats

Integrate returned coverage into `feedback_intake.md`, `drafting_panels.md`, and
`review_report.html`, plus `feedback_manager_context.md` when present. Make
blocked or partially addressed directives visible near the top of the report.

## Review Report

After each initial drafting pass or feedback revision, write or update
`review_report.html` as the sole human-review surface. It must be static and
self-contained with every reviewable PNG visually present.

Use this order:

1. Decision summary and reviewer action needed.
2. Directive-to-output status table.
3. Recommended final drafts.
4. Blockers and missing required outputs.
5. Conservative alternatives.
6. Exploratory alternatives.
7. Refined and initial subpanels, collapsed by default unless primary.
8. Full provenance appendix.

Include a prior-code fidelity table near the top report summary. Before
generating the report, write `report_manifest.csv` for every PNG under
`initial_subpanels/`, `refined_subpanels/`, and `final_figures/`.

Use `<skill_dir>/scripts/review_report_template.R` as the starting template.
Copy it into the local drafting package, usually as
`scripts/make_<package>_review_report.R`, adapt the CONFIG block, and run it with
the project-approved R runner. Preserve the directive-status table,
recommended-first order, visible blocker section, report coverage checks, and
collapsed raw galleries.

## Workflow

1. Resolve the requested package, prompt, plan, data, and available figure/text
   context from supplied paths and local artifacts.
2. Read `.agents/references/manuscript_figure_style.md`.
3. Identify planning source, ideation files, supplied feedback-manager handoff
   if any, and
   the most recent user-reviewed or polished version of each affected panel.
4. Write `feedback_intake.md` before scope decisions are settled. If a managed
   handoff was supplied, also write `feedback_manager_context.md`.
5. Write `prior_panel_disposition.csv` and `prior_code_fidelity.csv`.
6. For inherited redrafting, copy prior code locally, adapt it under `scripts/`,
   and generate code diffs.
7. Use scoped drafting subagents when useful and available.
8. Draft only options needed to answer feedback/prompt directives and approved exploratory
   ideas. Record skipped options in `not_drafted.md`.
9. Write local generation scripts, notes, manifests, legends, `drafting_panels.md`,
   revision-context notes, and `review_report.html`.
10. Visually inspect final assembled PNGs against the style reference, the last
    user-reviewed versions, and every feedback/prompt directive.
11. Review the generated HTML itself for missing images, hidden blockers, missing
    directive rows, or absent prior-code fidelity coverage.

## Constraints

- Run R scripts through the project-approved R runner.
- Do not start expensive refits, segmentation reruns, classifier retraining, or
  long jobs unless explicitly approved.
- Treat unrequested changes to data, filtering, statistical summary, visual
  encoding, layout, panel order, panel role, or interpretation as unacceptable
  unless the panel is explicitly novel or exploratory.
- In `drafting_panels.md`, explicitly assess visual-QC items from
  `.agents/references/manuscript_figure_style.md`, plus whether any summary
  display is inappropriate where raw data was requested.
