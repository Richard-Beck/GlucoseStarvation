---
name: archived-manuscript-figure-drafting
description: "Archived historical copy of manuscript-figure-drafting, superseded by manuscript-figure-workflow on 2026-07-06. Use only when explicitly auditing or comparing the pre-promotion figure drafting workflow."
---

# Manuscript Figure Drafting

## Purpose

Use this skill to turn an approved ideation proposal or approved redraft figure package into reviewable draft figure material.

Approval-only feedback files mean permission to draft; they do not mean the user carefully endorsed every agent interpretation. Drafts must continue to match the spirit of the most relevant direct user feedback and the scientific role of each requested panel.

For figure formatting and visual-QC expectations, read `.agents/references/manuscript_figure_style.md`.

## Inputs

Required:

- A work-package id such as `WP3` or a redraft figure-package id such as `figure_4_s8`.
- An approved ideation file, approved redraft plan, or explicit user instruction to draft.
- A planning source, usually `docs/manuscript_plan_v3.txt` or a redraft root containing `redraft_plan.md`.

Optional:

- Existing drafting feedback file.
- Existing draft outputs to diagnose or revise.
- Direct user-feedback files or notes relevant to the affected figure(s), including sibling package feedback when figure placement, abort/subsume decisions, or panel reuse are cross-package issues.

## Output Location

Use the most local review root implied by the request:

- Legacy work-package workflow: `agent-dev/manuscript_work_packages/<WP_ID>/drafting/`.
- Redraft figure-package workflow: `<redraft_root>/figure_generation/<package_id>/drafting/`.

Write:

- `initial_subpanels/<subpanel_id>/`: source-generated draft options plus `notes.md`.
- `refined_subpanels/`: selected/refined PNG subpanels.
- `final_figures/conservative/`: assembled draft figure PNGs rebuilt from local drafting scripts, plus a folder-level `legend.md`.
- `final_figures/recommended/`: assembled draft figure PNGs following the strongest synthesis from ideation and draft review, plus a folder-level `legend.md`.
- `final_figures/exploratory/`: optional assembled draft figure PNGs for more aggressive redesigns, plus a folder-level `legend.md` when present.
- `review_report.html`: self-contained human-review report with embedded PNGs, a directive-to-output status table, recommended drafts first, figure legends, panel-specific feedback history, review context, caveats, and provenance summary.
- `drafting_panels.md`: purpose, generation, selection rationale, provenance, caveats, and revision history.
- `feedback.md`: user-editable feedback file.
- `feedback_history/`: immutable snapshots of each actionable feedback/redrafting round, when there is more than one round or when revising existing draft outputs.
- `feedback_intake.md`: required intake ledger of direct feedback sources, package scope decisions, subagent assignments, and output/report coverage expectations.
- `report_manifest.csv`: required report coverage manifest for every reviewable PNG drafted under `initial_subpanels/`, `refined_subpanels/`, and `final_figures/`.
- `prior_panel_disposition.csv`: each relevant prior panel/subpanel marked `preserved`, `targeted_fix`, `moved`, `replaced`, or `dropped`, with rationale.
- `prior_code_fidelity.csv`: required provenance gate for every affected prior panel/subpanel, mapping prior reviewed PNG/code to copied local baseline code, active local code, diff artifacts, allowed user-feedback directives, and fidelity status.
- `prior_code/<panel_id>/`: required copied baseline code for every inherited panel when prior reviewed code exists.
- `code_diffs/<panel_id>.diff`: required code comparison for every inherited panel with copied baseline code and active local code.
- `not_drafted.md`: required if any ideation subpanel option is not drafted.

## Feedback Contract

Create `feedback.md` when absent. It should point reviewers to `review_report.html` and use final marker lines `IMPLEMENT:` or `APPROVED:`. Treat approval-only markers as workflow permission, not evidence that every interpretation was endorsed. Before acting on new feedback, snapshot it unchanged under `feedback_history/`.

## Direct Feedback Intake

Before launching any drafting subagent or generating figures, create `feedback_intake.md`.

Find all user feedback that bears on the target figure or subpanel, not only the current package's ideation feedback. Search the current redraft, earlier planning/feedback files in the same redraft, sibling figure-generation packages, and prior manuscript-version feedback or polishing/integration roots. For redraft figure-generation roots, also scan sibling `*/ideation/user-feedback.txt` files for cross-package instructions before deciding scope.

In `feedback_intake.md`, record:

- every direct feedback source read, with relative path and whether it contains substantive instructions, approval-only permission, abort/subsume instructions, investigation requests, or text/integration-only comments
- package scope decisions such as `draft`, `abort`, `subsume_into_<package>`, `defer_to_integration`, `defer_to_results_text`, or `blocked_pending_investigation`
- concrete user directives that must appear in subagent briefs, preserving the user's intent even when the wording is informal or self-correcting; give each directive a stable `directive_id`
- cross-version and cross-package implications, especially when feedback in one package or manuscript version changes another package's drafting scope
- each directive's planned drafted output, report section, responsible subagent, or explicit non-drafting disposition
- each directive's final report status: `addressed`, `partially_addressed`, `blocked`, `dropped`, or `deferred`, plus the output path, blocker explanation, or non-drafting disposition that justifies that status

Direct user feedback takes precedence over ideation report interpretations. If user feedback contradicts a candidate idea or decision checklist, follow the user feedback and document the disposition.

When reviewing or critiquing draft outputs, explicitly answer:

- Have the draft figure subpanels changed from the last user-reviewed version in any way that was not directly requested or clearly implied by user feedback? Treat unrequested changes as unacceptable unless the panel is explicitly marked exploratory or novel.
- Have the draft figures failed to address any direct user-feedback point, no matter how small? Treat omissions as unacceptable unless a blocker or explicit non-drafting disposition is documented.

## Prior Code Fidelity Gate

For every affected prior figure/subpanel, classify the drafting path before plotting:

- `inherited_preserve`: prior reviewed panel should be regenerated unchanged except for local-path relocation.
- `inherited_targeted_fix`: prior reviewed panel should be modified only to satisfy explicit user feedback.
- `inherited_move`: prior reviewed panel is reused in a different figure/subpanel position with its substantive content preserved.
- `inherited_replace`: prior reviewed panel is intentionally replaced because direct user feedback or the approved redraft plan requires it.
- `novel_no_prior`: no prior reviewed panel/code exists or the user explicitly requested a novel panel.
- `blocked_prior_missing`: the prior reviewed output or code should exist but cannot be identified.

Write `prior_code_fidelity.csv` before generating reviewable PNGs. Use these columns unless a package-specific reason requires an additional column:

```text
panel_id,inheritance_class,prior_png_path,prior_code_path,copied_baseline_code_path,
active_local_script_path,diff_path,allowed_change_directive_ids,fidelity_status,blocker
```

For `inherited_preserve`, `inherited_targeted_fix`, and `inherited_move`, prior reviewed code must be copied into the local drafting package before adaptation, usually under `prior_code/<panel_id>/`. The active local script under `scripts/` must then be copied/adapted from that local baseline, not recreated from scratch. Generate a plain-text diff under `code_diffs/<panel_id>.diff` comparing the copied baseline code to the active local script. Each substantive diff must be justified by a direct user-feedback directive id, a purely mechanical local-path/output-path relocation, or an explicit blocker.

For `inherited_replace`, record the prior reviewed PNG/code when identifiable and the directive ids or plan text that justify replacement. For `novel_no_prior`, explain why no prior-code copying is applicable. For `blocked_prior_missing`, stop treating the affected panel as ready: mark the blocker in `feedback_intake.md`, `prior_code_fidelity.csv`, `drafting_panels.md`, and the top blocker section of `review_report.html`.

Do not source, import, or assemble directly from old manuscript draft scripts, polishing scripts, prior draft PNGs, or ad-hoc figure-generation folders as a shortcut around this gate. Those files may be inspected and copied into the local package, but the local drafting package must contain the reproducible source used for every presented inherited panel.

## Subagent Brief Contract

Each drafting subagent must receive a short written brief that includes:

- the exact source paths for the relevant ideation files, direct feedback files, planning source, and prior panel material
- the inheritance class for each assigned panel from Prior Code Fidelity Gate
- the feedback excerpts or directive summary from `feedback_intake.md` that govern the assigned subpanel(s)
- the specific drafted outputs expected, including any requested A/B variants, broad option generation, investigation task, or abort/defer decision
- the local-script requirement from Local Regeneration And Prior-Version Fidelity and the expected output folder
- for inherited panels, the requirement to start from copied local prior code, generate a diff against the active local script, and map every substantive code change to a directive id, mechanical relocation, or blocker
- a request for a returned feedback-coverage statement mapping each assigned directive to generated output, skipped output, blocker, or caveat

The main drafter must integrate these returned feedback-coverage statements into `feedback_intake.md`, `drafting_panels.md`, and the corresponding `review_report.html` sections. If a directive remains blocked or only partially addressed, the blocker must be visible in the report's top directive-status table and not hidden only in an appendix.

## Final Figure Legends

Every populated `final_figures/<variant>/` folder must contain a `legend.md` review guide. For each assembled PNG, name the source script/subpanels, panel meanings, relevant feedback basis, caveats, and meaningful alternatives.

## Integrated Reviewer Report

After each initial drafting pass or feedback revision, write or update `review_report.html` as the sole human-review surface.

The report must be a static self-contained HTML file with every reviewable drafted PNG visually present, not merely linked. Present figures in narrative order with the relevant feedback history, caveats, source files, and visual-QC notes near the images.

The report must be organized as the review surface, not as a file dump. Use this order:

1. Decision summary and reviewer action needed.
2. Directive-to-output status table. Each direct user-feedback item needs a stable id, source, short directive, status, output or blocker disposition, and the first report section where it can be inspected.
3. Recommended final drafts.
4. Blockers and missing required outputs, if any.
5. Conservative alternatives.
6. Exploratory alternatives.
7. Refined and initial subpanels, collapsed by default unless a subpanel is itself the primary review item.
8. Full provenance appendix.

The report should let the reviewer answer "what should I approve, what changed because of my feedback, and what is still blocked?" without reading the appendix.

Include a prior-code fidelity table near the top report summary, before the raw provenance appendix. For each inherited, replaced, novel, or blocked panel, show the inheritance class, prior reviewed output/code source, copied local baseline path when applicable, active local script, diff artifact, allowed directive ids, fidelity status, and blocker if present.

Use `.agents/skills/manuscript-figure-drafting/scripts/review_report_template.R` as the starting template for report generation. Copy it into the local drafting package, usually as `scripts/make_<package>_review_report.R`, adapt the CONFIG block, and run it through `scripts/agentRrunner.sh`. Keep the generated HTML self-contained and regenerate it from the local script after figure or legend changes. A local report script may adapt the template, but it must preserve the required directive-status table, recommended-first order, visible blocker section, report-coverage checks, and collapsed raw galleries. Do not replace the template with a generic image gallery plus raw appendix.

Before generating the report, write `report_manifest.csv` for every PNG under `initial_subpanels/`, `refined_subpanels/`, and `final_figures/`. The report generation script must fail or flag a blocker if a reviewable PNG is absent from the report, if a direct feedback directive lacks a status-table row, if a `blocked` or `partially_addressed` directive is absent from the top report summary, if an inherited panel lacks a `prior_code_fidelity.csv` row/copy/diff without a blocker, or if recommended final drafts omit required panels without an explicit blocker/non-drafting disposition. Use labeled contact sheets only when every option remains visibly reviewable.

## Workflow

1. Read `docs/project_map.txt`.
2. Read `.agents/references/manuscript_figure_style.md`.
3. Identify the relevant planning source, ideation files, all direct user feedback, and the most recent user-reviewed or polished version of every affected panel.
4. Write `feedback_intake.md` before scope decisions are treated as settled. Explicitly classify abort, subsume, defer-to-text, defer-to-integration, investigation, preservation, and targeted-fix directives.
5. Write `prior_panel_disposition.csv` and `prior_code_fidelity.csv`; classify each affected panel as inherited, replaced, novel, or blocked.
6. For inherited redrafting work, copy the most recent accepted panel-generation code into `prior_code/<panel_id>/`, adapt the copied code into the active local script under `scripts/`, and generate `code_diffs/<panel_id>.diff`. Treat unrequested changes to data, filtering, layout, encodings, panel role, or numbering as likely unacceptable.
7. Use subagents for scoped drafting tasks when available. Each subagent should receive the governing feedback excerpts, prior-panel source paths, inheritance class, expected outputs, and a request for directive-by-directive coverage.
8. Draft only the options needed to answer direct feedback and approved exploratory ideas. Record skipped options in `not_drafted.md`.
9. Write local generation scripts, notes, manifests, legends, `drafting_panels.md`, `feedback.md`, and `review_report.html`.
10. Visually inspect final assembled PNGs and compare them against the last user-reviewed versions and every direct feedback point.
11. Review the generated `review_report.html` itself. Confirm that recommended final drafts appear before conservative/exploratory/raw galleries, each direct feedback item is traceable to an output or blocker, the prior-code fidelity table is visible near the top, and no required draft is hidden only in appendix text or a manifest.

For feedback revision, read the feedback snapshot, preserve it under `feedback_history/` if it is not already archived, and revise existing drafting outputs in place. Update `drafting_panels.md` so it records each feedback/redrafting round, what changed, source files, commands, and remaining caveats.

## Local Regeneration And Prior-Version Fidelity

Presented panels must be generated by scripts in the local drafting folder, usually under `scripts/`. For redrafting work, first find the most recent polished or otherwise user-reviewed version of the figure/subpanel and the code that produced it. Use that code as the starting point by copying/adapting it locally.

Fresh local generation is a reproducibility requirement, not permission to redesign. Any change to data source, filtering, statistical summary, visual encoding, layout, panel order, panel role, or interpretation that is not directly requested or clearly implied by primary user feedback is likely unacceptable. Mark genuinely novel or exploratory panels explicitly and keep them separate from targeted-fix recommendations.

For inherited panels, the copied baseline code, active local script, and diff artifact are part of the required deliverable. A panel with identifiable prior reviewed code is not ready for review if the local package contains only newly written code, only a sourced dependency on the old script, or only a copied PNG. When copying a large prior script that generated multiple panels, either copy the full script and identify the relevant panel-producing block in `prior_code_fidelity.csv`, or copy a minimal faithful extraction with a note explaining the extraction boundary.

Use existing repo material in these ways:

- Encourage inspection of existing code from across the repo to understand available analyses, plotting approaches, and data structures.
- Encourage direct use of shared code under `R/`.
- Encourage direct use of datasets and analysis outputs under `data/`.
- When useful code exists in a prior manuscript draft, polishing folder, work-package folder, or other ad-hoc location, copy the relevant logic into a local drafting script and adapt it there instead of sourcing or assembling from that location.
- Avoid direct dependencies on pre-existing manuscript draft scripts, polishing scripts, prior draft PNGs, or ad-hoc figure-generation locations.

Local drafting scripts should be the reproducible source for each presented figure. Their dependencies should generally be limited to `data/`, `R/`, standard package imports, and small local helper functions copied into the drafting folder when needed.

## Constraints

- Run R scripts through `scripts/agentRrunner.sh`.
- Do not start expensive refits, segmentation reruns, classifier retraining, or long jobs unless explicitly approved.
- Visually inspect final assembled PNGs before completion. Existence/readability checks are insufficient. In `drafting_panels.md`, explicitly assess the visual-QC items from `.agents/references/manuscript_figure_style.md`, plus whether any summary display is inappropriate where raw data was requested. If a final figure has a serious review-visible defect, generate an alternate or mark the defect as a blocker.
