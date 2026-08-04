# Figure Ideation Workflow

## Purpose

Create a reviewer-friendly ideation package that shows what is being revised,
which prompt, plan, or supplied feedback context is driving the work, and what
figure options could be generated next. The central output is
`ideation_report.html`; `ideas.md` is the concise source/sidecar for the
candidate options.

Assume figure-design judgment is fallible. Preserve supplied revision context,
show existing visual context, and generate a bounded
option set for important figure decisions. Do not draft figures, polish panels,
run major analyses, edit production code, or mark the proposal approved.

## Required Inputs

- A work-package or redraft figure-package id.
- A planning source, usually `docs/manuscript_plan_v3.txt`, `docs/manuscript_plan.txt`,
  or a redraft root containing `redraft_plan.md`.
- For revisions, the prompt or supplied feedback context and the existing
  figure or panel set being revised.
- Paths to existing integrated figures, polished figures, draft panels, prior
  subpanels, manifests, or legacy feedback files when they can be located.

## Output Location

Use the most local review root implied by the request:

- Legacy work-package workflow:
  `<work_package_root>/<WP_ID>/ideation/`
- Redraft figure-package workflow:
  `<redraft_root>/figure_generation/<package_id>/ideation/`

Write:

- `ideation_report.html`: primary self-contained review report with embedded
  existing figure/subpanel images when available, prompt or supplied feedback
  context, candidate ideas, and decision checklist.
- `ideas.md`: concise source for candidate outputs, rationale, preservation and
  change targets, and remaining choices.
- `feedback_manager_context.md`: optional pointer when a managed-feedback
  handoff was supplied.
- `existing_panel_disposition.csv`: required for revision or mixed work; mark
  each relevant existing figure/subpanel `preserve`, `targeted_fix`,
  `move_or_duplicate`, `replace`, `drop`, `uncertain`, or `not_applicable`.

## Workflow

1. Resolve the requested package, prompt, plan, data, and available figure/text
   context from supplied paths and local artifacts.
2. Inspect the existing figure or panel set being revised. Prefer rendered PNGs
   or integrated HTML over manifests alone.
3. Read only the planning and evidence files needed to understand the figure
   story, available inputs, and blocking dependencies.
4. Decide the scope:
   - `revision`: fix specific criticized parts while preserving the rest.
   - `greenfield`: propose new panels because no adequate figure exists.
   - `mixed`: combine a targeted revision with a clearly separated new evidence
     layer.
5. Write `existing_panel_disposition.csv` for revision or mixed work.
6. If a managed-feedback handoff was supplied, write
   `feedback_manager_context.md` with its concise pointer.
7. Write `ideas.md` with concrete alternatives for each important design
   decision. Group options so the user is not facing a flat catalog.
8. Write `ideation_report.html` as the main review surface.

## `ideation_report.html` Standard

The report must be a static self-contained HTML file. Embed relevant existing
PNGs directly when practical. If a source figure cannot be embedded cheaply,
include its path and a short explanation.

Use this report order:

1. Existing Figures And Panels.
2. Relevant Feedback or Prompt Context, grouped by figure, panel, or decision.
3. Panel Disposition for revision or mixed work.
4. Candidate Ideas tied to feedback and existing context.
5. Decision Checklist listing only choices needed before drafting.
6. Appendix with input context, source inventory, caveats, and
   deferred feedback questions.

Avoid making the report a directory index or a list of paths requiring
cross-reading.

## `ideas.md` Standard

Write `ideas.md` as a concise idea source for the HTML report. Start with the
figure decision being brainstormed and the feedback driving it. Translate
workflow ids, claim ids, model names, and package labels into plain language
immediately.

For revision work, preserve figure parts the feedback liked or left alone while
brainstorming alternatives for the criticized or missing evidence role. For
subpanel ideation, offer breadth within the affected role. For figure-level
ideation, offer a few coherent layout concepts rather than many interchangeable
parts.

For conditional analyses or unfinished outputs, reserve a slot without letting it
dominate the proposal. Keep detailed evidence audits in appendices or side notes
and summarize only decision-relevant conclusions in `ideas.md`.
