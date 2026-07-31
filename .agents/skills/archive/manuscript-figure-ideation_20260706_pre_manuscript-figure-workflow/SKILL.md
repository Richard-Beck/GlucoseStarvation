---
name: archived-manuscript-figure-ideation
description: "Archived historical copy of manuscript-figure-ideation, superseded by manuscript-figure-workflow on 2026-07-06. Use only when explicitly auditing or comparing the pre-promotion figure ideation workflow."
---

# Manuscript Figure Ideation

## Purpose

Create a reviewer-friendly ideation package that shows the user what is being revised, what the direct feedback said, and what figure options could be generated next. The package should help the user decide without forcing them to reconstruct context from scattered Markdown files, prior PNGs, and planning notes.

Assume the agent's figure-design judgment is unreliable. The skill compensates by preserving direct feedback, showing existing visual context, and generating a broad but bounded option set for each important figure decision. Do not collapse ideation into one overconfident design unless the user explicitly asks for a single recommended plan.

The central output is `ideation_report.html`, a self-contained review report. It should present, in this order:

1. Relevant existing figures and subpanels.
2. What the direct feedback said.
3. What the candidate ideas are and what choices remain.

Also write `ideas.md` as a concise source/sidecar for the idea options. It should answer, in plain language:

- What candidate outputs could be generated?
- Why is that the right next figure work?
- Which existing panels or ideas are starting points, liked context, or targets for revision?
- What choices does the user still need to make before drafting?

For revision work, direct user feedback is paramount. Raw feedback transcripts, user-authored notes, and feedback files outrank plan summaries and prior agent syntheses. Reread the relevant direct feedback yourself after any discovery step, interpret it explicitly, and ideate only the parts of the figure that the user appears to want changed. If the user liked or did not criticize a panel, preserve it unless the feedback clearly implies otherwise.

For greenfield work, start from the manuscript story and available evidence, then propose candidate panels that would make that story visible.

For mixed work, separate the targeted revision from genuinely new material. Do not let a small new evidence layer turn a targeted revision into a sprawling redesign.

Do not draft figures, polish panels, run major analyses, edit production code, or mark the proposal approved.

## Required Inputs

- A work-package or redraft figure-package id, such as `WP3` or `FG3_wp4_posterior`.
- The relevant planning source, usually `docs/manuscript_plan_v3.txt` or a redraft root containing `redraft_plan.md`.
- For revisions, the most relevant direct user feedback and the existing figure or panel set being revised.
- For revisions, paths to existing integrated figures, polished figures, draft panels, or prior subpanels when they can be located from the project map, manifest, plan, or feedback.

Use `docs/project_map.txt` first to orient to maintained roots and current workflow status.

## Output Location

Use the most local review root implied by the request:

- Legacy work-package workflow: `agent-dev/manuscript_work_packages/<WP_ID>/ideation/`.
- Redraft figure-package workflow: `<redraft_root>/figure_generation/<package_id>/ideation/`.

Write:

- `ideation_report.html`: the primary self-contained review report with embedded existing figure/subpanel images where available, relevant feedback, and candidate ideas.
- `ideas.md`: the concise idea source used by the report.
- `feedback_sources.md`: source paths and short notes for direct feedback, feedback summaries, plan files, and figure/manifest files inspected.
- `existing_panel_disposition.csv`: required for revision or mixed work; each relevant existing figure/subpanel marked `preserve`, `targeted_fix`, `move_or_duplicate`, `replace`, `drop`, `uncertain`, or `not_applicable`, with rationale.
- `feedback.md`: the user-editable feedback file.

## Feedback Contract

Create `feedback.md` if it does not exist. Include only this contract, adjusted for the package name:

```text
# Feedback for <package> ideation

Edit this file after reviewing ideation_report.html.
Primary review targets:
- ideation_report.html
- ideas.md
- feedback_sources.md
- existing_panel_disposition.csv, when present

When ready, append one final marker line:
- IMPLEMENT: revise ideation outputs from this feedback.
- APPROVED: permit drafting to proceed from this ideation package.

Do not leave either marker as the final nonblank line until you want the runner or agent to act.
```

Treat ideation outputs as tentative until `feedback.md` ends with `APPROVED` or the user explicitly approves the ideation in chat. An approval marker without substantive comments is permission to proceed, not proof that every idea is endorsed.

## Feedback Discovery

When subagents are available and active instructions permit their use, launch one feedback-scout subagent before writing ideas. Keep the task narrow:

- Identify direct user feedback relevant to the package, including raw transcripts, user-authored review notes, feedback files, and chat-derived notes.
- Return source paths, exact relevant sections or line ranges where possible, and a concise explanation of why each source matters.
- Distinguish direct feedback from agent summaries, planning records, and approval-only files.
- Do not propose figure designs or decide the package scope.

After the scout returns, the main agent must read the relevant direct feedback sources itself before writing `ideas.md` or `ideation_report.html`. If subagents are unavailable or not permitted, do the same discovery locally and record that no scout was used in `feedback_sources.md`.

## Workflow

1. Read `docs/project_map.txt`.
2. Run the feedback-discovery step above, then read the relevant direct feedback yourself. For revisions, direct feedback is more important than plans, registers, or previous agent summaries.
3. Inspect the existing figure or panel set being revised. Prefer the actual rendered PNGs or integrated HTML over manifests alone. Note what the feedback liked, disliked, questioned, or left alone.
4. Read only the planning and evidence files needed to understand the figure story, available inputs, and blocking dependencies.
5. Decide the real scope:
   - `revision`: fix specific criticized parts while preserving the rest.
   - `greenfield`: propose new panels because no adequate figure exists.
   - `mixed`: combine a targeted revision with a clearly separated new evidence layer.
6. Write `existing_panel_disposition.csv` for revision or mixed work.
7. Write `feedback_sources.md`, separating direct feedback from summaries/plans.
8. Write `ideas.md` as a concise source for what could be generated next. Give several concrete alternatives for each important design decision, but group them so the user is not facing a flat catalog.
9. Write `ideation_report.html` as the main review surface.
10. Write or refresh `feedback.md` without adding an approval marker.

## `ideation_report.html` Standard

The report must be a static self-contained HTML file. Embed relevant existing PNGs directly when practical; if a source figure is not a PNG or cannot be embedded cheaply, include a clear path and a short explanation instead of silently omitting it. Do not draft new data panels for ideation.

Use this report order:

1. **Existing Figures And Panels**: show the current figure/subpanel material being revised or used as a starting point, with short captions explaining why each item matters.
2. **Relevant Feedback**: summarize the direct feedback in reviewer-friendly language, grouped by figure/panel or decision. Preserve important wording from the user when it affects interpretation, and cite the source path.
3. **Panel Disposition**: for revision/mixed work, show what is preserved, fixed, moved, replaced, dropped, or uncertain.
4. **Candidate Ideas**: show the bounded option set from `ideas.md`, tied directly to feedback and existing visual context.
5. **Decision Checklist**: list only the decisions the user needs to make before drafting.
6. **Appendix**: include feedback-source inventory, plan/evidence files read, caveats, and any direct feedback that was intentionally deferred or narrowed.

The report should be the main artifact a user reviews. Avoid making it a directory index or a list of paths that requires cross-reading to understand the proposal.

## `ideas.md` Standard

Write `ideas.md` as the concise idea source for the HTML report. Avoid opaque IDs, workflow jargon, and internal bookkeeping. If claim IDs, feedback IDs, model names, or package labels matter, translate them immediately:

```markdown
The relevant feedback is about Figure 4's parameter-effect panel: uptake looks well supported, yield looks mixed, and the current figure does not show which models disagree.
```

Do not write:

```markdown
The plan maps this package to F14-F18 and C2/C6 and asks for posterior support audit outputs.
```

Start with candidate figure work, not the administrative scope. A useful `ideas.md` usually has:

- a brief statement of the figure decision being brainstormed
- a brief explanation of the user feedback driving it
- a preservation/change list for existing panels, summarized from `existing_panel_disposition.csv` when present
- a bounded option set for each panel or figure decision, each with "what it shows" and "why it helps"

Keep `ideas.md` focused on what could be generated, but do record deferred, out-of-scope, or intentionally narrowed user feedback in the report appendix and disposition files. If a deferred item changes how the user should judge the ideas, mention it briefly in `ideas.md`.

For revision work, preserve the figure parts that the feedback liked or left alone, but still brainstorm alternatives for the criticized or missing evidence role. Avoid turning one possible repair into a prescription. For example, do not say only "redraft the interval panel and emphasize support strength with labels." Instead, offer options such as:

- add text support labels
- encode support strength with color
- identify models using point shapes
- separate model disagreement into a side tile

For subpanel ideation, offer breadth within the criticized or missing evidence role. For figure-level ideation, offer a few coherent layout concepts rather than many interchangeable parts.

For conditional analyses or unfinished outputs, reserve a slot without letting it dominate the proposal:

```markdown
If the leave-one-line-out posterior summaries pass QC before drafting, one option is to add a compact robustness tile to the supplement. Another option is to fold the robustness result into the parameter-effect panel using model or subset markers.
```

Keep the file short enough to review in one sitting. If detailed evidence audits are needed, put them in a separate note or appendix and summarize only the decision-relevant conclusion in `ideas.md`.
