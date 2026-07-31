---
name: plan_redraft
description: "Plan the next manuscript redraft from a completed passing manuscript-integration output and user feedback. 
Use when Codex needs to intake a passing integration root or future integration contract, preserve raw feedback intent and 
directive strength through mandatory section-level subagent review, assign every task to Planner, Major analysis, Figure drafting workflow, or 
Manuscript integration buckets, create a new redraft root, and produce both internal planning artifacts and a single human-friendly HTML review checkpoint without changing integration outputs."
---

# Plan Redraft

## Purpose

Use this skill to turn a completed manuscript integration pass plus user feedback into a tentative plan for the next manuscript redraft. 
This is a planning and review adapter: it normalizes the integration handoff, preserves the user's concrete intent, defines executable work packages, 
and presents the proposed plan in a readable HTML checkpoint.

This skill is downstream of `.agents/skills/manuscript-integration/`. It must not regenerate figures, edit the claim graph, rewrite the manuscript, run major analyses, or change integration outputs unless the user explicitly expands scope.

## Core Rule

The main failure mode to avoid is a dense or abstract plan that the user approves because it is too painful to review. 
Treat the user-facing review experience as a primary deliverable, not a byproduct of the agent plan.

Never convert a high-confidence user directive into a vague default, optional question, or deferred item merely because the work is important. 
If the user clearly asked for an analysis, figure revision, or framing change, plan it as an action awaiting normal approval and scope control. 
Ask a question only when the user's intent, target artifact, or acceptable scope is genuinely unresolved.

Every generated redraft plan is tentative until the user explicitly approves it after reviewing the HTML checkpoint.
Do not treat `redraft_plan.md`, package tables, or prior approval text as permission to execute.

The workflow shape is:

```text
intake -> feedback fidelity -> manuscript structure -> feedback-understanding subagents -> planning subagents -> presentation subagents -> review-surface lint -> user HTML review -> approval or revision -> later execution
```

## Feedback Fidelity Rules

Raw feedback is authoritative. If raw transcripts and summarized feedback both exist for the same review session, 
read the raw transcript as the primary source and use summaries only as navigation aids. If only summaries exist, proceed with them and record that limitation.

When raw transcripts and summaries both exist, write:

```text
<redraft_root>/feedback_understanding/raw_summary_concordance.md
```

This file must list:

- feedback points present in both raw and summary sources;
- points weakened by summarization, such as "do this" becoming "consider";
- concrete examples, figure callouts, or wording proposals present in raw feedback but absent from summaries;
- uncertainty, frustration, praise, or preserve constraints that affect planning priority;
- any summary statement that contradicts or over-smooths the raw transcript.

Use the concordance to correct the plan before package decisions are finalized. Do not let a lossy summary override a direct raw-transcript instruction.

## Directive And Decision Rules

Classify each important feedback point before planning:

- `directive`: the user clearly asked for the change or analysis.
- `preference`: the user indicated a favored direction but left scope open.
- `question`: the user asked whether something is true or worth checking.
- `concern`: the user flagged a risk without specifying a remedy.
- `preserve`: the user liked or approved something that must survive revisions.

For `directive` with high or medium confidence, set the default action to execute after plan approval unless a feasibility audit shows it is impossible or much larger than the user likely understood. For `preference`, plan the favored path and identify any scope choice. For `question` or `concern`, plan the minimum audit needed to answer it.

Before asking the user to decide, check whether the decision was already settled in project documents, prior redraft plans, integration notes, known-issues files, or the raw feedback itself. If the decision is settled, state it as a settled premise and plan the necessary implementation or audit. Do not re-ask settled decisions as if they were open choices.

## Mandatory Subagent Rule

This skill requires subagent delegation. If subagent tooling is not available, ask the user to enable or authorize subagents. If subagents remain unavailable, abort normal redraft planning and write only a concise blocker note explaining that deep feedback interpretation cannot be done to standard.

Use three mandatory delegation rounds:

1. **Feedback understanding.** Launch one subagent per manuscript section or conceptual section. Each subagent must deeply interpret the relevant manuscript region, assets, and user feedback.
2. **Detailed planning.** Launch planning subagents to convert the interpreted feedback into work packages, feedback items, dependencies, and output contracts.
3. **User-presentation mapping.** Launch presentation subagents to convert the work packages into plain-English review statements for the user.

Do not skip a round because the main agent believes the answer is obvious. The planner may consolidate, correct, and arbitrate subagent outputs, but the core interpretation and planning passes must be delegated.

## Manuscript Structure Discovery

Before launching subagents, inspect the integration source and infer a review structure:

- If the manuscript is relatively developed, use manuscript sections such as Abstract, Introduction, Results section 1, Results section 2, Methods, Discussion, and figure-specific supplements.
- If the manuscript is still nascent, use conceptually coherent figure/analysis groups.
- Prefer 4-8 sections. Split a section if it contains unrelated scientific claims or asset families. Merge tiny sections if separate subagents would not improve understanding.
- Record the structure in `<redraft_root>/manuscript_section_map.md` with section title, manuscript location, figures/assets, claims, and candidate feedback scope.

The section map is an internal planning artifact. It should help subagents orient themselves; it is not the primary user review surface.

## Feedback Understanding Round

For each manuscript/conceptual section, assign a subagent this task:

```text
Read the assigned manuscript section and all relevant assets deeply. Inspect the user feedback transcript and identify every feedback point that may apply to this section, including direct requests, concerns, praise/preserve comments, uncertainty, contradictions, and implied scientific worries.

Do not plan implementation yet. Focus on understanding what the user means and why they likely mean it.

Return a structured Markdown or JSON object with:
- section_title
- assets_inspected
- current_manuscript_role: what this section is trying to accomplish
- feedback_points: one object per relevant point, with:
  - user_said: concise paraphrase, preserving intent
  - likely_reason: why the user probably raised it
  - affected_artifacts: manuscript/figure/analysis/prose targets
  - polarity: change | preserve | question | concern | context
  - confidence: high | medium | low
  - directive_type: directive | preference | question | concern | preserve
  - directive_strength: explicit | strong_inference | tentative | unclear
  - source_anchor: raw file path plus nearby quote/paraphrase location if available
  - specific_examples: concrete figure panels, conditions, phrases, or user-proposed wording
  - ambiguity_or_risk: none | ambiguous | contradictory | possibly_ill_advised | scientifically_consequential
  - simple_change: true/false
  - needs_major_analysis_consideration: true/false
  - settled_decision_candidate: true/false, with evidence if true
  - default_action: planned_execute | planned_audit | planned_defer | needs_user_scope
  - preserve_constraint: what should not be broken, if any
  - user_clarification_needed: plain-English question if needed
- section_level_interpretation: the deeper issue tying the feedback together
```

Subagents must inspect relevant figures, legends, integration reports, claim graph entries, manuscript text, and maintained analysis outputs as needed. They should not rely only on the feedback transcript. They must not edit files or run major analyses.

The main planner must merge these outputs into `<redraft_root>/feedback_understanding/section_feedback_synthesis.md`. Preserve uncertainties and disagreements rather than smoothing them away.

## Detailed Planning Round

After feedback is understood, launch planning subagents to turn interpreted feedback into work packages. Assign packages by coherent scientific or operational scope, not by convenience. Good package boundaries share a manuscript question, evidence source, execution bucket, target skill, and expected reviewer decision.

Every planned task must be assigned to exactly one execution bucket:

- `Planner`: active plan-state work, routing decisions, dependency tracking, and small checks/audits/computations that do not rise to major analysis and are needed to complete the plan or decide whether to hand work to another bucket.
- `Major analysis`: substantial analysis packages that need the major-analysis workflow before downstream figure or manuscript work.
- `Figure drafting workflow`: figure ideation, figure-specific data extraction or checks, draft PNG generation, figure provenance, visual QC, legends tied to draft figures, and polishing-gate handoff.
- `Manuscript integration`: claim/evidence reconciliation, manifest and numbering checks, manuscript-wide prose/rendering integration, final legends, and integration validation after figures and analyses stabilize.

Do not create standalone execution families such as `AUD`, `TXT`, or `INT` in the approved plan. If a task is an audit, classify it by who should execute it: `Planner` for small routing or evidence-readiness checks that the planner will complete during planning, `Figure drafting workflow` for figure-local checks/provenance that naturally belong inside a figure package, `Manuscript integration` for checks that belong in the integrated manuscript state, or `Major analysis` when the audit requires substantial new computation or method development.

The planner must not defer any task allocated to `Planner`. When planning or revising a redraft, the planner may use subagents to perform light checks, analyses, or audits needed to complete the plan, but those planner-owned checks must be resolved in the current planning pass and recorded in the plan. If the work cannot be completed in planning, assign it instead to `Figure drafting workflow`, `Major analysis`, or `Manuscript integration` with an explicit output contract. Prefer delegation to the appropriate figure drafting package when the check is figure-local or panel-provenance work. Assign to manuscript integration when the issue can be validated from the final manifest, claim graph, rendered manuscript, or integration report without blocking earlier work.

Planning subagents return:

```text
- package_title
- execution_bucket: Planner | Major analysis | Figure drafting workflow | Manuscript integration
- plain_english_goal
- feedback_points_addressed
- preserve_constraints
- blocked_by
- blocks
- can_run_in_parallel_with
- output_contract
- target_skill_or_workflow
- likely_compute_class
- major_analysis_decision: required | candidate | not_major_analysis
- planner_latitude: none | small_check_completed | delegate_to_figure_workflow | assign_to_integration_workflow
- planner_completion_status: completed_during_planning | delegated_to_workflow | not_applicable
- default_action_basis: why the proposed default preserves the directive strength
- prior_decision_check: settled | not_settled | not_applicable, with source if settled
- user_decision_needed: none | approval_only | true_unresolved_scope | true_unresolved_science
- risks_if_done_wrong
```

The main planner writes the internal package plan to `<redraft_root>/redraft_plan.md`. This file may use internal package identifiers and operational detail, but it is not the user-facing review checkpoint.

## User-Presentation Mapping Round

After work packages exist, launch presentation subagents to decide how to explain the plan to the user. Their task is to produce plain-English mappings from feedback to actions.

Rules for user-facing statements:

- Avoid internal identifiers such as `P1`, `MA1`, `FIG_WP4`, `MI2`, claim IDs, or workflow acronyms in the main presentation text.
- Avoid planning jargon such as "crosswalk", "scope matrix", "package", "artifact contract", "no-SUM", and unexplained model-set labels unless the user used those exact terms and they are necessary.
- Use simple statements of the form: "You said A, B, and C, so we are going to do Y."
- Aggregate repeated feedback that leads to the same action. Do not repeat "you said X so we will do Y" several times for the same Y.
- If the user explicitly asked for Y, write "we are going to do Y after approval" rather than "should we do Y?"
- When a real uncertainty remains, add: "We are not sure whether you meant Z" or "This may be scientifically risky because..."
- Preserve positive feedback explicitly: "You liked X, so we will keep that design feature while changing Y."
- Surface only true unresolved decisions early. Do not pad the first pass with settled actions.
- Use manuscript/figure locations the user can recognize, not internal package labels.
- Include enough concrete detail that the user can recognize their feedback without reopening the raw transcript, especially named panels, visible patterns, example phrases, and proposed analyses.

The main planner writes the consolidated mapping to `<redraft_root>/user_review_mapping.md` and uses it to build the HTML checkpoint.

## Major Analysis Definition

Use "major analysis package" for work that is likely to be time-consuming, compute-expensive, scientifically consequential, or operationally complex enough to need its own decision before figure generation. Examples include resegmenting many images, running new optimizations, developing new models, running NUTS, large simulation batches, or reworking a maintained model-fitting pipeline.

Do not classify ordinary figure-generation support as major analysis. Examples that usually belong under `Planner` or `Figure drafting workflow` include comparing existing NUTS QC panels, choosing a nonredundant QC visualization, legend placement, recomputing or exposing a small existing-export summary, checking a panel's provenance, and deciding whether an existing figure should be split into two panels.

If a lightweight computation from existing exports is scientifically decisive but quick, classify it as `Planner` only when it is needed to route the work before assigning downstream ownership. Classify it as `Figure drafting workflow` when it is figure-local support for a panel or legend. Classify it as `Manuscript integration` when it can be checked from final integration state. Promote it to `Major analysis` only if it needs new model fitting, large batch execution, new data generation, or substantial method development.

For each required or candidate major analysis package, make a tentative package-level decision:

- `tentatively_execute`: run after user approval.
- `tentatively_defer`: do not run in the next redraft unless the user overrides.
- `needs_user_scope`: the scientific question, acceptable runtime, or output target is ambiguous enough that user input is needed before choosing.

Use `tentatively_execute` when the user clearly requested the analysis and the target is specific, even if the work is scientifically consequential. Use `needs_user_scope` only when the request itself is ambiguous. Use `tentatively_defer` only when the user did not request the work, when it is clearly out of scope, or when a feasibility audit shows it is not currently practical. Do not use the major-analysis gate as a quiet deferral mechanism for work the user explicitly requested.

Planning itself must not launch major analyses. If a major-analysis candidate needs feasibility work, the detailed planning round should place the feasibility check in the `Planner` bucket only when the planner can complete it during the current planning pass by inspecting existing outputs, scripts, configs, manifests, and run records without executing the analysis. Otherwise, delegate the check to `Figure drafting workflow`, `Major analysis`, or `Manuscript integration` rather than leaving it as deferred planner work.

## Figure Generation Rule

Treat figure generation as a substantial workflow, not a cleanup bucket. Plan it as ideation or panel-choice review, drafting, polishing, legend/provenance updates, and visual QC.

When local skills are available, structure figure-generation packages so they can later be executed by `.agents/skills/manuscript-figure-ideation/`, `.agents/skills/manuscript-figure-drafting/`, and then the polishing-gate workflow if final polished outputs are needed. For each figure package, name the intended root, review files, approval markers, and output contracts in the internal plan.

Evidence checks that gate figure choices must be assigned either to `Planner` when they are small routing checks completed during planning, or to `Figure drafting workflow` when they are figure-local data/provenance checks that naturally belong inside the figure package. Do not hide a blocking check in prose without an owning bucket.

## Inputs

Required inputs:

- A completed passing integration source:
  - Preferred future form: an explicit integration contract if one exists.
  - Current compatibility form: an integration root such as `agent-dev/manuscript_integration/<run_id>/`.
- User feedback on that integrated manuscript spine or figure set.
  - Prefer raw feedback transcripts stored directly in the integration root.
  - If raw transcripts and cleaned summaries both exist, use the raw transcript as primary and record all feedback files in the intake contract.
  - If the user explicitly supplies feedback outside the integration root, proceed only if the intent is clear and record this as an external-feedback exception in the intake contract.
  - If no feedback path is supplied, search only the integration root for likely feedback files such as `raw_transcript`, `transcript`, `feedback.md`, `user_feedback.md`, `redraft_feedback.md`, `review_feedback.md`, `manuscript_feedback.md`, or `*feedback*.txt`.

If multiple plausible feedback files exist and their relationship is unclear, ask the user which one to use. Do not infer scientific intent from stale or unrelated notes. If the relationship is clear, such as raw transcript plus cleaned summary for the same session, use all of them with the raw file primary.

## Passing Integration Compatibility Check

When no formal integration contract exists, treat an integration root as a passing integration handoff only if these checks pass:

- Required files exist:
  - `figure_set_manifest.csv`
  - `claim_graph_integrated.json`
  - `integration_report.md`
  - `final_images/`
  - exactly one manuscript-mode output by default: `manuscript_bullet_draft.md`, `manuscript_update_recommendations.md`, or a generated manuscript HTML draft.
- `integration_report.md` records no blockers, or records only user-accepted exceptions.
- `final_images/` contains one whole-figure PNG for each unique `figure_id` in `figure_set_manifest.csv`.
- The graph validates, when the project validator is available:

```bash
python3 agent-dev/manuscript_claim_graph/validate_claim_graph.py <integration_root>/claim_graph_integrated.json
```

- The report's listed validation commands have either passed in the current session or are recorded as previously passing in the integration report.

If any check fails, do not create a normal redraft plan unless the user explicitly accepts the risk. Instead, stop with a concise blocker list, or create a clearly labeled blocker note only if the user asks for an artifact.

## Redraft Root

Create a new root for the next redraft under:

```text
agent-dev/manuscript_redrafts/<run_id>/
```

Default `run_id` rules:

- If the integration run id contains a version such as `v3`, infer the next redraft as `v4` unless the user specifies another version.
- Include the current date: for example `20260514_v4_redraft`.
- If the target exists, add a short numeric suffix rather than overwriting it.

Keep plan-redraft outputs inside this root. Do not copy large figures or data by default; reference integration outputs by path.

## Revising An Existing Plan

When the user reviews a tentative plan and asks for plan revisions, keep working in the current redraft root unless the user explicitly asks for a new root. Before overwriting the active planning documents, archive the prior review state in a subfolder of the current redraft root, for example:

```text
<redraft_root>/plan_revision_archive/<revision_id>/
```

At minimum, copy the old active `redraft_review.html`, old active `redraft_plan.md`, and the user feedback file or transcript that prompted the revision into that archive folder. Then update the active planning documents in place according to the user feedback. A separate delta artifact is not required unless the user asks for one.

## Output Contract

Write this file early:

```text
<redraft_root>/redraft_intake_contract.json
```

Use schema version `0.4` or later. Include at least:

```json
{
  "contract_type": "plan_redraft_intake",
  "schema_version": "0.4",
  "created_at": "ISO-8601 timestamp",
  "redraft_root": "agent-dev/manuscript_redrafts/<run_id>",
  "integration_source": {
    "type": "integration_root_compat",
    "path": "agent-dev/manuscript_integration/<run_id>",
    "run_id": "<run_id>",
    "required_artifacts_present": true,
    "blockers": "none",
    "claim_graph_validation": "passed|previously_reported_passed|not_run_with_reason",
    "manuscript_mode": "html_draft|bullet_draft|update_recommendations"
  },
  "integration_artifacts": {
    "figure_set_manifest": "<path>",
    "claim_graph": "<path>",
    "integration_report": "<path>",
    "final_images_dir": "<path>",
    "manuscript_source": "<path>",
    "figure_legends": "<path or NA>"
  },
  "feedback": {
    "primary_path": "<raw transcript, summary, or single feedback file>",
    "primary_type": "raw_transcript|summary|single_feedback_file",
    "all_feedback_paths": ["<path>"],
    "colocated_with_integration_root": true,
    "sha256": {"<path>": "sha256:<hash>"},
    "raw_summary_concordance_required": false,
    "external_feedback_exception": ""
  },
  "subagent_rounds": {
    "feedback_understanding": "completed|required_but_blocked",
    "detailed_planning": "completed|required_but_blocked",
    "user_presentation_mapping": "completed|required_but_blocked"
  },
  "review_artifacts": {
    "human_review_html": "agent-dev/manuscript_redrafts/<run_id>/redraft_review.html",
    "internal_plan": "agent-dev/manuscript_redrafts/<run_id>/redraft_plan.md",
    "raw_summary_concordance": "<path or NA>",
    "user_review_lint": "agent-dev/manuscript_redrafts/<run_id>/user_review_lint.md",
    "user_feedback_transcript": "agent-dev/manuscript_redrafts/<run_id>/plan_feedback_transcript.md",
    "review_status": "agent-dev/manuscript_redrafts/<run_id>/plan_review_status.json"
  },
  "plan_review": {
    "status": "tentative_pending_user_review",
    "approval_required_before_execution": true
  },
  "project_map_decision": ""
}
```

Detailed planning belongs in `redraft_plan.md`, not in the JSON contract.

## Internal Plan

Write:

```text
<redraft_root>/redraft_plan.md
```

This file is for execution agents. Keep it complete but do not optimize it as the user review surface.

Required sections:

- `Scope`: integration source, feedback source, inferred next version, and non-goals.
- `Input Checks`: compatibility checks run, validation commands, artifact status, and exceptions.
- `Manuscript Section Map`: section structure used for subagent delegation.
- `Subagent Rounds`: agents launched, assigned sections/packages, files inspected, and output paths.
- `Feedback Understanding Synthesis`: concise synthesis by section, including preserve constraints and ambiguous or risky feedback.
- `Feedback Fidelity Check`: raw/summary concordance, directive-strength preservation, and settled-decision checks.
- `Execution Bucket Assignment`: every planned task assigned to exactly one of `Planner`, `Major analysis`, `Figure drafting workflow`, or `Manuscript integration`, with the rationale for that assignment.
- `Major Analysis Gate`: required/candidate major-analysis packages, tentative decisions, blocked downstream outputs, and feasibility status.
- `Planner Work`: small routing checks/audits/computations the planner completed during planning or revision, plus checks explicitly delegated away from the planner. This section must not contain future deferred planner tasks.
- `Figure Drafting Workflow Packages`: ideation/drafting/polishing packages, figure-local checks, dependencies, roots, review files, and output contracts.
- `Manuscript Integration Packages`: prose/legend/claim-strength work, claim graph, manifest, figure numbering, render readiness, and final integration checks.
- `Package Dependency Map`: what can run in parallel and what is blocked.
- `Evidence Constraints`: claims, caveats, figure dependencies, and source limitations that must survive the redraft.
- `Execution Exit Criteria`: conditions for considering the redraft ready for another integration pass.

Use internal package identifiers here as needed. Do not ask the user to review this file as the main checkpoint.

## Human Review HTML

Write:

```text
<redraft_root>/redraft_review.html
```

This is the primary user-facing checkpoint. It must be a single self-contained HTML document with concise prose, recognizable manuscript locations, and clear review order. Optimize it for a user who may provide an unstructured feedback transcript.

The HTML must tell the user what to look at and in what order. If a prior integrated manuscript HTML with embedded figures exists, treat it as the only required companion artifact and link to it prominently. Do not require the user to open raw package tables, JSON, claim graphs, feedback registers, `redraft_plan.md`, or feedback-synthesis files. If internal files are linked, put them in a clearly optional technical appendix.

Required HTML structure:

1. **Review Instructions**
   - State that the plan is tentative.
   - Tell the user which manuscript/draft artifact to keep open.
   - Tell the user they may respond with free-form notes; no table editing is required.
   - Give a short suggested response pattern: "approve", "approve except...", "I meant...", "do not do...".
2. **First Pass: Decisions And Risks**
   - Cover only true unresolved choices, blocking risks, or major analyses that need explicit approval despite a clear user directive.
   - For each item, use plain English:
     - "Look at: <manuscript section or figure>."
     - "You said: <aggregated user feedback>."
     - "We plan to: <action>."
     - "We are unsure whether: <question>." when relevant.
     - "Please tell us: <decision needed>."
   - If the user already gave the decision, do not ask it again; state the action and any scope/approval needed.
3. **Second Pass: Complete Feedback Coverage**
   - Walk through the manuscript in user-recognizable order.
   - Use aggregated statements: "You said A, B, and C, so we will do Y."
   - Include preserve constraints: "You liked X, so we will keep X while changing Y."
   - Include concrete examples from feedback, such as named panels, cell lines, visible patterns, or wording requests, when abstraction would make the plan hard to verify.
   - Mention simple obvious fixes briefly, but do not replace specific feedback with generic bullets.
4. **What Will Happen After Approval**
   - Summarize execution phases in plain English.
   - Avoid internal identifiers in the main text.
   - Include optional collapsible technical details only if useful.
5. **How Approval Will Be Interpreted**
   - Explain that the next planner pass will read the user's free-form transcript and set the status to approved or requiring more feedback.
   - State that unresolved first-pass questions block execution.

Do not make the HTML a giant table. Use short sections, bullets, and grouped cards. The user should be able to understand the plan without knowing package IDs.

Before finalizing the HTML, write:

```text
<redraft_root>/user_review_lint.md
```

The lint must pass these checks:

- No internal package IDs, claim IDs, workflow acronyms, source roots, or file paths appear in the main explanatory text, except links in a clearly labeled technical appendix.
- No high-confidence directive is presented as "should we do this?" or default-deferred.
- Every first-pass question is a true unresolved question with a stated reason.
- At least the important second-pass sections contain user-recognizable specifics, not only generic "improve" or "clarify" bullets.
- Any settled prior decision is described as settled and tied to implementation work, not re-opened for approval.
- The approval language distinguishes normal approval to execute planned work from separate approval for genuinely ambiguous scope.

## Feedback Transcript And Review Status

Also create:

```text
<redraft_root>/plan_feedback_transcript.md
<redraft_root>/plan_review_status.json
```

`plan_feedback_transcript.md` is a blank or lightly prompted place where the user can paste free-form feedback after reading the HTML. It must not require structured tables.

`plan_review_status.json` starts as:

```json
{
  "status": "tentative_pending_user_review",
  "human_review_html": "<redraft_root>/redraft_review.html",
  "feedback_transcript": "<redraft_root>/plan_feedback_transcript.md",
  "blocking_questions": [],
  "approved_at": null,
  "requires_more_user_feedback": true
}
```

When the user later provides review feedback, the planning agent must parse the transcript and update this status:

- `approved_for_execution`: the user clearly approved the plan and no blocking first-pass questions remain.
- `requires_plan_revision`: the user requested changes that can be incorporated into a revised plan/review HTML.
- `requires_more_user_feedback`: ambiguity remains, or the user response does not resolve blocking decisions.

If status is not `approved_for_execution`, do not hand off to execution as though the plan is approved.

## Workflow

1. Read `docs/project_map.txt`.
2. Identify the integration source and feedback file.
3. Run the passing integration compatibility check.
4. Create the new redraft root.
5. Write the initial `redraft_intake_contract.json`.
6. Read the integration report, manuscript-mode output, manifest, claim graph, legends if present, and all relevant feedback, with raw transcripts primary.
7. If raw and summarized feedback both exist, write `feedback_understanding/raw_summary_concordance.md`.
8. Infer the manuscript/conceptual section map.
9. Launch the mandatory feedback-understanding subagents, one per section.
10. Consolidate their outputs into `feedback_understanding/section_feedback_synthesis.md`.
11. Check directive strength and prior settled decisions before package planning.
12. Launch the mandatory detailed-planning subagents.
13. Consolidate packages into `redraft_plan.md`.
14. Launch the mandatory user-presentation mapping subagents.
15. Write `user_review_mapping.md`.
16. Write `redraft_review.html`, optimized for human review.
17. Write `user_review_lint.md`; revise the mapping/HTML until it passes.
18. Write `plan_feedback_transcript.md` and `plan_review_status.json`.
19. Update `redraft_intake_contract.json` with actual artifact paths and subagent-round statuses.
20. Decide whether `docs/project_map.txt` needs a terse update. Update it only when the new redraft root becomes the current maintained manuscript-planning state, not for every scratch plan.

## Completion Standard

Finish only after the redraft root contains:

- `redraft_intake_contract.json`
- `manuscript_section_map.md`
- `feedback_understanding/raw_summary_concordance.md` when raw and summarized feedback both exist
- `feedback_understanding/section_feedback_synthesis.md`
- `redraft_plan.md`
- `user_review_mapping.md`
- `redraft_review.html`
- `user_review_lint.md`
- `plan_feedback_transcript.md`
- `plan_review_status.json`

Or finish after clearly reporting why the integration/feedback inputs or mandatory subagent delegation are insufficient.

The final response should list the redraft root, feedback file used, validation status, subagent-round status, human review HTML path, plan-review status, and whether `docs/project_map.txt` was updated.
