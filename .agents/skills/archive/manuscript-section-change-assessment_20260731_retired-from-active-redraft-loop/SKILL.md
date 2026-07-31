---
name: archived-manuscript-section-change-assessment
description: "Archived historical copy of manuscript-section-change-assessment, retired from the active redraft loop on 2026-07-31. Use only when explicitly auditing or comparing the former section-change assessment workflow."
---

# Manuscript Section Change Assessment

## Purpose

Use this skill after manuscript figures, data products, legends, semantic interpretations, feedback, or claim/evidence graph outputs have been regenerated or revised and a prior manuscript draft is available.

The goal is to determine, section by section, whether the existing manuscript text clearly needs updating to accommodate the current manuscript package. This skill returns the prior section text with comments/annotations on passages that appear stale, unsupported, contradicted, overstrong, missing necessary caveats, or mismatched to the current figures/data.

This skill does not build a manuscript spine, create section contracts, decide a new narrative order, write replacement Results prose, or state how the manuscript should change. It only marks places where prior text appears to demand an update and briefly explains why.

## Project Defaults

- Read `docs/project_map.txt` before assessment work.
- Prefer the newest explicitly requested manuscript integration root. If no root is specified, infer the latest coherent root under `agent-dev/manuscript_integration/`.
- Treat the prior manuscript draft or prior section sidecars as mandatory input. If no prior draft text exists, stop and report the blocker.
- Treat current figure-set outputs, current legends, current semantic interpretations, current claim graph, and current feedback-preflight outputs as context for evaluating whether the prior text is still compatible with the manuscript package.
- Keep reads narrow. Do not inspect raw figure-generation history, raw user feedback files, or analysis outputs unless required context is missing from the current integration outputs or feedback-preflight packet.
- Err toward `no_change` unless a concrete piece of prior text is clearly stale, contradicted, overstrong, unsupported, or missing a necessary qualification under the current manuscript package.

## Required Inputs

Establish these before launching section subagents:

- **Prior manuscript draft text:** rendered manuscript HTML, manuscript Markdown, section sidecars, or another source from which prior section text can be extracted. This is required.
- **Integration root or current manuscript package:** current figure-set integration and claim-graph integration outputs, or explicit equivalent paths.
- **Current figure set:** `figure_set_manifest.csv`, `figure_numbering_crosswalk.csv`, final figure paths, legends or legend paths, and accepted exceptions when present.
- **Current semantic interpretations:** `semantic_interpretation_index.csv` and the panel interpretation files relevant to each section.
- **Current claim graph:** `claim_graph_integrated.json`.
- **Claim/evidence audit:** `claim_evidence_update_table.csv`, `unresolved_claim_decisions.md`, and orphan/retired-evidence reports when present.
- **Feedback preflight packet:** a current `manuscript-feedback-preflight` output, including the feedback-sources list, file-level feedback audits, and any section/figure/claim routing summaries. If absent, run or request that preflight assembly before this skill proceeds.
- **Current figure legends, if applicable:** required for sections that discuss, cite, or depend on figures.
- **Output root:** usually a `section_change_assessment/` subfolder inside the current integration root unless the user requests another path.

## Subagent Requirement

Section change-assessment tasks must use subagents. Launch one subagent per prior manuscript section in scope, not per proposed future section.

If no callable subagent or multi-agent facility is available, stop before producing assessment outputs and ask the user for permission to proceed without the subagent requirement.

Each section subagent should receive only:

- the prior section heading, section type, and full prior section text;
- current legend text or legend paths for figures relevant to that section, if applicable;
- feedback-preflight excerpts routed to that section;
- relevant claim graph slice;
- relevant rows from claim/evidence audit and unresolved decisions;
- relevant figure/evidence rows, figure callout crosswalk rows, and accepted exceptions;
- relevant semantic interpretation files or excerpts for panels used by the section;
- any explicit user scope constraints for the assessment.

Each section subagent must return:

- the prior section text with comments/annotations inserted near passages that clearly demand an update;
- a brief rationale for each annotation, tied to current figures/data/legends/feedback/claims where possible;
- an overall assessment label: `no_change`, `targeted_changes`, `extensive_changes`, or `complete_redraft`;
- a short assessment rationale;
- blockers or missing context, if any.

Subagents must not write replacement prose, propose a new section structure, prescribe exact manuscript edits, or strengthen/soften claims beyond identifying why existing text may no longer fit.

## Annotation Standard

Annotated section outputs must preserve the prior section text as the main body. Annotations should be concise and local.

Use this inline format immediately after the affected sentence, clause, paragraph, heading, or figure callout:

```markdown
[ASSESSMENT: update_needed | reason: <brief rationale tied to current manuscript package>]
```

If a whole paragraph or section appears compatible with the current package, do not annotate it. If the whole section requires reassessment, add a short section-level note before the prior text and still mark any especially clear local triggers.

Allowed reasons include stale figure number, changed panel content, claim/evidence mismatch, overstrong wording, missing caveat, contradicted interpretation, feedback conflict, removed evidence, new required context, or unclear support.

Do not include replacement wording after an annotation.

## Workflow

1. **Establish scope**
   - Identify current integration root, prior draft source, feedback-preflight root, output root, and manuscript sections to assess.
   - Default to all Results sections plus any non-Results sections explicitly affected by revised figures, legends, feedback, or claim-graph changes.
   - Block if prior section text cannot be extracted.

2. **Verify preflight feedback assembly**
   - Confirm that a current `manuscript-feedback-preflight` packet exists for the same manuscript generation.
   - If the packet is absent, stale, or incompatible with the selected integration root, run or delegate preflight assembly before section assessment.
   - Route feedback excerpts to prior sections conservatively. Preserve file/source provenance in the routing table.

3. **Load current manuscript context**
   - Read current figure-set, legend, semantic-interpretation, claim-graph, and claim/evidence audit outputs.
   - Record unresolved scientific or narrative decisions as constraints, not as instructions to rewrite.
   - Do not invent new claims or new section purposes.

4. **Prepare section packets**
   - Extract prior section text and existing figure callouts.
   - Map prior figure callouts to current figures/panels using the crosswalk and manifest when possible.
   - Attach current legends and semantic interpretations for figures/panels used by the section.
   - Attach claim graph slices and feedback excerpts only when they are relevant to the prior section text.
   - Write or stage a `section_input_map.csv` before launching subagents.

5. **Launch section subagents**
   - Start one subagent per prior manuscript section in scope.
   - Ask each subagent to annotate the prior section text only where a clear update need exists.
   - Require conservative change-scope classification and rationale, not redraft instructions.

6. **Merge assessments**
   - Check that every section in scope has exactly one assessment output.
   - Normalize labels and reason categories without rewriting the subagents' substantive judgments.
   - Flag inconsistent or unsupported annotations for manual review rather than converting them into rewrite instructions.
   - Ensure section feedback, claim, legend, and semantic-interpretation sources are traceable.

7. **Classify change scope**
   - Use these labels:
     - `no_change`: no clear update needed.
     - `targeted_changes`: specific sentences, callouts, caveats, or small passages need attention, but the section mostly remains compatible.
     - `extensive_changes`: many local conflicts or changed dependencies make patching broad portions likely necessary.
     - `complete_redraft`: the prior section is broadly incompatible with the current manuscript package.
   - Assign labels by prior section. Do not assign a label by proposed future section.

8. **Produce assessment outputs**
   - Write the section input map, per-section annotated prior text files, summary change-scope table, feedback routing table, unresolved blocker list, and assessment report.
   - Keep rationale concise and evidence-linked. Avoid drafting, replacement language, or new manuscript architecture.

## Output Contract

Write outputs under the chosen output root:

```text
section_input_map.csv
section_assessments/
feedback_routing_by_section.csv
section_change_assessment.csv
unresolved_assessment_blockers.md
section_change_assessment_report.md
```

### section_input_map.csv

Required columns:

- `section_id`
- `section_type`: `abstract`, `introduction`, `results`, `discussion`, `methods`, or another explicit type.
- `prior_heading`
- `prior_text_path`
- `relevant_prior_callouts`
- `current_figure_callouts`
- `relevant_legend_paths`
- `relevant_semantic_interpretation_paths`
- `relevant_claim_ids`
- `relevant_evidence_ids`
- `relevant_feedback_ids`
- `notes`

### section_assessments/

Write one Markdown file per assessed prior section. Each file must include:

- section ID, prior heading, and prior text source;
- overall assessment label;
- short assessment rationale;
- blockers or missing context;
- annotated prior section text using the required annotation format.

### feedback_routing_by_section.csv

Required columns:

- `feedback_id`
- `feedback_source`
- `section_id`
- `prior_heading`
- `routing_basis`: `explicit_section`, `figure_callout`, `claim_topic`, `legend_topic`, `semantic_topic`, `global`, or `manual`.
- `feedback_summary`
- `directive_strength`: `must`, `should`, `could`, `question`, or `ambiguous`.
- `notes`

### section_change_assessment.csv

Required columns:

- `section_id`
- `section_type`
- `prior_heading`
- `assessment`: one of `no_change`, `targeted_changes`, `extensive_changes`, or `complete_redraft`.
- `annotation_count`
- `primary_update_drivers`
- `assessment_rationale`
- `blockers`
- `subagent_output_path`

## Validation Targets

A deterministic validation pass should be able to check:

- Required output files exist.
- A prior draft source is recorded.
- A compatible feedback-preflight packet is recorded.
- Every section in scope has one row in `section_input_map.csv`, one row in `section_change_assessment.csv`, and one Markdown file under `section_assessments/`.
- Every assessment label is one of `no_change`, `targeted_changes`, `extensive_changes`, or `complete_redraft`.
- Every annotation uses the required `[ASSESSMENT: update_needed | reason: ...]` format.
- Figure callouts, legend paths, semantic interpretation paths, claim IDs, evidence IDs, and feedback IDs are either resolvable or listed in `unresolved_assessment_blockers.md`.
- `section_change_assessment_report.md` records selected inputs, subagents used, feedback-preflight source, validation summary, assessment-label distribution, blockers, and downstream handoff notes.
- Outputs do not contain replacement Results prose, proposed section contracts, or a newly planned manuscript spine.

## Completion Standard

Finish only after prior draft text has been assessed for every section in scope, every section has an annotated prior-text output and an overall assessment label, preflight feedback has been routed by section, and unresolved blockers are documented.

The final response should list the output root, prior draft source, feedback-preflight source, number of assessed sections, assessment-label distribution, blockers, and whether the package is ready for a separate redrafting or rendering workflow.
