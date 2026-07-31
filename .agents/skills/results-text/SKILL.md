---
name: results-text
description: "Draft, revise, or coordinate clean, sharp manuscript Results prose from figure-description sources, with optional use of user-locked claim-graph claims. Use when Codex needs to write one or more Results sections from semantic figure interpretations or legends, delegate per-section Results drafting to subagents, or produce Results sidecar files for renderer injection."
---

# Results Text

## Purpose

Write journal-facing Results prose from figure-description evidence. The task is to discover the scientific argument already present in the assigned figure(s) and express it clearly while respecting locked claims.

This skill owns Results drafting, Results drafting delegation, and section-level Markdown sidecar files for renderer injection.

## Inputs

### Required

* A figure-description source describing the assigned figure(s) and subpanel(s). Prefer semantic interpretation files. Figure legends may substitute for or augment semantic interpretations when they describe the visible figure content accurately enough to support Results prose.

### Optional

* Claim graph excerpts for the section, especially locked claims.
* Approved prior Results text when the user or routing contract permits reuse.
* Figure manifest rows or legends selected by the coordinating agent when they
  are needed to identify panel order or visible terminology.

### Routing Metadata

* Section id/title, assigned output path, and target main/supplemental figures are routing metadata. They identify the drafting assignment but are not scientific evidence.

## Writer assignment packet

Before delegation, the coordinating agent must materialize one bounded packet
per Results section containing:

- section id, working title or target claim, assigned figures, and output path;
- the exact figure-description sources the writer may use;
- relevant locked-claim excerpts, when supplied;
- explicitly permitted prior Results text, when reuse is allowed;
- any selected manifest or legend excerpts needed for panel identity;
- the Results style guide; and
- user constraints or known source conflicts relevant to that section.

The packet is the writer's complete authorized context. Resolve paths and
context selection before delegation. Writers must not browse the repository,
search for additional evidence, inspect generating code or raw data, read other
manuscript sections, or substitute undeclared sources. If the packet is
insufficient, the writer must identify the missing fact in `Revision Notes`
rather than locate it independently.

## Style

Consult `references/results_style_exemplars.md` before drafting and again before finalizing. Treat it as the authority for Results style.

Match the guide’s emphasis on conclusion-first headings, direct rationale, evidence-linked claims, compact methodological context where needed, precise panel use, local caveats, and synthesis sentences that state what the evidence implies.

## Evidence Source Interpretation

Use the figure-description source as the primary scientific evidence.

Interpret sources as follows:

* **Semantic figure interpretations** are the preferred evidence source. Use them to determine what each figure and subpanel establishes.
* **Figure legends** may be used when they describe figure content, measurements, comparisons, conditions, or panel meaning. Do not treat generic legend boilerplate as a result.
* **Locked claims** are claim graph statements marked `user_fixed: true`. These represent user-approved scientific claims that should be preserved and stated as forcefully as the figure-description evidence allows.
* **Unlocked claim graph content** may help identify context, but it is not a priority for this skill.
* **Style exemplars** are style guidance only. They are never scientific evidence.

If sources conflict, or if an optional input asks for a change not supported by the figure-description source, preserve accurate figure-supported prose and flag the issue in `Revision Notes`.

## Locked Claims

When a locked claim is relevant to the assigned figures, try to state it directly and without unnecessary caveat. Prefer clear formulations such as:

* “These results show…”
* “These data indicate…”
* “Together, these findings support…”
* “This analysis identifies…”

Use softer formulations only when the figure-description evidence requires them, for example:

* “These data suggest…”
* “These results are consistent with…”
* “This pattern is compatible with…”

Do not dilute a locked claim reflexively. Caveat only when the figure-description source makes the limitation necessary for accurate interpretation.

## Delegation

Use subagent threads for Results writing when more than one Results section needs drafting or revision. Prefer one dedicated subagent per Results section, with a disjoint output file such as `manuscript_sections/results/<section_id>.md`.

Give each subagent only its writer assignment packet and the payloads explicitly
listed in that packet. Do not ask section writers to locate context.

Require each subagent to write the assigned sidecar file. The coordinating Results agent should inspect, merge only when necessary, and preserve sidecars as review artifacts.

If a callable subagent-thread tool is not available, use tool discovery if available. If subagents still cannot be made callable, ask the user whether to enable/use subagents or to proceed with an explicit non-delegated fallback. Do not proceed with main-agent solo multi-section Results drafting by default.

## Output File

Write the assigned sidecar file in Markdown. The completed file must have this structure:

```markdown
---
section_id: <id>
title: <conclusion-first title>
primary_claims: [<locked claim ids/text when provided, otherwise empty>]
main_figures: [<figure ids>]
supplemental_figures: [<figure ids>]
input_versions:
  figure_description: <semantic interpretations, legends, or both>
  claim_graph: <provided or not provided>
---

## Results Text

<Final journal-facing Results prose.>

## Revision Notes

1. <Brief note on unsupported locked claims, source conflicts, necessary caveats, or source limitations. If none, state "None.">
```

The renderer should inject `Results Text` by default. Keep `Revision Notes` as review/audit material unless explicitly requested for manuscript inclusion.

## Validation

Before finishing, confirm that:

* the Results prose is grounded in the figure-description source;
* `references/results_style_exemplars.md` was consulted before drafting and again before finalizing;
* relevant locked `user_fixed: true` claims are stated as forcefully as the evidence allows;
* unresolved source conflicts or unsupported locked claims are flagged in `Revision Notes`;
* the manuscript-facing text contains no internal workflow labels, file paths, generation bookkeeping, or audit language;
* panel citations support specific claims without forcing coverage of every assigned subpanel;
* delegated writers used only their declared assignment-packet inputs;
* the sidecar file exists at the assigned path and is ready for renderer injection.
