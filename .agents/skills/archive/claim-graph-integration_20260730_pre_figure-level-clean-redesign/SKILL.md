---
name: claim-graph-integration
description: "Update and validate a manuscript claim/evidence graph from a supplied current evidence package and an existing graph. Use when Codex needs to reconcile current evidence items with graph evidence nodes, update or create evidence nodes, preserve user-fixed claim contracts, identify claims affected by evidence changes, run bundled graph validation and plotting scripts, and produce claim/evidence audit outputs without owning figure generation, prose drafting, or broader manuscript coordination."
---

# Claim Graph Integration

## Purpose

Use this skill to refresh a claim/evidence graph from supplied current evidence
artifacts and an existing graph. The skill reconciles current evidence items
against graph evidence nodes, updates or creates evidence nodes, identifies
claims affected by each evidence item, and validates/renders the resulting
graph.

This skill treats current evidence artifacts as supplied public inputs. It
should not depend on producer-internal file layout, create a replacement
manifest for the evidence package, or prescribe how other consumers use the
refreshed graph. It may flag claim-support problems, but it must not silently
weaken, strengthen, remove, or rewrite claims to force agreement with evidence.

## Bundled Resources

- `references/claim_graph_spec.md`: graph schema and user-fixed claim contract.
  Read this before editing graph JSON, checking locked-claim rules, or resolving
  schema details.
- `scripts/validate_claim_graph.py`: graph validator. Run after graph edits.
- `scripts/plot_claim_graph.py`: graph visualization script. Run after graph or
  claim-edge edits unless the user explicitly declines visualization.
- `scripts/authenticate_claim_graph_inputs.py`: hash and compare only the exact
  evidence-package inputs and starting graph supplied to this skill.
- `references/input_authentication.md`: input-index, reuse, and audit-cache
  contract. Read before reusing a prior evidence-item audit or complete return.

Do not send agents back to the original graph-prototype directory for routine
schema context; use the bundled reference.

## Project Defaults

- Preserve any claim with `user_fixed: true` according to the bundled spec: do
  not change locked claim `text`, `tier`, or `status` unless the user explicitly
  requests that change.
- Keep reads narrow. Do not inspect manuscript drafts, raw user feedback, whole
  figure-generation history, or broad analysis outputs unless the user expands
  scope.
- Treat evidence-description files, figure/panel identifiers, provenance paths,
  legends, and accepted exceptions as optional evidence-package fields unless
  the user or supplied package contract marks them required.

## Required Inputs

Establish these before writing outputs:

- **Current evidence package:** one or more supplied files, tables, paths, or
  excerpts that define the current evidence items to reconcile. It should give
  each evidence item a stable current identifier plus enough source information
  to write or update a graph evidence node.
- **Evidence descriptions:** visual, semantic, quantitative, provenance, or
  legend/source descriptions for current evidence items when available. If an
  evidence item has no description, use only source fields that are explicit in
  the supplied package and record the limitation.
- **Starting claim graph:** old/current graph JSON.
- **Prior-linkage information:** optional mapping or lineage fields that help
  identify whether current evidence is stable, moved, changed, new, split,
  merged, removed, or replaced.
- **Output root:** the directory where refreshed graph and audit outputs should
  be written.
- **Prior successful authentication and return:** optional inputs used only for
  exact-input complete or item-level audit reuse.

If a current evidence item lacks enough description/source information to write
a defensible evidence node, mark it as blocked unless an accepted exception
explicitly says it is contextual, QC-only, intentionally not interpreted, or out
of scope.

## Subagent Requirement

Claim-graph integration tasks must use subagents for current evidence items that
lack an authenticated prior audit. First follow
`references/input_authentication.md`. Launch one subagent per cache-miss current
evidence item, not per prior graph evidence node. A valid complete-return reuse
requires no evidence-item subagents; an authenticated reused audit satisfies the
per-item audit requirement.

If any cache-miss item requires audit and no callable subagent or multi-agent
facility is available, stop before producing integration outputs and ask the
user for permission to proceed without the subagent requirement.

Each evidence-item subagent receives only:

- current evidence item row(s) or excerpts;
- linked evidence-description file(s) or excerpts, when available;
- relevant prior-linkage information, when available;
- the starting claim graph or a minimized excerpt/index of old evidence nodes
  and claims;
- the evidence-node guidelines from this skill, if needed.

Each subagent must perform exactly these tasks:

1. **Identify corresponding old evidence node(s).**
   - Match using explicit prior linkage first, then old evidence `source`, old
     evidence notes, source package, short content, provenance, figure/panel
     labels, and lineage status when those fields are available.
   - Classify the mapping as `same_node`, `moved_node`, `changed_node`,
     `split_from_old`, `merged_from_old`, `new_node`,
     `removed_or_replaced_prior_node`, or `no_prior_match`.

2. **Draft an updated or new evidence node.**
   - Preserve old evidence IDs when the current item is the same evidence after
     renumbering or movement.
   - Create a new evidence ID when the current item is genuinely new, split
     beyond one old node, or materially changed enough that preserving the old
     ID would be misleading.
   - Populate `kind`, `source`, and `notes` using the guidelines below.

3. **Identify claims modified by the evidence item.**
   - List claims previously linked to matching old evidence nodes.
   - List claims that should newly reference the evidence item.
   - List claims that should lose old evidence because the current evidence item
     changed, moved to a different claim role, or was removed/replaced.
   - List possible new or missing claims only as recommendations.
   - Do not rewrite claims or edit graph JSON.

The main agent merges subagent outputs, resolves collisions, updates graph JSON,
checks user-fixed claim rules, runs validation, and records unresolved
decisions.

## Evidence Node Guidelines

An evidence node is a compact, source-anchored statement of a key result,
pattern, or relationship observable in supplied evidence. It describes what the
evidence shows; it does not describe which claim the evidence supports or how
the graph should use it. Claim use belongs in claim evidence links,
claim-to-claim edges, and audit tables.

Use the graph schema fields unless the existing graph has an established
additional field.

- `id`: stable evidence identifier. Preserve when lineage is stable; create a
  normalized new ID when evidence is new or materially changed.
- `kind`: concise class such as `figure_panel`, `supplemental_panel`,
  `simulation_panel`, `model_diagnostic_panel`, `contextual_panel`, or another
  existing graph convention.
- `source`: normalized current source reference plus path or package pointer
  when available.
- `notes`: one short paragraph or compact bullets stating the key
  result/pattern/relationship shown by the evidence, plus only the
  source/provenance and caveats needed to understand that result.

Use evidence descriptions carefully:

- Use visual or semantic descriptions to ground the result in what the evidence
  visibly or explicitly shows.
- Use source clarification to refine what plotted or summarized values
  represent: data identity, encodings, units, transforms, denominators, and
  quantified checks when available.
- Use provenance to distinguish observed data, derived summaries, model
  diagnostics, posterior summaries, or simulations.
- Preserve caveats and not-assessed notes when they constrain what the evidence
  node may state.
- Do not infer evidence-node content from filenames, package names, legends
  alone, prior feedback, or manuscript prose.

Evidence-node notes should distinguish:

- result/pattern/relationship visible or explicit in the evidence;
- clarified plotted-data or source identity;
- provenance/data or model source;
- limitations, ambiguity, or not-assessed boundaries.

Do not include language such as "supports claim X", "qualifies claim Y",
"primary evidence for...", or "used to argue..." inside the evidence node. Put
that information in the evidence-item audit and `claim_evidence_update_table.csv`.

## Workflow

1. **Establish scope**
   - Identify current evidence package, starting graph, output root, current
     evidence item definition, and whether graph edits are allowed or only
     recommendations are requested.
   - Read `references/claim_graph_spec.md` only if graph edits, locked-claim
     rules, or schema interpretation are needed.

2. **Validate starting graph**
   - Run `scripts/validate_claim_graph.py <starting_graph>`.
   - If validation fails, either repair schema-level issues first or report a
     blocker.

3. **Authenticate exact inputs**
   - Build the minimal current-evidence-item index already needed by this
     workflow from the supplied evidence package. Do not discover or require
     producer artifacts outside the declared Claim Graph inputs.
   - Snapshot the starting graph, declared evidence-package paths, and canonical
     item payloads with `scripts/authenticate_claim_graph_inputs.py`.
   - Compare a prior passing `input_authentication.json` when supplied. Follow
     only the reported `reuse_complete`, `reuse_partial`, or `fresh_full` path.
   - Treat a missing or checksum-mismatched prior audit as a cache miss.

4. **Build minimal indexes**
   - Index current evidence items from the supplied package.
   - Index linked evidence descriptions and prior-linkage fields when supplied.
   - Index old evidence nodes by id, source, notes, source references, and
     linked claims.
   - Index claims by id, text, tier, status, user-fixed flag, evidence links,
     methods, assumptions, support edges, qualification edges, and undermining
     edges.

5. **Launch evidence-item subagents**
   - Start one subagent per current evidence item not covered by an authenticated
     reuse record.
   - Require each subagent to produce a structured evidence-item audit with
     old-node mapping, proposed evidence node, affected claims, confidence,
     blockers, and unresolved questions.
   - Save audits under `evidence_item_audits/` when producing files.

6. **Merge evidence-node updates**
   - Resolve duplicate proposed IDs, split/merge conflicts, and old-node
     collisions.
   - Preserve old evidence IDs for stable lineage.
   - Create new nodes for new or materially changed evidence.
   - Retire or leave unreferenced old evidence nodes only with an explicit audit
     note; do not silently delete evidence.

7. **Update claim evidence references**
   - Apply low-risk reference updates where the evidence item maps cleanly and
     claim meaning is unchanged.
   - Add or remove evidence links only when the evidence-item audit and global
     graph context support the change.
   - Preserve user-fixed claim text, tier, and status.
   - Add or adjust support/qualification edges only when the semantic meaning is
     unchanged or clearly required by the existing graph logic.

8. **Produce decision reports**
   - Flag claims whose current evidence is missing, weak, contextual only,
     contradictory, or materially changed.
   - Flag current evidence items that do not clearly attach to any claim.
   - Flag old evidence nodes that are removed, replaced, or no longer
     represented by the current evidence package.
   - Keep substantive claim wording/tier/status changes as unresolved decisions
     unless explicitly authorized.

9. **Validate, render, and finalize authentication**
   - Write `claim_graph_integrated.json`.
   - Run:

```bash
python3 .agents/skills/claim-graph-integration/scripts/validate_claim_graph.py \
  <output_root>/claim_graph_integrated.json
```

   - After graph or edge edits, run:

```bash
python3 .agents/skills/claim-graph-integration/scripts/plot_claim_graph.py \
  --graph <output_root>/claim_graph_integrated.json \
  --output <output_root>/claim_graph_integrated.png \
  --dot-output <output_root>/claim_graph_integrated.dot
```

   - If Graphviz is unavailable, keep the DOT output when possible and record
     the PNG rendering blocker.
   - Write `audit_index.json`, then final `input_authentication.json` including
     the integrated result-graph hash, authenticated audit receipts, and
     `input_reuse_report.json`. Deterministic graph validation is mandatory even
     for `reuse_complete`.

## Output Contract

Write outputs under the chosen output root:

```text
claim_graph_integrated.json
claim_graph_integrated.dot
claim_graph_integrated.png
claim_graph_validation_report.txt
claim_graph_integration_report.md
evidence_item_audits/
evidence_remapping.csv
evidence_node_updates.csv
claim_evidence_update_table.csv
orphan_or_unassigned_evidence.csv
retired_or_unrepresented_old_evidence.csv
unresolved_claim_decisions.md
current_evidence_item_index.json
audit_index.json
input_authentication.json
input_reuse_report.json
```

### evidence_remapping.csv

Required columns:

- `current_evidence_item_id`
- `current_source_ref`
- `old_evidence_ids`
- `mapping_status`
- `prior_linkage_status`
- `confidence`
- `mapping_evidence`
- `notes`

### evidence_node_updates.csv

Required columns:

- `evidence_id`
- `action`: `preserved`, `updated`, `created`, `retained_unreferenced`, or
  `retired_in_report_only`.
- `old_evidence_ids`
- `kind`
- `source`
- `evidence_description_path`
- `provenance_path`
- `notes_summary`
- `accepted_exception`

### claim_evidence_update_table.csv

Required columns:

- `claim_id`
- `claim_tier`
- `user_fixed`
- `old_evidence_ids`
- `new_evidence_ids`
- `change_type`: `unchanged`, `remapped`, `added`, `removed`,
  `needs_decision`, or `not_assessed`.
- `support_classification`: `direct`, `partial`, `contextual`, `QC`,
  `qualifying`, `contradictory`, `missing`, or `uncertain`.
- `recommended_action`
- `blocker`
- `notes`

## Validation Targets

A deterministic validation pass should be able to check:

- Required files exist.
- Starting and integrated graphs pass `scripts/validate_claim_graph.py`, or
  validation failures are recorded as blockers.
- User-fixed claim `text`, `tier`, and `status` match the graph contract unless
  explicitly authorized.
- Every current evidence item has an evidence-item audit.
- Every reused evidence-item audit has a valid prior path and SHA-256 recorded in
  `input_reuse_report.json`; every other item has a fresh subagent audit.
- `input_authentication.json` covers the exact declared current evidence-package
  inputs, starting graph, canonical item index, and integrated result graph.
- Every created/updated evidence node has `id`, `kind`, `source`, and `notes`.
- Every claim evidence reference points to an existing evidence node.
- Every current evidence item is mapped to a graph evidence node or listed in
  `orphan_or_unassigned_evidence.csv`.
- Every old evidence node made obsolete by current evidence-package lineage is
  represented in `retired_or_unrepresented_old_evidence.csv` or preserved with
  rationale.
- Plotting produced `.dot` and `.png`, or the rendering blocker is recorded.
- `claim_graph_integration_report.md` records selected inputs, authentication and
  reuse disposition, subagents used, graph edits made, validation commands and
  results, unresolved decisions, blockers, and accepted exceptions.

## Completion Standard

Finish only after every current evidence item has a fresh or authenticated reused
audit, current items are reconciled to graph evidence nodes, claim evidence links
are updated or flagged, user-fixed claim rules are preserved, graph validation
has run, graph visualization has been attempted, and unresolved scientific or
narrative decisions are clearly reported.

The final response should list the output root, integrated graph path, validator
status, plot status, number of evidence items audited, graph edits made, and
unresolved blockers or decisions.
