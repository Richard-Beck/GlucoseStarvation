---
name: method-table-provenance
description: "Build, verify, classify, analyze, outline, revise, or draft Methods-facing provenance from manuscript figures. Use when Codex must trace figure-generating scripts and analytical objects into a deduplicated recursive dependency graph, maintain a current locked provenance table with exact panel endpoints, or consume that table for Methods work."
---

# Method Table Provenance

Route Methods-facing work from manuscript figures through one current-snapshot
provenance graph.

## Inputs

Require a current figure or panel package with enough metadata to identify:

- every manuscript-visible panel;
- the current whole-figure artifact for each panel;
- figure-generating scripts or equivalent execution roots; and
- generating objects, inputs, manifests, or other trace evidence where known.

Accept manifests, package provenance, inventories, scripts, configs, saved
objects, legends, semantic notes, and explicit user corrections. Report gaps;
never invent provenance.

Existing outputs from this skill are optional inputs.

## Canonical Methods handoff

For a complete Methods-ready handoff, write singular current artifacts under a
user-specified or concise local root:

```text
<methods_root>/
  target_figure_set.tsv
  manuscript_endpoint_verification.md
  methods_text.md
  locked_provenance_table.md
  provenance_lock_verification.md
  node_classification.csv
  methods_spine.md
  final_methods_draft_audit.md
  methods_handoff.md
```

Keep construction evidence outside the canonical consumer contract unless the
user requests process reproduction:

```text
<methods_root>/
  figure_object_roots.tsv
  object_registry.tsv
  dependency_edges.tsv
  object_resolution_index.tsv
  object_fact_cards/
  unresolved_methods_issues.md
```

## Core construction policy

Keep the structural leaves of `locked_provenance_table.md` exactly equal to the
panels in `target_figure_set.tsv`. Keep every canonical row upstream-reachable
from those leaves.

Trace from deduplicated figure-generating scripts and generating objects, not
from independent per-panel chains. Preserve panel-to-object mappings so all
panel leaves can later attach to their last Methods-relevant object.

For each object, perform two adjacent operations:

1. Enumerate candidate direct scientific dependencies mechanically.
2. Resolve which candidates are actually consumed and what scientific
   transformation, if any, the object represents.

Do not postpone all semantic resolution until after a raw file graph is
complete. Do not treat every loaded file or sourced script as a scientific
parent. Exclude runtimes, packages, shell tools, and generic implementation
machinery, including project-local path/config readers and execution-only
status gates; this graph captures scientific lineage, not the software or
orchestration environment.

Use one router-owned object registry. Deduplicate before every worker wave.
Never reconcile by repeatedly merging prose chains.

Treat a prior locked table as:

- reusable when both lock and endpoint verification pass;
- a source of authenticated unaffected rows in a targeted update;
- a read-only diagnostic or later comparison baseline in a clean
  reconstruction; and
- never a source of unchecked inherited rows merely because those rows existed
  before.

## Provenance modes

Choose one mode before provenance construction:

- `verified_reuse`: both deterministic verifiers pass. Reuse the locked table
  without Object Resolvers. Reuse classifications, spine, and Methods prose only
  where their inputs remain authenticated.
- `targeted_update`: verification identifies a bounded changed scope. Retain
  only authenticated unaffected rows and queue new, changed, or unresolved
  objects on affected branches. Use Object Resolvers only for that queue.
- `clean_reconstruction`: the user requests a clean rebuild, the prior graph is
  unavailable, or invalidation is too broad to authenticate unaffected
  branches. Build from current figure/script/object roots and use one bounded
  Object Resolver per canonical unresolved object requiring semantic
  resolution. Materialize panel endpoints and other deterministic identities in
  the router.

Do not choose targeted update merely to reduce work. Escalate to clean
reconstruction when shared ancestors, aliases, or scope changes make unaffected
branch authentication ambiguous.

## Delegation

The coordinating role is `router`.

For `clean_reconstruction`, delegate one canonical unresolved object requiring
semantic resolution per Object Resolver. For `targeted_update`, apply the same
rule only to new, changed, or unresolved objects on affected branches.
`verified_reuse` uses no Object Resolvers. The router may materialize panel
endpoints, hashes, and other deterministic identities without delegation. Give
each resolver only its assigned object, bounded current evidence, and the
contract in `references/object_resolution.md`. A resolver:

- identifies only direct dependencies;
- resolves the local transformation and Methods role;
- does not recurse into newly discovered parents;
- does not edit shared registries or canonical outputs; and
- returns one fact card.

The router registers and deduplicates newly discovered objects before the next
wave. Bundle objects only when one function or action generates them
inseparably. Use bounded concurrency; default to at most eight simultaneous
workers unless the user specifies otherwise.

Use an efficient low-reasoning worker for clear extraction when available.
Escalate ambiguous dynamic paths, generator identity, or scientific
transformations to a stronger worker or the router.

Do not delegate deterministic manifest creation, hashing, registry validation,
canonical reconciliation, reachability pruning, or final endpoint validation.

## Router workflow

1. Read the request and define the requested canonical output.
2. Materialize one row per current panel in `target_figure_set.tsv`.
3. For a revision, run both deterministic verification workflows and select
   `verified_reuse`, `targeted_update`, or `clean_reconstruction`.
4. For targeted update, authenticate retained unaffected rows and seed the queue
   only with affected current roots and invalid/new objects. For clean
   reconstruction, materialize and deduplicate all current
   figure/script/object roots. Stop construction for verified reuse.
5. Initialize one object registry and direct-dependency edge registry for the
   selected construction scope.
6. Resolve queued objects in bounded worker waves until the registry reaches
   closure or explicit terminal/unresolved boundaries.
7. Reject fact cards that change assigned identity, invent schema values, or
   recurse through a newly discovered parent; reconcile valid cards centrally
   by canonical object identity and evidence.
8. Convert confirmed Methods-relevant objects and actions into one locked
   provenance table; omit unused and display-only dependencies.
9. Attach every panel leaf to its last Methods-relevant object.
10. Audit, prune unreachable rows, and run both bundled validators.
11. Propagate changes only to classifications, spine entries, or prose that the
    repaired provenance actually invalidates.

If Object Resolvers required by targeted update or clean reconstruction are
unavailable, stop and ask the user to enable them or narrow the construction
scope.

## Workflow routing

Read only the reference needed for the task:

- `references/provenance_table_verification.md`: verify locks and current panel
  endpoints.
- `references/provenance_table_construction.md`: build or revise the recursive
  object dependency graph and canonical locked table.
- `references/object_resolution.md`: prompt and output contract for one object
  worker.
- `references/node_classification.md`: classify nodes in a verified table.
- `references/graph_analysis.md`: inspect or reduce a verified graph.
- `references/methods_spine.md`: build a Methods outline from classified
  provenance.
- `references/methods_drafting.md`: draft Methods text from a verified spine.

Do not draft manuscript prose unless the user requests it.
