---
name: method-table-provenance
description: "Orient, validate, construct, classify, outline, revise, or draft Methods-facing provenance from manuscript figures. Use when Codex needs to trace current manuscript panels to Methods-relevant objects, maintain the canonical provenance table, or consume that table for Methods work."
---

# Method Table Provenance

Route Methods-facing provenance work from current manuscript panels. When
responding directly to user instructions or feedback, identify the coordinating
role as `router`.

## Required Input

Require a current manuscript figure or panel package with enough metadata to
identify every visible panel and trace it through generating artifacts and
Methods-relevant transformations. Manifests, provenance files, semantic notes,
scripts, configs, saved tables, and model outputs may supply those facts. Report
unsupported or contradictory edges instead of inventing them.

## Canonical Provenance Contract

Maintain one current-snapshot Markdown table with exactly these columns:

```text
id | parent | what | why | comment | lock_target | lock_selector | lock_kind | sha256 | hash_status
```

The table is a panel-rooted directed acyclic graph, not a historical archive:

- Use one structural leaf per current manifest panel, with an id such as
  `F2#panel_a`.
- Give every panel leaf its final whole-figure path and SHA256 plus its mandatory
  `#panel_<id>` selector. Panels in the same figure share the path and hash.
- Connect each panel leaf directly to its last Methods-relevant object or
  transformation. Do not represent rendering, styling, layout, composition,
  copying, or package assembly between them.
- Keep every row upstream-reachable from a current panel leaf. Exclude prior
  versions, overlays, disconnected components, and superseded nodes from the
  canonical table; preserve them outside the current bundle when history is
  needed.
- Use stable ids and current canonical repo-relative paths. Resolve duplicate
  rows that represent the same durable object or action.
- Allow repeated hashes for panels sharing a whole figure and for honest proxy,
  code, or representative locks. Repeated hashes alone do not make rows
  duplicates.

Materialize the current panel scope as `target_figure_set.tsv`:

```text
endpoint_id | artifact_path | lock_selector | sha256
```

Use one row per panel. Require `endpoint_id` to end with its selector, for
example `F2#panel_a` and `#panel_a`. The number of unique target hashes must
equal the number of current final figures.

## Canonical Methods Handoff

For a complete Methods-ready handoff, write one copy of each artifact under the
user-specified root:

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

Keep supporting fact cards, reviewer points, and unresolved-issue notes outside
the canonical consumer contract unless requested.

## Router Workflow

1. Identify the requested output and current panel package.
2. Materialize the panel-level `target_figure_set.tsv` from that package.
3. Run structural panel-leaf validation before reusing any existing provenance,
   classification, spine, or Methods text.
4. Run lock-target verification for the canonical table.
5. If either check fails, repair or reconstruct the affected current-snapshot
   provenance before propagating only necessary downstream changes.
6. Reuse downstream artifacts only when their upstream table passes both checks
   and their scope remains current.
7. Generate the minimum missing or invalid prerequisite chain and return the
   requested canonical artifact.

Use subagents only when a requested workflow genuinely requires distinct
tracing, writing, inspection, or review roles. Keep their assignments bounded
to the router-defined scope.

## Workflow Routing

Read only the relevant reference:

- `references/provenance_table_construction.md`: construct or replace the
  canonical current-snapshot table.
- `references/provenance_table_verification.md`: validate graph structure,
  panel leaves, and file locks.
- `references/node_classification.md`: classify nodes in a structurally valid
  canonical table.
- `references/methods_spine.md`: build a Methods outline from the validated
  table and classifications.
- `references/methods_drafting.md`: draft Methods text from the Methods spine.

Do not draft manuscript prose unless the user requests Methods drafting or a
routed workflow explicitly requires it.
